nextflow.enable.dsl = 2

include { PROCESS_READS } from './modules/2-process_reads.nf'
include { DECOMPRESS_READS } from './modules/2b-decompress.nf'
include { WRITE_LIBRARY_CONF } from './modules/2c-write_library_conf.nf'
include { MAKE_DATABASE_SRA } from './modules/3-make.database.sra.nf'
include { INDEX_BAITS_STOP } from './modules/3b-index_baits.nf'
include { MINING_SRA } from './modules/4-mining_sra.nf'
include { REFORMER_SRA } from './modules/5a-reformer_sra.nf'
include { EXONERATE_SRA } from './modules/5b-exonerate_sra.nf'
include { STITCHER_SRA }        from './modules/5c-ExAMiner_stitcher.sra.nf'
include { PREPARE_COORDINATES } from './modules/6a-write_exon_intron_list.nf'
include { MAPPER }               from './modules/6b-mapper.nf'
include { CONSENSUS }            from './modules/6c-consensus.nf'
include { PREPARE_ALIGNMENT_INPUT }            from './modules/7a-prepare_alignment_input.nf'
include { ALIGNER }            from './modules/7b-aligner.nf'
include { CLEANING } from './modules/8a-cleaning.nf'
include { FILTERING } from './modules/8b-filtering.nf'
include { PARTITIONFINDER } from './modules/9a-partitionfinder.nf'
include { MRBAYES }         from './modules/9b-mrbayes.nf'
include { GENERATE_DATASETS } from './modules/9c-generate_datasets.nf'

workflow {
    // ---------- sanity checks ----------

    def checkDir = { name, path ->
        def d = file(path)
        if (!d.isDirectory()) {
            error "Missing required directory '${name}': ${path}"
        }
    }

    def checkFile = { name, path ->
        def f = file(path)
        if (!f.isFile()) {
            error "Missing required file '${name}': ${path}"
        }
    }

    [
        infodir: params.infodir,
        datadir: params.datadir,
        baits_stop_dir: params.baits_stop_dir,
        baits_nostop_dir: params.baits_nostop_dir
    ].each { name, path ->
        checkDir(name, path)
    }

    [
        samples_list: params.samples_list,
        samples_info: params.samples_info,
        baits_list: params.baits_list,
        sra_conf: params.sra_conf
    ].each { name, path ->
        checkFile(name, path)
    }

    if (!file(params.baits_stop_dir).listFiles().any { it.name.endsWith('.fasta') }) {
        error "No FASTA files found in baits directory: ${params.baits_stop_dir}"
    }

    [
        samples_list: params.samples_list,
        samples_info: params.samples_info,
        baits_list: params.baits_list,
        sra_conf: params.sra_conf
    ].each { name, path ->
        if (!file(path).text.trim()) {
            error "Required input file is empty '${name}': ${path}"
        }
    }

    // Check sample and bait input files
    def expected_samples = file(params.samples_list).readLines()*.trim().findAll { it }
    def expected_baits_master = file(params.baits_list).readLines()*.trim().findAll { it }

    if (!expected_samples) {
        error "SAFEGUARD: params.samples_list (${params.samples_list}) resolved to zero sample IDs."
    }
    if (!expected_baits_master) {
        error "SAFEGUARD: params.baits_list (${params.baits_list}) resolved to zero bait IDs."
    }

    // ---------- samples ----------
    ch_samples = Channel
        .fromPath(params.samples_list)
        .splitText()
        .map { it.trim() }
        .filter { it }

    // ---------- STEP 2 - process reads ----------
	if (!params.skip_mapping || !params.skip_make_database_sra || !params.skip_mining_sra) {

		if (!params.skip_process_reads) {

			// Process raw reads
			ch_reads_raw = ch_samples.map { s ->
				tuple(
					s,
					file("${params.datadir}/${s}_1.fastq.gz", checkIfExists: true),
					file("${params.datadir}/${s}_2.fastq.gz", checkIfExists: true)
				)
			}

			PROCESS_READS(ch_reads_raw)

			ch_reads_proc = PROCESS_READS.out.reads

		} else {

			// Assume reads were already processed
			ch_reads_proc = ch_samples.map { s ->

				def r1 = file("${params.processed_reads_dir}/${s}_1.proc.fastq.gz")
				def r2 = file("${params.processed_reads_dir}/${s}_2.proc.fastq.gz")

				// fallback to datadir
				if (!r1.exists() || !r2.exists()) {
					r1 = file("${params.datadir}/${s}_1.proc.fastq.gz")
					r2 = file("${params.datadir}/${s}_2.proc.fastq.gz")
				}

				if (!r1.exists() || !r2.exists()) {
					error """
					Missing processed reads for sample: ${s}

					Looked in:

					1) ${params.processed_reads_dir}
						${s}_1.proc.fastq.gz
						${s}_2.proc.fastq.gz

					2) ${params.datadir}
						${s}_1.proc.fastq.gz
						${s}_2.proc.fastq.gz

					Either:
					- run with --skip_process_reads false to generate processed reads
					- provide existing processed reads in one of the locations above
					"""
				}

				tuple(s, r1, r2)
			}
		}

		// keep a copy of compressed reads for the mapper stage, before
		// ch_reads_proc gets overwritten by DECOMPRESS_READS below
		ch_reads_compressed = ch_reads_proc

		// Get insert sizes
		insert_sizes = file(params.samples_info)
			.readLines()
			.findAll { it.trim() }
			.collectEntries { line ->
				def fields = line.split(/\s+/)
				[(fields[0]): fields[1]]
			}
	}

    if (!params.skip_make_database_sra || !params.skip_mining_sra) {

        // ---------- 2b - decompress reads before SRA step ----------
        DECOMPRESS_READS(ch_reads_proc)
        ch_reads_proc = DECOMPRESS_READS.out.reads

        // ---------- 2c - library configuration files ----------

        ch_reads_for_conf = ch_reads_proc.map { sample, r1, r2 ->
            if (!insert_sizes.containsKey(sample)) {
                error """
                Missing insert size for sample:

                    ${sample}

                Not found in:
                    ${params.samples_info}
                """
            }
            tuple(sample, r1, r2, insert_sizes[sample])
        }

        WRITE_LIBRARY_CONF(ch_reads_for_conf)
        ch_library_conf = WRITE_LIBRARY_CONF.out.library_conf
    }

    // ---------- STEP 3 - make SRA database ----------

    if (!params.skip_make_database_sra) {

        ch_library_conf_for_mining = ch_library_conf.map { sample, r1, r2, conf ->
            tuple(sample, r1, r2, conf)
        }
      
        MAKE_DATABASE_SRA(ch_library_conf)

        ch_database = MAKE_DATABASE_SRA.out.database
            .join(ch_library_conf_for_mining, by: 0)


    } else if (!params.skip_mining_sra) {

        // Assume SRA databases already exist
        ch_database = ch_library_conf.map { sample, r1, r2, conf ->

            def db = file("${params.sra_db_dir}/${sample}")

            if (!db.exists()) {
                error """
                Missing SRA database for sample: ${sample}

                Expected directory:
                ${db}

                Either:
                1. run with --skip_make_database_sra false to build databases
                2. provide existing SRAssembler databases in:
                   ${params.sra_db_dir}
                """
            }

            tuple(sample, db, r1, r2, conf)
        }
    }

    // ---------- STEP 4 - mine SRA database ----------

    if (!params.skip_mining_sra) {

        INDEX_BAITS_STOP(Channel.value(file(params.baits_stop_dir)))
        ch_baits_stop_indexed = INDEX_BAITS_STOP.out.baits_stop_indexed

        ch_mining = ch_database
            .combine(ch_baits_stop_indexed)
            .map { sample, database, r1, r2, library_conf, baits_fasta, baits_fai, baits_id_map ->
                tuple(
                    sample, database, r1, r2,
                    file(params.baits_list),
                    baits_fasta, baits_fai, baits_id_map,
                    file(params.sra_conf),
                    library_conf
                )
            }

        MINING_SRA(ch_mining)

        ch_manifest = MINING_SRA.out.manifest

    } else if (!params.skip_reformat_sra) {

        log.info "Skipping SRA mining"

        def manifest = file("${params.sra_mining_dir}/manifest.tsv")

        if (!manifest.exists()) {
            error """
            Missing SRA mining manifest:

            ${manifest}

            Either:
            - run with --skip_mining_sra false to generate it
            - provide an existing MINING_SRA manifest
            """
        }

        // Reconstruct manifest channel from previous MINING_SRA output
        ch_manifest = Channel.fromPath(manifest)

    }


    // ---------- bait reference lookup (avoids staging entire 5000-file dirs per task) ----------
    ch_baits_stop_byid = Channel
        .fromPath("${params.baits_stop_dir}/*.fasta")
        .map { f ->
            def m = (f.name =~ /^.*_(.+)\.fasta$/)
            tuple(m[0][1], f)
        }

    ch_baits_nostop_byid = Channel
        .fromPath("${params.baits_nostop_dir}/*.fasta")
        .map { f ->
            def m = (f.name =~ /^.*_(.+)\.fasta$/)
            tuple(m[0][1], f)
        }


    // ---------- STEP 5 - reform per-bait fasta files ----------

    if (!params.skip_reformat_sra) {

        // manifest.tsv rows are "<bait>\t<sample>/<bait>/all_contigs.fasta",
        // written relative to the MINING_SRA task workdir the manifest itself
        // lives in — resolve against manifest_file's parent to get real paths.
        ch_bait_contigs = ch_manifest
            .flatMap { manifest_file ->
                    manifest_file.splitCsv(sep: '\t').collect { row ->
                        def bait = row[0]
                        def contig_file = manifest_file.getParent().resolve(row[1])
                        // row[1] = BMNH1045344/RZC33347.1/all_contigs.fasta
                        def sample = row[1].tokenize('/')[0]
                        tuple(bait, sample, file(contig_file))
                    }
                }
                .groupTuple(by: 0)   // -> tuple(bait, [samples], [contig_files])

        ch_reformer_ready = ch_bait_contigs
            .join(ch_baits_stop_byid)
            .join(ch_baits_nostop_byid)
            // -> tuple(bait, [samples], [contig_files], bait_stop_fasta, bait_nostop_fasta)

        REFORMER_SRA(ch_reformer_ready)

        ch_exonerate_input = REFORMER_SRA.out.exonerate_input
        ch_bait_ref = REFORMER_SRA.out.bait_ref

    } else {
        log.info "Skipping bait reformatting"

        def reform_dir = file(params.reformer_sra_dir)

        if (!reform_dir.isDirectory()) {
            error """
            Missing reformatted bait directory:
            ${reform_dir}

            Either:
            - run with --skip_reformat_sra false
            - provide previous REFORMER_SRA output
            """
        }

        ch_exonerate_input = Channel
            .fromPath("${reform_dir}/*", type: 'dir')
            .map { contigs_dir ->
                def bait = contigs_dir.name
                tuple(bait, contigs_dir)
            }


        ch_bait_ref = Channel.fromPath("${reform_dir}/*", type: 'dir')
            .filter { it.isDirectory() && !(it.name in [
                    'with.stop',
                    'without.stop',
                    '5.2-exonerate.output',
                    '5.3-stitcher.output'
                ])
            }
            .map { contigs_dir ->
                def bait = contigs_dir.name
                tuple(
                    bait,
                    file("${params.reformer_sra_dir}/with.stop/${bait}.fasta"),
                    file("${params.reformer_sra_dir}/without.stop/${bait}.fasta")
                )
            }
    }

    // ---------- STEP 5b - exonerate stitching ----------

    if (!params.skip_exonerate_sra) {

        ch_exonerate_ready = ch_exonerate_input
            .join(ch_bait_ref)
            // (bait, contigs_dir, stop_fasta, nostop_fasta) -> drop stop_fasta, EXONERATE_SRA doesn't use it
            .map { bait, contigs_dir, stop_fasta, nostop_fasta ->
                tuple(bait, contigs_dir, nostop_fasta)
            }

        EXONERATE_SRA(ch_exonerate_ready)

        ch_exonerate_results = EXONERATE_SRA.out.results

    } else if (!params.skip_stitcher_sra) {
        def exonerate_dir = file(params.exonerate_sra_dir)

        if (!exonerate_dir.isDirectory()) {
            error """
            Missing EXONERATE output directory:

            ${exonerate_dir}

            Either:
            - run with --skip_exonerate_sra false
            - provide existing EXONERATE_SRA output
            """
        }

        // reconstruct from previous EXONERATE_SRA publishDir
        ch_exonerate_results = Channel
            .fromPath("${exonerate_dir}/*/exonerate_out", type: 'dir')
            .map { exonerate_out ->
                def bait = exonerate_out.getParent().name
                tuple(bait, exonerate_out)
            }

        log.info "Skipping exonerate stitching"

    }

    // ---------- STEP 5c - stitching (ExAMiner) ----------

    if (!params.skip_stitcher_sra) {

        ch_stitcher_ready = ch_exonerate_results
            .join(ch_bait_ref)
            .map { bait, exonerate_dir, stop_fasta, nostop_fasta ->
                tuple(bait, exonerate_dir, nostop_fasta)
            }

        STITCHER_SRA(
            ch_stitcher_ready,
            Channel.value(file(params.samples_list)),
            Channel.value(file(params.outgroups_list)),
        )

        ch_tblout = STITCHER_SRA.out.tblout

        // manifest.tsv rows are "<sample>\t<mapping_out>/<sample>/<bait>.fasta",
        // written relative to the STITCHER_SRA task workdir the manifest itself
        // lives in - resolve against manifest_file's parent to get real paths.
        ch_sample_bait_fastas = STITCHER_SRA.out.manifest
            .flatMap { bait, manifest_file ->
                manifest_file.splitCsv(sep: '\t').collect { row ->
                    def sample = row[0]
                    def fasta_file = manifest_file.getParent().resolve(row[1])
                    tuple(sample, file(fasta_file))
                }
            }
            .groupTuple(by: 0)   // -> tuple(sample, [bait_fastas])

    } else if (
        params.results_from_third_party == null ||
        params.results_from_third_party == "" ||
        params.results_from_third_party == "stitcher"
    ) {

        // Assume stitching already ran; reconstruct from published output
        def stitched_dir = file("${params.stitcher_out_dir}")

        if (!stitched_dir.isDirectory()) {
            error """
            Missing stitcher output directory: ${stitched_dir}

            Either:
            - run with --skip_stitcher_sra false to regenerate it
            - provide existing STITCHER_SRA output in:
              ${params.stitcher_out_dir}
            """
        }

        ch_tblout = Channel.fromPath("${stitched_dir}/*/mapping_out/tblout/*.tsv")

        ch_sample_bait_fastas = Channel
            .fromPath("${stitched_dir}/*/mapping_out/*/*.fasta")
            .map { f -> tuple(f.getParent().getName(), f) }
            .groupTuple(by: 0)

        log.info "Skipping ExAMiner stitching"

    } else {
        // Stitcher is skipped and a third-party source will be used.
        log.info "Skipping ExAMiner stitching; using third-party source: ${params.results_from_third_party}"
    }

    // ---------- STEP 6a - exon/intron coordinates ----------

    if (!params.skip_consensus) {

        PREPARE_COORDINATES(ch_tblout.collect())
        ch_coordinates = PREPARE_COORDINATES.out.coordinates

    }

    // ---------- STEP 6b - mapping ----------

    if (!params.skip_mapping) {

        ch_mapper_input = ch_sample_bait_fastas
            .join(ch_reads_compressed)
            .map { sample, bait_fastas, r1, r2 ->
                if (!insert_sizes.containsKey(sample)) {
                    error "Missing insert size for sample: ${sample}"
                }
                tuple(sample, bait_fastas, r1, r2, insert_sizes[sample])
            }

        MAPPER(ch_mapper_input)

        ch_mapper_bam       = MAPPER.out.bam
        ch_mapper_reference = MAPPER.out.reference

    } else if (!params.skip_consensus) {

        // Assume mapping already ran; reconstruct from published output
        ch_mapper_bam = ch_samples.map { s ->
            def bam = file("${params.mapping_out_dir}/6.1-mapping/${s}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam")
            def bai = file("${bam}.bai")

            if (!bam.exists() || !bai.exists()) {
                error """
                Missing mapped BAM for sample: ${s}

                Expected:
                ${bam}
                ${bai}

                Either:
                - run with --skip_mapping false to regenerate it
                - provide existing MAPPER output in:
                  ${params.mapping_out_dir}/6.1-mapping
                """
            }

            tuple(s, bam, bai)
        }

        ch_mapper_reference = ch_samples.map { s ->
            def loci_fasta = file("${params.mapping_out_dir}/refs/${s}_loci.fasta")
            def loci_idx = file("${params.mapping_out_dir}/refs/${s}_loci.fasta.fai")

            if (!loci_fasta.exists() || !loci_idx.exists()) {
                error """
                Missing reference/index for sample: ${s}

                Expected:
                ${loci_fasta}
                ${params.mapping_out_dir}/refs/${s}_loci.fasta.*

                Either:
                - run with --skip_mapping false to regenerate it
                - provide existing MAPPER output in:
                  ${params.mapping_out_dir}/refs
                """
            }

            tuple(s, loci_fasta, loci_idx)
        }

        log.info "Skipping mapping"
    }

    // ---------- STEP 6c - consensus ----------

    if (!params.skip_consensus) {

        ch_consensus_input = ch_mapper_bam
            .join(ch_mapper_reference)
            // -> tuple(sample, bam, bai, loci_fasta, loci_fai_files)

        CONSENSUS(
            ch_consensus_input,
            Channel.value(file(params.baits_list)),
            ch_coordinates
        )

        ch_consensus_exons = CONSENSUS.out.exons

        } else if (!params.skip_alignment && !params.skip_mapping) {

        log.info "Skipping consensus calling"

        def consensus_dir = file(params.consensus_dir)

        if (!consensus_dir.isDirectory()) {
            error """
            Missing consensus output directory:

            ${consensus_dir}

            Either:
            - run with --skip_consensus false
            - provide existing CONSENSUS output
            """
        }
 
        ch_consensus_exons = Channel
            .fromPath("${params.consensus_dir}/${params.ucode}/2-stitched.exons/*.fasta")
            .map { f ->
                def sample = f.baseName.split('_')[-1]
                tuple(sample, f)
            }
            .groupTuple(by:0)

    }

    // ---------- STEP 7 - alignment (MACSE) ----------

    if (!params.skip_alignment) {

        if (params.skip_mapping && params.skip_consensus) {

        // Mapping and consensus are skipped.
        // Prepare alignment input directly from:
        // STITCHER (default), aTRAM, or Patchwork.

        log.info "Skipping mapping and consensus"
        def source = params.results_from_third_party ?: 'stitcher'

        if (source == 'stitcher' && !params.skip_stitcher_sra) {
            // STITCHER_SRA runs in this same execution — consume its
            // output channel directly, don't gate on disk.
            log.info "Preparing alignment input from STITCHER_SRA (this run)"

            ch_aligner_input = STITCHER_SRA.out.stitched
                .map { bait, stitched_dir ->
                    tuple(bait, [file("${stitched_dir}/${bait}.fasta")])
                }

        } else {
            log.info "Preparing alignment input from: ${source}"
            def alignment_source = source == 'atram'     ? file(params.atram_out_dir)
                                  : source == 'patchwork' ? file(params.patchwork_out_dir)
                                  :                          file(params.stitcher_out_dir)

            if (!alignment_source.isDirectory()) {
                error """
                Missing alignment input directory:

                ${alignment_source}

                Selected source:
                    ${params.results_from_third_party ?: 'stitcher'}
                """
            }

            /*
            * PREPARE_ALIGNMENT_INPUT converts:
            *
            * prepared/
            *   RZB38597.1.fasta
            *   RZB38599.1.fasta
            *   RZB38609.1.fasta
            *   ...
            *
            * into:
            * tuple(bait, [fasta]) = ALIGNER input
            */

            PREPARE_ALIGNMENT_INPUT(alignment_source)
            ch_aligner_input = PREPARE_ALIGNMENT_INPUT.out.alignment_input
                .flatten()
                .map { fasta -> tuple(fasta.baseName, [fasta]) }
        }

        ALIGNER(
            ch_aligner_input,
            Channel.value(file(params.samples_list))
            )
            
        ch_alignments = ALIGNER.out.alignment

        } else {
        // Normal pipeline:
        // STITCHER → MAPPER → CONSENSUS → ALIGNER

            ch_aligner_input = ch_consensus_exons
                .flatMap { sample, fasta_files ->
                    fasta_files.collect { fasta ->
                        def bait = fasta.baseName.replace("_${sample}", "")
                        tuple(bait, sample, fasta)
                    }
                }
                .groupTuple(by: 0)
                .map { bait, samples, exon_fastas ->
                    tuple(bait, exon_fastas)
                }

            ALIGNER(
                ch_aligner_input,
                Channel.value(file(params.samples_list))
            )

            ch_alignments = ALIGNER.out.alignment
        }

    } else if (!params.skip_cleaning) {

        // Assume alignment already ran; reconstruct from published output
        def align_out_dir = "${params.outdir}/7-aligning/${params.ucode}_N${params.n_percent}S${params.s_percent}/3-alignments"

        if (!file(align_out_dir).isDirectory()) {
            error """
            Missing alignment output directory: ${align_out_dir}

            Either:
            - run with --skip_alignment false to regenerate it
            - provide existing ALIGNER output in:
              ${align_out_dir}
            """
        }

        def nt_files = file("${align_out_dir}").listFiles().findAll { it.name.endsWith("_NT.fasta") }

        if (nt_files.isEmpty()) {
            error """
            Missing alignment FASTA files:

            ${align_out_dir}/*_NT.fasta

            Either:
            - run with --skip_alignment false
            - provide existing ALIGNER output
            """
        }

        ch_alignments = Channel
            .fromPath("${align_out_dir}/*_NT.fasta")
            .map { nt_fasta ->
                def bait = nt_fasta.baseName.replace('_NT', '')
                def aa_fasta = file("${align_out_dir}/${bait}_AA.fasta")

                if (!aa_fasta.exists()) {
                    error """
                    Missing AA fasta for bait '${bait}':
                    ${aa_fasta}

                    NT and AA outputs should always be published together by ALIGNER;
                    this suggests a partial/corrupted publishDir directory.
                    """
                }

                tuple(bait, nt_fasta, aa_fasta)
            }

        log.info "Skipping alignment"
    }

    // ---------- STEP 8a - cleaning (SEDA stop/frameshift filter + IQ-TREE gene trees) ----------

    if (!params.skip_cleaning) {

        ch_cleaning_input = ch_alignments
            .map { bait, nt_fasta, aa_fasta -> tuple(bait, nt_fasta) }

        CLEANING(
            ch_cleaning_input,
            Channel.value(file(params.outgroups_list))
        )

        ch_gene_trees = CLEANING.out.gene_tree
            .flatMap { bait, all_fasta, treefile, runtrees -> [all_fasta, treefile, runtrees] }
        ch_ingroup = CLEANING.out.ingroup
            .flatMap { bait, ingroup_fasta -> [ingroup_fasta] }

    } else if (!params.skip_filtering) {
        log.info "Skipping cleaning"

        // Assume cleaning already ran; reconstruct from published output
        def clean_dir = file("${params.outdir}/8-cleaning/${params.ucode}_N${params.n_percent}S${params.s_percent}_fSTP_fFRS/8.5-filtering.input")

        ch_gene_trees = Channel.fromPath([
            "${clean_dir}/*.all.fasta",
            "${clean_dir}/*.treefile",
            "${clean_dir}/*.runtrees"
        ])

        ch_ingroup = Channel.fromPath("${clean_dir}/*.ingroup.fasta")
    }

    // ---------- STEP 8b - filtering ----------
    if (!params.skip_filtering) {

        FILTERING(ch_gene_trees.mix(ch_ingroup).collect(), file(params.physsu_script))

    } else {

        log.info "Skipping filtering"

    }

    // ---------- STEP 9a - partition model selection (PartitionFinder) ----------

    if (!params.skip_partitionfinder) {

        if (!params.skip_filtering) {
            ch_pf_alignments = FILTERING.out.alignments_all
                .flatten()
                .map { fasta ->
                    tuple(fasta.baseName.replaceFirst(/\.all$/, ''), fasta)
                }

            ch_pf_stats = FILTERING.out.postfilter_stats

        } else {
            log.info "Reconstructing FILTERING output for PartitionFinder"

            def align_dir  = file("${params.final_out_dir}/alignments")
            def stats_file = file("${params.final_out_dir}/stats/postfilter_stats_all.txt")

            if (!align_dir.isDirectory() || !stats_file.isFile()) {
                error """
                Missing FILTERING output required for PartitionFinder:

                ${align_dir}
                ${stats_file}

                Either:
                - run with --skip_filtering false to regenerate it
                - provide existing FILTERING output
                """
            }

            ch_pf_alignments = Channel
                .fromPath("${align_dir}/*.all.fasta")
                .map { fasta ->
                    tuple(fasta.baseName.replaceFirst(/\.all$/, ''), fasta)
                }

            ch_pf_stats = Channel.value(stats_file)
        }

        ch_pf_ready = ch_pf_alignments
            .combine(ch_pf_stats)
            .combine(Channel.value(file(params.template_pf_cfg)))
            .map { bait, fasta, stats, template -> tuple(bait, fasta, stats, template) }

        PARTITIONFINDER(ch_pf_ready)

        ch_best_scheme = PARTITIONFINDER.out.best_scheme

    } else if (!params.skip_mrbayes) {

        log.info "Skipping PartitionFinder"

        def partfi_dir = file(params.partfi_dir)

        if (!partfi_dir.isDirectory()) {
            error """
            Missing PartitionFinder output directory:

            ${partfi_dir}

            Either:
            - run with --skip_partitionfinder false to regenerate it
            - provide existing PARTITIONFINDER output
            """
        }

        ch_best_scheme = Channel
            .fromPath("${partfi_dir}/*/analysis/best_scheme.txt")
            .map { f -> tuple(f.getParent().getParent().name, f) }
    }

    // ---------- STEP 9b - Bayesian gene trees (MrBayes) ----------

    if (!params.skip_mrbayes) {

        ch_mb_alignments = (!params.skip_filtering)
            ? FILTERING.out.alignments_all.flatten()   // one file per item
            : Channel.fromPath("${params.final_out_dir}/alignments/*.all.fasta")

        ch_mb_alignments = ch_mb_alignments
            .map { fasta ->
                tuple(fasta.baseName.replaceFirst(/\.all$/, ''), fasta)
            }

        ch_mb_ready = ch_mb_alignments.join(ch_best_scheme, by: 0)

        MRBAYES(
            ch_mb_ready,
            Channel.value(file(params.outgroups_list))
        )

    } else {
        log.info "Skipping MrBayes"
    }
    
    // ---------- STEP 9c - final datasets (crossval whitelist + BUSCO subset) ----------

    if (!params.skip_datasets) {

        // BI trees: staged from this run's MRBAYES tasks when they ran (this also makes
        // the step wait for every MrBayes task), otherwise reconstructed from disk
        if (!params.skip_mrbayes) {
            ch_bi_trees = MRBAYES.out.bi_tree.collect()
        } else {
            log.info "Reconstructing MRBAYES output for dataset generation"
            if (!file("${params.bi_trees_dir}/*.BI.treefile")) {
                error "No *.BI.treefile found in ${params.bi_trees_dir}; run with --skip_mrbayes false"
            }
            ch_bi_trees = Channel.fromPath("${params.bi_trees_dir}/*.BI.treefile").collect()
        }

        // FILTERING outputs must already be published if FILTERING is skipped
        if (params.skip_filtering) {
            checkDir('crossval_ml_trees_dir', params.crossval_ml_trees_dir)
            checkDir('crossval_alignments_dir', params.crossval_alignments_dir)
            checkFile('crossval_stats', params.crossval_stats)
        }

        GENERATE_DATASETS(
            ch_bi_trees,
            Channel.value(file(params.crossval_ml_trees_dir)),
            Channel.value(file(params.crossval_alignments_dir)),
            Channel.value(file(params.crossval_stats)),
            Channel.value(file(params.baits_list, checkIfExists: true)),
            Channel.value(file(params.outgroups_list, checkIfExists: true)),
            Channel.value(file(params.busco_table, checkIfExists: true)),
            Channel.value(file(params.datagene_script, checkIfExists: true))
        )

    } else {
        log.info "Skipping dataset generation"
    }
}
