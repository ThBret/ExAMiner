process PREPARE_COORDINATES {
    tag "coordinates_${params.ref_tax}"
    label "local_small"
    publishDir "${params.infodir}", mode: 'copy'

    input:
        path tblout_files

    output:
        path "list.coordinates.${params.ref_tax}.txt", emit: coordinates

    script:
    """
    # Combine stitcher tblout results into exon-intron coordinate file
    cat ${tblout_files} > list.coordinates.${params.ref_tax}.txt
    """
}
