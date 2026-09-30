from Bio.SeqRecord import SeqRecord
from Bio import SearchIO, SeqIO
from Bio.Seq import Seq
import numpy as np
import argparse
import logging
import utils
import csv
import sys
import os

CODON_SIZE = utils.CODON_SIZE

Segment = utils.Segment
Intron = utils.Intron
AlignmentReference = utils.AlignmentReference

# ============================================================================
# Discard/skip tracker
# ============================================================================
# Records why fragments, HSPs, or samples are skipped/discarded.
# `events` contains all individual events, while `sample_status` stores the
# final status of each sample. Overlap events are kept separately because
# they need more detailed sequence/coordinate information.
class DiscardTracker:
    def __init__(self, bait):
        self.bait = bait
        self.events = []          # every discard/skip/warning event, in order
        self.sample_status = {}   # sample -> 'kept' | 'discarded: <reason>'
        self.overlap_events = []  # dedicated overlap log

    def log(self, sample, stage, reason, detail=''):
        """Record an event without necessarily discarding the sample.
        Some events only affect one fragment/HSP or a sample's eligibility as a
        reference. The sample may still be successfully stitched later."""
        entry = {'bait': self.bait, 'sample': sample, 'stage': stage,
                 'reason': reason, 'detail': detail}
        self.events.append(entry)

    def discard(self, sample, stage, reason, detail=''):
        """Record an event and mark the sample as discarded (first reason
        wins if a sample is flagged more than once)."""
        self.log(sample, stage, reason, detail)
        if self.sample_status.get(sample, '').startswith('discarded'):
            return  # keep the first/root-cause reason
        self.sample_status[sample] = f'discarded: {stage} / {reason} ({detail})' if detail else f'discarded: {stage} / {reason}'

    def keep(self, sample, note=''):
        """Mark a sample as kept unless it was already discarded."""
        if not self.sample_status.get(sample, '').startswith('discarded'):
            self.sample_status[sample] = f'kept ({note})' if note else 'kept'

    def write_report(self, events_path, summary_path):
        os.makedirs(os.path.dirname(events_path), exist_ok=True)
        with open(events_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['bait', 'sample', 'stage', 'reason', 'detail'])
            writer.writeheader()
            for e in self.events:
                writer.writerow(e)

        os.makedirs(os.path.dirname(summary_path), exist_ok=True)
        with open(summary_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['bait', 'sample', 'status'])
            for sample, status in sorted(self.sample_status.items()):
                writer.writerow([self.bait, sample, status])

    def log_summary(self):
        kept = [s for s, st in self.sample_status.items() if st.startswith('kept')]
        discarded = [s for s, st in self.sample_status.items() if st.startswith('discarded')]
        logging.info(f'=== SUMMARY for bait {self.bait}: '
                     f'{len(kept)} kept, {len(discarded)} discarded '
                     f'(of {len(self.sample_status)} evaluated) ===')
        for sample in sorted(discarded):
            logging.info(f'  DISCARDED  {sample}: {self.sample_status[sample]}')

    def log_overlap(self, sample, segment1_idx, segment1_coord, segment2_idx, segment2_coord, overlap_aa,
                     resolution, overlap_seq_1_aa='', overlap_seq_2_aa='',
                     segment_1_seq_aa='', segment_2_seq_aa='', detail=''):
        """Record every segment-overlap check from STEP C-10, independent of the
        general discard/keep events, so overlaps can be assessed on their own
        regardless of whether they were trimmed or proved fatal to the sample.
        Sequence fields contain the overlap itself and the resulting full
        segment sequences after any trimming/masking."""
        entry = {
            'bait': self.bait, 'sample': sample,
            'segment_i': segment1_idx, 'segment_i_start': segment1_coord[0], 'segment_i_end': segment1_coord[1],
            'segment_j': segment2_idx, 'segment_j_start': segment2_coord[0], 'segment_j_end': segment2_coord[1],
            'overlap_aa': overlap_aa, 'resolution': resolution,
            'overlap_seq_i_aa': overlap_seq_1_aa, 'overlap_seq_j_aa': overlap_seq_2_aa,
            'segment_i_seq_aa': segment_1_seq_aa, 'segment_j_seq_aa': segment_2_seq_aa,
            'detail': detail,
        }
        self.overlap_events.append(entry)

    def write_overlap_report(self, path):
        """Write a human-readable report of all segment-overlap events."""
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            if not self.overlap_events:
                f.write(f'No segment overlaps recorded for bait {self.bait}.\n')
                return

            f.write(f'Overlap report for bait: {self.bait}\n')
            f.write(f'Total overlap events: {len(self.overlap_events)}\n')
            f.write('=' * 80 + '\n\n')

            for n, e in enumerate(self.overlap_events, start=1):
                f.write(f'[{n}] Sample: {e["sample"]}\n')
                f.write(f'    Segments: [{e["segment_i"]}] '
                        f'({e["segment_i_start"]}-{e["segment_i_end"]}) aa  <->  '
                        f'[{e["segment_j"]}] ({e["segment_j_start"]}-{e["segment_j_end"]}) aa\n')
                f.write(f'    Overlap length: {e["overlap_aa"]} aa\n')
                f.write(f'    Overlapping region, segment [{e["segment_i"]}] side: {e["overlap_seq_i_aa"]}\n')
                f.write(f'    Overlapping region, segment [{e["segment_j"]}] side: {e["overlap_seq_j_aa"]}\n')
                f.write(f'    Strategy applied: {e["resolution"]}\n')

                if e["detail"]:
                    f.write(f'    Detail: {e["detail"]}\n')

                if "dropped_fully_redundant" in e["resolution"]:
                    f.write(f'    Kept segment [{e["segment_i"]}] sequence:\n')
                    f.write(f'        {e["segment_i_seq_aa"]}\n')
                    f.write(f'    Removed segment [{e["segment_j"]}] sequence:\n')
                    f.write(f'        {e["segment_j_seq_aa"]}\n')
                elif "_masked" in e["resolution"]:
                    # Masking may fully consume either side (or both); an empty
                    # translated sequence for a side means it was dropped entirely.
                    if e["segment_i_seq_aa"]:
                        f.write(f'    Remaining segment [{e["segment_i"]}] sequence:\n')
                        f.write(f'        {e["segment_i_seq_aa"]}\n')
                    else:
                        f.write(f'    Masked (removed) segment [{e["segment_i"]}]: fully consumed by mask\n')

                    if e["segment_j_seq_aa"]:
                        f.write(f'    Remaining segment [{e["segment_j"]}] sequence:\n')
                        f.write(f'        {e["segment_j_seq_aa"]}\n')
                    else:
                        f.write(f'    Masked (removed) segment [{e["segment_j"]}]: fully consumed by mask\n')
                else:
                    f.write(f'    Resulting sequence, segment [{e["segment_i"]}]:\n')
                    f.write(f'        {e["segment_i_seq_aa"]}\n')
                    f.write(f'    Resulting sequence, segment [{e["segment_j"]}]:\n')
                    f.write(f'        {e["segment_j_seq_aa"]}\n')

                f.write('-' * 80 + '\n')

def _segment_ref_score(fragment, crop_coord, gap_pos):
    """Similarity score of a fragment's cropped region against the bait, for
    the same coordinate window used to build the Segment in STEP C-8.
    Mirrors the scoring already done in STEP C-4's HSP filtering, just
    applied per-segment instead of per-HSP."""
    ref_st = max(crop_coord[0], fragment.query_start)
    st = ref_st - fragment.query_start
    ref_end = min(crop_coord[1], fragment.query_end)
    end = ref_end - fragment.query_start
    st += sum([l for p, l in gap_pos if p <= st])
    end += sum([l for p, l in gap_pos if p < end])
    return utils.sim_score(fragment.query[st:end], fragment.hit[st:end])

def _resolve_overlap(segments, keep_idx, mod_idx, overlap_len, ov_seq_keep_aa, ov_seq_mod_aa,
                      sample, tracker, reason_tag, extra_detail=''):
    """Resolve an overlap by keeping one segment and trimming the other.
    If the overlap consumes the entire modified segment, that segment is
    removed. Otherwise only its overlapping end is trimmed.
    Returns: Index of a segment that should be removed, or None."""
    keep = segments[keep_idx]
    mod = segments[mod_idx]

    if overlap_len >= mod.length:
        tracker.log_overlap(sample, keep_idx, (keep.start, keep.end), mod_idx, (mod.start, mod.end),
                             overlap_len, resolution=f'{reason_tag}_dropped_fully_redundant',
                             overlap_seq_1_aa=ov_seq_keep_aa, overlap_seq_2_aa=ov_seq_mod_aa,
                             segment_1_seq_aa=str(Seq(keep.seq).translate()),
                             segment_2_seq_aa=str(Seq(mod.seq).translate()),
                             detail=f'segment[{mod_idx}] fully redundant with segment[{keep_idx}]; '
                                    f'dropped, no novel sequence. {extra_detail}')
        return mod_idx

    # Trim from whichever side actually overlaps `keep`
    if mod.start >= keep.start:
        mod.seq = mod.seq[overlap_len * CODON_SIZE:]
        mod.start += overlap_len
    else:
        mod.seq = mod.seq[:-(overlap_len * CODON_SIZE)]
        mod.end -= overlap_len
    mod.length = mod.end - mod.start

    tracker.log_overlap(sample, keep_idx, (keep.start, keep.end), mod_idx, (mod.start, mod.end),
                         overlap_len, resolution=f'{reason_tag}_trimmed',
                         overlap_seq_1_aa=ov_seq_keep_aa, overlap_seq_2_aa=ov_seq_mod_aa,
                         segment_1_seq_aa=str(Seq(keep.seq).translate()),
                         segment_2_seq_aa=str(Seq(mod.seq).translate()),
                         detail=f'trimmed {overlap_len}aa off segment[{mod_idx}]; '
                                f'new coords=({mod.start}, {mod.end}). {extra_detail}')
    return None

def _mask_overlap(segments, idx1, idx2, overlap_len, overlap_seq_1_aa, overlap_seq_2_aa,
                   sample, tracker, reason_tag, extra_detail=''):
    """Resolve an overlap by removing it from both segments.
    The resulting gap is filled with `unk_ch` during stitching. If the overlap
    consumes an entire segment, that segment is removed."""
    first = segments[idx1]
    second = segments[idx2]
    to_drop = set()

    if overlap_len >= first.length:
        to_drop.add(idx1)
    else:
        first.seq = first.seq[:-(overlap_len * CODON_SIZE)]
        first.end -= overlap_len
        first.length = first.end - first.start

    if overlap_len >= second.length:
        to_drop.add(idx2)
    else:
        second.seq = second.seq[overlap_len * CODON_SIZE:]
        second.start += overlap_len
        second.length = second.end - second.start

    tracker.log_overlap(
        sample, idx1, (first.start, first.end), idx2, (second.start, second.end),
        overlap_len, resolution=f'{reason_tag}_masked',
        overlap_seq_1_aa=overlap_seq_1_aa, overlap_seq_2_aa=overlap_seq_2_aa,
        segment_1_seq_aa='' if idx1 in to_drop else str(Seq(first.seq).translate()),
        segment_2_seq_aa='' if idx2 in to_drop else str(Seq(second.seq).translate()),
        detail=f'overlap masked with missing data instead of discarding the sample; '
               f'{overlap_len}aa removed from both segments (one dropped entirely if fully '
               f'consumed). {extra_detail}')
    return to_drop


def fragment_to_segment(fragment, crop_coord, gap_pos, sample, tracker):
    """Crop a fragment to a reference exon and return a Segment.
    Returns: Segment, or None if the fragment cannot be mapped to valid
    reference coordinates."""
    frag_seq, start, end = utils.crop_fragment_seq(fragment, crop_coord)
    detail = f'query_range={fragment.query_range}; exon={crop_coord}'

    if not frag_seq:
        logging.warning(f'Skipping fragment: no sequence after cropping ({detail})')
        tracker.log(sample, 'segment_creation', 'empty_cropped_sequence', detail)
        return None

    if start is None or end is None:
        logging.warning(f'Skipping fragment: undefined coordinates ({detail})')
        tracker.log(sample, 'segment_creation', 'undefined_coordinates', detail)
        return None

    segment = Segment(frag_seq, start, end)
    segment.score = _segment_ref_score(fragment, crop_coord, gap_pos)

    return segment


# ============================================================================
# Parser
# ============================================================================

parser = argparse.ArgumentParser()
parser.add_argument('bait', help='Bait fasta file path')
parser.add_argument('-t', '--threshold',
                    help='Percentage threshold to consider a partial reference',
                    type=float, default=0.8)
parser.add_argument('-b', '--nbuffer',
                    help='Length of the buffer of "N"s to add to both ends '
                    'in the output for mapping',
                    type=int,
                    default=1000)
parser.add_argument('-c', '--complete-only',
                    help='Do not consider partial references and '
                    'use complete references exclusively',
                    action='store_true')
parser.add_argument('-g', '--outgroups',
                    help='A file with a list of IDs of outgroup samples (one ID per line)',
                    required=False)
parser.add_argument('-s', '--samples',
                    help='A file with a list of IDs of samples, including possible outgroups '
                    '(one ID per line)',
                    required=True)
parser.add_argument('-e', '--exoneratepath',
                    help='Exonerate aligment path',
                    required=True)
parser.add_argument('-m', '--mappingout',
                    help='Directory for mapping output',
                    required=True)
parser.add_argument('-o', '--output',
                    help='Output directory',
                    required=True)
parser.add_argument('-u', '--unknown',
                    help='Character used to fill unknown or missing nucleotides '
                    'in the output alignments ("N" or "?", default: "N")',
                    choices=['?', 'N'],
                    default='N')
parser.add_argument('-S', '--scorediff',
                    help='Score difference to discard an alignment. An alignment '
                    'is discarded if it has a score that is less than the '
                    'reference score minus the value of this option. '
                    'Alignment scores range from 0 to 1 (default: 0.10)',
                    type=float,
                    default=0.10)
parser.add_argument('-M', '--overlap-mismatch-tolerance',
                    help='Maximum number of amino-acid mismatches within an overlapping '
                    'region between two segments that will still be resolved by trimming '
                    '(keeping the higher-scoring segment) rather than discarding the whole '
                    'sample (default: 0)',
                    type=int, default=0)
parser.add_argument('-C', '--overlap-conflict-resolution',
                    help='How to resolve a segment overlap that cannot be confidently '
                    'trimmed (conflicting AA sequences over the mismatch tolerance or '
                    'tied/unavailable segment scores). "mask" removes the overlapping '
                    'amino acids from both segments and lets the gap be filled with '
                    'missing-data characters, keeping the rest of the sample (default)'
                    'while "discard" drops the whole sample (more conservative).',
                    choices=['mask','discard'], default='mask')
# Optional custom path for the discard report; defaults to
# living alongside the existing tblout directory so nothing else changes.
parser.add_argument('-r', '--discardreport',
                    help='Directory to write the per-bait discard report TSVs '
                    '(default: <mappingout>/tblout)',
                    required=False, default=None)
args = parser.parse_args()

logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', level=logging.INFO)


# ============================================================================
# STEP A: IDENTIFICATION OF COMPLETE AND PARTIAL REFERENCES
# ============================================================================

# Minimum percentage of the bait that should be covered for an
# alignment to be considered as a partial reference
partref_threshold = args.threshold

# Remove ".fasta" extension from bait file name
bait_file = args.bait
bait = os.path.splitext(os.path.basename(bait_file))[0]

# Instantiate the tracker now that `bait` is known
tracker = DiscardTracker(bait)

# Exonerate output directory
exonerate_dir = args.exoneratepath

# Read sample IDs from sample file
samples = utils.list_from_file_lines(args.samples)
if args.outgroups:
    # Read outgroup IDs from outgroup file
    outgroups = utils.list_from_file_lines(args.outgroups)
else:
    outgroups = []

# Get sequence and length of the bait
bait_fasta = list(SeqIO.parse(bait_file, "fasta"))
if len(bait_fasta) != 1:
    raise ValueError(
        f"Expected exactly one bait sequence, found {len(bait_fasta)}"
    )

bait_seq = str(bait_fasta[0].seq)
bait_len = len(bait_fasta[0].seq)

# Dictionaries storing complete/partial bait alignment
# Key: sample ID, Value: AlignmentReference object
complete_refs = {}
partial_refs = {}

# Get reference coordinates from the non-outgroup samples
logging.info('Trying to get reference coordinates')
for sample in samples:
    if sample in outgroups:
        continue

    # Exonerate output files must follow a strict directory structure and naming convention:
    # <exonerate_dir>/<bait>/<sample>_exonerate.out
    fname = os.path.join(exonerate_dir, bait, f'{sample}_exonerate.out')

    # Parse exonerate output file
    try:
        query_results = list(SearchIO.parse(fname, 'exonerate-text'))
    except (ValueError, FileNotFoundError):
        logging.warning(f'Could not read exonerate output file ({bait}: {sample})')
        # This only affects candidacy as a *reference*, not
        # whether the sample gets stitched later on, so we log it as an
        # informational event rather than a final discard.
        tracker.log(sample, 'reference_search', 'no_exonerate_output_file', fname)
        continue

    if not query_results:
        logging.info(f'Empty exonerate output file ({bait}: {sample})')
        tracker.log(sample, 'reference_search', 'empty_exonerate_output_file', fname)
        continue

    # Assumes exactly one query (the bait) per file
    assert len(query_results) == 1
    exonerate_res = query_results[0]

    # Iterate over all hits for this query
    for hit in exonerate_res:
        # Iterate over all high-scoring alignments (HSPs) for this hit
        for hsp in hit:
            hsp_len = hsp.query_span

            # CASE 1: Full-length alignment (HSP covers the entire bait sequence)
            if hsp_len == bait_len:
                candidate_ref = AlignmentReference(hsp, bait_seq)
                if sample in complete_refs:
                    # Check if a complete reference is already present for this sample
                    logging.warning(f'A potential reference was already found '
                                    f'for sample {sample} in bait {bait}')
                    
                    # Keep the highest-scoring alignment
                    candidate_ref_score = candidate_ref.ref_aln_score()
                    current_ref_score = complete_refs[sample].ref_aln_score()
                    if candidate_ref_score > current_ref_score:
                        complete_refs[sample] = candidate_ref
                        logging.warning(f'Found another possible reference in the same exonerate output file with '
                                        f' a greater alignment score, substituting previous reference candidate '
                                        f'(sample {sample} in bait {bait})')
                else:
                    # First recorded full-length reference for this sample
                    complete_refs[sample] = candidate_ref

            # CASE 2: Biologically improbable alignment (greater length than bait sequence)
            elif hsp_len > bait_len:
                logging.warning('Somehow the alignment covers a length greater than the '
                                'length of the bait. This probably should not happen.')
                tracker.log(sample, 'reference_search', 'hsp_longer_than_bait',
                            f'hsp_len={hsp_len} bait_len={bait_len}')

            # Short HSPs below the partial-reference threshold are intentionally ignored.
            # A fragment should correspond to an exon or part of an exon, but should not cover
            # more than the length of an  exon assuming split sites were properly detected

            # CASE 3: Partial but high-coverage alignment (covering >= threshold% of bait but not its full length)
            elif hsp_len >= partref_threshold * bait_len:
                candidate_ref = AlignmentReference(hsp, bait_seq)
                if sample in partial_refs:
                    # Check if a partial reference is already present for this sample
                    logging.warning(f'A potential partial reference was already found '
                                    f'for sample {sample} in bait {bait}')
                    
                    # Keep the longest partial alignment
                    if hsp_len > partial_refs[sample].query_span:
                        partial_refs[sample] = candidate_ref
                else:
                    # First recorded partial reference for this sample
                    partial_refs[sample] = candidate_ref


# ============================================================================
# STEP B: REFERENCE SELECTION AND EXTRACTION OF EXON-INTRON COORDINATES
# ============================================================================

# Disregard partial references if '-c' or '--complete-only' is set to True
if args.complete_only:
    partial_refs = None

# Check for presence of at least a single reference
if not complete_refs and not partial_refs:
    logging.error(f'{bait}: No candidate samples with reference exons found.')
    sys.exit(1)

# Select potential references prioritising the complete ones
if complete_refs:
    potential_refs = complete_refs
elif partial_refs:
    logging.warning('Could not find any complete references for this bait.  '
                    'Using partial references.')
    potential_refs = partial_refs

# Get the consensus (most frequent) exon coordinates from potential references
# Example structure: ((1, 124), (125, 250), (250, 309), (309, 315), (315, 371), (372, 477))
ref_exon_coords = utils.get_ref_exon_coords(potential_refs)
# Filter references to only those matching the consensus exon structure
exonref = {k: r for k, r in potential_refs.items() if r.exon_coords == ref_exon_coords}

# Report reference candidates that were dropped for not matching the consensus exon structure (previously silent)
for k, r in potential_refs.items():
    if k not in exonref:
        tracker.log(k, 'reference_selection', 'exon_structure_mismatch_with_consensus',
                    f'sample_exon_coords={r.exon_coords} consensus={ref_exon_coords}')

# Calculate reference intron length as mean intron length from reference specimens sharing the consensus exon structure
intronref_lens = np.round(np.mean([r.intron_lens for _, r in exonref.items()], axis = 0)).astype(int)
# Store reference alignments for quality scoring
exonref_alns = {k: r.ref_aln for k, r in exonref.items()}

# Get split coordinates from one of the references 
# (should be the same in all of them since they have the same exon coordinates)
split_coords = next(iter(exonref.items()))[1].split_coords
split_pos_aa = [coord[0] for coord in split_coords]

# NOTE: Exon coordinates do not include the split aminoacid


# ============================================================================
# STEP C: FRAGMENT EVALUATION AND MAPPING TO THE EXON–INTRON MODEL
# ============================================================================

# Create mapping output directory
exintr_for_mapping_dir = args.mappingout
os.makedirs(os.path.join(exintr_for_mapping_dir, 'tblout'), exist_ok = True)

# Write log file to record exon/intron positions for each sample
exintr_tbl_file = os.path.join(exintr_for_mapping_dir, 'tblout', f'exon_intron_table_{bait}.tsv')

seq_records = []

# Run another iteration after establishing reference coordinates
for sample in samples:
    logging.info(f'Stitching sample {sample}')

    segments = []
    introns = []
    included_introns = []

    # If the sample already has a reference, just use that as its sequence.
    # This should avoid discarding it if there are more sequences that partially
    # align to the bait, producing an overlap that is never solved.
    if sample in exonref:
        segments = utils.segments_from_aln_ref(exonref[sample])
        introns = utils.introns_from_aln_ref(exonref[sample])
        included_introns = list(range(len(introns)))
        proper_sample = True
        tracker.keep(sample, note='used directly as consensus-structure reference')

    else:
        proper_sample = True

        # Read exonerate output for this sample
        fname = os.path.join(exonerate_dir, bait, f'{sample}_exonerate.out')

        try:
            query_results = list(SearchIO.parse(fname, 'exonerate-text'))
        except (ValueError, FileNotFoundError):
            logging.warning(f'Could not read exonerate output file ({bait}: {sample}).')
            tracker.discard(sample, 'stitching', 'no_exonerate_output_file', fname)
            continue

        if not query_results:
            tracker.discard(sample, 'stitching', 'empty_exonerate_output_file', fname)
            continue

        assert len(query_results) == 1

        exonerate_res = query_results[0]

        ### (1) Split HSPs into fragments
        # Get all high-scoring pairs (HSPs) from the alignment
        all_hsps = exonerate_res.hsps

        # Early filtering of bad alignments
        logging.info('Evaluating all fragments to discard poor alignments')
        hsp_idx_to_rm = []
        for hsp_idx, hsp in enumerate(all_hsps):
            for fragment_idx, fragment in enumerate(hsp):
                logging.info(f'{bait} {sample} {fragment.hit_id} {fragment.hit_description}, '
                             f'fragment {fragment_idx}')
                
                ### (2) Assign fragments to reference exons
                exon_overlaps = utils.exon_overlap(fragment.query_range, ref_exon_coords)

                # Stripping the 'X's of the split aminoacids is necessary because
                # they are not included in the coordinates `query_start' and `query_end'
                query_seq = str(fragment.query.seq).strip('X')
                gap_pos = utils.gap_pos_from_seq(query_seq)
                if sum([x > 0 for x in exon_overlaps]) > 1:
                    logging.warning('This fragment overlaps more than one reference exon')
                    tracker.log(sample, 'stitching', 'fragment_overlaps_multiple_exons_prefilter',
                                f'hsp_idx={hsp_idx} fragment_idx={fragment_idx} query_range={fragment.query_range} overlap={exon_overlaps}')

                ### (3) Crop fragments to exon boundaries
                # Get the exon with maximum overlap
                max_overlap_exon = exon_overlaps.index(max(exon_overlaps))
                crop_coord = ref_exon_coords[max_overlap_exon]
                ref_st = max(crop_coord[0], fragment.query_start)
                ref_end = min(crop_coord[1], fragment.query_end)

                ### (4) Compute similarity scores
                # Score fragment vs reference
                sc = _segment_ref_score(fragment, crop_coord, gap_pos)
                ref_aln_scores = []
                for _, aln in exonref_alns.items():
                    ref_aln_frg = utils.get_segment_from_aln(aln, (ref_st, ref_end))
                    ref_aln_scores.append(utils.sim_score(ref_aln_frg[0], ref_aln_frg[1]))

                # Reference alignment defines expected score
                # Fragment must not be worse than ref − scorediff
                ref_sc = np.mean(ref_aln_scores)
                threshold_sc = ref_sc - args.scorediff

                ### (5) Compare against reference mean − threshold
                # Discard bad HSPs (one bad fragment will lead to whole HSP being removed)
                if sc <= threshold_sc:
                    logging.info('Found fragment with low score')
                    logging.info(f'{fragment.hit_id}, fragment {fragment_idx} has score {sc}, which '
                                 f'is lower than {threshold_sc} (orig: {ref_sc}).')
                    hsp_idx_to_rm.append(hsp_idx)
                    # This removes the whole HSP, not just this fragment, so record which fragment triggered it
                    tracker.log(sample, 'stitching', 'hsp_removed_low_fragment_score',
                                f'hsp_idx={hsp_idx} fragment_idx={fragment_idx} query_range={fragment.query_range} '
                                f'score={sc:.4f} threshold={threshold_sc:.4f} ref_mean_score={ref_sc:.4f}')
                    break
                else:
                    logging.info(f'Fragment has good enough score ({sc} compared to {threshold_sc} '
                                 f'(orig: {ref_sc})).')

        # Keep only valid HSPs
        all_hsps = [h for i, h in enumerate(all_hsps) if i not in hsp_idx_to_rm]
        all_hsps = sorted(all_hsps, key = lambda x: x.query_start)

        # Avoid empty sequences full of 'N's
        if len(all_hsps) == 0:
            proper_sample = False
            n_removed = len(set(hsp_idx_to_rm))
            tracker.discard(sample, 'stitching', 'no_valid_hsps_remaining',
                            f'{n_removed} HSP(s) removed for low score; none passed filtering')
    
        # HSPs should be already ordered by start position
        for hsp_idx, hsp in enumerate(all_hsps):
            # These coordinates will be used to get the introns that will be included for this HSP 
            # after cropping it (excluding all of the introns outside the segment cropped after the exons)
            min_hsp_coord = np.inf
            max_hsp_coord = -1

            # Get introns in this HSP
            intron_pos = hsp.query_inter_ranges
            intron_lens = hsp.hit_inter_spans

            potential_introns = [Intron(l, p) for l, p in zip(intron_lens, intron_pos)]

            ### (6) Identify split codons
            # Get HSP gap positions
            hsp_query_seq = ''.join([str(f.query.seq.rstrip('X')) for f in hsp])
            hsp_gap_pos = utils.gap_pos_from_seq(hsp_query_seq)
            hsp_seq_nt = ''.join([''.join(ann['hit_annotation']) for ann in hsp.aln_annotation_all])

            for split_start, split_end in hsp.query_inter_ranges:
                if split_start in split_pos_aa:
                    logging.info('Found split site in the same position as the reference')
                    # Include the split site in the segments to stitch
                    n_gaps = sum([l for pos, l in hsp_gap_pos if pos <= (split_start - hsp.query_start)])
                    st = (split_start - hsp.query_start + n_gaps) * CODON_SIZE
                    split_seq = hsp_seq_nt[st:(st + (split_end - split_start) * CODON_SIZE)]
                    segments.append(Segment(split_seq, split_start, split_end))
                    segments[-1].score = None

            ### (7) Validate fragments against the reference exon structure
            # Check if fragments overlap a splice site
            proper_hsp = True
            for fragment_idx, fragment in enumerate(hsp):
                gap_pos = utils.gap_pos_from_seq(str(fragment.query.seq).strip('X'))
                overlaps_split_site = False
                for split in split_coords:
                    if utils.overlaps_point(fragment.query_range, split[0]):
                        logging.warning(f'This fragment overlaps a reference splice aminoacid in position {split}')
                        proper_hsp = False
                        overlaps_split_site = True
                        tracker.log(sample, 'stitching', 'fragment_overlaps_splice_site',
                                    f'hsp_idx={hsp_idx} fragment_idx={fragment_idx} query_range={fragment.query_range} splice_pos={split}')

                exon_overlaps = utils.exon_overlap(fragment.query_range, ref_exon_coords)

                # exon_overlaps example structure:
                #   [0, 35, 42, 0]
                #   exon 1:  0 aa
                #   exon 2: 35 aa
                #   exon 3: 42 aa
                #   exon 4:  0 aa

                ### (8) Handle fragments with significant overlap to multiple exons
                # Require at least 10 aa of overlap with each exon.
                significant_overlaps = [overlap >= 10 for overlap in exon_overlaps]

                # With example
                # --> [False, True, True, False]

                if sum(significant_overlaps) > 1:
                    # There might be a different exon-intron pattern from the reference,and an
                    # exon in this alignment might be overlapping multiple exons in the reference
                    logging.info(f'Fragment has significant overlap with multiple exons (sample={sample})')
                    for exon_idx in [i for i, x in enumerate(significant_overlaps) if x]:
                        crop_coord = ref_exon_coords[exon_idx]
                        segment = fragment_to_segment(fragment, crop_coord, gap_pos, sample, tracker)
                        if segment is not None:
                            segments.append(segment)
                            min_hsp_coord = min(min_hsp_coord, segment.start)
                            max_hsp_coord = max(max_hsp_coord, segment.end)
                else:
                    if sum(overlap > 0 for overlap in exon_overlaps) > 1 or overlaps_split_site:
                        logging.warning('This fragment overlaps more than one reference exon '
                                        'or a split aminoacid')
                        proper_hsp = False
                        tracker.log(sample, 'stitching', 'fragment_overlaps_multiple_exons_or_split',
                                    f'hsp_idx={hsp_idx} fragment_idx={fragment_idx} query_range={fragment.query_range} overlap={exon_overlaps}')

                    # If no significant overlap, assign the fragment to the exon it overlaps most
                    # Get the exon with maximum overlap
                    max_overlap_exon = exon_overlaps.index(max(exon_overlaps))
                    crop_coord = ref_exon_coords[max_overlap_exon]

                    # In partial references, some fragments might be outside the boundaries 
                    # of a reference exon without overlapping more than one reference exon
                    # or a split site, and they should be cropped as well.
                    segment = fragment_to_segment(fragment, crop_coord, gap_pos, sample, tracker)

                    if segment is not None:
                        segments.append(segment)
                        min_hsp_coord = min(min_hsp_coord, segment.start)
                        max_hsp_coord = max(max_hsp_coord, segment.end)

            ### (9) Intron reconstruction
            # This might happen if the HSP is outside the reference coordinates
            if min_hsp_coord is np.inf or max_hsp_coord == -1:
                logging.warning(f'HSP with coordinates {hsp.query_range} outside reference '
                                f'coordinates (partial reference)')
                # Note this is an intron-tracking gap for this HSP, not by itself a sample-level discard
                tracker.log(sample, 'stitching', 'hsp_outside_reference_coords_no_introns_assigned',
                            f'hsp_idx={hsp_idx} hsp_range={hsp.query_range}')
                min_hsp_coord, max_hsp_coord = 0, 0

            introns_in_hsp = [intr for intr in potential_introns
                           if intr.pos[0] >= min_hsp_coord and intr.pos[1] <= max_hsp_coord]
            # These are the introns that are already covered within the
            # coordinates of the added HSPs
            ref_introns_in_hsp = [i for i, pos in enumerate(split_coords)
                               if pos[0] > min_hsp_coord and pos[1] < max_hsp_coord]
            introns.extend(introns_in_hsp)
            included_introns.extend(ref_introns_in_hsp)

        ### (10) Resolve overlaps between cropped segments
        # Cropping each fragment independently can still leave segments sharing a few
        # amino acids. Compare all overlapping segment pairs and:
        #   - identical overlap       -> trim the lower-priority segment
        #   - minor mismatch          -> keep the higher-scoring segment
        #   - unresolved scores       -> discard or mask, depending on the CLI option
        #   - substantial mismatch    -> discard or mask, depending on the CLI option
        #
        # Segments are sorted first so that idx_first is normally the leftmost segment.
        logging.info('Checking for overlap after cropping exons')
        logging.info(f'No. segments: {len(segments)}')

        segments_to_drop = set()

        segments.sort(key=lambda s: s.start)   # sort before comparing
        for i in range(len(segments) - 1):
            for j in range(i + 1, len(segments)):
                if i in segments_to_drop or j in segments_to_drop:
                    continue

                segment1 = segments[i]
                segment2 = segments[j]
                segment1_coord = (segment1.start, segment1.end)
                segment2_coord = (segment2.start, segment2.end)
                overlap_len = utils.overlap(segment1_coord, segment2_coord)   # amino-acid units

                if overlap_len <= 0:
                    continue

                if segment1.start <= segment2.start:
                    idx_first, idx_second = i, j
                else:
                    idx_first, idx_second = j, i
                first, second = segments[idx_first], segments[idx_second]

                # Lengths can disagree between segments
                first_len_aa = len(first.seq) // CODON_SIZE
                second_len_aa = len(second.seq) // CODON_SIZE

                ov_seq_first_aa = str(Seq(first.seq[(first_len_aa - overlap_len) * CODON_SIZE:]).translate())
                ov_seq_second_aa = str(Seq(second.seq[:overlap_len * CODON_SIZE]).translate())

                n_mismatches = sum(a != b for a, b in zip(ov_seq_first_aa, ov_seq_second_aa))

                if n_mismatches == 0:
                    # Case 1: identical overlapping peptide -> trim (or drop if fully redundant)
                    dropped = _resolve_overlap(segments, idx_first, idx_second, overlap_len,
                                                ov_seq_first_aa, ov_seq_second_aa,
                                                sample, tracker, 'identical_overlap')
                    if dropped is not None:
                        segments_to_drop.add(dropped)

                elif n_mismatches <= args.overlap_mismatch_tolerance:
                    # Case 2: small mismatch -> keep the higher-scoring segment's copy
                    score_first = getattr(first, 'score', None)
                    score_second = getattr(second, 'score', None)
                    if score_first is not None and score_second is not None and score_first != score_second:
                        keep_idx, mod_idx = (idx_first, idx_second) if score_first >= score_second \
                                            else (idx_second, idx_first)

                        ov_seq_keep = ov_seq_first_aa if keep_idx == idx_first else ov_seq_second_aa
                        ov_seq_mod = ov_seq_second_aa if keep_idx == idx_first else ov_seq_first_aa
                        dropped = _resolve_overlap(segments, keep_idx, mod_idx, overlap_len, ov_seq_keep, ov_seq_mod,
                                                    sample, tracker, 'minor_mismatch',
                                                    extra_detail=f'{n_mismatches} aa mismatch(es); '
                                                                f'score_kept={getattr(segments[keep_idx], "score", None)} '
                                                                f'score_dropped={getattr(segments[mod_idx], "score", None)}')
                        if dropped is not None:
                            segments_to_drop.add(dropped)

                    else:
                        # Scores unavailable or tied - no principled way to pick which
                        # copy to keep, so treat this like an unresolved conflict rather
                        # than guessing (previously defaulted to keeping the earlier segment)
                        if args.overlap_conflict_resolution == 'mask':
                            dropped = _mask_overlap(
                                segments, idx_first, idx_second, overlap_len,
                                ov_seq_first_aa, ov_seq_second_aa,
                                sample, tracker, 'minor_mismatch_score_unavailable_or_tied',
                                extra_detail=f'{n_mismatches} aa mismatch(es); '
                                             f'score_first={score_first} score_second={score_second}')
                            segments_to_drop.update(dropped)
                            tracker.log(sample, 'stitching',
                                        'minor_mismatch_score_unavailable_or_tied_masked',
                                        f'segment[{idx_first}]=({first.start},{first.end}) overlaps '
                                        f'segment[{idx_second}]=({second.start},{second.end}) '
                                        f'overlap_aa={overlap_len} mismatches={n_mismatches}; masked, sample kept')
                        else:
                            logging.error(f'[{bait}] Overlapping segments in sample {sample} have a minor '
                                        f'mismatch ({n_mismatches} aa in {overlap_len}aa overlap) but scores are '
                                        f'unavailable or tied (score_first={score_first}, score_second={score_second}); '
                                        f'cannot decide which segment to keep')
                            proper_sample = False
                            tracker.discard(sample, 'stitching', 'minor_mismatch_score_unavailable_or_tied',
                                            f'segment[{idx_first}]=({first.start},{first.end}) overlaps '
                                            f'segment[{idx_second}]=({second.start},{second.end}) '
                                            f'overlap_aa={overlap_len} mismatches={n_mismatches} '
                                            f'score_first={score_first} score_second={score_second}')
                            tracker.log_overlap(sample, idx_first, (first.start, first.end),
                                                idx_second, (second.start, second.end), overlap_len,
                                                resolution='discarded_sample_score_unavailable_or_tied',
                                                overlap_seq_1_aa=ov_seq_first_aa, overlap_seq_2_aa=ov_seq_second_aa,
                                                segment_1_seq_aa=str(Seq(first.seq).translate()),
                                                segment_2_seq_aa=str(Seq(second.seq).translate()),
                                                detail=f'{n_mismatches} aa mismatch(es) within tolerance '
                                                    f'(--overlap-mismatch-tolerance={args.overlap_mismatch_tolerance}) '
                                                    f'but score_first={score_first} score_second={score_second} '
                                                    f'- cannot determine which segment to keep')

                else:
                    # Case 3: substantially different peptide sequence -> likely conflicting
                    # exon assignment or paralog
                    if args.overlap_conflict_resolution == 'mask':
                        logging.warning(f'[{bait}] Overlapping segments in sample {sample} have different '
                                    f'peptide sequences ({n_mismatches} aa mismatches in {overlap_len}aa overlap); '
                                    f'masking overlap instead of discarding sample')
                        dropped = _mask_overlap(
                            segments, idx_first, idx_second, overlap_len,
                            ov_seq_first_aa, ov_seq_second_aa,
                            sample, tracker, 'conflicting_segment_overlap',
                            extra_detail=f'{n_mismatches} aa mismatches exceeds '
                                         f'--overlap-mismatch-tolerance ({args.overlap_mismatch_tolerance})')
                        segments_to_drop.update(dropped)
                        tracker.log(sample, 'stitching', 'conflicting_segment_overlap_masked',
                                    f'segment[{idx_first}]=({first.start},{first.end}) overlaps '
                                    f'segment[{idx_second}]=({second.start},{second.end}) '
                                    f'overlap_aa={overlap_len} mismatches={n_mismatches}; masked, sample kept')
                    else:
                        logging.error(f'[{bait}] Overlapping segments in sample {sample} have different '
                                    f'peptide sequences ({n_mismatches} aa mismatches in {overlap_len}aa overlap)')
                        proper_sample = False
                        tracker.discard(sample, 'stitching', 'conflicting_segment_overlap',
                                        f'segment[{idx_first}]=({first.start},{first.end}) overlaps '
                                        f'segment[{idx_second}]=({second.start},{second.end}) '
                                        f'overlap_aa={overlap_len} mismatches={n_mismatches}')
                        tracker.log_overlap(sample, idx_first, (first.start, first.end),
                                            idx_second, (second.start, second.end), overlap_len,
                                            resolution='discarded_sample_conflicting_sequence',
                                            overlap_seq_1_aa=ov_seq_first_aa, overlap_seq_2_aa=ov_seq_second_aa,
                                            segment_1_seq_aa=str(Seq(first.seq).translate()),
                                            segment_2_seq_aa=str(Seq(second.seq).translate()),
                                            detail=f'{n_mismatches} aa mismatches exceeds '
                                                f'--overlap-mismatch-tolerance ({args.overlap_mismatch_tolerance})')

        if segments_to_drop:
            segments = [s for k, s in enumerate(segments) if k not in segments_to_drop]


# ============================================================================
# STEP D: FINAL SEQUENCE GENERATION
# ============================================================================

    # Sort segments again by start position
    segments = sorted(segments, key = lambda x: x.start)

    if proper_sample:
        # Stitch segments into the exon-only nucleotide sequence
        try:
            nt_seq_str = utils.stitch_segments(segments, bait_len, unk_ch = args.unknown)

        except AssertionError as error:
            proper_sample = False
            tracker.discard(sample, 'stitching', 'stitch_segments_assertion_failed',
                            f'{error}; bait_len={bait_len}; '
                            f'n_segments={len(segments)}; '
                            f'segments={[(s.start, s.end, len(s.seq)) for s in segments]}')

        if proper_sample:
            nt_seq_rec = SeqRecord(Seq(nt_seq_str), id = sample, description = '')
            seq_records.append(nt_seq_rec)
            segments_without_gaps = [Segment(s.seq.replace('-', ''), s.start, s.end)
                                    for s in segments]

            # Insert introns and mapping buffer
            try:
                seq_with_introns, exintr_tbl = utils.insert_introns_and_buffer(
                    segments_without_gaps, introns, bait_len, split_coords,
                    intronref_lens, included_introns, buf_len = args.nbuffer,
                    unk_ch = args.unknown)
            except AssertionError as error:
                proper_sample = False
                segment_inconsistencies = [
                    (s.start, s.end, (s.end - s.start) * CODON_SIZE, len(s.seq))
                    for s in segments_without_gaps
                    if len(s.seq) != (s.end - s.start) * CODON_SIZE
                ]
                tracker.discard(sample, 'stitching', 'insert_introns_assertion_failed', f'{error}; '
                                f'n_segments={len(segments_without_gaps)}; '
                                f'segments={[(s.start, s.end, len(s.seq)) for s in segments_without_gaps]}; '
                                f'n_introns={len(introns)}; '
                                f'introns={[(i.pos, i.length) for i in introns]}; '
                                f'n_included_introns={len(included_introns)}; '
                                f'included_introns={included_introns}; '
                                f'split_coords={split_coords}; '
                                f'n_segment_inconsistencies={len(segment_inconsistencies)}; '
                                f'segment_inconsistencies={segment_inconsistencies}; ')
            else:
                seq_with_introns_rec = SeqRecord(Seq(seq_with_introns), id = bait, description = '')
                outfile = os.path.join(exintr_for_mapping_dir, sample, f'{bait}.fasta')
                os.makedirs(os.path.dirname(outfile), exist_ok = True)
                SeqIO.write(seq_with_introns_rec, outfile, format = 'fasta')
                with open(exintr_tbl_file, 'a') as f:
                    for segment_name, start, end in exintr_tbl:
                        f.write(f'{sample}\t{bait}\t{segment_name}\t{start + 1}:{end}\n')
                tracker.keep(sample, note=f'{len(segments)} segment(s) stitched, {len(seq_with_introns)}bp written')

    # proper_sample is False but no explicit discard()
    # call fired above for this sample in this run (shouldn't normally
    # happen — every path that sets proper_sample=False also calls
    # tracker.discard — but this is a safety net so the summary/report
    # never silently omits a discarded sample).
    if not proper_sample:
        if not tracker.sample_status.get(sample, '').startswith('discarded'):
            tracker.discard(sample, 'stitching', 'proper_sample_false_unattributed',
                            'proper_sample was False but no specific event was logged - '
                            'check code paths for a missed tracker.discard() call')

# Create final output directory
outdir = args.output
os.makedirs(outdir, exist_ok = True)

# Write output stitched FASTA file
seq_out_file = os.path.join(outdir, f'{bait}.fasta')
SeqIO.write(seq_records, seq_out_file, format = 'fasta')


# ============================================================================
# Write discard report + human-readable summary to log
# ============================================================================
report_dir = args.discardreport if args.discardreport else os.path.join(exintr_for_mapping_dir, 'tblout')
events_path = os.path.join(report_dir, f'discard_events_{bait}.tsv')
summary_path = os.path.join(report_dir, f'discard_summary_{bait}.tsv')
overlap_path = os.path.join(report_dir, f'overlap_events_{bait}.txt')
tracker.write_report(events_path, summary_path)
tracker.write_overlap_report(overlap_path)
tracker.log_summary()
logging.info(f'Discard event log written to: {events_path}')
logging.info(f'Discard summary written to:   {summary_path}')
logging.info(f'Overlap event log written to: {overlap_path}')
