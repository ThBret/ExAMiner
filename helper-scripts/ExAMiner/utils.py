import Bio
from collections import Counter
import numpy as np
import re

CODON_SIZE = 3

BLOSUM62 = Bio.Align.substitution_matrices.load('BLOSUM62')
score_mat = np.array(BLOSUM62)
score_mat = np.append(score_mat, np.ones((score_mat.shape[0], 1)) * -1, axis=1)
score_mat = np.append(score_mat, np.ones((1, score_mat.shape[1])) * -1, axis=0)
BLOSUM62_WITH_GAPS = Bio.Align.substitution_matrices.Array(
    BLOSUM62.alphabet + '-', dims=2, data=score_mat, dtype=int)
## No penalty when two gaps match
BLOSUM62_WITH_GAPS['-', '-'] = 0

GAP_PAT = re.compile('-+')

## NOTE: make sure that the exonerate output files have three-letter
## abbreviations for amino acid names, not one-letter symbols.
def aln_score(seq1, seq2, score_mat=BLOSUM62_WITH_GAPS):
    """Calculate the sum of scores of the given sequences (it is
    assumed that they have been aligned).

    The scores are obtained from the provided substitution matrix.

    Args:
        seq1: first sequence (a string).
        seq2: second sequence (a string).
        score_mat: the substitution matrix (the BLOSUM62 matrix by default).
    Returns:
        the sum of scores of the given sequences.
    """
    assert len(seq1) == len(seq2)

    scores = [score_mat[res_pair[0]][res_pair[1]] for res_pair in zip(seq1, seq2)]

    return sum(scores)


def sim_score(seq1, seq2, score_mat=BLOSUM62_WITH_GAPS):
    """Calculate the similarity score of the given sequences (it is
    assumed that they have been aligned).

    The scores are obtained from the provided substitution matrix.

    Args:
        seq1: first sequence (a string).
        seq2: second sequence (a string).
        score_mat: the substitution matrix (the BLOSUM62 matrix by default).
    Returns:
        the similarity score of the given sequences.
    """
    assert len(seq1) == len(seq2)

    a12 = aln_score(seq1, seq2, score_mat)
    a11 = aln_score(seq1, seq1, score_mat)
    a22 = aln_score(seq2, seq2, score_mat)

    if a11 == 0 or a22 == 0 or a12 <= 0:
        return 0

    return a12 / np.sqrt(a11 * a22)


def overlap(coord1, coord2):
    """Calculate amount of overlap between two segments.

    The segment coordinates are given as a list or tuple of two
    elements (beginning, end).  The end coordinate is not included in
    the segment, as when indexing lists or strings.

    Args:
        coord1: first segment coordinates.  A list or tuple of two elements (beginning, end).
        coord2: second segment coordinates.  A list or tuple of two elements (beginning, end).
    Returns:
        the length of the overlap between the given segments (an integer).
    """
    maxstart = max(coord1[0], coord2[0])
    minend = min(coord1[1], coord2[1])
    return max(0, minend - maxstart)


def overlaps(coord1, coord2):
    """Return True if the provided segments overlap, False otherwise.

    Args:
        coord1: first segment coordinates.  A list or tuple of two integers (beginning, end).
        coord2: second segment coordinates.  A list or tuple of two integers (beginning, end).
    Returns:
        a boolean value.
    """
    return overlap(coord1, coord2) > 0


def overlaps_point(segment, point):
    """Return True if the provided point (position) is within the
    segment coordinates.

    Args:
        segment: segment coordinates.  A list or tuple of two integers (beginning, end).
        point: point coordinate (position). An integer.
    Returns:
        a boolean value.
    """
    ## Overlap does not count at one of the ends
    if point > segment[0] and point < segment[1]:
        return True
    else:
        return False


def exon_overlap(segm_coord, exon_coord):
    """Returns an array of the overlaps of the specified segment with
    each of the provided exons in the argument `exon_coord'.

    The overlap with the first exon is in index 0 of the returned
    array, and so on.

    Args:
        segm_coord: segment coordinates.  A list or tuple of two integers (beginning, end).
        exon_coord: exon coordinates.  A list of segments.
    Returns:
       an array of the overlaps of the specified segment with each of
       the provided exons in the argument `exon_coord'.
    """
    ov = []
    for i in range(len(exon_coord)):
        ov.append(overlap(segm_coord, exon_coord[i]))
    return ov


def stitch_segments(segments, bait_end, st=0, unk_ch='N'):
    """Stitch the provided segments in one sequence, filling in
    missing positions in the bait.
    
    Args:
        segments: a list of Segment objects.
        bait_end: the end coordinate of the bait, where the stitching ends (an integer).
        st: the start coordinate, where the stitching starts (an integer).
        unk_ch: character used to fill unknown or missing nucleotides
    Returns:
        a sequence with the stitched segments.
    """
    segments_sort = sorted(segments, key=lambda x: x.start)
    prev_end = st
    stitched_seq = ''
    for segm in segments_sort:
        curr_st = segm.start
        assert not prev_end > curr_st
        ## Fill with 'N's or '?'s
        stitched_seq += unk_ch * (curr_st - prev_end) * CODON_SIZE
        stitched_seq += segm.seq
        prev_end = segm.end
    assert not prev_end > bait_end
    stitched_seq += unk_ch * (bait_end - prev_end) * CODON_SIZE
    return stitched_seq


def normalise_gaps(gaps):
    """Normalise gap coordinates accounting for previous gaps.

    Args:
        gaps: a list of gaps, which are tuples (position, length).
    Returns:
        A list of gaps with normalised positions, i.e. the position of
        the amino acid of the bait in which they are inserted.
    """
    if not gaps:
        return []
    positions, lengths = zip(*gaps)
    positions = np.array(positions)
    lengths = np.array(lengths)
    cumulative_len = np.cumsum(lengths)
    ## Subtract the length of all previous gaps to normalise their positions
    positions[1:] -= cumulative_len[:-1]
    return list(zip(positions, lengths))


def denormalise_gaps(gaps):
    """Denormalise gap coordinates accounting for previous gaps.

    Args:
        gaps: a list of normalised gaps, which are tuples (position,
        length).
    Returns:
        A list of gaps with denormalised positions, i.e. their
        position in a sequence in which they are already inserted.
    """
    if not gaps:
        return []
    positions, lengths = zip(*gaps)
    positions = np.array(positions)
    lengths = np.array(lengths)
    cumulative_len = np.cumsum(lengths)
    positions[1:] += cumulative_len[:-1]
    return list(zip(positions, lengths))


def add_gaps_aa_seq(seq, gaps):
    """Adds gaps to an amino acid sequence in the specified positions.

    Args:
       seq: the original sequence (a string).
       gaps: a list of gaps, which are tuples (position, length).
    Returns:
       a string representing the sequence with the gaps inserted in
       the specified positions.
    """
    gaps = denormalise_gaps(gaps)
    seq_with_gaps = seq
    for p, l in gaps:
        ## The '=' in '<=' allows to add gaps at the very end,
        ## although this should not be necessary
        assert p <= len(seq_with_gaps)
        seq_with_gaps = seq_with_gaps[:p] + '-' * l + seq_with_gaps[p:]
    return seq_with_gaps


class Segment:
    """A class that represents a part of a sequence that is to be stitched.
    """
    def __init__(self, seq, start, end):
        assert start <= end
        self.seq = seq
        self.start = start
        self.end = end
        self.length = end - start

    def __repr__(self):
        return f'<Segment of length {len(self.seq)} nt, start: {self.start}, end: {self.end}>'

class Intron:
    """A class representing an intron."""
    def __init__(self, length, pos):
        self.length = length
        self.pos = pos
        self.split = pos[1] - pos[0] > 0

    def __repr__(self):
        if not self.split:
            spl_str = 'No split amino acid.'
        else:
            spl_str = f'Has a split amino acid.'
        return f'<Intron of length {self.length} in position {self.pos}. {spl_str}>'


class AlignmentReference:
    """A class representing a reference alignment."""
    def __init__(self, hsp, bait_seq):
        query_seq = ''.join([str(q.seq.rstrip('X')) for q in hsp.query_all])
        gap_pos = gap_pos_from_seq(query_seq)
        gap_pos = shift_gap_pos_by(gap_pos, hsp.query_range[0])
        query_seq_with_gaps = add_gaps_aa_seq(bait_seq, gap_pos)
        hit_nt_seq = Bio.Seq.Seq(''.join([''.join(ann['hit_annotation'])
                                          for ann in hsp.aln_annotation_all]))
        hit_aa_seq = str(hit_nt_seq.translate())

        if len(hit_aa_seq) < len(query_seq_with_gaps):
            hit_aa_seq = pad_aa_seq(hit_aa_seq, hsp.query_range[0],
                                    len(query_seq_with_gaps))

        assert len(hit_aa_seq) == len(query_seq_with_gaps)

        self.hsp = hsp
        self.query_span = hsp.query_span
        self.split_coords = tuple(hsp.query_inter_ranges)
        self.exon_coords = tuple(hsp.query_range_all)
        self.intron_lens = tuple(hsp.hit_inter_spans)
        self.ref_aln = [query_seq_with_gaps, hit_aa_seq]
        self.nt_seq = hit_nt_seq

    def __repr__(self):
        return f'<Alignment reference with exon coordinates {self.exon_coords}>'

    def ref_aln_score(self):
        return sim_score(self.ref_aln[0], self.ref_aln[1])


def get_segment_from_aln(aln, crd):
    """Get the segment of a given alignment between the given
    coordinates.

    Args:
        aln: the alignment.
        crd: the coordinates (beginning, end) of the part of the
        alignment to be obtained.  These coordinates correspond to
        amino acid positions in the bait.

    Returns:
        The segment of the alignment between the coordinates in `crd'.
    """
    query_seq = aln[0]
    gaps = gap_pos_from_seq(query_seq)

    st = crd[0]
    st += sum([l for p, l in gaps if p <= st])
    end = crd[1]
    end += sum([l for p, l in gaps if p < end])
    cropped_aln = [s[st:end] for s in aln]
    return cropped_aln


def gap_pos_from_seq(seq, pat=GAP_PAT):
    """Get the gaps from a given sequence.

    Args:
        seq: the sequence from which to get gap positions.
        pat: regexp indicating what is to be considered a gap.

    Returns:
        The position and length of gaps in sequence `seq'.
    """
    gaps = []
    for m in pat.finditer(seq):
        s = m.span()
        gaps.append((s[0], s[1] - s[0]))
    gaps = normalise_gaps(gaps)
    return gaps


def insert_introns_and_buffer(segments, introns, bait_end,
                              ref_spl_crd, ref_intr_len, included,
                              buf_len=1000, unk_ch='N'):
    """Generates the output for mapping.

    This function inserts introns of the appropriate length where they
    appear.  If the exonerate alignment does not cover the intron
    position, it inserts an intron with the mean length of introns in
    the corresponding position in the reference alignments.  It also
    adds buffer nucleotides at both ends of the sequence.

    Args:
        segments: a list of Segment objects.
        introns: a list of Intron objects.
        bait_end: the end coordinate of the bait (an integer).
        ref_spl_crd: reference split coordinates (exon-intron boundaries).
        ref_intr_len: mean length of introns in reference alignments (list of integers).
        included: which introns have already been included in the
            exonerate alignment (a list of indices).
        buf_len: length of the buffer that is to be added to both ends
            (an integer, default: 1000).
        unk_ch: character to use when a nucleotide is unknown, and
            in the buffer (default: 'N').
    Returns:
        - a string representing the sequence with the stitched segments,
        with introns inserted in the appropriate positions and buffers
        at both ends.
        - a table specifying the coordinates of the buffers, exons,
          and introns (1-based indexing, end coordinate included in
          the segment)
    """
    seq = ''
    tbl = []

    ref_intr_to_add = [Intron(ref_intr_len[i], pos)
                       for i, pos in enumerate(ref_spl_crd)
                       if i not in included]
    intr_to_add = introns + ref_intr_to_add
    intr_to_add = sorted(intr_to_add, key=lambda x: x.pos[0])

    cur_pos = 0
    # Start buffer
    buf = unk_ch * buf_len
    seq += buf
    tbl.append(['ftail', cur_pos, len(buf)])
    cur_pos += len(buf)
    last_start = 0
    for intr_idx, intron in enumerate(intr_to_add):
        intr_start = intron.pos[0]
        intr_end = intron.pos[1]
        ## This assumes that segments will not overlap introns (which
        ## should be true if things are cropped properly)
        segm_between = [s for s in segments if s.start >= last_start and s.end <= intr_start]
        exon_seq = stitch_segments(segm_between, intr_start, last_start, unk_ch)
        seq += exon_seq
        tbl.append([f'exon{intr_idx + 1}', cur_pos,
                    cur_pos + len(exon_seq)])
        cur_pos += len(exon_seq)

        if intron.split:
            ## NOTE: I had to give up on getting the correct split
            ## site because the coordinates in `query_split_codons' are
            ## inconsistent in some cases.  There might be a different
            ## way to get this information reliably.
            ##
            ## Add the split amino acid to the intron
            intr_len = intron.length + 3
        else:
            intr_len = intron.length

        seq += unk_ch * intr_len
        tbl.append([f'intron{intr_idx + 1}', cur_pos, cur_pos + intr_len])
        cur_pos += intr_len
            
        last_start = intron.pos[1]
    ## Get the remaining segments until the end
    remaining_segments = [s for s in segments if s.start >= last_start]
    exon_seq = stitch_segments(remaining_segments, bait_end, last_start, unk_ch)
    seq += exon_seq
    tbl.append([f'exon{len(intr_to_add) + 1}', cur_pos,
                cur_pos + len(exon_seq)])
    cur_pos += len(exon_seq)
    # End buffer
    seq += unk_ch * buf_len
    tbl.append(['etail', cur_pos, cur_pos + len(buf)])
    return seq, tbl


def get_count(ref_dict, attr):
    """Get the coordinate count from a certain type in a
    dictionary of possible references.

    Args:
        ref_dict: dictionary of possible references (of the
        AlignmentReference class), its keys being the names of the
        samples.
        attr: type of coordinates. Possible values are 'split_coords'
        or 'exon_coords'.
    Returns:
        The coordinate count in the provided 'ref_dict' array.
    """
    coords = [getattr(r, attr) for _, r in ref_dict.items()]
    return Counter(coords)



def shift_gap_pos_by(gap_pos, by):
    """Shift gap_pos by a number of positions.

    This is useful with partial references that start after position 0
    of the bait.

    Args:
        gap_pos: the original gap positions.  They consist of an array
        of (position, length) tuples.
        by: number of positions to shift
    Returns:
        An array with the modified positions of the gaps.
    """
    shift_gap_pos = [(pos + by, l) for pos, l in gap_pos]
    return shift_gap_pos


def pad_aa_seq(aa_seq, pos, tot_len):
    """Pad an amino acid sequence with 'X's to match the length of a
    longer amino acid sequence.

    Allows to align both sequences.

    Args:
        aa_seq: a string with the amino acid sequence to pad.
        pos: the start position of the alignment of 'aa_seq' to the
        longer sequence.
        tot_len: the length of the longer amino acid sequence.
    Returns:
        'aa_seq' padded with the necessary number of 'X's to match the
        longer sequence.
    """
    last_pos = pos + len(aa_seq)
    return 'X' * pos + aa_seq + 'X' * (tot_len - last_pos)
    

def crop_fragment_seq(frg, crop_crd):
    """Crop a fragment to the reference coordinates.

    Args:
        frg: the fragment to crop.
        crop_crd: the coordinates to which the fragment has to be cropped.
    Returns:
        The cropped nucleotide sequence of the fragment.
        The position of the amino acid in the bait where the cropped sequence starts.
        The position of the amino acid in the bait where the cropped sequence ends.
    """
    sgm_st = max(crop_crd[0], frg.query_start)
    sgm_end = min(crop_crd[1], frg.query_end)
    frg_st = sgm_st - frg.query_start
    frg_end = sgm_end - frg.query_start
    # Take bait gaps into account
    gap_pos = gap_pos_from_seq(str(frg.query.seq).strip('X'))
    frg_st += sum([l for pos, l in gap_pos if pos <= frg_st])
    frg_end += sum([l for pos, l in gap_pos if pos < frg_end])
    crop_st = frg_st * CODON_SIZE
    crop_end = frg_end * CODON_SIZE

    ## Get nucleotide sequence of the fragment
    codons = frg.aln_annotation['hit_annotation']
    # Remove nucleotides that belong to split amino acids
    codons = [cod for cod in codons if len(cod) == 3]
    frg_seq = ''.join(codons)

    nt_seq = frg_seq[crop_st:crop_end]

    ## This can happen if the fragment is outside the cropping coordinates
    if sgm_st > sgm_end:
        sgm_st, sgm_end = None, None

    return nt_seq, sgm_st, sgm_end


def segments_from_aln_ref(aln_ref):
    """Get a list of Segment objects from an alignment reference.

    Args:
        aln_ref: an AlignmentReference object.
    Returns:
        a list of Segment objects from the provided reference.
    """
    segments = []
    ref_hsp = aln_ref.hsp
    hsp_query_seq = ''.join([str(f.query.seq.rstrip('X')) for f in ref_hsp])
    hsp_gap_pos = gap_pos_from_seq(hsp_query_seq)
    hsp_seq_nt = ''.join([''.join(ann['hit_annotation'])
                          for ann in ref_hsp.aln_annotation_all])
    ## Split amino acids
    for spl_st, spl_end in aln_ref.split_coords:
        ngaps = sum([l for pos, l in hsp_gap_pos if pos <= (spl_st - ref_hsp.query_start)])
        st = (spl_st - ref_hsp.query_start + ngaps) * CODON_SIZE
        spl_seq = hsp_seq_nt[st:(st + (spl_end - spl_st) * CODON_SIZE)]
        segments.append(Segment(spl_seq, spl_st, spl_end))

    ## Sequence fragments
    for fragment in ref_hsp:
        codons = fragment.aln_annotation['hit_annotation']
        codons = [cod for cod in codons if len(cod) == 3]
        frg_seq = ''.join(codons)
        segments.append(Segment(frg_seq, fragment.query_start, fragment.query_end))

    return segments


def introns_from_aln_ref(aln_ref):
    """Get a list of Intron objects from an alignment reference.

    Args:
        aln_ref: an AlignmentReference object.
    Returns:
        a list of Intron objects from the provided reference.
    """
    return [Intron(l, p) for l, p in zip(aln_ref.intron_lens, aln_ref.split_coords)]


## IMPORTANT: spl_crd should be sorted before calling this function
def exon_crd_from_split_crd(spl_crd, beg, end):
    """Get a list of exon coordinates for a bait given the position of
    the splice sites and the start and end points of the positions
    covered by a contig.

    Args:
        spl_crd: the coordinates of the splice sites.  A list of (beginning, end) tuples.
        beg: the beginning coordinate (an integer).
        end: the end coordinate (an integer).

    Return:
        A list of exon coordinates.
    """
    if not spl_crd:
        return ((beg, end),)
    else:
        first_spl = spl_crd[0]
        assert first_spl[0] >= beg and first_spl[1] <= end
        remaining_exon_crd = exon_crd_from_split_crd(spl_crd[1:], first_spl[1], end)
        return ((beg, first_spl[0]),) + remaining_exon_crd


def get_ref_exon_coords(potential_refs):
    """Get the reference exon coordinates for the provided set of samples.

    Args:
        potential_refs: a dictionary of AlignmentReference objects,
        the keys being the sample names.
    Returns:
        a list of reference exon coordinates.
    """
    # Find most frequently occurring splice site coordinates
    spl_crd_cnt = get_count(potential_refs, 'split_coords')
    max_spl_cnt = max(spl_crd_cnt.values())
    most_common_crds = [k for k, v in spl_crd_cnt.items() if v == max_spl_cnt]

    # Filter specimens to only those with the most common splice coordinates
    common_refs = {k: r for k, r in potential_refs.items()
                   if r.split_coords in most_common_crds}
    
    # Find most frequently occurring exon coordinate structure
    exon_crd_cnt = get_count(common_refs, 'exon_coords')
    max_exon_cnt = max(exon_crd_cnt.values())
    poss_exon_crd = [k for k, v in exon_crd_cnt.items() if v == max_exon_cnt]
    assert len(poss_exon_crd) >= 1
    # Handle ties - select by longest query span
    if len(poss_exon_crd) > 1:
        ## Choose the longest reference (if there is a complete
        ## reference, all of them have the same length and it will
        ## choose one at random, which is what has been happening
        ## until now).  Complete references should only get here if
        ## there are multiple potential split coordinates appearing
        ## with the same frequency.
        max_len = max([r.query_span for r in common_refs.values()])
        refs_maxlen = {k: r for k, r in common_refs.items() if r.query_span == max_len}
        crd_cnt = get_count(refs_maxlen, 'exon_coords')
        return sorted(crd_cnt)[0]
    else:
        return poss_exon_crd[0]


def list_from_file_lines(infile):
    """Return a list whose elements are the lines of the specified
    file, without leading and trailing whitespace.

    Args:
        infile: the input file
    Returns:
        a list of strings
    """
    with open(infile) as f:
        arr = list(f)
    arr = [line.strip() for line in arr]
    return arr
