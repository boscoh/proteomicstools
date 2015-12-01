from __future__ import print_function
import os

import datafile


__doc__ = """

Sorts a bunch of peptides into clusters of overlapping residues
and identifies the minimum peptide of each cluster.


"""

def read_peptides_from_csv(fname):
    """
    Read peptide sequences, sorts by shortest first, and
    removes repeats.
    """
    peptides = []
    for entry in datafile.read_csv(fname):
        if 'Sequence' not in entry:
            continue
        peptide = {
            'sequence': entry['Sequence'],
            'overlaps': [],
            'subsets': [],
            'supersets': [],
            'groups': [],
            'attr': {}
        }
        for key in entry:
            if not key == 'Sequence':
                peptide['attr'][key] = entry[key]
        peptides.append(peptide)
    peptides.sort()
    peptides.sort(key=lambda p: len(p['sequence']))
    for i, peptide in enumerate(peptides):
        peptide['i_peptide'] = i
    return peptides


def read_peptides_from_txt(fname):
    """
    Read peptide sequences, sorts by shortest first, and
    removes repeats.
    """
    peptides = []
    for sequence in open(fname).read().split():
        peptides.append({
            'sequence': sequence,
            'overlaps': [],
            'subsets': [],
            'supersets': [],
            'groups': [],
            'attr': {}
        })
    peptides.sort()
    peptides.sort(key=lambda p: len(p['sequence']))
    for i, peptide in enumerate(peptides):
        peptide['i_peptide'] = i
    return peptides


def get_i_relative_to_ref(ref_seq, seq):
    """
    Identifies the position of seq relative to ref_seq.
    """
    if seq in ref_seq:
        return ref_seq.find(seq)
    
    if ref_seq in seq:
        return -seq.find(ref_seq)

    ref_left_overlap = get_overlap(ref_seq, seq)
    ref_right_overlap = get_overlap(seq, ref_seq)
    
    if ref_left_overlap and ref_left_overlap > ref_right_overlap:
        return len(ref_seq) - ref_left_overlap
    
    if ref_right_overlap and ref_right_overlap > ref_left_overlap:
        return -(len(seq) - ref_right_overlap)

    raise "No overlap"


def calc_alignment(group):
    """
    Creates an 'aligned_sequence' entry for each peptide
    in group that shows the alignment.
    """

    peptides = group['peptides']

    sequences = []
    for peptide in peptides:
        sequences.append(peptide['sequence'])

    ref_seq = sequences[0]
    if 'shortest' in group:
        ref_seq = group['shortest']['sequence']
    
    indices = [
        get_i_relative_to_ref(ref_seq, s) 
        for s in sequences
    ]

    # find max length
    left = -min(indices)
    max_len = 0
    for index, seq in zip(indices, sequences):
        str_len = left + index + len(seq)
        if str_len > max_len:
            max_len = str_len

    n_terminal_aligned = True
    c_terminal_aligned = True

    for index, peptide in zip(indices, peptides):
        seq = peptide['sequence']
        n_pre = left + index
        n_post = max_len - len(seq) - n_pre
        s = '.'*n_pre + seq + '.'*n_post
        peptide['aligned_sequence'] = s
        if s[0] == '.':
            n_terminal_aligned = False
        if s[-1] == '.':
            c_terminal_aligned = False

    if n_terminal_aligned and c_terminal_aligned:
        group['alignment'] = 'Single'
    elif n_terminal_aligned:
        group['alignment'] = 'N-aligned'
    elif c_terminal_aligned:
        group['alignment'] = 'C-aligned'
    else:
        group['alignment'] = 'Mixed-aligned'

    sequences = []
    for peptide in peptides:
        sequences.append(peptide['aligned_sequence'])
    n = len(sequences[0])

    # check sequenes
    for seq in sequences:
        if len(seq) != n:
            print('>>problem at calc_alignments')
            for seq2 in sequences:
                print(seq2)

    for peptide in peptides:
        assert len(peptide['aligned_sequence']) == max_len


def check_misalignment(group):
    """
    Check that 'aligned_sequence' entries are same length,
    and that every column has same character.
    """
    
    sequences = [
        p['aligned_sequence'] 
        for p in group['peptides']
    ]

    n = len(sequences[0])

    # check sequenes
    for seq in sequences:
        if len(seq) != n:
            print('>>problem at check_misalignment')
            for seq2 in sequences:
                print(seq2)

    misaligned = False
    for i in range(n):
        chars = set()
        for seq in sequences:
            if i < len(seq):
                if seq[i] != '.':
                    chars.add(seq[i])
        if len(chars) > 1:
            misaligned = True
            break

    group['misaligned'] = misaligned


def write_groups(groups, out_csv):

    headers = ['group', 'other_groups', 'align', 'consistent', 'sequence']

    titles = groups[0]['peptides'][0]['attr'].keys()

    headers.extend(titles)

    rows = [headers]

    for group in groups:

        calc_alignment(group)

        check_misalignment(group)

        for peptide in group['peptides']:
            row = []

            row.append(group['i'])

            other_groups = []
            for i_group in peptide['groups']:
                if i_group != group['i']:
                    other_groups.append(str(i_group))

            if other_groups:
                other_groups = ';'.join(other_groups)
            else:
                other_groups = ''

            row.append(other_groups)

            row.append(group['alignment'])
            row.append(not group['misaligned'])

            row.append(peptide['sequence'])

            for title in titles:
                if title in peptide['attr']:
                    row.append(peptide['attr'][title])
                else:
                    row.append('')

            row.append(' ')

            for c in peptide['aligned_sequence']:
                row.append(c)

            rows.append(row)

        rows.append([
            group['i'],
            '',
            group['alignment'],
            not group['misaligned'],
        ])

    datafile.write_csv(rows, out_csv)


##################
# min_overlap


def get_overlap(text1, text2):
    x = min(len(text1), len(text2))
    while x > 0:
        if text1[-x:] == text2[:x]:
            break
        x -= 1
    return x


def test_overlap(s1, s2, min_overlap):
    if s1 in s2:
        return True
    if s2 in s1:
        return True
    if get_overlap(s1, s2) >= min_overlap:
        return True
    if get_overlap(s2, s1) >= min_overlap:
        return True
    return False


def identify_overlaps(peptides, min_overlap):
    n = len(peptides)
    for i in range(n):
        min_seq = peptides[i]['sequence']
        if i % 50 == 0:
            print("Processed", i, "sequences...")
        for j in range(i + 1, n):
            test_seq = peptides[j]['sequence']
            if test_overlap(min_seq, test_seq, min_overlap):
                peptides[i]['overlaps'].append(j)
                peptides[j]['overlaps'].append(i)


def find_overlap_groups(peptides, min_overlap):

    groups = []

    def tag_group_of_peptide(i_peptide, i_group):
        peptide = peptides[i_peptide]
        seq = peptide['sequence']
        group = groups[i_group]

        if i_group in peptide['groups']:
            return

        for group_peptide in group['peptides']:
            group_seq = group_peptide['sequence']
            if not test_overlap(seq, group_seq, min_overlap):
                return

        peptide['groups'].append(i_group)
        group['peptides'].append(peptide)
        for i_overlap_peptide in peptide['overlaps']:
            tag_group_of_peptide(i_overlap_peptide, i_group)

    def new_group(i_peptide):
        i_group = len(groups)
        groups.append({
            'i': i_group,
            'shortest': peptides[i_peptide],
            'peptides': []
        })
        tag_group_of_peptide(i_peptide, i_group)

    n_peptide = len(peptides)
    for i_peptide in range(n_peptide):
        peptide = peptides[i_peptide]
        if len(peptide['groups']) == 0:
            new_group(i_peptide)

    return groups


def make_overlap_groups(peptides, groups_yaml, min_overlap):
    identify_overlaps(peptides, min_overlap)
    groups = find_overlap_groups(peptides, min_overlap)
    datafile.write_yaml(groups, groups_yaml)


def get_kernel(group):
    sequences = []
    for peptide in group['peptides']:
        sequences.append(peptide['sequence'])

    ref_seq = sequences[0]
    indices = [get_i_relative_to_ref(ref_seq, s) for s in sequences]

    def all_aligned(i):
        for i_seq_rel_ref, seq in zip(indices, sequences):
            i_rel_to_seq = i - i_seq_rel_ref
            if i_rel_to_seq < 0 or i_rel_to_seq >= len(seq):
                return False
        return True

    i_start = min(indices)
    while not all_aligned(i_start):
        i_start += 1

    i_end = i_start
    while all_aligned(i_end):
        i_end += 1

    return ref_seq[i_start:i_end]


def write_overlap_kernels(groups, out_csv):

    headers = ['group', 'size_group', 'sequence']

    rows = [headers]

    for group in groups:
        rows.append([
            group['i'],
            len(group['peptides']),
            get_kernel(group)
        ])

    datafile.write_csv(rows, out_csv)


def process_overlaps(in_f, min_overlap=6):

    print("Reading peptides from %s.." % in_f)
    if in_f.endswith('.csv'):
        peptides = read_peptides_from_csv(in_f)
    elif in_f.endswith('.txt'):
        peptides = read_peptides_from_txt(in_f)

    print("Generating groups...")
    make_overlap_groups(peptides, 'groups.yaml', min_overlap)

    groups = datafile.load_yaml('groups.yaml')
    base = os.path.splitext(in_f)[0]

    out_csv = base + '.cluster.csv'
    print("Writing groups %s" % out_csv)
    write_groups(groups, out_csv)

    out_csv = base + '.kernel.csv'
    print("Writing groups %s" % out_csv)
    write_overlap_kernels(groups, out_csv)


##################
# subset_groups


def identify_subsets(peptides):
    n = len(peptides)
    for i in range(n):
        peptide_i = peptides[i]['sequence']
        if i % 50 == 0:
            print("Processed", i, "sequences...")
        for j in range(n):
            if i == j:
                continue
            peptide_j = peptides[j]['sequence']
            if peptide_i in peptide_j:
                peptides[i]['supersets'].append(j)
                peptides[j]['subsets'].append(i)


def find_subset_groups(peptides):
    """
    This uses a recursive procedure to generate
    groups, as the first group is created, anything
    that doesn't fit is used to generate a new
    group.
    """

    groups = []

    def tag_group_of_peptide(i_peptide, i_group):
        peptide = peptides[i_peptide]
        seq = peptide['sequence']
        group = groups[i_group]

        if i_group in peptide['groups']:
            return

        min_seq = group['shortest']['sequence']
        if min_seq not in seq:
            return

        peptide['groups'].append(i_group)
        group['peptides'].append(peptide)

        for i_superset_peptide in peptide['supersets']:
            tag_group_of_peptide(i_superset_peptide, i_group)

    def new_group(i_peptide):
        i_group = len(groups)
        groups.append({
            'i': i_group,
            'shortest': peptides[i_peptide],
            'peptides': []
        })
        tag_group_of_peptide(i_peptide, i_group)

    n_peptide = len(peptides)
    for i_peptide in range(n_peptide):
        peptide = peptides[i_peptide]
        if len(peptide['groups']) == 0:
            new_group(i_peptide)

    return groups


def make_subset_groups(peptides, groups_yaml):
    identify_subsets(peptides)
    groups = find_subset_groups(peptides)
    datafile.write_yaml(groups, groups_yaml)


def write_subset_kernels(groups, out_csv):

    headers = ['group', 'size_group', 'alignment', 'sequence']

    rows = [headers]

    for group in groups:
        rows.append([
            group['i'],
            len(group['peptides']),
            group['alignment'],
            group['shortest']['sequence'],
        ])

    datafile.write_csv(rows, out_csv)


def process_subsets(in_f):
    """
    need to look for redundancies
    """

    print("Reading peptides from %s.." % in_f)
    if in_f.endswith('.csv'):
        peptides = read_peptides_from_csv(in_f)
    elif in_f.endswith('.txt'):
        peptides = read_peptides_from_txt(in_f)

    print("Generating groups...")
    make_subset_groups(peptides, 'groups.yaml')

    groups = datafile.load_yaml('groups.yaml')
    base = os.path.splitext(in_f)[0]

    out_csv = base + '.cluster.csv'
    print("Writing groups %s" % out_csv)
    write_groups(groups, out_csv)

    out_csv = base + '.kernel.csv'
    print("Writing groups %s" % out_csv)
    write_subset_kernels(groups, out_csv)


##################


# process_overlaps('peptides_ralf_aug15.txt')
# process_overlaps('IAb_B6_spleen.csv')
# process_subsets('pat_20150522_Collate_B57_clean.csv')
# process_subsets('peptides_ralf_aug15.txt')
# process_subsets('test.txt')

