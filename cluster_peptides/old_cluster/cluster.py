from __future__ import print_function

import datafile


MIN_OVERLAP = 6


def read_peptides_from_csv(fname):
    """
    Read peptide sequences, sorts by shortest first, and
    removes repeats.
    """
    peptides = []
    for entry in datafile.read_csv(fname):
        sequence = entry['Sequence']
        if not sequence:
            continue
        peptides.append({
            'sequence': sequence,
            'modifications': entry['Modifications'],
            'protein': entry['Protein ID'],
            'overlaps': [],
            'groups': [],
        })
    peptides.sort()
    peptides.sort(key=lambda p: len(p['sequence']))
    for i, peptide in enumerate(peptides):
        peptide['i_peptide'] = i
    return peptides


def get_overlap(text1, text2):
    x = min(len(text1), len(text2))
    while x > 0:
        if text1[-x:] == text2[:x]:
            break
        x -= 1
    return x


def test_overlap(s1, s2, MIN_OVERLAP):
    if s1 in s2:
        return True
    if s2 in s1:
        return True
    if get_overlap(s1, s2) >= MIN_OVERLAP:
        return True
    if get_overlap(s2, s1) >= MIN_OVERLAP:
        return True
    return False


def identify_overlaps(peptides):
    n = len(peptides)
    for i in range(n):
        min_seq = peptides[i]['sequence']
        if i % 50 == 0:
            print(i, min_seq)
        for j in range(i + 1, n):
            test_seq = peptides[j]['sequence']
            if test_overlap(min_seq, test_seq, MIN_OVERLAP):
                peptides[i]['overlaps'].append(j)
                peptides[j]['overlaps'].append(i)


def find_groups(peptides):

    groups = []

    def tag_group_of_peptide(i_peptide, i_group):
        peptide = peptides[i_peptide]
        seq = peptide['sequence']
        group = groups[i_group]

        if i_group in peptide['groups']:
            return

        for group_peptide in group['peptides']:
            group_seq = group_peptide['sequence']
            if not test_overlap(seq, group_seq, MIN_OVERLAP):
                return

        peptide['groups'].append(i_group)
        group['peptides'].append(peptide)
        for i_overlap_peptide in peptide['overlaps']:
            tag_group_of_peptide(i_overlap_peptide, i_group)

    def new_group(i_peptide):
        i_group = len(groups)
        groups.append({
            'i': i_group,
            'peptides': []
        })
        tag_group_of_peptide(i_peptide, i_group)

    n_peptide = len(peptides)
    for i_peptide in range(n_peptide):
        peptide = peptides[i_peptide]
        if len(peptide['groups']) == 0:
            new_group(i_peptide)

    return groups


def make_groups(peptides_csv, groups_yaml):
    peptides = read_peptides_from_csv(peptides_csv)
    identify_overlaps(peptides)
    groups = find_groups(peptides)
    datafile.write_yaml(groups, groups_yaml)


def is_misalign(s1, s2):
    misalignment = False
    n_char = max(map(len, [s1, s2]))
    for i_char in range(n_char):
        for s in [s1, s2]:
            if i_char == len(s):
                return False
        chars_in_col = set(s[i_char] for s in [s1, s2])
        if ' ' in chars_in_col:
            chars_in_col.remove(' ')
        if len(chars_in_col) > 1:
            return True
    return False


def get_i_relative_to_ref(ref_seq, seq):
    ref_left_overlap = get_overlap(ref_seq, seq)
    ref_right_overlap = get_overlap(seq, ref_seq)
    if seq in ref_seq:
        return ref_seq.find(seq)
    elif ref_seq in seq:
        return -seq.find(ref_seq)
    elif ref_left_overlap:
        return len(ref_seq) - ref_left_overlap
    elif ref_right_overlap:
        return -(len(seq) - ref_right_overlap)
    raise "No overlap"


print("Generating groups...")
make_groups('b57_clean.csv', 'groups.yaml')

print("Loading groups...")
groups = datafile.load_yaml('groups.yaml')

rows = [('i_group', 'sequence', 'modifications', 'protein')]

for group in groups:
    sequences = []
    for peptide in group['peptides']:
        sequences.append(peptide['sequence'])

    ref_seq = sequences[0]
    indices = [get_i_relative_to_ref(ref_seq, s) for s in sequences]

    # find max length
    left = -min(indices)
    max_len = 0
    for index, seq in zip(indices, sequences):
        str_len = left + index + len(seq)
        if str_len > max_len:
            max_len = str_len

    for index, peptide in zip(indices, group['peptides']):
        row = []

        row.append(group['i'])
        row.append(peptide['sequence'])
        row.append(peptide['modifications'])
        row.append(peptide['protein'])

        row.append(' ')

        seq = peptide['sequence']
        s = ' ' * (left + index) + seq
        s = s + ' ' * (max_len - len(s))
        for c in s:
            if c == ' ':
                row.append('.')
            else:
                row.append(c)
        rows.append(row)

    rows.append([])

datafile.write_csv(rows, 'cluster.csv')
