from __future__ import print_function
import csv
import uniprot
import re


def clean_seqid(seqid):
  pieces = seqid.split('|')
  return pieces[1]


def is_overlap(i1, seq1, i2, seq2):
    j1 = i1 + len(seq1)
    j2 = i2 + len(seq2)
    if i1 >= j2 or i2 >= j1:
        return False
    return True


def gap(i1, seq1, i2, seq2):
    if not is_overlap(i1, seq1, i2, seq2):
        j1 = i1 + len(seq1)
        j2 = i2 + len(seq2)
        return min(abs(i1 - j2), abs(i2 - j1))
    return None


fname = 'KL-A1_pep15.csv'

peptides = []

proteins = {}
for row in csv.DictReader(open(fname)):
    seqid = row['Protein']
    if seqid not in proteins:
        proteins[seqid] = []
    proteins[seqid].append(row)
    row['matches'] = []
    peptides.append(row)

seqids, fasta = uniprot.read_fasta('use.fasta')

for seqid, protein in proteins.items():
    bare_seqid = clean_seqid(seqid)
    full_sequence = fasta[bare_seqid]['sequence']
    for i_peptide1, peptide1 in enumerate(protein):
        seq1 = peptide1['Sequence']
        i1 = full_sequence.find(seq1)
        peptide1['Start'] = i1 + 1
        peptide1['End'] = i1 + len(seq1)
        for peptide2 in protein:
            if peptide1 == peptide2: 
                continue
            seq2 = peptide2['Sequence']
            i2 = full_sequence.find(seq2)
            match = {
                'overlap': is_overlap(i1, seq1, i2, seq2),
                'gap': gap(i1, seq1, i2, seq2),
                'seq1': seq1,
                'seq2': seq2,
                'i1': i1,
                'i2': i2,
            }
            peptide1['matches'].append(match)

f = open('KL-A1_pep15.out.csv', 'w')
writer = csv.writer(f)
in_keys = ['Sequence','Protein','Protein Description','Length']
headers = in_keys + ['Start', 'End', 'Overlap', '3AA', '4AA', '5AA']
writer.writerow(headers)
first_keys = in_keys + ['Start', 'End']
for peptide in peptides:
    row = [peptide[key] for key in first_keys]
    overlap = False
    aa3 = False
    aa4 = False
    aa5 = False
    for match in peptide['matches']:
        if match['overlap']:
            overlap = True
        else:
            gap = match['gap']
            if gap >= 3:
                aa3 = True
            if gap >= 4:
                aa4 = True
            if gap >= 5:
                aa5 = True
    for is_test in [overlap, aa3, aa4, aa5]:
      row.append('X' if is_test else '')
    writer.writerow(row)
f.close()


