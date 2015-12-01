import os
import csv


"""
Uses Nathan Crofts's clustering algorithm to generate a list
of smallest non-redundant peptides. Then aligns the
sequences using ClustalW2 mainly because Clustalw2 
allows a large gap penalty, which is our way of doing
gapless alignment.
"""


def make_reduced_fasta(sequences):
  """
  Given sequences which is a list of strings, will apply
  a set of filters, then cluster together peptides with a 5 residue
  common motif and return a nonhomologous set of shortes peptides, 
  one from each cluster. Designed by Nathan Croft.
  """

  def delete_same_but_longer():
    sequences.sort()
    n_sequence = len(sequences)
    # since deleting from a list, always traverse 
    # from the end first
    for i in reversed(range(1, n_sequence)):
      if sequences[i-1] in sequences[i]:
        del sequences[i]

  def reverse_sequences():
    n_sequence = len(sequences)
    for i in range(n_sequence):
      sequences[i] = sequences[i][::-1]
    sequences.sort()

  def delete_with_subpattern(i_position, subpattern_length):
    n_sequence = len(sequences)
    for i in reversed(range(n_sequence-1)):
      sequence = sequences[i]
      i_seq = i_position
      j_seq = i_position + subpattern_length
      if j_seq > len(sequence):
        continue
      test_seq = sequence[i_seq:j_seq]
      n_sequence_now = len(sequences)
      for j in reversed(range(i+1, n_sequence_now)):
        if test_seq in sequences[j]:
          del sequences[j]

  def delete_too_short(cutoff_length):
    n_sequence = len(sequences)
    for i in reversed(range(n_sequence)):
      if len(sequences[i]) < cutoff_length:
        del sequences[i]

  print "originally %d sequences" % len(sequences)

  delete_too_short(7)

  delete_same_but_longer()
  print "deleted contained sequences to %d" % len(sequences)

  reverse_sequences()
  delete_same_but_longer()
  reverse_sequences()
  print "deleted reverse contained sequences to %d" % len(sequences)

  for i_position in range(10):
    delete_with_subpattern(i_position, 5)
    print "deleted subpatterns starting at %d to %d" % (i_position, len(sequences))

  return sequences


def write_sequences_to_fasta(out_fasta, sequences):
  f = open(out_fasta, 'w')
  for i, seq in enumerate(sequences):
    title = ">seq%d \n" % (i)
    f.write(title)
    f.write(seq + "\n")
  f.close()


def read_csv_file(fname, parse_list):
  """
  A generic function to read a CSV file into a list of 
  properties. The parse_list describes the parsing of the CSV
  file for each row. The parse_list contains a list of 2-pules where
  the first entry gives the name of the property, and the
  second entry gives a function that parses the string in
  the CSV file. If it is to be read as a straight string,
  None should be put there.
  """

  result = []

  reader = csv.reader(open(fname, "Ur"))

  for i, words in enumerate(reader):
    if i == 0:
      continue
    result.append({})
    try:
      for i, (key, fn) in enumerate(parse_list):
        val = words[i]
        if fn is not None:
          val = fn(val)
        result[-1][key] = val
    except:
      print "problem reading line %d" % i
  return result


def read_sequences_from_fasta(fname):
  """
  Reads a bunch of sequences in from a FASTA file,
  and returns it as a list of strings.
  """
  sequences = []
  for l in open(fname):
    if l.startswith(">"):
      sequences.append("")
    else:
      sequences[-1] += l.strip()
  return sequences



def build_alignment_profile(sequences):
  """
  This produces an amino acid profile based on a set
  of sequences that has been aligned using clustalw2.
  The expected sequences includes "-" gaps in order
  from the left/right so that the main motif is in
  the middle. The results are written as a table
  in CSV form that can be imported into Excel.
  """

  AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
  n_aa = len(AMINO_ACIDS)
  n_res = len(sequences[0])

  alignment = []
  counts = []
  for i in range(n_res):
    alignment.append([0 for j in range(n_aa)])
    count = 0
    for sequence in sequences:
      aa = sequence[i]
      if aa not in AMINO_ACIDS:
        if aa != "-":
          print "Can't recognize character %s in %s" % (aa, sequence)
        continue
      i_aa = AMINO_ACIDS.index(aa)
      alignment[i][i_aa] += 1
      count += 1
    counts.append(count)

  for i in range(n_res):
    for i_aa in range(n_aa):
      count_aa = alignment[i][i_aa]
      alignment[i][i_aa] = int(count_aa/float(counts[i])*100)

  rows = []
  for i_aa in range(n_aa):
    aa = AMINO_ACIDS[i_aa]
    rows.append([aa] + [alignment[i][i_aa] for i in range(n_res)])
  rows.append([])
  rows.append(['counts'] + counts)

  return rows


def align_with_clustalw2(in_fasta, out_fasta):
  """
  Externally run with clustalw2 with very high penalties for
  gaps and gap extensions to make it gapless.
  """
  cmd = '"/Users/Ralfold/Monash/Python/clustalw2" '
  cmd += ' -infile=' + in_fasta
  cmd += ' -output=' + out_fasta
  cmd += ' -gapopen=10 -gapext=10 -output=fasta' 
  print cmd
  os.system(cmd)


def compute_motif(in_csv, parse_list, out_name):
  """
  The main execution loop. Basically, this takes a CSV file
  that stores the sequences, and generates 3 output files:
    1. fasta file with all the truncated sequences
    2. after externally running clustalw2, an aligned fasta
    3. a CSV with the alignment profile in percentages for 
       each amino acids.
  """

  peptides = read_csv_file(in_csv, parse_list)
  sequences = [p['sequence'] for p in peptides]

  sequences = make_reduced_fasta(sequences)

  fasta = out_name + '.fasta'
  write_sequences_to_fasta(fasta, sequences)

  clustalw2_fasta = out_name + '.clustalw2.fasta'  
  align_with_clustalw2(fasta, clustalw2_fasta)
  sequences = read_sequences_from_fasta(clustalw2_fasta)

  alignment_rows = build_alignment_profile(sequences)

  writer = csv.writer(open(out_name + '.csv', 'w'))
  for row in alignment_rows:
    writer.writerow(row)



compute_motif(
    'raw_peptides.csv', 
    [('sequence', None), ('gene', None), ('description', None),],
    'sorted_and_aligned_peptides')






