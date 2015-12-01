import os
import sys
import re

fastas = """

uniprot_sprot.fasta  uniprot_sprot_clean.fasta
uniprot_trembl.fasta  uniprot_trembl_clean.fasta


"""


"""
Do some name mangling for uniprot sequence ID's.
Makes the pure Uniprot ID the first word, or seqID of
the fasta file, and some escaping for HTML display
of the description.
"""


def transform(line):
  if not line.startswith('>'):
    return line
  if '|' not in line.split()[0]:
    return line
  line = re.sub('&', '&amp;', line)
  line = re.sub('\'', '&apos;', line)
  line = re.sub('\"', '&quot;', line)
  line = re.sub('<', '[', line)
  line = re.sub('>', ']', line)
  if line.startswith(']'):
    pieces = line[1:].split('|')
    line = '>{1} {0}|{1}|{2}'.format(*pieces)
  return line


def convert_fasta(fasta, out_fasta):
  if not os.path.isfile(fasta):
    print "Couldn't find {0}.".format(fasta)
    return
  print "Processing {0}".format(fasta)
  out_f = open(out_fasta, 'w')
  with open(fasta, 'Ur') as in_f:
    for i, line in enumerate(in_f):
      out_f.write(transform(line))
      if i % 100000 == 0 and i > 0:
        print 'Processed {0} lines'.format(i)


for line in fastas.splitlines():
  if line.strip():
    fasta, out_fasta = line.split()[:2]
    convert_fasta(fasta, out_fasta)



original_bash_script = """
# add escaping of special chars ', ", & and substitute '<' to '[' and '>' to ']' in def lines
# warning, substitution is expensive aka slow; comment out next four commands
# if you can't afford the time and deal with the special chars yourself
#
sed -e '/^>/s/&/\&amp;/g' \
    -e "/^>/s/'/\&apos;/g" \
    -e '/^>/s/"/\&quot;/g' \
    -e '/^>/s/</[/g' \
    -e '/^>/s/>/]/g' \
    -e 's/^]\(..\)|\(.*\)|\(.*\)/>\2 \1|\2|\3/' < uniprot_sprot.fasta > uniprot_sprot.fasta.tmp
mv -f uniprot_sprot.fasta.tmp uniprot_sprot.fasta

sed -e '/^>/s/&/\&amp;/g' \
    -e "/^>/s/'/\&apos;/g" \
    -e '/^>/s/"/\&quot;/g' \
    -e '/^>/s/</[/g' \
    -e '/^>/s/>/]/g' \
    -e 's/^]\(..\)|\(.*\)|\(.*\)/>\2 \1|\2|\3/' < uniprot_trembl.fasta > uniprot_trembl.fasta.tmp
mv -f uniprot_trembl.fasta.tmp uniprot_trembl.fasta
"""
