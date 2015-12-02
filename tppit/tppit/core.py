#
# (C) 2013 Bosco Ho
#
# A Python interface to the command-line tools of the TPP
# http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP


import sys
import os
import glob

import envoy
from jinja2 import Environment, FileSystemLoader


this_dir = os.path.dirname(__file__)
template_dir = os.path.join(this_dir, 'templates')
other_binary_dirs = [] # a common linux location


def check_files(*fnames):
  """
  Checks for existence of fnames. Raises error if not found.
  """
  for fname in fnames:
    if not os.path.isfile(fname):
      raise IOError("Error: can't find {}".format(fname))


def which(program):
  """
  Reproduces Unix 'which' and looks in other_binary_dirs
  """
  def is_binary(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
  fpath, fname = os.path.split(program)
  if fpath:
    if is_binary(program):
      return program
  else:
    binary_dirs = os.environ["PATH"].split(os.pathsep)
    binary_dirs.extend(other_binary_dirs)
    for path in binary_dirs:
      exe_file = os.path.join(path, program)
      if is_binary(exe_file):
        return exe_file
  return None


def change_ext(fname, new_ext):
  """
  Swaps extension to new_ext
  """
  if not new_ext:
    return fname
  base, ext = os.path.splitext(fname)
  if not new_ext.startswith('.'):
    new_ext = '.' + new_ext
  return base + new_ext


def render_jinja2(template_fname, env, out_fname):
  """
  Renders a file using a .jinja2 template with variables 
  piped in by env.
  """
  dirname, basename = os.path.split(template_fname)
  jinja2_env = Environment(loader=FileSystemLoader(dirname))
  template = jinja2_env.get_template(basename)
  open(out_fname, 'w').write(template.render(**env))


def quote(s):
  return '\'' + str(s) + '\''


def check_dir(d):
  if not os.path.isdir(d):
    os.makedirs(d)


def run(cmd, log_file=None):
  """
  Executes cmd with comprehensive directory search for
  location of binary (first word). In Windows, also 
  searches .exe extension
  """
  bin = cmd.split()[0]
  full_bin = bin
  if not which(full_bin):
    full_bin += '.exe'
    if not which(full_bin):
      raise ValueError("Couldn't find " + full_bin)
  cmd = cmd.replace(bin, full_bin)
  print cmd
  r = envoy.run(cmd)
  if log_file is not None:
    check_dir(os.path.dirname(log_file))
    open(log_file, 'w').write(r.std_out + '\n' + r.std_err)


def mascot2pepxml(mascot_dat, fasta_db, pepxml):
  check_files(mascot_dat, fasta_db)
  run('Mascot2XML ' + quote(mascot_dat) + ' -D' + quote(fasta_db) + ' -notgz -nodta', pepxml + '.log')
  default_pepxml = change_ext(mascot_dat, '.pep.xml')
  if default_pepxml != pepxml:
    os.rename(default_pepxml, pepxml)


def mzxml2mzml(mzxml):
  check_files(mzxml)
  output_dir = os.path.dirname(mzxml)
  run('msconvert ' + mzxml + ' -o "{}"'.format(output_dir), change_ext(mzxml, '.mzML') + '.log')


def tandem2pepxml(tandem, pepxml):
  check_files(tandem)
  run('Tandem2XML ' + quote(tandem) + ' ' + quote(pepxml))


def merge_pepxmls(pepxmls, merged_pepxml, min_pep_len=7):
  check_files(*pepxmls)
  args = ' ' + quote(merged_pepxml) + ' ' + ' '.join(map(quote, pepxmls))
  if min_pep_len is not None:
    args += ' -L' + str(min_pep_len)
  run('InteractParser ' + args, merged_pepxml + '.log')


def peptide_prophet(pepxml, prob_cutoff=0.05, decoy_prefix=None):
  """
  Calculates probabilities of peptide spectrum matches
  by fitting a two peak distribution to the frequency
  against the match scores.
  """
  check_files(pepxml)
  args = ' ' + quote(pepxml)
  if prob_cutoff is not None:
    args += ' MINPROB={}'.format(prob_cutoff)
  if decoy_prefix is not None:
    args += ' DECOY=' + decoy_prefix
  run('PeptideProphetParser ' + args, pepxml + '.peptideprophet.log')


def protein_prophet(pepxml, protxml, prob_cutoff=0.05):
  """
  Calculates probabilities of protein identification
  based on how many unique peptides are found in
  each protein group, or proteins indistinguishable
  from the peptide match.
  """
  check_files(pepxml)
  args = ' ' + quote(pepxml) + ' ' + quote(protxml)
  if prob_cutoff is not None:
    args += ' NOPLOT MINPROB={}'.format(prob_cutoff)
  run('ProteinProphet ' + args, protxml + '.log')


def resolve_proteins(pepxmls, protxml, prob_cutoff=None, decoy_prefix=None):
  """
  Uses the TPP stack to identify proteins given
  a set of peptide matches based on previous database
  searches.
  """
  interact_pepxml = protxml.replace('.prot.xml', '.pep.xml')
  merge_pepxmls(pepxmls, interact_pepxml)
  peptide_prophet(interact_pepxml, prob_cutoff=prob_cutoff, decoy_prefix=decoy_prefix)
  protein_prophet(interact_pepxml, protxml, prob_cutoff)


def xtandem_match_spectra(mzxml, fasta_db):
  check_files(mzxml)
  check_files(fasta_db)
  mzxml = os.path.abspath(mzxml)
  fasta_db = os.path.abspath(fasta_db)
  tandem = change_ext(mzxml, '.tandem')
  taxonomy_xml = change_ext(mzxml, '.taxonomy.xml')
  tandem_params = change_ext(mzxml, '.tandem.params')
  env = {
    'taxonomy_xml': taxonomy_xml,
    'database_id': 'mydatabase',
    'database': fasta_db,
    'input': mzxml,
    'output': tandem,
    'mass_error': 0.8
  }
  template = os.path.join(template_dir, 'taxonomy.xml.jinja2')
  render_jinja2(template, env, taxonomy_xml)
  template = os.path.join(template_dir, 'tandem.params.jinja2')
  render_jinja2(template, env, tandem_params)
  run('tandem ' + quote(tandem_params), tandem + '.log')


def make_decoy_fasta(fasta_db, decoy_db=None, decoy_prefix='DECOY_'):
  if decoy_db is None:
    decoy_db = fasta_db.replace('.fasta', '_decoy.fasta')
  run('decoyFasta ' + ' -t ' + decoy_prefix + ' ' + quote(fasta_db) + ' ' + quote(decoy_db))

 




