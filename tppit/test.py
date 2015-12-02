import glob
import tppit
import os



def test_mascot_to_protxml():
  tppit.other_binary_dirs.append('/usr/local/tpp/bin')
  mascots = ['example/mascot/F02204{}.dat'.format(i) for i in [3, 5]]
  fasta_db = 'example/mascot/HUMAN.fasta'
  pepxmls = [tppit.change_ext(m, '.pep.xml') for m in mascots]
  for mascot, pepxml in zip(mascots, pepxmls):
    tppit.mascot2pepxml(mascot, fasta_db, pepxml)
  tppit.resolve_proteins(pepxmls, 'example/mascot/interact.prot.xml')


def test_mzxml_to_protxml():
  tppit.other_binary_dirs.append('/usr/local/tpp/bin')
  fasta = 'example/xtandem/orf_trans.fasta'
  decoy_fasta = tppit.change_ext(fasta, '.decoy.fasta')
  tppit.make_decoy_fasta(fasta, decoy_fasta, 'DECOY_')
  mzxmls = ['example/xtandem/Seq2328{}_E1O1.mzXML'.format(i) for i in [2,3,4]]
  for mzxml in mzxmls:
    tppit.xtandem_match_spectra(mzxml, decoy_fasta)
  tandems = [tppit.change_ext(m, '.tandem') for m in mzxmls]
  pepxmls = [t + '.pep.xml' for t in tandems]
  for tandem, pepxml in zip(tandems, pepxmls):
    tppit.tandem2pepxml(tandem, pepxml)
  tppit.resolve_proteins(pepxmls, 'example/xtandem/interact.prot.xml', decoy_prefix='DECOY_')
  



if __name__ == "__main__":
  test_mascot_to_protxml()
  test_mzxml_to_protxml()
