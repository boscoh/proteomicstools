
# tppit

`tppit` is a Python module that wraps the Transatlantic Proteomics Pipeline (TPP).

The TPP combines a number of tools that convert shotgun massspec experimental data into peptides and protein identifications, and in some cases, even quantification. 

However, these tools have inconsistent interfaces and requires a bit of glue to tie them together. `tppit` wraps these command-line tools so that a consistent interface can be run in a automated workflow.

## Typical Protein Identification workflow

1. `msconvert` converts raw mass-spec data into the standard `mzXML`/`mzML` formats.

2. `tandem` is one (of many) search engines that matches massspec spectra to a database of candidate peptides (typically given as a fasta file of protein sequences with digestion enzyme).

3. `PeptideProphetParser` takes a set of spectrum-peptide matches with given expectations (in `.pep.xml` format) and assigns probabilities to the peptides by assuming a bi-modal frequency distribution of  expectations, corresponding to incorrect/correct matches.

4. `ProteinProphetParse` attempts to find a minimal grouping of proteins using parsimony based on the number of unique peptides in each proteins and the probabilities of each peptide.

## Why Python?

Doing this analysis means you will be processing a lot of files, and all of them will be interrelated in a lot of fiddly ways. Python is great at handling these details.

First load the module:

    import tppit

Now `tppit` will search for the TPP programs in the system path, but you can always specify it directly:

    tppit.other_binary_dirs.append('/usr/local/tpp/bin')

## Matching spectra to peptides

Let's say you have some MS/MS spectra in the form of `mzXML` files in the `example` directory: 

    mzxmls = ['example/xtandem/Seq2328{}_E1O1.mzXML'.format(i) for i in [2,3,4]]

You have some candidate proteins in: 

    fasta = 'example/xtandem/orf_trans.fasta'

Then run the search-engine provided by the TPP `x!tandem` with:

    for mzxml in mzxmls:
      tppit.xtandem_match_spectra(mzxml, fasta)

The resultant spectrum-peptide matches are saved in:

    tandems = [tppit.change_ext(m, '.tandem') for m in mzxmls]

## Converting to .pep.xml for further processing

Okay so you have a set of spectrum-peptide matches. For more TPP goodness, they need to be converted to `.pep.xml` format.

Let's say they are in the form of `.tandem` files as above. We need to convert them to `.pep.xml` files:

    for tandem, pepxml in zip(tandems, pepxmls):
      tppit.tandem2pepxml(tandem, pepxml)
    pepxmls = [t + '.pep.xml' for t in tandems]

As another example, there are some mascot `.dat` files in the examples:

    mascots = ['example/mascot/F02204{}.dat'.format(i) for i in [3, 5]]

Now the conversion of mascot files to `.pep.xml` will require the sequence database, so:

    fasta_db = '../db/HUMAN.fasta'

Now run the translator:

    for mascot, pepxml in zip(mascots, pepxmls):
      tppit.mascot2pepxml(mascot, fasta_db, pepxml)
    pepxmls = [tppit.change_ext(m, '.pep.xml') for m in mascots]

So finally, let's predict the proteins:

    protxml = 'example/xtandem/interact.prot.xml'
    tppit.resolve_proteins(pepxmls, protxml)
  
The function `resolve_proteins` will:

  1. combine all the `.pep.xml` files with `InteractParser`
  2. generate peptie probabilities with `PeptideProphetParser` 
  3. generate parsimony protein lists with `ProteinProphetParser`

  
