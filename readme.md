
# Proteomics Tools

This is a repository of the tools written at the Monash Proteomics Facility

- `best_seqid`: searches through a list of protein seqid's and looks up uniprot to find the longest and best-reviewed seqid's.

- `cluster_peptides`: a tool to cluster overlapping peptides together, and identify a kernel of shortest motifs

- `cluster_protein_sets`: for a given set of peptide-protein matches, groups together peptides that are found in the same protein

- `coverage_overlay`: this takes a series of peptides in the ProteinPilot 4 summary form and combines them together into one large file where a table shows in which file that a peptide was detected.

- `down_pilot`: a simple conversion tool for ProteinPilot 5 down to ProteinPilot 4

- `format_uniprot`: for the TPP package, fasta files need to be preprocessed to change the name entries. Although a perl script is provided, this does it in python.