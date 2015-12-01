

# collate_peptides.py 

collate_peptides.py is a script that collates peptide occurences across mulitple
experiments stored as A/B ProteinPilot peptide summaries.

collate_peptides.py creates a row for each peptide sequence found across
the files. A column is added for each source file, and an [X] is
marked if the peptide is found in that source file. The best match
across the source files is chosen, and it's properties are given in
the rest of the column.

## Usage - GUI

To use the collate.py in GUI mode:

    - in Windows, double-click on `collate_peptides.bat`
    - in Mac, double-click on `collate_peptides.command`
    - on the command-line, `python collate_peptides.py -i`

## Usage - command-line

1. copy `example.params` to a new file, say `job.params`

    For reference, it contains:

    {
        "output_dir": "example",
        "collate_file": "collate.csv",
        "reject_file": 'reject.csv',
        "delta_mass_cutoff": 0.05,
        "datasets": [
            ["example/DPrep1.txt", 97], 
            ["example/DPrep2.txt", 96.3],
            ["example/DPrep3.txt", 95.6],
            ["example/Pilot5_example.txt", 96],
        ],
        'contaminants': [
            ['example/B57-w632-contam.txt', 90],
        ]
    }

2. edit `job.params` for your directory and files
3. on the command-line in a terminal:
      
      python collate_peptides.py today.params

4. Look at the results in "collate.csv"



(C) 2015, Bosco K. Ho and Sri Ramararathinam

--

## Archive of results

## March 10, 2015 COLLATE UPDATE

Here are a few things that could be added/modified:
1) The delta mass cut off of 0.05 is applied in current form to reject peptides. can we also make the script reject the peptides that have negative -0.05 or lower. 
         reject peptides that are above 0.05 or below -0.05  (so +/- 0.05)
2) at the moment the collate only keeps the answers we need- i there any chance we can also include a separate file that we can call as “rejected" that has rest of the peptides that failed the criteria (confidence, deltamass)?
        i think at the moment script cleans and removes redundancy at individual file level and uses these cleaned versions for collating together. 
        - can we just take all the rejected peptides and compare them to put in a file called “rejected”?
            - the end point will be two files- "collate" and "rejected” 


3) Also can we add a column that shows the source of the sample from which the best peptide was sourced. that way each peptide will have the sample name attached to it. 

4) Pilot 5 has following columns would be great to add them in 
Pilot 5 Column headings to use for the BEST PEPTIDE in collate.txt and all peptides in rejected:
Conf
Sequence
Modifications
ProteinModifications
Cleavages
dMass
Obs MW
Obs m/z
Theor MW
Theor m/z
Theor z
Sc
Spectrum
Acq Time
Intensity (Peptide)
PrecursorIntensityAcquisition
Apex Time (Peptide)
Elution Peak Width (Peptide)
MS2Counts
Best Sample (Includes name of the file that the best peptide is derived from)

Thanks a a lot Bosco! Please let me know if you need more info or if you need any clarification. I have put example files in folder Collate_update

----

##  Jan2014FinalDataset (31/1/2014)

[done] 1) to create collate.txt like last time (see details below)


2) next step is to perform following analysis on three files: DQ_motif, DR_motif, DP_motif
Peptide sequences are in first column, corresponding protein names are in last column.
    folder DQ7sanitized (18/02/2014)
   
    1) correlation between number of peptides vs abundance of protein
    2) any correlation between number of peptides and protein size
    2) where do most peptides come from- N term, c-term or internal
    3) Is there any correlation between peptides presented and their location on 3D structure?
        we did this last time with a smaller number of peptides where you gave the following info- surface accessibility, helix, loop or sheet location etc .
    4) motif analysis
      you gave me the script for this - I can do this myself.


Files for processing and their relevant cut off values are in ConfCutoffValues.txt


## HOW TO GENERATE COLLATE.TXT

1. Compare all the files and put them in collate.txt
  - Clean up of peptide summaries 
    - Filter each file
      - Apply confidence cut offs for each file (see tab delimited file: ConfCutoffValues.txt)
    - Remove redundant peptides 
      - Keep peptide with best delta mass (column dmass)
      - Keep peptides with modification (ie. Same sequence but post translationally modified— determined by column: ‘Modifications’) - tha way modified and unmodified peptides with same sequence are counted as separate peptides.
  - Apply dmass cut off for each file - reject anything greater than 0.05 in dmass column (smaller than -0.05 or greater than +0.05)   [ONLY CHANGE from last version of your script]
2. Output in a sheet comparing peptides from all samples
  - Choose peptide with best dmass and show all relevant details (all the columns like retention time etc)

-- 

STAGE TWO:


1) correlation between number of peptides vs abundance of protein
2) any correlation between number of peptides and protein size
2) where do most peptides come from- N term, c-term or internal
3) Is there any correlation between peptides presented and their location on 3D structure?
    we did this last time with a smaller number of peptides where you gave the following info- surface accessibility, helix, loop or sheet location etc .
4) motif analysis
    you gave me the script for this - I can do this myself.


# miscellaneous

- apply 0.05 dmass cut off along with the confidence cut off while cleaning data.
- can we have peptide confidence column as well?
- perhaps if possible add one more column with length of peptide sequence - number of amino acids (can be easily done in excel by using ‘Len’ formula, but nice to have it done automatically )


