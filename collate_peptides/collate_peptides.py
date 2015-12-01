#!/usr/bin/env python


__version__ = "2"

# now handles Pilot5 .csv files as well as Pilot4 .txt files

__doc__ = """
collate.py v%s - collates and cleans-up AB/X ProteinPilot(4|5) files.

1. Removes matches with low probability and poor masses 
2. Finds a single best match for a peptide 
3. Collates peptide identifications across several files

(C) 2015, Bosco K. Ho and Sri Ramararathinam

Usage:
   python collate.py -i
   python collate.py param

param is a file with the following format:

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
""" % __version__


import sys
import os
import json
import glob
import csv
import subprocess
from collections import defaultdict


import tkform


def get_protein_pilot_headers():
    header4_str = """
    Sequence
    Modifications
    ProteinModifications
    dMass
    Conf
    Cleavages
    Prec MW
    Prec m/z
    Theor MW
    Theor m/z
    Theor z
    Sc
    Spectrum
    Time
    PrecursorSignal
    PrecursorElution
    Names
    Accessions
    """
    headers4 = [k.strip() for k in header4_str.splitlines()]

    header5_str = """
    Sequence
    Modifications
    ProteinModifications
    dMass
    Conf
    Cleavages
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
    Names
    Accessions
    """
    headers5 = [k.strip() for k in header5_str.splitlines()]

    headers = headers4
    for header in headers5:
        if header not in headers:
            headers.append(header)

    return headers


def logging(s, callback=None):
    sys.stdout.write(s)


def get_base(fname):
    return os.path.splitext(os.path.basename(fname))[0]


def get_ext(fname):
    return os.path.splitext(fname)[-1]


def guess_delimiter(fname):
    with open(fname, 'Ur') as f:
        line = f.readline()
        if '\t' in line:
            return '\t'
        elif ',' in line:
            return ','
        else:
            raise IOError("Can't recognize delimiter in first line of '%s'" % ext)


def get_headers(fname):
    delimiter = guess_delimiter(fname)
    with open(fname, 'Ur') as f:
        reader = csv.reader(f, delimiter=delimiter)
        return reader.next()


def read_csv(fname):
    delimiter = guess_delimiter(fname)
    with open(fname, 'Ur') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_csv(fname, rows):
    ext = get_ext(fname).lower()
    if ext == '.csv':
        delimiter = ','
    elif ext in ['.txt', '.tsv']:
        delimiter = '\t'
    else:
        raise IOError("Can't recognize fname extension '%s'" % ext)
    with open(fname, 'w') as f:
        writer = csv.writer(f, delimiter=delimiter)
        for row in rows:
            writer.writerow(row)


def open_file(filepath):
    if sys.platform.startswith('darwin'):
        subprocess.call(('open', filepath))
    elif os.name == 'nt':
        os.startfile(filepath)
    elif os.name == 'posix':
        subprocess.call(('xdg-open', filepath))


def mangle_fname(top_dir, fname, pre_ext):
    base, ext = os.path.splitext(os.path.basename(fname))
    return os.path.join(top_dir, base + pre_ext + ext)


def make_peptide_filter(confidence_cutoff, dmass_cutoff, contaminants):
    def is_good(peptide):
        sequence = peptide['Sequence'] + peptide['Modifications']
        if sequence in contaminants:
            return False
        if float(peptide['Conf']) < confidence_cutoff:
            return False
        if abs(float(peptide['dMass'])) > dmass_cutoff:
            return False
        return True
    return is_good


def write_clean_peptide_summary(
        fname, output_dir, peptide_filter):

    logging("Reading ProteinPilot '{}'\n".format(fname))
    headers = get_headers(fname)

    clean_peptides = []
    for peptide in read_csv(fname):
        if peptide_filter(peptide):
            clean_peptides.append(peptide)

    # filter out low confidence peptides and delta_mass are too big
    out_fname = mangle_fname(output_dir, fname, '.clean')
    logging("Writing '{}'\n".format(out_fname))
    rows = [headers]
    for peptide in clean_peptides:
        rows.append([peptide[h] for h in headers])
    write_csv(out_fname, rows)

    # only show one psm for each peptide
    unique_peptides = []
    peptides_by_seq = defaultdict(list)
    for peptide in clean_peptides:
        sequence = peptide['Sequence'] + peptide['Modifications']
        peptides_by_seq[sequence].append(peptide)
    for sequence, peptides in peptides_by_seq.items():
        peptides.sort(key=lambda p: abs(float(p['dMass'])))
        unique_peptides.append(peptides[0])
    unique_peptides.sort(key=lambda p: p['Sequence'])
    unique_peptides.sort(key=lambda p: -float(p['Conf']))

    out_fname = mangle_fname(output_dir, fname, '.unique')
    logging("Writing '{}'\n".format(out_fname))
    rows = [headers]
    for peptide in unique_peptides:
        rows.append([peptide[h] for h in headers])
    write_csv(out_fname, rows)


def get_grouped_peptides_by_seq(
        fnames, experiments, confidence_cutoffs, dmass_cutoff, contaminant_seqs):

    grouped_peptides_by_seq = {}
    rejected_peptides_by_seq = {}

    for experiment, in_fname, confidence_cutoff in zip(experiments, fnames, confidence_cutoffs):

        peptide_filter = make_peptide_filter(
            confidence_cutoff, dmass_cutoff, contaminant_seqs)

        for peptide in read_csv(in_fname):
            peptide['Source'] = experiment
            sequence = peptide['Sequence'] + peptide['Modifications']
            if peptide_filter(peptide):
                if sequence not in grouped_peptides_by_seq:
                    grouped_peptides_by_seq[sequence] = { 'experiments': {} }
                grouped_peptides_by_seq[sequence]['experiments'][experiment] = peptide
            else:
                if sequence not in rejected_peptides_by_seq:
                    rejected_peptides_by_seq[sequence] = { 'experiments': {} }
                rejected_peptides_by_seq[sequence]['experiments'][experiment] = peptide

    for group in grouped_peptides_by_seq.values():
        peptides = group['experiments'].values()
        peptides.sort(key=lambda e: abs(float(e['dMass'])))
        group['best_peptide'] = peptides[0]

    for group in rejected_peptides_by_seq.values():
        peptides = group['experiments'].values()
        peptides.sort(key=lambda e: abs(float(e['dMass'])))
        group['best_peptide'] = peptides[0]

    return grouped_peptides_by_seq, rejected_peptides_by_seq


def write_group_to_csv(
        grouped_peptides_by_seq, out_csv, experiments, print_headers, contaminant_seqs):

    rows = []

    found_headers = []
    for grouped_peptides in grouped_peptides_by_seq.values():
        for header in grouped_peptides['best_peptide']:
            if header in print_headers and header not in found_headers:
                found_headers.append(header)
    print_headers = [k for k in print_headers if k in found_headers]

    headers = []
    headers += experiments
    headers.append('N')
    headers.append('Contaminant')
    headers.append('Best Sample')
    headers += print_headers

    rows.append(headers)

    for seq in sorted(grouped_peptides_by_seq.keys()):
        row = []
        grouped_peptides = grouped_peptides_by_seq[seq]
        peptide = grouped_peptides['best_peptide']
        n = 0
        for experiment in experiments:
            if experiment in grouped_peptides['experiments']:
                n += 1
                row.append('X')
            else:
                row.append('')
        row.append(n)
        if seq in contaminant_seqs:
            row.append('X')
        else:
            row.append('')
        row.append(peptide['Source'])
        row += [peptide[header] if header in peptide else '' for header in print_headers]
        rows.append(row)

    i_move = headers.index('Modifications')
    for row in rows:
        sequence = row[i_move]
        del row[i_move]
        row.insert(0, sequence)

    i_move = headers.index('Sequence')
    for row in rows:
        sequence = row[i_move]
        del row[i_move]
        row.insert(0, sequence)

    write_csv(out_csv, rows)


def check_csv(fname):
    if not os.path.isfile(fname):
        raise IOError("Couldn't find " + fname)
    check_headers = ['Sequence', 'Modifications', 'Conf', 'dMass']
    headers = get_headers(fname)
    for check_header in check_headers:
        if check_header not in headers:
            raise IOError("Couldn't find '%s' column in %s" % (check_header, fname))


def collate_peptide_summaries(params):

    dmass_cutoff = float(params['delta_mass_cutoff'])

    output_dir = os.path.abspath(params['output_dir'])
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if len(params['datasets']) == 0:
        raise IOError("No files selected")

    # read contaminats
    logging(">> Loading contaminants:\n")
    contaminant_seqs = set()
    for dataset in params['contaminants']:
        fname = dataset[0]
        check_csv(fname)
        confidence_cutoff = 0.0
        peptide_filter = make_peptide_filter(
            confidence_cutoff, dmass_cutoff, [])
        logging("Reading ProteinPilot '{}'\n".format(fname))
        for peptide in read_csv(fname):
            if peptide_filter(peptide):
                seq =  peptide['Sequence'] + peptide['Modifications']
                contaminant_seqs.add(seq)

    # read experiments and filter them
    logging(">> Loading individual peptide summaries:\n")
    experiments = []
    fnames = []
    confidence_cutoffs = []
    for dataset in params['datasets']:
        fname = dataset[0]
        check_csv(fname)
        confidence_cutoff = float(dataset[1])
        peptide_filter = make_peptide_filter(
            confidence_cutoff, dmass_cutoff, contaminant_seqs)

        write_clean_peptide_summary(
            fname, output_dir, peptide_filter)

        experiments.append(get_base(fname))
        fnames.append(fname)
        confidence_cutoffs.append(confidence_cutoff)

    logging('>> Collating peptides ...\n')
    grouped_peptides_by_seq, rejected_peptides_by_seq = get_grouped_peptides_by_seq(
        fnames, experiments, confidence_cutoffs, dmass_cutoff, contaminant_seqs)

    reject_file = os.path.join(output_dir, params['reject_file'])
    logging("Writing rejected peptides to '{}'\n".format(reject_file))
    write_group_to_csv(rejected_peptides_by_seq, reject_file, experiments, get_protein_pilot_headers(), contaminant_seqs)

    collate_file = os.path.join(output_dir, params['collate_file'])
    logging("Writing collated peptides to ")
    logging("'{}'\n".format(collate_file), lambda: open_file(collate_file))
    write_group_to_csv(grouped_peptides_by_seq, collate_file, experiments, get_protein_pilot_headers(), contaminant_seqs)

    logging(">> Done\n")


class CollateForm(tkform.Form):

    def __init__(self, width=700, height=800):
        tkform.Form.__init__(self, 'Collate Peptides across Multiple Experiments', width, height)

        self.push_text("Collate Peptides across Multiple Experiments", 30)

        self.push_line(width=600)
        self.push_text("(C) 2015, Bosco K. Ho and Sri Ramararathinam")
        self.push_spacer()

        self.push_text("Load ProteinPilot peptide summary files (txt|csv)")
        self.push_text(u"Reorder the files by dragging \u2630 up and down")
        self.push_text("Reject peptides with the Confidence lower-cutoff")
        self.push_summary_list('datasets', '+ peptide summaries (txt|csv)')
        self.push_spacer()
        self.push_text("Load contaminant peptides")
        self.push_summary_list('contaminants', '+ peptide summaries (txt|csv)')
        self.push_spacer()
        self.push_labeled_param(
            'delta_mass_cutoff', 
            'Delta mass cutoff', 
            entry_text='0.05',
            width=4)
        self.push_labeled_param(
            'output_dir', 
            'Output directory', 
            os.getcwd(), 
            load_dir_text='select',
            width=50)
        self.push_labeled_param(
            'collate_file', 
            'Collated peptides file', 
            entry_text='collate.csv')
        self.push_labeled_param(
            'reject_file', 
            'Rejected peptides file', 
            entry_text='reject.csv')

        self.push_spacer()
        self.push_line(width=600)

        self.push_text("Output", 16)
        self.push_submit()
        self.push_output()

    def push_summary_list(
            self, param_id, load_file_text, is_label=True):
        file_list = tkform.ReorderableList(self.interior)

        def load_file():
            fnames = tkform.askopenfilenames(title=load_file_text)
            for fname in fnames:
                file_list.add_entry_label(fname, '95', 4)

        load_files_button = tkform.tk.Button(
            self.interior, text=load_file_text, command=load_file)
        self.push_row(load_files_button)

        self.push_row(file_list)
        self.mouse_widgets.append(file_list)
        self.param_entries[param_id] = file_list

    def run(self, params):
        self.clear_output()
        global logging
        logging = self.print_output
        collate_peptide_summaries(params)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print __doc__

    elif sys.argv[1] == '-i':
        # run the GUI
        CollateForm(800, -50).mainloop()

    elif sys.argv[1] == '-t':
        # testing for debugging purposes
        params = \
        {
            "output_dir": "example",
            "collate_file": "collate.csv",
            "reject_file": 'reject.csv',
            "delta_mass_cutoff": 0.05,
            "datasets": [
                ["example/DPrep1.txt", 97], 
                ["example/DPrep2.txt", 96.3],
                ["example/DPrep3.txt", 95.6],
                ["example/Pilot5_example.csv", 96],
            ],
            'contaminants': [
                ['example/B57-w632-contam.txt', 90],
            ]
        }
        collate_peptide_summaries(params)

    else:
        # run in command-line mode with a params file
        params_fname = sys.argv[1]
        logging('Loading ' + params_fname)
        params = eval(open(params_fname).read())
        collate_peptide_summaries(params)




