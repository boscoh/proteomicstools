#!/usr/bin/env python


__version__ = "2"

# now handles Pilot5 .csv files as well as Pilot4 .txt files

__doc__ = """
down_pilot %s

(C) 2015, Bosco K. Ho and Sri Ramararathinam

Usage:
   python down_pilot.py -i
   python down_pilot.py csv...

""" % __version__


import sys
import os
import json
import glob
import csv
import subprocess
from collections import defaultdict
from pprint import pprint

import tkform



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


transform = {
    'Acq Time': 'Time',
    'Apex Time (Peptide)': 'PrecursorElution',
    'Elution Peak Width (Peptide)': None,
    'Intensity (Peptide)': None,
    'MS2Counts': None,
    'Obs MW': 'Prec MW',
    'Obs m/z': 'Prec m/z',
    'PrecursorIntensityAcquisition': 'PrecursorSignal',
    'ProteinModifications': None
}


pilot4_headers = get_headers('example/DPrep1.txt')
pilot5_headers = get_headers('example/Pilot5_example.csv')



def convert(csv_fnames, logging=logging):
    this_dir = os.getcwd()
    for csv_fname in csv_fnames:
        csv_fname = os.path.relpath(csv_fname, this_dir)
        headers = get_headers(csv_fname)
        if headers == pilot4_headers:
            logging('Already pilot4, skipping ' + csv_fname + '\n')
        else:
            rows = []
            rows.append(pilot4_headers)
            for entry in read_csv(csv_fname):
                new_entry = {}
                for key, val in entry.items():
                    if key in transform:
                        new_key = transform[key]
                        if new_key is None:
                            continue
                        else:
                            new_entry[new_key] = val
                    else:
                        new_entry[key] = val

                row = []
                for key in pilot4_headers:
                    if key in new_entry:
                        row.append(new_entry[key])
                    elif key == "Specific":
                        row.append(1)
                rows.append(row)

            name, ext = os.path.splitext(csv_fname)
            new_csv_fname = name + '.pilot4.txt'
            write_csv(new_csv_fname, rows)
            logging('Made %s\n' % new_csv_fname)


class ConvertForm(tkform.Form):

    def __init__(self, width=700, height=800):
        tkform.Form.__init__(
            self, 
            'Convert ProteinPilot5 to ProteinPilot4', 
            width, 
            height)

        self.push_text(
            "Convert ProteinPilot5 to ProteinPilot4", 30)

        self.push_line(width=600)
        self.push_spacer()

        self.push_text(
            "Load ProteinPilot peptide summary files (txt|csv)")
        self.push_file_list_param(
            'datasets', 
            '+ peptide summaries (txt|csv)', 
            is_label=False)

        self.push_spacer()
        self.push_line(width=600)

        self.push_text("Output", 16)
        self.push_submit()
        self.push_output()

    def run(self, params):
        self.clear_output()
        convert([p[0] for p in params['datasets']], self.print_output)


def make_transform():

    transform_str = """\
    %Cov
    %Cov(50)
    %Cov(95)
    Accessions
    Acq Time [Time]
    Annotation
    Apex Time (Peptide) [PrecursorElution]
    Cleavages
    Conf
    Contrib
    Elution Peak Width (Peptide) []
    Intensity (Peptide) []
    MS2Counts []
    Modifications
    N
    Names
    Obs MW [Prec MW]
    Obs m/z [Prec m/z]
    PrecursorIntensityAcquisition [PrecursorSignal]
    ProteinModifications []
    Sc
    Sequence
    Spectrum
    Theor MW
    Theor m/z
    Theor z
    Total
    Unused
    Used
    dMass"""
    transform = {}
    for line in transform_str.splitlines():
        words = line.split()
        key_words = []
        new_words = []
        for word in words:
            if '[' not in word and ']' not in word:
                key_words.append(word)
            else:
                new_words.append(word.replace('[', '').replace(']', ''))
        key_word = ' '.join(key_words)
        new_word = ' '.join(new_words)
        if new_words:
            if new_word:
                transform[key_word] = new_word
            else:
                transform[key_word] = None

    return transform


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print __doc__
    elif sys.argv[1] == '-i':
        # run the GUI
        ConvertForm(800, -50).mainloop()
    elif sys.argv[1] == '-t':
        # testing for debugging purposes
        convert(
            ['example/DPrep1.txt', 
            'example/Pilot5_example.csv'])
    else:
        # run in command-line mode 
        convert(sys.argv[1:])



