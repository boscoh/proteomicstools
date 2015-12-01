#!/usr/bin/env python


__doc__ = """

best_seqid.py

Reads a csv|tsv|txt file and looks for uniprot accessions

Then fetches metadata from http://uniprot.org and chooses
the best uniprot accession based on:

1. whether the entry has been manually curated
2. longest sequence

Usage:
    
    python best_seqid.py -i 
    python best_seqid.py 

""" 


import sys
import os
import json
import glob
import csv
import subprocess
from collections import defaultdict

import requests
import uniprot
import tkform


def logging(s, dummy_callback=None):
    sys.stdout.write(s)


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


def open_file(fname):
    if sys.platform.startswith('darwin'):
        subprocess.call(('open', fname))
    elif os.name == 'nt':
        os.startfile(fname)
    elif os.name == 'posix':
        subprocess.call(('xdg-open', fname))


def mangle_fname(top_dir, fname, pre_ext):
    base, ext = os.path.splitext(os.path.basename(fname))
    return os.path.join(top_dir, base + pre_ext + ext)


class CacheSaver(object):

    def __init__(self, fname):
        self.fname = fname
        if not os.path.isfile(self.fname):
            self.data = {}
        else:
            with open(self.fname) as f:
                self.data = json.load(f)

    def dump(self):
        with open(self.fname, 'w') as f:
            json.dump(self.data, f)


def insert_best_seqid_column(params):
    csv = params['csv']
    uniprot_ids_header = params['uniprot_ids_header']
    delimiter = params['delimiter']
    out_csv = params['output_csv']
    uniprot_cache = params['cache_dir']

    if not csv:
        raise IOError('No file selected')

    if not os.path.isfile(csv):
        raise IOError(csv + ' not found')

    headers = get_headers(csv)

    if uniprot_ids_header not in headers:
        s = "Column header '%s' not found, available headers:\n" % uniprot_ids_header
        for header in headers:
            s += '   ' + header + '\n'
        raise IOError(s)

    logging('Reading %s\n' % csv)
    entries = read_csv(csv)

    all_seqids = []
    for entry in entries:
        tokens = entry[uniprot_ids_header].split(delimiter)
        entry['seqids'] = [s.strip() for s in tokens]
        all_seqids.extend(entry['seqids'])

    logging('Found %d potential Uniprot IDs\n' % len(all_seqids))

    uniprot_data = uniprot.batch_uniprot_metadata(
        all_seqids, uniprot_cache)

    for entry in entries:
        best_seqid = uniprot.sort_seqids_by_uniprot(
            entry['seqids'], uniprot_data)[0]
        entry['best_seqid'] = best_seqid
        entry['is_reviewed'] = False
        if best_seqid in uniprot_data:
            entry['is_reviewed'] = \
                uniprot_data[best_seqid]['is_reviewed']

    logging('Writing ')
    logging('%s\n' % os.path.abspath(out_csv), lambda: open_file(out_csv))
    headers = ['best_seqid', 'is_reviewed'] + get_headers(csv)
    rows = [headers]
    for entry in entries:
        rows.append([entry[h] for h in headers])
    write_csv(out_csv, rows)


class BestSeqidForm(tkform.Form):

    def __init__(self, width=700, height=800):
        tkform.Form.__init__(self, 'Best Uniprot ID Extractor', width, height)

        self.push_text("Best Uniprot ID Extractor", 30)

        self.push_line(width=600)
        self.push_text("(C) 2015 Bosco K Ho")

        text = """
        This program will search the Uniprot IDs in a table (csv|txt|tsv) against http://uniprot.org 
        and identify the best Uniprot ID - the longest protein that has been reviewed.

        The uniprot data will be cached (in case of interruption) and the results will be written
        to a table with two new columns: 'best_seqid' and 'is_reviewed'.
        """
        for line in text.splitlines():
          self.push_text(line.strip())

        self.push_labeled_param(
            'csv', 
            'Input table', 
            load_file_text='select (csv|txt|tsv)',
            width=50)
        self.push_labeled_param(
            'uniprot_ids_header', 
            'Column with Uniprot IDs',
            entry_text='Majority protein IDs',
            width=50)
        self.push_labeled_param(
            'delimiter', 
            'Delimiter between Uniprot IDs in table cell',
            entry_text=';',
            width=4)
        self.push_labeled_param(
            'cache_dir', 
            'Cache directory', 
            entry_text='uniprot.cache',
            load_dir_text='select',
            width=50)
        self.push_labeled_param(
            'output_csv', 
            'Output table', 
            entry_text='output.csv',
            load_file_text='select',
            width=50)

        self.push_spacer()
        self.push_line(width=600)

        self.push_text("Output", 16)
        self.push_submit()
        self.push_output()

    def run(self, params):
        self.clear_output()
        global logging
        logging = self.print_output
        uniprot.logging = self.print_output
        insert_best_seqid_column(params)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print __doc__

    elif sys.argv[1] == '-i':
        # run the GUI
        BestSeqidForm(800, -50).mainloop()

    elif sys.argv[1] == '-t':
        # test set of parameters
        params = {
            'csv': 'at3_150310.csv',
            'uniprot_ids_header': 'Majority protein IDs',
            'cache_dir': 'uniprot.cache',
            'delimiter': ';',
            'output_csv': 'at3_150310.best_seqid.csv'
        }
        insert_best_seqid_column(params)

    else:
        # run in command-line mode with a params file
        params_fname = sys.argv[1]
        logging('Loading params ' + params_fname)
        params = json.load(open(params_fname, 'Ur'))
        collate_peptide_summaries(params)







