import sys
import os
import csv
import json
from pprint import pprint
import glob
import os
import re
import shutil
import copy

from PIL import Image, ImageDraw
import datafile

import uniprot
import datafile


def transfer_files(in_dir, out_dir):
  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
  for src in glob.glob(os.path.join(in_dir, '*')):
    dst = os.path.join(out_dir, os.path.basename(src))
    shutil.copy(src, dst)


def default_protein():
    return {
        'sequence': '', 
        'sources': [
            { 
              'peptides': [], 
              'intervals':[], 
              'color': 'red',
            }
            for i in range(3)
        ],
    }


def make_csv(data, csv_fname):
    rows = []
    header = ['seqid', 'description' 'n_res', 'DP', 'DQ', 'DR', 'all']
    rows.append(header)
    for seqid in data:
        protein = data[seqid]
        n_res = protein['length']
        sources = protein['sources']
        n_source = len(sources)

        row = [seqid, protein['description']]
        coverage = [[0]*n_res for i in range(n_source+1)]
        for i_source in range(n_source):
            source = sources[i_source]
            for i, j in source['intervals']:
                for k in range(i, j):
                  coverage[i_source][k] = 1
                  coverage[-1][k] = 1
            row.append(sum(coverage[i_source])/float(n_res)*100.0)
        row.append(sum(coverage[-1])/float(n_res)*100.0)
        rows.append(row)  
    datafile.write_csv(rows, csv_fname)


def make_residues(n_res, sources, permutation):
    residues = [0 for i in range(n_res)]
    for i_source, source in enumerate(sources):
        source_residues = [0 for i in range(n_res)]
        if i_source not in permutation:
            continue
        for i_res, j_res in source['intervals']:
            for k_res in range(i_res, j_res):
                source_residues[k_res] = 1
        for i in range(n_res):
            residues[i] += source_residues[i]
    return residues


def make_intervals(residues):
    result = []
    i = 0
    curr_val = residues[i]
    n_res = len(residues)
    for j in range(1, n_res + 1):
        if j == n_res:
            result.append((i, j+1, curr_val))
            break
        if curr_val == residues[j]:
            continue
        result.append((i, j, curr_val))
        i = j
        curr_val = residues[i]
    return result


def make_png_from_data(data_json, out_png, permutation):
    data = datafile.load_json(data_json)
    n_protein = len(data)
    im = Image.new('RGB', (n_protein, n_protein), "white")

    draw = ImageDraw.Draw(im)

    for i_protein, (seqid, protein) in enumerate(data.items()):
        n_res = protein['length']
        residues = make_residues(n_res, protein['sources'], permutation)
        intervals = make_intervals(residues)
        for i_res, j_res, val in intervals:
            if val == 0:
                continue
            elif val == 1:
                color = "pink"
            elif val > 1:
                color = "blue"
            i = int(i_res/float(n_res)*n_protein)
            j = int(j_res/float(n_res)*n_protein)
            draw.line((i, i_protein, j, i_protein), fill=color)
            pass

    del draw

    im.save(out_png)


def make_protein_copy(proteins, source_indices):
    result = {}
    for seqid, protein in proteins.items():
        new_protein = {}
        for key, val in protein.items():
            if key == 'sources':
                new_protein['sources'] = []
                for i_source, source in enumerate(protein['sources']):
                    new_protein['sources'].append(copy.deepcopy(source))
                    if i_source not in source_indices:
                        new_protein['sources'][i_source]['peptides'] = []
                        new_protein['sources'][i_source]['intervals'] = []
                        print "Haha", i_source
            else:
                new_protein[key] = val
        result[seqid] = new_protein
    return result

    
for skip in ['', '.skip_dr4']:
    exp_sets = [
        ('DP', 'red', 0, list(range(1, 5))), 
        ('DQ', 'blue', 1, list(range(1, 6))), 
        ('DR', 'green', 2, list(range(1, 5)))]

    source_sets = []
    for cell_type, color, i_source, repeats in exp_sets:
        for i_repeat in repeats:
            if skip == '_skip_dr4':
                if cell_type =='DR' and i_repeat == 4:
                    continue
            source_sets.append(('%s%d' % (cell_type, i_repeat), color, i_source))
    print source_sets
    
    protein = {}
    for exp, color, i_source in source_sets:
        fname = 'Data_PeptideOverlay/%s_motif.csv' % exp
        for entry in datafile.read_csv(fname):
            seqid = uniprot.parse_fasta_header(entry['Accessions'])[0]
            if seqid not in protein:
                protein[seqid] = default_protein()
                protein[seqid]['description'] = entry['Names']
            source = protein[seqid]['sources'][i_source]
            source['color'] = color
            source['peptides'].append(entry['sequence'])

    seqids, fasta = uniprot.read_fasta('../db/uniprot_sprot.fasta')

    for seqid in protein:
        sequence = fasta[seqid]['sequence']
        protein[seqid]['sequence'] = sequence
        protein[seqid]['length'] = len(sequence)
        for source in protein[seqid]['sources']:
            for peptide in source['peptides']:
                i = sequence.index(peptide)
                j = i + len(peptide)
                source['intervals'].append([i, j])
            del source['peptides']

    print 'write overlay%s.csv' % skip
    make_csv(protein, 'overlay%s.csv' % skip)

    for permutation_group, permutation in [
        ('DP-DR', (0, 2)),
        ('DR-DQ', (1, 2)),
        ('DP-DQ', (0, 1)),
        ('common', (0, 1, 2))
    ]:
        out_dir = 'webapp.%s%s' % (permutation_group, skip)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        new_protein = make_protein_copy(protein, permutation)
        jsonp = os.path.join(out_dir, 'load_data.jsonp')
        print "Write json for webpage", jsonp
        datafile.write_jsonp(new_protein, jsonp, 'load_data')
        transfer_files('webapp.template', out_dir)

        data_json = 'out.%s%s.json' % (permutation_group, skip)
        out_png = 'out.%s%s.png' % (permutation_group, skip)
        print 'write %s' % out_png
        datafile.write_json(protein, data_json)
        make_png_from_data(data_json, out_png, permutation)


