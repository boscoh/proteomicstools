 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint

import math
import os
import json
import copy
import glob
import shutil
from pprint import pprint
import logging


import fasta

logger = logging.getLogger('proteins')


this_dir = os.path.abspath(os.path.dirname(__file__))


def new_protein(seqid):
  return {
    'attr': { 'seqid': seqid },
    'sources': [{'peptides': [], }],
  }


def new_peptide(peptide_sequence):
  return {
    'sequence': peptide_sequence,
    'intensity': 1,
    'attr': {}
  }


def peptide_iter(proteins):
  for seqid in proteins:
    protein = proteins[seqid]
    for source in protein['sources']:
      for peptide in source['peptides']:
        yield seqid, peptide


def determine_unique_peptides(proteins):
  seqids_by_seq = {}
  for seqid, peptide in peptide_iter(proteins):
    seq = peptide['sequence']
    if seq not in seqids_by_seq:
      seqids_by_seq[seq] = set()
    seqids_by_seq[seq].add(seqid)
  for seqid, peptide in peptide_iter(proteins):
    seq = peptide['sequence']
    n_seqid = len(seqids_by_seq[seq])
    peptide['attr']['is_unique'] = (n_seqid == 1)


def check_missing_fields(proteins):
  for seqid, peptide in peptide_iter(proteins):
    if 'modifications' not in peptide['attr']:
      peptide['attr']['modifications'] = []
    if 'intensity' not in peptide:
      peptide['intensity'] = 1
    if 'mask' not in peptide:
      peptide['mask'] = 0


def count_peptides(proteins, n_peptide_cutoff=1, is_skip_no_unique=False):
  for seqid in proteins.keys():
    protein = proteins[seqid]
    n_peptide = 0
    n_slice_populated = 0
    n_unique_peptide = 0
    for source in protein['sources']:
      peptides = source['peptides']
      unique_peptides = [p for p in peptides if p['attr']['is_unique']]
      n_unique_peptide += len(unique_peptides)
      n_peptide += len(peptides)
      if len(peptides) > 0:
        n_slice_populated += 1
    protein['attr']['n_slice_populated'] = n_slice_populated
    protein['attr']['n_unique_peptide'] = n_unique_peptide
    protein['attr']['n_peptide'] = n_peptide
    if n_peptide < n_peptide_cutoff:
      del proteins[seqid]
      continue
    if is_skip_no_unique:
      if protein['attr']['n_unique_peptide'] == 0:
        del proteins[seqid]


def find_peptide_proteins(proteins):
  n_source = len(proteins.values()[0]['sources'])
  seqids_by_sequence = {}
  for seqid in proteins:
    protein = proteins[seqid]
    for i_source in range(n_source):
      source = protein['sources'][i_source]
      peptides = source['peptides']
      for peptide in peptides:
        sequence = peptide['sequence']
        if sequence not in seqids_by_sequence:
          seqids_by_sequence[sequence] = set()
        seqids_by_sequence[sequence].add(seqid)
  for seqid in proteins:
    protein = proteins[seqid]
    for i_source in range(n_source):
      protein = proteins[seqid]
      source = protein['sources'][i_source]
      peptides = source['peptides']
      for peptide in peptides:
        sequence = peptide['sequence']
        for test_seqid in seqids_by_sequence[sequence]:
          if seqid != test_seqid:
            if 'other_seqids'not in peptide['attr']:
              peptide['attr']['other_seqids'] = []
            peptide['attr']['other_seqids'].append(test_seqid)
            logger.debug('{} in {} also found in {}'.format(sequence, seqid, test_seqid))


def change_seqids_in_proteins(proteins, clean_seqid):
  seqids = proteins.keys()
  for seqid in seqids:
    new_seqid = clean_seqid(seqid)
    if seqid != new_seqid:
      proteins[new_seqid] = proteins[seqid]
      del proteins[seqid]
      if 'attr' in proteins[new_seqid]:
        if 'seqid' in proteins[new_seqid]['attr']:
          proteins[new_seqid]['attr']['seqid'] = new_seqid


def load_fastas_into_proteins(
    proteins, fastas, clean_seqid=None, iso_leu_isomerism=False):
  if clean_seqid:
    change_seqids_in_proteins(proteins, clean_seqid)
    change_seqids_in_proteins(fastas, clean_seqid)
  for seqid in proteins.keys():
    protein = proteins[seqid]
    if seqid not in fastas:
      logger.debug("%s not found in fasta database" % seqid)
      del proteins[seqid]
      continue
    protein_sequence = fastas[seqid]['sequence']
    protein['description'] = fastas[seqid]['description']
    protein['sequence'] = protein_sequence
    if iso_leu_isomerism:
      protein_sequence = protein_sequence.replace("L", "I")
    n_peptide = 0
    for source in protein['sources']:
      peptides = source['peptides']
      for i_peptide in reversed(range(len(peptides))):
        peptide = peptides[i_peptide]
        peptide_sequence = peptide['sequence']
        if iso_leu_isomerism:
          peptide_sequence = peptide_sequence.replace("L", "I")
        i = protein_sequence.find(peptide_sequence)
        if i < 0:
          logger.debug("'{}' not found in {}".format(peptide_sequence, seqid))
          del peptides[i_peptide]
          continue
        peptide['i'] = i 
        # peptide['j'] = i + len(peptide_sequence)
      n_peptide += len(peptides)
    if n_peptide == 0:
      del proteins[seqid]


def load_fasta_db_into_proteins(
    proteins, fasta_db, clean_seqid=None, iso_leu_isomerism=False):
  seqids, fastas = fasta.read_fasta(fasta_db)
  load_fastas_into_proteins(proteins, fastas, clean_seqid, iso_leu_isomerism)


def merge_two_proteins(proteins1, proteins2):
  """
  Merges two proteins structures. In particular, it grafts the
  'sources' together, treating the sources in each proteins as 
  distinct, and maintaing the order.
  """
  if len(proteins1) == 0:
    return proteins2
  seqid = proteins1.keys()[0]
  n_source1 = len(proteins1[seqid]['sources'])
  seqid = proteins2.keys()[0]
  n_source2 = len(proteins2[seqid]['sources'])
  for seqid in proteins2:
    protein2 = proteins2[seqid]
    sources2 = protein2['sources']
    if seqid in proteins1:
      protein1 = proteins1[seqid]
      protein1['sources'].extend(sources2)
    else:
      sources1 = [{'peptides': []} for i in range(n_source1)]
      combined_sources = sources1 + sources2
      proteins1[seqid] = protein2
      proteins1[seqid]['sources'] = combined_sources
  for seqid in proteins1:
    peptides = proteins1[seqid]['sources']
    if len(peptides) == n_source1:
      peptides.extend([{'peptides': []} for i in range(n_source2)])
  return proteins1


def save_data_js(data, js_fname):
  f = open(js_fname, 'w')
  f.write('var data = \n')
  f.write(json.dumps(data, indent=None))
  f.close()


def transfer_newer_files(in_dir, out_dir):
  for src in glob.glob(os.path.join(in_dir, '*')):
    dst = os.path.join(out_dir, os.path.basename(src))
    if not os.path.isfile(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
      shutil.copy(src, dst)


def make_peptograph_directory(data, out_dir):
  # sanity checks
  proteins = data['proteins']
  determine_unique_peptides(proteins)
  count_peptides(proteins)
  check_missing_fields(proteins)
  find_peptide_proteins(proteins)
  for seqid, protein in proteins.items():
    for source in protein['sources']:
      peptides = source['peptides']
      peptides.sort(key=lambda peptide: len(peptide['sequence']))
      peptides.sort(key=lambda peptide: peptide['i'])
  if 'source_labels' not in data:
    data['source_labels'] = []
  if 'color_names' not in data:
    data['color_names'] = ['', '', '']
  if 'mask_labels' not in data:
    data['mask_labels'] = []

  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

  save_data_js(data, os.path.join(out_dir, 'data.js'))
  transfer_newer_files(os.path.join(this_dir, 'templates/peptograph'), out_dir)
  index_html = os.path.abspath(os.path.join(out_dir, 'index.html'))
  logger.info("Made peptograph in " + index_html)


def make_proteins_directory(data, out_dir):
  # sanity checks
  proteins = data['proteins']
  determine_unique_peptides(proteins)
  count_peptides(proteins)
  check_missing_fields(proteins)
  find_peptide_proteins(proteins)
  for seqid, protein in proteins.items():
    for source in protein['sources']:
      peptides = source['peptides']
      peptides.sort(key=lambda peptide: len(peptide['sequence']))
      peptides.sort(key=lambda peptide: peptide['i'])
  if 'source_labels' not in data:
    data['source_labels'] = []
  if 'color_names' not in data:
    data['color_names'] = ['', '', '']
  if 'mask_labels' not in data:
    data['mask_labels'] = []

  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

  save_data_js(data, os.path.join(out_dir, 'data.js'))
  transfer_newer_files(os.path.join(this_dir, 'templates/proteins'), out_dir)
  index_html = os.path.abspath(os.path.join(out_dir, 'index.html'))
  logger.info("Made peptograph in " + index_html)





