import csv
import collections

import datafile


def read_peptides_from_csv(fname):
  peptides = {}
  for entry in datafile.read_csv(fname):
    if int(entry['Uni']) == 0:
      continue
    n = entry['Uni']
    seqids = entry['protein'].split(',')
    sequence = entry['Sequence'].split('.')[1]
    if sequence.startswith('n'):
      sequence = sequence[1:]
    if sequence in peptides:
      existing_seqids = set(peptides[sequence]['seqids'])
      if existing_seqids != seqids:
        peptides[sequence]['seqids'] = existing_seqids.union(seqids)
    else:
      peptides[sequence] = { 
        'sequence': sequence,
        'seqids': seqids,
      }
  return peptides


def assign_peptides_to_proteins(peptides):
  proteins = {}
  for sequence in sorted(peptides.keys()):
    for seqid in peptides[sequence]['seqids']:
      if seqid not in proteins:
        proteins[seqid] = {
          'peptide_set': set()
        }
      proteins[seqid]['peptide_set'].add(sequence)
  return proteins


def generate_groups_with_similar_peptides(proteins):
  groups = []
  def get_i_group(peptides):
    for i_group, group in enumerate(groups):
      if group['peptide_set'] == peptides:
        return i_group
    return None
  for seqid in proteins:
    peptides = proteins[seqid]['peptide_set']
    i_group = get_i_group(peptides)
    if i_group is None:
      groups.append({
        'peptide_set': peptides,
        'seqids': [seqid]
        })
    else:
      groups[i_group]['seqids'].append(seqid)
  return groups


def collect_clusters_with_intersecting_peptides(groups):
  i_cluster = None
  clusters = []
  n_group = len(groups)
  for i_group1 in range(n_group):
    group1 = groups[i_group1]
    if 'i_cluster' not in group1:
      clusters.append(set([i_group1]))
      group1['i_cluster'] = len(clusters)-1
    i_cluster = group1['i_cluster']
    for i_group2 in range(i_group1+1, n_group):
      group2 = groups[i_group2]
      for i_group_in_cluster in clusters[i_cluster]:
        group = groups[i_group_in_cluster]
        intersect = group['peptide_set'].intersection(group2['peptide_set'])
        if len(intersect) > 0:
          group2['i_cluster'] = i_cluster
          clusters[i_cluster].add(i_group2)
          break
  # determine if a sibling is a subset or superset of another sibling
  for i_cluster, cluster in enumerate(clusters):
    for i_group1 in cluster:
      group = groups[i_group1]
      peptides1 = group['peptide_set']
      group['superset'] = []
      group['subset'] = []
      for i_group2 in cluster:
        if i_group1 == i_group2:
          continue
        peptides2 = groups[i_group2]['peptide_set']
        if peptides1.issubset(peptides2):
          group['superset'].append(i_group2)
        if peptides1.issuperset(peptides2):
          group['subset'].append(i_group2)
    for i_sibling in range(len(cluster)):
      i_group = list(cluster)[i_sibling]
      group = groups[i_group]
      superset = group['superset']
      subset = group['subset']
      print i_cluster+1, i_sibling+1, i_group, superset, subset, ' '.join(group['peptide_set'])
  return clusters


def extract_primary_seqid(clusters):
  for i_cluster, cluster in enumerate(clusters):
    supersets = []
    subsets = []
    i_group1 = 0
    for i_group1 in cluster:
      group = groups[i_group1]
      seqids = list(group['seqids'])
      seqids.sort(key=lambda s: len(s))
      group['seqid'] = seqids[0]
      group['other_seqids'] = seqids[1:]


def write_out_clusters(fname, groups, clusters):
  f = open(fname, 'w')
  writer = csv.writer(f)
  writer.writerow(['group', 'sibling', 'seqid', 'other_seqids', 'peptides', 'subset_seqids'])
  for i_cluster, cluster in enumerate(clusters):
    for i_sibling, i_group in enumerate(cluster):
      group = groups[i_group]
      if len(group['superset']) > 0:
        continue
      subset = group['subset']
      subset_seqids = []
      for i in subset:
        subset_seqids.extend(groups[i]['seqids'])
      row = [
        i_cluster+1, 
        i_sibling+1, 
        group['seqid'], 
        ','.join(group['other_seqids']), 
        ','.join(group['peptide_set']), 
        ','.join(subset_seqids)]  
      writer.writerow(row)


peptides = read_peptides_from_csv('peptide_list.csv')
proteins = assign_peptides_to_proteins(peptides)
groups = generate_groups_with_similar_peptides(proteins)
clusters = collect_clusters_with_intersecting_peptides(groups)
extract_primary_seqid(clusters)
write_out_clusters('proteins.csv', groups, clusters)

peptides = []
for i_cluster, cluster in enumerate(clusters):
  for i_sibling, i_group in enumerate(cluster):
    group = groups[i_group]
    if len(group['superset']) > 0:
      continue
    sequences = group['peptide_set']
    for sequence in sequences:
      peptide = {
        'sequence': sequence,
        'seqid': group['seqid'],
        'other_seqids': group['other_seqids'],
        'secondary_seqids': []
      }
      peptides.append(peptide)
      if len(group['subset']) > 0:
        for i_group_subset in group['subset']:
          subset_group = groups[i_group_subset]
          if sequence in subset_group['peptide_set']:
            peptide['secondary_seqids'].extend(list(subset_group['seqids']))


fname = 'peptide_sorted.csv'
f = open(fname, 'w')
writer = csv.writer(f)
writer.writerow(['sequence', 'seqid', 'other_seqids', 'secondary_seqids'])
peptides.sort(key=lambda p:p['sequence'])
for p in peptides:
  row = [
    p['sequence'], 
    p['seqid'], 
    ':'.join(p['other_seqids']), 
    ':'.join(p['secondary_seqids'])]
  writer.writerow(row)






