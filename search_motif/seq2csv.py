
aa_list = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
    'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
    'T', 'V', 'W', 'Y']


def read_sequences_from_aln(fname):
    sequences = []
    for line in open(fname):
        if not line.startswith('Seq'):
            continue
        piece = line.split()[1]
        chars = [c for c in piece if c.isupper()]
        sequences.append("".join(chars))
    return sequences


def convert_sequences_to_motfi(sequences):
    n = len(sequences[0])

    total_aa = len(sequences)*n

    motif = {}
    for i in range(n):
        p = "P{}".format(i+1)
        motif[p] = {}
        for aa in aa_list:
            motif[p][aa] = 0
        for sequence in sequences:
            aa = sequence[i]
            motif[p][aa] += 1
        for aa in motif[p]:
            motif[p][aa] *= 100.0/float(len(sequences))

    return motif


def write_motif_to_csv(fname, motif):
    titles = ['Amino Acid'] + list(motif.keys())
    f = open(fname, 'w')
    f.write(','.join(titles) + '\n')
    for aa in aa_list:
        row = [aa]
        for p in titles[1:]:
            row.append(str(motif[p][aa]))
        f.write(','.join(row) + '\n')


def make_motif_csv_from_aln(aln, csv):
    sequences = read_sequences_from_aln(aln)
    motif = convert_sequences_to_motfi(sequences)
    write_motif_to_csv(csv, motif)


def make_motif_csv_from_txt(txt, csv):
    sequences = open(txt, 'Ur').read().split()
    seq_lens = map(len, sequences)
    seq_lens = list(set(seq_lens))
    for n in [8, 9, 10]:
        mer = str(n) + 'mer'
        sub_sequences = [s for s in sequences if len(s) == n]    
        sub_csv = csv.replace('.csv', '.' + mer + '.csv')
        motif = convert_sequences_to_motfi(sub_sequences)
        write_motif_to_csv(sub_csv, motif)
        print len(sub_sequences), sub_csv


# make_motif_csv_from_aln(
#     'human/gibbscluster/gibbs.1g.aln', 
#     'gibbs.pssm.csv')


make_motif_csv_from_txt(
    'bat_mhc_alleles/M1_classIpeps.txt', 
    'bat_mhc_alleles/M1_classIpeps.csv')