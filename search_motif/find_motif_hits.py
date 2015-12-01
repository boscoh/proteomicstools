

import copy
import math
import os
import matplotlib.pylab as pylab

import datafile
import uniprot
import glob


"""
Calculates the amino acid propensities from a fasta file.
"""



def print_matrix(matrix, header_fmt=" %5s", val_fmt=" % 5.2f"):

    keys1 = matrix.keys()
    keys2 = matrix[keys1[0]].keys()

    header = "   "
    for key1 in keys1:
        header += header_fmt % key1
    print header

    for key2 in keys2:
        s = key2 + "  "
        for key1 in keys1:
            s += val_fmt % matrix[key1][key2]
        print s



def count_aa_in_fasta(fasta, counts_txt):

    if os.path.isfile(counts_txt):
        with open(counts_txt) as f:
            return eval(f.read())

    seqids, fastas = uniprot.read_fasta(fasta)

    counts = {}
    n = 0
    for entry in fastas.values():
        seq = entry['sequence']
        for c in seq:
            if c not in counts:
                counts[c] = 0
            counts[c] += 1
        n += len(seq)

    with open(counts_txt, 'w') as f:
        f.write(repr(counts))

    return counts



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
        open(sub_csv + '.count', 'w').write(
            str(len(sub_sequences)))


def read_blosum_txt(blosum_txt):
    blosum = {}
    titles = None
    for line in open(blosum_txt):
        if line.startswith('#'):
            continue
        if titles is None:
            titles = line.split()
            for aa1 in titles:
                blosum[aa1] = {}
                for aa2 in titles:
                    blosum[aa1][aa2] = 0.0
        else:
            tokens = line.split()
            aa1 = tokens[0]
            for aa2, val in zip(titles, tokens[1:]):
                blosum[aa1][aa2] = float(val)
    return blosum



def read_freqs_csv(freqs_csv):
    freqs = {}
    for entry in datafile.read_csv(freqs_csv):
        aa = entry["Amino Acid"]
        if aa == "Total":
            continue
        for key, val in entry.items():
            if key.startswith("P"):
                i = int(key[1:])
                if i not in freqs:
                    freqs[i] = {}
                freqs[i][aa] = float(val)/100.
    return freqs



def get_motif_log_odds(
    motif_csv, aa_probs, 
    alpha=1.0, beta=0.01, variant="fixed"):

    raw_freqs = read_freqs_csv(motif_csv)
    # print_matrix(raw_freqs)
    aa_list = raw_freqs[raw_freqs.keys()[0]].keys()

    pseudo_freqs = copy.deepcopy(raw_freqs)

    if variant == "const":
        for p in pseudo_freqs:
            for aa in pseudo_freqs[p]:
                pseudo_freqs[p][aa] = 1.0
    elif variant == "prob":
        for p in pseudo_freqs:
            for aa in pseudo_freqs[p]:
                pseudo_freqs[p][aa] = aa_probs[aa]
    elif variant == "subs":
        # http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        blosum = read_blosum_txt('blosum62.txt')
        # print_matrix(blosum, header_fmt=" %4s", val_fmt=" %4d")

        # pseudo_freqs calculates the "background" freqs 
        # from aa substitution and will be used to avoid f=0 for 
        # log-odds calculation
        for p in pseudo_freqs:
            for aa1 in pseudo_freqs[p]:
                substitution_prob = 0.0
                for aa2 in aa_list:
                    odds_of_substitution = math.pow(2, 0.5 * blosum[aa1][aa2])
                    substitution_prob += raw_freqs[p][aa2] * odds_of_substitution
                pseudo_freqs[p][aa1] = aa_probs[aa1] * substitution_prob
        # print_matrix(pseudo_freqs)
    else:
        raise Exception

    freqs = copy.deepcopy(raw_freqs)
    for p in raw_freqs:
        for aa in raw_freqs[p]:
            freq = raw_freqs[p][aa]
            pseudo_freq = pseudo_freqs[p][aa]
            corrected_freq = (alpha*freq + beta*pseudo_freq)/(alpha + beta)
            freqs[p][aa] = corrected_freq
    # print_matrix(freqs)

    log_odds = copy.deepcopy(freqs)
    for p in freqs:
        for aa in freqs[p]:
            odds = freqs[p][aa]/aa_probs[aa]
            log_odds[p][aa] = math.log(odds)
    # print_matrix(log_odds)

    return log_odds



def score_seq(log_odds, seq):
    n_position = len(log_odds.keys())
    score = 0
    for j in range(0, n_position):
        p = j+1
        score += log_odds[p][seq[j]]
    return score


def score_full_sequence(log_odds_pssm, test_seq):
    n_position = len(log_odds_pssm.keys())
    scores = []
    n_test_seq = len(test_seq) - n_position
    for i in range(n_test_seq):
        seq = test_seq[i:i+n_position]
        score = score_seq(log_odds_pssm, seq)
        scores.append(score)
    return scores


def make_html_for_seq_scores(test_seq, n_position, scores):

    max_score = max(scores)
    min_score = min(scores)

    half_width = max(abs(max_score), abs(min_score))
    delta = 1.0/half_width*100

    def td(width, text, color=""):
        if color:
            color = "background-color:%s" % color
        text = "<td style='width:%dpx;%s'>%s</td>" % (width, color, text)
        return text

    piece = ""
    for i in range(0, len(test_seq) - n_position):
        piece += "<table style='table-layout:fixed;'>"
        seq = test_seq[i:i+n_position]
        s = scores[i]
        piece += td(40, i+1)
        piece += td(160, seq)
        piece += td(50, "%.2f" % s)
        if s < 0:
            piece += td((half_width + s)*delta, "", "lightgrey")
            piece += td(-s*delta, "", "blue")
            piece += td(half_width*delta, "", "lightgrey")
        else:
            piece += td(half_width*delta, "", "lightgrey")
            piece += td(s*delta, "", "red")
            piece += td((half_width - s)*delta, "", "lightgrey")
        piece += "</table>"

    return piece


def write_log_odds(log_odds_pssm, out_txt):
    aa_list = log_odds_pssm[log_odds_pssm.keys()[0]].keys()
    with open(out_txt, 'w') as f:
        f.write('Last position-specific scoring matrix computed\n')
        for aa in aa_list:
            f.write('\t %s' % aa)
        f.write('\n')
        for p in log_odds_pssm:
            f.write("%d A" % p)
            for aa in log_odds_pssm[p]:
                f.write("\t%.3f" % log_odds_pssm[p][aa])
            f.write('\n')



def get_aa_probs_from_fasta(background_fasta):
    cache = background_fasta + '.probs.txt'
    aa_counts = count_aa_in_fasta(background_fasta, cache)
    aa_str = "ACEDGFIHKMLNQPSRTWVY"
    aa_list = [c for c in aa_str]
    n = sum(aa_counts[aa] for aa in aa_list)
    aa_probs = { aa: float(aa_counts[aa])/n for aa in aa_list }
    # for aa, p in aa_probs.items():
    #     print "%s: %.f%%" % (aa, p*100)
    return aa_probs


def get_aa_probs_from_csv(freqs_csv):
    freqs = {}
    for line in open(freqs_csv, "Ur"):
        if line.startswith("Amino acid"):
            continue
        words = line.split(",")
        freqs[words[0]] = float(words[1])

    total = float(sum(freqs.values()))
    for aa in freqs:
        freqs[aa] = freqs[aa] / total

    return freqs




def run(
    motif_csv, test_fasta, out_html, aa_probs,
    variant='fixed', 
    alpha=1.0, beta=0.01):

    log_odds_pssm = get_motif_log_odds(
        motif_csv, aa_probs, 
        alpha=alpha, beta=beta, 
        variant=variant)

    write_log_odds(log_odds_pssm, motif_csv.replace('.pssm', '').replace('.csv', '.pssm'))

    n_position = len(log_odds_pssm.keys())

    seqids, fastas = uniprot.read_fasta(test_fasta)

    saved_scores = []

    html = "<body>"
    for seqid, entry in fastas.items():
        # print "Sequence", seqid
        test_seq = entry['sequence']
        scores = score_full_sequence(log_odds_pssm, test_seq)
        html += "<h1>%s</h1>" % entry['description']
        html += make_html_for_seq_scores(test_seq, n_position, scores)
        saved_scores.append({
            'seqid': seqid,
            'description': entry['description'],
            'scores': scores,
            'seq': test_seq
        })
    html += "</body>"
    open(out_html, 'w').write(html)

    titles = ['seqid', 'description', 'position', 'seq', 'score']
    rows = [titles]
    for entry in saved_scores:
        n_seq = len(entry['seq'])
        for i in range(0, n_seq - n_position):
            score = entry['scores'][i]
            pos = i + 1
            if score > 0:
                row = [
                    entry['seqid'], 
                    entry['description'],
                    pos, 
                    score, 
                    entry['seq'][i:i+n_position]
                ]
                rows.append(row)
    datafile.write_csv(rows, out_html + '.csv')


# make_motif_csv_from_aln(
#     'human/gibbscluster/gibbs.1g.aln', 
#     'gibbs.pssm.csv')

# params = """
# gibbs 3050
# 8mer 29
# 9mer 1166
# 10mer 730
# 11mer 311
# 12mer 343
# """

# if not os.path.isdir('human/result'):
#     os.makedirs('human/result')

# for variant in ['const', 'prob', 'subs']:
#     for name, n in pairs:
#         for beta_type in ['beta_0.1', 'beta_0.01', 'beta_n']:
#             if beta_type == 'beta_0.01':
#                 beta = 0.01
#             elif beta_type == 'beta_0.1':
#                 beta = 0.1
#             else:
#                 beta = 1./float(n)
#             csv = 'human/{}.pssm.csv'.format(name)
#             html = 'human/result/{}.{}.{}.html'.format(name, variant, beta_type)
#             run(csv, 'human/test.fasta', html, variant, beta=beta)



# pairs = [
#     ("AW_10mer", 36),
#     ("AW_9mer", 82),
#     ("AW_8mer", 21),
# ]

# if not os.path.isdir('bat/result'):
#     os.makedirs('bat/result')

# # aa_probs = get_aa_probs_from_fasta('../db/uniprot_sprot.fasta')
# aa_probs = get_aa_probs_from_csv('bat/1_BatProteome_AAfreq_AW05102015.csv')

# for variant in ['const', 'prob', 'subs']:
#     for name, n in pairs:
#         for beta_type in ['beta_0.1', 'beta_0.01', 'beta_n']:
#             if beta_type == 'beta_0.01':
#                 beta = 0.01
#             elif beta_type == 'beta_0.1':
#                 beta = 0.1
#             else:
#                 beta = 1./float(n)
#             csv = 'bat/{}.csv'.format(name)
#             html = 'bat/result/{}.{}.{}.html'.format(name, variant, beta_type)
#             run(csv, 'bat/test.fasta', html, aa_probs, variant, beta=beta)



if not os.path.isdir('bat_mhc_alleles/result'):
    os.makedirs('bat_mhc_alleles/result')

for exp in ["M1_classIpeps", "M2_classIpeps", "UK_classIpeps"]:

    make_motif_csv_from_txt(
        'bat_mhc_alleles/%s.txt' % exp, 
        'bat_mhc_alleles/%s.csv' % exp)

    aa_probs = get_aa_probs_from_csv(
        'bat/1_BatProteome_AAfreq_AW05102015.csv')

    for mer in ['8mer', '9mer', '10mer']:
        name = exp + '.' + mer
        for beta_type in ['beta_0.01', 'beta_n']:
            csv = 'bat_mhc_alleles/{}.csv'.format(name)
            count = csv + '.count'
            if beta_type == 'beta_0.01':
                beta = 0.01
            else:
                n = int(open(count).read())
                beta = 1./float(n)
            html = 'bat_mhc_alleles/result/{}.{}.html'.format(name, beta_type)
            run(csv, 'bat_mhc_alleles/2_HendraVirusProteome.fasta', html, aa_probs, 'prob', beta=beta)



