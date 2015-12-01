seqs = open('human/sequences.txt').read().split()
seq_lens = map(len, seqs)

for i in ['gibbs', 8, 9, 10, 11, 12]:
    print i,
    if i == 'gibbs':
        print len(seqs)
    else:
        print seq_lens.count(i)
