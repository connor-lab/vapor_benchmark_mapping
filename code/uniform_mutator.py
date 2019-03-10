import sys
import random
from Bio import SeqIO

# Parse a fasta file of sequences, an index for which one to retrieve, and a mutation probability p
seqf, pid = sys.argv[1:]
p = 1-float(pid)/100.
seqr = [s for s in SeqIO.parse(seqf, "fasta")][0]
seqh = seqr.description
seq = str(seqr.seq)
mutseq = ""
for c in seq:
    roll = random.uniform(0,1)
    if roll < p:
        alt = [b for b in "ATCG" if b != c]
        mutseq += random.choice(alt)
    else:
        mutseq += c

print(">"+seqh)
print(mutseq)
