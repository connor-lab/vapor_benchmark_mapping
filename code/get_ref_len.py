from Bio import SeqIO
import sys

seq = [str(r.seq) for r in SeqIO.parse(sys.argv[1], "fasta")][0]
print(len(seq))
