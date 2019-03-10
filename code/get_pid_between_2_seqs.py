from Bio import pairwise2

def get_pid(q, r):
	aln = pairwise2.align.globalms(q, r, 1,-1,-1,-1, one_alignment_only=True)
	s1 = aln[0][0]
	s2 = aln[0][1]
	c = 0
	for i in range(len(s1)):
		if s1[i] != s2[i]:
			c += 1
	return 100.0 - 100 * float(c)/len(s1)	

if __name__ == "__main__":
	import sys
	f1name = sys.argv[1]
	f2name = sys.argv[2]
	# Will be garbage collected anyway, simple script
	f1seq = "".join([l.strip() for l in open(f1name)][1:])
	f2seq = "".join([l.strip() for l in open(f2name)][1:])
	print(get_pid(f1seq, f2seq))
