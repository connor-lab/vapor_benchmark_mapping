""" Retrieves the nth sequence from a fasta """

if __name__ == "__main__":
    import sys
    from Bio import SeqIO
    fname = sys.argv[1]
    n = int(sys.argv[2])
    for i, r in enumerate(SeqIO.parse(fname, "fasta")):
        if i+1 == n:
            print(">"+r.description)
            print(str(r.seq))
            sys.stderr.write("Fetching " + str(n) + r.description+"\n")
            sys.exit(0)
    sys.exit("n < the number of sequences")

        
