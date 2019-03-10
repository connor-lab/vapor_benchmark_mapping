""" Retrieves a sequence given its gid from a file of fastas, headers are >gi|X| """

if __name__ == "__main__":
	import random
	import sys
	sequence_fname = sys.argv[1]
	sequences = {}
	with open(sequence_fname) as f:
		temp_header = ""
		temp_seq = ""
		lines = [l for l in f]
		header_indices = [i for i in range(len(lines)) if lines[i][0] == ">"]
		rando = random.choice(header_indices)
		switch = 0
		for i in range(len(lines)):
			line = lines[i]	
			if line[0] == ">":
				if i == rando:
					switch = 1
					print(line.strip())
				else:
					if switch == 1:
						break
			else:
				if switch == 1:
					print(line.strip())


