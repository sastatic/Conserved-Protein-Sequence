import os
import time
from collections import defaultdict
from Bio import SeqIO

TEST_PATH = "./global_aligned_proteins/"
MSA_PATH = "./msa_db2/"

def get_subject_sequence_set():
	names = os.listdir(TEST_PATH)
	nset = defaultdict(list)

	for name in names:
		count = 0
		for record in SeqIO.parse(TEST_PATH + name, "fasta"):
			nset[count].append(record)
			count += 1

	for i, records in iter(nset.items()):
		PATH = MSA_PATH + str('{0:03d}'.format(i+1)) + '.fasta'
		with open(PATH, 'a') as outFile:
			for record in records:
				SeqIO.write(record, outFile, 'fasta')

if __name__ == '__main__':
	t0 = time.time()
	get_subject_sequence_set()
	t1 = time.time()
	print(t1-t0)
