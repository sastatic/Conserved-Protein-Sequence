import os
from Bio import SeqIO
import time
from Bio import Align
import multiprocessing
from numpy import array_split
from Bio.SubsMat.MatrixInfo import blosum62

TEST_PATH = "./non_redundent_protein/"
ALN_PATH = "./global_aligned_proteins/"

alnr = Align.PairwiseAligner()
alnr.mode = 'global'
alnr.extend_gap_score = -1
alnr.open_gap_score = -3
alnr.match = 1
alnr.mismatch = -1

def minimum_seq_file():
	names = os.listdir(TEST_PATH)
	nd = dict()
	for name in names:
		val = int(name.split('-')[0])
		nd[val] = name
	nd = sorted(nd.items())
	names = []
	for _, v in iter(nd):
		names.append(v)
	nl = []
	for a in names[1:]:
		nl.append([names[0], a])
	return nl


def function_aligning(x):
	print ("Process Start.")
	fl1 = x[0]
	fl2 = x[1]
	output_file_name = ALN_PATH + fl2.split('-')[1]
	t3 = time.time()
	with open(output_file_name, 'w') as outFile:
		for seq1 in SeqIO.parse(TEST_PATH + fl1, 'fasta'):
			mxscr = -999999999
			for seq2 in SeqIO.parse(TEST_PATH + fl2, 'fasta'):
				score = alnr.score(seq1.seq, seq2.seq)
				if (mxscr < score):
					mxscr = score
					res = seq2
			SeqIO.write(res, outFile, 'fasta')
	t4 = time.time()
	print ("Took " + str(t4-t3) + " seconds to ", end='\t')
	print ("complete 'GLOBAL' alignment for: " + fl2.split('-')[-1])


if __name__ == '__main__':
	t0 = time.time()

	nl = minimum_seq_file()
	p2 = multiprocessing.Pool(multiprocessing.cpu_count())
	p2.map(function_aligning, nl)

	t1 = time.time()
	print(t1-t0)