import os
import time
import multiprocessing
from shutil import copyfile
from Bio import AlignIO
from Bio import SeqIO
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62

ALN_PATH = "./msa_db2/"
# ALN_PATH = "./aligned_protein/"

alnr = Align.PairwiseAligner()
alnr.mode = 'local'
alnr.substitution_matrix = blosum62
alnr.extend_gap_score = -0.123
alnr.open_gap_score = -1.53

class scoreobj():
	def __init__(self, name, score = 0):
		self.name = name
		self.score = score



def alignment_score(obj):
	name = obj.name
	s = []
	for seq in SeqIO.parse(ALN_PATH + name, 'fasta'):
		s.append(seq)
	score = 0
	for i in range(len(s)):
		for j in range(i+1,len(s)):
			score += alnr.score(s[i], s[j])	
	obj.score = score

def main():
	t0 = time.time()
	names = os.listdir(ALN_PATH)

	lst = []
	for name in names:
		obj = scoreobj(name)
		lst.append(obj)

	p = multiprocessing.Pool(multiprocessing.cpu_count())
	p.map(alignment_score, lst)

	n = int(input("Enter Number Of Desired Sequence: "))

	lst.sort(key = lambda x: x.score, reverse=True)
	for i, a in enumerate(lst):
		if n == 0:
			break
		n -= 1
		copyfile(ALN_PATH + a.name, './score_wise/' + str(i+1) + '.fasta')
		copyfile(ALN_PATH + a.name, './best_n/' + a.name)

	t1 = time.time()

if __name__ == '__main__':
	main()
