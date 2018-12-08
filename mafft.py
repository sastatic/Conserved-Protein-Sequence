import os
import time
import multiprocessing
from Bio.Align.Applications import MafftCommandline

TEST_PATH = './msa_db2/'
ALN_PATH = "./aligned_protein/"


def function_to_align(name):
	mafft_cline = MafftCommandline(input=TEST_PATH + name)
	stdout, stderr=mafft_cline()

	with open(ALN_PATH + name, 'w') as ifh:
		ifh.write(stdout)


def main():
	t0 = time.time()
	names = os.listdir(TEST_PATH)

	p = multiprocessing.Pool(multiprocessing.cpu_count())
	p.map(function_to_align, names)

	t1 = time.time()
	print (t1-t0)


if __name__ == '__main__':
	main()