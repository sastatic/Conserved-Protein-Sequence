import os
import time
from shutil import copyfile

SRC_DATA = './aligned_protein/'

def main():
	nd = dict()
	with open('scores', 'r') as fl:
		for line in fl:
			[file, score] = line.split(',')
			score = int(score)
			nd[score] = file
	nd = sorted(nd.items(), reverse=True)
	n = int(input("Enter Number Of Desired Sequence: "))
	for score, file in iter(nd):
		if n == 0:
			break
		copyfile(SRC_DATA + file, './score_wise/' + str(score) + '.fasta')
		copyfile(SRC_DATA + file, './best_n/' + file)
		n -= 1



if __name__ == '__main__':
	main()