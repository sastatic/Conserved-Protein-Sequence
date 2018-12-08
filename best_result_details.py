import csv
import os
from Bio import SeqIO


SRC_DATA = './best_n/'

def main():
	files = os.listdir(SRC_DATA)
	lst = []

	for file in files:
		ls = []
		nm = []
		for seq in SeqIO.parse(SRC_DATA + file, 'fasta'):
			s = seq.description
			s = s.split('|')[-1]
			o = s.find('OS')
			name = s[:o]
			ox = s.find('OX')
			org = s[o+3:ox].split('(')[0]
			nm.append(name)
			ls.append(org)
		lst.append(nm)
	with open('result.csv', mode='w') as csv_file:
		writer = csv.DictWriter(csv_file, fieldnames=ls)
		writer.writeheader()
		for proteins in lst:
			d = dict()
			for i, protein in enumerate(proteins):
				d[ls[i]] = protein
			writer.writerow(d)

if __name__ == '__main__':
	main()