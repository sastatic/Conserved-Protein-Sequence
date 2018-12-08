#############################
# import all libraries used #
#############################

import os
from Bio import SeqIO
import time
import multiprocessing
from numpy import array_split

#########################################
# setting up data path and destination  #
# path for non redundent proteins		#
#########################################

DATA_PATH = "./proteins/"
TEST_PATH = "./non_redundent_protein/"
Ignore_record = ['HYPOTHETICAL', 'PREDICTED', 'UNCHARACTERIZED', 'PUTATIVE']

def remove_redundent_sequence(input_file_name, output_file_name):
	#############################################################
	# This function is used to remove similar protein sequences #
	#############################################################
	seq_count_before = 0
	seq_count_after = 0
	with open(output_file_name, 'w') as outFile:
		records = set()
		for record in SeqIO.parse(input_file_name, "fasta"):
			seq_count_before += 1
			x = record.seq.__str__()
			flag1 = x not in records
			s = record.description.upper()
			flag2 = not any(x in s for x in Ignore_record)
			if flag1 and flag2:
				seq_count_after += 1
				records.add(x)
				SeqIO.write(record, outFile, 'fasta')
	return (seq_count_before, seq_count_after)


def functionX(names):
	nums = names.__len__()
	for i in range(nums):
		input_file_name = DATA_PATH + names[i]
		output_file_name = TEST_PATH + names[i]
		[seq_count_before, seq_count_after] = remove_redundent_sequence(input_file_name, output_file_name)
		renamed = TEST_PATH + str(seq_count_after) + '-' + names[i]
		os.rename(output_file_name, renamed)


def prepare_non_redundent_data(_path):
	###############################################################
	# This function is used to prepare for non redundent proteins #
	###############################################################
	names = os.listdir(_path)
	names = sorted(names)
	p = multiprocessing.Pool(multiprocessing.cpu_count())
	p.map(functionX, array_split(names, multiprocessing.cpu_count()))


if __name__ == "__main__":
	t0 = time.time()
	prepare_non_redundent_data(DATA_PATH)
	t1 = time.time()
	print("Time taken to complete task is : ", end='')
	print(t1-t0)
	names = os.listdir(TEST_PATH)
	for name in names:
		print(name)
