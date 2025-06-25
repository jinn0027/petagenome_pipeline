def batch_iterator(iterator, args, f_table) :
	entry = True #Make sure we loop once
	i, j = 0, 0

	while entry :
		try :
			#entry = iterator.next()
			entry = next(iterator)
		except StopIteration :
			entry = None
		if entry is None :
			#End of file
			break

		i += 1 # counts of all entries
		if len(entry.seq) >= args.min_length :
			j += 1 # counts of remained entries
			org_name = entry.id
			if args.rename:
				new_name = args.prefix + format(j, '09d')
				entry.id = new_name
				entry.description = ''
			else:
				new_name = org_name
			# write correspondence table between original and new names
			if f_table:
				f_table.write(new_name+'\t'+org_name+'\n')
			yield entry

import sys
import os
import re
import argparse # to use this module python version must be higher than 2.7
# parse input arguments
parser = argparse.ArgumentParser(description='This script is to filter fasta or fastq file. \
Copyright (C) 2015 Yasumasa Kimura, Division of Systems Immunology, Institute of Medical Science, The University of Tokyo, Tokyo, Japan.')
parser.add_argument('f_in', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input file name')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 2018.04.30')
parser.add_argument('-t', '--file_type', action='store', dest='file_type', type=str, default="fasta", help='input file type [fasta/fastq]')
parser.add_argument('-min', '--min_length', action='store', dest='min_length', type=int, default=500, help='minimum sequence length to output')
parser.add_argument('-o', '--output', action='store', dest='out_', default=None, help='output file name')
parser.add_argument('--rename', action='store_true', default=False, help='rename entry id, default:False')
parser.add_argument('--prefix', action='store', dest='prefix', default='Contig', help='prefix of contig name default:Contig')
parser.add_argument('--table', action='store', dest='table_', type=str, default=None, help='file of correspondence table between original and new names')
args = parser.parse_args()

if args.out_:
	f_out = open(args.out_, "w")
else:
	f_out = sys.stdout

if args.table_:
	f_table = open(args.table_, "w")
else:
	f_table = None

from Bio import SeqIO
record_iter = SeqIO.parse(args.f_in, args.file_type)
for i, batch in enumerate(batch_iterator(record_iter, args, f_table)) :
	count = SeqIO.write(batch, f_out, args.file_type)

f_out.close()
if args.table_:
	f_table.close()
