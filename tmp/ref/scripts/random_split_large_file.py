
import sys
import os
import re
import random
import argparse # to use this module python version must be higher than 2.7
# parse input arguments
parser = argparse.ArgumentParser(description='This script is to split large file into the specified number of files. \
Copyright (C) 2015 Yasumasa Kimura, Division of Systems Immunology, Institute of Medical Science, The University of Tokyo, Tokyo, Japan.')
parser.add_argument('input_fasta', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input fasta file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 2016.4.22')
parser.add_argument('-o', '--output', action='store', dest='out_base_', default=None, help='basename of output files')
parser.add_argument('-n', '--number_of_split', action='store', dest='n_split', type=int, default=None, help='the number of output files splitted')
args = parser.parse_args()


if args.out_base_:
	out_ = args.out_base_
elif args.input_fasta.name != '<stdin>':
	out_ = os.path.basename(args.input_fasta.name)
else:
	out_ = "noname"

a_f_out = []
if args.n_split is None:
	print 'Option -n is required. Please specifiy the number of output files splitted.'
else:
	for i in range(0, args.n_split):
		f_out_ = open("%s.group_%i.fa" % (out_, i+1), "w")
		a_f_out.append(f_out_)

print out_

from Bio import SeqIO
record_iter = SeqIO.parse(args.input_fasta, "fasta")
for entry in record_iter:
	i = random.randint(1,args.n_split)
	a_f_out[i - 1].write(">%s\n%s\n" % (entry.id, entry.seq))

# file close
for i in range(0, args.n_split):
	a_f_out[i].close()

