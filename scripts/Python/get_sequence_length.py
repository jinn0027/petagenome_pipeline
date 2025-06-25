import sys
import os
import argparse # to use this module python version must be higher than 2.7

# parse input arguments
parser = argparse.ArgumentParser(description='This script is to output sequence length in fasta or fastq file. \
Copyright (C) 2015 Yasumasa Kimura, Division of Systems Immunology, Institute of Medical Science, The University of Tokyo, Tokyo, Japan.')
parser.add_argument('input_', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input file name')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 2015.06.16')
parser.add_argument('-t', '--file_type', action='store', dest='file_type', type=str, default="fasta", help='input file type [ fasta | fastq ]. [default=%(default)s]')
args = parser.parse_args()

from Bio import SeqIO
record_iter = SeqIO.parse(args.input_, args.file_type)
for entry in record_iter:
	id = entry.id
	length = len(entry.seq)
	print( "%s\t%d" % (id, length))

