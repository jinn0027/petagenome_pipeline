import sys
import os
import re
import argparse # to use this module python version must be higher than 2.7
# parse input arguments
parser = argparse.ArgumentParser(description='This script is to filter fasta file. \
Copyright (C) 2015 Yasumasa Kimura, Division of Systems Immunology, Institute of Medical Science, The University of Tokyo, Tokyo, Japan.')
parser.add_argument('f_in', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input fasta file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 2016.5.24')
parser.add_argument('-o', '--output', action='store', dest='out_', default=None, help='output file name')
parser.add_argument('-e', '--output_excluded', action='store', dest='out_ex_', default=None, help='output file name that contains excluded reads')
parser.add_argument('-f', '--filter_out', action='store', dest='filter_', default=None, help='file contains id to filter out')
parser.add_argument('-i', '--info_column', action='store', dest='col_info', type=int, default=None, help='column number in filter_out file. information in the column will be added to fasta file header.')
parser.add_argument('-d', '--delimiter', action='store', dest='delim', default=None, help='delimiter for id')
args = parser.parse_args()

if args.filter_ is None:
	sys.stderr.write("File containing id to filter out must be specifiec by -f option.")
	exit

d_filter_out = {}
f_filter_out = open(args.filter_)

# make dictionary of id to filter out
for line in f_filter_out:
	line = line.rstrip()
	a_all = line.split("\t")
	# split id field
	if args.delim is None:
		id = a_all[0]
	else:
		id = (a_all[0]).split(args.delim)[0]

	# id info
	if args.col_info:
		d_filter_out[id] = a_all[args.col_info]
	else:
		d_filter_out[id] = 1

f_filter_out.close()
sys.stderr.write("Number of ids to filter out : %d\n" % len(d_filter_out.keys()))


if args.out_:
	f_out = open(args.out_, "w")
else:
	f_out = sys.stdout

if args.out_ex_:
	f_out_ex = open(args.out_ex_, "w")


count = 0
i = 0
header = ""
sequence = ""
for line in args.f_in:
	line = line.rstrip()
	i += 1

	# check if the array start with id line
	if len(line) > 1:
		if line[0] == ">":
			# check previous header
			if len(header) > 1:
				# split id field
				if args.delim is None:
					id = header[1:]
				else:
					id = (header[1:]).split(args.delim)[0]

				if d_filter_out.has_key(id):
					count += 1
					if args.out_ex_:
						if args.col_info:
							# add information to fasta file header
							f_out_ex.write("%s %s\n%s\n" % (header, d_filter_out[id], sequence))
						else:
							f_out_ex.write("%s\n%s\n" % (header, sequence))
				else:
					f_out.write("%s\n%s\n" % (header, sequence))

			# store current header
			header = line
			# initialize sequence
			sequence = ""
		else:
			# store sequences
			sequence += line

#--------------------
# for the final fasta
# split id field
if args.delim is None:
	id = header[1:]
else:
	id = (header[1:]).split(args.delim)[0]

if d_filter_out.has_key(id):
	count += 1
	if args.out_ex_:
		if args.col_info:
			# add information to fasta file header
			f_out_ex.write("%s %s\n%s\n" % (header, d_filter_out[id], sequence))
		else:
			f_out_ex.write("%s\n%s\n" % (header, sequence))
else:
	f_out.write("%s\n%s\n" % (header, sequence))

# close file handlers
if args.out_:
	f_out.close()

if args.out_ex_:
	f_out_ex.close()

sys.stderr.write("Number of ids filtered out : %d\n" % count)

