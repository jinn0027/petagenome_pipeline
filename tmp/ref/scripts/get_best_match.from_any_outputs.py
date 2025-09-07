#!/usr/bin/env python

import sys
import glob
import os
import re

def run(in_, f_out, args):
	d_best_in ={}

	for in_ in glob.glob(in_): # wild card can be used for the file name
		f_in = open(in_, "r")

		for line in f_in:
			a_all = line[:-1].split("\t")
			queryId = a_all[args.col_id]
			score = float(a_all[args.col_score])

			if d_best_in.has_key(queryId):
				# best hit so far
				best_score = (d_best_in[queryId])['score']
				# update the best hit if current data has higher score
				if best_score < score:
					d_best_in[queryId] = {'line':line, 'score':score}
			else:
				# initalize with first data
				d_best_in[queryId] = {'line':line, 'score':score}

		f_in.close()

	# output
	for d_value in d_best_in.values():
		f_out.write(d_value['line'])


def main(argv):
	import argparse # to use this module python version must be higher than 2.7
	# parse input arguments
	parser = argparse.ArgumentParser(description='This script is to obtain best aligned hit from multiple blast tabular output files mapped to different referece sequences. \
Copyright (C) 2015 Yasumasa Kimura, Division of Systems Immunology, Institute of Medical Science, The University of Tokyo, Tokyo, Japan.',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_', help='input file (wild card can be used for the file name)')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 2016.07.25')
	parser.add_argument('-o', '--output', action='store', dest='out_', help='output file name', default=None)
	parser.add_argument('-i', '--col_id', action='store', dest='col_id', type=int, default=0, help='column that is summarised for.')
	parser.add_argument('-s', '--col_score', action='store', dest='col_score', type=int, default=None, help='column that has scores to evaluate.')
	args = parser.parse_args()

	if args.col_score is None:
		sys.stderr.write("Column of scores must be specifiec by -s or --col_score option.")
		exit

	# open output files
	if args.out_:
		f_out  = open(args.out_, 'w')
	else:
		f_out = sys.stdout

	# run
	run(args.input_, f_out, args)
	if args.out_:
		f_out.close()

if __name__ == "__main__":
	main(sys.argv[1:])

