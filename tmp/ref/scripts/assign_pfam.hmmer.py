#!/usr/bin/env python

import sys
import glob
import os
import re


def read_pfam_info(pfam_info_):
	import gzip
	d_pfam_info = {}
	f_pfam_info = open(pfam_info_)

	# make dictionary of Pfam id to general function
	for line in f_pfam_info:
		a_all = line[:-1].split("\t")
		pfamid = a_all[0]
		info = a_all[1]
		d_pfam_info[pfamid] = info

	f_pfam_info.close()
	return d_pfam_info


def run(in_, thre_score, thre_eval, f_out, pfam_func_, pfam_clans_, pfam_desc_, dir_data_pfam_):
	f_in = open(in_, "r")
	# Pfam dirctory
	if dir_data_pfam_ == None:
		dir_data_pfam_ = "/home/user/data/Pfam"

	# read Pfam information
	if pfam_func_ == None:
		pfam_func_ = "%s/Pfam-A.general_function.txt" % dir_data_pfam_
	d_pfam_func = read_pfam_info(pfam_func_)

	if pfam_clans_ == None:
		pfam_clans_ = "%s/Pfam-A.clans.txt" % dir_data_pfam_
	d_pfam_clans = read_pfam_info(pfam_clans_)

	if pfam_desc_ == None:
		pfam_desc_ = "%s/Pfam-A.description.txt" % dir_data_pfam_
	d_pfam_desc = read_pfam_info(pfam_desc_)

	# read input file and filter by threshold
	p = re.compile(r'\s+')
	for line in f_in:
		try:
			if line[0] == "#":
				continue
			line = line[:-1]
			a_all = p.split(line)
			target_name = a_all[0]
			pfamid = a_all[1].split(".")[0]
			query_name = a_all[2]
			eval = a_all[4]
			score = a_all[5]
			bias = a_all[6]

			if float(eval) > float(thre_eval):
				continue

			if float(score) < thre_score:
				continue

			if d_pfam_func.has_key(pfamid):
				pfam_func = d_pfam_func[pfamid]
			else:
				pfam_func = "Unassigned"

			if d_pfam_clans.has_key(pfamid):
				pfam_clans = d_pfam_clans[pfamid]
			else:
				pfam_clans = "Unassigned"

			if d_pfam_desc.has_key(pfamid):
				pfam_desc = d_pfam_desc[pfamid]
			else:
				pfam_desc = "Unassigned"

			# swap 1st and 3rd elements
			a_all[0] = query_name
			a_all[2] = target_name
			f_out.write("%s\t%s\t%s\t%s\t%s\n" % ("\t".join(a_all[0:18]), pfamid, pfam_func, pfam_clans, pfam_desc))

		except Exception, e:
			sys.stderr.write("Error occurred: %s for %s\n" % (e, line))

	f_in.close()


def main(argv):
	import argparse # to use this module python version must be higher than 2.7
	# parse input arguments
	parser = argparse.ArgumentParser(description='This script is to parse hmmscan results of Resfams. \
Copyright (C) 2015 Yasumasa Kimura, Division of Systems Immunology, Institute of Medical Science, The University of Tokyo, Tokyo, Japan.')
	parser.add_argument('input_', help='input tab-delimited hmmscan result data')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 2016.2.5')
	parser.add_argument('-o', '--output', action='store', dest='out_', help='output file name', default=None)
	parser.add_argument('-p', '--pfam_dir', action='store', dest='pfam_dir_', help='pfam directory', default=None)
	parser.add_argument('-f', '--func', action='store', dest='pfam_func_', help='pfam information to assign', default=None)
	parser.add_argument('-c', '--clans', action='store', dest='pfam_clans_', help='pfam clans to assign', default=None)
	parser.add_argument('-d', '--desc', action='store', dest='pfam_desc_', help='pfam description to assign', default=None)
	parser.add_argument('-t', '--thre_bitScore', action='store', dest='thre_score', type=int, default=0, help='threshold for threshold for bitScore. if the bitScore is higher than this value the read will be counted.')
	parser.add_argument('-e', '--thre_eVal', action='store', dest='thre_eval', type=float, default=1.0, help='threshold for threshold for eVal. if the eVal is lower than this value the read will be counted.')
	args = parser.parse_args()

	# output file
	if args.out_:
		f_out = open(args.out_, 'w')
	else:
		f_out = sys.stdout

	run(args.input_, args.thre_score, args.thre_eval, f_out, args.pfam_func_, args.pfam_clans_, args.pfam_desc_, args.pfam_dir_)
	if args.out_:
		f_out.close()

if __name__ == "__main__":
	main(sys.argv[1:])

