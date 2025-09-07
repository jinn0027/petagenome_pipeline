#!/usr/bin/env python

import sys

def run(f_in, f_out, f_out_ko, f_out_pathway, f_out_module, db_date):
	d_gene2ko = {}
	d_ko2pathway = {}
	d_ko2module = {}
	#if db_date is None:
	ko_genes_ = "/home/user/kegg.20180916/genes/ko/ko_genes.list"
	ko_pathway_ = "/home/user/kegg.20180916/genes/ko/ko_pathway.list"
	ko_module_ = "/home/user/kegg.20180916/genes/ko/ko_module.list"

	# make dictionary of gene to ko
	f_ko_genes = open(ko_genes_)
	for line in f_ko_genes:
		a_all = line[:-1].split("\t")
		_ko = a_all[0]
		gene = a_all[1]
		if _ko[0:3] != "ko:":
			sys.stderr.write("KO in different format in ko_genes %s" % _ko)
			continue
		ko = _ko[3:]
		# store gene - ko information
		# some genes have more than one KO
		if d_gene2ko.has_key(gene):
			d_ko = d_gene2ko[gene]
			d_ko[ko] = 1
		else:
			d_gene2ko[gene] = {ko:1}
	f_ko_genes.close()

	# make dictionary of ko to pathway
	f_ko_pathway = open(ko_pathway_)
	for line in f_ko_pathway:
		a_all = line[:-1].split("\t")
		_ko = a_all[0]
		_pathway = a_all[1]
		if _ko[0:3] != "ko:":
			sys.stderr.write("KO in different format in ko_pathway %s" % _ko)
			continue
		if _pathway[0:8] != "path:map":
			continue
		ko = _ko[3:]
		pathway = _pathway[8:]
		# store ko - pathway information
		# KO involves in several pathways
		if d_ko2pathway.has_key(ko):
			d_pathway = d_ko2pathway[ko]
			d_pathway[pathway] = 1
		else:
			d_ko2pathway[ko] = {pathway:1}
	f_ko_pathway.close()

	# make dictionary of ko to module
	f_ko_module = open(ko_module_)
	for line in f_ko_module:
		a_all = line[:-1].split("\t")
		_ko = a_all[0]
		_module = a_all[1]
		if _ko[0:3] != "ko:":
			sys.stderr.write("KO in different format in ko_module %s" % _ko)
			continue
		if _module[0:3] != "md:":
			continue
		ko = _ko[3:]
		module = _module[3:]
		# store ko - module information
		# KO involves in several modules
		if d_ko2module.has_key(ko):
			d_module = d_ko2module[ko]
			d_module[module] = 1
		else:
			d_ko2module[ko] = {module:1}
	f_ko_module.close()

	# read input file and put KEGG annotation
	for line in f_in:
		a_all = line[:-1].split("\t")
		query = a_all[0]
		subj = a_all[1]
		gene = subj.split(" ")[0]

		out_kos = ""
		out_pathways = ""
		out_modules = ""
		d_ko = d_gene2ko.get(gene)
		if d_ko is not None:
			a_out_ko = []
			d_out_pathway = {}
			d_out_module = {}
			for ko in d_ko.keys():
				if ko == "":
					continue
				a_out_ko.append(ko)
				# related pathway
				d_pathway = d_ko2pathway.get(ko)
				if d_pathway is not None:
					for pathway in d_pathway.keys():
						d_out_pathway[pathway] = 1
					# output pathway
					for pathway in d_out_pathway.keys():
						f_out_pathway.write("%s\t%s\t%s\n" % (pathway, ko, query))

				# related module
				d_module = d_ko2module.get(ko)
				if d_module is not None:
					for module in d_module.keys():
						d_out_module[module] = 1
					# output module
					for module in d_out_module.keys():
						f_out_module.write("%s\t%s\t%s\n" % (module, ko, query))


			if len(a_out_ko) > 0:
				out_kos = ",".join(a_out_ko)
				# output ko
				for ko in a_out_ko:
					f_out_ko.write("%s\t%s\n" % (ko, query))

			if len(d_out_pathway) > 0:
				out_pathways = ",".join(d_out_pathway.keys())

			if len(d_out_module) > 0:
				out_modules = ",".join(d_out_module.keys())

		# output gene with KEGG annotation
		f_out.write("%s\t%s\t%s\t%s\n" % (line[:-1], out_kos, out_pathways, out_modules))



import argparse # to use this module python version must be higher than 2.7
# parse input arguments
parser = argparse.ArgumentParser(description='This script is to put KEGG KO and pathway IDs. \
Copyright (C) 2018 Yasumasa Kimura, Division of Systems Immunology, Institute of Medical Science, The University of Tokyo, Tokyo, Japan.')
parser.add_argument('f_in', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input blast-tab-delimited file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 2018.06.18')
parser.add_argument('-d', '--db_date', action='store', dest='db_date', help='date of KEGG data', default=None)
parser.add_argument('-o', '--output', action='store', dest='out', help='output file basename', default=None)
args = parser.parse_args()

# output file
if args.out:
	f_out = open("%s.info.txt" % args.out, 'w')
	f_out_ko = open("%s.ko.txt" % args.out, 'w')
	f_out_pathway = open("%s.pathway.txt" % args.out, 'w')
	f_out_module = open("%s.module.txt" % args.out, 'w')
else:
	sys.stderr.write("output file basename is required using option -o")
	sys.exit()

# main
run(args.f_in, f_out, f_out_ko, f_out_pathway, f_out_module, args.db_date)
f_out.close()
f_out_ko.close()
f_out_pathway.close()
f_out_module.close()


