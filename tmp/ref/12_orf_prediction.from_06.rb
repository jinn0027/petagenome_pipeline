
#--------------------------------------------------
# Program: 12_orf_prediction.from_06.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'assembly_label:SPAdes_meta_c', 'l:', 'q:')
assembly_label = h_params["assembly_label"]
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
profile_ = "#{dir_home_}/pipeline/profile"

l_thre = "5000.c_1500"
c_thre = "1500"

"---------------------"
prodigal_dir = ""
"---------------------"

predictor = "prodigal"

# loading setting file
abort "Abort: './setting.cfg' does not exist." if ! File.exist?('./setting.cfg')
inifile = IniFile.load('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	dir_project_ = inifile[section]['dir_project_']
	merge_label = inifile[section]['merge_label']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil

	# input directories
	dir_05_ = "#{dir_project_}/d05_circular_contigs/#{merge_label}"
	dir_06_ = "#{dir_project_}/d06_decontamination/#{merge_label}"
	# output directory
	dir_out_ = "#{dir_project_}/d12_orf_prediction/#{merge_label}"
	# log directory
	dir_log_ = "#{dir_project_}/d12_orf_prediction/log"

	query_label = "#{merge_label}_#{assembly_label}"
	# input
	contig_c_len_ = "#{dir_05_}/contig.circular.c_#{c_thre}.length.txt"
	contig_ex_fa_ = "#{dir_06_}/contig.decontaminated.#{query_label}.#{l_thre}.extended.fa"

	# output
	out_ex = "#{dir_out_}/#{predictor}.#{query_label}.#{l_thre}.extended.predicted"
	out = "#{dir_out_}/#{predictor}.#{query_label}.#{l_thre}.predicted"

	# check
	abort "Abort: 'contig_c_len_':#{contig_c_len_} does not exist." if ! File.exist?(contig_c_len_)
	abort "Abort: 'contig_ex_fa_':#{contig_ex_fa_} does not exist." if ! File.exist?(contig_ex_fa_)
	if Dir.exist?(dir_out_) then
		abort "Abort: 'output directory':#{dir_out_} already exists."
	else
		FileUtils.mkdir_p(dir_out_)
	end

	# move to log directory
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	Dir.chdir(dir_log_)

	# qsub file
	puts "start #{merge_label}"
	qsub_ = "p12_#{predictor}.#{project}.#{query_label}.#{l_thre}.sh"
	f_qsub = open(qsub_, "w")
	f_qsub.puts "\#!/bin/bash"
	f_qsub.puts "\#$ -S /bin/bash"
	f_qsub.puts "\#$ -l s_vmem=16G,mem_req=16G"
	f_qsub.puts "\#$ -cwd"
	f_qsub.puts "\# version #{version}"
	f_qsub.puts "source #{profile_}"
	f_qsub.puts "PRODIGAL_=\"#{prodigal_dir}/prodigal\""
	f_qsub.puts "SCRIPT_JOIN_=#{dir_scripts_}/join_with_tab.rb"
	f_qsub.puts "SCRIPT_FILTER_FA_=\"#{dir_scripts_}/filter_fasta_by_id.py\""
	f_qsub.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""

	# run prodigal
	f_qsub.puts "${PRODIGAL_} -a #{out_ex}.original.faa  -d #{out_ex}.original.ffn  -p meta  -i #{contig_ex_fa_} -f gbk  -o #{out_ex}.original.gbk"
	# convert ORF name
	f_qsub.puts "cat #{out_ex}.original.faa | perl -pe \"s/_\\d+ # /__/\" | perl -pe \"s/ # /,/\" | perl -pe \"s/ # /_/\" | perl -pe \"s/ # .+$//\" > #{out_ex}.faa"
	f_qsub.puts "cat #{out_ex}.original.ffn | perl -pe \"s/_\\d+ # /__/\" | perl -pe \"s/ # /,/\" | perl -pe \"s/ # /_/\" | perl -pe \"s/ # .+$//\" > #{out_ex}.ffn"
	# length of ORFs
	f_qsub.puts "python ${SCRIPT_LENGTH_} #{out_ex}.ffn > #{out_ex}.ffn.length.txt"
	# position of ORFs
	f_qsub.puts "grep \">\" #{out_ex}.original.ffn | perl -pe \"s/_\\d+ # /\\t/\" | perl -pe \"s/ # /\\t/\" | perl -pe \"s/ # /\\t/\" | perl -pe \"s/ # .+(partial=\\d{2};).+$/\\t\\1/\" | perl -pe \"s/^>//\" | paste #{out_ex}.ffn.length.txt - > #{out_ex}.ffn.position.txt"

	# excluding ORFs in redundant region
	f_qsub.puts "ruby ${SCRIPT_JOIN_} #{contig_c_len_} 1 #{out_ex}.ffn.position.txt 3 | awk -F \"\\t\" '{OFS=\"\\t\"}  {if($8 == \"partial=10;\" || $5 > $2) print $3,$4,$1,$2,$5,$6,$7,$8}' > #{out_ex}.exclude.txt"
	f_qsub.puts "python ${SCRIPT_FILTER_FA_}  -f #{out_ex}.exclude.txt  -o #{out}.ffn  #{out_ex}.ffn"
	f_qsub.puts "python ${SCRIPT_FILTER_FA_}  -f #{out_ex}.exclude.txt  -o #{out}.faa  #{out_ex}.faa"

	# length of ORFs
	f_qsub.puts "python ${SCRIPT_LENGTH_} #{out}.ffn > #{out}.ffn.length.txt"
	f_qsub.puts "python ${SCRIPT_LENGTH_} #{out}.faa > #{out}.faa.length.txt"
	# add contig id
	f_qsub.puts "perl -pe \"s/^(.+\\.\\d{9})__/\\1\\t\\1__/\" #{out}.ffn.length.txt | cut -f 1 | paste #{out}.ffn.length.txt -  > #{out}.ffn.length.contig.txt"
	f_qsub.puts "perl -pe \"s/^(.+\\.\\d{9})__/\\1\\t\\1__/\" #{out}.faa.length.txt | cut -f 1 | paste #{out}.faa.length.txt -  > #{out}.faa.length.contig.txt"

	f_qsub.close
	if q_opt != nil then
		system("qsub -q #{q_opt} #{qsub_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} #{qsub_}")
	else
		system("qsub #{qsub_}")
	end
end
