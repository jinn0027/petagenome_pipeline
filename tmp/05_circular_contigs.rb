
#--------------------------------------------------
# Program: 05_circular_contigs.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'n_threads:6', 'assembly_label:SPAdes_meta_c', 'l:', 'q:')
n_threads = h_params["n_threads"] # the number of threads used in BLASTN search
assembly_label = h_params["assembly_label"]
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
profile_ = "#{dir_home_}/pipeline/profile"

"---------------------"
blast_dir = ""
"---------------------"

# E-value cutoff for circular formation
e_thre = "1e-10"
# % identity for circular formation
pi_thre = "100"
# minimum alignment length for circular formation
al_thre = "50"

# % identity for removal of redundancy
pi_thre_rd = "95"
# % query coverage for removal of redundancy
qc_thre_rd = "95"

# minimum length for linear contigs
len_l = "5000"
# minimum length for circular contigs
len_c = "1500"

# loading setting file
inifile = IniFile.load('./setting.cfg')
abort "Abort: './setting.cfg' does not exist." if ! File.exist?('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	type = inifile[section]['type']
	dir_project_ = inifile[section]['dir_project_']
	merge_label = inifile[section]['merge_label']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'type' in ./setting.cfg" if type == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil

	# check project type
	if type == "virome" then
		abort "Abort: 'type' is #{type} but '#{$0}' was run." if File.dirname($0) !~ /virome$/
		l_thre = "1000"
	elsif type == "bacteriome" then
		abort "Abort: 'type' is #{type} but '#{$0}' was run." if File.dirname($0) !~ /bacteriome$/
		l_thre = "5000"
	else
		abort "Abort: 'type' in ./setting.cfg must be 'virome' or 'bacteriome'."
	end

	query_label = "#{merge_label}_#{assembly_label}"
	# input directory
	dir_in_ = "#{dir_project_}/d04_pool_contigs/#{merge_label}"
	# output directory
	dir_out_ = "#{dir_project_}/d05_circular_contigs/#{merge_label}"
	# log directory
	dir_log_ = "#{dir_project_}/d05_circular_contigs/log"

	# input
	query_fa_ = "#{dir_in_}/contig.#{query_label}.#{l_thre}.fa"

	# output
	out = "#{dir_out_}/contig.updated.#{query_label}.#{len_l}.c_#{len_c}"

	# check
	abort "Abort: 'query_fa_':#{query_fa_} does not exist." if ! File.exist?(query_fa_)
	abort "Abort: 'output directory':#{dir_out_} already exists." if Dir.exist?(dir_out_)

	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	# move to log directory
	Dir::chdir(dir_log_)

	# qsub file
	puts "start #{merge_label}"
	qsub_ = "p05_circular_contigs.#{project}.#{query_label}.sh"
	f_qsub = open(qsub_, "w")
	f_qsub.puts "\#!/bin/bash"
	f_qsub.puts "\#$ -S /bin/bash"
	f_qsub.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G -pe def_slot #{n_threads}"
	f_qsub.puts "\#$ -cwd"
	f_qsub.puts "\# version #{version}"
	f_qsub.puts "source #{profile_}"
	f_qsub.puts "SCRIPT_SELF_=\"#{dir_scripts_}/run_blastn.self.contig.3.rb\""
	f_qsub.puts "MAKEBLASTDB_=\"#{blast_dir}/bin/makeblastdb\""

	# main script for exploring circular contigs
	f_qsub.puts "ruby ${SCRIPT_SELF_}  -d #{dir_out_}  -i #{query_fa_}  -n #{n_threads}  --len_l #{len_l} --len_c #{len_c} --pi_rd #{pi_thre_rd} --qc_rd #{qc_thre_rd} "
	# copy results to adjust to the existing scripts
	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.fa  #{out}.fa"
	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.length.txt  #{out}.length.txt"
	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.extended.fa  #{out}.extended.fa"
	f_qsub.puts "cp  #{dir_out_}/contig.updated.#{len_l}.c_#{len_c}.extended.length.txt  #{out}.extended.length.txt"
	# blastdb
	f_qsub.puts "${MAKEBLASTDB_} -in #{out}.fa -out #{out} -dbtype nucl -parse_seqids"

	# name information
	# updated (final) contigs
	f_qsub.puts "cut -f 1 #{out}.length.txt > #{out}.txt"
	f_qsub.puts "cat #{out}.txt | perl -pe \"s/l(\\.\\d+)$/n\\1/\" | perl -pe \"s/c(\\.\\d+)$/n\\1/\" | paste - #{out}.txt > #{dir_out_}/name.updated.#{query_label}.#{len_l}.c_#{len_c}.txt"
	# redundant (excluded) contigs
	f_qsub.puts "cut -f 1,6 #{dir_out_}/contig.redundant.c_#{len_c}.length.txt | perl -pe \"s/l(\\.\\d+)\\t/n\\1\\t/\" | perl -pe \"s/c(\\.\\d+)\\t/n\\1\\t/\" > #{dir_out_}/name.redundant.#{query_label}.#{len_l}.c_#{len_c}.txt"
	# all contigs
	f_qsub.puts "cat #{dir_out_}/name.updated.#{query_label}.#{len_l}.c_#{len_c}.txt #{dir_out_}/name.redundant.#{query_label}.#{len_l}.c_#{len_c}.txt > #{dir_out_}/name.all.#{query_label}.#{len_l}.c_#{len_c}.txt"

	f_qsub.close
	if q_opt != nil then
		system("qsub -q #{q_opt} #{qsub_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} #{qsub_}")
	else
		system("qsub #{qsub_}")
	end
end

