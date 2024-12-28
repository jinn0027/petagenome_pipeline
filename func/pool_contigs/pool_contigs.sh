#!/bin/bash

exit


#--------------------------------------------------
# Program: 04_pool_contigs.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'mem:128', 'l:', 'q:', 'pre_assembly_label:SPAdes_meta', 'assembly_label:SPAdes_meta_c')
mem = h_params["mem"] # memory limit
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q
pre_assembly_label = h_params["pre_assembly_label"]
assembly_label = h_params["assembly_label"]

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
#dir_R_ = "/usr/local/package/r/3.4.1/bin"
profile_ = "#{dir_home_}/pipeline/profile"

"---------------------"
blast_dir = ""
cdhit_dir = ""
"---------------------"

# memory usage
if mem == "16" then
	s_vmem = "5.3"
	n_threads = "3"
elsif mem == "32" then
	s_vmem = "5.3"
	n_threads = "6"
elsif mem == "64" then
	s_vmem = "5.3"
	n_threads = "12"
elsif mem == "128"
	s_vmem = "5.3"
	n_threads = "24"
elsif mem == "256"
	# this setting will use lmem.q
	s_vmem = "41.7"
	n_threads = "6"
else
	abort "please specify --mem as 16, 32, 64, 128 or 256."
end

# loading setting file
abort "Abort: ./setting.cfg does not exist." if ! File.exist?('./setting.cfg')
inifile = IniFile.load('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	type = inifile[section]['type']
	dir_project_ = inifile[section]['dir_project_']
	samples = inifile[section]['samples']
	merge_label = inifile[section]['merge_label']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'type' in ./setting.cfg" if type == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	a_samples = samples.gsub(" ","").split(',')

	# check project type
	if type == "virome" then
		prefix = "vPContig"
		abort "Abort: 'type' is #{type} but '#{$0}' was run." if File.dirname($0) !~ /virome$/
		l_thre = "1000"
	elsif type == "bacteriome" then
		prefix = "bPContig"
		abort "Abort: 'type' is #{type} but '#{$0}' was run." if File.dirname($0) !~ /bacteriome$/
		l_thre = "5000"
	else
		abort "Abort: 'type' in ./setting.cfg must be 'virome' or 'bacteriome'."
	end

	# check input fasta files
	h_samples = {}
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		dir_in_ = "#{dir_project_}/d03_assembly/#{sample}"
		input_fa_ = "#{dir_in_}/contig.#{sample}_#{pre_assembly_label}.#{l_thre}.fa"
		if File.exist?(input_fa_) then
			abort "Abort: duplicated sample name #{sample}" if h_samples.key?(sample)
			h_samples[sample] = input_fa_
			puts "Check OK: 'input_fa_':#{input_fa_} exists."
		else
			abort "Abort: 'input_fa_':#{input_fa_} does not exist."
		end
	end
	puts "#samples: #{h_samples.length}"

	# output directory
	dir_out_ = "#{dir_project_}/d04_pool_contigs/#{merge_label}"
	if Dir.exist?(dir_out_) then
		abort "Abort: directory for merge_label '#{merge_label}':#{dir_out_} already exists."
	else
		FileUtils.mkdir_p(dir_out_)
	end

	# log directory
	dir_log_ = "#{dir_project_}/d04_pool_contigs/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	# make a copy of 'setting.cfg' to keep a record
	system("cp #{dir_project_}/setting.cfg #{dir_project_}/setting.#{merge_label}.cfg")

	# output
	out = "#{dir_out_}/contig.#{merge_label}_#{assembly_label}"
	merged_fa_ = "#{dir_out_}/merged_contig.#{pre_assembly_label}.#{merge_label}.#{l_thre}.fa"
	rep_fa_ = "#{dir_out_}/rep_contig.#{pre_assembly_label}.#{merge_label}.#{l_thre}.fa"
	out_table_name_ = "#{dir_out_}/name.#{merge_label}_#{assembly_label}.#{l_thre}.txt"

	#------------------------------
	puts "start #{merge_label}"
	qsub_ = "p04_pool_contigs.#{project}.#{merge_label}_#{assembly_label}.sh"
	f_qsub = open(qsub_, "w")
	f_qsub.puts "\#!/bin/bash"
	f_qsub.puts "\#$ -S /bin/bash"
	f_qsub.puts "\#$ -l s_vmem=#{s_vmem}G,mem_req=#{s_vmem}G -pe def_slot #{n_threads}"
	f_qsub.puts "\#$ -cwd"
	f_qsub.puts "\# version #{version}"
	f_qsub.puts "source #{profile_}"
	f_qsub.puts "CDHIT_=\"#{cdhit_dir}/cd-hit-est\""
	f_qsub.puts "SCRIPT_RENAME_=\"#{dir_scripts_}/filter_contig.rename.py\""
	f_qsub.puts "SCRIPT_CLSTR_=\"#{dir_scripts_}/parse.cdhit_clstr.rb\""
	f_qsub.puts "SCRIPT_JOIN_=\"#{dir_scripts_}/join_with_tab.rb\""
	f_qsub.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""
	f_qsub.puts "SCRIPT_STATS_=\"#{dir_scripts_}/stats.assembly.R\""
	f_qsub.puts "MAKEBLASTDB_=\"#{blast_dir}/bin/makeblastdb\""

	#------------------------------
	# concatenate assemblies
	f_qsub.puts "echo \"start concatenation\""
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		dir_in_ = "#{dir_project_}/d03_assembly/#{sample}"
		input_fa_ = "#{dir_in_}/contig.#{sample}_#{pre_assembly_label}.#{l_thre}.fa"
		if l == 0 then
			f_qsub.puts "cat #{input_fa_} >  #{merged_fa_}"
		else
			f_qsub.puts "cat #{input_fa_} >> #{merged_fa_}"
		end
	end

	#------------------------------
	# removing redundancy by CD-HIT-EST
	f_qsub.puts "echo \"start cd-hit-est\""
	f_qsub.puts "${CDHIT_} -c 0.95 -G 1 -mask NX -d 150 -n 10 -T #{n_threads} -M #{mem}000 -i #{merged_fa_} -o #{rep_fa_}"

	# rename and filter (L>=#{l_thre}) contigs
	# this contig name format is nessecary for later procedures.
	# the corresponding name table will be output in #{out}.#{l_thre}.name.txt.
	f_qsub.puts "python ${SCRIPT_RENAME_} -min #{l_thre} --rename --prefix #{prefix}.#{project}.#{merge_label}.n. --table #{out}.#{l_thre}.name.txt #{rep_fa_} > #{out}.#{l_thre}.fa"

	# summaize contig name ( contig in each sample / representative contig / renamed representative contig )
	f_qsub.puts "ruby ${SCRIPT_CLSTR_} -i #{rep_fa_}.clstr --include_rep > #{rep_fa_}.name"
	f_qsub.puts "ruby ${SCRIPT_JOIN_} #{rep_fa_}.name 2 #{out}.#{l_thre}.name.txt 2 | awk -F '\\t' '{OFS=\"\\t\"} {print $2,$1,$3}' > #{out_table_name_}"
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		out_table_name_each_ = "#{dir_out_}/name.#{sample}_#{pre_assembly_label}.#{merge_label}.#{l_thre}.txt"
		f_qsub.puts "awk -F \"\\t\" '{OFS=\"\\t\"} { if ($1 ~ /#{sample}/) print $0 }' #{out_table_name_} > #{out_table_name_each_}"
	end

	# get length of contigs
	f_qsub.puts "python ${SCRIPT_LENGTH_} #{out}.#{l_thre}.fa > #{out}.#{l_thre}.length.txt"
	# stats of assemblies
	f_qsub.puts "R --vanilla --slave --args #{out}.#{l_thre}.length.txt 2 < ${SCRIPT_STATS_} > #{out}.#{l_thre}.stats.txt"
	# blastdb
	f_qsub.puts "${MAKEBLASTDB_} -in #{out}.#{l_thre}.fa -out #{out}.#{l_thre} -dbtype nucl -parse_seqids"

	f_qsub.close
	if q_opt != nil then
		system("qsub -q #{q_opt} #{qsub_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} #{qsub_}")
	else
		system("qsub #{qsub_}")
	end
end
