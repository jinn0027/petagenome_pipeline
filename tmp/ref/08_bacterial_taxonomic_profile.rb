
#--------------------------------------------------
# Program: 08_bacterial_taxonomic_profile.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'n_threads:6', 'l:', 'q:')
n_threads = h_params["n_threads"]
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
profile_ = "#{dir_home_}/pipeline/profile"

"-----------------"
metaphlan_dir = ""
"-----------------"

db = "./metaphlan_db_v20/mpa_v20_m200"
db_label = "metaphlan2_v20"


# load setting file
abort "Abort: './setting.cfg' does not exist." if ! File.exist?('./setting.cfg')
inifile = IniFile.load('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	dir_project_ = inifile[section]['dir_project_']
	samples = inifile[section]['samples']
	merge_label = inifile[section]['merge_label']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	a_samples = samples.gsub(" ","").split(',')

	# output directory
	dir_out_base_ = "#{dir_project_}/d08_bacterial_taxonomic_profile"
	FileUtils.mkdir_p(dir_out_base_) unless Dir.exist?(dir_out_base_)
	# log directory
	dir_log_ = "#{dir_out_base_}/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	a_jid = [] # job names to be finished before concatenation
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		# input
		dir_in_ = "#{dir_project_}/d02_error_correction/#{sample}"
		query1_ = "#{dir_in_}/corrected/cutadapt_prinseq.#{sample}_1.00.0_0.cor.fastq.gz"
		query2_ = "#{dir_in_}/corrected/cutadapt_prinseq.#{sample}_2.00.0_0.cor.fastq.gz"

		# output
		dir_out_ = "#{dir_out_base_}/#{sample}"
		out_bwt2_ = "#{dir_out_}/bowtie2.#{sample}.#{db_label}.sam.bz2"
		out_profile_ = "#{dir_out_}/profile.#{sample}.#{db_label}.txt"
		out_parsed_ = "#{dir_out_}/profile.#{sample}.#{db_label}.parsed.txt"

		# check
		abort "Abort: 'query1_':#{query1_} does not exist." if ! File.exist?(query1_)
		abort "Abort: 'query2_':#{query2_} does not exist." if ! File.exist?(query2_)
		if Dir.exist?(dir_out_) then
			abort "Abort: 'output directory':#{dir_out_} already exists."
		else
			FileUtils.mkdir_p(dir_out_)
		end

		# qsub file
		puts "start #{sample}"
		qsub_ = "p08_metaphlan2.#{project}.#{sample}.sh"
		f_qsub = open(qsub_, "w")
		f_qsub.puts "\#!/bin/bash"
		f_qsub.puts "\#$ -S /bin/bash"
		f_qsub.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G -pe def_slot #{n_threads}"
		f_qsub.puts "\#$ -cwd"
		f_qsub.puts "\# version #{version}"
		f_qsub.puts "source #{profile_}"
		f_qsub.puts "METAPHLAN2_=\"#{metaphlan_dir}/metaphlan2.py\""
		f_qsub.puts "SCRIPT_PARSE_=\"#{dir_scripts_}/parse.mpa-report.py\""
		f_qsub.puts "DB=\"#{db}\""
		f_qsub.puts "OUT_BWT2_=\"#{out_bwt2_}\""
		f_qsub.puts "OUT_PROFILE_=\"#{out_profile_}\""
		f_qsub.puts "OUT_PARSED_=\"#{out_parsed_}\""
		# run metaphlan
		f_qsub.puts "zcat #{query1_} #{query2_} | ${METAPHLAN2_} --bowtie2db ${DB} --input_type fastq --bowtie2out ${OUT_BWT2_} --nproc #{n_threads} > ${OUT_PROFILE_}"
		# parse profile
		f_qsub.puts "python ${SCRIPT_PARSE_}  -o ${OUT_PARSED_}  ${OUT_PROFILE_}"
		f_qsub.close
		if q_opt != nil then
			system("qsub -q #{q_opt} #{qsub_}")
		elsif l_opt != nil then
			system("qsub -l #{l_opt} #{qsub_}")
		else
			system("qsub #{qsub_}")
		end
		a_jid.push("#{qsub_}")
	end

	#------------------------------
	# join profiles
	#------------------------------
	# output
	dir_out_merge_ = "#{dir_out_base_}/#{merge_label}"
	FileUtils.mkdir_p(dir_out_merge_) unless Dir.exist?(dir_out_merge_)
	out_merged_profile_ = "#{dir_out_merge_}/profile.#{merge_label}.#{db_label}.txt"
	out_merged_parsed_ = "#{dir_out_merge_}/profile.#{merge_label}.#{db_label}.parsed.txt"

	# temporary files
	temp_merged_ = "#{dir_out_merge_}/temp.merged.coverage.txt"
	header_ = "#{dir_out_merge_}/temp.header.txt"

	# qsub for join
	qsub_join_ = "p08_join_metaphlan2.#{project}.#{merge_label}.sh"
	f_qsub_join = open(qsub_join_, "w")
	f_qsub_join.puts "\#!/bin/bash"
	f_qsub_join.puts "\#$ -S /bin/bash"
	f_qsub_join.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G"
	f_qsub_join.puts "\#$ -cwd"
	f_qsub_join.puts "\# version #{version}"
	f_qsub_join.puts "SCRIPT_JOIN_=#{dir_scripts_}/join_with_tab.option.rb"
	f_qsub_join.puts "SCRIPT_PARSE_=\"#{dir_scripts_}/parse.mpa-report.py\""
	f_qsub_join.puts "echo \"taxonomy\t#{a_samples.join("\t")}\" > #{header_}"
	# initial output file
	f_qsub_join.puts "if [ -f #{temp_merged_} ]; then\n  rm #{temp_merged_} \nfi"
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		in_profile_ = "#{dir_out_base_}/#{sample}/profile.#{sample}.#{db_label}.txt"
		# concatenate all results to get uniq taxonomy names
		f_qsub_join.puts "grep -v \"^#\" #{in_profile_} >> #{temp_merged_}"
	end
	# get all uniq taxonomy names
	f_qsub_join.puts "cut -f 1 #{temp_merged_} | sort | uniq > #{out_merged_profile_}"

	# join profile in each sample
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		# initial output file
		in_profile_ = "#{dir_out_base_}/#{sample}/profile.#{sample}.#{db_label}.txt"
		# add coverages of current sample to final output file
		f_qsub_join.puts "ruby ${SCRIPT_JOIN_}  #{out_merged_profile_} 1 #{in_profile_} 1  > #{temp_merged_}"
		# concatenate missing contig as count 0
		f_qsub_join.puts "ruby ${SCRIPT_JOIN_}  #{out_merged_profile_} 1 #{in_profile_} 1 \"-v 1\" | grep -v \"^#\" | sed -e \"s/$/\t0/\" | cat  >> #{temp_merged_}"
		# sort
		f_qsub_join.puts "sort -k1,1 #{temp_merged_} > #{out_merged_profile_}"
	end
	# final result
	f_qsub_join.puts "sort -k1,1 #{temp_merged_} | cat #{header_} - > #{out_merged_profile_}"
	f_qsub_join.puts "wait; rm #{temp_merged_} #{header_}"
	# parse profile
	f_qsub_join.puts "python ${SCRIPT_PARSE_}  -o #{out_merged_parsed_} #{out_merged_profile_}"
	f_qsub_join.close
	# qsub a job that will run after all the jobs above

	f_qsub_join.close
	if q_opt != nil then
		system("qsub -q #{q_opt} -hold_jid #{a_jid.join(",")} #{qsub_join_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} -hold_jid #{a_jid.join(",")} #{qsub_join_}")
	else
		system("qsub -hold_jid #{a_jid.join(",")} #{qsub_join_}")
	end
end

