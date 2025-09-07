
#--------------------------------------------------
# Program: 16_orf_profile.01_coverage.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"

# contig info
assembly_label = "SPAdes_meta_c"
l_thre = "5000.c_1500"
task = "bbmap"
predictor = "prodigal"

# scripts
script_join_ = "#{dir_scripts_}/join_with_tab.option.rb"
script_normalize_ ="#{dir_scripts_}/normalize.contig_profile_table.R"

# loading setting file
abort "Abort: ./setting.cfg does not exist." if ! File.exist?('./setting.cfg')
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
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')

	puts "start #{merge_label}"
	db_label = "#{merge_label}_#{assembly_label}"
	# input directories
	dir_12_ = "#{dir_project_}/d12_orf_prediction/#{merge_label}"
	dir_10_ = "#{dir_project_}/d10_contig_coverage/#{merge_label}"

	# output directory
	dir_out_ = "#{dir_project_}/d16_orf_profile_table/#{merge_label}"
	# log directory
	dir_log_ = "#{dir_project_}/d16_orf_profile_table/log"

	# input
	in_orf_len_ = "#{dir_12_}/#{predictor}.#{db_label}.#{l_thre}.extended.predicted.ffn.length.txt"

	# check
	abort "Abort: 'in_orf_len_':#{in_orf_len_} does not exist." if ! File.exist?(in_orf_len_)
	if Dir.exist?(dir_out_) then
		abort "Abort: 'output directory':#{dir_out_} already exists."
	else
		FileUtils.mkdir_p(dir_out_)
	end

	FileUtils.mkdir_p(dir_log_) unless FileTest.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	# temporary files
	temp_merged_ = "#{dir_out_}/temp.merged.coverage.txt"
	temp_coverage_ = "#{dir_out_}/temp.coverage.txt"
	header_ = "#{dir_out_}/temp.header.txt"

	qsub_ = "p16_orf_profile.#{project}.#{merge_label}.sh"
	f_qsub = open(qsub_, "w")
	f_qsub.puts "\#!/bin/bash"
	f_qsub.puts "\#$ -S /bin/bash"
	f_qsub.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G"
	f_qsub.puts "\#$ -cwd"
	f_qsub.puts "\# version #{version}"

	f_qsub.puts "echo \"orf\tlength\t#{a_samples.join("\t")}\" > #{header_}"

	#------------------------------
	# merge read covered length files
	#------------------------------
	f_qsub.puts "echo \"start covered_length\""
	# output
	label_out = "covered_length.orf.#{predictor}.#{db_label}.#{l_thre}"
	out_coverage_ = "#{dir_out_}/table.#{label_out}.txt"
	# initial output file
	f_qsub.puts "cp #{in_orf_len_} #{out_coverage_}"
	for q in 0..(a_samples.length - 1) do
		sample = a_samples[q]
		# input file of read coverages of contigs
		in_coverage_ = "#{dir_10_}/#{sample}/#{task}-#{db_label}.#{l_thre}.#{sample}.orf.covered_length.txt"
		if File.exist?(in_coverage_) then
			puts "Check OK: 'in_coverage_':#{in_coverage_} exists."
		else
			abort "Abort: 'in_coverage_':#{in_coverage_} does not exist."
		end

		f_qsub.puts "echo \"start #{sample}\""
		# covered length
		f_qsub.puts "cut -f 1,3 #{in_coverage_} > #{temp_coverage_}"
		# add coverages of current sample to final output file
		f_qsub.puts "ruby #{script_join_}  #{out_coverage_} 1  #{temp_coverage_} 1  > #{temp_merged_}"
		# concatenate missing contig as count 0
		f_qsub.puts "ruby #{script_join_}  #{out_coverage_} 1  #{temp_coverage_} 1 \"-v 1\" | grep -v orf[[:space:]]length | sed -e \"s/$/\t0/\" >> #{temp_merged_}"
		# sort by contig id
		f_qsub.puts "sort -k1,1 #{temp_merged_} > #{out_coverage_}"
	end
	# final result
	f_qsub.puts "sort -k1,1 #{temp_merged_} | cat #{header_} - > #{out_coverage_}"

	#------------------------------
	# merge read coverage files
	#------------------------------
	f_qsub.puts "echo \"start coverage\""
	# output
	label_out = "coverage.orf.#{predictor}.#{db_label}.#{l_thre}"
	out_coverage_ = "#{dir_out_}/table.#{label_out}.txt"
	# initial output file
	f_qsub.puts "cp #{in_orf_len_} #{out_coverage_}"
	for q in 0..(a_samples.length - 1) do
		sample = a_samples[q]
		# input file of read coverages of contigs
		in_coverage_ = "#{dir_10_}/#{sample}/#{task}-#{db_label}.#{l_thre}.#{sample}.orf.coverage.txt"
		if File.exist?(in_coverage_) then
			puts "Check OK: 'in_coverage_':#{in_coverage_} exists."
		else
			abort "Abort: 'in_coverage_':#{in_coverage_} does not exist."
		end

		f_qsub.puts "echo \"start #{sample}\""
		# coverage
		f_qsub.puts "awk -F '\\t' '{OFS=\"\\t\"} {print $4, $7/$5}' #{in_coverage_} > #{temp_coverage_}"
		# add coverages of current sample to final output file
		f_qsub.puts "ruby #{script_join_}  #{out_coverage_} 1  #{temp_coverage_} 1  > #{temp_merged_}"
		# concatenate missing contig as count 0
		f_qsub.puts "ruby #{script_join_}  #{out_coverage_} 1  #{temp_coverage_} 1 \"-v 1\" | grep -v orf[[:space:]]length | sed -e \"s/$/\t0/\" >> #{temp_merged_}"
		# sort by contig id
		f_qsub.puts "sort -k1,1 #{temp_merged_} > #{out_coverage_}"
	end
	# final result
	f_qsub.puts "sort -k1,1 #{temp_merged_} | cat #{header_} - > #{out_coverage_}"
	f_qsub.puts "rm #{dir_out_}/temp.*"
	f_qsub.close
	#system("qsub #{qsub_}")
	system("sh #{qsub_}")
end
