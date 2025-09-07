
#--------------------------------------------------
# Program: 01_trim_qc.HiSeq.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'l:', 'q:')
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
profile_ = "#{dir_home_}/pipeline/profile"

"---------------------"
cutadapt_dir = ""
prinseq_dir = ""
fastqc_dir = ""
"---------------------"

cmd_cutadapt = "#{cutadapt_dir}/bin/cutadapt --minimum-length 50"
cmd_prinseq = "#{prinseq_dir}/bin/prinseq-lite.pl -trim_right 10 -trim_left 10 -trim_qual_right 20 -trim_qual_left 20 -trim_qual_window 20 -min_len 75 -derep 1 -lc_method dust -lc_threshold 7 -trim_ns_right 1 -ns_max_n 0"

ADPT_FWD="AATGATACGGCGACCACCGAGAUCTACAC"
ADPT_REV="CAAGCAGAAGACGGCATACGAGAT"

# loading raw data information
abort "Abort: './info.raw_data.cfg' does not exist." if ! File.exist?('./info.raw_data.cfg')
inifile = IniFile.load('./info.raw_data.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	dir_project_ = inifile[section]['dir_project_']
	dir_run_ = inifile[section]['dir_run_']
	samples = inifile[section]['samples']
	# check settings
	abort "Abort: no 'project' in ./info.raw_data.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./info.raw_data.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'dir_run_' in ./info.raw_data.cfg" if dir_run_ == nil
	abort "Abort: 'dir_run_':#{dir_run_} does not exist." if ! Dir.exist?(dir_run_)
	abort "Abort: no 'samples' in ./info.raw_data.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')

	# output directory
	dir_out_base_ = "#{dir_project_}/d01_trim_qc"
	FileUtils.mkdir_p(dir_out_base_) unless FileTest.exist?(dir_out_base_)
	dir_fastqc_base_ = "#{dir_project_}/d00_fastqc"
	# log directory
	dir_log_ = "#{dir_out_base_}/log"
	FileUtils.mkdir_p(dir_log_) unless FileTest.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		puts "start #{sample}"
		abort "Abort: '#{dir_run_}/Sample_#{sample}' does not exist." if ! Dir.exist?("#{dir_run_}/Sample_#{sample}")

		# fastqc directory
		dir_fastqc_ = "#{dir_fastqc_base_}/trim_qc/#{sample}"
		FileUtils.mkdir_p(dir_fastqc_) unless Dir.exist?(dir_fastqc_)
		# output directory of each sample
		dir_out_ = "#{dir_out_base_}/#{sample}"

		if Dir.exist?(dir_out_) then
			abort "Abort: 'output directory':#{dir_out_} already exists."
		else
			FileUtils.mkdir_p(dir_out_)
		end

		# arrays to store information for concatenation
		a_hold_jid = []
		a_fwd_qc_fq_gz_ = []
		a_rev_qc_fq_gz_ = []

		# run cutadapt and prinseq for each fastq file
		Dir::glob("#{dir_run_}/Sample_#{sample}/*_R1_???.fastq.gz").each do |query1_gz_|
			query1_ = query1_gz_.sub(".fastq.gz",".fastq")
			query2_gz_ = query1_gz_.dup
			query2_gz_[/(_R1_)\d{3}.fastq.gz/, 1] = "_R2_"
			query2_ = query2_gz_.sub(".fastq.gz",".fastq")
			abort "Abort: '#{query2_gz_}' does not exist." if ! File.exist?(query2_gz_)
			abort "Abort: 'query1_gz_':#{query1_gz_} and 'query2_gz_':#{query2_gz_} are the same." if query1_gz_ == query2_gz_

			# fastq.gz files are uncompressed in the output directory
			original_dir_ = File.dirname("#{query1_}")
			query1_ = query1_.sub(original_dir_, dir_out_)
			query2_ = query2_.sub(original_dir_, dir_out_)

			label = File.basename("#{query1_gz_}",".fastq.gz").sub(/^(#{sample})_[a-zA-Z0-9\-]+_(L\d+)_R1_(\d+)$/,'\1_\2_\3')

			# output file base names
			out_cutadapt = "#{dir_out_}/cutadapt.#{label}"
			out_prinseq = "#{dir_out_}/cutadapt_prinseq.#{label}"

			puts "start #{label}"
			qsub_ = "p01_qc_trim.#{project}.#{label}.sh"
			f_qsub = open(qsub_, "w")
			f_qsub.puts "\#!/bin/bash"
			f_qsub.puts "\#$ -S /bin/bash"
			f_qsub.puts "\#$ -l s_vmem=16G,mem_req=16G"
			f_qsub.puts "\#$ -cwd"
			f_qsub.puts "\# version #{version}"
			f_qsub.puts "source #{profile_}"
			f_qsub.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""
			f_qsub.puts "FASTQC_=\"#{fastqc_dir}/bin/fastqc\""

			f_qsub.puts "if [ ! -f #{query1_} ]; then\n  gunzip -c #{query1_gz_} > #{query1_} \nfi"
			f_qsub.puts "if [ ! -f #{query2_} ]; then\n  gunzip -c #{query2_gz_} > #{query2_} \nfi"
			f_qsub.puts "QUERY1_=\"#{query1_}\""
			f_qsub.puts "QUERY2_=\"#{query2_}\""
			f_qsub.puts "ADPT_FWD=\"#{adapter_fwd}\""
			f_qsub.puts "ADPT_REV=\"#{adapter_rev}\""
			f_qsub.puts "CMD_CUTADAPT=\"#{cmd_cutadapt}\""
			f_qsub.puts "CMD_PRINSEQ=\"#{cmd_prinseq}\""

			# cutadapt
			f_qsub.puts "echo \"start #{cmd_cutadapt}\""
			f_qsub.puts "${CMD_CUTADAPT} -a ${ADPT_FWD} -o #{out_cutadapt}_1.fastq  ${QUERY1_}"
			f_qsub.puts "${CMD_CUTADAPT} -a ${ADPT_REV} -o #{out_cutadapt}_2.fastq  ${QUERY2_}"

			# FASTQC
			f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} #{out_cutadapt}_1.fastq"
			f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} #{out_cutadapt}_2.fastq"

			# prinseq
			f_qsub.puts "echo \"start #{cmd_prinseq}\""
			f_qsub.puts "${CMD_PRINSEQ} -fastq #{out_cutadapt}_1.fastq -fastq2 #{out_cutadapt}_2.fastq -out_good #{out_prinseq} -out_bad #{out_prinseq}_bad"

			# length
			f_qsub.puts "python ${SCRIPT_LENGTH_}  -t fastq  #{out_cutadapt}_1.fastq > #{out_cutadapt}_1.length.txt"
			f_qsub.puts "python ${SCRIPT_LENGTH_}  -t fastq  #{out_cutadapt}_2.fastq > #{out_cutadapt}_2.length.txt"
			f_qsub.puts "python ${SCRIPT_LENGTH_}  -t fastq  #{out_prinseq}_1.fastq > #{out_prinseq}_1.length.txt"
			f_qsub.puts "python ${SCRIPT_LENGTH_}  -t fastq  #{out_prinseq}_2.fastq > #{out_prinseq}_2.length.txt"

			# FASTQC
			f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} #{out_prinseq}_1.fastq"
			f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} #{out_prinseq}_2.fastq"
			f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} #{out_prinseq}_1_singletons.fastq"
			f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} #{out_prinseq}_2_singletons.fastq"

			# remove or compress files
			f_qsub.puts "if [ -f #{query1_} ]; then\n  rm #{query1_} \nfi"
			f_qsub.puts "if [ -f #{query2_} ]; then\n  rm #{query2_} \nfi"
			f_qsub.puts "gzip -f #{out_cutadapt}_1.length.txt"
			f_qsub.puts "gzip -f #{out_cutadapt}_2.length.txt"
			f_qsub.puts "rm #{out_cutadapt}_1.fastq"
			f_qsub.puts "rm #{out_cutadapt}_2.fastq"
			f_qsub.puts "gzip -f #{out_prinseq}_1.fastq"
			f_qsub.puts "gzip -f #{out_prinseq}_2.fastq"
			f_qsub.puts "gzip -f #{out_prinseq}_1_singletons.fastq"
			f_qsub.puts "gzip -f #{out_prinseq}_2_singletons.fastq"
			f_qsub.puts "gzip -f #{out_prinseq}_bad_1.fastq"
			f_qsub.puts "gzip -f #{out_prinseq}_bad_2.fastq"
			f_qsub.puts "gzip -f #{out_prinseq}_1.length.txt"
			f_qsub.puts "gzip -f #{out_prinseq}_2.length.txt"
			f_qsub.close
			if q_opt != nil then
				system("qsub -q #{q_opt} -N #{qsub_} #{qsub_}")
			elsif l_opt != nil then
				system("qsub -l #{l_opt} -N #{qsub_} #{qsub_}")
			else
				system("qsub -N #{qsub_} #{qsub_}")
			end

			# store quality filtered fastq file names for concatenation
			a_hold_jid.push(qsub_)
			a_fwd_qc_fq_gz_.push("#{out_prinseq}_1.fastq.gz")
			a_rev_qc_fq_gz_.push("#{out_prinseq}_2.fastq.gz")
		end

		# concatenate the quality filtered fastq files
		qsub_cat_ = "p01_qc_trim.cat.#{project}.#{sample}.sh"
		output_fwd = "#{dir_out_}/cutadapt_prinseq.#{sample}_1"
		output_rev = "#{dir_out_}/cutadapt_prinseq.#{sample}_2"
		f_qsub_cat = open(qsub_cat_, "w")
		f_qsub_cat.puts "\#!/bin/bash"
		f_qsub_cat.puts "\#$ -S /bin/bash"
		f_qsub_cat.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G"
		f_qsub_cat.puts "\#$ -cwd"
		f_qsub_cat.puts "\# version #{version}"
		f_qsub_cat.puts "source #{profile_}"
		f_qsub_cat.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""
		# concatenate
		f_qsub_cat.puts "zcat #{a_fwd_qc_fq_gz_.join(" ")} | gzip -c > #{output_fwd}.fastq.gz"
		f_qsub_cat.puts "zcat #{a_rev_qc_fq_gz_.join(" ")} | gzip -c > #{output_rev}.fastq.gz"
		# length
		f_qsub_cat.puts "zcat #{output_fwd}.fastq.gz | python ${SCRIPT_LENGTH_} -t fastq | gzip -c > #{output_fwd}.length.txt.gz"
		f_qsub_cat.puts "zcat #{output_rev}.fastq.gz | python ${SCRIPT_LENGTH_} -t fastq | gzip -c > #{output_rev}.length.txt.gz"
		f_qsub_cat.close
		# qsub a job that will run after all the jobs above
		if q_opt != nil then
			system("qsub -q #{q_opt} -hold_jid #{a_hold_jid.join(",")} #{qsub_cat_}")
		elsif l_opt != nil then
			system("qsub -l #{l_opt} -hold_jid #{a_hold_jid.join(",")} #{qsub_cat_}")
		else
			system("qsub -hold_jid #{a_hold_jid.join(",")} #{qsub_cat_}")
		end
	end
end

