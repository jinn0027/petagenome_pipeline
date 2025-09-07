
#--------------------------------------------------
# Program: 10_contig_coverage.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'n_threads:6', 'assembly_label:SPAdes_meta_c', 'l:', 'q:')
n_threads = h_params["n_threads"]
assembly_label = h_params["assembly_label"]
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"

l_thre = "5000.c_1500"

"-----------------"
bbmap_dir = ""
samtools_dir = ""
"-----------------"

task = "bbmap"

# number of threads
if n_threads.to_i.is_a?(Integer) then
	if n_threads.to_i < 1 or n_threads.to_i > 24 then
		abort "please specify --n_threads as integer 1-24."
	end
else
	abort "please specify --n_threads as integer 1-24."
end
n_threads = n_threads.to_i.to_s
mem = (5.3 * n_threads.to_i * 0.85).round.to_s

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

	# database (contig sequences)
	db_label = "#{merge_label}_#{assembly_label}"
	dir_contig_ = "#{dir_project_}/d05_circular_contigs/#{merge_label}"
	contig_ = "#{dir_contig_}/contig.updated.#{db_label}.#{l_thre}.fa"
	abort "Abort: 'contig_':#{contig_} does not exist." if ! File.exist?(contig_)

	# output directory
	dir_out_base_ = "#{dir_project_}/d10_contig_coverage"
	FileUtils.mkdir_p(dir_out_base_) unless Dir.exist?(dir_out_base_)
	# log directory
	dir_log_ = "#{dir_out_base_}/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	#------------------------------
	# qusb for making BBMap index
	#------------------------------
	puts "start #{merge_label}"
	qsub_ref_ = "p10_#{task}.index.#{project}.#{db_label}.#{l_thre}.sh"
	f_qsub_ref = open(qsub_ref_, "w")
	f_qsub_ref.puts "\#!/bin/bash"
	f_qsub_ref.puts "\#$ -S /bin/bash"
	f_qsub_ref.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G -pe def_slot #{n_threads}"
	f_qsub_ref.puts "\#$ -cwd"
	f_qsub_ref.puts "\# version #{version}"
	f_qsub_ref.puts "BBMAP_=\"#{bbmap_dir}/bbmap/#{task}.sh\""
	# remove existing index 
	# ! To be sure no bbmap run for this index.
	f_qsub_ref.puts "if [ -e ./ref ]; then"
	f_qsub_ref.puts "  rm -R ./ref" # remove existing index
	f_qsub_ref.puts "fi"
	f_qsub_ref.puts "${BBMAP_} -Xmx#{mem}g threads=#{n_threads} ref=#{contig_}"
	f_qsub_ref.close
	if q_opt != nil then
		system("qsub -q #{q_opt} #{qsub_ref_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} #{qsub_ref_}")
	else
		system("qsub #{qsub_ref_}")
	end

	#------------------------------
	# qusb for read mapping
	#------------------------------
	for q in 0..(a_samples.length - 1) do
		sample = a_samples[q]
		# reads after error correction
		dir_read_ = "#{dir_project_}/d02_error_correction/#{sample}"
		query1_ = "#{dir_read_}/corrected/cutadapt_prinseq.#{sample}_1.00.0_0.cor.fastq.gz"
		query2_ = "#{dir_read_}/corrected/cutadapt_prinseq.#{sample}_2.00.0_0.cor.fastq.gz"
		abort "Abort: 'query1_':#{query1_} does not exist." if ! File.exist?(query1_)
		abort "Abort: 'query2_':#{query2_} does not exist." if ! File.exist?(query2_)

		dir_out_ = "#{dir_out_base_}/#{merge_label}/#{sample}"
		out = "#{dir_out_}/#{task}-#{db_label}.#{l_thre}.#{sample}"
		if Dir.exist?(dir_out_) then
			abort "Abort: 'output directory':#{dir_out_} already exists."
		else
			FileUtils.mkdir_p(dir_out_)
		end

		qsub_ = "p10_#{task}.#{sample}.#{project}.#{db_label}.#{l_thre}.sh"
		f_qsub = open(qsub_, "w")
		f_qsub.puts "\#!/bin/bash"
		f_qsub.puts "\#$ -S /bin/bash"
		f_qsub.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G -pe def_slot #{n_threads}"
		f_qsub.puts "\#$ -cwd"
		f_qsub.puts "\# version #{version}"
		f_qsub.puts "BBMAP_=\"#{bbmap_dir}/bbmap/#{task}.sh\""

		# mapping reads by bbmap
		f_qsub.puts "${BBMAP_} -Xmx#{mem}g threads=#{n_threads} in=#{query1_} in2=#{query2_}  ambiguous=random minid=0.95 pairlen=1500 out=#{out}.sam scafstats=#{out}.scafstats"
		# sort while converting to bam
		f_qsub.puts "#{samtools_dir}/samtools sort -@ #{n_threads} -O bam -o #{out}.sort.bam #{out}.sam"

		# total nucleotides covering each contig
		f_qsub.puts "echo \"contig\tlength\ttotal_nucleotides\taverage_coverage\tcovered_length\" > #{out}.coverage.txt"
		# through depth profile
		f_qsub.puts "#{samtools_dir}/samtools depth -m 10000000 -a #{out}.sort.bam | awk -F '\\t' '{OFS=\"\\t\"} {a[$1]+=$3; l[$1]+=1; if($3>0) b[$1]+=1;} END{for(k in a)print k,l[k],a[k],a[k]/l[k],b[k]/l[k];}' >> #{out}.coverage.txt"
		# covered length > 0.75
		f_qsub.puts "cat #{out}.coverage.txt | awk -F '\\t' '{OFS=\"\\t\"} {if($5>0.75) print $0}' > #{out}.coverage.0.75.txt"

		f_qsub.puts "if [ -f #{out}.sam ]; then\n  rm #{out}.sam \nfi"
		f_qsub.close
		if q_opt != nil then
			system("qsub -q #{q_opt} -hold_jid #{qsub_ref_} #{qsub_}")
		elsif l_opt != nil then
			system("qsub -l #{l_opt} -hold_jid #{qsub_ref_} #{qsub_}")
		else
			system("qsub -hold_jid #{qsub_ref_} #{qsub_}")
		end
	end
end

