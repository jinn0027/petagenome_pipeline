
#--------------------------------------------------
# Program: 02_error_correction.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'mem:128', 'l:', 'q:')
mem = h_params["mem"] # memory limit
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
profile_ = "#{dir_home_}/pipeline/profile"

"---------------------"
spades_dir = ""
fastqc_dir = ""
"---------------------"

cmd_bayeshammer = "#{spades_dir}/bin/spades.py --only-error-correction "

# memory usage
if mem == "32" then
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
	abort "please specify --mem as 32, 64, 128 or 256."
end

# load setting file
abort "Abort: './setting.cfg' does not exist." if ! File.exist?('./setting.cfg')
inifile = IniFile.load('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	dir_project_ = inifile[section]['dir_project_']
	samples = inifile[section]['samples']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')

	# output directory
	dir_out_base_ = "#{dir_project_}/d02_error_correction"
	FileUtils.mkdir_p(dir_out_base_) unless Dir.exist?(dir_out_base_)
	dir_fastqc_base_ = "#{dir_project_}/d00_fastqc"
	# log directory
	dir_log_ = "#{dir_out_base_}/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		# input
		dir_in_ = "#{dir_project_}/d01_trim_qc/#{sample}"
		query1_ = "#{dir_in_}/cutadapt_prinseq.#{sample}_1.fastq.gz"
		query2_ = "#{dir_in_}/cutadapt_prinseq.#{sample}_2.fastq.gz"
		query1_base = File.basename("#{query1_}",".fastq.gz")
		query2_base = File.basename("#{query2_}",".fastq.gz")
		abort "Abort: 'query1_':#{query1_} does not exist." if ! File.exist?(query1_)
		abort "Abort: 'query2_':#{query2_} does not exist." if ! File.exist?(query2_)

		# output
		dir_out_ = "#{dir_out_base_}/#{sample}"
		abort "Abort: 'output directory':#{dir_out_} already exists." if Dir.exist?(dir_out_)
		# fastqc directory
		dir_fastqc_ = "#{dir_fastqc_base_}/error_correction/#{sample}"
		FileUtils.mkdir_p(dir_fastqc_) unless File.exist?(dir_fastqc_)

		puts "start #{sample}"
		qsub_ = "p02_error_correction.#{project}.#{sample}.sh"
		f_qsub = open(qsub_, "w")
		f_qsub.puts "\#!/bin/bash"
		f_qsub.puts "\#$ -S /bin/bash"
		f_qsub.puts "\#$ -l s_vmem=#{s_vmem}G,mem_req=#{s_vmem}G -pe def_slot #{n_threads}"
		f_qsub.puts "\#$ -cwd"
		f_qsub.puts "\# version #{version}"
		f_qsub.puts "source #{profile_}"
		f_qsub.puts "QUERY1_=\"#{query1_}\""
		f_qsub.puts "QUERY2_=\"#{query2_}\""
		f_qsub.puts "DIR_OUT_=\"#{dir_out_}\""
		f_qsub.puts "CMD_BH=\"#{cmd_bayeshammer}\""
		f_qsub.puts "FASTQC_=\"#{fastqc_dir}/bin/fastqc\""
		f_qsub.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""

		# run Bayeshammer
		f_qsub.puts "${CMD_BH} --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}"

		# rename
		f_qsub.puts "mv ${DIR_OUT_}/corrected/#{File.basename("#{query1_}",".gz")}.00.0_0.cor.fastq.gz ${DIR_OUT_}/corrected/#{query1_base}.00.0_0.cor.fastq.gz"
		f_qsub.puts "mv ${DIR_OUT_}/corrected/#{File.basename("#{query2_}",".gz")}.00.0_0.cor.fastq.gz ${DIR_OUT_}/corrected/#{query2_base}.00.0_0.cor.fastq.gz"

		# FastQC
		f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} ${DIR_OUT_}/corrected/#{query1_base}.00.0_0.cor.fastq.gz"
		f_qsub.puts "${FASTQC_} -o #{dir_fastqc_} ${DIR_OUT_}/corrected/#{query2_base}.00.0_0.cor.fastq.gz"

		# length
		f_qsub.puts "zcat ${DIR_OUT_}/corrected/#{query1_base}.00.0_0.cor.fastq.gz | python ${SCRIPT_LENGTH_}  -t fastq  | gzip -c >  ${DIR_OUT_}/corrected/#{query1_base}.00.0_0.cor.length.txt.gz"
		f_qsub.puts "zcat ${DIR_OUT_}/corrected/#{query2_base}.00.0_0.cor.fastq.gz | python ${SCRIPT_LENGTH_}  -t fastq  | gzip -c >  ${DIR_OUT_}/corrected/#{query2_base}.00.0_0.cor.length.txt.gz"

		f_qsub.close
		if q_opt != nil then
			system("qsub -q #{q_opt} #{qsub_}")
		elsif l_opt != nil then
			system("qsub -l #{l_opt} #{qsub_}")
		else
			system("qsub #{qsub_}")
		end
	end
end
