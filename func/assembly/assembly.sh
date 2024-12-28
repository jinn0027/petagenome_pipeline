#!/bin/bash

exit


#--------------------------------------------------
# Program: 03_assembly.rb
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
#dir_R_ = "/usr/local/package/r/3.4.1/bin"
profile_ = "#{dir_home_}/pipeline/profile

"---------------------"
spades_dir = ""
blast_dir = ""
"---------------------"

cmd_spades = "#{spades_dir}/bin/spades.py --only-assembler "
assembly_label = "SPAdes_meta"

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
abort "Abort: ./setting.cfg does not exist." if ! File.exist?('./setting.cfg')
inifile = IniFile.load('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	type = inifile[section]['type']
	dir_project_ = inifile[section]['dir_project_']
	samples = inifile[section]['samples']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'type' in ./setting.cfg" if type == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')

	# check project type
	if type == "virome" then
		prefix = "v"
		abort "Abort: 'type' is #{type} but '#{$0}' was run." if File.dirname($0) !~ /virome$/
	elsif type == "bacteriome" then
		prefix = "b"
		abort "Abort: 'type' is #{type} but '#{$0}' was run." if File.dirname($0) !~ /bacteriome$/
	else
		abort "Abort: 'type' in ./setting.cfg must be 'virome' or 'bacteriome'."
	end

	# output directory
	dir_out_base_ = "#{dir_project_}/d03_assembly"
	FileUtils.mkdir_p(dir_out_base_) unless Dir.exist?(dir_out_base_)
	# log directory
	dir_log_ = "#{dir_out_base_}/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		# input
		dir_in_ = "#{dir_project_}/d02_error_correction/#{sample}"
		query1_ = "#{dir_in_}/corrected/cutadapt_prinseq.#{sample}_1.00.0_0.cor.fastq.gz"
		query2_ = "#{dir_in_}/corrected/cutadapt_prinseq.#{sample}_2.00.0_0.cor.fastq.gz"
		abort "Abort: 'query1_':#{query1_} does not exist." if ! File.exist?(query1_)
		abort "Abort: 'query2_':#{query2_} does not exist." if ! File.exist?(query2_)

		# output
		dir_out_ = "#{dir_out_base_}/#{sample}"
		abort "Abort: 'output directory':#{dir_out_} already exists." if Dir.exist?(dir_out_)
		out = "#{dir_out_}/contig.#{sample}_#{assembly_label}"

		puts "start #{sample}"
		qsub_ = "p03_assembly.#{project}.#{sample}_#{assembly_label}.sh"
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
		f_qsub.puts "CMD_SPADES_=\"#{cmd_spades}\""
		f_qsub.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""
		f_qsub.puts "SCRIPT_RENAME_=\"#{dir_scripts_}/filter_contig.rename.py\""
		f_qsub.puts "SCRIPT_STATS_=\"#{dir_scripts_}/stats.assembly.R\""
		f_qsub.puts "MAKEBLASTDB_=\"#{blast_dir}/makeblastdb\""
		# SPAdes_meta
		f_qsub.puts "${CMD_SPADES_} --pe1-1 ${QUERY1_} --pe1-2 ${QUERY2_} --meta --memory #{mem} --threads #{n_threads} -o ${DIR_OUT_}"

		# length filter
		f_qsub.puts "python ${SCRIPT_RENAME_} -min 1000 --rename --prefix #{prefix}.#{sample}.n. --table #{out}.1000.name.txt ${DIR_OUT_}/scaffolds.fasta > #{out}.1000.fa"
		f_qsub.puts "python ${SCRIPT_RENAME_} -min 5000 #{out}.1000.fa > #{out}.5000.fa"
		f_qsub.puts "python ${SCRIPT_RENAME_} -min 10000 #{out}.1000.fa > #{out}.10000.fa"

		# get length of contigs
		f_qsub.puts "python ${SCRIPT_LENGTH_} #{out}.1000.fa > #{out}.1000.length.txt"
		f_qsub.puts "python ${SCRIPT_LENGTH_} #{out}.5000.fa > #{out}.5000.length.txt"
		f_qsub.puts "python ${SCRIPT_LENGTH_} #{out}.10000.fa > #{out}.10000.length.txt"

		# stats of assemblies
		f_qsub.puts "R --vanilla --slave --args #{out}.1000.length.txt 2 < ${SCRIPT_STATS_} > #{out}.1000.stats.txt"
		f_qsub.puts "R --vanilla --slave --args #{out}.5000.length.txt 2 < ${SCRIPT_STATS_} > #{out}.5000.stats.txt"
		f_qsub.puts "R --vanilla --slave --args #{out}.10000.length.txt 2 < ${SCRIPT_STATS_} > #{out}.10000.stats.txt"

		# blastdb
		f_qsub.puts "${MAKEBLASTDB_} -in #{out}.1000.fa -out #{out}.1000 -dbtype nucl -parse_seqids"
		f_qsub.puts "${MAKEBLASTDB_} -in #{out}.5000.fa -out #{out}.5000 -dbtype nucl -parse_seqids"
		f_qsub.puts "${MAKEBLASTDB_} -in #{out}.10000.fa -out #{out}.10000 -dbtype nucl -parse_seqids"

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

