
#--------------------------------------------------
# Program: 14_kegg_search.split.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'no_split', 'db_date:20180916', 'assembly_label:SPAdes_meta_c', 'l:', 'q:')
no_split = h_params["no_split"] # without spliting input fasta when it is already split
db_date = h_params["db_date"] # date of KEGG data
assembly_label = h_params["assembly_label"]
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"
profile_ = "#{dir_home_}/pipeline/profile"


"---------------------"
openmpi_dir = ""
ghostmp_dir = ""
"---------------------"


l_thre = "5000.c_1500"
predictor = "prodigal"

db_label = "kegg_prokaryotes"
#db = "#{dir_data_}/ghostmp_db/ghostmp.kegg.prokaryotes"
db = "./kegg.20180916.genes/ghostmp.kegg.prokaryotes"
# check db
if Dir::glob("#{db}.*").length == 0 then
	abort "Abort: 'db':#{db} does not exist."
else
	puts "Using 'db':#{db}"
end

n_al = "1"
e_thre = "1e-5"
b_thre = 50

# loading setting file
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
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')

	query_label = "#{merge_label}_#{assembly_label}"
	# input directory
	dir_in_ = "#{dir_project_}/d12_orf_prediction/#{merge_label}"
	# output directory
	dir_out_ = "#{dir_project_}/d14_kegg_search/#{merge_label}"
	# log directory
	dir_log_ = "#{dir_project_}/d14_kegg_search/log"

	# input
	query_ = "#{dir_in_}/#{predictor}.#{query_label}.#{l_thre}.predicted.faa"
	# output
	output = "#{dir_out_}/ghostmp-#{db_label}.#{predictor}.#{query_label}.#{l_thre}.n_#{n_al}"

	# check
	abort "Abort: 'query_':#{query_} does not exist." if ! File.exist?(query_)
	if Dir.exist?(dir_out_) then
		abort "Abort: 'output directory':#{dir_out_} already exists."
	else
		FileUtils.mkdir_p(dir_out_)
	end
	# check the existence of the GHOST-MP output file
	# since it causes an error in GHOST-MP
	abort "Abort: ghostmp output file already exists: #{output}.txt" if File.exist?("#{output}.txt")

	# move to log directory
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	Dir.chdir(dir_log_)

	puts "start #{merge_label}"

	# split input file
	n_split = a_samples.length
	if no_split then
		for i in 1..n_split do
			query_split_ = "#{query_}.group_#{i}.fa"
			if File.exist?(query_split_) then
				puts "Check OK: 'query_split_':#{query_split_} exists."
			else
				abort "Abort: 'query_split_':#{query_split_} does not exist."
			end
		end
	else
		system("source #{profile_}; python #{dir_scripts_}/random_split_large_file.py -n #{n_split.to_s} -o #{query_} #{query_}")
	end


	a_jid = [] # job names to be finished before concatenation
	for i in 1..n_split do
		split_label = "#{query_label}.group_#{i}"

		query_split_ = "#{query_}.group_#{i}.fa"
		output_split = "#{dir_out_}/ghostmp-#{db_label}.#{predictor}.#{split_label}.#{l_thre}.n_#{n_al}"

		# qsub file
		qsub_ = "p14_ghostmp-#{db_label}.#{predictor}.#{project}.#{split_label}.#{l_thre}.sh"
		f_qsub = open(qsub_, "w")
		f_qsub.puts "\#!/bin/bash"
		f_qsub.puts "\#$ -S /bin/bash"
		f_qsub.puts "\#$ -l s_vmem=42.4G,mem_req=42.4G"
		f_qsub.puts "\#$ -pe mpi-fillup 8"
		f_qsub.puts "\#$ -cwd"
		f_qsub.puts "\# version #{version}"
		f_qsub.puts "export OMP_NUM_THREADS=1"
		#f_qsub.puts "export MPI_HOME=/usr/local/package/openmpi/3.1.0-icc"
		f_qsub.puts "export PATH=${openmpi_dir}/bin:${PATH}"
		f_qsub.puts "export LD_LIBRARY_PATH=${openmpi_dir}/lib:${LD_LIBRARY_PATH}"
		f_qsub.puts "BIN_=\"#{ghostmp_dir}/src/ghostmp_search\""
		f_qsub.puts "DB=\"#{db}\""
		f_qsub.puts "QUERY_=\"#{query_split_}\""
		f_qsub.puts "OUTPUT=\"#{output_split}\""
		# run GHOST-MP
		f_qsub.puts "${openmpi_dir}/bin/mpiexec -np ${NSLOTS} -machinefile ${TMPDIR}/machines ${BIN_} -a ${OMP_NUM_THREADS} -q p -t p -b #{n_al} -i ${QUERY_} -d ${DB} -o ${OUTPUT}.txt"
		f_qsub.close
		if q_opt != nil then
			system("qsub -q #{q_opt} -N #{qsub_} #{qsub_}")
		elsif l_opt != nil then
			system("qsub -l #{l_opt} -N #{qsub_} #{qsub_}")
		else
			system("qsub -N #{qsub_} #{qsub_}")
		end
		a_jid.push("#{qsub_}")
	end

	# qsub for concatenation
	split_label_all = "#{query_label}.group_*"
	output_split_all = "#{dir_out_}/ghostmp-#{db_label}.#{predictor}.#{split_label_all}.#{l_thre}.n_#{n_al}"

	qsub_cat_ = "p14_cat_ghostmp-#{db_label}.#{predictor}.#{project}.#{query_label}.#{l_thre}.sh"
	f_qsub_cat = open(qsub_cat_, "w")
	f_qsub_cat.puts "\#!/bin/bash"
	f_qsub_cat.puts "\#$ -S /bin/bash"
	f_qsub_cat.puts "\#$ -l s_vmem=21.2G,mem_req=21.2G"
	f_qsub_cat.puts "\#$ -cwd"
	f_qsub_cat.puts "OUTPUT=\"#{output}\""
	f_qsub_cat.puts "SCRIPT_FILTR_=\"#{dir_scripts_}/filter.blast_output.py\""
	f_qsub_cat.puts "SCRIPT_INFO_=\"#{dir_scripts_}/assign_kegg.blast.py\""
	f_qsub_cat.puts "\# version #{version}"
	# concatenate results
	f_qsub_cat.puts "ls -lh #{output_split_all}.txt"
	f_qsub_cat.puts "n_done=$(ls -lh #{output_split_all}.txt | wc -l | cut -d' ' -f 1 )"
	f_qsub_cat.puts "if [ ${n_done} -eq #{n_split} ]; then"
	f_qsub_cat.puts "  cat #{output_split_all}.txt > ${OUTPUT}.txt"
	f_qsub_cat.puts "  wait"
	#f_qsub_cat.puts "  rm #{output_split_all}.txt"
	# filter by E-value and bitscore
	f_qsub_cat.puts "  python ${SCRIPT_FILTR_} -ev #{e_thre} -bs #{b_thre} -o ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.txt  ${OUTPUT}.txt"
	# assign KEGG information
	# using KEGG information file in kegg.#{db_date}
	f_qsub_cat.puts "  python ${SCRIPT_INFO_} -o ${OUTPUT}.e_#{e_thre}.b_#{b_thre} -d #{db_date} ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.txt"
	# put contig name at the 1st column
	f_qsub_cat.puts "  perl -pe \"s/^(.+\\.\\d{9})__/\\1\\t\\1__/\" ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.info.txt > ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.info.contig.txt"
	f_qsub_cat.puts "  cut -f 2 ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.ko.txt | perl -pe \"s/^(.+\\.\\d{9})__.+$/\\1/\" | paste ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.ko.txt - > ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.ko.contig.txt"
	f_qsub_cat.puts "  cut -f 3 ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.pathway.txt | perl -pe \"s/^(.+\\.\\d{9})__.+$/\\1/\" | paste ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.pathway.txt - > ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.pathway.contig.txt"
	f_qsub_cat.puts "  cut -f 3 ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.module.txt | perl -pe \"s/^(.+\\.\\d{9})__.+$/\\1/\" | paste ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.module.txt - > ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.module.contig.txt"
	f_qsub_cat.puts "else"
	f_qsub_cat.puts "  echo \"All the split files (#{n_split}) are not finished. ${n_done} files are done.\""
	f_qsub_cat.puts "fi"

	f_qsub_cat.close
	# qsub a job that will run after all the jobs above
	if q_opt != nil then
		system("qsub -q #{q_opt} -hold_jid #{a_jid.join(",")} #{qsub_cat_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} -hold_jid #{a_jid.join(",")} #{qsub_cat_}")
	else
		system("qsub -hold_jid #{a_jid.join(",")} #{qsub_cat_}")
	end
end
