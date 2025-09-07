
#--------------------------------------------------
# Program: 13_viral_protein_search.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"
require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'n_threads:8', 'assembly_label:SPAdes_meta_c', 'db_date:20181130', 'l:', 'q:')
n_threads = h_params["n_threads"]
assembly_label = h_params["assembly_label"]
db_date = h_params["db_date"] # date of RefSeq data
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
#dir_R_ = "/usr/local/package/r/3.4.1/bin"
dir_data_ = "#{dir_home_}/data"
profile_ = "#{dir_home_}/pipeline/profile"

"---------------------"
openmpi_dir = ""
ghostmp_dir = ""
"---------------------"

l_thre = "5000.c_1500"
predictor = "prodigal"

db_label = "viral_refseq_all"
#db = "#{dir_data_}/ghostmp_db/ghostmp.RefSeq_Proteins.all_viruses.id"
#refseq_lineage_ = "#{dir_data_}/ViralRefSeq.20171101/id_lineage.RefSeq_Proteins.all_viruses.txt"
db = "./virus_ghostmp_db/Refseq_all_viruses_ghostmp_20200505"
refseq_lineage_ = "./virus_ghostmp_db/id_lineage.RefSeq_Proteins.viral_refseq_202005.txt"
# check db
if Dir::glob("#{db}.*").length == 0 then
	abort "Abort: 'db':#{db} does not exist."
else
	puts "Using 'db':#{db}"
end

n_al = "3"
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
	merge_label = inifile[section]['merge_label']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil

	query_label = "#{merge_label}_#{assembly_label}"
	# input directory
	dir_in_ = "#{dir_project_}/d12_orf_prediction/#{merge_label}"
	# output directory
	dir_out_ = "#{dir_project_}/d13_viral_protein_search/#{merge_label}"
	# log directory
	dir_log_ = "#{dir_project_}/d13_viral_protein_search/log"

	# input
	query_ = "#{dir_in_}/#{predictor}.#{query_label}.#{l_thre}.predicted.faa"
	# output
	output = "#{dir_out_}/ghostmp-#{db_label}.#{predictor}.#{query_label}.#{l_thre}.n_#{n_al}"
	count_orf_ = "#{dir_out_}/count.contig.orf_viral_family.#{db_label}.#{query_label}.#{l_thre}.txt"
	table_orf_ = "#{dir_out_}/table.contig.orf_viral_family.#{db_label}.#{query_label}.#{l_thre}.txt"

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

	# qsub file
	puts "start #{merge_label}"
	qsub_ = "p13_ghostmp-#{db_label}.#{predictor}.#{project}.#{query_label}.#{l_thre}.sh"
	f_qsub = open(qsub_, "w")
	f_qsub.puts "\#!/bin/bash"
	f_qsub.puts "\#$ -S /bin/bash"
	f_qsub.puts "\#$ -l s_vmem=21.2G,mem_req=21.2G"
	f_qsub.puts "\#$ -pe mpi_4 #{n_threads}"
	#f_qsub.puts "\#$ -pe mpi-fillup 24"
	f_qsub.puts "\#$ -cwd"
	f_qsub.puts "\# version #{version}"
	f_qsub.puts "source #{profile_}"
	f_qsub.puts "export OMP_NUM_THREADS=1"
	#f_qsub.puts "export MPI_HOME=/usr/local/package/openmpi/1.8.4-static-gcc"
	f_qsub.puts "export PATH=${openmpi_dir}/bin:${PATH}"
	f_qsub.puts "export LD_LIBRARY_PATH=${MPI_HOME}/lib:${LD_LIBRARY_PATH}"
	f_qsub.puts "BIN_=\"#{ghostmp_dir}/src/ghostmp_search\""
	f_qsub.puts "DB=\"#{db}\""
	f_qsub.puts "QUERY_=\"#{query_}\""
	f_qsub.puts "OUTPUT=\"#{output}\""
	f_qsub.puts "SCRIPT_FILTR_=\"#{dir_scripts_}/filter.blast_output.py\""
	f_qsub.puts "SCRIPT_JOIN_=#{dir_scripts_}/join_with_tab.option.rb"
	# run GHOST-MP
	f_qsub.puts "${MPI_HOME}/bin/mpiexec -np ${NSLOTS} -machinefile ${TMPDIR}/machines ${BIN_} -a ${OMP_NUM_THREADS} -q p -t p -b #{n_al} -i ${QUERY_} -d ${DB} -o ${OUTPUT}.txt"
	# filter by E-value and bitscore
	f_qsub.puts "python ${SCRIPT_FILTR_} -ev #{e_thre} -bs #{b_thre} -o ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.txt  ${OUTPUT}.txt"
	# put taxonomic lineage
	f_qsub.puts "ruby ${SCRIPT_JOIN_}  ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.txt 2  #{refseq_lineage_} 1  >  ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.with_lineage.txt"
	# put LCA lineage
	f_qsub.puts "SCRIPT_LCA_=#{dir_scripts_}/set_lca.from_lineage.majority.R"
	f_qsub.puts "R --vanilla --slave --args  ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.with_lineage.txt  ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.lca  2  kingdom,phylum,class,order,family,subfamily,genus,species  17,18,19,20,21,22,23,24  0.5  taxid  <  ${SCRIPT_LCA_}"
	# put contig name at the 1st column
	f_qsub.puts "perl -pe \"s/^(.+\\.\\d{9})__/\\1\\t\\1__/\" ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.lca.uc_.txt > ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.lca.uc_.contig.txt"
	f_qsub.puts "perl -pe \"s/^(.+\\.\\d{9})__/\\1\\t\\1__/\" ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.lca.txt > ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.lca.contig.txt"

	# count ORFs with affiliated viral families
	f_qsub.puts "cat ${OUTPUT}.e_#{e_thre}.b_#{b_thre}.lca.uc_.contig.txt | awk -F '\\t' '{OFS=\"\\t\"} {a[$1\"\\t\"$9]+=1} END{for(k in a)print k,a[k];}' | sort -k1,2 > #{count_orf_}"
	f_qsub.puts "R --vanilla --slave --args #{count_orf_} #{table_orf_} contig < #{dir_scripts_}/make.spread_table.R"

	f_qsub.close
	if q_opt != nil then
		system("qsub -q #{q_opt} #{qsub_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} #{qsub_}")
	else
		system("qsub #{qsub_}")
	end
end


