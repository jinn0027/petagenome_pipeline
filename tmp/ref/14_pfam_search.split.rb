
#--------------------------------------------------
# Program: 14_pfam_search.split.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'n_threads:6', 'assembly_label:SPAdes_meta_c', 'pfam_version:32.0', 'no_split', 'n_files:', 'l:', 'q:')
n_threads = h_params["n_threads"] # the number of threads used in BLASTN search
assembly_label = h_params["assembly_label"]
pfam_ver = h_params["pfam_version"] # pfam version
no_split = h_params["no_split"] # without spliting input fasta when it is already split
n_files = h_params["n_files"] # the number of files after split
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"
profile_ = "#{dir_home_}/pipeline/profile"

"---------------------"
prodigal_dir = ""
hmmscan_dir = ""
"---------------------"

l_thre = "5000.c_1500"
predictor = "prodigal"

db_label = "pfam_A"
#dir_pfam_ = "#{dir_data_}/Pfam.#{pfam_ver}"
db = "./Pfam.32.0/Pfam-A.hmm"
# check db
if Dir::glob("#{db}.*").length == 0 then
	abort "Abort: 'db':#{db} does not exist."
else
	puts "Using 'db':#{db}"
end

e_thre = "1e-5"

# loading setting file
inifile = IniFile.load('./setting.cfg')
abort "Abort: './setting.cfg' does not exist." if ! File.exist?('./setting.cfg')

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

	query_label = "#{merge_label}_#{assembly_label}"
	# input directory
	dir_in_ = "#{dir_project_}/d12_orf_prediction/#{merge_label}"
	# output directory
	dir_out_ = "#{dir_project_}/d14_pfam_search/#{merge_label}"
	FileUtils.mkdir_p(dir_out_) unless FileTest.exist?(dir_out_)
	# log directory
	dir_log_ = "#{dir_project_}/d14_pfam_search/log"
	FileUtils.mkdir_p(dir_log_) unless FileTest.exist?(dir_log_)
	# move to log directory
	Dir::chdir(dir_log_)

	# input
	query_ = "#{dir_in_}/#{predictor}.#{query_label}.#{l_thre}.predicted.faa"
	abort "Abort: 'query_':#{query_} does not exist." if ! File.exist?(query_)

	# output
	output = "#{dir_out_}/hmmer-#{db_label}.#{predictor}.#{query_label}.#{l_thre}"

	puts "start #{merge_label}"

	# number of split
	if n_files == nil then
		abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
		a_samples = samples.gsub(" ","").split(',')
		n_files = a_samples.length
	end
	# split input file
	if no_split then
		for i in 1..n_files.to_i do
			query_split_ = "#{query_}.group_#{i}.fa"
			if File.exist?(query_split_) then
				puts "Check OK: 'query_split_':#{query_split_} exists."
			else
				abort "Abort: 'query_split_':#{query_split_} does not exist."
			end
		end
	else
		system("source #{profile_}; python #{dir_scripts_}/random_split_large_file.py -n #{n_files.to_s} -o #{query_} #{query_}")
	end

	a_jid = [] # job names to be finished before concatenation
	for i in 1..n_files.to_i do
		split_label = "#{query_label}.group_#{i}"

		query_split_ = "#{query_}.group_#{i}.fa"
		output_split = "#{dir_out_}/hmmer-#{db_label}.#{predictor}.#{split_label}.#{l_thre}"

		# qsub file
		qsub_ = "p14_hmmer-#{db_label}.#{predictor}.#{project}.#{split_label}.#{l_thre}.sh"
		f_qsub = open(qsub_, "w")
		f_qsub.puts "\#!/bin/bash"
		f_qsub.puts "\#$ -S /bin/bash"
		f_qsub.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G -pe def_slot #{n_threads}"
		f_qsub.puts "\#$ -cwd"
		f_qsub.puts "HMMSCAN_=\"#{hmmscan_dir}/bin/hmmscan\""
		f_qsub.puts "DB=\"#{db}\""
		f_qsub.puts "QUERY_=\"#{query_split_}\""
		f_qsub.puts "OUT=\"#{output_split}\""
		f_qsub.puts "\# version #{version}"
		# run hmmscan
		f_qsub.puts "${HMMSCAN_} --cpu #{n_threads}  -E 1e-5  --tblout ${OUT}.txt  -o ${OUT}.stdout  ${DB}  ${QUERY_}"
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
	output_split_all = "#{dir_out_}/hmmer-#{db_label}.#{predictor}.#{split_label_all}.#{l_thre}"

	qsub_cat_ = "p14_cat_hmmer-#{db_label}.#{predictor}.#{project}.#{query_label}.#{l_thre}.sh"
	f_qsub_cat = open(qsub_cat_, "w")
	f_qsub_cat.puts "\#!/bin/bash"
	f_qsub_cat.puts "\#$ -S /bin/bash"
	f_qsub_cat.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G"
	f_qsub_cat.puts "\#$ -cwd"
	f_qsub_cat.puts "OUT=\"#{output}\""
	f_qsub_cat.puts "SCRIPT_INFO_=#{dir_scripts_}/assign_pfam.hmmer.py"
	f_qsub_cat.puts "SCRIPT_BEST_=#{dir_scripts_}/get_best_match.from_any_outputs.py"
	f_qsub_cat.puts "SCRIPT_CATEGORY_=#{dir_scripts_}/get_best_pfam_category.py"
	f_qsub_cat.puts "\# version #{version}"
	# concatenate results
	f_qsub_cat.puts "ls -lh #{output_split_all}.txt"
	f_qsub_cat.puts "n_done=$(ls -lh #{output_split_all}.txt | wc -l | cut -d' ' -f 1 )"
	f_qsub_cat.puts "if [ ${n_done} -eq #{n_files} ]; then"
	f_qsub_cat.puts "  cat #{output_split_all}.txt > ${OUT}.txt"
	f_qsub_cat.puts "  cat #{output_split_all}.stdout | gzip -c > ${OUT}.stdout.gz"
	f_qsub_cat.puts "  wait"
	f_qsub_cat.puts "  rm #{output_split_all}.txt"
	f_qsub_cat.puts "  rm #{output_split_all}.stdout"
	# assign Pfam information
	# using pfam information file in #{dir_pfam_}
	f_qsub_cat.puts "  python  ${SCRIPT_INFO_}  -p #{dir_pfam_}  -e #{e_thre}  -o ${OUT}.e_#{e_thre}.info.txt  ${OUT}.txt"
	f_qsub_cat.puts "  perl -pe \"s/^(.+\\.\\d{9})__/\\1\\t\\1__/\" ${OUT}.e_#{e_thre}.info.txt > ${OUT}.e_#{e_thre}.info.contig.txt"
	# extract most significant hit for each ORF
	f_qsub_cat.puts "  python  ${SCRIPT_BEST_}  -i 0  -s 5  ${OUT}.e_#{e_thre}.info.txt > ${OUT}.e_#{e_thre}.info.best.txt"
	# extract most significant hit having functional category
	f_qsub_cat.puts "  python  ${SCRIPT_CATEGORY_} ${OUT}.e_#{e_thre}.info.txt > ${OUT}.e_#{e_thre}.info.best_category.txt"
	f_qsub_cat.puts "else"
	f_qsub_cat.puts "  echo \"All the split files (#{n_files}) are not finished. ${n_done} files are done.\""
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


