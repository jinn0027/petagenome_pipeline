
#--------------------------------------------------
# Program: 06_prophage_detection.each_sample.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'n_threads:6', 'assembly_label:SPAdes_meta', 'l:', 'q:')
n_threads = h_params["n_threads"]
assembly_label = h_params["assembly_label"]
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"
profile_ = "#{dir_home_}/pipeline/profile"

"---------------------"
virsorter_dir = ""
"---------------------"

l_thre = "5000"

# virsorter database
dir_data_vs_ = "#{dir_data_}/virsorter-data"
a_vs_db_ids = ["1", "2"]
h_vs_db = {"1" => "refseqdb", "2" => "viromedb"}

# use microbiome mode
opt_virome = ""
opt_label = ""

# number of threads
if n_threads.to_i.is_a?(Integer) then
	if n_threads.to_i < 1 or n_threads.to_i > 24 then
		abort "please specify --n_threads as integer 1-24."
	end
else
	abort "please specify --n_threads as integer 1-24."
end
n_threads = n_threads.to_i.to_s

# loading setting file
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
	dir_out_base_ = "#{dir_project_}/d06_prophage_detection"
	FileUtils.mkdir_p(dir_out_base_) unless Dir.exist?(dir_out_base_)
	# log directory
	dir_log_ = "#{dir_out_base_}/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)


	# move to log directory
	Dir.chdir(dir_log_)

	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		dir_in_ = "#{dir_project_}/d03_assembly/#{sample}"
		query_label = "#{sample}_#{assembly_label}"
		# input
		query_fa_ = "#{dir_in_}/contig.#{query_label}.#{l_thre}.fa"
		abort "Abort: 'query_fa_':#{query_fa_} does not exist." if ! File.exist?(query_fa_)

		for vs_db_id in a_vs_db_ids
			vs_db_label = h_vs_db[vs_db_id]
			# output
			dir_out_ = "#{dir_out_base_}/#{sample}/virsorter#{opt_label}-#{vs_db_label}.#{query_label}.#{l_thre}"
			if Dir.exist?(dir_out_) then
				abort "Abort: 'output directory':#{dir_out_} already exists."
			else
				FileUtils.mkdir_p(dir_out_)
			end
	
			query_org_id_ = "#{dir_out_}/input_sequences.original_id.txt"
			query_cnv_id_ = "#{dir_out_}/input_sequences.converted_id.txt"
			query_id_table_ = "#{dir_out_}/input_sequences.id.table.txt"
	
			# qsub file
			puts "start #{sample}"
			qsub_ = "p06_01_virsorter#{opt_label}-#{vs_db_label}.#{project}.#{query_label}.#{l_thre}.sh"
			f_qsub = open(qsub_, "w")
			f_qsub.puts "\#!/bin/bash"
			f_qsub.puts "\#$ -S /bin/bash"
			f_qsub.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G -pe def_slot #{n_threads}"
			f_qsub.puts "\#$ -cwd"
			f_qsub.puts "\# version #{version}"
			f_qsub.puts "source #{profile_}"
			f_qsub.puts "VIRSORTER_=\"#{virsorter_dir}/wrapper_phage_contigs_sorter_iPlant.pl\""
			f_qsub.puts "SCRIPT_JOIN_=#{dir_scripts_}/join_with_tab.rb"
	
			# run VirSorter
			f_qsub.puts "perl ${VIRSORTER_} #{opt_virome} --db #{vs_db_id} -data-dir #{dir_data_vs_} -ncpu #{n_threads} -wdir #{dir_out_} -f #{query_fa_}"
			# parse main result
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.csv | perl -pe \"s/,/\\t/g\" | perl -pe \"s/VIRSorter_//\" | perl -pe \"s/-circular//\" > #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv"
	
			# row number of each category header
			f_qsub.puts "N_cat_2=$(grep -e \"## 2\" -n #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv | sed -e 's/:.*//')"
			f_qsub.puts "N_cat_4=$(grep -e \"## 4\" -n #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv | sed -e 's/:.*//')"
			f_qsub.puts "N_cat_3=$(grep -e \"## 3\" -n #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv | sed -e 's/:.*//')"
			f_qsub.puts "N_cat_5=$(grep -e \"## 5\" -n #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv | sed -e 's/:.*//')"
			f_qsub.puts "N_cat_6=$(grep -e \"## 6\" -n #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv | sed -e 's/:.*//')"
			# extract rows belonging to each category
			f_qsub.puts "sed -e \"${N_cat_4},\\$d\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_1_2_3.tsv"
			f_qsub.puts "sed -e \"${N_cat_3},\\$d\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_1_2.tsv"
			f_qsub.puts "sed -e \"${N_cat_2},\\$d\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_1.tsv"
			f_qsub.puts "sed -n -e \"${N_cat_2},${N_cat_3}p\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_2.tsv"
			f_qsub.puts "sed -n -e \"${N_cat_3},${N_cat_4}p\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_3.tsv"
			f_qsub.puts "sed -n -e \"${N_cat_4},${N_cat_5}p\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_4.tsv"
			f_qsub.puts "sed -n -e \"${N_cat_5},${N_cat_6}p\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_5.tsv"
			f_qsub.puts "sed -n -e \"${N_cat_6},\\$p\" #{dir_out_}/VIRSorter_global-phage-signal.parsed.tsv > #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_6.tsv"
	
			# conversion of contig names
			f_qsub.puts "grep \">\" #{query_fa_} | perl -pe \"s/^>//\" > #{query_org_id_}"
			f_qsub.puts "perl -pe \"s/\\./_/g\" #{query_org_id_} > #{query_cnv_id_}"
			f_qsub.puts "paste #{query_org_id_} #{query_cnv_id_} > #{query_id_table_}"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_1_2_3.tsv 1 | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_1_2_3.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_1_2.tsv 1   | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_1_2.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_1.tsv 1     | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_1.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_2.tsv 1     | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_2.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_3.tsv 1     | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_3.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_4.tsv 1     | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_4.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_5.tsv 1     | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_5.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{query_id_table_} 2 #{dir_out_}/VIRSorter_global-phage-signal.parsed.cat_6.tsv 1     | cut -f 2 > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_6.txt"
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_4.txt #{dir_out_}/VIRSorter_global-phage-signal.id.cat_5.txt > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_4_5.txt"
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_4.txt #{dir_out_}/VIRSorter_global-phage-signal.id.cat_5.txt #{dir_out_}/VIRSorter_global-phage-signal.id.cat_6.txt > #{dir_out_}/VIRSorter_global-phage-signal.id.cat_4_5_6.txt"
			# summary : contig - virsorter category
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_1.txt | sed -e \"s/$/\t1/\" > #{dir_out_}/contig.VIRSorter_global-phage-signal.txt"
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_2.txt | sed -e \"s/$/\t2/\" >> #{dir_out_}/contig.VIRSorter_global-phage-signal.txt"
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_3.txt | sed -e \"s/$/\t3/\" >> #{dir_out_}/contig.VIRSorter_global-phage-signal.txt"
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_4.txt | sed -e \"s/$/\t4/\" >> #{dir_out_}/contig.VIRSorter_global-phage-signal.txt"
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_5.txt | sed -e \"s/$/\t5/\" >> #{dir_out_}/contig.VIRSorter_global-phage-signal.txt"
			f_qsub.puts "cat #{dir_out_}/VIRSorter_global-phage-signal.id.cat_6.txt | sed -e \"s/$/\t6/\" >> #{dir_out_}/contig.VIRSorter_global-phage-signal.txt"
	
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
end

