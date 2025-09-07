
#--------------------------------------------------
# Program: 16_orf_profile.03_integration.kegg.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'kegg_db:kegg_prokaryotes')
db_label = h_params["kegg_db"] # kegg db used

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"

assembly_label = "SPAdes_meta_c"
l_thre = "5000.c_1500"
predictor = "prodigal"

# scripts
script_join_ = "#{dir_scripts_}/join_with_tab.rb"

# loading setting file
inifile = IniFile.load('./setting.cfg')
abort "Abort: ./setting.cfg does not exist." if ! File.exist?('./setting.cfg')

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
	samples = a_samples.join(",")

	puts "start #{merge_label}"
	query_label = "#{merge_label}_#{assembly_label}"

	# contig annotation
	a_coverage_type = ["covered_length","coverage"]
	#a_coverage_type = ["coverage"]
	for c in 0..(a_coverage_type.length - 1) do
		coverage_type = a_coverage_type[c]
		puts "start #{coverage_type}"
		# input
		dir_16_ = "#{dir_project_}/d16_orf_profile_table/#{merge_label}"
		in_orf_prof_ = "#{dir_16_}/table.#{coverage_type}.orf.#{predictor}.#{query_label}.#{l_thre}.txt"
		in_orf_ann_ = "#{dir_16_}/orf_annotation.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
		system("cut -f 2- #{in_orf_ann_} > #{in_orf_ann_}.no_contig")
	
		#--------------------------------------------------
		# KEGG
		n_al = "1"
		e_thre = "1e-5"
		b_thre = "50"
	
		dir_14_ = "#{dir_project_}/d14_kegg_search/#{merge_label}"
	
		# KEGG Pathway & KO
		puts "start KEGG pathway & KO"
		# input
		in_pathway = "#{dir_14_}/ghostmp-#{db_label}.#{predictor}.#{query_label}.#{l_thre}.n_#{n_al}.e_#{e_thre}.b_#{b_thre}.pathway"
		# output
		out_pathway_ko_cov_ = "#{dir_16_}/orf.#{db_label}.pathway.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
		out_pathway_ko_cov_ann_ = "#{dir_16_}/orf_annotation.#{db_label}.pathway.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
	
		system("head -n 1 #{in_orf_prof_} | perl -pe \"s/orf/orf\\tpathway\\tKO/\" > #{out_pathway_ko_cov_}")
		system("ruby #{script_join_} #{in_pathway}.txt 3 #{in_orf_prof_} 1 >> #{out_pathway_ko_cov_}")
		system("head -n 1 #{in_orf_ann_}.no_contig | perl -pe \"s/orf/orf\\tpathway\\tKO/\" > #{out_pathway_ko_cov_ann_}")
		system("ruby #{script_join_} #{in_pathway}.txt 3 #{in_orf_ann_}.no_contig 1 >> #{out_pathway_ko_cov_ann_}")
	
	
		# KEGG Module & KO
		puts "start KEGG Module & KO"
		# input
		in_module = "#{dir_14_}/ghostmp-#{db_label}.#{predictor}.#{query_label}.#{l_thre}.n_#{n_al}.e_#{e_thre}.b_#{b_thre}.module"
		# output
		out_module_ko_cov_ = "#{dir_16_}/orf.#{db_label}.module.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
		out_module_ko_cov_ann_ = "#{dir_16_}/orf_annotation.#{db_label}.module.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
	
		system("head -n 1 #{in_orf_prof_} | perl -pe \"s/orf/orf\\tmodule\\tKO/\" > #{out_module_ko_cov_}")
		system("ruby #{script_join_} #{in_module}.txt 3 #{in_orf_prof_} 1 >> #{out_module_ko_cov_}")
		system("head -n 1 #{in_orf_ann_}.no_contig | perl -pe \"s/orf/orf\\tmodule\\tKO/\" > #{out_module_ko_cov_ann_}")
		system("ruby #{script_join_} #{in_module}.txt 3 #{in_orf_ann_}.no_contig 1 >> #{out_module_ko_cov_ann_}")
	
	
		# KEGG KO
		puts "start KEGG KO"
		# input
		in_ko = "#{dir_14_}/ghostmp-#{db_label}.#{predictor}.#{query_label}.#{l_thre}.n_#{n_al}.e_#{e_thre}.b_#{b_thre}.ko"
		# output
		out_ko_cov_ = "#{dir_16_}/orf.#{db_label}.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
		out_ko_cov_ann_ = "#{dir_16_}/orf_annotation.#{db_label}.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
	
		system("head -n 1 #{in_orf_prof_} | perl -pe \"s/orf/orf\\tKO/\" > #{out_ko_cov_}")
		system("ruby #{script_join_} #{in_ko}.txt 2 #{in_orf_prof_} 1 >> #{out_ko_cov_}")
		system("head -n 1 #{in_orf_ann_}.no_contig | perl -pe \"s/orf/orf\\tKO/\" > #{out_ko_cov_ann_}")
		system("ruby #{script_join_} #{in_ko}.txt 2 #{in_orf_ann_}.no_contig 1 >> #{out_ko_cov_ann_}")
	
		system("wait; rm #{in_orf_ann_}.no_contig")

	end
end
