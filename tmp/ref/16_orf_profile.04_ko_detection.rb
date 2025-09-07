
#--------------------------------------------------
# Program: 16_orf_profile.04_ko_detection.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'date:20180916', 'kegg_db:kegg_prokaryotes', 'thresholds:0.85')
date = h_params["date"] # date of KEGG data
db_label = h_params["kegg_db"] # kegg db used
thresholds = h_params["thresholds"] # thresholds for ORF coverage

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"

assembly_label = "SPAdes_meta_c"
l_thre = "5000.c_1500"
predictor = "prodigal"
coverage_type = "covered_length"

# scripts
script_pathway_ko_ = "#{dir_scripts_}/kegg_pathway_ko_detection.orf_coverage.R"
script_module_ko_ = "#{dir_scripts_}/kegg_module_ko_detection.orf_coverage.R"

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
	# input file
	dir_16_ = "#{dir_project_}/d16_orf_profile_table/#{merge_label}"
	dir_out_ = dir_16_

	#------------------------------
	# Pathway - KO
	puts "start Pathway - KO"
	in_pathway_ = "#{dir_16_}/orf.#{db_label}.pathway.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
	# run script
	system("R --vanilla --slave --args #{in_pathway_} #{dir_out_} #{samples} #{date} #{db_label} #{thresholds} <  #{script_pathway_ko_}")

	#------------------------------
	# Module - KO
	puts "start Module - KO"
	in_module_ = "#{dir_16_}/orf.#{db_label}.module.ko.#{coverage_type}.#{predictor}.#{query_label}.#{l_thre}.txt"
	# run script
	system("R --vanilla --slave --args #{in_module_} #{dir_out_} #{samples} #{date} #{db_label} #{thresholds} <  #{script_module_ko_}")

end
