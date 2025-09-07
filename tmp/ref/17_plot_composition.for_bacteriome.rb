
#--------------------------------------------------
# Program: 17_plot_composition.for_bacteriome.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'host:human', 'max_taxa:29', 'base_size:8')
host = h_params["host"] # host organism
max_taxa = h_params["max_taxa"] # max number of taxa to show in plot
base_size = h_params["base_size"] # base_size in ggplot2

dir_home_ = ""
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
#dir_R_ = "/usr/local/package/r/3.4.1/bin"
profile_ = "#{dir_home_}/pipeline/profile"

# annotation info
db_label = "metaphlan2_v20"
tax_ann = "read.#{db_label}"

script_composition_ = "#{dir_scripts_}/plot_composition.bacteriome.R"

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
	org_merge_label = inifile[section]['org_merge_label']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	org_merge_label = merge_label if org_merge_label == nil
	a_samples = samples.gsub(" ","").split(',')
	samples = a_samples.join(",")

	# input directories
	dir_in_ = "#{dir_project_}/d08_bacterial_taxonomic_profile/#{org_merge_label}"
	# output directory
	dir_out_ = "#{dir_project_}/d17_plot_composition/#{merge_label}"
	FileUtils.mkdir_p(dir_out_) unless FileTest.exist?(dir_out_)

	# input
	in_ = "#{dir_in_}/profile.#{org_merge_label}.#{db_label}.parsed.txt"
	abort "Abort: 'in_':#{in_} does not exist." if ! File.exist?(in_)

	# output
	out_pdf_ = "#{dir_out_}/composition.#{merge_label}.#{tax_ann}.pdf"

	puts "start #{merge_label}"
	system("source #{profile_}; R --vanilla --slave --args #{in_} #{dir_out_} #{out_pdf_} #{host} #{samples} #{max_taxa} #{base_size} < #{script_composition_}")
end
