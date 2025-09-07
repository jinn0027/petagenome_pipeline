
#--------------------------------------------------
# Program: 16_orf_profile.02_annotation.for_bacteriome.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"

# contig info
assembly_label = "SPAdes_meta_c"
l_thre = "5000.c_1500"
predictor = "prodigal"
tax_ann = "ppsp"

# loading setting file
inifile = IniFile.load('./setting.cfg')
abort "Abort: ./setting.cfg does not exist." if ! File.exist?('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	type = inifile[section]['type']
	dir_project_ = inifile[section]['dir_project_']
	samples = inifile[section]['samples']
	merge_label = inifile[section]['merge_label']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'type' in ./setting.cfg" if type == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')

	# checking project type
	if type != "bacteriome" then
		abort "Abort: this script is for virome. 'type' in ./setting.cfg must be 'bacteriome'."
	end

	puts "start #{merge_label}"
	# input directories
	dir_09_ = "#{dir_project_}/d09_bacterial_contig_taxonomy/#{merge_label}"
	dir_out_ = "#{dir_project_}/d16_orf_profile_table/#{merge_label}"

	# coverage type
	a_coverage_type = ["covered_length","coverage"]
	# header for merged contigs
	a_header = ["contig","orf","orf_length",a_samples,"contig_length","lca","rank","phylum_name","class_name","order_name","family_name","subfamily_name","genus_name","species_name","method"].flatten
	for c in 0..(a_coverage_type.length - 1) do
		coverage_type = a_coverage_type[c]
		puts "start #{coverage_type}"
		# input
		in_cov_ = "#{dir_out_}/table.#{coverage_type}.orf.#{predictor}.#{merge_label}_#{assembly_label}.#{l_thre}.txt"
		in_taxa_ = "#{dir_09_}/taxonomy.#{tax_ann}.#{merge_label}_#{assembly_label}.#{l_thre}.length.uc_.parsed.txt"
		# output
		out_ = "#{dir_out_}/orf_annotation.#{coverage_type}.#{predictor}.#{merge_label}_#{assembly_label}.#{l_thre}.txt"
		# check
		abort "Abort: 'in_cov_':#{in_cov_} does not exist." if ! File.exist?(in_cov_)
		abort "Abort: 'in_taxa_':#{in_taxa_} does not exist." if ! File.exist?(in_taxa_)

		# put contig name at the 1st column
		in_cov_contig_ = "#{dir_out_}/table.#{coverage_type}.orf.#{predictor}.#{merge_label}_#{assembly_label}.#{l_thre}.contig.txt"
		system("perl -pe \"s/^(.+\\.\\d{9})__/\\1\\t\\1__/\" #{in_cov_} > #{in_cov_contig_}")
		system("echo \"#{a_header.join("\t")}\" > #{out_}")
		system("ruby #{dir_scripts_}/join_with_tab.rb #{in_cov_contig_} 1 #{in_taxa_} 1 >> #{out_}")
		system("wait; rm #{in_cov_contig_}")
	end
end

