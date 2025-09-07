
#--------------------------------------------------
# Program: 16_contig_profile_table.01_coverage.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'assembly_label:SPAdes_meta_c')
assembly_label = h_params["assembly_label"]

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"

# contig info
l_thre = "5000.c_1500"
task = "bbmap"

script_join_ = "#{dir_scripts_}/join_with_tab.option.rb"
script_normalize_ ="#{dir_scripts_}/normalize.contig_profile_table.R"

# loading setting file
abort "Abort: ./setting.cfg does not exist." if ! File.exist?('./setting.cfg')
inifile = IniFile.load('./setting.cfg')

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

	# check project type
	if type != "bacteriome" then
		abort "Abort: this script is for virome. 'type' in ./setting.cfg must be 'bacteriome'."
	end

	puts "start #{merge_label}"
	db_label = "#{merge_label}_#{assembly_label}"
	# input directories
	dir_05_ = "#{dir_project_}/d05_circular_contigs/#{merge_label}"
	dir_10_ = "#{dir_project_}/d10_contig_coverage/#{merge_label}"

	# output directory
	dir_out_ = "#{dir_project_}/d16_contig_profile_table/#{merge_label}"

	# input
	in_contig_len_ = "#{dir_05_}/contig.updated.#{db_label}.#{l_thre}.length.txt"

	# check
	abort "Abort: 'in_contig_len_':#{in_contig_len_} does not exist." if ! File.exist?(in_contig_len_)
	if Dir.exist?(dir_out_) then
		abort "Abort: 'output directory':#{dir_out_} already exists."
	else
		FileUtils.mkdir_p(dir_out_)
	end

	# read coverages of contigs
	a_cov_thre = ["", ".0.75"]

	for c in 0..(a_cov_thre.length - 1) do
		cov_thre = a_cov_thre[c]
		puts "coverage threshold : #{cov_thre}"
		# output
		label_out = "coverage#{cov_thre}.contig.updated.#{db_label}.#{l_thre}"
		out_coverage_ = "#{dir_out_}/table.#{label_out}.txt"

		# temporary files
		temp_merged_ = "#{dir_out_}/temp.merged.coverage.txt"
		temp_coverage_ = "#{dir_out_}/temp.coverage.txt"
		header_ = "#{dir_out_}/temp.header.txt"
		system("echo \"contig\tlength\t#{a_samples.join("\t")}\" > #{header_}")
		# initial output file
		system("cp #{in_contig_len_} #{out_coverage_}")

		#------------------------------
		# merge read coverage files
		#------------------------------
		for q in 0..(a_samples.length - 1) do
			sample = a_samples[q]
			puts "start #{sample}"
			# input file of read coverages of contigs
			in_coverage_ = "#{dir_10_}/#{sample}/#{task}-#{db_label}.#{l_thre}.#{sample}.coverage#{cov_thre}.txt"
			abort "Abort: 'in_coverage_':#{in_coverage_} does not exist." if ! File.exist?(in_coverage_)
			system("cut -f 1,3 #{in_coverage_} > #{temp_coverage_}")

			# add coverages of current sample to final output file
			system("ruby #{script_join_}  #{out_coverage_} 1  #{temp_coverage_} 1  > #{temp_merged_}")
			# concatenate missing contig as count 0
			system("ruby #{script_join_}  #{out_coverage_} 1 #{temp_coverage_} 1 \"-v 1\" | grep -v contig[[:space:]]length | sed -e \"s/$/\t0/\" | cat  >> #{temp_merged_}")
			# sort by contig id
			system("sort -k1,1 #{temp_merged_} > #{out_coverage_}")
		end
		# final result
		system("sort -k1,1 #{temp_merged_} | cat #{header_} - > #{out_coverage_}")
		# normalization
		system("R --vanilla --slave --args #{dir_out_} #{out_coverage_} #{label_out} <  #{script_normalize_}")

		system("rm #{dir_out_}/temp.*")
		puts "\n"
	end
end
