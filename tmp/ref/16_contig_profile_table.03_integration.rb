
#--------------------------------------------------
# Program: 16_contig_profile_table.03_integration.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('kcp', 'base_size:7', 'assembly_label:SPAdes_meta_c')
do_kegg = h_params["k"] # do kegg integration
do_crispr = h_params["c"] # do crispr integration
do_prophage = h_params["p"] # do prophage integration
base_size = h_params["base_size"] # max number of taxa to show in plot
assembly_label = h_params["assembly_label"]

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"

l_thre = "5000.c_1500"
predictor = "prodigal"
tax_ann = "ppsp"

# annotation info
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
	dir_16_ = "#{dir_project_}/d16_contig_profile_table/#{merge_label}"

	a_coverage_type = ["coverage.length","coverage.rate"]
	for c in 0..(a_coverage_type.length - 1) do

		coverage_type = a_coverage_type[c]
		puts coverage_type
		in_ann_ = "#{dir_16_}/contig_annotation.#{query_label}.#{l_thre}.#{tax_ann}.#{coverage_type}.updated.txt"

		if do_kegg then
			#--------------------------------------------------
			# KEGG
			puts "start KEGG"
			db_label = "kegg_prokaryotes"
			n_al = "1"
			e_thre = "1e-5"
			b_thre = "50"

			# input
			dir_14_ = "#{dir_project_}/d14_kegg_search/#{merge_label}"
			in_pathway_ = "#{dir_14_}/ghostmp-#{db_label}.#{predictor}.#{query_label}.#{l_thre}.n_#{n_al}.e_#{e_thre}.b_#{b_thre}.pathway.contig.txt"
			# output
			out_pathway_cov_ = "#{dir_16_}/contig_annotation.kegg_pathway.#{query_label}.#{l_thre}.#{tax_ann}.#{coverage_type}.updated.txt"
			out_pathway_ko_cov_ = "#{dir_16_}/contig_annotation.kegg_pathway.ko.#{query_label}.#{l_thre}.#{tax_ann}.#{coverage_type}.updated.txt"

			# header
			system("head -n 1 #{in_ann_} | perl -pe \"s/contig/contig\\tpathway\\torf/\" > #{out_pathway_cov_}")
			# body
			system("cut -f 1,3,4 #{in_pathway_} | sort | uniq > #{in_pathway_}.uniq")
			system("ruby #{script_join_} #{in_pathway_}.uniq 3 #{in_ann_} 1 >> #{out_pathway_cov_}")
			system("rm #{in_pathway_}.uniq")

			# kegg pathway & KO
			# header
			system("head -n 1 #{in_ann_} | perl -pe \"s/contig/contig\\tpathway\\tKO\\torf/\" > #{out_pathway_ko_cov_}")
			# body
			system("ruby #{script_join_} #{in_pathway_} 4 #{in_ann_} 1 >> #{out_pathway_ko_cov_}")

		end

		if do_crispr then
			#--------------------------------------------------
			# CRISPR spacer
			#spacer_type = "merged_spacer.crisprdetect"
			spacer_type = "each_spacer.crisprdetect"
			puts "start CRISPR spacer: #{spacer_type}"

			# input
			dir_07_ = "#{dir_project_}/d07_crispr_detection/#{merge_label}"
			in_count_crispr_ = "#{dir_07_}/count.#{spacer_type}.#{query_label}.#{l_thre}.txt"
			# output
			out_crispr_ann_ = "#{dir_16_}/contig_annotation.#{spacer_type}.#{query_label}.#{l_thre}.#{tax_ann}.#{coverage_type}.updated.txt"

			# header
			system("head -n 1 #{in_ann_} | perl -pe \"s/contig/contig\\tcount_spacer/\" > #{out_crispr_ann_}")
			# body
			system("ruby #{script_join_} #{in_count_crispr_} 1 #{in_ann_} 1 >> #{out_crispr_ann_}")
		end

		if do_prophage then
			#--------------------------------------------------
			# Prophage
			puts "start Prophage"
			prophage_type = "prophage.virsorter_cat_4_5"

			# input
			dir_06_ = "#{dir_project_}/d06_virsorter"
			in_count_prophage_ = "#{dir_06_}/count.#{merge_label}.VIRSorter_global-phage-signal.id.cat_4_5.rename.txt"
			# output
			out_prophage_ann_ = "#{dir_16_}/contig_annotation.#{prophage_type}.#{query_label}.#{l_thre}.#{tax_ann}.#{coverage_type}.updated.txt"

			# header
			system("head -n 1 #{in_ann_} | perl -pe \"s/contig/contig\\tcount_prophage/\" > #{out_prophage_ann_}")
			# body
			system("ruby #{script_join_} #{in_count_prophage_} 1 #{in_ann_} 1 >> #{out_prophage_ann_}")
		end

	end

end
