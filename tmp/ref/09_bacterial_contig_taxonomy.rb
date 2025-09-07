
#--------------------------------------------------
# Program: 09_bacterial_contig_taxonomy.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'mem:32', 'l:', 'q:', 'maxLeafClades:500', 'minPercentInLeaf:0.05', 'db_date:20181130', 'assembly_label:SPAdes_meta_c')
mem = h_params["mem"] # memory
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q
db_date = h_params["db_date"] # date of RefSeq data
assembly_label = h_params["assembly_label"]


# basic settings in PyloPythiaS+
rankIdCut = "6" # 6:species level
maxLeafClades = h_params["maxLeafClades"] # "500"
minPercentInLeaf = h_params["minPercentInLeaf"] # "0.05"
minSeqLen = "1000"

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_configs_ = "#{dir_home_}/pipeline/configs"
dir_data_ = "#{dir_home_}/data"
profile_ = "#{dir_home_}/pipeline/profile"

"-----------------"
ppsp_dir = ""
"-----------------"

l_thre = "5000.c_1500"

tax_map_ = "./ppsp/root_map_all.txt"
abort "Abort: 'tax_map_':#{tax_map_} does not exist." if ! File.exist?(tax_map_)

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

	# input directory
	dir_in_ = "#{dir_project_}/d05_circular_contigs/#{merge_label}"
	# output directory
	dir_out_ = "#{dir_project_}/d09_bacterial_contig_taxonomy/#{merge_label}"
	# log directory
	dir_log_ = "#{dir_project_}/d09_bacterial_contig_taxonomy/log"

	query_label = "#{merge_label}_#{assembly_label}"
	# input
	query_fa_ = "#{dir_in_}/contig.updated.#{query_label}.#{l_thre}.fa"
	query_len_ = "#{dir_in_}/contig.updated.#{query_label}.#{l_thre}.length.txt"
	config_temp_ = "./ppsp/config_ppsp_vm_refNCBI201502.template.cfg"

	# output
	config_ = "#{dir_out_}/config_ppsp_vm_refNCBI201502.#{query_label}.#{l_thre}.cfg"
	out_ppsp_ = "#{dir_out_}/output/contig.updated.#{query_label}.#{l_thre}.fa.pOUT"
	out_tax = "#{dir_out_}/taxonomy.ppsp.#{query_label}.#{l_thre}"

	# check
	abort "Abort: 'query_fa_':#{query_fa_} does not exist." if ! File.exist?(query_fa_)
	abort "Abort: 'query_len_':#{query_len_} does not exist." if ! File.exist?(query_len_)
	if Dir.exist?(dir_out_) then
		abort "Abort: 'output directory':#{dir_out_} already exists."
	else
		FileUtils.mkdir_p(dir_out_)
	end

	# move to log directory
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	Dir.chdir(dir_log_)

	# generation of config file
	f_temp = open(config_temp_)
	text_config = f_temp.read
	f_temp.close
	# write settings
	text_config.sub!("[dir_out_]",dir_out_)
	text_config.sub!("[query_fa_]",query_fa_)
	text_config.sub!("[rankIdCut]",rankIdCut)
	text_config.sub!("[maxLeafClades]",maxLeafClades)
	text_config.sub!("[minPercentInLeaf]",minPercentInLeaf)
	text_config.sub!("[minSeqLen]",minSeqLen)
	f_config = open(config_, "w")
	f_config.puts text_config
	f_config.close

	# qsub file
	puts "start #{merge_label}"
	qsub_ = "p09_ppsp.#{project}.#{query_label}.#{l_thre}.sh"
	f_qsub = open(qsub_, "w")
	f_qsub.puts "\#!/bin/bash"
	f_qsub.puts "\#$ -S /bin/bash"
	f_qsub.puts "\#$ -l s_vmem=#{mem}G,mem_req=#{mem}G"
	f_qsub.puts "\#$ -cwd"
	f_qsub.puts "\# version #{version}"
	f_qsub.puts "source #{profile_}"
	f_qsub.puts "PPSP_=#{ppsp_dir}/pps/tools/ppsp"
	f_qsub.puts "SCRIPT_JOIN_=#{dir_scripts_}/join_with_tab.option.rb"

	# run ppsp
	f_qsub.puts "${PPSP_} -c #{config_} -n -g -o s16 mg -t -p c -r -s"

	# add taxonomy information
	f_qsub.puts "ruby ${SCRIPT_JOIN_} #{out_ppsp_} 2 #{tax_map_} 1  >  #{out_ppsp_}.with_lineage.txt"
	# parse taxonomy information
	f_qsub.puts "cat #{out_ppsp_}.with_lineage.txt | perl -pe \"s/\\d+://g\" | awk -F '\\t' '{OFS=\"\\t\"} {print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,\"ppsp\"}' > #{out_tax}.txt"

	# set uc_ (uncharacterized) to lower rank taxa
	f_qsub.print "cat #{out_tax}.txt "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} { if($3==\"superkingdom\") print $1,$2,$3,$4,$5,\"uc_\"$5,\"uc_\"$5,\"uc_\"$5,\"uc_\"$5,\"uc_\"$5,\"uc_\"$5,\"uc_\"$5,$13; else print $0; }' "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} { if($3==\"phylum\")       print $1,$2,$3,$4,$5,$6,\"uc_\"$6,\"uc_\"$6,\"uc_\"$6,\"uc_\"$6,\"uc_\"$6,\"uc_\"$6,$13; else print $0;  }' "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} { if($3==\"class\")        print $1,$2,$3,$4,$5,$6,$7,\"uc_\"$7,\"uc_\"$7,\"uc_\"$7,\"uc_\"$7,\"uc_\"$7,$13; else print $0;  }' "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} { if($3==\"order\")        print $1,$2,$3,$4,$5,$6,$7,$8,\"uc_\"$8,\"uc_\"$8,\"uc_\"$8,\"uc_\"$8,$13; else print $0;  }' "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} { if($3==\"family\")       print $1,$2,$3,$4,$5,$6,$7,$8,$9,\"uc_\"$9,\"uc_\"$9,\"uc_\"$9,$13; else print $0;  }' "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} { if($3==\"subfamily\")    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,\"uc_\"$10,\"uc_\"$10,$13; else print $0;  }' "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} { if($3==\"genus\")        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,\"uc_\"$11,$13; else print $0;  }' "
	f_qsub.print "> #{out_tax}.uc_.txt \n"

	# put prefix
	f_qsub.puts "cat #{out_tax}.uc_.txt | awk -F '\\t' '{OFS=\"\\t\"} {print $1,$2,$3,$4,\"k__\"$5,\"p__\"$6,\"c__\"$7,\"o__\"$8,\"f__\"$9,\"sf__\"$10,\"g__\"$11,\"s__\"$12,$13}' > #{out_tax}.uc_.parsed.txt"

	# put length and remove kingdom for final output
	f_qsub.print "ruby ${SCRIPT_JOIN_} #{query_len_} 1 #{out_tax}.uc_.parsed.txt 1 "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} {print $1,$2,$5,$4,$7,$8,$9,$10,$11,$12,$13,$14}' "
	f_qsub.print "> #{out_tax}.length.uc_.parsed.txt \n"

	# add unclassified contigs
	a_unknown = ["unknown","no_rank","p__unknown","c__unknown","o__unknown","f__unknown","sf__unknown","g__unknown","s__unknown","ppsp"]
	f_qsub.print "ruby ${SCRIPT_JOIN_} #{query_len_} 1 #{out_tax}.uc_.parsed.txt 1 \"-v 1\" "
	f_qsub.print "| awk -F '\\t' '{OFS=\"\\t\"} {print $1,$2,\"#{a_unknown.join("\",\"")}\"}' "
	f_qsub.print ">> #{out_tax}.length.uc_.parsed.txt \n"

	# to save disk space
	f_qsub.puts "tar -czf  #{dir_out_}/working/projectDir.tar.gz  #{dir_out_}/working/projectDir  --remove-files"

	f_qsub.close
	if q_opt != nil then
		system("qsub -q #{q_opt} #{qsub_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} #{qsub_}")
	else
		system("qsub #{qsub_}")
	end
end
