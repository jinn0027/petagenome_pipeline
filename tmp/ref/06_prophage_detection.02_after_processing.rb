
#--------------------------------------------------
# Program: 06_prophage_detection.02_after_processing.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'n_threads:3', 'l:', 'q:', 'assembly_label_each:SPAdes_meta', 'assembly_label:SPAdes_meta_c', 'l_thre_each:5000', 'l_thre:5000.c_1500')
n_threads = h_params["n_threads"]
l_opt = h_params["l"] # type of queue such as ljob
q_opt = h_params["q"] # qruop queue such as sysimm.q, hgc0951.q
assembly_label_each = h_params["assembly_label_each"]
assembly_label = h_params["assembly_label"]
l_thre_each = h_params["l_thre_each"]
l_thre = h_params["l_thre"]

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_data_ = "#{dir_home_}/data"
profile_ = "#{dir_home_}/pipeline/profile"

blast_dir = ""
cdhit_dir = ""
bedtools_dir = ""
seqkit_dir = ""

prefix = "bPContig"

# virsorter database
dir_data_vs_ = "#{dir_data_}/virsorter-data"
a_vs_db_ids = ["1", "2"]
h_vs_db = {"1" => "refseqdb", "2" => "viromedb"}

# use microbiome mode
opt_virome = ""
opt_label = ""

task = "megablast"
n_al = "1"
pi_thre = "95"
e_thre = "1e-100"

category = "cat_4_5"
#category = "cat_4_5_6"

# loading setting file
abort "Abort: './setting.cfg' does not exist." if ! File.exist?('./setting.cfg')
inifile = IniFile.load('./setting.cfg')

inifile.each_section do |section|
	puts "start #{section}"
	# store settings
	project = inifile[section]['project']
	dir_project_ = inifile[section]['dir_project_']
	merge_label = inifile[section]['merge_label']
	org_merge_label = inifile[section]['org_merge_label']
	samples = inifile[section]['samples']
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')
	# org_merge_label is used when merge_label is different from original
	if org_merge_label == nil then
		org_merge_label = merge_label
	end

	dir_04_ = "#{dir_project_}/d04_pool_contigs/#{org_merge_label}"
	dir_05_ = "#{dir_project_}/d05_circular_contigs/#{org_merge_label}"
	dir_09_ = "#{dir_project_}/d09_bacterial_contig_taxonomy/#{org_merge_label}"
	# output directory
	dir_out_base_ = "#{dir_project_}/d06_prophage_detection"
	dir_out_ = "#{dir_out_base_}/#{merge_label}_proph"
	FileUtils.mkdir_p(dir_out_) unless Dir.exist?(dir_out_)
	# log directory
	dir_log_ = "#{dir_out_base_}/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)

	# move to log directory
	Dir.chdir(dir_log_)

	# input
	contig_ = "#{dir_05_}/contig.updated.#{org_merge_label}_#{assembly_label}.#{l_thre}.fa"
	table_name_05_ = "#{dir_05_}/name.all.#{org_merge_label}_#{assembly_label}.#{l_thre}.txt"
	bac_taxa_ = "#{dir_09_}/taxonomy.ppsp.#{org_merge_label}_#{assembly_label}.#{l_thre}.length.uc_.parsed.txt"

	# output
	out_name_ = "#{dir_out_}/name.#{merge_label}.VIRSorter_global-phage-signal.id.#{category}.txt"
	out_count = "#{dir_out_}/count.#{merge_label}.VIRSorter_global-phage-signal.id.#{category}"
	out_count_merged = "#{dir_out_}/count.merged_prophage.virsorter.#{category}.#{merge_label}_#{assembly_label}.#{l_thre}"

	prophage = "#{dir_out_}/prophage.#{merge_label}.VIRSorter_global-phage-signal.id.#{category}"

	contig_proph = "#{dir_out_}/contig.updated.#{merge_label}_#{assembly_label}.#{l_thre}.virsorter.#{category}"
	blast_out = "#{dir_out_}/#{task}-#{merge_label}_#{assembly_label}.#{l_thre}.virsorter.#{category}.pi_#{pi_thre}"
	merged_prophage = "#{dir_out_}/merged_prophage.virsorter.#{category}.#{merge_label}_proph_#{assembly_label}.#{l_thre}"

	# check
	abort "Abort: 'contig_':#{contig_} does not exist." if ! File.exist?(contig_)

	# qsub file
	qsub_ = "p06_01_after.virsorter#{opt_label}.#{category}.#{project}.#{merge_label}.#{l_thre}.bash"
	f_qsub = open(qsub_, "w")
	f_qsub.puts "\#!/bin/bash"
	f_qsub.puts "\#$ -S /bin/bash"
	f_qsub.puts "\#$ -l s_vmem=5.3G,mem_req=5.3G -pe def_slot #{n_threads}"
	f_qsub.puts "\#$ -cwd"
	f_qsub.puts "\# version #{version}"
	f_qsub.puts "source #{profile_}"
	f_qsub.puts "SCRIPT_JOIN_=#{dir_scripts_}/join_with_tab.rb"
	f_qsub.puts "MAKEBLASTDB_=\"#{blast_dir}/makeblastdb\""
	f_qsub.puts "BLASTN_=\"#{blast_dir}/bin/blastn\""
	f_qsub.puts "CDHIT_=\"#{cdhit_dir}/bin/cd-hit-est\""
	f_qsub.puts "BEDTOOLS_=\"#{bedtools_dir}/bin/bedtools\""
	f_qsub.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""
	f_qsub.puts "SCRIPT_RENAME_=\"#{dir_scripts_}/filter_contig.rename.py\""

	#------------------------------
	# extract contigs containing prophages
	#------------------------------
	File.unlink(out_name_) if File.exist?(out_name_)
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		query_label = "#{sample}_#{assembly_label_each}"

		# input
		table_name_ = "#{dir_04_}/name.#{query_label}.#{org_merge_label}.#{l_thre_each}.txt"
		abort "Abort: 'table_name_':#{table_name_} does not exist." if ! File.exist?(table_name_)

		for vs_db_id in a_vs_db_ids
			vs_db_label = h_vs_db[vs_db_id]
			# virsorter output directory
			dir_each_ = "#{dir_out_base_}/#{sample}/virsorter#{opt_label}-#{vs_db_label}.#{query_label}.#{l_thre_each}"
			# contig names in category 4 & 5
			f_qsub.puts "echo \"start #{sample} #{vs_db_label}\""
			f_qsub.puts "awk -F '\\t' '{OFS=\"\\t\"} {a[$1]+=1} END{for(k in a)print k,a[k];}' #{dir_each_}/VIRSorter_global-phage-signal.id.#{category}.txt > #{dir_each_}/count.VIRSorter_global-phage-signal.id.#{category}.txt"
			f_qsub.puts "ruby ${SCRIPT_JOIN_} #{dir_each_}/count.VIRSorter_global-phage-signal.id.#{category}.txt 1 #{table_name_} 1 > #{dir_each_}/name.#{merge_label}.VIRSorter_global-phage-signal.id.#{category}.txt"
			f_qsub.puts "cat #{dir_each_}/name.#{merge_label}.VIRSorter_global-phage-signal.id.#{category}.txt >> #{out_name_}"
		end
	end
	f_qsub.puts "awk -F '\\t' '{OFS=\"\\t\"} {if(a[$4]<$2) a[$4]=$2} END{for(k in a)print k,a[k];}' #{out_name_} | ruby ${SCRIPT_JOIN_}  #{table_name_05_} 1 - 1 | cut -f 2,3 > #{out_count}.txt"


	#------------------------------
	# extract prophage sequences
	#------------------------------
	File.unlink("#{prophage}.fa") if File.exist?("#{prophage}.fa")
	for l in 0..(a_samples.length - 1) do
		sample = a_samples[l]
		query_label = "#{sample}_#{assembly_label_each}"

		for vs_db_id in a_vs_db_ids
			vs_db_label = h_vs_db[vs_db_id]
			# virsorter output directory
			dir_each_ = "#{dir_out_base_}/#{sample}/virsorter#{opt_label}-#{vs_db_label}.#{query_label}.#{l_thre_each}"
			dir_seq_ = "#{dir_each_}/Predicted_viral_sequences"
			query_id_table_ = "#{dir_each_}/input_sequences.id.table.txt"
			# contig names in category 4 & 5
			f_qsub.puts "echo \"start #{sample}\""

			# convert fasta id of predicted prophage sequences
			f_qsub.print "cat #{dir_seq_}/VIRSorter_prophages_cat-4.fasta | perl -pe \"s/VIRSorter_//\" | perl -pe \"s/-circular//\" "
			f_qsub.print " | #{dir_tools_}/bin/seqkit replace --pattern \"^(\\w+)(_gene_\\d+_gene_\\d+-)\" --replacement \"{kv}:${2}\" "
			f_qsub.print "   --kv-file <(awk '{OFS=\"\\t\"} {print $2,$1}' #{query_id_table_}) > #{dir_seq_}/VIRSorter_prophages_cat-4.fa\n"
			f_qsub.print "cat #{dir_seq_}/VIRSorter_prophages_cat-5.fasta | perl -pe \"s/VIRSorter_//\" | perl -pe \"s/-circular//\" "
			f_qsub.print " | #{dir_tools_}/bin/seqkit replace --pattern \"^(\\w+)(_gene_\\d+_gene_\\d+-)\" --replacement \"{kv}:${2}\" "
			f_qsub.print "   --kv-file <(awk '{OFS=\"\\t\"} {print $2,$1}' #{query_id_table_}) > #{dir_seq_}/VIRSorter_prophages_cat-5.fa\n"
			f_qsub.print "cat #{dir_seq_}/VIRSorter_prophages_cat-6.fasta | perl -pe \"s/VIRSorter_//\" | perl -pe \"s/-circular//\" "
			f_qsub.print " | #{dir_tools_}/bin/seqkit replace --pattern \"^(\\w+)(_gene_\\d+_gene_\\d+-)\" --replacement \"{kv}:${2}\" "
			f_qsub.print "   --kv-file <(awk '{OFS=\"\\t\"} {print $2,$1}' #{query_id_table_}) > #{dir_seq_}/VIRSorter_prophages_cat-6.fa\n"

			# concatenate predicted prophage sequences
			if category == "cat_4_5" then
				f_qsub.puts "cat #{dir_seq_}/VIRSorter_prophages_cat-4.fa #{dir_seq_}/VIRSorter_prophages_cat-5.fa >> #{prophage}.fa"
			elsif category == "cat_4_5_6" then
				f_qsub.puts "cat #{dir_seq_}/VIRSorter_prophages_cat-4.fa #{dir_seq_}/VIRSorter_prophages_cat-5.fa #{dir_seq_}/VIRSorter_prophages_cat-6.fa >> #{prophage}.fa"
			end
		end
	end

	# make blastdb of contigs where prophages were
	f_qsub.puts "cat #{contig_} | #{seqkit_dir}/bin/seqkit grep -n --pattern-file <(cut -f 1 #{out_count}.txt) -o #{contig_proph}.fa"
	f_qsub.puts "${MAKEBLASTDB_} -in #{contig_proph}.fa -out #{contig_proph} -dbtype nucl -parse_seqids"
	f_qsub.puts "${BLASTN_} -task #{task} -num_threads #{n_threads} -query #{prophage}.fa -db #{contig_proph} -perc_identity #{pi_thre} -evalue #{e_thre} -outfmt 6 -num_alignments #{n_al} -out #{blast_out}.txt"

	# get prophage length
	f_qsub.puts "python2.7 ${SCRIPT_LENGTH_} #{contig_proph}.fa > #{contig_proph}.length.txt"
	f_qsub.puts "python2.7 ${SCRIPT_LENGTH_} #{prophage}.fa > #{prophage}.length.txt"
	f_qsub.puts "perl -pe \"s/:/\\t/\" #{prophage}.length.txt | cut -f 2 | awk -F '-' '{if($1<$2) print $2-$1; else print $1-$2;}' | paste #{prophage}.length.txt - > #{prophage}.length.proph_len.txt"
	f_qsub.puts "ruby ${SCRIPT_JOIN_} #{blast_out}.txt 1 #{prophage}.length.proph_len.txt 1 > #{blast_out}.with_proph_len.txt"

	#*** make each virsorter-predicted region (since alignments can be separated by a gap and a predicted region sometimese exceeds a contig region) ***
	# bed file of pooled_contig::individual_prophage
	f_qsub.puts "cat #{blast_out}.with_proph_len.txt | awk -F '\\t' '{OFS=\"\\t\"} {if($9<$10) print $2\"::\"$1,$9-1,$10,$1,$14,\"+\",$13,$7,$8; else print $2\"::\"$1,$10-1,$9,$1,$14,\"-\",$13,$7,$8;}' | sort -k1,1 -k2,2n > #{blast_out}.bed"
	# merged alignments separated by a gap (Ns)
	f_qsub.puts "${BEDTOOLS_} merge -i #{blast_out}.bed -d 150 -s -c 4,5,6,7,8,9 -o distinct,distinct,distinct,distinct,min,max | perl -pe \"s/::/\\t/\" | cut -f 1,3- | sort -k1,1 -k2,2n > #{blast_out}.merged.tmp1.bed"
	# extend a region if the whole sequence aligned and the predicted region is longer than the sequence
	f_qsub.puts "cat #{blast_out}.merged.tmp1.bed | awk -F '\\t' '{OFS=\"\\t\"} {if($5>$7 && $8==1 && $9==$7) {if($6==\"+\") $3+=$5-$7; else $2-=$5-$7;}; if($2<0) $2=0; print $0;}' | sort -k1,1 -k2,2n > #{blast_out}.merged.tmp2.bed"
	#*** merge ovrelapped virsorter-predicted regions ***
	f_qsub.puts "${BEDTOOLS_} merge -i #{blast_out}.merged.tmp2.bed -d 0 -c 4,5,6 -o distinct,distinct,distinct > #{blast_out}.merged.tmp3.bed"
	# correct if a region exceed the length of the aligned pooled contig
	f_qsub.puts "ruby ${SCRIPT_JOIN_} #{blast_out}.merged.tmp3.bed 1 #{contig_proph}.length.txt 1 | awk -F '\\t' '{OFS=\"\\t\"} {if($3>$7) $3=$7; print $0}' > #{blast_out}.merged.tmp4.bed"
	#*** extract prophages > 3kb ***
	f_qsub.puts "awk -F '\\t' '{OFS=\"\\t\"} {if($3 - $2 > 3000) print $0}' #{blast_out}.merged.tmp4.bed | cut -f 1-4,6 > #{blast_out}.merged.bed"

	# get count of prophages at merged contigs
	f_qsub.puts "awk -F '\\t' '{OFS=\"\\t\"} {a[$1]+=1} END{for(k in a)print k,a[k];}' #{blast_out}.merged.bed > #{out_count_merged}.txt"
	f_qsub.puts "ruby ${SCRIPT_JOIN_} #{out_count_merged}.txt 1 #{bac_taxa_} 1 > #{out_count_merged}.bacterial_taxonomy.txt"

	# get merged prophage fasta sequence from merged contigs
	f_qsub.puts "${BEDTOOLS_} getfasta -fi #{contig_proph}.fa -bed #{blast_out}.merged.bed -fo #{merged_prophage}.original.fa"
	f_qsub.puts "python2.7 ${SCRIPT_RENAME_} -min 1500 --rename --prefix #{prefix}.#{project}.#{merge_label}_proph.l. --table #{merged_prophage}.name.txt #{merged_prophage}.original.fa > #{merged_prophage}.fa"
	# get length
	f_qsub.puts "python2.7 ${SCRIPT_LENGTH_} -t fasta #{merged_prophage}.fa > #{merged_prophage}.length.txt"

	f_qsub.close
	if q_opt != nil then
		system("qsub -q #{q_opt} #{qsub_}")
	elsif l_opt != nil then
		system("qsub -l #{l_opt} #{qsub_}")
	else
		system("qsub #{qsub_}")
	end
	#system("bash #{qsub_}")
end

