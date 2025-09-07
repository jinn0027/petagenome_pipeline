
#--------------------------------------------------
# Program: 10_orf_coverage.rb
# Author: Yasumasa Kimura
#--------------------------------------------------
version = "1.0.0"

require 'fileutils'
require 'inifile'
require 'optparse'
h_params = ARGV.getopts('', 'mem:15.9')
mem = h_params["mem"]

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"

assembly_label = "SPAdes_meta_c"
l_thre = "5000.c_1500"


"---------------------"
bbmap_dir = ""
bedtools_dir = ""
samtools_dir = ""
"---------------------"


task = "bbmap"
predictor = "prodigal"

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
	# check settings
	abort "Abort: no 'project' in ./setting.cfg" if project == nil
	abort "Abort: no 'dir_project_' in ./setting.cfg" if dir_project_ == nil
	abort "Abort: 'dir_project_':#{dir_project_} does not exist." if ! Dir.exist?(dir_project_)
	abort "Abort: no 'merge_label' in ./setting.cfg" if merge_label == nil
	abort "Abort: no 'samples' in ./setting.cfg" if samples == nil
	a_samples = samples.gsub(" ","").split(',')

	# database (contig sequences)
	db_label = "#{merge_label}_#{assembly_label}"
	dir_coverage_ = "#{dir_project_}/d10_contig_coverage/#{merge_label}"
	dir_orf_ = "#{dir_project_}/d12_orf_prediction/#{merge_label}"

	# log directory
	dir_log_ = "#{dir_project_}/d10_contig_coverage/log"
	FileUtils.mkdir_p(dir_log_) unless Dir.exist?(dir_log_)
	# move to log directory
	Dir.chdir(dir_log_)

	# positions of ORFs
	orf_position = "#{dir_orf_}/#{predictor}.#{db_label}.#{l_thre}.extended.predicted.ffn.position"
	if ! File.exist?("#{orf_position}.bed") then
		puts "make #{orf_position}.bed"
		abort "Abort: 'orf_position':#{orf_position}.txt does not exist." if ! File.exist?("#{orf_position}.txt")
		system("awk -F '\\t' '{OFS=\"\\t\"} { if($6==1) print $3,$4-1,$5,$1,$2,\"+\"; else print $3,$4-1,$5,$1,$2,\"-\"; }' #{orf_position}.txt > #{orf_position}.bed")
	end

	#------------------------------
	# qusb for read mapping
	#------------------------------
	for q in 0..(a_samples.length - 1) do
		sample = a_samples[q]
		# coverages of contigs
		read_cov = "#{dir_coverage_}/#{sample}/#{task}-#{db_label}.#{l_thre}.#{sample}"
		out_orf = "#{dir_coverage_}/#{sample}/#{task}-#{db_label}.#{l_thre}.#{sample}.orf"
		abort "Abort: 'read_cov':#{read_cov}.sort.bam does not exist." if ! File.exist?("#{read_cov}.sort.bam")

		qsub_ = "p10_orf.#{task}.#{sample}.#{project}.#{db_label}.#{l_thre}.sh"
		f_qsub = open(qsub_, "w")
		f_qsub.puts "\#!/bin/bash"
		f_qsub.puts "\#$ -S /bin/bash"
		f_qsub.puts "\#$ -l s_vmem=#{mem}G,mem_req=#{mem}G"
		f_qsub.puts "\#$ -cwd"
		f_qsub.puts "\# version #{version}"
		f_qsub.puts "BEDTOOLS_=\"#{bedtools_dir}/bin/bedtools\""
		f_qsub.puts "SCRIPT_COV_=\"#{dir_scripts_}/get_covered_region_bed.from_depth.rb\""

		# Bed file of regions covered by reads
		f_qsub.puts "if [ -f #{read_cov}.depth.gz ]; then"
		f_qsub.puts "  zcat #{read_cov}.depth.gz | ruby ${SCRIPT_COV_} > #{read_cov}.region.bed"
		f_qsub.puts "else"
		f_qsub.puts "  #{samtools_dir}/samtools depth -m 10000000 -a #{read_cov}.sort.bam | ruby ${SCRIPT_COV_} > #{read_cov}.region.bed"
		f_qsub.puts "fi"
		# Intersections of ORF reagions and regions covered by reads
		f_qsub.puts "${BEDTOOLS_} intersect -a #{orf_position}.bed -b #{read_cov}.region.bed -wo > #{out_orf}.intersect.bed"
		# ORF length covered by reads
		f_qsub.puts "echo \"orf\tlength\tcovered_length\" > #{out_orf}.covered_length.txt"
		f_qsub.puts "awk -F '\\t' '{OFS=\"\\t\"} { a[$4]+=$10; l[$4]=$5; } END{for(k in a) print k,l[k],a[k]/l[k];}' #{out_orf}.intersect.bed | sort -k1,1 >> #{out_orf}.covered_length.txt"

		# ORF coverage by reads
		f_qsub.puts "#{samtools_dir}/samtools index #{read_cov}.sort.bam"
		f_qsub.puts "#{samtools_dir}/samtools bedcov #{orf_position}.bed #{read_cov}.sort.bam | awk -F '\\t' '{OFS=\"\\t\"} {if($7!=0) print $0}' > #{out_orf}.coverage.txt"

		f_qsub.close
		system("qsub #{qsub_}")
	end
end

