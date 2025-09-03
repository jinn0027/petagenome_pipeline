
require 'fileutils'
require 'bio'

dir_home_ = ""
dir_tools_ = "#{dir_home_}/tools"
dir_scripts_ = "#{dir_home_}/pipeline/scripts"
dir_R_ = "/usr/local/package/r/current3/bin"

require 'optparse'
h_params = ARGV.getopts('', 'dir_out_:', 'in_contig_:', 'n_threads:1', 'pi_self:100', 'al_self:50', 'pi_rd:98', 'qc_rd:98', 'len_l:5000', 'len_c:1500')
dir_out_ = h_params["dir_out_"] # directory for outputs
in_contig_ = h_params["in_contig_"] # fasta file of contigs
n_threads = h_params["n_threads"] # number of threads used in blast search
# parameters for circular contig formation
pi_self = h_params["pi_self"] # % identity for circular formation
al_self = h_params["al_self"] # alignment length for circular formation
len_l = h_params["len_l"] # minimum length of linear contigs
len_c = h_params["len_c"] # minimum length of circular contigs
# parameters for redundant contig detection
pi_rd = h_params["pi_rd"] # % identity for removal of redundancy
qc_rd = h_params["qc_rd"] # % query coverage for removal of redundancy

# sh file
sh_ = "sh.self.contig.#{$$}.sh"

task = "megablast"
e_thre = "1e-10"

dir_in_each_ = "#{dir_out_}/each_contig"
dir_out_each_ = "#{dir_out_}/each_blast"
dir_out_circular_ = "#{dir_out_}/each_circular"
FileUtils.mkdir_p(dir_in_each_) unless Dir.exist?(dir_in_each_)
FileUtils.mkdir_p(dir_out_) unless Dir.exist?(dir_out_)
FileUtils.mkdir_p(dir_out_each_) unless Dir.exist?(dir_out_each_)
FileUtils.mkdir_p(dir_out_circular_) unless Dir.exist?(dir_out_circular_)

# output files
out_blast_all_ = "#{dir_out_}/#{task}.all_vs_all.#{e_thre}.pi_#{pi_self}.txt"
out_blast_self_ = "#{dir_out_}/#{task}.self.#{e_thre}.pi_#{pi_self}.txt"
out_blast_circular_ = "#{dir_out_}/#{task}.circular.#{e_thre}.pi_#{pi_self}.al_#{al_self}.txt"
out_linear = "#{dir_out_}/contig.linear"
out_circular = "#{dir_out_}/contig.circular"
out_updated = "#{dir_out_}/contig.updated"

# initialize output files
if FileTest.exist?("#{out_blast_circular_}") then File.unlink("#{out_blast_circular_}") end
puts "start #{sh_}"
f_sh = open(sh_, "w")
f_sh.puts "OUT_=\"#{out_blast_circular_}\""
f_sh.puts "BLASTN_=\"#{dir_tools_}/ncbi-blast-2.7.1+/bin/blastn\""
f_sh.puts "MAKEBLASTDB_=\"#{dir_tools_}/ncbi-blast-2.7.1+/bin/makeblastdb\""
f_sh.puts "SCRIPT_LENGTH_=\"#{dir_scripts_}/get_sequence_length.py\""
f_sh.puts "SCRIPT_REDUNDANCY_=\"#{dir_scripts_}/extract_contig_redundancy.3.rb\""
f_sh.puts "SCRIPT_FILTER_FA_=\"#{dir_scripts_}/filter_fasta_by_id.py\""
f_sh.puts "SCRIPT_RENAME_=\"#{dir_scripts_}/filter_contig.rename.py\""
# make blastdb if it does not exist.
in_contig = "#{File.dirname(in_contig_)}/#{File.basename(in_contig_,'.*')}"
if ! File.exist?("#{in_contig}.nsq") and ! File.exist?("#{in_contig}.nal")
	f_sh.puts "${MAKEBLASTDB_} -in #{in_contig_}  -out #{in_contig}  -dbtype nucl  -parse_seqids"
end
f_sh.puts "${BLASTN_} -task #{task} -num_threads #{n_threads} -query #{in_contig_} -db #{in_contig} -evalue #{e_thre} -perc_identity #{pi_self} -outfmt 6 -num_alignments 5 -out #{out_blast_all_}"
f_sh.puts "awk -F \"\\t\" '{OFS=\"\\t\"}  { if ($1 == $2) print $0 }' #{out_blast_all_} > #{out_blast_self_}"

f_in = Bio::FlatFile.open(Bio::FastaFormat, in_contig_)
#------------------------------
# detection of circular contigs
#------------------------------
f_in.each_entry do |entry|
	# make fasta file of each contig
	query_each_ = "#{dir_in_each_}/#{entry.definition}.fa"
	f_each = open(query_each_, "w")
	f_each.puts ">#{entry.definition}"
	f_each.puts "#{entry.naseq.upcase}"
	f_each.close

	# megablast search of contig against itself
	out_blast_each_ = "#{dir_out_each_}/megablast.#{entry.definition}.txt"
	f_sh.puts "if [ $(grep -c #{entry.definition} #{out_blast_self_}) -gt 1 ]; then"
	f_sh.puts "  grep #{entry.definition} #{out_blast_self_} > #{out_blast_each_}"
	# get contig length (BLAST position is 1-based)
	f_sh.puts "  LENGTH=$(cat #{out_blast_each_} | head -n1 | cut -f 4)"
	# extract putative circular contigs
	f_sh.puts "  cat #{out_blast_each_} | tail -n +2 | grep -w ${LENGTH} | awk -F \"\\t\" '{OFS=\"\\t\"} {if($4>=#{al_self} && $9==1) print}' >> ${OUT_}"

	# get position of the overlap between 3' and 5' end
	f_sh.puts "  POS_DUP=$(cat #{out_blast_each_} | sort -t$'\\t' -k4 -n -r | tail -n +2 | grep -w ${LENGTH} | awk -F \"\\t\" '{OFS=\"\\t\"} {if($4>=#{al_self} && $9==1) print $7}' | head -n 1)"
	f_sh.puts "  if [ -n \"${POS_DUP}\" ]; then"
	f_sh.puts "    POS_END=`expr ${POS_DUP} - 1`"
	f_sh.puts "  else"
	f_sh.puts "    POS_END=0"
	f_sh.puts "  fi"

	# get circlular contig
	f_sh.puts "  if [ ${POS_END} -ge #{len_c} ]; then"
	f_sh.puts "    POS_END=`expr ${POS_DUP} - 1`"
	# write circular contig in fasta format
	f_sh.puts "    SEQ=$(sed -n 2p #{query_each_} | cut -c 1-${POS_END})"
	f_sh.puts "    echo \">#{entry.definition}\" > #{dir_out_circular_}/#{entry.definition}.cut.fa"
	f_sh.puts "    echo \"${SEQ}\" >> #{dir_out_circular_}/#{entry.definition}.cut.fa"
	# extend circular contig with minimum query length bp (the first 1,500 nucleotides are duplicated and added at the contig’s end.)
	#f_sh.puts "    SEQ_EX=$(sed -n 2p #{query_each_} | cut -c 1-#{len_c})"
	# extend circular contig with whole contig sequence (duplicated sequences are concatenated.)
	f_sh.puts "    SEQ_EX=$(sed -n 2p #{query_each_})"
	f_sh.puts "    echo \">#{entry.definition}\" > #{dir_out_circular_}/#{entry.definition}.extended.fa"
	f_sh.puts "    echo \"${SEQ}${SEQ_EX}\" >> #{dir_out_circular_}/#{entry.definition}.extended.fa"
	f_sh.puts "  fi"
	# removal of each files
	f_sh.puts "  rm #{out_blast_each_}"
	f_sh.puts "fi"
	f_sh.puts "rm #{query_each_}"
end

# concatenate circular contigs
f_sh.puts "cat  #{dir_out_circular_}/*.cut.fa > #{out_circular}.c_#{len_c}.all.fa"
f_sh.puts "cat  #{dir_out_circular_}/*.extended.fa > #{out_circular}.c_#{len_c}.all.extended.fa"
# get linear contigs
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_circular}.c_#{len_c}.all.fa > #{out_circular}.c_#{len_c}.all.length.txt"
f_sh.puts "wait"
f_sh.puts "python ${SCRIPT_FILTER_FA_} -f #{out_circular}.c_#{len_c}.all.length.txt #{in_contig_} > #{out_linear}.all.fa"

#------------------------------
# rename contigs:  circular:c  linear:l
#------------------------------
f_sh.puts "wait"
f_sh.puts "cat #{dir_out_circular_}/*.cut.fa | perl -pe \"s/n(\\.\\d+)$/c\\1/\" > #{out_circular}.c_#{len_c}.all.fa"
f_sh.puts "cat  #{dir_out_circular_}/*.extended.fa | perl -pe \"s/n(\\.\\d+)$/c\\1/\" > #{out_circular}.c_#{len_c}.all.extended.fa"
f_sh.puts "python ${SCRIPT_FILTER_FA_} -f #{out_circular}.c_#{len_c}.all.length.txt #{in_contig_} | perl -pe \"s/n(\\.\\d+)$/l\\1/\" > #{out_linear}.all.fa"

# get length of contigs
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_circular}.c_#{len_c}.all.fa > #{out_circular}.c_#{len_c}.all.length.txt"
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_circular}.c_#{len_c}.all.extended.fa > #{out_circular}.c_#{len_c}.all.extended.length.txt"
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_linear}.all.fa > #{out_linear}.all.length.txt"

# updated contigs
f_sh.puts "cat  #{out_circular}.c_#{len_c}.all.fa  #{out_linear}.all.fa > #{out_updated}.c_#{len_c}.all.fa"
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_updated}.c_#{len_c}.all.fa > #{out_updated}.c_#{len_c}.all.length.txt"

#------------------------------
# removal of redundancy in contigs
#------------------------------
out_blast_rd = "#{dir_out_}/#{task}.redundancy.updated.c_#{len_c}.#{e_thre}.pi_#{pi_rd}"
out_rd_info_ = "#{dir_out_}/info.redundancy.updated.c_#{len_c}.#{e_thre}.pi_#{pi_rd}.qc_#{qc_rd}.txt"
out_ex_contig_ = "#{dir_out_}/contig.redundant.c_#{len_c}.length.txt"

# blast search for extracting redundancy
f_sh.puts "${MAKEBLASTDB_} -in #{out_updated}.c_#{len_c}.all.fa  -out #{out_updated}.c_#{len_c}.all  -dbtype nucl  -parse_seqids"
f_sh.puts "${BLASTN_} -task #{task} -num_threads #{n_threads} -query #{out_updated}.c_#{len_c}.all.fa -db #{out_updated}.c_#{len_c}.all -evalue #{e_thre} -perc_identity #{pi_rd} -outfmt 6 -num_alignments 50 -out #{out_blast_rd}.all.txt"
f_sh.puts "awk -F \"\\t\" '{OFS=\"\\t\"}  { if ($1 != $2) print $0 }' #{out_blast_rd}.all.txt > #{out_blast_rd}.txt"

# extract redundant contigs
f_sh.puts "ruby ${SCRIPT_REDUNDANCY_}  -b #{out_blast_rd}.txt  -l #{out_updated}.c_#{len_c}.all.length.txt  -c #{qc_rd}  -d 6  -i #{out_rd_info_}  -o #{out_ex_contig_}"

# get deduplicated contigs
f_sh.puts "python ${SCRIPT_FILTER_FA_} -f #{out_ex_contig_} #{out_linear}.all.fa > #{out_linear}.fa"
f_sh.puts "python ${SCRIPT_FILTER_FA_} -f #{out_ex_contig_} #{out_circular}.c_#{len_c}.all.fa > #{out_circular}.c_#{len_c}.fa"
f_sh.puts "python ${SCRIPT_FILTER_FA_} -f #{out_ex_contig_} #{out_circular}.c_#{len_c}.all.extended.fa > #{out_circular}.c_#{len_c}.extended.fa"
f_sh.puts "python ${SCRIPT_FILTER_FA_} -f #{out_ex_contig_} #{out_updated}.c_#{len_c}.all.fa > #{out_updated}.c_#{len_c}.fa"

f_sh.puts "python ${SCRIPT_LENGTH_} #{out_linear}.fa > #{out_linear}.length.txt"
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_circular}.c_#{len_c}.fa > #{out_circular}.c_#{len_c}.length.txt"
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_circular}.c_#{len_c}.extended.fa > #{out_circular}.c_#{len_c}.extended.length.txt"
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_updated}.c_#{len_c}.fa > #{out_updated}.c_#{len_c}.length.txt"

#------------------------------
# extract large contigs
#------------------------------
f_sh.puts "python ${SCRIPT_RENAME_} -min #{len_l} #{out_linear}.fa > #{out_linear}.#{len_l}.fa"
f_sh.puts "python ${SCRIPT_LENGTH_} #{out_linear}.#{len_l}.fa > #{out_linear}.#{len_l}.length.txt"

# updated contigs
f_sh.puts "cat  #{out_circular}.c_#{len_c}.fa  #{out_linear}.#{len_l}.fa > #{out_updated}.#{len_l}.c_#{len_c}.fa"
f_sh.puts "cat  #{out_circular}.c_#{len_c}.extended.fa  #{out_linear}.#{len_l}.fa > #{out_updated}.#{len_l}.c_#{len_c}.extended.fa"
f_sh.puts "cat  #{out_circular}.c_#{len_c}.length.txt  #{out_linear}.#{len_l}.length.txt | sort -k 1,1  > #{out_updated}.#{len_l}.c_#{len_c}.length.txt"
f_sh.puts "cat  #{out_circular}.c_#{len_c}.extended.length.txt  #{out_linear}.#{len_l}.length.txt | sort -k 1,1  > #{out_updated}.#{len_l}.c_#{len_c}.extended.length.txt"

f_sh.close
system("sh #{sh_}")

