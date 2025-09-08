
require 'fileutils'

# Inputs
#in_blast_ = ARGV[0] # blast result file
#contig_length_ = ARGV[1] # contig length file
#cov_thre = ARGV[2].to_i # threshold for coverage of query sequence
#d_merge = ARGV[3].to_i # minimum distance for merging regions
#out_info_ = ARGV[4] # output file of summary information
#out_rd_ = ARGV[5] # output file of redundant contigs

require 'optparse'
h_params = ARGV.getopts('', 'blast_file:', 'length_file:', 'cov:95', 'd_merge:6', 'info_file:', 'out_file:')
in_blast_ = h_params["blast_file"] # blast result file
contig_length_ = h_params["length_file"] # contig length file
cov_thre = h_params["cov"].to_i # threshold for coverage of query sequence
d_merge = h_params["d_merge"].to_i # minimum distance for merging regions
out_info_ = h_params["info_file"] # output file of summary information
out_rd_ = h_params["out_file"] # output file of redundant contigs

	# load contig length
	h_contig_length = {}
	f_length = open(contig_length_)
	while line = f_length.gets
		a_all = line.split("\t")
		contig = a_all[0]
		length = a_all[1].to_i
		h_contig_length[contig] = length
	end
	f_length.close()

	# summary information of blast hits
	h_pair_info = {}
	f_in = open("| sort -t $'\\t' -k1,1 -k2,2 -k7,7n #{in_blast_}")
	pre_pair = ""
	pre_q_end = ""
	pre_s_end = ""
	strand = ""
	al_length = ""
	flag = ""
	query_region = ""
	subject_region = ""

	#------------------------------ 
	# load blast hits and judge if an aligned contig is contiguous
	#------------------------------ 
	while line = f_in.gets
		line.chomp!
		a_all = line.split("\t")
		query = a_all[0]
		subject = a_all[1]
		al_len = a_all[3].to_i
		q_start = a_all[6].to_i
		q_end = a_all[7].to_i
		s_start = a_all[8].to_i
		s_end = a_all[9].to_i
		pair = "#{query}:#{subject}"

		if pair != pre_pair then
		# 1st hit record
			if s_start < s_end then
				strand = "+"
			else
				strand = "-"
			end
			flag = ""
			q_start_whole = q_start
			q_end_whole = q_end
			al_length = al_len
			query_region = "#{q_start}-#{q_end}"
			subject_region = "#{s_start}-#{s_end}"

			pre_pair = pair
			pre_q_end = q_end
			pre_s_end = s_end
			h_pair_info[pair] = [q_start_whole, q_end_whole, strand, al_length, flag, query_region, subject_region]
		else
		# 2nd or later hit record
			a_pair_info = h_pair_info[pair]
			strand = a_pair_info[2]
			flag = a_pair_info[4]
			ignore = false

			# judge if an aligned contig is contiguous
			s_contig_length = h_contig_length[subject]
			if q_start <= pre_q_end + d_merge and q_start >= pre_q_end - d_merge then
				if strand == "+" then
					if s_start > s_end then
						flag += "subject strand is not consistent. "
					elsif s_start <= pre_s_end + d_merge and s_start >= pre_s_end - d_merge then
					elsif s_start <= d_merge and pre_s_end + d_merge >= s_contig_length then
						# border of circular contig
					else
						flag += "subject is not contiguous. "
					end
				else # strand == "-"
					if s_start < s_end then
						flag += "subject strand is not consistent. "
					elsif s_start <= pre_s_end + d_merge and s_start >= pre_s_end - d_merge then
					elsif s_start + d_merge >= s_contig_length and pre_s_end <= d_merge then
						# border of circular contig
					else
						flag += "subject is not contiguous. "
					end
				end
			else
				if q_end <= pre_q_end + d_merge then
					# sub sequence of the previous query
					# --> ignore
					ignore = true
				else
					flag += "query is not contiguous. "
				end
			end

			if !ignore then
				q_start_whole = a_pair_info[0]
				q_end_whole = q_end
				al_length = a_pair_info[3] + al_len
				query_region = h_pair_info[pair][5] + ",#{q_start}-#{q_end}"
				subject_region = h_pair_info[pair][6] + ",#{s_start}-#{s_end}"
				pre_pair = pair
				pre_q_end = q_end
				pre_s_end = s_end
				h_pair_info[pair] = [q_start_whole, q_end_whole, strand, al_length, flag, query_region, subject_region]
			end

		end

	end
	f_in.close()

	# output total hit length
	f_out_info = open("#{out_info_}", "w")
	h_common = {}
	h_subj2common = {}
	h_subj2query = {}
	h_query2common = {}
	id_common = 0
	h_pair_info.keys.each do |pair|
		query = pair.split(":")[0]
		subject = pair.split(":")[1]
		a_info = h_pair_info[pair]
		al_length = a_info[3]
		flag = a_info[4]
		# length of contigs
		q_contig_length = h_contig_length[query]
		s_contig_length = h_contig_length[subject]
		# coverage of query
		#coverage_query = (a_info[1] - a_info[0]).to_f / q_contig_length.to_f
		coverage_query = al_length.to_f / q_contig_length.to_f
		coverage_subject = al_length.to_f / s_contig_length.to_f

		# ouput of redandant contigs
		
		if coverage_query * 100 >= cov_thre.to_f and flag == "" then
			f_out_info.puts [pair, "redundant", query, q_contig_length, subject, s_contig_length, coverage_query, coverage_subject, a_info].flatten.join("\t")

			if h_subj2common.key?(subject) then
				# subject is already registered as common
				common = h_subj2common[subject]
				# 1. no need to register subject as common
				# 2. add query to its query list
				a_query = h_subj2query[subject]
				a_query.push(query)
				# 3. register query info
				h_query2common[query] = common

			elsif h_subj2common.key?(query) then
				# query is already registered as common
				common = h_subj2common[query]
				# 1. add subject as common
				h_subj_info = h_common[common]
				h_subj_info[subject] = s_contig_length
				h_subj2common[subject] = common
				# 2. make new query list
				h_subj2query[subject] = [query]
				# 3. register query info
				h_query2common[query] = common

			elsif h_query2common.key?(subject) then
				# subject is linked to common as a query
				common = h_query2common[subject]
				# 1. add subject as common
				h_subj_info = h_common[common]
				h_subj_info[subject] = s_contig_length
				h_subj2common[subject] = common
				# 2. make new query list
				h_subj2query[subject] = [query]
				# 3. register query info
				h_query2common[query] = common

			elsif h_query2common.key?(query) then
				# query is linked to common as a query
				common = h_query2common[query]
				# 1. add subject as common
				h_subj_info = h_common[common]
				h_subj_info[subject] = s_contig_length
				h_subj2common[subject] = common
				# 2. make new query list
				h_subj2query[subject] = [query]
				# 3. no need to register query info

			else
				# new common contig
				id_common += 1
				# 1. register new subject as common
				h_common[id_common] = {subject => s_contig_length}
				h_subj2common[subject] = id_common
				# 2. make new query list
				h_subj2query[subject] = [query]
				# 3. register query info
				h_query2common[query] = id_common

			end

		else
			f_out_info.puts [pair, "", query, q_contig_length, subject, s_contig_length, coverage_query, coverage_subject, a_info].flatten.join("\t")
		end

	end
	f_out_info.close()

	# output information of common contigs
	f_out_com = open("#{out_rd_}", "w")
	h_common.each_pair do |common, h_subj_info|
		# selection of representative conitg
		a_subject = []
		a_length = []
		a_query_count = []
		# sort by contig length in decreasing order. if contigh lengths are equal, sort by contig name in decreasing lexicographical order
		h_subj_info.to_a.sort{|a, b| (b[1] <=> a[1]) * 2 + (b[0] <=> a[0])}.each do |subject, s_contig_length|
			a_subject.push(subject)
			a_length.push(s_contig_length)
		end
		# representative conitg
		rep_subject = a_subject[0]
		# redundant contigs
		a_query_redundant = h_subj2query[rep_subject]
		# output
		a_query_redundant.each do |query_redundant|
			q_contig_length = h_contig_length[query_redundant]
			f_out_com.puts [query_redundant, q_contig_length, common, a_subject.join(","), a_length.join(","), rep_subject].flatten.join("\t")
		end
	end
	f_out_com.close()

