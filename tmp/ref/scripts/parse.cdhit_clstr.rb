
require 'fileutils'
require 'optparse'
h_params = ARGV.getopts('i:', 'include_rep')
# input clstr file
in_clstr_ = h_params["i"]
include_rep = h_params["include_rep"]

f_in = open(in_clstr_)
rep_contig = ""
a_contigs = []
# loop CD-HIT clstr file
while line = f_in.gets
	line.chomp!
	if line[0] == ">" then
		if rep_contig != "" then
			# output
			a_contigs.each do |contig|
				puts [contig, rep_contig].join("\t")
			end
			rep_contig = ""
			a_contigs = []
		end
	elsif line[-1] == "*" then
		# representative contig
		rep_contig = line.split(">")[1][0..-6]
		# output representative contig
		if include_rep then
			puts [rep_contig, rep_contig].join("\t")
		end
	else
		# clustered contig
		temp = line.split(">")[1]
		contig = temp[0..(temp.index("... at") - 1)]
		a_contigs.push(contig)
	end
end
f_in.close()

# for last cluster
if rep_contig != "" then
	# output
	a_contigs.each do |contig|
		puts [contig, rep_contig].join("\t")
	end
end
