
dir_ = File.dirname($0)

puts "--- 01 make contig coverage profile ---"
system("ruby #{dir_}/16_contig_profile_table.01_coverage.for_bacteriome.rb")
puts ""
puts "--- 02 combine profile with taxonomic annotation ---"
system("ruby #{dir_}/16_contig_profile_table.02_annotation.for_bacteriome.rb")
