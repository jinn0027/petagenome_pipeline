
require 'fileutils'

in_1_ = ARGV[0]
col_1 = ARGV[1]
in_2_ = ARGV[2]
col_2 = ARGV[3]
option = ARGV[4]

dir_ = Dir.getwd
# sort file
in_1_sort_ = "#{dir_}/1_sort_#{col_1}.#{$$}"
system("sort -t$'\\t' -k#{col_1},#{col_1} #{in_1_} > #{in_1_sort_}")
# sort file
in_2_sort_ = "#{dir_}/2_sort_#{col_2}.#{$$}"
system("sort -t$'\\t' -k#{col_2},#{col_2} #{in_2_} > #{in_2_sort_}")

system("join #{option} -1 #{col_1} -2 #{col_2} -t \"$(printf '\011')\"  #{in_1_sort_}  #{in_2_sort_}")

# unlink sort file
File.unlink("#{in_1_sort_}")
File.unlink("#{in_2_sort_}")

