
start = nil
pre_id = ""
pre_pos = 0
while line = STDIN.gets
	line.chomp!
	a_all = line.split("\t")
	id = a_all[0]
	pos = a_all[1].to_i
	d = a_all[2].to_i

	if id != pre_id then
		if start != nil then
			puts [pre_id, start, pre_pos].join("\t")
		end
		if d > 0 then
			start = pos - 1
		else
			start = nil
		end
	else
		if start != nil then
			if d==0 then
				puts [pre_id, start, pre_pos].join("\t")
				start = nil
			end
		else
			if d > 0 then
				start = pre_pos
			end
		end
	end
	pre_id = id
	pre_pos = pos
end

