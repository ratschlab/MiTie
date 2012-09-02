function exons = segments2exons(segments, transcript)


exons = segments(transcript(1), :); 

for j = 2:length(transcript)
	seg = segments(transcript(j), :);
	if seg(1)==exons(end, 2)+1
		exons(end, 2) = seg(2);
	else
		exons(end+1, :) = seg;
	end
end
ilen = exons(2:end, 1) - exons(1:end-1, 2);
assert(all(ilen>1))
