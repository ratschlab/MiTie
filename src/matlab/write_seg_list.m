function write_seg_list(res_dir, out_file)

fd = fopen(out_file, 'w');
if fd<0
	error('could not open file: %i for writing\n', out_file)
end
l = load_pred(res_dir);

for j = 1:length(l)

	if l(j).solvetime<0
		continue
	end

	seg = l(j).segments;
	for s = 1:size(seg, 1)
		weight = 0;
		for t = 1:length(l(j).transcripts)
			if ismember(s, l(j).transcripts{t})
				weight = weight+l(j).weights(end, t);
			end
		end
		fprintf(fd, 's\t%s\t%s\t%i\t%i\t%.2f\n', l(j).gene.chr, l(j).gene.strand, seg(s, 1), seg(s, 2), weight*l(j).cov_scale(end));
	end
	sa = l(j).seg_admat(:,:,end);
	[s1, s2] = find(sa>=0);
	for i=1:length(s1)
		if s1(i)>=s2(i)
			continue;
		end
		weight = 0;
		for t = 1:length(l(j).transcripts)
			if ismember(s1(i), l(j).transcripts{t}) && ismember(s2(i), l(j).transcripts{t}) && ~any(ismember(s1(i)+1:s2(i)-1, l(j).transcripts{t}))
				
				weight = weight+l(j).weights(end, t);
			end
		end
		fprintf(fd, 'i\t%s\t%s\t%i\t%i\t%.2f\n', l(j).gene.chr, l(j).gene.strand, seg(s1(i), 2), seg(s2(i), 1), weight*l(j).cov_scale(end));
		
	end
end
fclose(fd);
