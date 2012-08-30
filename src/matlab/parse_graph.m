function genes = parse_graph(fn_graph)

fd = fopen(fn_graph, 'r');
max_num = 1e6;
cnt = 0;
while ~feof(fd)
	line = fgetl(fd);
	if (mod(cnt, 100)==0)
		fprintf('\r parsed %i graphs    ', cnt);
	end
	if (strfind(line, 'region')==1)
		items = separate(line);
		chr = items{2};
		strand = items{3};
		start = str2num(items{4});
		stop = str2num(items{5});

		% parse segment header line
		line = fgetl(fd);
		assert(strfind(line, 'segments')==1)
		items = separate(line);
		num_seg = str2num(items{2});

		% parse segments
		segments = zeros(num_seg, 2);
		for j=1:num_seg
			line = fgetl(fd);
			[a b] = sscanf(line, '%i\t%i\t%*f');
			assert(length(a)==2);
			start = a(1);
			stop = a(2);
			segments(j, :) = [start, stop];
		end

		admat = zeros(num_seg+2, num_seg+2);
		while 1
			line = fgetl(fd);
			if (strfind(line, 'end')==1)
				break;
			end
			[a b] = sscanf(line, '%i\t%i');
			assert(length(a)==2);
			seg1 = a(1);
			seg2 = a(2);
			admat(seg1, seg2) = 1;
		end
		genes(max_num-cnt).chr = chr;
		genes(max_num-cnt).strand = strand;
		genes(max_num-cnt).start = start;
		genes(max_num-cnt).stop = stop;
		genes(max_num-cnt).segments = segments;
		% first and last node are source and sink
		genes(max_num-cnt).seg_admat = admat(2:end-1, 2:end-1);
		genes(max_num-cnt).initial = admat(1, 2:end-1);
		genes(max_num-cnt).terminal = admat(end, 2:end-1);
		cnt = cnt+1;
	end
end
genes(1:max_num-cnt) = [];
genes = genes(end:-1:1);
for j = 1:length(genes)
	genes(j).id = j;
end
fclose(fd);


return
