function genes = load_graph_bin(fn_graph, num, fn_bam)

if nargin>2
	if iscell(fn_bam)
		genes = load_regions_bin(fn_graph, num, fn_bam{:});
	else
		genes = load_regions_bin(fn_graph, num, fn_bam);
	end
else
	genes = load_regions_bin(fn_graph, num);
end
if length(genes)==0 %|| gene_cnt>100  
	% close input stream
	close_flag = -1;
	load_regions_bin(fn_graph, close_flag);
end
for j = 1:length(genes)
	genes(j).id = j;
	genes(j).initial   = double(genes(j).seg_admat(1, 2:end-1)>-2);
	genes(j).terminal  = double(genes(j).seg_admat(end, 2:end-1)>-2);
	genes(j).seg_admat = genes(j).seg_admat(2:end-1, 2:end-1);
	for s = 1:size(genes(j).segments, 2)-1
		if genes(j).segments(2, s)+1==genes(j).segments(1, s+1)
			genes(j).seg_admat(s, s+1) = -1;
			genes(j).seg_admat(s+1, s) = -1;
		end
	end


	% create pair list
	if ~isempty(genes(j).pair_mat)
		[a b] = find(genes(j).pair_mat>0);
		cnt = 0;
		all_pl = zeros(length(a), 3);
		for k = 1:length(a)
			if a(k)<b(k)
				cnt = cnt+1;
				all_pl(cnt, :) = [a(k), b(k), genes(j).pair_mat(a(k), b(k))];
			end
		end
		all_pl(cnt+1:end, :) = [];
		genes(j).pair_list = all_pl;
		%if sample==1
		%	pair_list = all_pl;
		%else
		%	for k = 1:size(all_pl)
		%		[tmp idx] = ismember(all_pl(k, 1:2), pair_list(:, 1:2), 'rows');
		%		if tmp>0
		%			pair_list(idx, 3) = pair_list(idx, 3) + all_pl(k, 3);
		%		end
		%	end
		%	[tmp, idx] = ismember(all_pl(:, 1:2), pair_list(:, 1:2), 'rows');
		%	pair_list = [pair_list; all_pl(idx==0, :)];
		%end
	end

end

