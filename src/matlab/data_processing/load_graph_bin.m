function genes = load_graph_bin(fn_graph, fn_bam, gene_idx)

genes = [];
for gene_id=gene_idx
	if nargin>2 && ~isempty(fn_bam)
		if iscell(fn_bam)
			gg = load_regions_bin(fn_graph, gene_id, fn_bam{:});
		else
			gg = load_regions_bin(fn_graph, gene_id, fn_bam);
		end
	else
		gg = load_regions_bin(fn_graph, gene_id);
	end
	if isempty(genes)
		genes = gg;
	else
		genes = [genes, gg];
	end
end
for j = 1:length(genes)
	genes(j).id = j;
	% make admat symmetric
	for k = 1:size(genes(j).seg_admat, 1)
		for l = k+1:size(genes(j).seg_admat, 2) 
			genes(j).seg_admat(k,l) = max(genes(j).seg_admat(k,l), genes(j).seg_admat(l, k));
			genes(j).seg_admat(l,k) = max(genes(j).seg_admat(k,l), genes(j).seg_admat(l, k));
		end
	end
	genes(j).initial   = double(genes(j).seg_admat(1, 2:end-1, 1)>-2);
	genes(j).terminal  = double(genes(j).seg_admat(end, 2:end-1, 1)>-2);
	genes(j).seg_admat = genes(j).seg_admat(2:end-1, 2:end-1, :);
	for s = 1:size(genes(j).segments, 2)-1
		if genes(j).segments(2, s)+1==genes(j).segments(1, s+1)
			genes(j).seg_admat(s, s+1, :) = -1;
			genes(j).seg_admat(s+1, s, :) = -1;
		end
	end

	% loop over samples and make sure that all connections 
	% can be found in each sample
	all_connections = sum(genes(j).seg_admat>-1, 3)>0; % all connections that have RNA-seq evidence in any sample
	for s = 1:size(genes(j).seg_admat, 3)
		genes(j).seg_admat(:,:,s) = max(genes(j).seg_admat(:,:,s), all_connections*2-2); 
	end

	% create pair list
	if ~isempty(genes(j).pair_mat)
		[a b] = find(sum(genes(j).pair_mat, 3)>0);
		cnt = 0;
		all_pl = zeros(length(a), 3);
		for k = 1:length(a)
			if a(k)<b(k)
				cnt = cnt+1;
				all_pl(cnt, 1:2) = [a(k), b(k)];
				for s = 1:size(genes(j).pair_mat, 3)
					all_pl(cnt, 2+s) =  genes(j).pair_mat(a(k), b(k), s);
				end
			end
		end
		all_pl(cnt+1:end, :) = [];
		genes(j).pair_list = all_pl;
	end

end

