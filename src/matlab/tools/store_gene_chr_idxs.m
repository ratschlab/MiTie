function all_idxs = store_gene_chr_idxs(genes, chr_nums);

	max_chr_num = max(chr_nums);
	all_idxs = cell(max_chr_num);
	all_num = zeros(1, max_chr_num);
	
	% allocate some mem
	for j = 1:length(all_idxs)
		all_idxs{j} = zeros(1, 1000);
	end

	for j = 1:length(genes)
		chr_num = genes(j).chr_num;
		if all_num(chr_num)+1>length(all_idxs{chr_num})
			% double allocatet mem
			all_idxs{chr_num} = [all_idxs{chr_num}, zeros(1, length(all_idxs{chr_num}))];
		end
		all_idxs{chr_num}(all_num(chr_num)+1) = j;
		all_num(chr_num) = all_num(chr_num)+1; 
	end
	all_genes_ids = [];
	for j = 1:length(all_idxs)
		all_idxs{j}(all_num(j)+1:end) = [];
		all_genes_ids = [all_genes_ids, all_idxs{j}];
		assert(all(all_idxs{j}>0))
	end
	assert(length(all_genes_ids)==length(genes))
	assert(length(all_genes_ids)==length(unique(all_genes_ids)))
return

