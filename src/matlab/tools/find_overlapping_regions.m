function idx = find_overlapping_regions(regions1, regions2, consider_strands)
%%

if isempty(regions1) & isempty(regions2)
  idx = [];
  return;
end
  
if isempty(regions2)&&(nargin==2 || consider_strands)
	%% find overlapping pairs within one set of regions
	chr_nums = unique([regions1.chr_num]);
	all_gene_chr_idxs1 = store_gene_chr_idxs(regions1, chr_nums);%fast if there are many contigs
	idx = [];
	for chr=chr_nums
		for s='+-'
			%fprintf('\rchr%i strand%s',chr,s)
			chr_genes_idx1 = all_gene_chr_idxs1{chr}; 
			str_genes_idx1 = find([regions1(chr_genes_idx1).strand]==s);
			r1_chr_idx = chr_genes_idx1(str_genes_idx1);

			if isempty(r1_chr_idx)
				continue
			end 
			idx = compute_overlap(regions1, r1_chr_idx, regions1, r1_chr_idx, idx);
		end
	end
	idx(idx(:,1)==idx(:,2),:) = [];
	idx = sort(idx')';
	idx = unique(idx, 'rows');
	return
end

if ~isfield(regions1, 'chr_num') || ~isfield(regions2, 'chr_num') || length([[regions1.chr_num] [regions2.chr_num]]) ~= length(regions1)+length(regions2)
	chr = unique([{regions1.chr} {regions2.chr}]); 	
	regions1 = set_chr_num(regions1, chr);	
	regions2 = set_chr_num(regions2, chr);	
end

chr_nums = unique([[regions1.chr_num] [regions2.chr_num]]);
all_gene_chr_idxs1 = store_gene_chr_idxs(regions1, chr_nums);%fast if there are many contigs
all_gene_chr_idxs2 = store_gene_chr_idxs(regions2, chr_nums);

if nargin==2 || consider_strands
	idx = [];
	for chr=chr_nums
		for s='+-'
			%fprintf('\rchr%i strand%s',chr,s)
			%r1_chr_idx = find([regions1.chr_num]==chr&[regions1.strand]==s);
			%r2_chr_idx = find([regions2.chr_num]==chr&[regions2.strand]==s);
			chr_genes_idx1 = all_gene_chr_idxs1{chr}; % precomputed to make it faster for many contigs
			str_genes_idx1 = find([regions1(chr_genes_idx1).strand]==s);
			r1_chr_idx = chr_genes_idx1(str_genes_idx1);

			chr_genes_idx2 = all_gene_chr_idxs2{chr}; 
			str_genes_idx2 = find([regions2(chr_genes_idx2).strand]==s);
			r2_chr_idx = chr_genes_idx2(str_genes_idx2);

			if isempty(r1_chr_idx)||isempty(r2_chr_idx)
				continue
			end 
			idx = compute_overlap(regions1, r1_chr_idx, regions2, r2_chr_idx, idx);
		end
	end
else
	idx = [];
	for chr=chr_nums
		%fprintf('\rchr%i',chr)
		r1_chr_idx = all_gene_chr_idxs1{chr};
		r2_chr_idx = all_gene_chr_idxs2{chr};

		if isempty(r1_chr_idx)||isempty(r2_chr_idx)
			continue
		end 
		idx = compute_overlap(regions1, r1_chr_idx, regions2, r2_chr_idx, idx);
	end
end
return

function idx = compute_overlap(regions1, r1_idx, regions2, r2_idx, idx)
		starts1 = [regions1(r1_idx).start];
		stops1  = [regions1(r1_idx).stop];
		starts2 = [regions2(r2_idx).start];
		stops2  = [regions2(r2_idx).stop];
		[a b] = interval_overlap(starts1, stops1, starts2, stops2);
		idx = [idx; [r1_idx(a)' r2_idx(b)']];
return
