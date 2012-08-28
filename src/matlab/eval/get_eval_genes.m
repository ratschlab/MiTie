function genes = get_eval_genes(mode)
assert(isequal(mode, 'train') || isequal(mode, 'test'))

fn_genes = sprintf('/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/eval_%s.mat', mode);

if ~fexist(fn_genes)
	fn_genes2 = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr1.mat';
	
	load(fn_genes2, 'genes')
	
	% training set
	if isequal(mode, 'train')
		genes = genes(1:500);
	else
		genes = genes(501:1000);
	end

	fn_save='/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr_all_expr_orig.mat';
	xx = load(fn_save);
	ids = [xx.genes.id];
	% remove transcripts with zero expression in all samples
	for j = 1:length(genes)
		gidx = find(ids==genes(j).id);
		assert(length(gidx)==1)
		expr = sum(xx.genes(gidx).expr_orig(2:end, :)); % row one is the true hidden expression value
		rank = -ones(1, length(genes(j).exons));
		for k = 1:length(genes(j).exons)
			rank(k) = sum(genes(j).expr_orig>genes(j).expr_orig(k))+1;
		end
		idx = find(expr>0);
		assert(any(rank(idx)==1))
		genes(j).exons 		= genes(j).exons(idx);
		genes(j).transcripts= genes(j).transcripts(idx);
		genes(j).expr_orig 	= expr(idx);
		genes(j).rank = rank(idx);
	end
	save(fn_genes, 'genes')
else
	load(fn_genes, 'genes')
end
