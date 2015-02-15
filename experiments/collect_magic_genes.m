function genes = collect_magic_genes()

base='/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/MiTie/pred_'

genes = [];
for chr = {'Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'}
	for strand ='+-'
		fname=sprintf('%s%s%s/res_genes.mat', base, chr{1}, strand); 
		l = load(fname, 'genes');
		if isempty(genes)
			genes = l.genes;
		else
			genes = [genes l.genes];
		end
	end
end

num = zeros(1,length(genes));
for j = 1:length(genes); for k = 1:length(genes(j).transcripts), if strfind(genes(j).transcripts{k}, 'new_trans'), num(j) = num(j)+1; end, end, end
fprintf('%i genes have new transcripts predicted (%.2f%%)\n', sum(num), mean(num)*100)
