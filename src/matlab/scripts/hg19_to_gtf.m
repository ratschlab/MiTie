

for j = 1:5

	fn_mat = sprintf('/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr%i.mat', j)
	fn_gtf = sprintf('/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr%i.gtf', j)

	clear genes
	load(fn_mat, 'genes')
	write_gtf(genes, fn_gtf, 'annotation')
end
