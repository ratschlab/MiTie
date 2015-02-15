#!/bin/bash
MAT="/fml/ag-raetsch/share/software/matlab-7.6/bin/matlab -nojvm -nodesktop -nosplash"
addpaths="addpath matlab; mip_paths; dbstop error; "

for chr in Chr1 Chr2 Chr3 Chr4 Chr5; do
	for s in + -; do
		echo start with $chr$s
		${MAT} -r "$addpaths transcript_predictions_MAGIC('$chr', '$s'); exit"
	done
done

#sleep 1000
base=/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/MiTie
for chr in Chr1 Chr2 Chr3 Chr4 Chr5; do
	for s in + -; do
		echo start with $chr$s
		dir_=$base/pred_$chr$s/
		fn_save=$dir_/res_genes.mat
		add_weights=1
		mmr=1
		write_gtf_flag=1
		echo ${MAT} -r "$addpaths collect_results('$dir_', '$fn_save', $add_weights, $mmr, $write_gtf_flag); exit"
		${MAT} -r "$addpaths collect_results('$dir_', '$fn_save', $add_weights, $mmr, $write_gtf_flag); exit"
	done
done

${MAT} -r "$addpaths collect_magic_genes(); exit"

# debug code
#fn_graph = '/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/MiTie/graphs/seg_graph30_filtered_Chr1+'
#fn_bam = '/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/alignments_magic/bam_mmr/all_merged.WTCHG_34033_162.sorted.mmr.bam'
#genes = load_regions_bin(fn_graph, 10000, fn_bam)
