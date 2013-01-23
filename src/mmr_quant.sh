#!/bin/bash

MAT="/fml/ag-raetsch/share/software/matlab-7.6/bin/matlab -nojvm -nodesktop -nosplash"
addpaths="addpath matlab; "


base_dir=/fml/ag-raetsch/nobackup2/projects/mmr/human_simulation/5000_genes_20000000_reads
mitie_dir=$base_dir/MiTie
fn_gtf=/fml/ag-raetsch/nobackup2/projects/mmr/human_simulation/annotation/hg19_subsample_5000_genes.nochr.gtf

for mmr in "" mmr0. mmr1.; do
	for noise in "" noise0.01. noise0.02. noise0.03.; do
		fn_bam=$base_dir/hg19_subsample_5000_genes.gtf.${noise}fastq.gz.mapped.2.${mmr}sorted.bam
		#ls -l $fn_bam

		if [ -z $mmr ] && [ -z $noise ]; then 
			out_dir=$mitie_dir/unfiltered
		elif [ -z $mmr ]; then 
			out_dir=$mitie_dir/unfiltered.${noise%.}
		elif [ -z $noise ]; then 
			out_dir=$mitie_dir/${mmr%.}
		else
			out_dir=$mitie_dir/${mmr}${noise%.}
		fi
		echo $out_dir
		mkdir -p $out_dir

		fn_graph=$out_dir/graph.bin
		opts="--gtf-offset 0"
		#./generate_segment_graph ${fn_graph}.tmp $opts --gtf $fn_gtf $fn_bam 
		mv ${fn_graph}.tmp ${fn_graph}

		mip_dir=$out_dir/pred
		mkdir -p $mip_dir
		quantify=1
		idx="1:100000" # predict on all genes
		${MAT} -r "dbstop error; $addpaths mip_paths; transcript_predictions('$fn_graph', '$fn_bam', '$mip_dir', [], $idx , $quantify); exit"
		exit 1;
	done
done
