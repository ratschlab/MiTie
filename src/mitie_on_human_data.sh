#!/bin/bash
#

fn_bam_brain=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.mmr.sorted.bam
fn_bam_uhr=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.mmr.sorted.bam
fn_bam_brain_all=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.sorted.bam
fn_bam_uhr_all=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.sorted.bam

out_dir=/fml/ag-raetsch/nobackup/projects/mip/chris_mason/MiTie
fn_regions=$out_dir/regions.flat

dir=`dirname $0`

if [ ! -f $fn_regions ]
then
	echo run define_regions => $fn_regions
	#valgrind --tool=cachegrind $dir/define_regions $fn_bam_all -o $fn_regions
	#valgrind --leak-check=full $dir/define_regions $fn_bam_all -o $fn_regions
	$dir/define_regions $fn_bam_brain $fn_bam_uhr -o $fn_regions --shrink --cut-regions --max-len 100000
fi

# create segment graph and store in file
fn_graph=${out_dir}/graph

if [ ! -f $fn_graph ]
then
	echo
	echo generate graph on bam file $fn_bam_all
	echo
	#opts="--mismatches 3 --min-exonic-len 5"
	echo $dir/generate_segment_graph $fn_graph.tmp $opts --split-chr --regions $fn_regions --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001 $fn_bam_brain_all $fn_bam_uhr_all
	$dir/generate_segment_graph $fn_graph.tmp $opts --split-chr --regions $fn_regions --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001 $fn_bam_brain_all $fn_bam_uhr_all
	
	if [ ! "$?" -eq "0" ]; then
		echo generate_segment_graph return: $?
		exit 0
	fi
	mv ${fn_graph}.tmp $fn_graph
fi
exit 0


##############################	
##############################	
mip_dir=$out_dir/MiTie_pred/
mkdir -p $mip_dir
MAT="/fml/ag-raetsch/share/software/matlab-7.6/bin/matlab -nojvm -nodesktop -nosplash"
addpaths="addpath matlab; "
${MAT} -r "dbstop error; $addpaths mip_paths; denovo('$fn_graph', {'$fn_bam_brain', '$fn_bam_uhr'}, '$mip_dir'); exit"


