#!/bin/bash

out_dir=/cbio/grlab/nobackup/projects/sequencing_runs/F_graminearum/mitie/

fn_bam_all="/cbio/grlab/nobackup/projects/sequencing_runs/F_graminearum/alignments_tophat/F_graminearum_1.bam /cbio/grlab/nobackup/projects/sequencing_runs/F_graminearum/alignments_tophat/F_graminearum_2.bam"

mkdir -p $out_dir

dir=`dirname $0`
# create covered regions
fn_regions=${out_dir}/regions.flat
#fn_regions=~/tmp/sample_regin.flat
if [ ! -f $fn_regions ]
then
	echo run define_regions => $fn_regions
	#valgrind --tool=cachegrind $dir/define_regions $fn_bam_all -o $fn_regions
	#valgrind --leak-check=full $dir/define_regions $fn_bam_all -o $fn_regions
	$dir/define_regions $fn_bam_all -o $fn_regions --shrink --cut-regions --max-len 100000
fi

# create segment graph and store in file
fn_graph=${out_dir}/graph_gtf.bin

#fn_gtf=~/tmp/sample.gtf
if [ ! -f $fn_graph ]
then
	echo
	echo generate graph on bam file $fn_bam_all
	echo
	$dir/generate_segment_graph $fn_graph.tmp $opt --regions $fn_regions --max-junk 200000 --seg-filter 0.05 --region-filter 50 --tss-tts-pval 0.0001 $fn_bam_all

	if [ ! "$?" -eq "0" ]; then
		echo generate_segment_graph return: $?
		exit 0
	fi
	mv ${fn_graph}.tmp $fn_graph
fi

num_trans=0 # number of additional transcripts 
lambda=3
eta1=1.4
eta2=0.4
./transcript_prediction $fn_graph $fn_bam_all --max-num-trans $num_trans --param-eta1 $eta1 --param-eta2 $eta2 --param-lambda $lambda 

