#!/bin/bash
#


if [ -z $1 ]; then
	echo usage $0 '<max num mismatches>'
	exit 0
else
	num_missmatches=$1;
fi
eta1=1.00
eta2=0.1
lambda=0

out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_quant_exon_eta1_${eta1}_eta2_${eta2}_lambda_${lambda}_mm$num_missmatches
#out_dir=~/tmp

fn_bam_all=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired_200000_4_$num_missmatches.bam


mkdir -p $out_dir

dir=`dirname $0`

# create segment graph and store in file
fn_graph=${out_dir}/graph_gtf.h5

#fn_gtf=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr_max_trans.gtf
fn_gtf=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr1.gtf
if [ ! -f $fn_graph ]
then
	echo
	echo generate graph on bam file $fn_bam_all
	echo
	$dir/generate_segment_graph $fn_graph.tmp --gtf $fn_gtf

	if [ ! "$?" -eq "0" ]; then
		echo generate_segment_graph return: $?
		exit 0
	fi
	mv ${fn_graph}.tmp $fn_graph
fi

#valgrind --leak-check=full $dir/generate_segment_graph $fn_regions $fn_graph --seg-filter 0.05 --region-filter 0 --gtf $fn_gtf $fn_bam_all


##############################	
# mip 
# generate expected quantification values for each segment
# and expected intron counts for each intron;
# treat bam files as separate samples
##############################	
fn_bam_iter=$fn_bam_all
mip_dir=$out_dir
MAT="/cbio/grlab/share/software/matlab/matlab_R2012b/bin/matlab -nojvm -nodesktop -nosplash"
addpaths="addpath matlab; "
quantify=1

${MAT} -r "dbstop error; $addpaths; mip_paths; C.num_transcripts = 0; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir', C, 1:10 , $quantify, $eta1, $eta2, $lambda); exit"
exit 0;

for x in `seq 1 1010`; do
	fn_res=$mip_dir/gene$x.mat
	if [ ! -f $fn_res ]; then
		${MAT} -r "dbstop error; $addpaths; mip_paths; C.num_transcripts = 0; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir', C, $x , $quantify, $eta1, $eta2, $lambda); exit"
		sleep 3
	fi
done
#${MAT} -r "dbstop error; $addpaths; mip_paths; C.num_transcripts = 0; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir', C, [], $quantify, $eta1, $eta2, $lambda); exit"
#echo "debug_on_error(1); $addpaths; mip_paths; C.num_transcripts = 0; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir', C, [], $quantify, $eta1, $eta2, $lambda); exit" | octave

##############################	
# eval mip
##############################		
#${MAT} -r "dbstop error; $addpaths mip_paths; eval_dir('$mip_dir', 'train'); exit"

