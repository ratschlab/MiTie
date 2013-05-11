#!/bin/bash
#
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_gtf_sample${sample} /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_sample${sample} /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_fixes_sample${sample} /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_filter_sample${sample} /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_no_shrink_sample${sample} /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
#
#

if [ -z $1 ]; then
	sample=2
else
	sample=$1
fi
eta1=1.00
eta2=0.42
lambda=3

#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_release_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_tss0.01_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_enum5_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_enum5_1e4_repeat_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_seg_filter0.01_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions40000_no_junc_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_GP_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_GP_sample${sample}_filter
out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_quant_sample${sample}_eta1_${eta1}_eta2_${eta1}_lambda_${lambda}
#out_dir=~/tmp
#fn_bam_all=/cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam

#for s in `seq 1 $sample`; do
#	fn_bam_all="$fn_bam_all /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${s}_merged_err_1.no_junc.sorted.paired.bam"
#	#fn_bam_all="$fn_bam_all /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam"
#done
fn_bam_all=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired_200000_4_5.bam 


mkdir -p $out_dir

dir=`dirname $0`

# create segment graph and store in file
fn_graph=${out_dir}/graph_gtf.bin

#fn_gtf=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr_max_trans.gtf
fn_gtf=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr1.gtf
#fn_gtf=~/tmp/sample.gtf
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
echo "dbstop error; $addpaths; mip_paths; C.num_transcripts = 0; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir', C, [], $quantify, $eta1, $eta2, $lambda); exit"
${MAT} -r "dbstop error; $addpaths; mip_paths; C.num_transcripts = 0; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir', C, [], $quantify, $eta1, $eta2, $lambda); exit"

##############################	
# eval mip
##############################		
#${MAT} -r "dbstop error; $addpaths mip_paths; eval_dir('$mip_dir', 'train'); exit"

