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
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_release_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_tss0.01_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_enum5_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_enum5_1e4_repeat_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_seg_filter0.01_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions40000_no_junc_sample${sample}
#out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_GP_sample${sample}
out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/mip_mmr_GP_sample${sample}_filter
#out_dir=~/tmp
#fn_bam_all=/cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam

#for s in `seq 1 $sample`; do
#	fn_bam_all="$fn_bam_all /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${s}_merged_err_1.no_junc.sorted.paired.bam"
#	#fn_bam_all="$fn_bam_all /cbio/grlab/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam"
#done
fn_bam_all=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired_200000_4_5.bam 


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

fn_gtf=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr_max_trans.gtf
#fn_gtf=~/tmp/sample.gtf
if [ ! -f $fn_graph ]
then
	echo
	echo generate graph on bam file $fn_bam_all
	echo
	#opts="--mismatches 3 --min-exonic-len 5"
	$dir/generate_segment_graph $fn_graph.tmp $opt --regions $fn_regions --max-junk 200000 --gtf-offset 50000 --seg-filter 0.05 --region-filter 50 --tss-tts-pval 0.0001 --gtf $fn_gtf $fn_bam_all
	#valgrind --leak-check=full --show-reachable=yes $dir/generate_segment_graph $fn_graph.tmp $opt --regions $fn_regions --max-junk 200000 --gtf-offset 50000 --seg-filter 0.05 --region-filter 50 --tss-tts-pval 0.0001 $fn_bam_all
	#exit -1
	#$dir/generate_segment_graph $fn_graph.tmp $opt --regions $fn_regions --gtf-offset 50000 --seg-filter 0.05 --region-filter 50 --tss-tts-pval 0.0001 $fn_bam_all
	#$dir/generate_segment_graph $fn_graph.tmp $opt --regions $fn_regions --gtf-offset 50000 --seg-filter 0.05 --region-filter 50 --tss-tts-pval 0.0001 --gtf $fn_gtf $fn_bam_all
	#gdb $dir/generate_segment_graph
	#$dir/generate_segment_graph ${fn_graph}.tmp $opts --few-regions --gtf-offset 40000 --seg-filter 0.01 --region-filter 50 --tss-tts-pval 0.0001 --gtf $fn_gtf $fn_bam_all 
	

	if [ ! "$?" -eq "0" ]; then
		echo generate_segment_graph return: $?
		exit 0
	fi
	mv ${fn_graph}.tmp $fn_graph
fi

#valgrind --leak-check=full $dir/generate_segment_graph $fn_regions $fn_graph --seg-filter 0.05 --region-filter 0 --gtf $fn_gtf $fn_bam_all

#$dir/generate_segment_graph $fn_graph --regions $fn_regions --seg-filter 0.05 --region-filter 0 --gtf $fn_gtf $fn_bam_all
#$dir/generate_segment_graph $fn_graph --few-regions --regions $fn_regions --seg-filter 0.05 --region-filter 100 $fn_bam_all

for iter in `seq 1`
do


	##############################	
	# mip 
	# generate expected quantification values for each segment
	# and expected intron counts for each intron;
	# treat bam files as separate samples
	##############################	
	fn_bam_iter=$fn_bam_all
	mip_dir=$out_dir/res_iter${iter}_gtf/
	mkdir -p $mip_dir
	#MAT="/cbio/grlab/share/software/matlab-7.6/bin/matlab -nojvm -nodesktop -nosplash"
	#MAT="matlab -nojvm -nodesktop -nosplash"
	MAT="/cbio/grlab/share/software/matlab/matlab_R2012b/bin/matlab -nojvm -nodesktop -nosplash"
	addpaths="addpath matlab; "
	echo "dbstop error; $addpaths; mip_paths; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir'); exit"
	${MAT} -r "dbstop error; $addpaths mip_paths; transcript_predictions('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir'); exit"
	
	##############################	
	# eval mip
	##############################		
	#${MAT} -r "dbstop error; $addpaths mip_paths; eval_dir('$mip_dir', 'train'); exit"

done
