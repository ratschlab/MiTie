#!/bin/bash
#
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_gtf_sample${sample} /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_sample${sample} /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_fixes_sample${sample} /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_filter_sample${sample} /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
# e.g. export sample=1 && ./mmr_mip/mip_mmr.sh /fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_no_shrink_sample${sample} /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam $sample
#
#

if [ -z $1 ]; then
	sample=2
else
	sample=$1
fi
#out_dir=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_release_sample${sample}
#out_dir=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_tss0.01_sample${sample}
#out_dir=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_enum5_sample${sample}
#out_dir=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_enum5_1e4_repeat_sample${sample}
#out_dir=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_seg_filter0.01_sample${sample}
#out_dir=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions40000_no_junc_sample${sample}
out_dir=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_GP_sample${sample}
#out_dir=~/tmp
#fn_bam_all=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam

for s in `seq 1 $sample`; do
	fn_bam_all="$fn_bam_all /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${s}_merged_err_1.no_junc.sorted.paired.bam"
	#fn_bam_all="$fn_bam_all /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.new.sorted.paired.bam"
done

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

fn_gtf=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr_max_trans.gtf
#fn_gtf=~/tmp/sample.gtf
if [ ! -f $fn_graph ]
then
	echo
	echo generate graph on bam file $fn_bam_all
	echo
	#opts="--mismatches 3 --min-exonic-len 5"
	$dir/generate_segment_graph $fn_graph.tmp $opt --split-chr --regions $fn_regions --max-junk 200000 --gtf-offset 50000 --seg-filter 0.05 --region-filter 50 --tss-tts-pval 0.0001 --gtf $fn_gtf $fn_bam_all
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

for iter in `seq 1 2`
do
	echo $iter	
	fn_bam_iter=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.no_junc.mmr.iter$iter.bam

	##############################	
	# mmr
	# generate bam file for iteration 
	##############################	
	if [ "$iter" -gt "1" ]; then
		mmr=$HOME/svn/projects/rnageeq/mm_resolve/threaded_oop_mip/mmr
		INPUT=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${sample}_merged_err_1.no_junc.ID_sorted.paired.bam
		seg_list_iter=$out_dir/segment_list_$(($iter - 1)).txt
		LOSS_FILE=/fml/ag-raetsch/home/akahles/git/software/RNAgeeq/mm_resolve/threaded_oop_mip/poisson_3.flat
		THREADS=6
		OUTFILE=$HOME/tmp/mmr_iter$iter.bam
		ITER=1
		$mmr -o $OUTFILE -t $THREADS -z -S -I $ITER -f -F 1 -p -b -m -s $seg_list_iter -l $LOSS_FILE -r 75 -v $INPUT || exit -1;
		
		samtools sort $OUTFILE ${fn_bam_iter%.bam} && samtools index $fn_bam_iter
	fi

	if [ "$iter" -lt "2" ]; then
		fn_bam_iter=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.no_junc.sorted.paired.best.bam
		for s in `seq 2 $sample`; do
			fn_bam_iter="$fn_bam_iter /fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias${s}_merged_err_1.no_junc.sorted.paired.best.bam"
		done
	fi
	if [ ! -f $fn_bam_iter ]; then
		exit -1;
	fi
	echo $fn_bam_iter

	##############################	
	# mip 
	# generate expected quantification values for each segment
	# and expected intron counts for each intron;
	# treat bam files as separate samples
	##############################	
	mip_dir=$out_dir/res_iter${iter}_gtf/
	mkdir -p $mip_dir
	MAT="/fml/ag-raetsch/share/software/matlab-7.6/bin/matlab -nojvm -nodesktop -nosplash"
	#MAT="matlab -nojvm -nodesktop -nosplash"
	addpaths="addpath matlab; "
	echo "dbstop error; $addpaths; denovo('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir'); exit"
	${MAT} -r "dbstop error; $addpaths mip_paths; denovo('$fn_graph', {'`echo $fn_bam_iter | sed "s/ /','/g"`'}, '$mip_dir'); exit"
	#echo "dbstop error; $addpaths mip_paths; denovo('$fn_graph', '$fn_bam_iter', '$mip_dir'); exit"
	#${MAT} -r "dbstop error; $addpaths mip_paths; denovo('$fn_graph', '$fn_bam_iter', '$mip_dir'); exit"
	
	##############################	
	# eval mip
	##############################		
	${MAT} -r "dbstop error; $addpaths mip_paths; eval_dir('$mip_dir', 'train'); exit"

	##############################		
	# write segment and intron list	
	##############################		
	seg_list_iter=$out_dir/segment_list_$iter.txt
	${MAT} -r "dbstop error; $addpaths mip_paths; write_seg_list('$mip_dir', '$seg_list_iter'); exit"
done
