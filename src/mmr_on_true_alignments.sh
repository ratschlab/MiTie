#!/bin/bash
	
seg_list_iter=/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_quant_mm0_nb_true/segment_list.txt
fn_bam_iter=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.no_junc.mmr_anno_true_new2.bam
#fn_bam_iter=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.no_junc.mmr.bam

mmr=$HOME/svn/projects/rnageeq/mm_resolve/threaded_oop_mip/mmr
INPUT=/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.no_junc.ID_sorted.paired.bam
LOSS_FILE=/fml/ag-raetsch/home/akahles/git/software/RNAgeeq/mm_resolve/threaded_oop_mip/poisson_3.flat
THREADS=6
OUTFILE=$HOME/tmp/mmr_$RANDOM.bam
ITER=5
ZERO_SEGMENTS="-z"
OPTS="-t $THREADS $ZERO_SEGMENTS --strand-specific -I $ITER --pre-filter --filter-dist 0 -p -b --insert-size 100000 --insert-dev 2.0"
#$mmr -o $OUTFILE $OPTS -v $INPUT
$mmr -o $OUTFILE $OPTS -m -s $seg_list_iter -l $LOSS_FILE -r 75 -v $INPUT
	
samtools sort $OUTFILE ${fn_bam_iter%.bam} && samtools index $fn_bam_iter

rm $OUTFILE
