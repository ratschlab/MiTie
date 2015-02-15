#!/bin/bash
#
<<<<<<< HEAD

#fn_bam_brain=/cbio/grlab/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.mmr.sorted.bam
#fn_bam_uhr=/cbio/grlab/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.mmr.sorted.bam
#fn_bam_brain_all=/cbio/grlab/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.sorted.bam
#fn_bam_uhr_all=/cbio/grlab/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.sorted.bam

#out_dir=/cbio/grlab/nobackup/projects/mip/chris_mason/MiTie
=======
source config.sh
>>>>>>> ea98fc4555a5b6a523e8e025b582a384f09c815a

fn_bam="testdata/*toy.bam"
out_dir=testdata/results
fn_regions=$out_dir/regions.flat
mkdir -p $out_dir
rm $out_dir/*

dir=`dirname $0`

# define regions for prediction
##############################	
if [ ! -f $fn_regions ]
then
	echo run define_regions => $fn_regions
	echo $dir/define_regions $fn_bam_few -o $fn_regions --shrink --cut-regions --max-len 100000
	$dir/define_regions $fn_bam -o $fn_regions --shrink --cut-regions --max-len 100000
fi

# create segment graph and store in file
##############################	
fn_graph=${out_dir}/graphs.h5
opts="--few-regions --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001"

echo
echo generate graph on bam file $fn_bam without annotation
echo

$dir/generate_segment_graph ${fn_graph}.tmp $opts --regions $fn_regions $fn_bam
mv ${fn_graph}.tmp $fn_graph
num_graphs=`h5dump --dataset=Graph_meta_info $fn_graph  | grep -A 5 num_graphs | tail -n 1`
echo saved $num_graphs graphs to file


echo
echo generate graph on bam file $fn_bam with annotation
echo
fn_graph_gtf=${out_dir}/graphs_gtf.h5
fn_gtf=testdata/Homo_sapiens.GRCh37.68.chr20.gtf

$dir/generate_segment_graph ${fn_graph_gtf}.tmp $opts --regions $fn_regions --gtf $fn_gtf $fn_bam
mv ${fn_graph_gtf}.tmp $fn_graph_gtf
num_graphs_gtf=`h5dump --dataset=Graph_meta_info $fn_graph_gtf  | grep -A 5 num_graphs | tail -n 1`
echo saved $num_graphs_gtf graphs to file


# perform transcript predictions
##############################	

# without annotation
fn_result=$out_dir/result.gtf
fn_quant=$out_dir/quant.txt
eta1=1.00
eta2=0.10
lambda=0
quantify=0
order=1 # order of the polynom to approximate the loss function (1 or 2)
num_trans=1

echo ./transcript_prediction $fn_graph $fn_bam --max-num-trans $num_trans --param-eta1 $eta1 --param-eta2 $eta2 --param-lambda $lambda --C-intron 100.0 --C-exon 1.0 --C-num-trans 10.0 --fn-quant $fn_quant --fn-out $fn_result --order $order
#gdb --args 
time ./transcript_prediction $fn_graph $fn_bam --max-num-trans $num_trans --param-eta1 $eta1 --param-eta2 $eta2 --param-lambda $lambda --C-intron 100.0 --C-exon 1.0 --C-num-trans 10.0 --fn-quant $fn_quant --fn-out $fn_result --order $order

# with annotation
fn_result_gtf=$out_dir/result_gtf.gtf
fn_quant=$out_dir/quant_gtf.txt
./transcript_prediction $fn_graph_gtf $fn_bam --max-num-trans $num_trans --param-eta1 $eta1 --param-eta2 $eta2 --param-lambda $lambda --C-intron 100.0 --C-exon 1.0 --C-num-trans 10.0 --fn-quant $fn_quant --fn-out $fn_result_gtf --order $order


