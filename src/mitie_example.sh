#!/bin/bash
#

#fn_bam_brain=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.mmr.sorted.bam
#fn_bam_uhr=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.mmr.sorted.bam
#fn_bam_brain_all=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.sorted.bam
#fn_bam_uhr_all=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.sorted.bam

#out_dir=/fml/ag-raetsch/nobackup/projects/mip/chris_mason/MiTie

fn_bam="testdata/*toy.bam"

out_dir=testdata/results
mkdir -p $out_dir
fn_regions=$out_dir/regions.flat

dir=`dirname $0`

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

echo 0

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
source config.sh
echo running matlab: $MATLAB_BIN
MAT="$MATLAB_BIN -nojvm -nodesktop -nosplash"
addpaths="addpath matlab; "

# without annotation
mip_dir=$out_dir/MiTie_pred/
mkdir -p $mip_dir
eta1=1.50
eta2=0.00
lambda=0
quantify=0

${MAT} -r "dbstop error; $addpaths; mip_paths; C.num_transcripts = 5; transcript_predictions('$fn_graph', {'`echo $fn_bam | sed "s/ /','/g"`'}, '$mip_dir', C, 1:$num_graphs , $quantify, $eta1, $eta2, $lambda); exit"

# with annotation
mip_dir_gtf=$out_dir/MiTie_pred_gtf/
mkdir -p $mip_dir_gtf
${MAT} -r "dbstop error; $addpaths; mip_paths; C.num_transcripts = 5; transcript_predictions('$fn_graph_gtf', {'`echo $fn_bam | sed "s/ /','/g"`'}, '$mip_dir_gtf', C, 1:$num_graphs , $quantify, $eta1, $eta2, $lambda); exit"

# collect predictions and write gtf file
##############################	
add_weights=0;
mmr=1;
write_gtf=1;

# without annotation
fn_genes_mat="$mip_dir/res_genes.mat";
${MAT} -r "dbstop error; $addpaths mip_paths; collect_results('$mip_dir', '$fn_genes_mat', $add_weights, $mmr, $write_gtf); exit"

# with annotation
fn_genes_mat="$mip_dir_gtf/res_genes.mat";
${MAT} -r "dbstop error; $addpaths mip_paths; collect_results('$mip_dir_gtf', '$fn_genes_mat', $add_weights, $mmr, $write_gtf); exit"

echo you can find the resulting transcript prediction in $mip_dir/res_genes.gtf and $mip_dir_gtf/res_genes.gtf

