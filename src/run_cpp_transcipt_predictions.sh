


if [ -z $1 ]; then
	sample=1
else
	sample=$1
fi
eta1=1.00
eta2=0.10
lambda=3

fn_bam_all=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired_200000_4_5.bam 
out_dir=/cbio/grlab/nobackup/projects/mip/human_sim/cpp_mip_quant_sample${sample}_eta1_${eta1}_eta2_${eta2}_lambda_${lambda}
mkdir -p $out_dir

dir=`dirname $0`

# create segment graph and store in file
fn_graph=${out_dir}/graph_gtf.h5

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

h5dump --dataset=Graph_meta_info  $fn_graph

#valgrind --leak-check=full $dir/generate_segment_graph $fn_regions $fn_graph --seg-filter 0.05 --region-filter 0 --gtf $fn_gtf $fn_bam_all


##############################	
# mip 
# generate expected quantification values for each segment
# and expected intron counts for each intron;
# treat bam files as separate samples
##############################	
fn_quant=$out_dir/quant.txt
fn_gff=$out_dir/transcripts.gtf

if [ ! -f $fn_quant.asdf ]
then 
	num_trans=0 # number of additional transcripts 
	order=1
	./transcript_prediction $fn_graph $fn_bam_all $mip_dir --max-num-trans $num_trans --param-eta1 $eta1 --param-eta2 $eta2 --param-lambda $lambda --C-intron 10.0 --C-num-trans 100.0 --fn-quant $fn_quant --fn-out $fn_gff --order $order
else
	echo $fn_quant exists
fi

fn_quant_anno=/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr1_quant.txt
./eval_quant  $fn_quant_anno $fn_quant


