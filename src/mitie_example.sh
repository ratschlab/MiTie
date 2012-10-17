#!/bin/bash
#

#fn_bam_brain=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.mmr.sorted.bam
#fn_bam_uhr=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.mmr.sorted.bam
#fn_bam_brain_all=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.Brain.sorted.bam
#fn_bam_uhr_all=/fml/ag-raetsch/nobackup2/projects/mip_spladder/alignments/human/ILM_S.UHR.sorted.bam

#out_dir=/fml/ag-raetsch/nobackup/projects/mip/chris_mason/MiTie

orig_data=/fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE/
#fn_bam="$orig_data/*merge.chr20.bam"
fn_bam="/fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqA549CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqAg04450CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqBjCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqCd20CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqGm12878CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqGm12878CytosolPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqGm12878NucleusPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqH1hescCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHelas3CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHelas3CytosolPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHelas3NucleusPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHepg2CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHepg2CytosolPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHepg2NucleusPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHmecCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHsmmCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqHuvecCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqK562CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqK562CytosolPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqMcf7CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqMonocd14CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqNhekCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqNhlfCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqSknshraCellPapAln.merge.chr20.bam"
#fn_bam="$orig_data/wgEncodeCshlLongRnaSeqSknshraCellPapAln.merge.chr20.bam"

fn_bam_few="/fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqA549CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqAg04450CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqBjCellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqCd20CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqGm12878CellPapAln.merge.chr20.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE//wgEncodeCshlLongRnaSeqGm12878CytosolPapAln.merge.chr20.bam"

out_dir=/fml/ag-raetsch/nobackup/projects/mip/ENCODE/MiTie_all_max1e6
mkdir -p $out_dir
fn_regions=$out_dir/regions.flat

dir=`dirname $0`

if [ ! -f $fn_regions ]
then
	echo run define_regions => $fn_regions
	#valgrind --tool=cachegrind $dir/define_regions $fn_bam_all -o $fn_regions
	#valgrind --leak-check=full $dir/define_regions $fn_bam_all -o $fn_regions
	echo $dir/define_regions $fn_bam_few -o $fn_regions --shrink --cut-regions --max-len 100000
	gdb  $dir/define_regions
fi

# create segment graph and store in file
dir_graph=${out_dir}/graphs
dir_log=${out_dir}/logs
dir_reg=${out_dir}/regions
mkdir -p $dir_graph
mkdir -p $dir_log
mkdir -p $dir_reg

fn_graph=${out_dir}/graphs_xxx
fn_graph=

if [ ! -f $fn_graph ]
then
	echo
	echo generate graph on bam file $fn_bam
	echo
	#opts="--mismatches 3 --min-exonic-len 5"
	#echo $dir/generate_segment_graph $fn_graph.tmp $opts --split-chr --regions $fn_regions --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001 $fn_bam
	#$dir/generate_segment_graph $fn_graph.tmp $opts --split-chr --regions $fn_regions --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001 $fn_bam
	echo 

	cnt=0
	while read line         
	do
		mem_req=5
		cnt=$(($cnt+1))
		fn_reg=$dir_reg/region$cnt.flat
		fn_gr=$dir_graph/graph$cnt.bin
		fn_log=$dir_log/log$cnt
		if [ ! -f $fn_gr ]; then
			echo $line > ${fn_reg}
			echo "$dir/generate_segment_graph $fn_gr $opts --regions $fn_reg --few-regions --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001 $fn_bam" | qsub -o $fn_log -cwd -l h_vmem=${mem_req}G -j y -p 31 -N gen_graph
		fi
		num_jobs=`qstat | wc -l`
		while [ "$num_jobs" -gt "2500" ]; do
			echo "too many jobs in queue ($num_jobs), waiting..."
			sleep 120
			num_jobs=`qstat | wc -l`
		done
	done < $fn_regions

	#valgrind --leak-check=full $dir/generate_segment_graph $fn_graph.tmp --few-regions $opts --regions $fn_regions --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001 $fn_bam
	#echo $dir/generate_segment_graph $fn_graph.tmp $opts --split-chr --regions $fn_regions --max-junk 200000 --seg-filter 0.05 --region-filter 100 --tss-tts-pval 0.0001 $fn_bam
	#gdb $dir/generate_segment_graph


	if [ ! "$?" -eq "0" ]; then
		echo generate_segment_graph return: $?
		exit 0
	fi
	mv ${fn_graph}.tmp $fn_graph
else
	echo $fn_graph exists
fi


##############################	
##############################	
mip_dir=$out_dir/MiTie_pred/
mkdir -p $mip_dir
MAT="/fml/ag-raetsch/share/software/matlab-7.6/bin/matlab -nojvm -nodesktop -nosplash"
addpaths="addpath matlab; "
for fn_gr in $dir_graph/graph*.bin; do
	if [ -s $fn_gr ]; then 
		echo ls -l $fn_gr
		${MAT} -r "dbstop error; $addpaths mip_paths; denovo_encode('$fn_gr', {'`echo $fn_bam | sed "s/ /','/g"`'}, '$mip_dir'); exit"
	fi
done


