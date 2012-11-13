#!/bin/bash

mkdir testdata
fn_gtf_orig=/fml/ag-raetsch/nobackup/projects/mip/human_annotation/Homo_sapiens.GRCh37.68.gtf
fn_gtf=testdata/Homo_sapiens.GRCh37.68.chr20.gtf
grep ^20 $fn_gtf_orig | sed 's/^20/chr20/' | grep ENSG00000125841 | grep "[[:space:]]exon[[:space:]]" > $fn_gtf

#region: RAS 
REG=chr20:19,915,743-19,983,103
REG=chr20:327527-334147


bam_orig="/fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE/wgEncodeCshlLongRnaSeqGm12878CellPapAlnRep1.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE/wgEncodeCshlLongRnaSeqK562CellPapAlnRep1.bam /fml/ag-raetsch/nobackup/projects/sequencing_runs/human/ENCODE/wgEncodeCshlLongRnaSeqHepg2CellLongnonpolyaAlnRep2.bam"

for f in $bam_orig; do
	res=`basename $f`
	res=testdata/${res%.bam}.toy.bam
	echo $res
	samtools view $f $REG -h | samtools view -Sb - > $res 
	samtools index $res
done
