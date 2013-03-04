#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <math.h>
#include <algorithm>
#include <assert.h>
#include "bam.h"
#include "region.h"
#include "get_reads_direct.h"
#include <fstream>
#include "gtf_tools.h"
#include "math_tools.h"
#include "tools.h"
#include "file_stats.h"

//includes for samtools
////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "razf.h"
#include "bgzf.h"
#include "razf.h"
////////////////////////////////////////////////////////////////////////////////////////////////////
struct Config
{
	vector<char*> bam_files;
	char* fn_gtf;
	char* fn_regions;
	char* fn_out;
	float seg_filter;
	float tss_pval;
	bool split_chr;
	bool strand_specific;
	bool reads_by_chr;
	bool mm_filter;
	int intron_len_filter;
	int filter_mismatch;
	int exon_len_filter;
	int region_filter;
	int gtf_offset;
	int max_junk;
};

int parse_args(int argc, char** argv,  Config* c)
{
	if (argc<2)
	{
			fprintf(stdout, "Usage: %s <fn_out> [options] <fn_bam1> <fn_bam2> ...\n", argv[0]);
			fprintf(stdout, "\n");
			fprintf(stdout, "options:\n");
			fprintf(stdout, "\t--regions\t\t(file name) specify regions flat file (e.g. output of define_regions)\n");
			fprintf(stdout, "\t--few-regions \t\t(flag) load data not for whole chromosome, but for each region separate \n");
			fprintf(stdout, "\t--max-junk \t\t(default 2e7) if flag --few-regions not set, load reads for junk at once\n");
			fprintf(stdout, "\t\t\t\t -> more efficient if there are only a few region to predict on\n");
			fprintf(stdout, "\t--region-filter\t\t(default 1000) discard regions with less reads \n");
			fprintf(stdout, "\t--gtf\t\t\t(file name) specify gtf file\n");
			fprintf(stdout, "\t--gtf-offset\t\t(default 10000) region arround annotated gene to be included in search for alternative transcripts\n");
			fprintf(stdout, "\t--tss-tts-pval\t\t(default 0.0001) p-value cutoff for transcription start and transcription termination site discovery\n");
			fprintf(stdout, "\t--split-chr\t\t(flag) generate one output file for each chromosom\n");
			fprintf(stdout, "\t\n");
			fprintf(stdout, "read related options:\n");
			fprintf(stdout, "\t--min-exonic-len\t(default 0) minimal number of aligned positions on each side of an intron\n");
			fprintf(stdout, "\t--mismatches\t\t(default 10) maximal number of mismatches\n");
			fprintf(stdout, "\t--best\t\t\t(flag) filter for best alignment\n");
			fprintf(stdout, "\t--no-strand \t\t(flag) data are not strand specific \n");
			fprintf(stdout, "\t\t\n");
			return -1;
	}

	// defaults
	c->fn_gtf = NULL;
	c->fn_regions = NULL;	
	c->fn_out = argv[1];	
	c->gtf_offset = 10000;
	c->strand_specific = true;
	c->reads_by_chr = true;
	c->intron_len_filter = 200000;
	c->filter_mismatch = 10;
	c->exon_len_filter = 0;
	c->mm_filter = false;
	c->region_filter = 1000;
	c->seg_filter = 0.05;
	c->tss_pval = 1e-4;
	c->split_chr=false;
	c->max_junk = 2e7;

    for (int i = 2; i < argc; i++)  
    {
    	if (strcmp(argv[i], "--no-strand") == 0)
		{
			c->strand_specific = false;
		}
		else if (strcmp(argv[i], "--split-chr") == 0)
        {
			c->split_chr = true;
        }
	    else if (strcmp(argv[i], "--best") == 0)
        {
			c->mm_filter = true;
        }
	    else if (strcmp(argv[i], "--max-junk") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --max-junk\n") ;
                return -1;
            }
            i++;
			c->max_junk = atoi(argv[i]);

        }
	    else if (strcmp(argv[i], "--gtf-offset") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --gtf-offset\n") ;
                return -1;
            }
            i++;
			c->gtf_offset = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--min-exonic-len") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --min-exonic-len\n") ;
                return -1;
            }
            i++;
			c->exon_len_filter = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--mismatches") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --mismatches\n") ;
                return -1;
            }
            i++;
			c->filter_mismatch = atoi(argv[i]);
        }
    	else if (strcmp(argv[i], "--few-regions") == 0)
		{
			c->reads_by_chr = false;
		}
	    else if (strcmp(argv[i], "--region-filter") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --region-filter\n") ;
                return -1;
            }
            i++;
			c->region_filter = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--regions") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --regions\n") ;
                return -1;
            }
            i++;
			c->fn_regions = argv[i];
        }
	    else if (strcmp(argv[i], "--gtf") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --gtf\n") ;
                return -1;
            }
            i++;
			c->fn_gtf = argv[i];
        }
	    else if (strcmp(argv[i], "--seg-filter") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --seg-filter\n") ;
                return -1;
            }
            i++;
			c->seg_filter = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--tss-tts-pval") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --tss_tts_pval\n") ;
                return -1;
            }
            i++;
			c->tss_pval = atof(argv[i]);
        }
		else
		{
			c->bam_files.push_back(argv[i]);
		}	
	}
	//assert(c->bam_files.size()>0);
	if (!c->fn_regions && !c->fn_gtf)
	{
		printf("either --gtf, or --fn-regions, or both have to be specified\n");
	}
	return 0;
}

void find_tss_and_cleave(Region* region, Config* c)
{
	if (region->transcripts.size()==0)
	{
		printf("no transcript found\n");
		return;
	}

	// collect read starts and end positions
	vector<int> rstarts;
	vector<int> rends;
	for (uint i=0; i<region->reads.size(); i++)
	{
		rstarts.push_back(region->reads[i]->start_pos);
		rends.push_back(region->reads[i]->get_last_position());
	}
	sort(rstarts.begin(), rstarts.end());
	sort(rends.begin(), rends.end());

	int max_len = 500;
	int win=30;
	float best_pval1 = 1.0;
	int best_tss = -1;
	float best_pval2 = 1.0;
	int best_stop = -1;
	for (int i=0; i<region->transcripts.size(); i++)
	{
		// find start
		int start = region->transcripts[i].front().first;
		for (int pos=start-1; pos>region->start && pos>start-max_len; pos-=10)
		{
			int val1 = 0;
			int val2 = 0;
			for (uint s=0; s<rstarts.size(); s++)
			{
				if (rstarts[s]>pos && rstarts[s]<=pos+win)
					val1++;
				else if (rstarts[s]<=pos && rstarts[s]>pos-win)
					val2++;
			}
			float p_val = bino_cfd(val1, val1+val2, 0.25);
			if (p_val<best_pval1)
			{
				best_pval1 = p_val;
				best_tss = pos; 
				printf("tss: pos:%i, p_val: %.16f,  log(pval): %4f val1 %i, val2 %i\n", pos, p_val , log10(p_val), val1, val2);
			}
			if (val1 == 0 || val2 == 0)
			{
				if (best_pval1>c->tss_pval);
					best_tss = pos;
				break;
			}
		}
		// find end of coverage if p_val too large
		if (best_pval1>c->tss_pval)
		{
			printf("tss: did not find significant drop in coverage (pval: %.4f)\n", best_pval1);
		}

		
		// find stop
		int stop = region->transcripts[i].back().second;
		for (int pos=stop+1; pos<region->stop && pos<stop+max_len; pos+=10)
		{
			int val1 = 0;
			int val2 = 0;
			for (uint s=0; s<rends.size(); s++)
			{
				if (rends[s]>pos && rends[s]<=pos+win)
					val1++;
				else if (rends[s]<=pos && rends[s]>pos-win)
					val2++;
			}
			float p_val = bino_cfd(val1, val1+val2, 0.25);
			printf("stop: pos:%i, p_val: %.16f,  log(pval): %4f val1 %i, val2 %i\n", pos, p_val , log10(p_val), val1, val2);
			if (p_val<best_pval2)
			{
				best_pval2 = p_val;
				best_stop = pos; 
				printf("take\n");
			}
			if (val1 == 0 || val2 == 0)
			{
				if (best_pval2>c->tss_pval);
					best_stop = pos;
				break;
			}
		}
		// find end of coverage if p_val too large
		if (best_pval2>c->tss_pval)
		{
			printf("stop: did not find significant drop in coverage (pval: %.4f)\n", best_pval2);
		}

		// add new segments to transcript
		if (region->strand=='+')
		{
			segment* UTR5 = new segment(best_tss, start-1);
			segment* UTR3 = new segment(stop+1, best_stop); 
			region->transcripts[i].push_back(*UTR3);
			region->coding_flag[i].push_back(3);
			region->transcripts[i].insert(region->transcripts[i].begin(), *UTR5);
			region->coding_flag[i].insert(region->coding_flag[i].begin(), 5);
		}
		else
		{
			segment* UTR5 = new segment(stop+1, best_stop); 
			segment* UTR3 = new segment(best_tss, start-1);
			region->transcripts[i].push_back(*UTR5);
			region->coding_flag[i].push_back(5);
			region->transcripts[i].insert(region->transcripts[i].begin(), *UTR3);
			region->coding_flag[i].insert(region->coding_flag[i].begin(), 3);
		}
	}
}
bool compare_chr_and_strand(const Region* reg1, const Region* reg2)
{
	if (reg1->chr_num<reg2->chr_num)
		return true; 
	if (reg1->chr_num>reg2->chr_num)
		return false;

	if (reg1->strand=='+' && reg2->strand=='-')
		return true;
	if (reg1->strand=='-' && reg2->strand=='+')
		return false;
	
	if (reg1->start<reg2->start)
		return true;
	else 
		return false;
}

int main(int argc, char* argv[])
{
	Config c;
	int ret = parse_args(argc, argv, &c);
	if (ret!=0)
		return ret;

	bamFile fd;
	bam_header_t* header;
	if (c.bam_files.size()>0)
	{
		fd = bam_open(c.bam_files[0], "r");
		if (fd==0)
		{
			printf("failed to open bam file: %s\n", c.bam_files[0]);
			return -2;
		}
		header = bam_header_read(fd);
		if (header == 0)
		{
			fprintf(stderr, "[%s] Invalid BAM header.\n", argv[0]);
			return -2;
		}
	}

	FILE* fd_null = fopen("/dev/null", "w");
	FILE* fd_out = fopen(c.fn_out, "w");

	vector<Region*> gtf_regions;
	char strand_prev;
	vector<CRead*>::iterator curr; 
	Region* reg = NULL;
	int last_stop=0;
	char* chr_prev = (char*) "xxx";

	const char* format = determine_format(c.fn_gtf);
	printf("loading regions from %s file: %s\n", format, c.fn_gtf);
	if (strcmp(format, "gtf")==0)
		gtf_regions = parse_gtf(c.fn_gtf);
	else if (strcmp(format, "gff3")==0)
		gtf_regions = parse_gff(c.fn_gtf);
	else
	{
		printf("could not determine format of annotation file: %s\n", c.fn_gtf);
		exit(-1);
	}
	printf("number of regions from gtf file: %i\n", (int) gtf_regions.size());
	printf("config: reads_by_chr: %i\n", c.reads_by_chr);

	if (c.reads_by_chr)
	{
		for (int i=0; i<gtf_regions.size(); i++)
		{
			set_chr_num(gtf_regions[i], header);
		}
		sort(gtf_regions.begin(), gtf_regions.end(), compare_chr_and_strand);
	}

	printf("process gtf regions ... \n");
	for (int i=0; i<gtf_regions.size(); i++)
	{
		printf("\rprocess gtf region %i (%i)", i, (int) gtf_regions.size());
		int chr_len = 1e7;
		if (c.bam_files.size()>0)
		{
			set_chr_num(gtf_regions[i], header);
			chr_len = header->target_len[gtf_regions[i]->chr_num];
		}

		gtf_regions[i]->start = std::max(0, gtf_regions[i]->start-c.gtf_offset);
		gtf_regions[i]->stop = std::min(chr_len-1, gtf_regions[i]->stop+c.gtf_offset);

		// shrink start and stop according to read coverage	
		int num_bam = c.bam_files.size();
		if (!c.strand_specific)
			gtf_regions[i]->strand = '0';
		
		gtf_regions[i]->fd_out = fd_null;
		if (c.reads_by_chr)// more fast if many regions are considered
		{
			bool get_reads = false;
			if (strcmp(gtf_regions[i]->chr, chr_prev)!=0 || (gtf_regions[i]->strand != strand_prev && c.strand_specific) || reg->stop<gtf_regions[i]->stop)
			{
				get_reads = true;
				chr_prev = gtf_regions[i]->chr;
				strand_prev = gtf_regions[i]->strand;
				printf("starting with chr: %s%c\n", chr_prev, strand_prev);
				delete reg;
				reg = new Region(gtf_regions[i]->start, gtf_regions[i]->start, chr_prev, strand_prev);
				set_chr_num(reg, header);
				reg->stop = std::min(gtf_regions[i]->start+c.max_junk, (int) header->target_len[reg->chr_num]);
			}

			if (get_reads)
			{
				printf("get reads from %i bam files for region %s:%i->%i\n", (int)c.bam_files.size(), reg->chr, reg->start, reg->stop);

				// get reads for the large region
				int num_bam = c.bam_files.size();
				if (!c.strand_specific)
					reg->strand = '0';
		
				reg->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);

				printf("sort reads by start position ... \n");
				sort(reg->reads.begin(), reg->reads.end(), CRead::compare_by_start_pos);
				printf("done\n");

				curr = reg->reads.begin();
				last_stop=0;
			}

			// get reads from chromosom region
			// if regions overlapp, decrement the read pointer accordingly
			if (gtf_regions[i]->start<last_stop)
			{
				if (curr == reg->reads.end() && curr != reg->reads.begin())
					curr--;
				while (curr != reg->reads.begin() && (*curr)->start_pos>=gtf_regions[i]->start)
					curr--;
			}
			last_stop = gtf_regions[i]->stop;

			//printf("add reads to gtf_region(%s%c:%i-%i)\n", gtf_regions[i]->chr, gtf_regions[i]->strand, gtf_regions[i]->start, gtf_regions[i]->stop);
			int num_reads = 0;
			while (curr != reg->reads.end())
			{
				if ((*curr)->start_pos>=gtf_regions[i]->stop)
					break;

				if ((*curr)->start_pos>=gtf_regions[i]->start && (*curr)->get_last_position()<gtf_regions[i]->stop)
				{
					gtf_regions[i]->reads.push_back(*curr);
					num_reads++;
				}
				curr++;
			}
			//printf("added %i reads\n", num_reads);
		}
		else if (c.bam_files.size()>0)
		{
			gtf_regions[i]->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);
		}

		// find tss and cleavage site and 
		// add new UTR exons to each transcript
		find_tss_and_cleave(gtf_regions[i], &c);

		// write to gff/gtf file
		write_gff(fd_out, gtf_regions[i], "MiTie");

		if (c.reads_by_chr)
			gtf_regions[i]->reads.clear();
		else
			gtf_regions[i]->clear_reads();
		delete[] gtf_regions[i]->coverage;
	}
	printf("done\n");

	// cleanup
	for (int i=0; i<gtf_regions.size(); i++)
		delete gtf_regions[i];
	
	fclose(fd_null);
	delete reg;
	bam_close(fd);
	bam_header_destroy(header);

	return 0;
}
