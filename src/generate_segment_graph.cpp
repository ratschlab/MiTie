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
#include "tools.h"

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
	int gtf_offset;
	bool strand_specific;
	char* fn_regions;
	char* fn_out;
	bool reads_by_chr;
	int intron_len_filter;
	int filter_mismatch;
	int exon_len_filter;
	bool mm_filter;
	int region_filter;
	float seg_filter;
	float tss_pval;

};

int parse_args(int argc, char** argv,  Config* c)
{
	if (argc<4)
	{
			fprintf(stdout, "Usage: %s <fn_out> [options] <fn_bam1> <fn_bam2> ...\n", argv[0]);
			fprintf(stdout, "options:\n");
			fprintf(stdout, "--no-strand \t\t (flag) data are not strand specific \n");
			fprintf(stdout, "--few-regions \t\t (flag) load data not for whole chromosome, but for each region separate \n");
			fprintf(stdout, "\t\t\t\t -> (flag) more efficient if there are only a few region to predict on\n");
			fprintf(stdout, "--region-filter \t\t (default 1000) discard regions with less reads \n");
			fprintf(stdout, "--regions \t\t (file name) specify regions flat file (e.g. output of define_regions)\n");
			fprintf(stdout, "--gtf \t\t (file name) specify gtf file\n");
			fprintf(stdout, "--seg-filter \t\t (default 0.05) segments are filtered out if fraction of covered nucleotides is below\n");
			fprintf(stdout, "--tss-tts-pval \t\t (default 0.0001) p-value cutoff for transcription start and transcription termination site discovery\n");
			fprintf(stdout, "--min-exonic-len\t\t\n");
			fprintf(stdout, "--mismatches\t\t\n");
			fprintf(stdout, "--gtf-offset\t\t\n");
			fprintf(stdout, "\t\t\n");
			return -1;
	}

	// defaults
	c->fn_gtf = NULL;
	c->gtf_offset = 0;
	c->fn_regions = NULL;	
	c->fn_out = argv[1];	
	c->strand_specific = true;
	c->reads_by_chr = true;
	c->intron_len_filter = 200000;
	c->filter_mismatch = 10;
	c->exon_len_filter = 0;
	c->mm_filter = false;
	c->region_filter = 1000;
	c->seg_filter = 0.05;
	c->tss_pval = 1e-4;

    for (int i = 2; i < argc; i++)  
    {
    	if (strcmp(argv[i], "--no-strand") == 0)
		{
			c->strand_specific = false;
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
	assert(c->bam_files.size()>0);
	if (!c->fn_regions && !c->fn_gtf)
	{
		printf("either --gtf, or --fn-regions, or both have to be specified\n");
	}
	return 0;
}
int main(int argc, char* argv[])
{
	Config c;
	int ret = parse_args(argc, argv, &c);
	if (ret!=0)
		return ret;

	std::ofstream* ofs = new std::ofstream();
	ofs->open(c.fn_out, std::ios::binary);

	if (!ofs)
	{
		fprintf(stderr, "[%s] Could not open file: %s for writing\n", argv[0], c.fn_out);
		return -2;
	}

	bamFile fd = bam_open(c.bam_files[0], "r");
	if (fd==0)
	{
		printf("failed to open bam file: %s\n", c.bam_files[0]);
		return -2;
	}
	bam_header_t* header = bam_header_read(fd);
	if (header == 0)
	{
		fprintf(stderr, "[%s] Invalid BAM header.\n", argv[0]);
		return NULL;
	}

	vector<Region*> regions;
	if (c.fn_regions)
	{
		printf("loading regions from flat file: %s\n", c.fn_regions);
		regions = parse_regions(c.fn_regions);
		printf("number of regions from flat file: %i\n", (int) regions.size());
	}
	int num_reg = regions.size();

	vector<Region*> gtf_regions;
	if (c.fn_gtf)
	{
		printf("loading regions form gtf file: %s\n", c.fn_gtf);
		gtf_regions = parse_gtf(c.fn_gtf);
		for (int i=0; i<gtf_regions.size(); i++)
		{
			gtf_regions[i]->start = std::max(0, gtf_regions[i]->start-c.gtf_offset);
			gtf_regions[i]->stop = gtf_regions[i]->stop+c.gtf_offset;
		}
		printf("number of regions from gtf file: %i\n", (int) gtf_regions.size());
	}
	int num_gtf = gtf_regions.size();
	// add gtf regions
	for (int i=0; i<gtf_regions.size(); i++)
	{
		regions.push_back(gtf_regions[i]);
	}

	if (num_reg>0 && num_gtf>0)
	{
		printf("merging %i flat file regions with %i gtf regions\n", num_reg, num_gtf);
		// self overlap
		bool change = true;
		int iter = 0;
		while (change)
		{
			printf("merge iteration: %i size:%i\n", iter++, (int) regions.size());
			if (iter >10)
				break;
			change = false;	
			vector<vector<int> > ov_list = region_overlap(regions, regions);
			for (int i=0; i<regions.size(); i++)
			{
				for (int j=0; j<ov_list[i].size(); j++)
				{
					// self overlap
					if (ov_list[i][j]==i)
						continue;

					if (strcmp(regions[i]->chr, regions[ov_list[i][j]]->chr)!=0 || regions[i]->strand != regions[ov_list[i][j]]->strand)
						continue;

					if (regions[ov_list[i][j]]->start==-1 || regions[i]->start==-1)
					{
						continue;
					}
					change = true;

					//printf("reg: %s%c:%i->%i\n", regions[i]->chr, regions[i]->strand, regions[i]->start, regions[i]->stop);
					//printf("ov: %s%c:%i->%i\n", regions[ov_list[i][j]]->chr, regions[ov_list[i][j]]->strand, regions[ov_list[i][j]]->start, regions[ov_list[i][j]]->stop);
					// merge regions
					regions[i]->start = std::min(regions[i]->start, regions[ov_list[i][j]]->start);
					regions[i]->stop = std::max(regions[i]->stop, regions[ov_list[i][j]]->stop);
					for (int k=0; k<regions[ov_list[i][j]]->transcripts.size(); k++)
						regions[i]->transcripts.push_back(regions[ov_list[i][j]]->transcripts[k]);

					regions[ov_list[i][j]]->start = -1;

				}
			}
			// remove merged regions
			vector<Region*> tmp;
			for (int i=0; i<regions.size(); i++)
				if (regions[i]->start>-1)
					tmp.push_back(regions[i]);

			regions = tmp;
		}
	}

	FILE* fd_null = fopen("/dev/null", "w");

	vector<int> bias_vector(100, 0);
	//filter regions
	int cnt = 0;
	char* chr_prev = (char*) "xxx";
	char strand_prev;
	vector<CRead*>::iterator curr; 
	Region* reg = NULL;
	for (int i=0; i<regions.size(); i++)
	{
		if (c.reads_by_chr)// more fast if many regions are considered
		{
			if (strcmp(regions[i]->chr, chr_prev)!=0 || (regions[i]->strand != strand_prev && c.strand_specific))
			{
				chr_prev = regions[i]->chr;
				strand_prev = regions[i]->strand;
				printf("starting with chr: %s%c\n", chr_prev, strand_prev);
				delete reg;
				reg = new Region(regions[i]);
				reg->start = 1;
				set_chr_num(reg, header);
				reg->stop = header->target_len[reg->chr_num];
				printf("get reads for region %s:%i->%i\n", reg->chr, reg->start, reg->stop);

				// get reads for the large region
				int num_bam = c.bam_files.size();
				if (!c.strand_specific)
					reg->strand = '0';
			
				reg->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);
				curr = reg->reads.begin();
			}

			// get reads from chromosom region
			while (curr<reg->reads.end())
			{
				//printf("%p, %i read start\n", (*curr), (*curr)->start_pos);
				if ((*curr)->start_pos>=regions[i]->stop)
					break;
				if ((*curr)->get_last_position()>regions[i]->start)
					regions[i]->reads.push_back(*curr);
				curr++;
			}
		}
		else
		{
			int num_bam = c.bam_files.size();
			if (!c.strand_specific)
				regions[i]->strand = '0';
			
			regions[i]->fd_out = fd_null;
			regions[i]->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);
		}
		//printf("%i reads in region: %s:%i-%i\n", (int)regions[i]->reads.size(), regions[i]->chr, regions[i]->start, regions[i]->stop);
		printf("\r region %i(%i)", i, (int)regions.size());

		int num_reads=0;
		for (int r=0; r<regions[i]->reads.size(); r++)
			if (regions[i]->reads[r]->mismatches<1)
				num_reads++;
		
		bool is_annotated = regions[i]->transcripts.size()>0;
		if (num_reads<=c.region_filter && !is_annotated)
		{
			if (c.reads_by_chr)
				regions[i]->reads.clear();
			else
				regions[i]->clear_reads();
			continue;
		}
		cnt++;	
		//regions[i]->fd_out = stdout;
		regions[i]->fd_out = fd_null;
		regions[i]->generate_segment_graph(c.seg_filter, c.tss_pval);
		//regions[i]->add_bias_counts(&bias_vector);

		// write region in binary file
		regions[i]->write_binary(ofs);

		// cleanup
		if (c.reads_by_chr)
			regions[i]->reads.clear();// dont delete the reads
		else
			regions[i]->clear_reads();

		//if (cnt>200)
		//	break;
	}
	printf("read %i regions from file\n", (int)regions.size());
	printf("%i regions after filtering\n", cnt);

	// cleanup
	for (int i=0; i<regions.size(); i++)
		delete regions[i];
	
	delete reg;
	//fclose(fd_out);
	ofs->close();
	bam_close(fd);
	bam_header_destroy(header);
}
