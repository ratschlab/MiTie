#include <stdio.h>
#include <mex.h>
#include <assert.h>
#include <vector>
  using std::vector;
#include "region.h"
#include "get_var.h"
#include <fstream>

/*
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs<2)
	{
		mexErrMsgTxt("Usage: load_regions_bin(fname, num_genes, fn_bam1, fn_bam2, ...)");	
	}

	char* fname = get_string(prhs[0]);
	
	int num_return = get_int(prhs[1]);
	static std::ifstream ifs;
	if (num_return<0 && ifs && ifs.is_open())
	{
		printf("load_regions_bin: closing binary file\n"); 
		ifs.close();
		return;
	}

	vector<char*> bam_files;
	for (int i=2; i<nrhs; i++)
		bam_files.push_back(get_string(prhs[i]));
		
	printf("load_regions_bin: found %i bam files\n", (int)bam_files.size());
	for (int i=0; i<bam_files.size(); i++)
		printf("\t%s\n", bam_files[i]);

	if (!ifs.is_open())
	{
		printf("load_regions_bin: opening binary file\n"); 
		ifs.open(fname, std::ios::binary);
	}
	FILE* fd_null = fopen("/dev/null", "w");

	vector<Region*> regions; 
	int cnt = 0;
	//while (!ifs.eof()&& cnt<num_return)
	while (ifs.good()&& cnt<num_return)
	{
		if ((cnt++)%100==0)
			printf("\r %i", cnt);
		Region* r = new Region();
		int ret = r->read_binary(&ifs);
		// TODO load reads and compute coverage and intron counts
		if (!ifs.eof() && ret==0)
			regions.push_back(r);
		else
		{
			printf("could not read from file %s, %i\n", fname, ret);
			delete r;
			break;
		}

		if (bam_files.size()>0)
		{
			//printf("load_regions_bin: bam_files[0]: %s\n", bam_files[0]);
			int num_bam_files = bam_files.size();
			int intron_len_filter = 200000;
			int filter_mismatch = 10;
			int exon_len_filter = 0;
			bool mm_filter = false;
			//printf("chr: %s\n", r->chr);
			r->fd_out = fd_null;
			r->get_reads(&bam_files[0], num_bam_files, intron_len_filter, filter_mismatch, exon_len_filter, mm_filter);
			//printf("admat.size():%i  segments.size():%i\n", (int) r->segments.size(), (int) r->admat.size());
			r->update_coverage_information();
			r->compute_seg_cov();
			r->compute_pair_mat();
			//printf("admat.size():%i  segments.size():%i\n", (int) r->segments.size(), (int) r->admat.size());
		}

	}
	//ifs.close();
	//r->write_segment_graph(stdout);

	const char *field_names[] = {"chr", "strand", "start", "stop", "segments", "coverage", "seg_admat", "transcripts", "pair_mat"};
	int nof_fields = 9;

	int num = regions.size();
	printf("\nparsed %i regions from file\n", num);

	mwSize dims[2] = {1, num};
	plhs[0] = mxCreateStructArray(2, dims, nof_fields, field_names);

	// get index of fields
	int chr_field = mxGetFieldNumber(plhs[0], "chr");
	int strand_field = mxGetFieldNumber(plhs[0], "strand");
	int start_field = mxGetFieldNumber(plhs[0], "start");
	int stop_field = mxGetFieldNumber(plhs[0], "stop");
	int segments_field = mxGetFieldNumber(plhs[0], "segments");
	int coverage_field = mxGetFieldNumber(plhs[0], "coverage");
	int seg_admat_field = mxGetFieldNumber(plhs[0], "seg_admat");
	int transcripts_field = mxGetFieldNumber(plhs[0], "transcripts");
	int pair_field = mxGetFieldNumber(plhs[0], "pair_mat");


	cnt = 0;
	for (int j=0; j<num; j++)
	{
		if ((cnt++)%100==0)
			printf("\r %i\t(%i)", cnt, num);

		// chromosome
		mxSetFieldByNumber(plhs[0], j, chr_field, mxCreateString(regions[j]->chr));

		// strand
		char strand[2];
		strand[0] = regions[j]->strand;
		strand[1] = '\0';
		mxSetFieldByNumber(plhs[0], j, strand_field, mxCreateString(strand));

		// start 
		mxArray* start = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(start) = regions[j]->start;
		mxSetFieldByNumber(plhs[0], j, start_field, start);

		// stop
		mxArray* stop = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(stop) = regions[j]->stop;
		mxSetFieldByNumber(plhs[0], j, stop_field, stop);

		// segments
		int num_seg = regions[j]->segments.size();
		mxArray* segments = mxCreateDoubleMatrix(2, num_seg, mxREAL);
		double* seg = mxGetPr(segments);
		for (int i=0; i<num_seg; i++)
		{
			seg[2*i] = (double) regions[j]->segments[i].first;
			seg[2*i+1] = (double) regions[j]->segments[i].second;
		}
		mxSetFieldByNumber(plhs[0], j, segments_field, segments);

		// coverage
		assert(regions[j]->seg_cov.size()==num_seg);
		mxArray* coverage = mxCreateDoubleMatrix(1, num_seg, mxREAL);
		double* cov = mxGetPr(coverage);
		for (int i=0; i<num_seg; i++)
			cov[i] = regions[j]->seg_cov[i];
		mxSetFieldByNumber(plhs[0], j, coverage_field, coverage);

		// segment adjacency matrix
		int len = regions[j]->admat.size();
		//printf("admat.size(): %i, %i\n", len, seg_admat_field);
		mxArray* seg_admat = mxCreateDoubleMatrix(len, len, mxREAL);
		double* mat = mxGetPr(seg_admat);
		for (int i=0; i<len; i++)
		{
			for (int k=0; k<len; k++)
			{
				mat[len*i+k] = regions[j]->admat[i][k];
			}
		}
		mxSetFieldByNumber(plhs[0], j, seg_admat_field, seg_admat);

		//transcripts
		int num_trans = regions[j]->transcript_paths.size();
		mxArray* trans_mat = mxCreateDoubleMatrix(num_trans, num_seg, mxREAL);
		double* trm = mxGetPr(trans_mat);
		for (int i=0; i<num_trans; i++)
		{
			for (int k=0; k<regions[j]->transcript_paths[i].size(); k++)
			{
				int p = regions[j]->transcript_paths[i][k];
				if (p>=num_seg || p<0)
				{
					printf("Error: trans: %i, segment: %i num_seg: %i\n", i, p, num_seg);
					continue;
				}
				trm[i*num_seg+p] = 1;
			}
		}
		mxSetFieldByNumber(plhs[0], j, transcripts_field, trans_mat);

		// segment adjacency matrix
		len = regions[j]->pair_mat.size();
		//printf("admat.size(): %i, %i\n", len, seg_admat_field);
		mxArray* pair_mat = mxCreateDoubleMatrix(len, len, mxREAL);
		double* pmat = mxGetPr(pair_mat);
		for (int i=0; i<len; i++)
		{
			for (int k=0; k<len; k++)
			{
				pmat[len*i+k] = regions[j]->pair_mat[i][k];
			}
		}
		mxSetFieldByNumber(plhs[0], j, pair_field, pair_mat);
	}

	// cleanup
	for (int j=0; j<num; j++)
		delete regions[j];
}
