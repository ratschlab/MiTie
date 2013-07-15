#include <stdio.h>
#include <mex.h>
#include <assert.h>
#include <vector>
  using std::vector;
#include "region.h"
#include "get_var.h"
#include <algorithm>
	using std::max;
#include <fstream>
#include "tools.h"

/*
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs<2)
	{
		mexErrMsgTxt("Usage: load_regions_bin(fname, num_genes, fn_bam1, fn_bam2, ...)");	
	}

	char* fname = get_string(prhs[0]);
	
	int gene_id = get_int(prhs[1]);

	vector<char*> bam_files;
	for (int i=2; i<nrhs; i++)
		bam_files.push_back(get_string(prhs[i]));
		
	//printf("load_regions_bin: found %i samples\n", (int)bam_files.size());
	vector<vector<char*> > samples;
	for (int i=0; i<bam_files.size(); i++)
	{
		//printf("Sample %i:\n", i+1);
		vector<char*> bams = separate(bam_files[i], ',');
		//for (int j = 0; j<bams.size(); j++)
		//	printf("\t%s\n", bams[j]);

		samples.push_back(bams);
	}
	FILE* fd_null = fopen("/dev/null", "w");

	vector<Bam_Region*> regions; 
	int cnt = 0;
	{
		Bam_Region* r = new Bam_Region();
		int ret = r->read_HDF5(fname, gene_id);
		if (ret==0)
			regions.push_back(r);
		else
		{
			printf("could not read gene %i from file: %s\n", gene_id, fname);
			delete r;
			exit(-1);
		}
	}
	//r->write_segment_graph(stdout);

	const char *field_names[] = {"chr", "strand", "start", "stop", "segments", "coverage", "seg_admat", "transcripts", "transcript_names", "pair_mat"};
	int nof_fields = 10;

	int num = regions.size();
	//printf("\nparsed %i regions from file\n", num);

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
	int transcript_names_field = mxGetFieldNumber(plhs[0], "transcript_names");
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
		mxArray* coverage = mxCreateDoubleMatrix(samples.size(), num_seg, mxREAL);
		double* cov = mxGetPr(coverage);
		mxSetFieldByNumber(plhs[0], j, coverage_field, coverage);

		// segment adjacency matrix
		int seg_admat_size = regions[j]->admat.size();
		//printf("admat.size(): %i, %i\n", len, seg_admat_field);
		//mxArray* seg_admat = mxCreateDoubleMatrix(len, len, mxREAL);
		int dims[3];
		dims[0] = seg_admat_size;
		dims[1] = seg_admat_size;
		dims[2] = samples.size();
		mxArray* seg_admat = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		double* seg_mat = mxGetPr(seg_admat);
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
				trm[p*num_trans+i] = 1;
			}
		}
		mxSetFieldByNumber(plhs[0], j, transcripts_field, trans_mat);
		
		//transcript names
		mxArray* names_cell = mxCreateCellArray(1, &num_trans);
		for (int i=0; i<regions[j]->transcript_names.size(); i++)
		{
			const char* name = regions[j]->transcript_names[i].c_str();
			mxArray* sname = mxCreateString(name);
			mxSetCell(names_cell, i, sname);
		}
		mxSetFieldByNumber(plhs[0], j, transcript_names_field, names_cell);

		// segment adjacency matrix
		int len = regions[j]->pair_mat.size();
		dims[0] = len;
		dims[1] = len;
		dims[2] = samples.size();
		mxArray* pair_mat = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		//printf("admat.size(): %i, %i\n", len, seg_admat_field);
		//mxArray* pair_mat = mxCreateDoubleMatrix(len, len, mxREAL);
		mxSetFieldByNumber(plhs[0], j, pair_field, pair_mat);


		for (int k = 0; k<samples.size(); k++)
		{
			if (bam_files.size()>k)
			{
				// each string in bam_files may be a comma separated list of files
				vector<char*> bams = samples[k];

				int intron_len_filter = 200000;
				int filter_mismatch = 10;
				int exon_len_filter = 0;
				bool mm_filter = false;
				//printf("chr: %s\n", regions[j]->chr);
				regions[j]->fd_out = fd_null;
				regions[j]->clear_reads();
				regions[j]->get_reads(&bams[0], bams.size(), intron_len_filter, filter_mismatch, exon_len_filter, mm_filter);
				//printf("admat.size():%i  segments.size():%i\n", (int) regions[j]->segments.size(), (int) regions[j]->admat.size());
				regions[j]->update_coverage_information();
				//printf("num_reads: %i\n", (int) regions[j]->reads.size());
				regions[j]->compute_coverage();
				regions[j]->compute_seg_cov();
				regions[j]->compute_pair_mat();
				//printf("admat.size():%i  segments.size():%i\n", (int) regions[j]->segments.size(), (int) regions[j]->admat.size());
			}

			for (int i=0; i<num_seg; i++)
				cov[i*samples.size()+k] = regions[j]->seg_cov[i];

			for (int i=0; i<seg_admat_size; i++)
			{
				for (int l=0; l<seg_admat_size; l++)
				{
					seg_mat[k*seg_admat_size*seg_admat_size+seg_admat_size*i+l] = regions[j]->admat[i][l];
				}
			}

			double* pmat = mxGetPr(pair_mat);
			for (int i=0; i<len; i++)
			{
				for (int l=0; l<len; l++)
				{
					pmat[k*len*len+len*i+l] = regions[j]->pair_mat[i][l];
				}
			}

		}
	}


	// cleanup
	for (int j=0; j<num; j++)
		delete regions[j];
	
	fclose(fd_null);
}
