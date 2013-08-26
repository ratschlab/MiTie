#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <fstream>
#include "bam.h"
#include "region.h"
#include "bam_region.h"
#include "get_reads_direct.h"
#include "gtf_tools.h"
#include "tools.h"
#include "file_stats.h"
#include "graph_tools.h"
#include "basic_tools.h"
#include "QP.h"
#include "vector_op.h"
#include "create_loss_param.h"
#include <time.h>
#include "transcript_prediction.h"
#include "config.h"

#ifdef USE_CPLEX
#include "solve_qp_cplex.h"
#endif

#ifdef USE_GLPK
#include "solve_lp_glpk.h"
#endif

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


#define USE_HDF
//#define DEBUG
//#define DEBUG2


int parse_args(int argc, char** argv,  Config* c)
{
	if (argc<4)
	{
			fprintf(stdout, "Usage: %s <fn_graph> [options] <fn_bam1,fn_bam2> <fn_bam3> ...\n", argv[0]);
			fprintf(stdout, "(this corresponds to two samples where sample 1 has two bam files and sample two has one)\n");
			fprintf(stdout, "general options:\n");
			fprintf(stdout, "\t--fn-quant\t\t(default NULL) file name to write quantification values\n");
			fprintf(stdout, "\t--fn-out\t\t(default NULL) file name to write transcript structures\n");
			fprintf(stdout, "\t--num-trans\t\t(default 5) maximal number of transcripts predicted in addition to the annotated transcripts\n");
			fprintf(stdout, "\t\t\n");

			fprintf(stdout, "loss function:\n");
			fprintf(stdout, "\t--order\t\t(default 2) order of the polynom to fit the loss function 1 (lp) or 2 (qp)\n");
			fprintf(stdout, "\t--param-eta1\t\t(default 1.0)\n");
			fprintf(stdout, "\t--param-eta2\t\t(default 0.1)\n");
			fprintf(stdout, "\t--param-lambda\t\t(default 3)\n");
			fprintf(stdout, "\t\t\n");

			fprintf(stdout, "read related options:\n");
			fprintf(stdout, "\t--min-exonic-len\t(default 0) minimal number of aligned positions on each side of an intron\n");
			fprintf(stdout, "\t--mismatches\t\t(default 10) maximal number of mismatches\n");
			fprintf(stdout, "\t--best\t\t\t(flag) filter for best alignment\n");
			fprintf(stdout, "\t\t\n");

			fprintf(stdout, "regularization related options:\n");
			fprintf(stdout, "\t--C-exon\t\t(default 1.0) loss param exon read count\n");
			fprintf(stdout, "\t--C-intron\t\t(default 10.0) loss param intron read count\n");
			fprintf(stdout, "\t--C-pair\t\t(default 1.0) unexplained pair count penalty\n");
			fprintf(stdout, "\t--C-num-trans\t\t(default 100.0) l_0 norm penalty for newly predicted transcripts\n");
			fprintf(stdout, "\t--C-num-trans-predef\t(default 1.0) l_0 norm penalty for annotated transcripts\n");
			fprintf(stdout, "\t\t\n");
			return -1;
	}

	// defaults
	c->fn_graph = argv[1];	
	c->fn_quant = NULL;
	c->fn_gtf = NULL;
	c->intron_len_filter = 200000;
	c->filter_mismatch = 10;
	c->exon_len_filter = 0;
	c->mm_filter = false;
	c->max_num_trans = 5;
	c->min_trans_len = 100;
	c->use_pair = false;

	c->order = 2;
	c->eta1 = 1.0;
	c->eta2 = 0.1;
	c->lambda = 3;

	c->C_exon = 1.0;
	c->C_intron = 10.0;
	c->C_pair = 1.0;
	c->C_num_trans = 100.0;
	c->C_num_trans_predef = 1.0;

    for (int i = 2; i < argc; i++)  
    {
	   	if (strcmp(argv[i], "--min-exonic-len") == 0)
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
	    else if (strcmp(argv[i], "--max-num-trans") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --max-num-trans\n") ;
                return -1;
            }
            i++;
			c->max_num_trans = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--order") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --order\n") ;
                return -1;
            }
            i++;
			c->order = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--param-eta1") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --param-eta1\n") ;
                return -1;
            }
            i++;
			c->eta1 = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--param-eta2") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --param-eta2\n") ;
                return -1;
            }
            i++;
			c->eta2 = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--param-lambda") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --param-lambda\n") ;
                return -1;
            }
            i++;
			c->lambda = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--C-exon") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --C-exon\n") ;
                return -1;
            }
            i++;
			c->C_exon = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--C-intron") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --C-intron\n") ;
                return -1;
            }
            i++;
			c->C_intron = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--C-pair") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --C-pair\n") ;
                return -1;
            }
            i++;
			c->C_pair = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--C-num-trans") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --C-num-trans\n") ;
                return -1;
            }
            i++;
			c->C_num_trans = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--C-num-trans-predef") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --C-num-trans-predef\n") ;
                return -1;
            }
            i++;
			c->C_num_trans_predef = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--fn-quant") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --fn-quant\n") ;
                return -1;
            }
            i++;
			c->fn_quant = argv[i];
        }
	    else if (strcmp(argv[i], "--fn-out") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --fn-out\n") ;
                return -1;
            }
            i++;
			c->fn_gtf = argv[i];
        }
		else
		{
			c->bam_files.push_back(argv[i]);
		}	
	}
	return 0;
}
int connect_neighbors(vector<vector<vector<float> > >* all_admat, vector<segment>* segments)
{
	int num_samples = all_admat->size();
	for (int j = 1; j<segments->size()-1; j++)
	{
		if (segments->at(j).second+1 == segments->at(j+1).first)
		{
			assert(all_admat->at(0).size()>j+2);
			//printf("segment[%i], %i,%i %.2f\n", j, segments->at(j).first, segments->at(j).second, all_admat->at(0)[j+1][j+2]);
			//printf("segment[%i], %i,%i\n", j, segments->at(j).first, segments->at(j).second);
			//printf("\n");

			for (int s=0; s<num_samples; s++)
			{
				assert(all_admat->at(s).size()>j+2);
				all_admat->at(s)[j+1][j+2] = -1;
			}
		}
	}

}

int union_connections(vector<vector<vector<float> > >* all_admat)
{
	// assuming that the first and the last segment are artificial source/sink nodes
	assert(all_admat->at(0).size()>0);
	assert(all_admat->at(0).size()==all_admat->at(0)[0].size());
	int num_nodes = all_admat->at(0).size();
	int num_samples = all_admat->size();

	for (int j=1; j<num_nodes-1; j++)
	{
		for (int k=1; k<num_nodes-1; k++)
		{
			float max = -2;
			for (int s=0; s<num_samples; s++)
			{
				if (all_admat->at(s)[j][k]>max)
					max = all_admat->at(s)[j][k];
			}
			if (max<0)
				continue;

			for (int s=0; s<num_samples; s++)
			{
				if (all_admat->at(s)[j][k]<0)
					all_admat->at(s)[j][k] = 0;
			}
		}
	}
	return 0;
}

semi_sparse_3d_matrix<float> compute_intron_list(vector<vector<vector<float> > >* all_admat)
{
	semi_sparse_3d_matrix<float> list; 
	int r = all_admat->size();
    int s = all_admat->at(0).size();;
    for (int j=1; j<s-1; j++)
	{
        for (int k=j+1; k<s-1; k++)
		{
            if (all_admat->at(0)[j][k]>=0)
            {
				// valid intron not just a neighboring segment
				for (int sample=0; sample<r; sample++)
				{
					//seg.features.push_back(all_admat->at(sample)[j][k]);
					list.add(j-1, k-1, all_admat->at(sample)[j][k]);
				}
			}
		}
	}
	return list;
}

semi_sparse_3d_matrix<float> compute_pair_list(vector<vector<vector<int> > >* all_pair_mat)
{
	semi_sparse_3d_matrix<float> list; 
	int r = all_pair_mat->size();
    int s = all_pair_mat->at(0).size();;
    for (int j=1; j<s-1; j++)
	{
        for (int k=j+1; k<s-1; k++)
		{
            if (all_pair_mat->at(0)[j][k]>=0)
            {
				for (int sample=0; sample<r; sample++)
				{
					list.add(j, k, (float) all_pair_mat->at(sample)[j][k]);
				}
			}
		}
	}
	return list;
}

void Tr_Pred::get_loss_coef(float (&ret)[4], float obs)
{
	// find bin
	int j=0;
	for (;j<loss_param->size(); j++)
	{
		if (obs<=loss_param->at(j)[0])
			break;
	}
	if (j==0 || j==loss_param->size()-1)
	{
		ret[0] = loss_param->at(j)[1]; //left_l
		ret[1] = loss_param->at(j)[2]; //left_q
		ret[2] = loss_param->at(j)[3]; //right_l
		ret[3] = loss_param->at(j)[4]; //right_q
		return ;
	}
	double left_pos = loss_param->at(j-1)[0];
	double right_pos = loss_param->at(j)[0];
	for (int i=0; i<4; i++)
	{
		// linear interpolation between parameters
		ret[i] = ((obs-left_pos)*loss_param->at(j-1)[i+1] + (right_pos-obs)*loss_param->at(j)[i+1])/(right_pos-left_pos);

		//printf("%i, %.5f, %.12f, %.12f, %.12f\n", i, obs, loss_param->at(j-1)[i+1], ret[i], loss_param->at(j)[i+1]);
		float eps = 1e-10;
		if (loss_param->at(j-1)[i+1]<loss_param->at(j)[i+1])
		{
			assert(ret[i]>=loss_param->at(j-1)[i+1]-eps);
			assert(ret[i]<=loss_param->at(j)[i+1]+eps);
		}
		else
		{
			assert(ret[i]<=loss_param->at(j-1)[i+1]+eps);
			assert(ret[i]>=loss_param->at(j)[i+1]-eps);
		}
	}
}

void Tr_Pred::make_qp()
{

	int len = 0;
	int all_len[graph->segments.size()];
	for (int i=0; i<graph->segments.size(); i++)
	{
		all_len[i] = graph->segments[i].second-graph->segments[i].first;
		len += all_len[i];
	}
	printf("length of all segments: %i\n", len);
	int min_trans_len = min(len, config->min_trans_len);

	semi_sparse_3d_matrix<float> intron_list = compute_intron_list(all_admat);
	semi_sparse_3d_matrix<float> pair_list = compute_pair_list(all_pair_mat);

	printf("intron_list.size(): %lu\n", intron_list.size());
	printf("pair_list.size(): %lu\n", pair_list.size());

	int num_introns = intron_list.size();
	int num_samples = all_admat->size();
	vector<int> binary_var_index_all;

	int num_annotated_trans = graph->transcript_paths.size();
	int max_num_trans = config->max_num_trans + graph->transcript_paths.size();




	// compute the number of variables
	//
	// U_st : usage of segment s in transcript t                                     (binary)
	// I_t  : indicator if transcript t has weight>=0 in any sample                  (binary)
	// E_str : expected coverage of segment s in transcript t in sample r                 
	// W_tr  : weight of transcript t in sample r
	// L_sr  : loss of segment s (deviation from the observed coverage) in sample r
	// C_ctr : Spliced read expected coverage (c is a subset of sxs) in sample r
	// D_cr  : deviation of expected intron coverage from observed intron cov in sample r
	//
	// var = [ U_st I_t E_str W_tr L_sr C_ctr D_cr]
	int c = num_introns;
	int t = max_num_trans;
	int s = graph->segments.size();
	int r = num_samples;
	//int nc = num_clusters;
	printf("#introns: %i, max #transcripts: %i(%lu+%i), #segments: %i, #samples: %i\n", c, t, graph->transcript_paths.size(), t-graph->transcript_paths.size(),  s, r);
	

	// compute the number of variables of the optimization problem
	int num_var = 0;
	// indices of variables in the parameter vector
	vector<int> U_idx = range(0, s*t-1); 					num_var += U_idx.size();
	vector<int> I_idx = range(num_var, num_var+t-1);        num_var += I_idx.size();
	vector<int> E_idx = range(num_var, num_var+s*t*r-1);    num_var += E_idx.size();
	vector<int> W_idx = range(num_var, num_var+t*r-1);      num_var += W_idx.size();
	vector<int> L_idx = range(num_var, num_var+s*r*2-1);    num_var += L_idx.size();
	vector<int> C_idx = range(num_var, num_var+c*t*r-1);    num_var += C_idx.size();
	vector<int> D_idx = range(num_var, num_var+c*r*2-1);    num_var += D_idx.size();
	
#ifdef DEBUG2
	printf("U_idx: %i-%i\n", U_idx[0], U_idx[U_idx.size()-1]);
	printf("I_idx: %i-%i\n", I_idx[0], I_idx[I_idx.size()-1]);
	printf("E_idx: %i-%i\n", E_idx[0], E_idx[E_idx.size()-1]);
	printf("W_idx: %i-%i\n", W_idx[0], W_idx[W_idx.size()-1]);
	printf("L_idx: %i-%i\n", L_idx[0], L_idx[L_idx.size()-1]);
	printf("C_idx: %i-%i\n", C_idx[0], C_idx[C_idx.size()-1]);
	printf("D_idx: %i-%i\n", D_idx[0], D_idx[D_idx.size()-1]);
#endif
	//if use_cluster
	//    M_idx = num_var+1:num_var+r*nc; num_var = num_var+length(M_idx); % cluster assignment
	//    m_idx = num_var+1:num_var+t*nc; num_var = num_var+length(m_idx); % cluster centroids
	//    K_idx = num_var+1:num_var+nc; num_var = num_var+length(K_idx);   % cluster usage
	//end
	if (config->use_pair)
	{
		//TODO implement
		printf("pair feature not yet implemented in cpp\n");
		return;
	    int np = pair_list.size();
	    vector<int>  P_idx = range(num_var, num_var+t*np-1); num_var += P_idx.size();
	    vector<int> SP_idx = range(num_var, num_var+np-1);   num_var += SP_idx.size(); // slacks for pair observation
		printf("P_idx: %i-%i\n", P_idx[0], P_idx[P_idx.size()-1]);
		printf("SP_idx: %i-%i\n", SP_idx[0], SP_idx[SP_idx.size()-1]);

	}

	printf("QP has %i variables\n", num_var);

	// compute scale factor for all counts
	float cov_scale[r];
	for (int i=0; i<num_samples; i++)
	{
		float max_val = max<float>(all_seg_cov->at(i));
		for (int j=0; j<all_admat->at(i).size(); j++)
		{
			float max_row = max<float>(&(all_admat->at(i)[j]));
			float min_row = min<float>(&(all_admat->at(i)[j]));
			assert(min_row>=-2);
			//printf("max_val:%.2f, max_row:%.2f\n", max_val, max_row);
			if (max_row>max_val)
				max_val = max_row;
		}
		cov_scale[i] = max_val;

		if (max_val<=0)
			continue;

		mult(all_seg_cov->at(i), 1/max_val);

		//for (int j=0; j<all_admat->at(i).size(); j++)
		//{
		//	mult(&all_admat->at(i)[j], 1/max_val);
		//}
	}

	QP* qp = new QP(num_var);

	// set lower and upper bounds
	///////////////////////////////////////////////////////////////
	qp->lb = vector<float>(num_var, -1e20);
	qp->ub = vector<float>(num_var, 1e20);
	qp->binary_idx = vector<int>(num_var, 0);
	for (int i=0; i<U_idx.size(); i++)
	{
		qp->lb[U_idx[i]] = 0;
		qp->ub[U_idx[i]] = 1;
		qp->binary_idx[U_idx[i]] = 1;
	}
	for (int i=0; i<E_idx.size(); i++)
	{
		qp->lb[E_idx[i]] = 0;
		qp->ub[E_idx[i]] = 1;
	}
	for (int i=0; i<W_idx.size(); i++)
	{
		qp->lb[W_idx[i]] = 0;
		qp->ub[W_idx[i]] = 1;
	}
	for (int i=0; i<C_idx.size(); i++)
	{
		qp->lb[C_idx[i]] = 0;
		qp->ub[C_idx[i]] = 1;
	}
	for (int i=0; i<L_idx.size(); i++)
	{
		qp->lb[L_idx[i]] = 0;
	}
	for (int i=0; i<D_idx.size(); i++)
	{
		qp->lb[D_idx[i]] = 0;
	}

	// fix solution for known transcripts 
	//for (int i=0; i<num_annotated_trans; i++)
	int debug_trans = 0;
	for (int i=0; i<num_annotated_trans; i++)
	{
		int k = 0;
		for (int j=0; j<graph->transcript_paths[i].size(); j++)
		{
			//transcript path is zero based
			while (k<graph->transcript_paths[i][j])
			{
				qp->ub[U_idx[k*t+i]] = 0;
#ifdef DEBUG
				printf("k%i:U[%i]==0\n", k, k*t+i);
#endif
				k++;
			}
			assert(k==graph->transcript_paths[i][j]);
			qp->lb[U_idx[k*t+i]] = 1;
#ifdef DEBUG
			printf("k%i:U[%i]==1\n", k, k*t+i);
#endif
			k++;
		}
		while (k<s)
		{
			qp->ub[U_idx[k*t+i]] = 0;
#ifdef DEBUG
			printf("k%i:U[%i]==0\n", k, k*t+i);
#endif
			k++;
		}
	}

	// create objective function:
	// min_x x'Qx + f'x
	///////////////////////////////////////////////////////////////
	assert(all_seg_cov->size()==r);

	qp->F = vector<float>(num_var, 0.0);
	
	// l_0 norm on transcript quantification
	for (int i=0; i<num_annotated_trans; i++)
	{
		qp->F[I_idx[i]] = config->C_num_trans_predef;
	}
	for (int i=num_annotated_trans; i<t; i++)
	{
		qp->F[I_idx[i]] = config->C_num_trans;
	}

	// set objective coefficients for exon counts
	for (int i=0; i<r; i++)
	{
		for (int j=0; j<s; j++)
		{
			int readlen = 100;
			int x1 = L_idx[i*s+j]; // loss right hand side
			int x2 = L_idx[i*s+j+s*r]; // loss left hand side
			float val = all_seg_cov->at(i)->at(j) * cov_scale[i];

			//get linear and quadratic coefficients of the loss function
			// [ll lq rl rq]
			float coef[4];
			get_loss_coef(coef, val);

#ifdef DEBUG
			printf("val[%i][%i]:%.4f (left)coef[1]:%.10f, (right)coef[3]:%.10f\n", i, j, val, coef[1], coef[3]);
			printf("reg_parts: %.10f %.10f %.10f %.10f %.10f\n", 1.0/r, coef[3], config->C_exon, (float) all_len[j]/readlen, cov_scale[i]); 
			printf("Q[%i][%i]: %.10f\n", x1, x1, 1.0/r*coef[3]*config->C_exon*pow(((float) all_len[j])/readlen*cov_scale[i], 2)); 
			printf("Q[%i][%i]: %.10f\n", x2, x2, 1.0/r*coef[1]*config->C_exon*pow(((float) all_len[j])/readlen*cov_scale[i], 2));
#endif
			qp->Q.set(x1,x1, 1.0/r*coef[3]*config->C_exon*pow(((float) all_len[j])/readlen*cov_scale[i], 2));
			qp->Q.set(x2,x2, 1.0/r*coef[1]*config->C_exon*pow(((float) all_len[j])/readlen*cov_scale[i], 2));
			qp->F[x1] = 1.0/r*coef[2]*config->C_exon*all_len[j]/readlen*cov_scale[i];
			qp->F[x2] = 1.0/r*coef[0]*config->C_exon*all_len[j]/readlen*cov_scale[i];
		}
	}
	// set objective coefficients for intron counts
	for (int i=0; i<r; i++)
	{
		intron_list.reset_it();
		int j=0;
		while (c>0)//dont do this if there are no introns
		{
			int x1 = D_idx[i*c+j];
			int x2 = D_idx[i*c+j+c*r];
			int tmp1 = 0;
			int tmp2 = 0;
			vector<float>* conf = intron_list.next(&tmp1, &tmp2);

			if (tmp1==-1)
				break;

			float val = conf->at(i);
			//printf("val[%i][%i]:%.4f\n", i, j, val);
			float coef[4];
			get_loss_coef(coef, val);

			//get linear and quadratic coefficients of the loss function
			qp->Q.set(x1,x1, 1.0/r*coef[3]*config->C_intron*pow(cov_scale[i], 2));
			qp->Q.set(x2,x2, 1.0/r*coef[1]*config->C_intron*pow(cov_scale[i], 2));
			qp->F[x1] = 1.0/r*coef[2]*config->C_intron*cov_scale[i];
			qp->F[x2] = 1.0/r*coef[0]*config->C_intron*cov_scale[i];
			j++;
		}
	}

	// create constraints
	// Ax <= b
	// for constraints with their index in eq_idx:
	// Ax = b

	bool A1 = true;
	bool A2 = true;
	bool A3 = true;
	bool A4 = true;
	bool A5 = true;
	bool A6 = true;
	bool A8 = true;
	bool A9 = true;
	bool A10 = true;
	bool A11 = true;
	bool A12 = true;
	bool A13 = true;
	bool A14 = true;
	bool A15 = true;
	bool A18 = true;
	// loss
	// L_sr1: loss right hand side
	// L_sr2: loss left hand side
	// sum_t E_str -L_sr1 + L_sr2= O_sr
	int cc = 0; // constraint count
	if (A1)
	{
		for (int i=0; i<r; i++)// loop over samples
		{
			for (int j=0; j<s; j++)
			{
				for (int k=0; k<t; k++)
				{
					int idx = i*s*t + j*t + k;
					qp->A.set(cc, E_idx[idx], 1); // this is where we included Reginas profiles 
					//if (i==0)
					//	printf("E[%i]+ ", idx); 
				}
				//if (i==0)
				//	printf("L[%i] - L[%i] = %.2f\n", i*s+j+s*r, i*s+j, all_seg_cov->at(i)->at(j));
				qp->A.set(cc, L_idx[i*s+j], -1); 
				qp->A.set(cc, L_idx[i*s+j+s*r], 1); 
				qp->b.push_back(all_seg_cov->at(i)->at(j));
				qp->eq_idx.push_back(1);
				cc++;
			}
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());
#ifdef DEBUG
	printf("A1: %i\n", cc);
#endif

	// if no segments are selected, W_t has to be 0
	// W_x<=\sum_j U_{jx}
	if (A2)
	{
		for (int i=0; i<r; i++)
		{
			for (int x=0; x<t; x++)
			{
				qp->A.set(cc, W_idx[i*t+x], 1); // W_xi
				for (int j=0; j<s; j++)
				{
					qp->A.set(cc, U_idx[j*t+x], -1); // -U_jx
				}
			}
			cc++;
			qp->b.push_back(1);
			qp->eq_idx.push_back(0);
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());
#ifdef DEBUG
	printf("A2: %i\n", cc);
#endif


	// this takes the connectivity of the splice graph and makes 
	// sure that segment j is followed by any of the segments 
	// it is connected to in the splice graph G=(N,E)
	//
	// (1) U_{jx} <= \sum_{k \in \{k | (j,k)\in E\}} U_{cx} 
	// (2) U_{jx} <= \sum_{k \in \{k | (k,j)\in E\}} U_{cx} 
	//
	// for initial (and terminal) segments make sure that there are no
	// downstream (upstream) segments used
	//
	// terminal:
	// (3) \sum_{k=j+1}^s U_kx <= (s-j) - (s-j)*U_jx
	// <=> \sum_{k=j+1}^s U_kx + (s-j)*U_jx <= (s-j)
	//
	// initial:
	// (4) \sum_{k=1}^{j-1} U_jx <= (j-1) - (j-1)*U_jx
	// <=> \sum_{k=1}^{j-1} U_jx + (j-1)*U_jx <= (j-1)

	if (A3)
	{
#ifdef DEBUG2
		for (int x=0; x<t; x++)
#else
		for (int x=num_annotated_trans; x<t; x++)
#endif
		{
			for (int j=0; j<s; j++)
			{
				// get all connected nodes
				vector<int> children; 
				for (int i=j+1; i<s; i++)// state 0 is the source and state s+1 is the sink state
				{
					if (all_admat->at(0)[j+1][i+1]>=-1)
						children.push_back(i); // zero based
				}
				vector<int> parents;
				for (int i=0; i<j; i++)
				{
					if (all_admat->at(0)[i+1][j+1]>=-1)
						parents.push_back(i);// zero based
				}
				// constraints for children
				if (!children.empty() && !graph->is_terminal(j+1))
				{

					qp->A.set(cc, U_idx[j*t+x], 1); 
#ifdef DEBUG
					if (x==debug_trans)
						 printf("U[%i] ", j*t+x); 
#endif
					for (int k=0; k<children.size(); k++)
					{
						qp->A.set(cc, U_idx[(children[k])*t+x], -1); 
#ifdef DEBUG
						if (x==debug_trans)
							printf("-U[%i] ", (children[k])*t+x);
#endif
					}
#ifdef DEBUG
					if (x==debug_trans)
						printf(" <= 0\n");
#endif
					cc++;
					qp->b.push_back(0);
					qp->eq_idx.push_back(0);
				}
				else if (false && children.empty())
				{
					//make sure there is no downstream segment used if 
					//U_jx is used
					qp->A.set(cc, U_idx.at(j*t+x), s-j); // (s-j)U_jx
					for (int k=j; k<s; k++)
						qp->A.set(cc, U_idx.at(k*t+x), 1); // U_kx
					cc++;
					qp->b.push_back(s-j);
					qp->eq_idx.push_back(0);
				}

				// constraints for parents
				if (!parents.empty() && ! graph->is_initial(j+1))
				{
					qp->A.set(cc, U_idx[j*t+x], 1); // U_jx
#ifdef DEBUG
					if (x==debug_trans)
						printf("U[%i] ", j*t+x);
#endif
					for (int k=0; k<parents.size(); k++)
					{
						qp->A.set(cc, U_idx.at((parents[k])*t+x), -1); // -U_kx

#ifdef DEBUG
						if (x==debug_trans)
							printf("-U[%i] ", (parents[k])*t+x);
#endif
					}
#ifdef DEBUG
					if (x==debug_trans)
						printf(" <= 0\n");
#endif
					cc++;
					qp->b.push_back(0);
					qp->eq_idx.push_back(0);
				}
				else if (parents.empty() && j>0)
				{
					// make sure there are no upsteam segments if 
					// U_jx is used
					qp->A.set(cc, U_idx[j*t+x], j); // U_jx
#ifdef DEBUG
					if (x==debug_trans)
						printf("%i*U[%i]", j, j*t+x);
#endif
					for (int k=0; k<j; k++)
					{
						qp->A.set(cc, U_idx[k*t+x], 1); // U_kx
#ifdef DEBUG
						if (x==debug_trans)
							printf(" + U[%i]", k*t+x);
#endif
					}
					cc++;
#ifdef DEBUG
					if (x==debug_trans)
						printf(" <= %i\n", j);
#endif
					qp->b.push_back(j);
					qp->eq_idx.push_back(0);
				}
			}
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());
#ifdef DEBUG
	printf("A3: %i\n", cc);
#endif

	// E_str
	// define helper variables computing the expected coverage for each segment
	// The expected coverage is given by 
	// E_str = U_st * W_tr
	// bound them in [0, 1]

	// set E_str == W_tr if 
	// segment s is used in transcript t
	for (int i=0; i<r; i++)
	{
		for (int j=0; j<s; j++)
		{
			for (int k=0; k<t; k++)
			{
				int pos = i*s*t+j*t+k;
				if (A4)
				{
					// E_str - U_st <= 0
					qp->A.set(cc, E_idx[pos], 1);
					qp->A.set(cc, U_idx[j*t+k], -1);
					cc++;
					qp->b.push_back(0);
					qp->eq_idx.push_back(0);
				}
				if (A5)
				{
					// E_str + U_st - W_tr <= 1
					qp->A.set(cc, E_idx[pos], 1);
					qp->A.set(cc, U_idx[j*t+k], 1);
					qp->A.set(cc, W_idx[i*t+k], -1);
					cc++;
					qp->b.push_back(1);
					qp->eq_idx.push_back(0);
				}
				if (A6)
				{
					// -E_str + U_st + W_tr <= 1
					qp->A.set(cc, E_idx[pos], -1);
					qp->A.set(cc, U_idx[j*t+k], 1);
					qp->A.set(cc, W_idx[i*t+k], 1);
					cc++;
					qp->b.push_back(1);
					qp->eq_idx.push_back(0);
				}
			}
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());
#ifdef DEBUG
	printf("A4 A5 A6: %i\n", cc);
#endif


	// I_t
	// indicator for transcripts with weight >0 in any sample
	// sum_r W_tr - r*I_t <=0 
	// -I_1 <= -1 predict at least one transcript
	if (A9)
	{
		for (int x=0; x<t; x++)
		{
			for (int i=0; i<r; i++)
			{
				qp->A.set(cc, W_idx[i*t+x], 1); // W_tr
			}
			qp->A.set(cc, I_idx[x], -r); // -r*I_t
			cc++;
			qp->b.push_back(0);
			qp->eq_idx.push_back(0);
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());
#ifdef DEBUG
	printf("A9: %i\n", cc);
#endif


	//
	// S_ss
	// spliced read penalty
	//
	// C_jkt = W_t*U_jt*prod_i=j+1^k-1{1-U_it}*U_kt
	//
	// (1) C_jkt<=U_jt
	// (2) C_jkt<=U_kt
	// (3) C_jkt<=1-U_it for all j<i<k
	// (4) C_jkt<=W_t-U_jt+1-U_kt+1+ (sum_i U_it) 
	// (5) C_jkt>= W_t+U_jt-1+U_kt-1- (sum_i U_it)
	//
	// compute expected intron coverage - observed intron coverage
	// (6) D_jk = sum_t C_jkt - O_jk 
	//
	// for introns (j,k) not in splicegraph
	// do not allow the usage of this intron:
	// (7) U_{jt}+U_{kt} <= \sum_{i=j+1}^{k-1} U_{it} +1 
	//


	// (1) C_jkt-U_jt<=0				-> A10
	// (2) C_jkt-U_kt<=0 				-> A11
	// (3) C_jkt+U_it<=1 for all j<i<k	-> A12

	for (int i=0; i<r; i++)
	{
		for (int x=0; x<t; x++)
		{
			intron_list.reset_it();
			int xx=0;
			while (true)
			{
				int j = 0;
				int k = 0;
				intron_list.next(&j, &k);

				if (j==-1)
					break;

				assert(j<s);
				assert(k<s);

				int cnt = i*t*c+x*c+xx ;  
				if (A10)
				{
					qp->A.set(cc, U_idx[j*t+x], -1);
					qp->A.set(cc, C_idx[cnt], 1);
					cc++;
					qp->b.push_back(0);
					qp->eq_idx.push_back(0);
				}
				if (A11)
				{
					qp->A.set(cc, U_idx[k*t+x], -1);
					qp->A.set(cc, C_idx[cnt], 1);
					cc++;
					qp->b.push_back(0);
					qp->eq_idx.push_back(0);
				}
				if (A12)
				{
					for (int l=j+1; l<k; l++)
					{
						qp->A.set(cc, U_idx[l*t+x], 1);
						qp->A.set(cc, C_idx[cnt], 1);
						cc++;
						qp->b.push_back(1);
						qp->eq_idx.push_back(0);
					}
				}
				xx++;
			}
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());
#ifdef DEBUG
	printf("A10 A11 A12: %i\n", cc);
#endif


	// A13 
	// (4) C_jkt-W_t+U_jt+U_kt-(sum_i U_it) <= 2
	if (A13)
	{
		for (int i=0; i<r; i++)
		{
			for (int x=0; x<t; x++)
			{
				intron_list.reset_it();
				int xx=0;
				while (true)
				{
					int j = 0;
					int k = 0;
					intron_list.next(&j, &k);

					if (j==-1)
						break;

					int cnt = i*t*c+x*c+xx;  
					qp->A.set(cc, U_idx[j*t+x], 1);
					qp->A.set(cc, U_idx[k*t+x], 1);
					qp->A.set(cc, C_idx[cnt], 1);
					qp->A.set(cc, W_idx[i*t+x], -1);
					//if (i==0 && x==0)
					//	printf("C[%i] - W[%i] + U[%i] + U[%i]", cnt, i*t+x, k*t+x, j*t+x);
					for (int l=j+1; l<k; l++)
					{
						qp->A.set(cc, U_idx[l*t+x], -1);
					//	if (i==0 && x==0)
					//		printf("-U[%i] ", l*t+x);
					}
					//if (i==0 && x==0)
					//	printf(" <= 2\n");
					cc++;
					qp->b.push_back(2);
					qp->eq_idx.push_back(0);
					xx++;
				}
			}
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());
#ifdef DEBUG
	printf("A13: %i\n", cc);
#endif

	// (5) -C_jkt+W_t+U_jt+U_kt-(sum_i U_it) <= 2
	if (A14)
	{
		for (int i=0; i<r; i++)
		{
			for (int x=0; x<t; x++)
			{
				intron_list.reset_it();
				int xx=0;
				while (true)
				{
					int j = 0;
					int k = 0;
					intron_list.next(&j, &k);

					if (j==-1)
						break;

					int cnt = i*t*c+x*c+xx;  
					qp->A.set(cc, U_idx[j*t+x], 1);
					qp->A.set(cc, U_idx[k*t+x], 1);
					qp->A.set(cc, C_idx[cnt], -1);
					qp->A.set(cc, W_idx[i*t+x], 1);
					for (int l=j+1; l<k; l++)
					{
						qp->A.set(cc, U_idx[l*t+x], -1);
					}
					cc++;
					qp->b.push_back(2);
					qp->eq_idx.push_back(0);
					xx++;
				}
			}
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());


	// (6) D_jk1 + D_jk2 = sum_t C_jkt - O_jk
	// (6) sum_t C_jkt - D_jk1 + D_jk2  =  O_jk
	if (A15)
	{
		for (int i=0; i<r; i++)
		{
			intron_list.reset_it();
			int xx=0;
			while (c>0)
			{
				int j = 0;
				int k = 0;
				vector<float>* conf = intron_list.next(&j, &k);

				if (j==-1)
					break;

				assert(conf->size()==r);
				for (int x=0; x<t; x++)
				{
					int cnt = i*t*c+x*c+xx;
					qp->A.set(cc, C_idx[cnt], 1); 
				}
				qp->A.set(cc, D_idx[i*c+xx], -1); 
				qp->A.set(cc, D_idx[i*c+xx+c*r], 1); 
		
				cc++;
				qp->b.push_back(conf->at(i)/cov_scale[i]);
				qp->eq_idx.push_back(1);
				xx++;
			}
		}
	}

	// for introns (j,k) not in splicegraph
	// do not allow the usage of this intron:
	// (7) U_{jt}+U_{kt} <= \sum_{i=j+1}^{k-1} U_{it} +1 
	// (7) U_{jt}+U_{kt}-\sum_{i=j+1}^{k-1} U_{it} <= 1 
	//
	// negative formulation
	// all not known introns are forbidden
	if (A18)
	{
#ifdef DEBUG2
		for (int x=0; x<t; x++)
#else
		for (int x=num_annotated_trans; x<t; x++)
#endif
		{
			for (int j=0; j<s; j++)
			{
				 // find the last segment that is connected to segment j
				 // constraint matrix A3 will make sure that one of them is 
				 // used if the segment is not terminal
				 // => if the segment is not terminal introns larger than to 
				 // the last connected segment are excluded any way and 
				 // we do not need to do this here
				bool no_neighbors = false;
				vector<int> children; 
				for (int i=j+1; i<s; i++)
				{
					if (all_admat->at(0)[j+1][i+1]>=-1)
						children.push_back(i);// zero based
				}

				int maxs = 0;
				if (graph->is_terminal(j+1))
				{
					maxs = s;
				}
				else if (children.size()==0)
				{
					// this case is handled by A3
					// it makes sure that no downstream segment is used
					// thus there is also no intron
					continue;
				}
				else
				{
					maxs = max<int>(&children);
				}
				assert(maxs<=s);

				for (int k=j+1; k<maxs; k++)
				{
					// check if j->k is a valid intron or neighboring segment
					if (all_admat->at(0)[j+1][k+1]>=-1)
						continue;
				
					qp->A.set(cc, U_idx[j*t+x], 1); 
					qp->A.set(cc, U_idx[k*t+x], 1); 
#ifdef DEBUG
					if (x==debug_trans)
						printf("U[%i] + U[%i] ", j*t+x, k*t+x);
#endif
					for (int l=j+1; l<k; l++)
					{
						qp->A.set(cc, U_idx[l*t+x], -1);
#ifdef DEBUG
						if (x==debug_trans)
							printf("-U[%i] ", l*t+x);
#endif
					}
					cc++;
					qp->b.push_back(1);
					qp->eq_idx.push_back(0);
#ifdef DEBUG
					if (x==debug_trans)
						printf(" <= 1\n");
#endif
				}
			}
		}
	}
	assert(cc==qp->b.size());
	assert(cc==qp->eq_idx.size());



	bool success = true;
#ifdef USE_CPLEX
	printf("solve qp using cplex\n");
	qp->result = solve_qp_cplex(qp, success);
//#else
#ifdef USE_GLPK
	printf("solve lp using glpk\n");
	qp->result = solve_lp_glpk(qp, success);
#else
	printf("no solver available; stopping here\n");
	exit(-1);
#endif
#endif

	assert(qp->result.size() == qp->num_var);
	printf("res[E_idx]\n");
	// sum_t E_str -L_sr = O_sr
	for (int i=0; i<r; i++)// loop over samples
	{
		for (int j=0; j<s; j++)
		{
			float exp_sr = 0.0;
			for (int k=0; k<t; k++)
			{
				int idx = i*s*t + j*t + k;
				exp_sr += qp->result[E_idx[idx]];
			}
			printf("%.2f ", exp_sr);
		}
		printf("\n");
	}

	printf("res[L_idx]\n");
	for (int i=0; i<r; i++)
	{
		printf("left\t");
		for (int j=0; j<s; j++)
		{
			printf("%.2f ", qp->result[L_idx[i*s+j+ s*r]]);
		}
		printf("\n");
		printf("right\t");
		for (int j=0; j<s; j++)
		{
			printf("%.2f ", qp->result[L_idx[i*s+j]]);
		}
		printf("\n");
	}

	printf("res[W_idx]: res[U_idx]\n");
	for (int i=0; i<t; i++)
	{
		for (int j=0; j<r; j++)
			printf(" %.2f", qp->result[W_idx[j*t+i]]);
		printf(": "); 
		for (int j=0; j<s; j++)
		{
			printf("%i ", qp->result[U_idx[j*t+i]]>0.5);
		}
		printf("\n");
	}



	// prepare output
	//
	
	for (int i=num_annotated_trans; i<t; i++)
	{
		if (qp->result[I_idx[i]]<0.5)// check if transcript is expressed (using the binary indicator variable)
			continue;

		vector<int> new_path;
		for (int j=0; j<s; j++)
		{
			if (qp->result[U_idx[j*t+i]]>0.5)
				new_path.push_back(j);
		}
		assert(new_path.size()>0);
		graph->transcript_paths.push_back(new_path);

		char name[1000];
		sprintf(name, "mitie_%i_%i", graph->id, i);
		graph->transcript_names.push_back(string(name));
		if (graph->gene_names.size()==i)
		{
			graph->gene_names.push_back(graph->gene_names[i-1]);
		}
	}

	if (config->fn_gtf)
	{
		//printf("%i transcripts\n", graph->transcripts.size());
		graph->transcripts.clear();
		graph->compute_transcripts_from_paths();
		//printf("after: %i transcripts\n", graph->transcripts.size());
		//for (int i=0; i<graph->transcripts.size(); i++)
		//	printf("trans%i: %i\n", i, graph->transcripts[i].size());
		FILE* fd_gtf = fopen(config->fn_gtf, "a");
		if (!fd_gtf)
		{
			printf("could not open file %s for writing\n", config->fn_gtf);
			exit(-1);
		}
		write_gtf(fd_gtf, graph, "mitie");
		fclose(fd_gtf);
	}
	
	// write quantification values to file
	FILE* fd_quant; 
	if (config->fn_quant && (fd_quant = fopen(config->fn_quant, "a")))
	{
		for (int i=0; i<t; i++)
		{
			if (i<graph->transcript_names.size())
			{
				fprintf(fd_quant, "%s", graph->transcript_names[i].c_str());
				for (int j=0; j<r; j++)
				{
					fprintf(fd_quant, "\t%.5f", qp->result[W_idx[j*t+i]]*cov_scale[j]); 
				}
				fprintf(fd_quant, "\n"); 
			}
		}
		fclose(fd_quant);
	}
	else if (config->fn_quant)
	{
		printf("could not open file: %s for writing\n", config->fn_quant);
		exit(-1);
	}
}

int main(int argc, char* argv[])
{
	Config c;
	int ret = parse_args(argc, argv, &c);
	if (ret!=0)
		return ret;


	time_t t = time(NULL);
	vector<vector<double> > loss_param = create_loss_parameters(c.eta1, c.eta2, c.lambda, c.order);
	printf("time diff: %lu\n", time(NULL)-t);


	printf("loading graphs from file: %s\n", c.fn_graph);
	
#ifndef USE_HDF
	static std::ifstream ifs;
	ifs.open(c.fn_graph, std::ios::binary);
#endif
	FILE* fd_null = fopen("/dev/null", "w");


	printf("found %i samples\n", (int)c.bam_files.size());
	vector<vector<char*> > samples;
	for (int i=0; i<c.bam_files.size(); i++)
	{
		printf("Sample %i:\n", i+1);
		vector<char*> bams = separate(c.bam_files[i], ',');
		for (int j = 0; j<bams.size(); j++)
			printf("\t%s\n", bams[j]);

		samples.push_back(bams);
	}

	// load graphs from binary file
	vector<Bam_Region*> graphs; 
	int cnt = 0;
	//int num_return = 1e6;// run on all
	int num_return = 3;
#ifdef USE_HDF
	while (cnt<num_return)
#else
	while (ifs.good()&& cnt<num_return)
#endif
	{
		if ((cnt++)%100==0)
			printf("\r %i", cnt);
		Bam_Region* r = new Bam_Region();
#ifdef USE_HDF		
		int ret = r->read_HDF5(c.fn_graph, cnt);
		r->id = cnt;
		if (ret==0)
#else
		int ret = r->read_binary(&ifs);
		if (!ifs.eof() && ret==0)
#endif
			graphs.push_back(r);
		else
		{
			printf("read from file %s, %i, read %lu graphs\n", c.fn_graph, ret, graphs.size());
			delete r;
			break;
		}
	}
	printf("\n");

	printf("obtained %i graphs\n", graphs.size());

	if (true)
	{
		for (uint j=0; j<graphs.size(); j++)
		{
			int num_paths = graphs[j]->compute_num_paths();
			printf("num_paths\t%i\n", num_paths);
		}
	}

	for (uint j=0; j<graphs.size(); j++)
	{
		// store coverage informaton for all samples
		vector<vector<vector<float> > > all_admat;
		vector<vector<vector<int> > > all_pair_mat;
		vector<vector<float>* > all_seg_cov;

		for (uint k = 0; k<samples.size(); k++)
		{
			// each string in bam_files may be a comma separated list of files
			vector<char*> bams = samples[k];

			int intron_len_filter = 200000;
			int filter_mismatch = 10;
			int exon_len_filter = 0;
			bool mm_filter = false;//multi mappers
			//printf("chr: %s\n", graphs[j]->chr);
			graphs[j]->fd_out = fd_null;
			graphs[j]->clear_reads();
			graphs[j]->get_reads(&bams[0], bams.size(), intron_len_filter, filter_mismatch, exon_len_filter, mm_filter);
			printf("admat.size():%i  segments.size():%i\n", (int) graphs[j]->segments.size(), (int) graphs[j]->admat.size());
			graphs[j]->update_coverage_information();
			printf("num_reads: %i\n", (int) graphs[j]->reads.size());
			graphs[j]->compute_coverage();
			graphs[j]->compute_seg_cov();
			graphs[j]->compute_pair_mat();
			printf("admat.size():%i  segments.size():%i\n", (int) graphs[j]->segments.size(), (int) graphs[j]->admat.size());
			all_admat.push_back(graphs[j]->admat);
			//print_mat(&graphs[j]->admat, "%.2f ");
			all_pair_mat.push_back(graphs[j]->pair_mat);
			all_seg_cov.push_back(new vector<float>(graphs[j]->seg_cov));
		}

		// process the connectivity matrix of the graph
		// add implicit connections between neighboring segments
		connect_neighbors(&all_admat, &graphs[j]->segments);

		// make sure connections are valid in each sample if 
		// they have evidence in one sample
		union_connections(&all_admat);

		Tr_Pred* tr_pred = new Tr_Pred();
		tr_pred->graph = graphs[j];
		tr_pred->all_admat = &all_admat;
		tr_pred->all_pair_mat = &all_pair_mat;
		tr_pred->all_seg_cov = &all_seg_cov; 
		tr_pred->config = &c;
		tr_pred->loss_param = &loss_param;
		
		tr_pred->make_qp(); 

		//printf("number of constraints: %lu\n", tr_pred->qp->b.size());

		// solve QP
		//qp->result = vector<double>(qp->b.size(), 1);
		//qp->result = solve_qp_cplex(qp);
		//qp->line_search();

		//printf("obj = %.4f\n", tr_pred->qp->compute_obj());

		for (uint k = 0; k<samples.size(); k++)
			delete all_seg_cov[k];
		delete tr_pred;
	}

	for (uint j=0; j<graphs.size(); j++)
	{
		delete graphs[j];
	}
#ifndef USE_HDF
	ifs.close();
#endif
}


