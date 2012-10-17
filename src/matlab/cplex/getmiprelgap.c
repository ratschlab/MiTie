/* return the relative gap for mixed integer problems 
	relgap = (bestinteger - bestobjective) / (1e-10 + |bestobjective|) 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"

#include "cplex.h"

#define myMalloc CPXmalloc 
#define myFree CPXfree
#define STD_OUT stderr

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
	long *lpenv=NULL, *p_lp=NULL;
	CPXENVptr     env = NULL;
	CPXLPptr      lp = NULL;
	if (nrhs > 2 || nrhs < 2) 
	{
		mexErrMsgTxt("Usage: [relgap] = getmiprelgap(env, lp)");
		return;
	}
	if (1 != mxGetM(prhs[1]))
	{
		mexErrMsgTxt("Dimension error (arg 2).");
		return;
	}
    p_lp = (long*) mxGetPr(prhs[1]);
    
	if (1 != mxGetM(prhs[0])) 
	{
		mexErrMsgTxt("Dimension error (arg 1).");
		return;
	}
	lpenv = (long*) mxGetPr(prhs[0]);
  
	if (nlhs > 1 || nlhs < 1) 
	{
		mexErrMsgTxt("Usage: [relgap] = getmiprelgap(env, lp)");
		return;
	}

	/* Initialize the CPLEX environment */
	env = (CPXENVptr) lpenv[0] ;
	lp=(CPXLPptr)p_lp[0] ;

	int status;
	double gap=0.0;
	double obj_best=0.0;
	double obj_mip=0.0;
	/*int status = CPXgetmiprelgap (env, lp, &gap); this only exists in cplex11 or later*/
	/* get best relaxed solution*/
	status = CPXgetbestobjval(env, lp, &obj_best);
	if ( status ) 
	{
    	fprintf (STD_OUT, "CPXgetbestobjval failed.\n");
		char  errmsg[1024];
		CPXgeterrorstring (env, status, errmsg);
		fprintf (STD_OUT, "%s", errmsg);
  	}

	/* get best integer*/
	status = CPXgetmipobjval(env, lp, &obj_mip);
	if ( status ) 
	{
    	fprintf (STD_OUT, "CPXgetmipobjval failed.\n");
		char  errmsg[1024];
		CPXgeterrorstring (env, status, errmsg);
		fprintf (STD_OUT, "%s", errmsg);
  	}
	gap = (obj_mip - obj_best) / (1e-10 + obj_best); 

	fprintf(stdout, "Best: %.4f, MIP_obj: %.4f, gap:%.4f", obj_best, obj_mip, gap);

	/* return relative gap */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* ret = (double*) mxGetData(plhs[0]);
	ret[0] = gap;

	return ;
} 
