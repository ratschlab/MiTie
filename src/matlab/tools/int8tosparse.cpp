/*
 * =============================================================
 * int8tosparse.c
 * This example demonstrates how to populate a sparse
 * matrix. For the purpose of this example, you must pass in a
 * non-sparse 2-dimensional argument of type int8.
 *
 * Mathworks fulltosparse.c modified by J. Tursa for int8.
 *
=============================================================
 */

#include <stdio.h>
#include <math.h> /* Needed for the ceil() prototype. */
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Declare variables. */
	mwSize *irs,*jcs;
	double *sr;
	signed char *pr;
	
	/* Check for proper number of input and output arguments. */ 
	if (nrhs != 1) {
		mexErrMsgTxt("One input argument required.");
	} 
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}
	
	/* Check data type of input argument. */
	if (!(mxIsInt8(prhs[0]))) {
		mexErrMsgTxt("Input argument must be of type int8.");
	} 
	
	if (mxGetNumberOfDimensions(prhs[0]) != 2) {
		mexErrMsgTxt("Input argument must be two dimensional\n");
	}
	
	/* Get the size and pointers to input data. */
	int m = mxGetM(prhs[0]);
	int n = mxGetN(prhs[0]);
	pr = (signed char*) mxGetPr(prhs[0]);
	
	/* Allocate space for sparse matrix. */
	int nzmax = 1;
	signed char* pr2 = pr;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<m; j++)
		{
			if (pr2[j])
				nzmax++;
		}
		pr2+=m;
	}

	plhs[0] = mxCreateSparse(m,n,nzmax,mxREAL);
	sr = mxGetPr(plhs[0]);
	irs = mxGetIr(plhs[0]);
	jcs = mxGetJc(plhs[0]);
    
	/* Copy nonzeros. */
	int k = 0; 
	for (int j=0; j<n; j++)
	{
		jcs[j] = k;
		for (int i=0; i<m; i++)
		{
			if (pr[i])
			{
				sr[k] = pr[i];
				irs[k] = i;
				k++;
			}
		}
		pr += m;
	}
	jcs[n] = k;
}
