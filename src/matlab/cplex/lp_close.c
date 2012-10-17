/* lpex3.c, example of using CPXaddrows to solve a problem */

#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

#include "cplex.h"

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
  long *lpenv, *p_lp;
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status ;
  
  if (nlhs !=0 || nrhs != 2) {
    mexErrMsgTxt("Usage: lp_close(lpenv, p_lp)");
    return;
  }
  
  if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
      ||  mxIsSparse(prhs[1])
      || !mxIsDouble(prhs[1]) 
      ||  mxGetN(prhs[1])!=1 ) {
    mexErrMsgTxt("2nd argument (p_lp) must be "
		 "a column vector.");
    return;
  }
  if (1 != mxGetM(prhs[1])) {
    mexErrMsgTxt("Dimension error (arg 2).");
    return;
  }
  p_lp = (long*) mxGetPr(prhs[1]);

  if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
      ||  mxIsSparse(prhs[0])
      || !mxIsDouble(prhs[0]) 
      ||  mxGetN(prhs[0])!=1 ) {
    mexErrMsgTxt("1st argument (lpenv) must be "
		 "a column vector.");
    return;
  }
  if (1 != mxGetM(prhs[0])) {
    mexErrMsgTxt("Dimension error (arg 1).");
    return;
  }
  lpenv = (long*) mxGetPr(prhs[0]);
  
  /* Initialize the CPLEX environment */
  env = (CPXENVptr) lpenv[0] ;
  lp=(CPXLPptr)p_lp[0] ;
  
  if ( lp != NULL ) {
    status = CPXfreeprob (env, &lp);
    if ( status ) {
      fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
    }
  }
  
  if ( status ) {
    char  errmsg[1024];
    fprintf(stderr, "Could not close CPLEX environment.\n");
    CPXgeterrorstring(env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
  }
}     




