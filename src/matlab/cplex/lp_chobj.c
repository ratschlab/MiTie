/* lpex3.c, example of using CPXaddrows to solve a problem */

#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

#include "cplex.h"

#define myMalloc CPXmalloc 
#define myFree CPXfree
#define STD_OUT stderr

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
  double *values=NULL, *idxd=NULL ;
  int *idx=NULL, n=0, display=0, i=0;
  long *lpenv=NULL, *p_lp=NULL;
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status, lpstat;
  double        objval;
  char **names ;
  
  if (nrhs > 7 || nrhs < 1) {
    mexErrMsgTxt("Usage: [how] "
		 "= lp_chobj(lpenv,p_lp,idxs,values,disp)");
    return;
  }
  switch (nrhs) {
  case 5:
    if (mxGetM(prhs[4]) != 0 || mxGetN(prhs[4]) != 0) {
      if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) 
	  ||  mxIsSparse(prhs[4])
	  || !(mxGetM(prhs[4])==1 && mxGetN(prhs[4])==1)) {
	mexErrMsgTxt("5th argument (display) must be "
		     "an integer scalar.");
	return;
      }
      display = *mxGetPr(prhs[4]);
    }
  case 4:
    if (mxGetM(prhs[3]) != 0 || mxGetN(prhs[3]) != 0) {
      if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) 
	  || mxIsSparse(prhs[3]) || mxGetN(prhs[2])!=1) {
	mexErrMsgTxt("4th argument (values) must be "
		     "a column vector.");
	return;
      }
      n = mxGetN(prhs[3]);
      values = mxGetPr(prhs[3]);
    }
  case 3:
    if (mxGetM(prhs[2]) != 0 || mxGetN(prhs[2]) != 0) {
      if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
	  ||  mxIsSparse(prhs[2])
	  || !mxIsDouble(prhs[2]) 
	  ||  mxGetN(prhs[2])!=1 ) {
	mexErrMsgTxt("3st argument (idx) must be "
		     "a column vector.");
	return;
      }
      if (n != 0 && n != mxGetM(prhs[2])) {
	mexErrMsgTxt("Dimension error (arg 3 and later).");
	return;
      }
      idxd = mxGetPr(prhs[2]);
      n = mxGetM(prhs[2]);
    }
  case 2:
    if (mxGetM(prhs[1]) != 0 || mxGetN(prhs[1]) != 0) {
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
    }
  case 1:
    if (mxGetM(prhs[0]) != 0 || mxGetN(prhs[0]) != 0) {
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
    }
  }
  
  if (nlhs > 1 || nlhs < 1) {
    mexErrMsgTxt("Usage: [how] "
		 "= lp_chobj(lpenv,p_lp,idxs,values,disp)");
    return;
  }
  if (display>3) fprintf(STD_OUT, "(n=%i) \n", n) ;
  if (display>2) fprintf(STD_OUT, "argument processing finished\n") ;

  /* Initialize the CPLEX environment */
  env = (CPXENVptr) lpenv[0] ;
  lp=(CPXLPptr)p_lp[0] ;

  if (display>2) 
    fprintf(STD_OUT, "calling CPXchgobj\n") ;
  idx=CPXmalloc(n*sizeof(int)) ;
  for (i=0; i<n; i++) 
      idx[i]=(int) idxd[i];
  status = CPXchgobj(env, lp, n, idx, values);
  if ( status ) {
    fprintf (STD_OUT, "CPXchgobj failed.\n");
    goto TERMINATE;
  }

 TERMINATE:
  if (status) {
    char  errmsg[1024];
    CPXgeterrorstring (env, status, errmsg);
    fprintf (STD_OUT, "%s", errmsg);
    if (nlhs >= 1) 
      plhs[0] = mxCreateString(errmsg) ;
  } else
    if (nlhs >= 1) 
      plhs[0] = mxCreateString("OK") ;
  ;
  return ;
} 
