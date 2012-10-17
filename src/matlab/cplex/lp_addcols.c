/* lpex3.c, example of using CPXaddrows to solve a problem */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  double *c=NULL, *A=NULL, 
    *l=NULL, *u=NULL, *x=NULL, *lambda=NULL ;
  int *iA=NULL, *kA=NULL, nnz=0, m=0, n=0, display=0, i=0, num_cols=0;
  long *lpenv=NULL, *p_lp=NULL;
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status, lpstat;
  double        objval;
  char **names ;
#ifndef MX_COMPAT_32
  long *iA_=NULL, *kA_=NULL ;
#endif
  
  if (nrhs > 7 || nrhs < 1) {
    mexErrMsgTxt("Usage: [how] "
		 "= lp_addcols(lpenv,p_lp,c,A,l,u,disp)");
    return;
  }
  switch (nrhs) {
  case 7:
    if (mxGetM(prhs[6]) != 0 || mxGetN(prhs[6]) != 0) {
      if (!mxIsNumeric(prhs[6]) || mxIsComplex(prhs[6]) 
	  ||  mxIsSparse(prhs[6])
	  || !(mxGetM(prhs[6])==1 && mxGetN(prhs[6])==1)) {
	mexErrMsgTxt("7th argument (display) must be "
		     "an integer scalar.");
	return;
      }
      display = *mxGetPr(prhs[6]);
    }
  case 6:
    if (mxGetM(prhs[5]) != 0 || mxGetN(prhs[5]) != 0) {
      if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5]) 
	  ||  mxIsSparse(prhs[5])
	  || !mxIsDouble(prhs[5]) 
	  ||  mxGetN(prhs[5])!=1 ) {
	mexErrMsgTxt("6th argument (u) must be "
		     "a column vector.");
	return;
      }
      u = mxGetPr(prhs[5]);
      n = mxGetM(prhs[5]);
    }
  case 5:
    if (mxGetM(prhs[4]) != 0 || mxGetN(prhs[4]) != 0) {
      if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) 
	  ||  mxIsSparse(prhs[4])
	  || !mxIsDouble(prhs[4]) 
	  ||  mxGetN(prhs[4])!=1 ) {
	mexErrMsgTxt("5th argument (l) must be "
		     "a column vector.");
	return;
      }
      if (n != 0 && n != mxGetM(prhs[4])) {
	mexErrMsgTxt("Dimension error (arg 5 and later).");
	return;
      }
      l = mxGetPr(prhs[4]);
      n = mxGetM(prhs[4]);
    }
  case 4:
    if (mxGetM(prhs[3]) != 0 || mxGetN(prhs[3]) != 0) {
      if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) 
	  || !mxIsSparse(prhs[3]) ) {
	mexErrMsgTxt("4n argument (A) must be "
		     "a sparse matrix.");
	return;
      }
      if (m != 0 && m != mxGetM(prhs[3])) {
	mexErrMsgTxt("Dimension error (arg 4 and later).");
	return;
      }
      if (n != 0 && n != mxGetN(prhs[3])) {
	mexErrMsgTxt("Dimension error (arg 4 and later).");
	return;
      }
      m = mxGetM(prhs[3]);
      n = mxGetN(prhs[3]);
      
      A = mxGetPr(prhs[3]);
#ifdef MX_COMPAT_32
      iA = mxGetIr(prhs[3]);
      kA = mxGetJc(prhs[3]);
#else
      iA_ = mxGetIr(prhs[3]);
      kA_ = mxGetJc(prhs[3]);

	  iA = myMalloc(mxGetNzmax(prhs[3])*sizeof(int)) ;
	  for (i=0; i<mxGetNzmax(prhs[3]); i++)
		  iA[i]=iA_[i] ;
	  
	  kA = myMalloc((n+1)*sizeof(int)) ;
	  for (i=0; i<n+1; i++)
		  kA[i]=kA_[i] ;
#endif
      nnz=mxGetNzmax(prhs[3]); 
    }
  case 3:
    if (mxGetM(prhs[2]) != 0 || mxGetN(prhs[2]) != 0) {
      if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
	  ||  mxIsSparse(prhs[2])
	  || !mxIsDouble(prhs[2]) 
	  ||  mxGetN(prhs[2])!=1 ) {
	mexErrMsgTxt("3st argument (c) must be "
		     "a column vector.");
	return;
      }
      if (n != 0 && n != mxGetM(prhs[2])) {
	mexErrMsgTxt("Dimension error (arg 3 and later).");
	return;
      }
      c = mxGetPr(prhs[2]);
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
		 "= lp_addcols(lpenv,p_lp,c,A,l,u,disp)");
    return;
  }
  if (display>3) fprintf(STD_OUT, "(m=%i, n=%i, nnz=%i) \n", m, n, nnz) ;
  if (display>3) fprintf(STD_OUT, "(kA[0]=%i, iA[0]=%i, A[0]=%f) \n", kA[0], iA[0], A[0]) ;
  if (display>3) fprintf(STD_OUT, "(l[0]=%f, u[0]=%f) \n", l[0], u[0]) ;
  if (display>2) fprintf(STD_OUT, "argument processing finished\n") ;

  /* Initialize the CPLEX environment */
  env = (CPXENVptr) lpenv[0] ;
  lp=(CPXLPptr)p_lp[0] ;

  if (display>2) 
    fprintf(STD_OUT, "calling CPXaddcols\n") ;
  names=CPXmalloc(n*sizeof(char*)) ;
  num_cols=CPXgetnumcols(env, lp)+1 ;
  for (i=0; i<n; i++) 
    {
      char buff[20] ;
      sprintf(buff, "x%i", i+num_cols) ;
      names[i]=CPXmalloc((strlen(buff)+1)*sizeof(char)) ;
      strcpy(names[i], buff) ;
    } ;
  status = CPXaddcols(env, lp, n, nnz, c, kA, iA, A, l, u, names);
  if ( status ) {
    fprintf (STD_OUT, "CPXaddcols failed.\n");
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

#ifndef MX_COMPAT_32
  if (iA) myFree(iA) ;
  if (kA) myFree(kA) ;
#endif

  return ;
} 
