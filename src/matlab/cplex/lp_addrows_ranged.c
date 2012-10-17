/* lpex3.c, example of using CPXaddrows to solve a problem */

#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
/*#include "mcc.h"*/

#include "cplex.h"

#define myMalloc CPXmalloc 
#define myFree CPXfree
#define STD_OUT stderr

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
  int i, j, nrow;
  double *b=NULL, *A=NULL, 
    *l=NULL, *u=NULL, *x=NULL, *lambda=NULL, *r=NULL;
  int *iA=NULL, *kA=NULL, nnz=0, m=0, n=0, display=0, *ind;
  long *lpenv=NULL, *p_lp=NULL;
  char *Sense=NULL ;
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status, lpstat;
  double        objval;
#ifndef MX_COMPAT_32
  long *iA_=NULL, *kA_=NULL ;
#endif
  
  if (nrhs > 6 || nrhs < 1) {
    mexErrMsgTxt("Usage: [how] "
		 "= lp_addrows_ranged(env,lp,A,b,r,disp) (b <= A x <= b+r)");
    return;
  }

  switch (nrhs) {
  case 6:
    if (mxGetM(prhs[5]) != 0 || mxGetN(prhs[5]) != 0) {
      if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5]) 
	  ||  mxIsSparse(prhs[5])
	  || !(mxGetM(prhs[5])==1 && mxGetN(prhs[5])==1)) {
	mexErrMsgTxt("6th argument (display) must be "
		     "an integer scalar.");
	return;
      }
      display = *mxGetPr(prhs[5]);
    }
  case 5:
    if (mxGetM(prhs[4]) != 0 || mxGetN(prhs[4]) != 0) {
      if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) 
	  ||  mxIsSparse(prhs[4])
	  || !mxIsDouble(prhs[4]) 
	  ||  mxGetN(prhs[4])!=1 ) {
	mexErrMsgTxt("5th argument (r) must be "
		     "a column vector.");
	return;
      }
    }
  case 4:
    if (mxGetM(prhs[3]) != 0 || mxGetN(prhs[3]) != 0) {
      if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) 
	  ||  mxIsSparse(prhs[3])
	  || !mxIsDouble(prhs[3]) 
	  ||  mxGetN(prhs[3])!=1 ) {
	mexErrMsgTxt("4th argument (b) must be "
		     "a column vector.");
	return;
      }
      if (m != 0 && m != mxGetM(prhs[3])) {
	mexErrMsgTxt("Dimension error (arg 4 and later).");
	return;
      }
      b = mxGetPr(prhs[3]);
      r = mxGetPr(prhs[4]);
      m = mxGetM(prhs[3]);
    }
  case 3:
    if (mxGetM(prhs[2]) != 0 || mxGetN(prhs[2]) != 0) {
      if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
	  || !mxIsSparse(prhs[2]) ) {
	mexErrMsgTxt("3n argument (A) must be "
		     "a sparse matrix.");
	return;
      }
      if (m != 0 && m != mxGetN(prhs[2])) {
	mexErrMsgTxt("Dimension error (arg 3 and later).");
	return;
      }
      if (n != 0 && n != mxGetM(prhs[2])) {
	mexErrMsgTxt("Dimension error (arg 3 and later).");
	return;
      }
      m = mxGetN(prhs[2]);
      n = mxGetM(prhs[2]);
      
      A = mxGetPr(prhs[2]);
#ifdef MX_COMPAT_32
      iA = mxGetIr(prhs[2]);
      kA = mxGetJc(prhs[2]);
#else
      iA_ = mxGetIr(prhs[2]);
      kA_ = mxGetJc(prhs[2]);

	  iA = myMalloc(mxGetNzmax(prhs[2])*sizeof(int)) ;
	  for (i=0; i<mxGetNzmax(prhs[2]); i++)
		  iA[i]=iA_[i] ;
	  
	  kA = myMalloc((n+1)*sizeof(int)) ;
	  for (i=0; i<n+1; i++)
		  kA[i]=kA_[i] ;
#endif

      /*{
	int k=0 ;
	mxArray* a = mxCreateDoubleMatrix(1, 1, mxREAL);
	for (; k<28; k++)
	  printf("%i,", iA[k]) ;
	printf("\n") ;
	for (k=0; k<28; k++)
	  printf("%i,", kA[k]) ;
	printf("\n") ;
	for (k=0; k<28; k++)
	  printf("%f,", A[k]) ;
	printf("\n") ;
	}*/
      /*nnz=mxGetNzmax(prhs[2]); */
      nnz=kA[m] ;
      if (display>3)
	fprintf(STD_OUT, "nnz=%i\n", nnz) ;

      Sense=myMalloc((m+1)*sizeof(char)) ;
      for (i=0; i<m; i++)
	Sense[i] = 'R';

      
      Sense[m]=0 ;
      if (display>3)
	fprintf(STD_OUT, "Sense=%s\n", Sense) ;
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
		 "= lp_addrows_range(lpenv,p_lp,A,b,r,disp)");
    return;
  }
  if (display>3) fprintf(STD_OUT, "(m=%i, n=%i, nnz=%i) \n", m, n, nnz) ;

  /* Initialize the CPLEX environment */
  env = (CPXENVptr) lpenv[0] ;
  lp=(CPXLPptr)p_lp[0] ;


  nrow = CPXgetnumrows(env, lp);
  ind = (int *) myMalloc(nrow * sizeof(int));
  
  for (i=0; i < m; i++)
    ind[i] = nrow + i;
  
  if (display>2) 
    fprintf(STD_OUT, "calling CPXaddrows \n") ;
  status = CPXaddrows (env, lp, 0, m, nnz, b, Sense, kA, iA, A, NULL, NULL);
  if ( status ) {
    fprintf (STD_OUT,"CPXaddrows failed.\n");
    goto TERMINATE;
  }
  
  status = CPXchgrngval (env, lp, m, ind, r);
  if ( status ) {
    fprintf (STD_OUT,"CPXchgrngval failed.\n");
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
  /*  if (Sense) myFree(Sense) ;*/
#ifndef MX_COMPAT_32
  if (iA) myFree(iA) ;
  if (kA) myFree(kA) ;
#endif

  return ;
}     
