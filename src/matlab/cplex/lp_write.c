/* lpex3.c, example of using CPXaddrows to solve a problem */

#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

#include "cplex.h"

#define myMalloc CPXmalloc 
#define myFree CPXfree
#define STD_OUT stderr

static char *err_str[] = {
  "OK",
  "Optimal solution found\n", /* 1 */
  "Problem infeasible\n",
  "Problem unbounded\n",
  "Objective limit exceeded in Phase II\n",
  "Iteration limit exceeded in Phase II\n",
  "Iteration limit exceeded in Phase I\n",
  "Time limit exceeded in Phase II\n",
  "Time limit exceeded in Phase I\n",
  "Problem non-optimal, singularities in Phase II\n",
  "Problem non-optimal, singularities in Phase I\n",
  "Optimal solution found, unscaled infeasibilities\n",
  "Aborted in Phase II\n",
  "Aborted in Phase I\n",
  "Aborted in barrier, dual infeasible \n",
  "Aborted in barrier, primal infeasible\n",
  "Aborted in barrier, primal and dual infeasible\n",
  "Aborted in barrier, primal and dual feasible\n",
  "Aborted in crossover  \n",
  "Infeasible or unbounded\n",
  "","","","","","","","","","","","",
  "Converged, dual feasible, primal infeasible\n", /* 32 */
  "Converged, primal feasible, dual infeasible\n",
  "Converged, primal and dual infeasible\n",
  "Primal objective limit reached\n",
  "Dual objective limit reached\n",
  "Primal has unbounded optimal face\n",
  "Non-optimal solution found, primal-dual feasible\n",
  "Non-optimal solution found, primal infeasible\n",
  "Non-optimal solution found, dual infeasible\n",
  "Non-optimal solution found, primal-dual infeasible\n",
  "Non-optimal solution found, numerical difficulties\n",
  "Barrier found inconsistent constraints\n",
  "\n",
};

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
  int i, j;
  int m=0, n=0;
  long *lpenv=NULL, *p_lp=NULL;
  char fname[1024]; 
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status, lpstat;
  double        objval;
  
  if (nrhs != 3) {
    mexErrMsgTxt("Usage: [how] "
		 "= lp_write(lpenv, p_lp, fname)");
    return;
  }
  switch (nrhs) {
  case 3:
    if (mxGetM(prhs[2]) != 0 || mxGetN(prhs[2]) != 0) {
      if (mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
	  ||  mxIsSparse(prhs[2]) || !mxIsChar(prhs[2])
	  || !(mxGetM(prhs[2])==1 && mxGetN(prhs[2])>=1)) {
	mexErrMsgTxt("3th argument (fname) must be "
		     "a string.");
	return;
      }
      mxGetString(prhs[2], fname, 1024);
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
  /* Initialize the CPLEX environment */
  env = (CPXENVptr) lpenv[0] ;
  lp=(CPXLPptr)p_lp[0] ;
  
  if (nlhs !=1) {
    mexErrMsgTxt("Usage: [how] "
		 "= lp_write(lpenv,p_lp,fname)");
    return;
  }

  fprintf(STD_OUT, "saving to %s\n", fname) ;
  status=CPXwriteprob (env, lp, fname, NULL) ;
  
 TERMINATE:
  {
    char  errmsg[1024];
    CPXgeterrorstring (env, status, errmsg);
    if (status) 
      fprintf (STD_OUT, "%s", errmsg);
    if (nlhs >= 1)
      if (status)
	plhs[0] = mxCreateString(errmsg) ;
      else
	plhs[0] = mxCreateString("OK") ;
    return ;
  }
}     
 
