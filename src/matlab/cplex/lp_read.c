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
  char fname[1024] ;
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status, lpstat;
  double        objval;
  
  if (nrhs !=2) {
    mexErrMsgTxt("Usage: [p_lp, how] "
		 "= lp_read(lpenv, fname)");
    return;
  }
  switch (nrhs) {
  case 2:
    if (mxGetM(prhs[1]) != 0 || mxGetN(prhs[1]) != 0) {
      if (mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
	  ||  mxIsSparse(prhs[1]) || !mxIsChar(prhs[1])
	  || !(mxGetM(prhs[1])==1 && mxGetN(prhs[1])>=1)) {
	mexErrMsgTxt("2nd argument (fname) must be "
		     "a string.");
	return;
      }
      mxGetString(prhs[1], fname, 1024);
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
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  p_lp = (long*) mxGetPr(plhs[0]);

  /* Initialize the CPLEX environment */
  env = (CPXENVptr) lpenv[0] ;

  /* Create the problem */    
  fprintf(STD_OUT, "calling CPXcreateprob \n") ;
  lp = CPXcreateprob (env, &status, "xxx");
  if ( lp == NULL ) {
    fprintf (STD_OUT,"Failed to create subproblem\n");
    status = 1;
    goto TERMINATE;
  } 
  if (p_lp) *p_lp=(long) lp ;
  
  if (nlhs > 2 || nlhs < 1) {
    mexErrMsgTxt("Usage: [p_lp, how] "
		 "= lp_read(lpenv,fname)");
    return;
  }
  fprintf(STD_OUT, "loading from %s\n", fname) ;
  status=CPXreadcopyprob (env, lp, fname, NULL) ;
  
 TERMINATE:
  {
    char  errmsg[1024];
    CPXgeterrorstring (env, status, errmsg);
    if (status) 
      fprintf (STD_OUT, "%s", errmsg);
    if (nlhs >= 1)
      if (status)
	plhs[1] = mxCreateString(errmsg) ;
      else
	plhs[1] = mxCreateString("OK") ;
    return ;
  }
}     
 

