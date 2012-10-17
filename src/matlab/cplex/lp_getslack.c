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
  "Optimal integer solution found\n",
  "Optimal sol. within epgap or epagap tolerance found\n",
  "Integer infeasible\n",
  "Mixed integer solutions limit exceeded\n",
  "Node limit exceeded, integer solution exists\n",
  "Node limit exceeded, no integer solution\n",
  "Time limit exceeded, integer solution exists\n",  
  "Time limit exceeded, no integer solution\n",
  "Error termination, integer solution exists\n",
  "Error termination, no integer solution\n",
  "Treememory limit, integer solution exists\n", 
  "Treememory limit, no integer solution exists\n", 
  "Aborted, integer solution exists\n",
  "Aborted, no integer solution\n",  
  "Problem optimal with unscaled infeasibilities\n",
  "Out of memory, no tree, integer solution exists\n",
  "Out of memory, no tree, no integer solution\n",  
  "Node file size limit, integer solution exists\n",  
  "Node file size limit, no integer solution\n",  
  "\n",
};

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
  int i, j;
  double *slack=NULL;
  int display=0, m=0, n=0;
  long *lpenv=NULL, *p_lp=NULL;
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status;
  
  if (nrhs != 2 || nrhs < 1) {
    mexErrMsgTxt("Usage: [slack,how] "
		 "= lp_getslack(lpenv, p_lp)");
    return;
  }
  
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

  /* Initialize the CPLEX environment */
  env = (CPXENVptr) lpenv[0] ;
  lp=(CPXLPptr)p_lp[0] ;
  
  /* Turn on output to the screen */
  if (display>0)
    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
  else
    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
  if ( status ) {
    fprintf (STD_OUT, 
	     "Failure to turn on screen indicator, error %d.\n", status);
    goto TERMINATE;
  }
  if (display>2)
    status = CPXsetintparam (env, CPX_PARAM_SIMDISPLAY, 2);
  else 
    status = CPXsetintparam (env, CPX_PARAM_SIMDISPLAY, display);
  
  if ( status ) {
    fprintf (STD_OUT,"Failed to turn up simplex display level.\n");
    goto TERMINATE;
  }
  
  n=CPXgetnumcols(env, lp);
  m=CPXgetnumrows(env, lp);


  plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
  slack = mxGetPr(plhs[0]);
 
  status = CPXgetslack (env, lp, slack, 0, m-1);

 if (status) {
    char  errmsg[1024];
    CPXgeterrorstring (env, status, errmsg);
    fprintf (STD_OUT, "%s", errmsg);
    if (nlhs >= 1) 
      plhs[1] = mxCreateString(errmsg) ;
  } else
    if (nlhs >= 1) 
      plhs[1] = mxCreateString("OK") ;
  ;

 TERMINATE:
  return ;
}     
 












