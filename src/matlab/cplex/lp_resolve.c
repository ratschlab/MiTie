/* lpex3.c, example of using CPXaddrows to solve a problem */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"

#include "cplex.h"

#define myMalloc CPXmalloc 
#define myFree CPXfree
#define STD_OUT stderr

static err_str_number = 61 ;

static char *err_str[] = {
  "OK",
  "CPX_STAT_OPTIMAL",
  "CPX_STAT_UNBOUNDED",
  "CPX_STAT_INFEASIBLE",
  "CPX_STAT_INForUNBD",
  "CPX_STAT_OPTIMAL_INFEAS",
  "CPX_STAT_NUM_BEST",
  "CPX_STAT_FEASIBLE_RELAXED",
  "CPX_STAT_OPTIMAL_RELAXED",
  "?",
  "CPX_STAT_ABORT_IT_LIM",
  "CPX_STAT_ABORT_TIME_LIM",
  "CPX_STAT_ABORT_OBJ_LIM",
  "CPX_STAT_ABORT_USER",
  "Aborted in barrier, dual infeasible ",
  "Aborted in barrier, primal infeasible",
  "Aborted in barrier, primal and dual infeasible",
  "Aborted in barrier, primal and dual feasible",
  "Aborted in crossover",
  "Infeasible or unbounded",
  "","","","","","","","","","","","",
  "Converged, dual feasible, primal infeasible", /* 32 */
  "Converged, primal feasible, dual infeasible",
  "Converged, primal and dual infeasible",
  "Primal objective limit reached",
  "Dual objective limit reached",
  "Primal has unbounded optimal face",
  "Non-optimal solution found, primal-dual feasible",
  "Non-optimal solution found, primal infeasible",
  "Non-optimal solution found, dual infeasible",
  "Non-optimal solution found, primal-dual infeasible",
  "Non-optimal solution found, numerical difficulties",
  "Barrier found inconsistent constraints",
  "",
  "Optimal soluton with the tolerance defined by epgap or epagap has been found",
  "",
  "The limit on mixed integer solutions has been reached",
  "",
  "",
  "Time limit exceeded, but integer solution exists", /*CPXMIP_TIME_LIM_FEAS */
  "Time limit exceeded; no integer solution", /* CPXMIP_TIME_LIM_INFEAS */
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "Problem has an unbounded ray",
};

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
  int i, j;
  double *x=NULL, *lambda=NULL ;
  int display=0, m=0, n=0;
  long *lpenv=NULL, *p_lp=NULL;
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;
  int           status, lpstat;
  double        objval;
  char          opt_method[128]="auto" ;
  
  if (nrhs > 4 || nrhs < 1) {
    mexErrMsgTxt("Usage: [x,lambda,how] "
		 "= lp_solve(lpenv, p_lp, disp, method)");
    return;
  }
  switch (nrhs) {
  case 4:
    if (mxGetM(prhs[3]) != 0 || mxGetN(prhs[3]) != 0) {
      if (mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || !mxIsChar(prhs[3]) 
	  ||  mxIsSparse(prhs[3])
	  || !(mxGetM(prhs[3])==1 && mxGetN(prhs[3])>=1)) {
	mexErrMsgTxt("4th argument (method) must be "
		     "a string.");
	return;
      }
      mxGetString(prhs[3], opt_method, 128) ;
    }
  case 3:
    if (mxGetM(prhs[2]) != 0 || mxGetN(prhs[2]) != 0) {
      if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
	  ||  mxIsSparse(prhs[2])
	  || !(mxGetM(prhs[2])==1 && mxGetN(prhs[2])==1)) {
	mexErrMsgTxt("3rd argument (display) must be "
		     "an integer scalar.");
	return;
      }
      display = *mxGetPr(prhs[2]);
    }
  case 2:
    if (display>3)
      fprintf(STD_OUT,"argument #2\n") ;
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
    if (display>3)
      fprintf(STD_OUT,"argument #1\n") ;
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
  if (display>3)
    fprintf(STD_OUT,"input argument processing finished") ;

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

  status= CPXsetstrparam  (env, CPX_PARAM_WORKDIR, "/tmp/") ;
  if ( status ) {
    fprintf (STD_OUT,"Failed to set OOC work directory.\n");
    goto TERMINATE;
  }

  n=CPXgetnumcols(env, lp);
  m=CPXgetnumrows(env, lp);

  if (nlhs > 3 || nlhs < 1) {
    mexErrMsgTxt("Usage: [x,lambda,how] "
		 "= lp_solve(lpenv,p_lp,disp,method)");
    return;
  }

  if (strcmp(opt_method, "primal") &&  
	  strcmp(opt_method, "dual") &&
      strcmp(opt_method, "mip") && 
      strcmp(opt_method, "bar") && 
	  strcmp(opt_method, "hybbar") &&
	  strcmp(opt_method, "hybbar-d") &&
	  strcmp(opt_method, "hybbar-p") &&
      strcmp(opt_method, "hybnet") && 
      strcmp(opt_method, "hybnet-d") && 
      strcmp(opt_method, "hybnet-p") && 
      strcmp(opt_method, "sift") && 
      strcmp(opt_method, "lp") && 
	  strcmp(opt_method, "qp-primal") &&  
	  strcmp(opt_method, "qp-dual") &&
      strcmp(opt_method, "qp-net") && 
	  strcmp(opt_method, "qp-bar") &&
	  strcmp(opt_method, "qp-sift") &&
	  strcmp(opt_method, "qp-con") &&
	  strcmp(opt_method, "qp-auto") &&
	  strcmp(opt_method, "auto"))
	  mxErrMsgTxt("method \\in " 
				  "{'lp','primal','dual','bar','hybbar','hybnet','hybnet-p','hybnet-d','hybbar-p',hybbar-d','sift','mip','auto','qp-primal','qp-dual','qp-net','qp-sift','qp-bar','qp-con','qp-auto'}\n") ;
  
  switch (nlhs) {
  case 3:
  case 2:
    plhs[1] = mxCreateDoubleMatrix(m, 1, mxREAL);
    lambda = mxGetPr(plhs[1]);
  case 1:
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x = mxGetPr(plhs[0]);
    break;
  }
  if (display>2) 
    fprintf(STD_OUT, "output argument processing finished\n") ;

  /*fprintf(STD_OUT, "saving to /tmp/lp_cplex.txt\n") ;
  status=CPXwriteprob (env, lp, "/tmp/lp_cplex.txt", "LP") ;
  fprintf(STD_OUT, "CPXwriteprob=%i\n", status) ; */

  if (display>2) 
    fprintf(STD_OUT, "calling optimizer '%s'\n", opt_method) ;

  if (!strcmp(opt_method, "primal"))
    status = CPXprimopt (env, lp);
  else if (!strcmp(opt_method, "dual"))
    status = CPXdualopt (env, lp);
  else if (!strcmp(opt_method, "bar"))
    status = CPXbaropt (env, lp);
  else if (!strcmp(opt_method, "lp"))
	  status = CPXbaropt (env, lp);
  else if (strcmp(opt_method, "qp-primal")==0)
  {
	  status = CPXsetintparam (env, CPX_PARAM_QPMETHOD, 1);
	  status = CPXqpopt (env, lp);
  }
  else if (strcmp(opt_method, "qp-dual")==0)
  {
	  status = CPXsetintparam (env, CPX_PARAM_QPMETHOD, 2);
	  status = CPXqpopt (env, lp);
  }
  else if (strcmp(opt_method, "qp-net")==0)
  {
	  status = CPXsetintparam (env, CPX_PARAM_QPMETHOD, 3);
	  status = CPXqpopt (env, lp);
  }
  else if (strcmp(opt_method, "qp-bar")==0)
  {
	  status = CPXsetintparam (env, CPX_PARAM_QPMETHOD, 4);
	  status = CPXqpopt (env, lp);
  }
  else if (strcmp(opt_method, "qp-sift")==0)
  {	  
	  status = CPXsetintparam (env, CPX_PARAM_QPMETHOD, 5);
	  status = CPXqpopt (env, lp);
	  status = CPXqpopt (env, lp);
  }
  else if (strcmp(opt_method, "qp-con")==0)
  {
	  status = CPXsetintparam (env, CPX_PARAM_QPMETHOD, 6);
	  status = CPXqpopt (env, lp);
  }
  else if (strcmp(opt_method, "qp-auto")==0)
  {
	  status = CPXsetintparam (env, CPX_PARAM_QPMETHOD, 0);
	  status = CPXqpopt (env, lp);
  }
  else if (!strcmp(opt_method, "mip"))
	  status = CPXmipopt (env, lp);
  else if (!strcmp(opt_method, "hybbar"))
	  status = CPXhybbaropt (env, lp, CPX_ALG_NONE);
  else if (!strcmp(opt_method, "hybnet"))
	  status = CPXhybnetopt (env, lp, CPX_ALG_NONE);
  else if (!strcmp(opt_method, "hybbar-p"))
	  status = CPXhybbaropt (env, lp, CPX_ALG_PRIMAL);
  else if (!strcmp(opt_method, "hybnet-p"))
	  status = CPXhybnetopt (env, lp, CPX_ALG_PRIMAL);
  else if (!strcmp(opt_method, "hybbar-d"))
	  status = CPXhybbaropt (env, lp, CPX_ALG_DUAL);
  else if (!strcmp(opt_method, "hybnet-d"))
	  status = CPXhybnetopt (env, lp, CPX_ALG_DUAL);
  else if (!strcmp(opt_method, "sift"))
	  status = CPXsiftopt (env, lp);
  else if (!strcmp(opt_method, "auto")) {
	  status = CPXdualopt (env, lp);
	  if ( !status ) 
		  status = CPXsolution (env, lp, &lpstat, &objval, x, lambda, NULL, NULL);
	  if (status || (lpstat!=1)) {
		  if (display>1) 
			  fprintf (STD_OUT,"CPXdualopt failed (%i:%i).\n", status, lpstat);
		  status = CPXprimopt (env, lp);
		  if ( !status ) 
			  status=CPXsolution (env, lp, &lpstat, &objval, x, lambda, NULL, NULL);
		  if (status || (lpstat!=1)) {
			  if (display>1) 
				  fprintf (STD_OUT,"CPXprimopt failed (%i:%i).\n", status, lpstat);
			  status = CPXbaropt (env, lp);
			  strcpy(opt_method, "bar") ;
		  } ;
	  } ;
  }
  else 
	  fprintf(STD_OUT, "unknown method\n") ;
  
  if (display>3)
	  fprintf(STD_OUT, "CPX%sopt=%i\n", opt_method, status) ;
  if ( status ) {
	  fprintf (STD_OUT,"CPX%sopt failed.\n", opt_method);
	  goto TERMINATE;
  }
  
  if (strcmp(opt_method,"mip"))
  { /* not a MIP program */
	  if (display>2)
		  fprintf(STD_OUT, "calling CPXsolution:\n") ;
	  status = CPXsolution (env, lp, &lpstat, &objval, x, lambda, NULL, NULL);
	  
	  if (lpstat>100) lpstat-=100-43 ;
	  if ( status ) {
		  fprintf (STD_OUT,"CPXsolution=%i, lpstat=%i, terminating.\n",status, lpstat);
		  goto TERMINATE;
	  }
	  if (display>1 && lpstat>=0 && lpstat<err_str_number)
		  fprintf (STD_OUT, "Solution status: %s\n", err_str[lpstat]);
	  if (display>2)
		  fprintf (STD_OUT, "Objective value %g\n", objval);
	  
	  if (nlhs >= 3) 
		  if (lpstat==1)
			  plhs[2] = mxCreateString(err_str[0]) ;
		  else if (lpstat==44)
			  plhs[2] = mxCreateString(err_str[0]) ;
		  else if  (lpstat>=0 && lpstat<err_str_number)
			  plhs[2] = mxCreateString(err_str[lpstat]) ;
		  else
		  {
			  char buf[100] ;
			  sprintf(buf, "error code %i", lpstat) ;
			  plhs[2] = mxCreateString(buf) ;
		  }
  } else
  { /* MIP program */
	  int numit=0,cur_numrows=0,cur_numcols=0 ;
	  if (display>2)
		  fprintf(STD_OUT, "calling CPXgetmipobjval:\n") ;
	  status = CPXgetmipobjval(env, lp, &objval);
	  /*fprintf(STD_OUT, "status = %i, objval=%f\n",  status, objval) ;*/
	  
	  if ( status ) {
		  fprintf (STD_OUT,"CPXgetmipobjval=%i\n",status, lpstat);
		  goto TERMINATE;
	  }
	  numit = CPXgetmipitcnt(env, lp);
	  cur_numrows = CPXgetnumrows (env, lp);
	  cur_numcols = CPXgetnumcols (env, lp);
	  fprintf (STD_OUT,"CPXgetnumrows=%i  CPXgetnumcols=%i\n",cur_numrows,cur_numcols);
	  status = CPXgetmipx(env, lp, x, 0, cur_numcols-1);
	  if ( status ) {
		  fprintf (STD_OUT,"CPXgetmipx=%i\n",status, lpstat);
		  goto TERMINATE;
	  }
	  status = CPXgetmipslack(env, lp, lambda, 0, cur_numrows-1);
	  if ( status ) {
		  fprintf (STD_OUT,"CPXgetmipslack=%i\n",status, lpstat);
		  goto TERMINATE;
	  }
	  lpstat = CPXgetstat(env, lp);
	  if (lpstat>100) lpstat-=100-43 ;
	  
	  if ( status ) {
		  fprintf (STD_OUT,"CPXsolution=%i, lpstat=%i.\n",status, lpstat);
		  goto TERMINATE;
	  }
	  
	  if (display>1  && lpstat>=0 && lpstat<err_str_number)
		  fprintf (STD_OUT, "Solution status: %s\n", err_str[lpstat]);
	  if (display>2) {
		  fprintf (STD_OUT, "Objective value %g\n", objval);
		  fprintf (STD_OUT, "Number of Simplex Iterations: %i\n", numit);
	  } ;
	  
	  if (nlhs >= 3) 
		  if (lpstat==1)
			  plhs[2] = mxCreateString(err_str[0]) ;
		  else if (lpstat==44)
			  plhs[2] = mxCreateString(err_str[0]) ;
		  else if  (lpstat>=0 && lpstat<err_str_number)
			  plhs[2] = mxCreateString(err_str[lpstat]) ;
		  else
		  {
			  char buf[100] ;
			  sprintf(buf, "error code %i", lpstat) ;
			  plhs[2] = mxCreateString(buf) ;
		  }
  } ;
TERMINATE:
  if (status) {
	  char  errmsg[1024];
	  CPXgeterrorstring (env, status, errmsg);
	  fprintf (STD_OUT, "%s", errmsg);
	  if (nlhs >= 3) 
		  plhs[2] = mxCreateString(errmsg) ;
  } ;
  return ;
}     

