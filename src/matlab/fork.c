/* lpex3.c, example of using CPXaddrows to solve a problem */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include <unistd.h>

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
  int num=0 ;
  
  if (nrhs != 1 || nlhs !=1) {
    mexErrMsgTxt("Usage: [id] "
		 "= fork(num)");
    return;
  }
  if (mxGetM(prhs[0]) != 0 || mxGetN(prhs[0]) != 0) 
  {
      if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		  ||  mxIsSparse(prhs[0])
		  || !(mxGetM(prhs[0])==1 && mxGetN(prhs[0])==1)) {
		  mexErrMsgTxt("1st argument (num) must be "
					   "an integer scalar.");
		  return;
      }
      num = *mxGetPr(prhs[0]);
	  if (num<0 || num>50)
		  mexErrMsgTxt("num too small or too large");
  }
  fprintf(stdout, "forking into %i processes\n", num) ;
  {
	  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	  double *x = mxGetPr(plhs[0]);

	  int i=0 ;
	  for (i=1; i<num; i++)
	  {
		  pid_t pid=fork() ;
		  if (pid<0)
			  mexErrMsgTxt("fork failed") ;
		  
		  if (pid==0)
		  {
			  *x = i ;
			  return ;
		  }
	  }
	  *x=0.0 ;
  }

  return ;
}     

