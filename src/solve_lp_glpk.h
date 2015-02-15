#ifndef _SOLVE_LP_GLPK_H__
#define _SOLVE_LP_GLPK_H__

#include "QP.h"
#include <vector>
	using std::vector;

void simplify(QP* qp);	
vector<double> solve_lp_glpk(QP* qp, bool& success);	

#endif
