
#include <stdio.h>
#include <glpk.h>
#include "solve_lp_glpk.h"
#include "QP.h"

vector<double> solve_lp_glpk(QP* qp, bool& success)
{ 

		vector<double> res(qp->num_var, 0);

		/*declarevariables*/
		glp_prob*lp;
		int num_constraints=qp->b.size();
		double ar[num_constraints+1];
		int ia[num_constraints+1];
		int ja[num_constraints+1];

		/*createproblem*/
		lp=glp_create_prob();
		glp_set_prob_name(lp,"mitie_lp");
		glp_set_obj_dir(lp,GLP_MIN);
		/*fillproblem*/
		glp_add_rows(lp,num_constraints);
		//glp_set_row_name(lp,1,"p1");
		for(int i=0;i<num_constraints;i++)
		{
			
			if (qp->eq_idx[i]>0.5)// equality constraints
				glp_set_row_bnds(lp,i+1,GLP_FX, qp->b[i] , qp->b[i]);
			else
				glp_set_row_bnds(lp,i+1,GLP_UP, 0.0 , qp->b[i]);

			// fill constraint matrix
			vector<int>* idx = NULL;
			vector<float>* coef = NULL;
			int ret = qp->A.get(&idx, &coef, i);
			assert(idx!=NULL);
			assert(coef!=NULL);
			assert(ret>0);
			assert(idx->size()==coef->size());

			for (int j=0; j<idx->size(); j++)
			{
				ia[i] = i+1;
				ja[i] = idx->at(j)+1;
				ar[i] = coef->at(j);
			}
		}

		glp_add_cols(lp, qp->num_var);
		for(int i=0;i<qp->num_var;i++)
		{
			// set box constraints for variables
			if (qp->lb[i]>-1e10)
			{
				glp_set_col_bnds(lp,i+1,GLP_LO,qp->lb[i],0.0);
			}
			if (qp->ub[i]<1e10)
			{
				glp_set_col_bnds(lp,i+1,GLP_LO,0.0,qp->ub[i]);
			}
			glp_set_obj_coef(lp,i+1,qp->F[i]);
		}

		printf("num_constraints: %i\n", num_constraints);
		//glp_load_matrix(lp,num_constraints-1,ia,ja,ar);

		/*solveproblem*/
		glp_simplex(lp,NULL);

		/*recoveranddisplayresults*/
		double z = glp_get_obj_val(lp);
		printf("glpk obj: %.3f\n", z);

		for(int nt=0;nt<qp->num_var;nt++)
		{
			res[nt]=glp_get_col_prim(lp,nt+1);
		}
		/*housekeeping*/
		glp_delete_prob(lp);
		glp_free_env();

		return res;
}
