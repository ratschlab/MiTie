
#include <ilcplex/ilocplex.h>
#include "QP.h"
#include "solve_qp_cplex.h"

ILOSTLBEGIN

vector<double> solve_qp_cplex(QP* qp, bool& success)
{
	vector<double> res(qp->num_var, 0);

	IloEnv env;

	try {

		IloModel model(env);

		// variables 
		IloNumVarArray var(env);

		for (int i=0; i<qp->num_var; i++)
		{
			if (qp->binary_idx[i]==0)
				var.add(IloNumVar(env, qp->lb[i], qp->ub[i]));
			else
				var.add(IloNumVar(env, qp->lb[i], qp->ub[i], ILOINT));
		}

		// constraints
		IloRangeArray con(env);
		//TODO remove these fake constraints
		for (int i=0; i<qp->num_var; i++)
		{
			con.add( var[i]<= 1e10);
		}

		// iterate over constraint matrix
		for (int i=0; i<qp->b.size(); i++)
		{
			vector<int>* idx = NULL;
			vector<float>* coef = NULL;
			int ret = qp->A.get(&idx, &coef, i);
			assert(idx!=NULL);
			assert(coef!=NULL);
			assert(ret>0);
			assert(idx->size()==coef->size());

			IloNumExpr c_expr(env);
			for (int j=0; j<idx->size(); j++)
			{
				//printf("j:%i %.3f %i (num_var: %i)\n", j, coef->at(j), idx->at(j), qp->num_var);
				assert(idx->at(j)<qp->num_var);
				c_expr+= coef->at(j)*var[idx->at(j)];
			}
			if (qp->eq_idx[i]>0)
			{
				con.add(c_expr == qp->b[i]);
				//env.out() << i << c_expr << "==" << qp->b[i] << endl;
			}
			else
			{
				con.add(c_expr <= qp->b[i]);
				//env.out() << i << c_expr << "<=" << qp->b[i] << endl;
			}
		}

		// objective
		IloNumExpr obj_expr(env);
		qp->Q.reset_it();
		while (true)
		{
			int i=0;
			int j=0;
			float val = qp->Q.next(&i, &j);
			//printf("%i %i %.2f\n", i, j, val);
			if (i==-1)
				break;
			assert(i<qp->num_var && j<qp->num_var);
			obj_expr+= var[i] * var[j] * val;
		}

		for (int i=0; i<qp->F.size(); i++)
		{
			if (qp->F[i]!=0)
				obj_expr+= var[i] * qp->F[i];
		}

		IloObjective obj(env, IloMinimize(env, obj_expr));
		//model.add(IloMinimize(env, var[0] + 2 * var[1] + 3 * var[2] + var[3]));

		//env.out() << "objective:  " << IloMinimize(env, obj_expr) << endl;
		//env.out() << "constraints:  " << con << endl;
		// create model
		//model.add(var);
		model.add(con);
		model.add(obj);

		// solve...
		IloCplex cplex(model);
		cplex.solve();

		IloNumArray vals(env);
		cplex.getValues(vals, var);
		//env.out() << "Values        = " << vals << endl;
		//cplex.getSlacks(vals, con);
		//env.out() << "Slacks        = " << vals << endl;

		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << endl;


		for (int i=0; i<qp->num_var; i++)
			res[i] = vals[i];
	}
	catch (IloException& e) 
	{
		cerr << "Concert exception caught: " << e << endl;
		success = false;
	}
	catch (...)
	{
		cerr << "Unknown exception caught" << endl;
		success = false;
	}
	env.end();
	return res;
}
