

CPLEX_LIB ?= /fml/ag-raetsch/share/software/ilog/cplex91/lib/x86-64_RHEL3.0_3.2/static_pic/
CPLEX_INC ?= /fml/ag-raetsch/share/software/ilog/cplex91/include/ilcplex
MEX ?= /fml/ag-raetsch/share/software/matlab/bin/mex

all: cplex_init_quit.mexa64 lp_gen.mexa64 qp_gen.mexa64 lp_solve.mexa64 qp_solve.mexa64 qp_solve_start.mexa64 lp_resolve.mexa64 lp_addrows.mexa64 lp_addrows_ranged.mexa64 lp_addcols.mexa64 lp_close.mexa64 lp_write.mexa64 lp_read.mexa64 lp_set_param.mexa64 lp_delcols.mexa64 lp_delrows.mexa64 lp_chobj.mexa64 lp_chgprobtype.mexa64 lp_chgctype.mexa64 lp_cpmipstart.mexa64 getmiprelgap.mexa64

cplex_init_quit.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims cplex_init_quit.c
lp_gen.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_gen.c
qp_gen.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims qp_gen.c
lp_solve.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_solve.c
qp_solve.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims qp_solve.c
qp_solve_start.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims qp_solve_start.c
lp_resolve.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_resolve.c
lp_addrows.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_addrows.c
lp_addrows_ranged.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_addrows_ranged.c
lp_addcols.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_addcols.c
lp_close.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_close.c
lp_write.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_write.c
lp_read.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_read.c
lp_set_param.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_set_param.c
lp_delcols.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_delcols.c
lp_delrows.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_delrows.c
lp_chobj.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_chobj.c
lp_chgprobtype.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_chgprobtype.c
lp_chgctype.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_chgctype.c
lp_cpmipstart.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_cpmipstart.c
getmiprelgap.mexa64:
	${MEX} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims getmiprelgap.c

