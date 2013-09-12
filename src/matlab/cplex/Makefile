

CPLEXLIBDIR ?= /cbio/grlab/share/software/ilog/CPLEX_Studio124/cplex/lib/x86-64_sles10_4.1/static_pic/
CPLEXINCDIR ?= /cbio/grlab/share/software/ilog/CPLEX_Studio124/cplex/include

CPLEX_LIB = ${CPLEXLIBDIR}
CPLEX_INC = ${CPLEXINCDIR}/ilcplex/
MEX_BIN ?= /cbio/grlab/share/software/matlab/matlab_R2012b/bin/mex
FLAGS = -g -DMX_NEW_MATLAB_CPLEX

all: cplex_init_quit.mexa64 lp_gen.mexa64 qp_gen.mexa64 lp_solve.mexa64 qp_solve.mexa64 qp_solve_start.mexa64 lp_resolve.mexa64 lp_addrows.mexa64 lp_addrows_ranged.mexa64 lp_addcols.mexa64 lp_close.mexa64 lp_write.mexa64 lp_read.mexa64 lp_set_param.mexa64 lp_delcols.mexa64 lp_delrows.mexa64 lp_chobj.mexa64 lp_chgprobtype.mexa64 lp_chgctype.mexa64 lp_cpmipstart.mexa64 getmiprelgap.mexa64

cplex_init_quit.mexa64: cplex_init_quit.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims cplex_init_quit.c
lp_gen.mexa64: lp_gen.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_gen.c
qp_gen.mexa64: qp_gen.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims qp_gen.c
lp_solve.mexa64: lp_solve.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_solve.c
qp_solve.mexa64: qp_solve.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims qp_solve.c
qp_solve_start.mexa64: qp_solve_start.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims qp_solve_start.c
lp_resolve.mexa64: lp_resolve.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_resolve.c
lp_addrows.mexa64: lp_addrows.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_addrows.c
lp_addrows_ranged.mexa64: lp_addrows_ranged.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_addrows_ranged.c
lp_addcols.mexa64: lp_addcols.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_addcols.c
lp_close.mexa64: lp_close.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_close.c
lp_write.mexa64: lp_write.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_write.c
lp_read.mexa64: lp_read.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_read.c
lp_set_param.mexa64: lp_set_param.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_set_param.c
lp_delcols.mexa64: lp_delcols.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_delcols.c
lp_delrows.mexa64: lp_delrows.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_delrows.c
lp_chobj.mexa64: lp_chobj.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_chobj.c
lp_chgprobtype.mexa64: lp_chgprobtype.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_chgprobtype.c
lp_chgctype.mexa64: lp_chgctype.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_chgctype.c
lp_cpmipstart.mexa64: lp_cpmipstart.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims lp_cpmipstart.c
getmiprelgap.mexa64: getmiprelgap.c
	${MEX_BIN} ${FLAGS} -I${CPLEX_INC} -L${CPLEX_LIB} -lcplex -lpthread -largeArrayDims getmiprelgap.c

