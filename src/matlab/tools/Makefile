MEX_BIN ?= /fml/ag-raetsch/share/software/matlab-7.6/bin/mex
MATLAB_INCDIR ?= /fml/ag-raetsch/share/software/matlab/extern/include
MKOCT_BIN ?= mkoctfile
mex_incl = "-I$(MATLAB_INCDIR)"
all: int16tosparse.mexa64 int8tosparse.mex int8tosparse.mexa64 

int8tosparse.mexa64: int8tosparse.cpp
	${MEX_BIN} -g -O -largeArrayDims int8tosparse.cpp ${mex_incl}
int16tosparse.mexa64: int16tosparse.cpp
	${MEX_BIN} -g -O -largeArrayDims int16tosparse.cpp ${mex_incl}
#
int8tosparse.mex: int8tosparse.cpp
	${MKOCT_BIN} -g --mex int8tosparse.cpp

