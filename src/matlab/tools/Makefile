MEX_BIN ?= /fml/ag-raetsch/share/software/matlab-7.6/bin/mex
MATLAB_INCDIR ?= /fml/ag-raetsch/share/software/matlab/extern/include
MKOCT_BIN ?= mkoctfile
mex_incl = "-I$(MATLAB_INCDIR)"

all: octave matlab
matlab: factln.mexa64 int16tosparse.mexa64 int8tosparse.mexa64 interval_overlap.mexa64 
octave: factln.mex int16tosparse.mex int8tosparse.mex interval_overlap.mex

int8tosparse.mexa64: int8tosparse.cpp
	${MEX_BIN} -g -O -largeArrayDims int8tosparse.cpp ${mex_incl}
int16tosparse.mexa64: int16tosparse.cpp
	${MEX_BIN} -g -O -largeArrayDims int16tosparse.cpp ${mex_incl}
interval_overlap.mexa64: interval_overlap.cpp
	${MEX_BIN} -g -O interval_overlap.cpp ${mex_incl}
factln.mexa64: factln.cpp
	${MEX_BIN} -g -O factln.cpp ${mex_incl}

#
int16tosparse.mex: int16tosparse.cpp
	${MKOCT_BIN} -g --mex int16tosparse.cpp
int8tosparse.mex: int8tosparse.cpp
	${MKOCT_BIN} -g --mex int8tosparse.cpp
interval_overlap.mex: interval_overlap.cpp
	${MKOCT_BIN} -g --mex interval_overlap.cpp
factln.mex: factln.cpp
	${MKOCT_BIN} -g --mex factln.cpp

