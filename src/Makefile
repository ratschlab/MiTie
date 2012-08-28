
CXX = g++ -O9
#CXX = g++ -g
NGS = ../../../../tools/ngs
NGS_TOOLS = ${NGS}/get_reads_direct.cpp ${NGS}/read.cpp 
SAMDIR = /fml/ag-raetsch/share/software/samtools.svn_up
SAMINCL = -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses
SKIMMER = ../../../transcript_skimmer/src
GENOME = ../../../mGene_core/genome/
TOOLS = .
INCL = -I${NGS} -I${SKIMMER} -I${GENOME} -I${TOOLS} ${SAMINCL}

#### Mex stuff
MEX_BIN ?= /fml/ag-raetsch/share/software/matlab-7.6/bin/mex
MATLAB_INCDIR ?= /fml/ag-raetsch/share/software/matlab/extern/include
MEX_UTILS=${HOME}/svn/projects/mGene_core/tools/



all: ath_count_reads_in_region define_regions generate_segment_graph load_regions_bin.mexa64

define_regions: define_regions.cpp ${SKIMMER}/region.cpp ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${NGS_TOOLS}
	$(CXX) define_regions.cpp ${SKIMMER}/region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${INCL} -o define_regions

generate_segment_graph: gtf_tools.h generate_segment_graph.cpp ${SKIMMER}/region.cpp ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${NGS_TOOLS}
	$(CXX) generate_segment_graph.cpp ${SKIMMER}/region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${INCL} -o generate_segment_graph

ath_count_reads_in_region: ath_count_reads_in_region.cpp ${SKIMMER}/region.cpp ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${NGS_TOOLS}
	$(CXX) ath_count_reads_in_region.cpp ${SKIMMER}/region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${INCL} -o ath_count_reads_in_region

load_regions_bin.mexa64: load_regions_bin.cpp ${SKIMMER}/region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp
	${MEX_BIN} -g -O load_regions_bin.cpp ${SKIMMER}/region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp -I${MEX_UTILS} -I${MATLAB_INCDIR} ${INCL}
