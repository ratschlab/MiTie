
CXX = g++ -O9
#CXX = g++ -g
NGS = tools/ngs
NGS_TOOLS = ${NGS}/get_reads_direct.cpp ${NGS}/read.cpp 
SAMDIR = tools/samtools
SAMINCL = -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses
GENOME = tools/genome/
TOOLS = tools/
INCL = -I${NGS} -I${GENOME} -I${TOOLS} ${SAMINCL}

#### Mex stuff
MEX_BIN ?= /fml/ag-raetsch/share/software/matlab-7.6/bin/mex
MATLAB_INCDIR ?= /fml/ag-raetsch/share/software/matlab/extern/include
MEX_UTILS = tools/mex/



all: subfolders define_regions generate_segment_graph ath_count_reads_in_region load_regions_bin.mexa64

subfolders:
	@echo "##################################################" 
	@echo "# compiling samtools"
	@echo "##################################################" 
	$(MAKE) -C tools/samtools
	@echo "##################################################" 
	@echo "# compiling mex tools"
	@echo "##################################################" 
	$(MAKE) -C matlab/tools

define_regions: define_regions.cpp region.cpp ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${NGS_TOOLS}
	$(CXX) define_regions.cpp region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${INCL} -o define_regions

generate_segment_graph: generate_segment_graph.cpp region.cpp ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${NGS_TOOLS}
	$(CXX) generate_segment_graph.cpp region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${INCL} -o generate_segment_graph

ath_count_reads_in_region: ath_count_reads_in_region.cpp region.cpp ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${NGS_TOOLS}
	$(CXX) ath_count_reads_in_region.cpp region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp ${INCL} -o ath_count_reads_in_region

load_regions_bin.mexa64: load_regions_bin.cpp region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp
	${MEX_BIN} -g -O load_regions_bin.cpp region.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp -I${MEX_UTILS} -I${MATLAB_INCDIR} ${INCL}
