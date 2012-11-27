
#### please adjust these paths or define the variables accordingly before starting make
export CPLEX_LIB ?= /fml/ag-raetsch/share/software/ilog/cplex91/lib/x86-64_RHEL3.0_3.2/static_pic/
export CPLEX_INC ?= /fml/ag-raetsch/share/software/ilog/cplex91/include/ilcplex
export MEX_BIN ?= /fml/ag-raetsch/share/software/matlab-7.6/bin/mex


#CXX = g++ -O9
CXX = g++ -g
NGS = tools/ngs
NGS_TOOLS = ${NGS}/get_reads_direct.cpp ${NGS}/read.cpp 
SAMDIR = tools/samtools
SAMINCL = -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses
GENOME = tools/genome/
GENOME_TOOLS = ${GENOME}/genome.cpp ${GENOME}/genome_tools.cpp
TOOLS_DIR = tools/
TOOLS = ${TOOLS_DIR}file_stats.cpp ${TOOLS_DIR}gtf_tools.cpp ${TOOLS_DIR}tools.cpp 
INCL = -I. -I${NGS} -I${GENOME} -I${TOOLS_DIR} ${SAMINCL}

#### Mex stuff
MATLAB_INCDIR ?= /fml/ag-raetsch/share/software/matlab/extern/include
MEX_UTILS = tools/mex/




all: subfolders define_regions generate_segment_graph ath_count_reads_in_region eval_pacbio load_regions_bin.mexa64

subfolders:
	@echo "##################################################" 
	@echo "# compiling samtools"
	@echo "##################################################" 
	$(MAKE) -C tools/samtools
	@echo "##################################################" 
	@echo "# compiling mex tools"
	@echo "##################################################" 
	$(MAKE) -C matlab/tools
	@echo "##################################################" 
	@echo "# compiling cplex interface"
	@echo "##################################################" 
	$(MAKE) -C matlab/cplex

define_regions: define_regions.cpp region.cpp ${TOOLS} ${GENOME_TOOLS} ${NGS_TOOLS}
	$(CXX) define_regions.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS} ${INCL} -o define_regions

generate_segment_graph: generate_segment_graph.cpp region.cpp ${TOOLS} ${GENOME_TOOLS} ${NGS_TOOLS}
	$(CXX) generate_segment_graph.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS} ${INCL} -o generate_segment_graph

ath_count_reads_in_region: ath_count_reads_in_region.cpp region.cpp ${TOOLS} ${GENOME_TOOLS} ${NGS_TOOLS}
	$(CXX) ath_count_reads_in_region.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS} ${INCL} -o ath_count_reads_in_region

eval_pacbio: eval_pacbio.cpp region.cpp ${TOOLS} ${GENOME_TOOLS} ${NGS_TOOLS} 
	$(CXX) eval_pacbio.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS} ${INCL} -o eval_pacbio

load_regions_bin.mexa64: load_regions_bin.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS}
	${MEX_BIN} -g -O load_regions_bin.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS} -I${MEX_UTILS} -I${MATLAB_INCDIR} ${INCL}

load_regions_bin_MAGIC.mexa64: load_regions_bin_MAGIC.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS}
	${MEX_BIN} -g -O load_regions_bin_MAGIC.cpp region.cpp ${TOOLS} ${NGS_TOOLS} ${GENOME_TOOLS} -I${MEX_UTILS} -I${MATLAB_INCDIR} ${INCL}
