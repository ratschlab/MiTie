#ifndef _TOOLS_H__
#define _TOOLS_H__
#include <assert.h>
#include <math.h>
#include <limits>
#include <assert.h>
#include <vector>
  using std::vector;
#include <algorithm>
  using std::sort;
  using std::min;
  using std::max;
//#ifdef __USE_BAM__
#ifdef BAM_BAM_H
	#include <bam_region.h>
#endif
//	#include "bam.h"
//#endif

bool compare_second(segment intr1, segment intr2);

vector<char*> separate(char* str, char sep);

void write_regions(vector<Region*> regions, FILE* fd);

// interval overlapp code
////////////////////////////////////////////////////////////////////////////////
typedef struct {
    int start;
    int stop;
	int idx;
	int set_id;
} interval_t;

bool compare (interval_t i, interval_t j); 

bool overlaps(interval_t a, interval_t b);

bool leftOf(interval_t a, interval_t b);

void scan(interval_t f, vector<interval_t>* Wf, interval_t g, vector<interval_t>* Wg, vector<int>* overlap);

vector<int> interval_overlap(vector<int> starts1, vector<int> stops1, vector<int> starts2, vector<int>stops2);

vector<vector<int> > region_overlap(vector<Region*> regions1, vector<Region*> regions2);

#ifdef BAM_BAM_H
	vector<vector<int> > region_overlap(vector<Bam_Region*> regions1, vector<Bam_Region*> regions2);

	vector<Bam_Region*> parse_regions(char* fn_regions);

	void set_chr_num(Region* reg, bam_header_t* header);
#endif

#endif 
