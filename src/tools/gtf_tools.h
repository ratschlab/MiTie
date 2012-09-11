#ifndef _GTF_TOOLS_H__
#define _GTF_TOOLS_H__
#include <assert.h>
#include "region.h"

char* get_attribute(char* line, const char* tag);

vector<char*> get_fields(char* line);

vector<Region*> parse_gtf(char* gtf_file);

bool compare_second(segment s1, segment s2);

#endif
