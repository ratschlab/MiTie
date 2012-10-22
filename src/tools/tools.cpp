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
#include <region.h>
#include "bam.h"

vector<char*> separate(char* str, char sep)
{
	vector<char*> ret;
	int i=0;
	ret.push_back(str);
	while (str[i]!=0)
	{
		if (str[i]==sep)
		{
			str[i]=0;
			ret.push_back(str+i+1);
		}
		i++;
	}
	return ret;
}

bool compare_second(segment intr1, segment intr2)
{
	if (intr1.first<intr2.first)
		return true;
	if (intr1.first>intr2.first)
		return false;
	return (intr1.second<intr2.second);
}

vector<Region*> parse_regions(char* fn_regions)
{
	FILE* fd = fopen(fn_regions, "r");
	if (!fd)
	{
		fprintf(stderr, "Could not open file: %s for reading\n", fn_regions);
	
	}
	int num_bytes = 1000;

	vector<Region*> regions;
	while (!feof(fd))
	{
		char line[num_bytes]; 
		if (!fgets(line, num_bytes, fd))
		{
			break;
		}
		if (line[0]=='%' || line[0]=='#')
			continue;
		char* chr = new char[100];
		char* strand = new char[100];
		int start = 0;
		int stop = 0;
		int num_read = sscanf(line, "%s\t%s\t%i\t%i", chr, strand, &start, &stop);
		if (num_read!=4)
		{
			fprintf(stderr, "Error parsing line: %s\n", line);
			exit(-2);
		}
		Region* reg = new Region(start, stop, chr, strand[0]);
		regions.push_back(reg);
		delete[] strand; 
		delete[] chr;
	}

	fclose(fd);

	return regions;
}
void write_regions(vector<Region*> regions, FILE* fd)
{
	for (int i=0; i<regions.size(); i++)
	{
		fprintf(fd, "%s\t%c\t%i\t%i\n", regions[i]->chr, regions[i]->strand, regions[i]->start, regions[i]->stop);
	}
}
void set_chr_num(Region* reg, bam_header_t* header)
{
	for (int j=0; j<header->n_targets; j++)
	{
		if (strcmp(reg->chr, header->target_name[j])==0)
		{
			reg->chr_num = j;
			return;
		}
	}
	fprintf(stderr, "Did not find chr name in header: %s\n", reg->chr);
	exit(-2);
}

// interval overlapp code
////////////////////////////////////////////////////////////////////////////////
typedef struct {
    int start;
    int stop;
	int idx;
	int set_id;
} interval_t;

bool compare (interval_t i, interval_t j) 
{ 
	return (i.start<j.start); 
}

bool overlaps(interval_t a, interval_t b)
{
	int v = min(a.stop,b.stop) - max(a.start,b.start) + 1;
	return (v >= 1);
}
bool leftOf(interval_t a, interval_t b)
{
	return (a.stop < b.start);
}

void scan(interval_t f, vector<interval_t>* Wf, interval_t g, vector<interval_t>* Wg, vector<int>* overlap)
{
	vector<interval_t>::iterator i;
	i=Wg->begin();
	while (i<Wg->end())
	{
		interval_t g2 = *i;
		if (leftOf(g2,f))
		{
			Wg->erase(i);// inefficient if Wg is large
			// this moves all elements, therefore i is not incremented
		}
		else if (overlaps(g2,f))
		{
			if (g2.set_id==1)
			{
				//printf("overlap: [%i | %i, %i] [%i | %i, %i]\n", g2.idx, g2.start, g2.stop, f.idx, f.start, f.stop);
				overlap->push_back(g2.idx);
				overlap->push_back(f.idx);
			}
			else if (f.set_id==1)
			{
				//printf("overlap: [%i | %i, %i] [%i | %i, %i]\n", f.idx, f.start, f.stop, g2.idx, g2.start, g2.stop);
				overlap->push_back(f.idx);
				overlap->push_back(g2.idx);
			}	
			i++;
		}
		else
		{
			printf("never happens??\n");
			i++;
		}
	}
	if (!leftOf(f, g))
	{
		Wf->push_back(f);
		//printf("push: [%i, %i] size:%i\n", f.start, f.stop, Wf->size());
	}
}

vector<int> interval_overlap(vector<int> starts1, vector<int> stops1, vector<int> starts2, vector<int>stops2)
{	
	int num_intervals1 = starts1.size();
	assert(stops1.size()==num_intervals1);
	int num_intervals2 = starts2.size();
	assert(stops2.size()==num_intervals2);

	vector<interval_t> intervals1;
	for (int i=0; i<num_intervals1; i++)
	{
		interval_t interval;
		interval.start = starts1[i];
		interval.stop = stops1[i];
		interval.set_id = 1;
		interval.idx = i;
		intervals1.push_back(interval);
		//printf("int1: [%i, %i] \n",intervals1[i].start, intervals1[i].stop);
	}
	interval_t i;
	i.start = std::numeric_limits<int>::max();
	i.stop = std::numeric_limits<int>::max();
	i.set_id = std::numeric_limits<int>::max();
	i.idx = std::numeric_limits<int>::max();
	intervals1.push_back(i);

	//printf("num_intervals1: %i\n", intervals1.size());
	vector<interval_t> intervals2;
	for (int i=0; i<num_intervals2; i++)
	{
		interval_t interval;
		interval.start = starts2[i];
		interval.stop = stops2[i];
		interval.set_id = 2;
		interval.idx = i;
		intervals2.push_back(interval);
		//printf("int2: [%i, %i] \n",intervals2[i].start, intervals2[i].stop);
	}
	intervals2.push_back(i);
	//printf("num_intervals2: %i\n", intervals2.size());

	sort(intervals1.begin(), intervals1.end(), compare);	
	sort(intervals2.begin(), intervals2.end(), compare);	


	vector<int> overlap;
	vector<interval_t> Wx;
	vector<interval_t> Wy;
	vector<interval_t>::iterator x = intervals1.begin();
	vector<interval_t>::iterator y = intervals2.begin();
	while (x<intervals1.end() && y<intervals2.end())
	{
		//vector<interval_t>::iterator x;
		//vector<interval_t>::iterator y;
		//if (it1>intervals1.end())
		//	x = inf_interval();
		//else
		//	x = it1;
		//if (it2>intervals2.end())
		//	y = inf_interval();
		//else
		//	y=it2;

		if (x->start <= y->start)
		{
				scan(*x, &Wx, *y, &Wy, &overlap);
				x++;
		}
		else
		{
			if (y<=intervals2.end())
			{
				scan(*y, &Wy, *x, &Wx, &overlap);
				y++;
			}
		}
	}
	return overlap;
}

vector<vector<int> > region_overlap(vector<Region*> regions1, vector<Region*> regions2)
{
	// compute overlap
	vector<int> starts1;
	vector<int> stops1;
	vector<int> starts2;
	vector<int> stops2;

	for (int i=0; i<regions1.size(); i++)
	{
		starts1.push_back(regions1[i]->start);
		stops1.push_back(regions1[i]->stop);
	}
	for (int i=0; i<regions2.size(); i++)
	{
		starts2.push_back(regions2[i]->start);
		stops2.push_back(regions2[i]->stop);
	}

	vector<int> ov = interval_overlap(starts1, stops1, starts2, stops2);

	//printf("overlap.size(): %i\n", (int) ov.size());
	vector<vector<int> > ov_list(regions1.size());
	for (int i=0; i<ov.size(); i+=2)
	{
		//printf("(%i,%i) %i->%i, %i->%i\n", ov[i], ov[i+1], regions1[ov[i]]->start, regions1[ov[i]]->stop, regions2[ov[i+1]]->start, regions2[ov[i+1]]->stop);
		ov_list[ov[i]].push_back(ov[i+1]);
	}
	return ov_list;
}

