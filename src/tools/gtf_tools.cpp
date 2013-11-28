#include <assert.h>
#include <string>
	using std::string;
#include <string.h>
#include <stdio.h>
#include <vector>
	using std::vector;
#include <map>
	using std::map;
#include "gtf_tools.h"
#include <algorithm>
#include "region.h"
#include "tools.h"
#include "basic_tools.h"
#include <limits> 

char* get_attribute(char* line, const char* tag)
{
	int len=strlen(tag);
	char* pos = strstr(line, tag);
	if (pos==NULL)
		return NULL;

	char* start = pos+len;

	while (*start!='"' && *start !=';' && start<line+strlen(line))
	{
		start++;
	}
	start++;
	char* result = new char[101];
	
	for (int i=0; i<100; i++)
	{
		if (start+i>=line+strlen(line))
		{
			printf("%p, %p\n", start+i, line+strlen(line));
			break;
		}

		if (start[i]=='"' || start[i]==';')
			break;

		result[i] = start[i];
		result[i+1] = '\0';
	}
	//printf("\n");
	return result;
}

const char* determine_format(char* filename)
{
	return determine_format(filename, "Parent");
}
const char* determine_format(char* filename, const char* gff_link_tag)
{
	FILE* fd = fopen(filename, "r");
	if (!fd)
	{
		printf("could not open file: %s\n", filename);
		exit(-1);
	}

	int cnt = 0;
	int trid_cnt = 0;
	int parent_cnt = 0;
	const char* ret = "unknown";

	while (~feof(fd))
	{
		char line[1000];
		if (fgets(line, 1000, fd)==NULL) break;

		cnt++;

		vector<char*> fields = get_fields(line);

		if (fields.size()<9)
			continue;

		char* tr_id = strstr(fields[8], "transcript_id");
		if (tr_id)
		{
			ret = "gtf";
			trid_cnt++;
		}

		char* parent = strstr(fields[8], gff_link_tag);
		if (parent)
		{
			ret = "gff3";
			parent_cnt++;
		}

		if (parent_cnt>10 && trid_cnt==0)
		{
			ret = "gff3";
			break;
		}

		if (parent_cnt==0 && trid_cnt>10)
		{
			ret = "gtf"; 
			break;
		}

		if (cnt>100000)
			break;
	}
	return ret;
}

char* get_gff_attribute(char* line, const char* tag)
{
	int len=strlen(tag);
	char* pos = strstr(line, tag);
	if (pos==NULL)
		return NULL;

	char* start = pos+len-1;

	while (*start!='=' && start<line+strlen(line))
	{
		start++;
	}
	start++;
	char* result = new char[101];
	
	for (int i=0; i<100; i++)
	{
		if (start+i>=line+strlen(line))
		{
			//printf("%p, %p\n", start+i, line+strlen(line));
			break;
		}

		if (start[i]=='"' || start[i]==';')
			break;

		result[i] = start[i];
		result[i+1] = '\0';
	}
	//printf("\n");
	return result;
}

vector<Region*> parse_gff(char* gtf_file)
{
	return parse_gff(gtf_file, "Parent");
}

vector<Region*> parse_gff(char* gtf_file, const char* link_tag)
{
	FILE* fd = fopen(gtf_file, "r");
	if (!fd)
	{
		printf("Could not open file: %s\n", gtf_file);
		exit(-1);
	}
	int cnt = 0;
	map<string, Region*> transcripts;
	while (~feof(fd))
	{
		char line[1000];
		if (fgets(line, 1000, fd)==NULL) break;

		if (++cnt%1000==0)
			printf("\rreading line %i (%i transcripts)", cnt, (int) transcripts.size());

		vector<char*> fields = get_fields(line);
		if (fields.size()!=9)
			continue;

		char* type = fields[2];
		char strand = fields[6][0];
		char* chr = fields[0];

		//printf("%s\n", fields[8]);

		// parse exons
		if (strcmp(type, "five_prime_UTR")==0 || strcmp(type, "three_prime_UTR")==0 || strcmp(type, "CDS")==0)
		{
			char* tr_id = get_gff_attribute(fields[8], link_tag);
			if (!tr_id)
			{
				printf("Could not find ID: %s\n", fields[8]);
				exit(-1);
			}
			string transcript_id(tr_id);
			delete[] tr_id;

			int start = atoi(fields[3]);
			int stop = atoi(fields[4]);
			//printf("exons: %i->%i\n", start, stop);	
			segment* seg = new segment;
			seg->first = start;
			seg->second = stop;
			if (strcmp(type, "five_prime_UTR")==0)
				seg->flag = 5;
			else if (strcmp(type, "three_prime_UTR")==0)
				seg->flag = 3;
			else if (strcmp(type, "CDS")==0)
				seg->flag = 4;
			else	
			{
				printf("Invalid exon type: %s\n", type);
				exit(-1);
			}

			Region* reg = transcripts[transcript_id];
			if (!reg)
			{
				reg = new Region(start, stop, chr, strand);
				transcripts[transcript_id] = reg;
				vector<segment> vec;
				reg->transcripts.push_back(vec);
			}
			reg->transcripts[0].push_back(*seg);

			delete seg;
		}
	}

	vector<Region*> regions = regions_from_map(transcripts);

	regions = merge_overlapping_regions(regions);

	return regions;
}

void write_gtf(FILE* fd, Region* region, const char* source)
{
	int start = std::numeric_limits<int>::max();
	int end = 0;
	for (uint i=0; i<region->transcripts.size(); i++)
	{
		assert(region->transcripts[i].size()>0);
		if (region->transcripts[i].front().first<start)
			start = region->transcripts[i].front().first;
		if (region->transcripts[i].back().second>end)
			end = region->transcripts[i].back().second;
	}

	for (uint i=0; i<region->transcripts.size(); i++)
	{
		assert(i<region->transcript_names.size());
		start = region->transcripts[i].front().first;
		end = region->transcripts[i].back().second;
		if (region->gene_names.size()>i)
			fprintf(fd, "%s\t%s\ttranscript\t%i\t%i\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"\n", region->chr, source, start, end, region->strand, region->gene_names[i].c_str(), region->transcript_names[i].c_str());
		else
			fprintf(fd, "%s\t%s\ttranscript\t%i\t%i\t.\t%c\t.\ttranscript_id \"%s\"\n", region->chr, source, start, end, region->strand, region->transcript_names[i].c_str());

		for (uint j=0; j<region->transcripts[i].size(); j++)
		{
			start = region->transcripts[i][j].first;
			end = region->transcripts[i][j].second;
			fprintf(fd, "%s\t%s\texon\t%i\t%i\t.\t%c\t.\ttranscript_id \"%s\"\n", region->chr, source, start, end, region->strand, region->transcript_names[i].c_str());
		}
	}
}

void write_gff(FILE* fd, Region* region, const char* source)
{
	int start = std::numeric_limits<int>::max();
	int end = 0;
	for (uint i=0; i<region->transcripts.size(); i++)
	{
		assert(region->transcripts[i].size()>0);
		if (region->transcripts[i].front().first<start)
			start = region->transcripts[i].front().first;
		if (region->transcripts[i].back().second>end)
			end = region->transcripts[i].back().second;
	}
	for (uint i=0; i<region->transcripts.size(); i++)
	{
		assert(i<region->transcript_names.size());
		start = region->transcripts[i].front().first;
		end = region->transcripts[i].back().second;
		fprintf(fd, "%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=%s\n", region->chr, source, start, end, region->strand, region->transcript_names[i].c_str());
		for (uint j=0; j<region->transcripts[i].size(); j++)
		{
			char tag[1000];
			if (region->transcripts[i][j].flag==3)
				sprintf(tag, "three_prime_UTR"); 
			else if (region->transcripts[i][j].flag==4)
				sprintf(tag, "CDS"); 
			else if (region->transcripts[i][j].flag==5)
				sprintf(tag, "five_prime_UTR"); 
			else
				sprintf(tag, "exon");

			start = region->transcripts[i][j].first;
			end = region->transcripts[i][j].second;
			fprintf(fd, "%s\t%s\t%s\t%i\t%i\t.\t%c\t.\tParent=%s\n", region->chr, source, tag, start, end, region->strand, region->transcript_names[i].c_str());
		}
	}
}

vector<Region*> regions_from_map(map<string, Region*> transcripts)
{
	vector<Region*> regions;
	map<string, Region*>::iterator it;
	for (it=transcripts.begin(); it!=transcripts.end(); it++)
	{
		Region* reg = it->second;
		string transcript_id = it->first;
		
		reg->transcript_names.push_back(transcript_id);

		// sort exons
		sort(reg->transcripts[0].begin(), reg->transcripts[0].end(), compare_second);

		// warn if exons overlapp
		for (unsigned int i=1; i<reg->transcripts[0].size(); i++)
			if (reg->transcripts[0][i-1].second>=reg->transcripts[0][i].first && reg->transcripts[0][i-1].flag==reg->transcripts[0][i].flag)
				printf("\nWarning: found overlapping exons in transcript %s: (%i %i) (%i %i)\n", transcript_id.c_str(), reg->transcripts[0][i-1].first, reg->transcripts[0][i-1].second, reg->transcripts[0][i].first, reg->transcripts[0][i].second);

		// addjust start and stop
		reg->start = reg->transcripts[0].front().first;
		reg->stop = reg->transcripts[0].back().second;

		regions.push_back(reg);
	}
	return regions;
}

vector<Region*> merge_overlapping_regions(vector<Region*> regions)
{
	// merge overlapping transcripts into genes
	bool change = true;
	int iter = 0;
	printf("\n");
	while (change)
	{
		printf("parse_gtf: merge iteration: %i size:%i\n", iter++, (int) regions.size());
		if (iter>10)
			break;
		change = false;	

		vector<vector<int> > ov_list = region_overlap(regions, regions);
		for (uint i=0; i<regions.size(); i++)
		{
			for (uint j=0; j<ov_list[i].size(); j++)
			{
				// self overlap
				if (ov_list[i][j]==(int) i)
					continue;
				if (strcmp(regions[i]->chr, regions[ov_list[i][j]]->chr)!=0)
					continue;
				if (regions[i]->strand != regions[ov_list[i][j]]->strand)
					continue;
				if (regions[ov_list[i][j]]->start==-1 || regions[i]->start==-1)
					continue;

				Region* r1 = regions[i];
				Region* r2 = regions[ov_list[i][j]];
				change = true;
				// append transcripts of r2 to r1
				r1->transcripts.insert(r1->transcripts.end(), r2->transcripts.begin(), r2->transcripts.end()); 
				r1->transcript_names.insert(r1->transcript_names.end(), r2->transcript_names.begin(), r2->transcript_names.end()); 
				r1->gene_names.insert(r1->gene_names.end(), r2->gene_names.begin(), r2->gene_names.end()); 

				r1->start = std::min(r1->start, r2->start);
				r1->stop = std::max(r1->stop, r2->stop);
				r2->start = -1;
			}
		}
		// remove merged regions
		vector<Region*> tmp;
		for (uint i=0; i<regions.size(); i++)
			if (regions[i]->start>-1)
				tmp.push_back(regions[i]);
			else
				delete regions[i];

		regions = tmp;
	}
	return regions;
}

vector<Region*> parse_gtf(char* gtf_file)
{
	return parse_gtf(gtf_file, NULL);
}

vector<Region*> parse_gtf(char* gtf_file, char* gene_name)
{
	FILE* fd = fopen(gtf_file, "r");
	if (!fd)
	{
		printf("could not open file: %s\n", gtf_file);
		exit(-1);
	}
	int cnt = 0;
	map<string, Region*> transcripts;
	while (~feof(fd))
	{
		char line[1000];
		if (fgets(line, 1000, fd)==NULL) break;

		if (++cnt%1000==0)
			printf("\rreading line %i (%i transcripts)", cnt, (int) transcripts.size());

		if (line[0] == '#') // disregard comments
			continue;

		vector<char*> fields = get_fields(line);

		if (fields.size()!=9)
		{
			printf("wrong number of fields (%i): skip over line\n", (int) fields.size());
			for (int i=0; i<(int)fields.size(); i++)
			{
				printf("%s\t", fields[i]);
			}
			printf("\n");
			continue;
		}

		char* type = fields[2];
		char strand = fields[6][0];
		char* chr = fields[0];

		// parse exons
		if (strcmp(type, "exon")==0 || strcmp(type, "CDS")==0)
		{
			char* tr_id = get_attribute(fields[8], "transcript_id");
			char* gene_id = get_attribute(fields[8], "gene_id");
			char* c_gene_name = get_attribute(fields[8], "gene_name");
			if (!tr_id)
			{
				printf("Could not find transcript_id: %s\n", fields[8]);
				exit(-1);
			}

			if (gene_name)
			{
				bool skip=true;
				if (gene_id && strcmp(gene_name, gene_id)==0)
					skip=false;
				if (c_gene_name && strcmp(gene_name, c_gene_name)==0)
					skip=false;
				if (strcmp(gene_name, tr_id)==0)
					skip=false;
				if (skip)
					continue;
			}

			string transcript_id(tr_id);
			delete[] tr_id;

			int start = atoi(fields[3]);
			int stop = atoi(fields[4]);
			//printf("exons: %i->%i\n", start, stop);	
			segment* seg = new segment;
			seg->first = start;
			seg->second = stop;
			if (strcmp(type, "exon")==0)
				seg->flag = 5;
			else if (strcmp(type, "CDS")==0)
				seg->flag = 4;
			else
				assert(false);

			Region* reg = transcripts[transcript_id];
			if (!reg)
			{
				reg = new Region(start, stop, chr, strand);
				transcripts[transcript_id] = reg;
				vector<segment> vec;
				reg->transcripts.push_back(vec);
				if (gene_id)
					reg->gene_names.push_back(string(gene_id));
				else 
					reg->gene_names.push_back(string("-"));
			}
			reg->transcripts[0].push_back(*seg);

			transcript_id.clear();
			delete seg;
		}
		//else
		//{
		//	printf("skip over line: %s\n", line);
		//}
	}

	// remove duplicated exons ("exon" and "CDS")
	map<string, Region*>::iterator it;
	cnt = 0;
	for (it=transcripts.begin(); it!=transcripts.end(); it++)
	{
		cnt++;
		Region* reg = it->second;
		// sort exons
		sort(reg->transcripts[0].begin(), reg->transcripts[0].end(), compare_second);

		// discard duplicated exons and split CDS/UTR exons

		vector<segment>::iterator tit = reg->transcripts[0].begin();
		tit++;
		for (; tit != reg->transcripts[0].end(); tit++)
		{
			if ((tit-1)->first==tit->first)
			{
				if (!((tit-1)->flag==4 || tit->flag==4))
				{
					// one should be a CDS exon
					printf("%s (cnt:%i): (%i %i %i)(%i %i %i)\n", it->first.c_str(), cnt, (tit-1)->first, (tit-1)->second, (tit-1)->flag, tit->first, tit->second, tit->flag);
					exit(-1);
				}
				if ((tit-1)->second==tit->second)
				{
					tit->flag = -1;//erase
					(tit-1)->flag=4;
				}
				else if ((tit-1)->flag==4)
				{
					assert((tit-1)->second<tit->second);
					tit->first = (tit-1)->second+1;
				}
				else if (tit->flag==4)
				{
					assert(tit->second<(tit-1)->second);
					(tit-1)->first = tit->second+1;
				}
			}
			else if ((tit-1)->second==tit->second)
			{
				// first is not identical
				assert(((tit-1)+1)->flag = 4);
				(tit-1)->second = tit->first-1;
			}
			else if ((tit-1)->second>tit->second)
			{
				// cds is included in a single exon

				if (tit->flag != 4)
				{
					printf("\ntit-1: %i %i %i\n", (tit-1)->first, (tit-1)->second, (tit-1)->flag);
					printf("tit:   %i %i %i\n", tit->first, tit->second, tit->flag);
				}
				assert(tit->flag==4);
				{
					// create three exons
					segment* seg = new segment();
					seg->first = tit->second+1;
					seg->second = (tit-1)->second;
					seg->flag = 5;
					(tit-1)->second = tit->first-1;

					tit = reg->transcripts[0].insert(tit+1, *seg);

					//printf("insert segment: %i %i %i\n", seg->first, seg->second, seg->flag);
					delete seg;
					//for (int jj=0; jj<reg->transcripts[0].size(); jj++)
					//	printf("seg: %i %i (%i)\n", reg->transcripts[0][jj].first, reg->transcripts[0][jj].second, reg->transcripts[0][jj].flag);
					//exit(-1);
				}
			}
		}

		// erase exons flagged with -1
		// and fix coding flags
		vector<segment> trans;
		int last_flag = -1;
		for (tit = reg->transcripts[0].begin(); tit != reg->transcripts[0].end(); tit++)
		{
			if (tit->flag != -1)
				trans.push_back(*tit);

			// gtf files do not distinguish between 5' and 3' UTRs 
			if (reg->strand == '+')
			{
				if ((last_flag == 4 || last_flag == 3) && tit->flag == 5)
				{
					tit->flag = 3; 
				}
			}
			if (reg->strand == '-')
			{
				if ((last_flag == -1 || last_flag == 3) && tit->flag == 5)
				{
					tit->flag = 3; 
				}
			}
			last_flag = tit->flag;
		}
		//reg->transcripts[0].clear();
		reg->transcripts[0] = trans;
	}
	
	vector<Region*> regions = regions_from_map(transcripts);
	
	regions = merge_overlapping_regions(regions);

	return regions;
}
