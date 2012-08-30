#ifndef _GRAPH_TOOLS_H__
#define _GRAPH_TOOLS_H__
#include <assert.h>

char* get_attribute(char* line, char* tag)
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
	printf("\n");
	return result;
};

vector<char*> get_fields(char* line)
{
	vector<char*> ret;
	int i=0;
	ret.push_back(line);
	while (line[i]!=0)
	{
		if (line[i]=='\t')
		{
			line[i]=0;
			ret.push_back(line+i+1);
		}
		i++;
	}
	return ret;
};

vector<Region*> parse_gtf(char* gtf_file)
{
	
	vector<Region*> regions;
	FILE* fd = fopen(gtf_file, "r");
	if (!fd)
	{
		printf("could not open file: %s\n", gtf_file);
		return regions;
	}
	vector<segment> transcript;
	while (~feof(fd))
	{
		char line[1000];
		if (fgets(line, 1000, fd)==NULL) break;
		//printf("line: %s\n", line);

		vector<char*> fields = get_fields(line);
		//for (int i=0; i<fields.size(); i++)
		//	printf("%s\n", fields[i]);	

		if (fields.size()!=9)
		{
			printf("wrong number of fields: skip over line\n");
			continue;
		}
		int start = atoi(fields[3]);
		int stop = atoi(fields[4]);
		char strand = fields[6][0];
		char* chr = new char[100];
		strcpy(chr, fields[0]);

		if (strcmp(fields[2], "transcript")==0)
		{
			// create new region for current transcript
			//printf("%s%c:%i..%i\n", chr, strand, start, stop);
			
			if (regions.size()>0)
			{
				assert(transcript.size()>0);
				assert(regions.back()->start = transcript.front().first);
				assert(regions.back()->stop = transcript.back().second);
				regions.back()->transcripts.push_back(transcript);
				transcript.clear();
			}
			Region* reg = new Region(start, stop, chr, strand);
			regions.push_back(reg);
		}
		// parse exons
		else if (strcmp(fields[2], "exon")==0)
		{
			int start = atoi(fields[3]);
			int stop = atoi(fields[4]);
			//printf("exons: %i->%i\n", start, stop);	
			segment* seg = new segment;
			seg->first = start;
			seg->second = stop;
			transcript.push_back(*seg);
		}
		else
		{
			printf("skip over line: %s\n", line);
		}
	//	char* tr_id = get_attribute(line, "transcript_id");
	//	printf("tr_id: %s\n", tr_id);
	//	delete[] tr_id;
	}
	regions.back()->transcripts.push_back(transcript);

	return regions;
}

#endif
