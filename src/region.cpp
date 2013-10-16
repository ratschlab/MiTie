#include "assert.h"
#include "region.h"
#include "genome.h"
#include <vector>
	using std::vector;
#include <fstream>
#include <string.h>
	using std::string;

/** default constructor*/
Region::Region()
{
	start = -1; 
	stop = -1;
	strand = '\0';
	chr_num = -1;
	chr = NULL;
	seq = NULL;
	gio = NULL;
}


/** constructor*/
Region::Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = pchr_num;
	chr = NULL;
	seq = NULL;
	fd_out = stdout;

	// initialize genome information object
	gio = new Genome(); 
	int ret = gio->init_genome((char*) gio_fname);
	if (ret<0)
	{   
		fprintf(stderr, "error reading genome info file: %s\n", gio_fname);
		return;
	}
}
/** constructor*/
Region::Region(int pstart, int pstop, int pchr_num, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = pchr_num;
	chr = NULL;
	seq = NULL;
	fd_out = stdout;
	gio = NULL;
}
/** constructor*/
Region::Region(int pstart, int pstop, char* pchr, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = -1;
	chr = new char[strlen(pchr)+1];
	strcpy(chr, pchr);
	seq = NULL;
	fd_out = stdout;
	gio = NULL;
}
/** constructor*/
Region::Region(Region* reg)
{
	start = reg->start; 
	stop = reg->stop;
	strand = reg->strand;
	chr_num = reg->chr_num;
	if (reg->chr)
	{
		chr = new char[strlen(reg->chr)+1];
		strcpy(chr, reg->chr);
	}
	else
	{
		chr = NULL;
	}
	seq = NULL;
	fd_out = reg->fd_out;
	gio = NULL;
	segments = reg->segments;
	transcripts = reg->transcripts;
	transcript_names = reg->transcript_names;
	transcript_paths = reg->transcript_paths;
}

Region::~Region()
{
	delete[] seq;	
	delete[] chr;
	intron_list.clear();
	unique_introns.clear();
}

void Region::set_gio(Genome* pgio)
{
	gio = pgio;
}

int Region::check_TIS(int i)
{
	assert(i<transcripts.size());
	assert(i>=0);

	int pos=0;
	for (unsigned int j=0; j<transcripts[i].size(); j++)
	{
		int flag = transcripts[i][j].flag;
		if (flag==4)//CDS
		{
			if (strand=='+')
			{
				// tis is the start of the first coding exon
				pos = transcripts[i][j].first-start;
				break;
			}
			else
			{
				// tis is the end of the last coding exon
				pos = transcripts[i][j].second-start;
			}
		}
	}
	if (pos==0)
		return -1;

	char* seq = get_sequence();
	int seq_len = strlen(seq);
	if (pos<0)
	{
		printf("check_TIS: Warning: found pos: %i<0\n", pos);
		return -2;
	}
	assert(pos<seq_len);


	if (strand=='+')
	{
		//printf("TIS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-3], seq[pos-2], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2], seq[pos+3], seq[pos+4], seq[pos+5]);
		return (seq[pos]=='A' || seq[pos]=='a') && (seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='G' || seq[pos+2]=='g');
	}
	else
	{
		//printf("TIS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-5], seq[pos-4], seq[pos-3], seq[pos-2], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2], seq[pos+3]);
		return (seq[pos]=='T' || seq[pos]=='t') && (seq[pos-1]=='A' || seq[pos-1]=='a') && (seq[pos-2]=='C' || seq[pos-2]=='c');
	}
}

int Region::check_TTS(int i)
{
	assert(i<transcripts.size());
	assert(i>=0);

	int pos=0;
	for (unsigned int j=0; j<transcripts[i].size(); j++)
	{
		int flag = transcripts[i][j].flag;
		if (flag==4)//CDS
		{
			if (strand=='+')
			{
				// tis is the start of the first coding exon
				pos = transcripts[i][j].second-start;
			}
			else
			{
				// tis is the end of the last coding exon
				pos = transcripts[i][j].first-start;
				break;
			}
		}
	}
	if (pos==0)
		return -1;

	char* seq = get_sequence();
	int seq_len = strlen(seq);
	if (pos<0)
	{
		printf("check_TTS: Warning: found pos: %i<0\n", pos);
		return -2;
	}

	assert(pos<seq_len);


	if (strand=='+')
	{
		//printf("TTS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-3], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2], seq[pos+3], seq[pos+4], seq[pos+5], seq[pos+6]);

		bool found = false;
		found = found || ((seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='A' || seq[pos+2]=='a') && (seq[pos+3]=='A' || seq[pos+3]=='a'));
		found = found || ((seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='G' || seq[pos+2]=='g') && (seq[pos+3]=='A' || seq[pos+3]=='a'));
		found = found || ((seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='A' || seq[pos+2]=='a') && (seq[pos+3]=='G' || seq[pos+3]=='g'));
		return found;
	}
	else
	{
		//printf("TTS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-6], seq[pos-5], seq[pos-4], seq[pos-3], seq[pos-2], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2]);
		bool found = false;
		found = found || ((seq[pos-3]=='T' || seq[pos-3]=='t') && (seq[pos-2]=='T' || seq[pos-2]=='t') && (seq[pos-1]=='A' || seq[pos-1]=='a'));
		found = found || ((seq[pos-3]=='T' || seq[pos-3]=='t') && (seq[pos-2]=='C' || seq[pos-2]=='c') && (seq[pos-1]=='A' || seq[pos-1]=='a'));
		found = found || ((seq[pos-3]=='C' || seq[pos-3]=='c') && (seq[pos-2]=='T' || seq[pos-2]=='t') && (seq[pos-1]=='A' || seq[pos-1]=='a'));
		return found;

	}
	return 0;
}

void Region::load_genomic_sequence()
{
	if (!check_region())
	{
		printf("load_genomic_sequence: check_region failed\n");
		print(stderr);
		exit(-1);
	}
	seq = gio->read_flat_file(chr_num, start-1, stop);
}

void Region::print(_IO_FILE*& fd)
{
	fprintf(fd, "region %s\n", get_region_str());
	fprintf(fd, "region start:\t%i\n", start);
	fprintf(fd, "region stop:\t%i\n", stop);
	fprintf(fd, "region strand:\t%c\n", strand);
	if (chr_num>0)
		fprintf(fd, "region chr_num:\t%i\n", chr_num);
	if (gio && chr_num>0)
		fprintf(fd, "region chr:\t%s\n", gio->get_contig_name(chr_num));
	else if (chr)
		fprintf(fd, "region chr:\t%s\n", chr);
}

char* Region::get_region_str()
{
	char* reg_str = new char[1000];
	if (chr)
		sprintf(reg_str, "%s:%i-%i", chr, start, stop);
	else if (gio && chr_num>0)
		sprintf(reg_str, "%s:%i-%i", gio->get_contig_name(chr_num), start, stop);
	else
	{
		fprintf(stderr, "genome information object not set");
		exit(-1);
	}
	return reg_str;
}



