
#ifndef _REGION_H__
#define _REGION_H__

//#define READ_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "genome.h"
#include <utility>
	using std::pair;
#include <vector>
	using std::vector;
#include <fstream>
#include <string>
	using std::string;

//typedef pair<int, int> segment;
class segment{
	public:
		segment(){first = 0; second = 0; flag = -1;};
		segment(int f, int s){first = f; second = s; flag=-1;};
		segment(int f, int s, int fl){first = f; second = s; flag=fl;};
		int first;
		int second;
		int flag;//4:CDS, 3: 3'UTR 5:5'UTR

		friend bool operator== (const segment& lhs, const segment& rhs)
		{return lhs.first==rhs.first && lhs.second==rhs.second;}

		friend bool operator<  (const segment& lhs, const segment& rhs)
		{return lhs.first<rhs.first || (!(rhs.first<lhs.first) && lhs.second<rhs.second);}

};
#define NO_CONNECTION -2
#define NEIGHBOR -1
#define CONNECTION 0

//#define INNER=0
//#define INITIAL=1
//#define TERMINAL=2

class Region
{
	public:
		int id;
		int start; 
		int stop;
		int chr_num;
		char* chr;
		char strand; 
		vector<segment> intron_list;
		vector<segment> unique_introns;
		vector<segment> segments;
		vector<vector<segment> > transcripts;
		vector<string> transcript_names;
		vector<string> gene_names;
		vector<vector<int> > transcript_paths;
		FILE* fd_out;

		Genome* gio; 
    	
			
		/** default constructor*/	
		Region();

		/** constructor*/
		Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname);
		Region(int pstart, int pstop, int pchr_num, char pstrand);
		Region(int pstart, int pstop, char* chr, char pstrand);
		Region(Region* reg);

		/** destructor*/
		~Region();


		char* get_sequence()
		{
			if (!seq)
				load_genomic_sequence();
			return seq;
		};

		bool check_region()
		{
			if (!(start>=0 && stop<gio->contig_len(chr_num)))
			{
				fprintf(stderr, "chr_num: %i, chr_len: %i\n", chr_num, gio->contig_len(chr_num));
			}
			return (start>=0 && stop<gio->contig_len(chr_num));	
		};


		/** get the triplet codon for pos in transcript i */
		char* get_triplet(int i, int pos, int* offset);

		/** find translation initiation site for transcript i*/
		int get_TIS(int i); 
		/** check translation initiation consensus for transcript i*/
		int check_TIS(int i); 

		/** find translation termination site for transcript i */
		int get_TTS(int i); 
		/** check translation termination consensus for transcript i */
		int check_TTS(int i); 

		void set_gio(Genome* pgio);
	
		void load_genomic_sequence(); 			

		char* get_region_str();

		virtual void print(FILE*& fd);

	private: 
		char* seq;

};

#endif
