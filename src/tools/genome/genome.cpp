#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>


//#include "common.h"
#include "genome.h"


Genome::Genome()
{
	num_contigs=0 ;
	
	contig_fnames=NULL ;
	contig_fasta_fnames=NULL ;
	contig_names=NULL ;
	alph= (char*) "acgt" ;
	contig_files=NULL ;
	
	num_est_files=0 ;
	est_fnames=NULL ;
	
	num_prot_files=0 ;
	prot_fnames=NULL ;
	
	num_cdna_files=0 ;
	cdna_fnames=NULL ;
	
	num_annotation_files=0 ;
	annotation_fnames=NULL ;
	
	basedir=NULL ;
	
	//genome defaults
	MAX_DOUBLE_QGAP    = 5 ; 
	MAX_DOUBLE_TGAP    = 5 ;
	MAX_CDNA_DELETION  = 3 ;
	MAX_CDNA_INSERTION = 3 ;
	MIN_INTRON_LEN    = 20 ;
	MAX_INTRON_LEN    = 50000 ;
	TERMINAL_EXON_END_TOL = 3 ;
	MERGE_EST_TRANSCRIPTS = 1 ;
	MAX_PERFECT_EST_MATCH_LEN = 100000 ;
	BLAT_BEST_HIT_ONLY = 0 ;
	MIN_SSQUALITY_SCORE = 100 ;
	MIN_EXONQUALITY_SCORE = 95 ;
	
	BLAT_BEST_HIT_MARGIN = 0.05 ;
	MIN_EST_COVER_FRAC = 0.95 ;
	MIN_PROT_COVER_FRAC = 0.95 ;
	MIN_CDNA_COVER_FRAC = 0.99 ;

}

Genome::~Genome()
{
	free(basedir) ;
	int i;
	for (i=0; i<num_contigs; i++)
	{
		free(contig_fnames[i]) ;
		free(contig_fasta_fnames[i]) ;
		free(contig_names[i]) ;
	}
	free(contig_fnames) ;
	free(contig_fasta_fnames) ;
	free(contig_names) ;

	for (i=0; i<num_contigs; i++)
		if (contig_files[i])
			fclose(contig_files[i]) ;
	free(contig_files) ;
	
	for (i=0; i<num_est_files; i++)
		free(est_fnames[i]) ;
	free(est_fnames) ;

	for (i=0; i<num_prot_files; i++)
		free(prot_fnames[i]) ;
	free(prot_fnames) ;

	for (i=0; i<num_cdna_files; i++)
		free(cdna_fnames[i]) ;
	free(cdna_fnames) ;
	cdna_fnames=NULL ;
	num_cdna_files=0 ;

	for (i=0; i<num_annotation_files; i++)
		free(annotation_fnames[i]) ;
	free(annotation_fnames) ;

	free(alph) ;
}

char* Genome::read_flat_file(int chr_num)
{
	if (chr_num>=num_contigs)
	{
		fprintf(stderr, "contig num exceeds number of contigs\n") ;
		exit(-1);
	}

	int filesize = this->contig_len(chr_num);
	//fprintf(stdout, "contig_fnames[%i]: %s %i\n",chr_num, this->contig_fnames[chr_num], filesize);	
	char* flat_fname = this->contig_fnames[chr_num];

	FILE* file = fopen(flat_fname, "r"); 
	char* str = new char[filesize+1];
	if (!str)
    {
        fprintf(stderr, "Error allocating mem: %i\n", filesize+1);
        exit(-1);
    }

	return this->read_line(file, str, filesize);
}
char* Genome::read_flat_file(int chr_num, int start, int stop)
{
	if (start>stop)
	{
		fprintf(stderr, "Genome::read_flat_file: bad interval [%i, %i]\n", start, stop);
		exit(-1);
	}
	if (chr_num>=num_contigs)
	{
		fprintf(stderr, "contig num exceeds number of contigs\n") ;
		exit(-1);
	}
	int filesize = this->contig_len(chr_num);

	if (start>=filesize || start<0)
	{
		fprintf(stderr, "start index (%i) out of bounds [0, %i]\n", start, filesize-1);
		exit(-1);
	}
	if (stop>=filesize || stop<0)
	{
		fprintf(stderr, "stop index (%i) out of bounds [0, %i]\n", stop, filesize-1);
		exit(-1);
	}
	int interval_size = stop-start+1; 
	//fprintf(stdout, "contig_fnames[%i]: %s %i\n",chr_num, this->contig_fnames[chr_num], filesize);	
	char* flat_fname = this->contig_fnames[chr_num];

	FILE* file = fopen(flat_fname, "r"); 
	char* str = new char[interval_size+1];
	if (!str)
	{
		fprintf(stderr, "Error allocating mem: %i\n", interval_size+1);
		exit(-1);
	}
	fseek(file, start, SEEK_SET);
	int num_read = fread(str, sizeof(char), interval_size, file);
	if (num_read!=interval_size)
	{
		fprintf(stderr, "Error; read %i bytes, interval_size: %i\n", num_read, interval_size);
		exit(-1);
	}
	fclose(file);
	return str;
}
int Genome::contig_len(int chr_num)
{
	if (chr_num>=num_contigs)
	{
		fprintf(stderr, "contig num exceeds number of contigs\n") ;
		exit(-1);
	}

	char* flat_fname = this->contig_fnames[chr_num];

	struct stat filestatus;
	stat( flat_fname, &filestatus );
	return filestatus.st_size;
}

/* read a line from a file */
char* Genome::read_line(FILE*infile, char* line, size_t max_len) 
{
	size_t i;
	int c ;
	
	i=0 ;
	while (!feof(infile))
    {
		c=fgetc(infile) ;
		if ((c==EOF) || (c==13) || (c==10))
			break ;
		line[i]=c ;
		
		i++ ;
		if (i>max_len)
		{
			fprintf(stderr, "line buffer not large enough\n") ;
			exit(-1) ;
		} ;
    } ;
	line[i]=0 ;
	
	if (i==0)
		strcpy(line,"") ;
	return line ;
}

char* Genome::path_with_basedir(char* b1, char* basedir)
{
	char buf[1000] ;
	if (b1[0]!='/')
	{
		if (basedir[strlen(basedir)-1]=='/')
			sprintf(buf, "%s%s", basedir, b1);
		else
			sprintf(buf, "%s/%s", basedir, b1);
		return strdup(buf) ;
	} else
		return strdup(b1) ;
}



int Genome::init_genome(char *fname)
{
	if (fname==NULL)
		fname= (char*) "genome.conf" ;
	//fprintf(stderr, "%s\n", fname) ;
	
	FILE * fd=fopen(fname, "r") ;
	if (fd==NULL)
		return -1 ;
	int contigs_found = 0 ;
	int ests_found = 0 ;
	int cdnas_found = 0 ;
	int alph_found=0 ;
	int basedir_found=0 ;
	int annotation_found=0 ;
	basedir = (char*) malloc(1000*sizeof(char)) ;
   
	while (!feof(fd))
	{
		char line[1000];
		char *ret=read_line(fd, line, 1000) ;
		if (ret==NULL)
			break ;
		//fprintf(stdout, "line:\t %s\n", line) ;
		if (line[0]=='#')
			continue ;
		
		if (strncmp(line,"BASEDIR",7)==0)
		{
			basedir_found = 1 ;
			int args=sscanf(line, "BASEDIR %s", basedir) ;
			if (args!=1)
			{
				fprintf(stderr, "string argument expected after BASEDIR\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MIN_EST_COVER_FRAC", 18)==0)
		{
			int args=sscanf(line, "MIN_EST_COVER_FRAC %lf", &MIN_EST_COVER_FRAC) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MIN_EST_COVER_FRAC\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MIN_PROT_COVER_FRAC", 19)==0)
		{
			int args=sscanf(line, "MIN_PROT_COVER_FRAC %lf", &MIN_PROT_COVER_FRAC) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MIN_PROT_COVER_FRAC\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MIN_CDNA_COVER_FRAC", 19)==0)
		{
			int args=sscanf(line, "MIN_CDNA_COVER_FRAC %lf", &MIN_CDNA_COVER_FRAC) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MIN_CDNA_COVER_FRAC\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MAX_DOUBLE_QGAP", strlen("MAX_DOUBLE_QGAP"))==0)
		{
			int args=sscanf(line, "MAX_DOUBLE_QGAP %i", &MAX_DOUBLE_QGAP) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MAX_DOUBLE_QGAP\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MAX_DOUBLE_TGAP", strlen("MAX_DOUBLE_TGAP"))==0)
		{
			int args=sscanf(line, "MAX_DOUBLE_TGAP %i", &MAX_DOUBLE_TGAP) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MAX_DOUBLE_TGAP\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MAX_CDNA_DELETION", strlen("MAX_CDNA_DELETION"))==0)
		{
			int args=sscanf(line, "MAX_CDNA_DELETION %i", &MAX_CDNA_DELETION) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MAX_CDNA_DELETION\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MAX_CDNA_INSERTION", strlen("MAX_CDNA_INSERTION"))==0)
		{
			int args=sscanf(line, "MAX_CDNA_INSERTION %i", &MAX_CDNA_INSERTION) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MAX_CDNA_INSERTION\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MIN_INTRON_LEN", strlen("MIN_INTRON_LEN"))==0)
		{
			int args=sscanf(line, "MIN_INTRON_LEN %i", &MIN_INTRON_LEN) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after MIN_INTRON_LEN\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MAX_INTRON_LEN", 14)==0)
		{
			int args=sscanf(line, "MAX_INTRON_LEN %i", &MAX_INTRON_LEN) ;
			if (args!=1)
			{
				fprintf(stderr, "integer argument expected after MAX_INTRON_LEN\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"TERMINAL_EXON_END_TOL", strlen("TERMINAL_EXON_END_TOL"))==0)
		{
			int args=sscanf(line, "TERMINAL_EXON_END_TOL %i", &TERMINAL_EXON_END_TOL) ;
			if (args!=1)
			{
				fprintf(stderr, "real argument expected after TERMINAL_EXON_END_TOL\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MERGE_EST_TRANSCRIPTS", strlen("MERGE_EST_TRANSCRIPTS"))==0)
		{
			int args=sscanf(line, "MERGE_EST_TRANSCRIPTS %i", &MERGE_EST_TRANSCRIPTS) ;
			if (args!=1 || ((MERGE_EST_TRANSCRIPTS!=0) && (MERGE_EST_TRANSCRIPTS!=1)))
			{
				fprintf(stderr, "0 or 1 expected after MERGE_EST_TRANSCRIPTS\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MAX_PERFECT_EST_MATCH_LEN", strlen("MAX_PERFECT_EST_MATCH_LEN"))==0)
		{
			int args=sscanf(line, "MAX_PERFECT_EST_MATCH_LEN %i", &MAX_PERFECT_EST_MATCH_LEN) ;
			if (args!=1)
			{
				fprintf(stderr, "integer expected after MAX_PERFECT_EST_MATCH_LEN\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"BLAT_BEST_HIT_ONLY", strlen("BLAT_BEST_HIT_ONLY"))==0)
		{
			int args=sscanf(line, "BLAT_BEST_HIT_ONLY %i", &BLAT_BEST_HIT_ONLY) ;
			if (args!=1 || (BLAT_BEST_HIT_ONLY!=0 && BLAT_BEST_HIT_ONLY!=1))
			{
				fprintf(stderr, "integer expected after BLAT_BEST_HIT_ONLY\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"BLAT_BEST_HIT_MARGIN", strlen("BLAT_BEST_HIT_MARGIN"))==0)
		{
			int args=sscanf(line, "BLAT_BEST_HIT_MARGIN %lf", &BLAT_BEST_HIT_MARGIN) ;
			if (args!=1)
			{
				fprintf(stderr, "integer expected after BLAT_BEST_HIT_MARGIN\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MIN_SSQUALITY_SCORE", strlen("MIN_SSQUALITY_SCORE"))==0)
		{
			int args=sscanf(line, "MIN_SSQUALITY_SCORE %i", &MIN_SSQUALITY_SCORE) ;
			if (args!=1)
			{
				fprintf(stderr, "integer expected after MIN_SSQUALITY_SCORE\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"MIN_EXONQUALITY_SCORE", strlen("MIN_EXONQUALITY_SCORE"))==0)
		{
			int args=sscanf(line, "MIN_EXONQUALITY_SCORE %i", &MIN_EXONQUALITY_SCORE) ;
			if (args!=1)
			{
				fprintf(stderr, "integer expected after MIN_EXONQUALITY_SCORE\n") ;
				fclose(fd) ;
				return -1 ;
			}
		}
		else if (strncmp(line,"CONTIGS",7)==0)
		{
			if (strlen(basedir)==0)
			{
				fprintf(stderr, "BASEDIR line has to appear before all other config lines\n") ;
				fclose(fd) ;
				return -1 ;
			}
			
			contigs_found = 1 ;
			int args=sscanf(line, "CONTIGS %i", &num_contigs) ;
			if (args!=1)
			{
				fprintf(stderr, "integer argument expected after CONTIGS\n") ;
				fclose(fd) ;
				return -1 ;
			}
			//fprintf(stderr, "Using %i contigs\n", num_contigs) ;
			contig_fnames=(char**)malloc(num_contigs*sizeof(char*)) ;
			contig_fasta_fnames=(char**)malloc(num_contigs*sizeof(char*)) ;
			contig_names=(char**)malloc(num_contigs*sizeof(char*)) ;
			contig_files=(FILE**)malloc(num_contigs*sizeof(FILE*)) ;
			
			int i ;
			for (i=0; i<num_contigs; i++)
			{
				contig_fnames[i]=NULL ;
				contig_fasta_fnames[i]=NULL ;
				contig_names[i]=NULL ;
				contig_files[i]=NULL ;
				
				char b1[1000], b2[1000], b3[1000] ;
				char *ret=read_line(fd, line, 1000) ;
				assert(ret!=NULL) ;
				int args=sscanf(line, "%s %s %s", b1, b2, b3) ;
				assert(args==3) ;
				contig_names[i]=strdup(b1) ;
				contig_fnames[i]=path_with_basedir(b2,basedir) ;
				contig_fasta_fnames[i]=path_with_basedir(b3,basedir) ;
				//if (num_contigs<30)
				//fprintf(stderr, "%s in file %s\n", contig_names[i], contig_fnames[i]) ;
			}
		}
		else if (strncmp(line,"ESTFILES",8)==0)
		{
			if (strlen(basedir)==0)
			{
				fprintf(stderr, "BASEDIR line has to appear before all other config lines\n") ;
				fclose(fd) ;
				return -1 ;
			}
			ests_found = 1 ;
			int args=sscanf(line, "ESTFILES %i", &num_est_files) ;
			if (args!=1)
			{
				fprintf(stderr, "integer argument expected after ESTFILES\n") ;
				fclose(fd) ;
				return -1 ;
			}
			est_fnames=(char**)malloc(num_est_files*sizeof(char*)) ;
			
			int i ;
			for (i=0; i<num_est_files; i++)
			{
				est_fnames[i]=NULL ;
				
				char b1[1000] ;
				char *ret=read_line(fd, line, 1000) ;
				assert(ret!=NULL) ;
				int args=sscanf(line, "%s", b1) ;
				assert(args==1) ;
				est_fnames[i]=path_with_basedir(b1,basedir) ;
			}
		}
		else if (strncmp(line,"PROTFILES",9)==0)
		{
			if (strlen(basedir)==0)
			{
				fprintf(stderr, "BASEDIR line has to appear before all other config lines\n") ;
				fclose(fd) ;
				return -1 ;
			}
			int args=sscanf(line, "PROTFILES %i", &num_prot_files) ;
			if (args!=1)
			{
				fprintf(stderr, "integer argument expected after PROTFILES\n") ;
				fclose(fd) ;
				return -1 ;
			}
			prot_fnames=(char**)malloc(num_prot_files*sizeof(char*)) ;
			
			int i ;
			for (i=0; i<num_prot_files; i++)
			{
				prot_fnames[i]=NULL ;
				
				char b1[1000] ;
				char *ret=read_line(fd, line, 1000) ;
				assert(ret!=NULL) ;
				int args=sscanf(line, "%s", b1) ;
				assert(args==1) ;
				prot_fnames[i]=path_with_basedir(b1,basedir) ;
			}
		}
		else if ((strncmp(line,"CDNAFILES",9)==0) || (strncmp(line,"FLCDNAFILES",11)==0))
		{
			if (strlen(basedir)==0)
			{
				fprintf(stderr, "BASEDIR line has to appear before all other config lines\n") ;
				fclose(fd) ;
				return -1 ;
			}
			cdnas_found = 1 ;
			int args = 0 ;
			int old_num_cdna_files = num_cdna_files ;
			if (strncmp(line,"CDNAFILES",9)==0)
				args = sscanf(line, "CDNAFILES %i", &num_cdna_files) ;
			else
				args = sscanf(line, "FLCDNAFILES %i", &num_cdna_files) ;
			if (args!=1)
			{
				fprintf(stderr, "integer argument expected after [FL]CDNAFILES\n") ;
				fclose(fd) ;
				return -1 ;
			}
			if (num_cdna_files>0)
			{
				if (cdna_fnames)
				{
					assert(old_num_cdna_files!=0) ;
					cdna_fnames = (char**)realloc(cdna_fnames, (num_cdna_files+old_num_cdna_files)*sizeof(char*)) ;
				}
				else
					cdna_fnames = (char**)malloc(num_cdna_files*sizeof(char*)) ;
				
				int i ;
				for (i=0; i<num_cdna_files; i++)
				{
					cdna_fnames[i+old_num_cdna_files]=NULL ;
					
					char b1[1000] ;
					char *ret=read_line(fd, line, 1000) ;
					assert(ret!=NULL) ;
					int args=sscanf(line, "%s", b1) ;
					assert(args==1) ;
					cdna_fnames[i+old_num_cdna_files] = path_with_basedir(b1,basedir) ;
				}
				num_cdna_files += old_num_cdna_files ;
			}
		}
		else if (strncmp(line,"ANNOTATIONFILES",15)==0)
		{
			if (strlen(basedir)==0)
			{
				fprintf(stderr, "BASEDIR line has to appear before all other config lines\n") ;
				fclose(fd) ;
				return -1 ;
			}
			annotation_found = 1 ;
			int args=sscanf(line, "ANNOTATIONFILES %i", &num_annotation_files) ;
			if (args!=1)
			{
				fprintf(stderr, "integer argument expected after ANNOTATIONFILES\n") ;
				fclose(fd) ;
				return -1 ;
			}
			annotation_fnames=(char**)malloc(num_annotation_files*sizeof(char*)) ;
			
			int i ;
			for (i=0; i<num_annotation_files; i++)
			{
				annotation_fnames[i]=NULL ;
				
				char b1[1000] ;
				char *ret=read_line(fd, line, 1000) ;
				assert(ret!=NULL) ;
				int args=sscanf(line, "%s", b1) ;
				assert(args==1) ;
				annotation_fnames[i]=path_with_basedir(b1,basedir) ;
			}
		}
		else if (strncmp(line,"ALPHABET",8)==0)
		{
			alph_found=1 ;
			char buf[1000] ;
			int args=sscanf(line, "ALPHABET %s", buf) ;
			assert(args==1) ;
			assert(strlen(buf)==4) ;
			alph=strdup(buf) ;
//			fprintf(stderr, "using alphabet '%s'\n", alph) ;
		}
		else if (strlen(line)>0)
		{
			fprintf(stderr, "unrecognized line '%s'\n", line) ;
			fclose(fd) ;
			return -1 ;
		}
	}
	
	assert(alph_found==1) ;
	assert(contigs_found==1) ;
	assert(ests_found==1) ;
	assert(cdnas_found==1) ;
	assert(basedir_found==1) ;
	assert(annotation_found==1) ;

	fclose(fd) ;
	return 0 ;
}

char* Genome::get_contig_name(int chr_num)
{
	if (chr_num<num_contigs && chr_num>=0)
		return contig_names[chr_num];
	
	fprintf(stderr, "chr_num %d out of range (num_contigs: %d)", chr_num, num_contigs);
	return NULL;
}


/* returns the contig idx from the string */
int Genome::get_contig_num(char *contig_name)
{
	int i, contig_num=-1 ;
	
	for (i=0; i<num_contigs; i++)
		if (strcmp(contig_name, contig_names[i])==0)
			contig_num=i ;
	if (contig_num==-1)
    {
		fprintf(stderr, "contig name not found (%s)\n", contig_name) ;
		exit(-1) ;
    } 
	return contig_num ;
} 
void Genome::complement(char* str, const int len)
{   
    int i;
    for (i=0; i<len-1; i++)
    {   
        switch (str[i])
        {   
            case 'A': str[i] = 'T'; break;
            case 'a': str[i] = 't'; break;
            case 'C': str[i] = 'G'; break;
            case 'c': str[i] = 'g'; break;
            case 'G': str[i] = 'C'; break;
            case 'g': str[i] = 'c'; break;
            case 'T': str[i] = 'A'; break;
            case 't': str[i] = 'a'; break;
            case 'n': str[i] = 'n'; break;
            case 'N': str[i] = 'N'; break;
            default : str[i] = str[i]; break;
        }

    }
    return;
}

