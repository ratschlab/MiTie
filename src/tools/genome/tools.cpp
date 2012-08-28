#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
using namespace std;
#define MAXLINE 1000

void save_score_pos(char* path_prefix, char** score_names, float** scores, int num_scores, int* pos, int num_pos)
{
	
	// open output files
	char filename[MAXLINE]; 
	fstream score_fds[num_scores];
	fstream pos_fd; 
    sprintf(filename, "%s.pos", path_prefix);
    pos_fd.open(filename, fstream::out);
	
    for (int i = 0; i < num_scores; i++) 
	{
        sprintf(filename, "%s.%s", path_prefix, score_names[i]);
        score_fds[i].open(filename, fstream::out);
    }

    fprintf(stdout, "found path prefix: %s\n", path_prefix);
    fprintf(stdout, "writing pos file and %d score files, each %d bytes...\n", num_scores, (int)(num_pos * sizeof(float)));

	// write pos file
	unsigned *pos_buf = (unsigned*) malloc(num_pos*sizeof(unsigned)) ;
    for (int i = 0; i < num_pos; i++) {
        unsigned curr_pos = (unsigned) pos[i];
        if (i && pos[i] < pos[i - 1]) 
		{
            fprintf(stderr, "pos[%d] = %u < %u = pos[%d]", i, pos[i], pos[i - 1], i - 1);
            exit(-1);
        }
		pos_buf[i] = curr_pos ;
    }
    pos_fd.write((char*) pos_buf, (int) sizeof(unsigned)*num_pos);
	pos_fd.close();
    free(pos_buf) ;
	
	// write score files
	float *score_buf = (float*) malloc(num_pos*sizeof(float)) ;
    for (int j = 0; j < num_scores; j++) 
	{
    	for (int i = 0; i < num_pos; i++) 
		{
        	float score = (float) scores[j][i];
        	score_buf[i]=score ;
    	}
	
      	score_fds[j].write((char*) score_buf, (int) sizeof(float)*num_pos);
		score_fds[j].close();
    }
    free(score_buf) ;
}

