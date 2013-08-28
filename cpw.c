/* cpw.c */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>

static int max_length (const char * filename, int window) {
	FILE  *fp;
	char   line[32];
	int    i, loc, max;
	
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "file open error\n");
		exit(1);
	}
		
	max = 0;
	while (fgets(line, sizeof(line), fp) != NULL) {
		if (sscanf(line, "%d", &i) != 1) {
			fprintf(stderr, "file read error\n");
			exit(1);
		}
		loc = (float)i / (float)window;
		if (loc > max) max = loc;
	}
	
	return max +1;
}

static int* read_sequence (const char * filename, int window, int max) {
	FILE  *fp;
	int    i;
	char   line[32];
	double loc;
	int   *seq;
	
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "file open error\n");
		exit(1);
	}
	
	
	seq = malloc(max * sizeof(int));
	for (i = 0; i < max; i++) seq[i] = 0;
	
	while (fgets(line, sizeof(line), fp) != NULL) {
		if (sscanf(line, "%d", &i) != 1) {
			fprintf(stderr, "file read error\n");
			exit(1);
		}
		
		loc = (float)i / (float)window; 
		if (loc >= max) {
			fprintf(stderr, "maximum length exceeded\n");
		}
		seq[(int)loc]++;
	}
	
	return seq;
}


static char * usage = "\
usage: cpw <file> <W>\n\
  file = one line per read position\n\
  W    = window size [int] 0 < W < 1000\n\
";

int main (int argc, char ** argv) {
	int window, slen, i;
	int *seq;

	if (argc != 3) {
		fprintf(stderr, "%s", usage);
		exit(1);
	}
	
	window = atoi(argv[2]);
	
	slen = max_length(argv[1], window);
	seq  = read_sequence(argv[1], window, slen);
	
	for (i = 0; i < slen; i++) {
		printf("%d\n", seq[i]);
	}

	return 0;
}
