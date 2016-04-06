#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static void error(char *s) {
    fprintf(stderr, "error: %s\n", s);
    exit(-1);
}

int main(int argc, char **argv) {
    int		    i, nrow, ncol, npln, ntot;
    char	    *fc, *temp;
    FILE	    *F;
    float	    *A, *B, max;
    if (argc < 4) {
	fprintf(stderr, "expect: <input file A> <input file B> <output file>\n");
	fprintf(stderr, "input files have a name like <name>_<row>_<col>_<plain>\n");
	exit(-1);
    }
/* process first intput file */
    fc = strdup(argv[1]);
    F = fopen(fc,"rb");
    if (! F) {
	fprintf(stderr, "unable to open: %s\n", fc);
	exit(-1);
    }
    temp = strtok(fc,"_.");
    temp = strtok(0,"_."); nrow = atoi(temp);
    temp = strtok(0,"_."); ncol = atoi(temp);
    temp = strtok(0,"_."); npln = atoi(temp);
    fprintf(stderr, "assuming file %s is %d x %d x %d\n", argv[1], nrow, ncol, npln);
    ntot = nrow*ncol*npln;
    A = malloc(ntot*sizeof(float));
    i = (int) fread(A, sizeof(float), ntot, F);
    if (i != ntot) {
	fprintf(stderr, "found only %d floats in file %s\n", i, argv[1]);
	exit(-1);
    }
    fclose(F);
/* process first intput file */
    fc = strdup(argv[2]);
    F = fopen(fc,"rb");
    if (! F) {
	fprintf(stderr, "unable to open: %s\n", fc);
	exit(-1);
    }
    temp = strtok(fc,"_.");
    temp = strtok(0,"_."); if (atoi(temp) != nrow) error("number of rows does not match");
    temp = strtok(0,"_."); if (atoi(temp) != ncol) error("number of columns does not match");
    temp = strtok(0,"_."); if (atoi(temp) != npln) error("number of planes does not match");
    fprintf(stderr, "assuming file %s is %d x %d x %d\n", argv[2], nrow, ncol, npln);
    B = malloc(ntot*sizeof(float));
    i = (int) fread(B, sizeof(float), ntot, F);
    if (i != ntot) {
	fprintf(stderr, "found only %d floats in file %s\n", i, argv[2]);
	exit(-1);
    }
    fclose(F);
/* add them */
    max = 0;
    for (i=0; i<ntot; i++) {
	A[i] += B[i];
	if (A[i] > max) max = A[i];
    }

/* write output file */
    fc = argv[3];
    F = fopen(fc,"wb");
    if (! F) {
	fprintf(stderr, "unable to open: %s\n", fc);
	exit(-1);
    }
    fwrite(A, sizeof(float), ntot, F);
    fclose(F);
    return(0);
}
