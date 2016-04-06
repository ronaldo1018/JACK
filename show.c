#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

int main(int argc, char **argv) {
    int		    i, x, y, z, nrow, ncol, npln, ntot;
    char	    *fc, *temp;
    FILE	    *F;
    float	    *A, max, v;
    if (argc < 2) {
	fprintf(stderr, "expect: <input file>\n");
	fprintf(stderr, "input files have a name like <name>_<row>_<col>_<plane>.raw\n");
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
/* get the maximum value */
    max = 0.0;
    for (i=0; i<ntot; i++) {
	if (A[i] > max) max = A[i];
    }
    if (max == 0.0) {
	printf("file is all zeros!\n");
	return(0);
    }
/* print it out */
    for (z=0; z<npln; z++) {
	printf("plane %d/%d\n", z, npln);
	for (y=0; y<ncol; y++) {
	    for (x=0; x<nrow; x++) {
		v = A[x + ncol*(y + nrow*z)];
		printf("%1d", (int)(9.0*v/max));
	    }
	    printf("\n");
	}
    }
    return(0);
}
