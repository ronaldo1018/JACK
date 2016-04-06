#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "coord.h"
#define	OUTSIDE	    0
#define	INSIDE	    1
#define	SURFACE	    2
#define	SQUARE(x)   ((x)*(x))
/* this program reads a cube of specified size, adds an ellipsoid at a
 * specified location with a specified size, and dumps out the combined file */

/* center and radius are global */
static float	    rx, ry, rz, cx, cy, cz;

int on_sphere(int x, int y, int z) {
    float    d;
    d = sqrtf(SQUARE((x-cx)/rx) + SQUARE((y-cy)/ry) + SQUARE((z-cz)/rz));
    if (d < 0.95) {
	return(INSIDE);
    } else if (d < 1.05) {
	return(SURFACE);
    }
    return(OUTSIDE);
}

void print_pixel(FILE *F, short value, int flag, int max) {
    int	    x;
    if (value < 0) {
	fprintf(F, "#000 ");
    } else if (flag > 1) {
	fprintf(F, "#F00 ");
    } else if (flag < -1) {
	fprintf(F, "#00F ");
    } else {
	x = (value*16)/max;
	fprintf(F, "#%01x%01x%01x ", x, x, x);
    }
}

int main(int argc, char **argv) {
    int		    i, j, k, idx, nrow, ncol, npln, ntot, min, max, n;
/* tumor */
    int		    t_xmin, t_xmax, t_ymin, t_ymax, t_zmin, t_zmax;
/* exclusion */
    int		    e_xmin, e_xmax, e_ymin, e_ymax, e_zmin, e_zmax;
    char	    *fc, *temp;
    FILE	    *F;
    short	    *v, *zero;
    if (argc < 6) {
	fprintf(stderr, "expect: <input file> <location file> <0/1 file> <x/y file> <x/z file> <y/z file\n");
	fprintf(stderr, "<input file> with a name like <name>_<row>_<col>_<plain>\n");
	fprintf(stderr, "<location file> has the tumor location\n");
	exit(-1);
    }
/* process intput file */
    fc = argv[1];
    F = fopen(fc,"rb");
    if (! F) {
	fprintf(stderr, "unable to open: %s\n", fc);
	exit(-1);
    }
    temp = strtok(fc,"_.");
    temp = strtok(0,"_."); nrow = atoi(temp);
    temp = strtok(0,"_."); ncol = atoi(temp);
    temp = strtok(0,"_."); npln = atoi(temp);
    fprintf(stderr, "assuming file is %d x %d x %d\n", nrow, ncol, npln);

    ntot = nrow*ncol*npln;
    v = malloc(ntot*sizeof(short));
    zero = malloc(ntot*sizeof(short));
    if (! v) {
	fprintf(stderr, "failed to malloc %d shorts\n", ntot);
	exit(-1);
    }
    i = (int) fread(v, sizeof(short), ntot, F);
    if (i != ntot) {
	fprintf(stderr, "found only %d shorts in file\n", i);
	exit(-1);
    }
    fclose(F);
/* process location file */
    fc = argv[2];
    F = fopen(fc,"r");
    if (! F) {
	fprintf(stderr, "unable to open: %s\n", fc);
	exit(-1);
    }
    if (fscanf(F, "%d %d %d %d %d %d", &t_xmin, &t_xmax, &t_ymin, &t_ymax, &t_zmin, &t_zmax) != 6) {
	fprintf(stderr, "could not read %s\n", fc);
	exit(-1);
    }
    if (fscanf(F, "%d %d %d %d %d %d", &e_xmin, &e_xmax, &e_ymin, &e_ymax, &e_zmin, &e_zmax) != 6) {
	fprintf(stderr, "could not read %s\n", fc);
	exit(-1);
    }
    fclose(F);
/* calculate center and radius for hyper-ellipse */
    rx = (t_xmax - t_xmin) * 0.5f;
    ry = (t_ymax - t_ymin) * 0.5f;
    rz = (t_zmax - t_zmin) * 0.5f;
    cx = (t_xmax + t_xmin) * 0.5f;
    cy = (t_ymax + t_ymin) * 0.5f;
    cz = (t_zmax + t_zmin) * 0.5f;
    fprintf(stderr, "tumor x [%d - %d] -> c = %f r = %f\n", t_xmin, t_xmax, cx, rx);
    fprintf(stderr, "tumor y [%d - %d] -> c = %f r = %f\n", t_ymin, t_ymax, cy, ry);
    fprintf(stderr, "tumor z [%d - %d] -> c = %f r = %f\n", t_zmin, t_zmax, cz, rz);
    fprintf(stderr, "exclude x [%d - %d]\n", e_xmin, e_xmax);
    fprintf(stderr, "exclude y [%d - %d]\n", e_ymin, e_ymax);
    fprintf(stderr, "exclude z [%d - %d]\n", e_zmin, e_zmax);
/* calculate min/max */
    min = 100000;
    max = -100000;
    for (i=0; i<ntot; i++) {
	zero[i] = 0;
	if (v[i] < 0) continue;
	if (v[i] < min) min = v[i];
	if (v[i] > max) max = v[i];
    }
    fprintf(stderr, "range of input: %d - %d\n", min, max);
/* insert the cancer region */
    n = 0;
    for (i=t_xmin; i<=t_xmax; i++) {
	for (j=t_ymin; j<=t_ymax; j++) {
	    for (k=t_zmin; k<=t_zmax; k++) {
		idx = nrow*ncol*k + ncol*j + i;
		switch (on_sphere(i, j, k)) {
		    case OUTSIDE :
			break;
		    case SURFACE :
			zero[idx] = 2;
			++n;
			break;
		    case INSIDE :
			zero[idx] = 1;
			++n;
			break;
		}
	    }
	}
    }
    fprintf(stderr, "found %d points within tumor\n", n);    
/* insert the exclusion region */
    n = 0;
    for (i=e_xmin; i<=e_xmax; i++) {
	for (j=e_ymin; j<=e_ymax; j++) {
	    for (k=e_zmin; k<=e_zmax; k++) {
		idx = nrow*ncol*k + ncol*j + i;
		if (idx < 0 || idx > ntot) {
		    printf("i=%d j=%d k=%d -> idx=%d > ntot=%d!\n",
			   i, j, k, idx, ntot);
		    continue;
		}
		if (i == e_xmin || i == e_xmax ||
		    j == e_ymin || j == e_ymax ||
		    k == e_zmin || k == e_zmax) {
		    zero[idx] = -2;
		} else {
		    zero[idx] = -1;
		}
		++n;
	    }
	}
    }
    fprintf(stderr, "found %d points in exlusion region \n", n);    
/* write the zero/one location file */
    fc = argv[3];
    F = fopen(fc,"w");
    if (! F) {
	fprintf(stderr, "unable to open: %s\n", fc);
	exit(-1);
    }
    for (i=0; i<npln; i++) {		/* Z */
	for (j=0; j<nrow; j++) {	/* Y */
	    for (k=0; k < ncol; k++) {	/* X */
		idx = i*nrow*ncol + j*nrow + k;
		if (zero[idx]) {
		    fprintf(F, "%d %d %d %d\n", k, j, i, zero[idx]);
		}
	    }
	}
    }
    fclose(F);
    fprintf(stderr, "done writing zero/one file to %s \n", fc);    
/* write the XY file */
    fc = argv[4];
    F = fopen(fc,"w");
    for (i=0; i<npln; i++) {		/* Z */
	for (j=0; j<nrow; j++) {	/* Y */
	    for (k=0; k < ncol; k++) {	/* X */
		idx = i*nrow*ncol + j*nrow + k;
		print_pixel(F, v[idx], zero[idx], max);
	    }
	    fprintf(F, "\n");
	}
    }
    fclose(F);
    fprintf(stderr, "done writing XY file to %s \n", fc);    
/* write the XZ file */
    fc = argv[5];
    F = fopen(fc,"w");
    for (j=0; j<nrow; j++) {		/* Y */
	for (i=0; i<npln; i++) {	/* Z */
	    for (k=0; k < ncol; k++) {	/* X */
		idx = i*nrow*ncol + j*nrow + k;
		print_pixel(F, v[idx], zero[idx], max);
	    }
	    fprintf(F, "\n");
	}
    }
    fclose(F);
    fprintf(stderr, "done writing XZ file to %s \n", fc);    
/* write the YZ file, note pics are rotated by 90 degrees */
    fc = argv[6];
    F = fopen(fc,"w");
    for (k=0; k < ncol; k++) {		/* X */
	for (j=0; j<nrow; j++) {	/* Y */
	    for (i=0; i<npln; i++) {	/* Z */
		idx = i*nrow*ncol + j*nrow + k;
		print_pixel(F, v[idx], zero[idx], max);
	    }
	    fprintf(F, "\n");
	}
    }
    fclose(F);
    fprintf(stderr, "done writing YZ file to %s \n", fc);    
    return(0);
}
