#include    "jack.h"

table *read_table_from_file(char *fname) {
    FILE		    *F;
    table		    *T;
    int			    i, size, np;
    REAL		    val, min, max;
    char		    line[MAXLIN], data[MAXLIN];
    enum material_type	    M;
    enum particle_type	    P;
    F = fopen(fname, "r");
    if (! F) {
	    return(NULL);
    }
/* read initial header */
    while (fgets(line, MAXLIN, F) && line[0] == '#') { };
/* read the number of point, minimun and maximum */
#if USE_FLOAT
    if (sscanf(line, "%d %f %f", &np, &min, &max) != 3) {
#else
    if (sscanf(line, "%d %lf %lf", &np, &min, &max) != 3) {
#endif
	fprintf(stderr, "expected: NP MIN MAX on\n\t\"%s\"\n", line);
	return(NULL);
    }
    if (np < 2 || max < min) {
	fprintf(stderr, "illogical: NP MIN MAX on\n\t\"%s\"\n", line);
	return(NULL);
    }
/* read and decode the material */
    if (! fgets(line, MAXLIN, F)) {
	fprintf(stderr, "ran out of data from %s\n", fname);
	return(NULL);
    }
    sscanf(line, "%s", data);
    M = UNKNOWN_M;
    for (i=0; i<NMATERIAL; i++) {
	if (! strcmp(data, material_name[i])) M = (enum material_type) i;
    }
    if (M == UNKNOWN_M) {
	fprintf(stderr, "unknown material %s\n", data);
	return(NULL);
    }
/* read and decode the particle */
    if (! fgets(line, MAXLIN, F)) {
	fprintf(stderr, "ran out of data from %s\n", fname);
	return(NULL);
    }
    sscanf(line, "%s", data);
    P = UNKNOWN_P;
    for (i=0; i<NPARTICLE; i++) {
	if (! strcmp(data, particle_name[i])) P = (enum particle_type) i;
    }
    if (P == UNKNOWN_P) {
	fprintf(stderr, "unknown particle %s\n", data);
	return(NULL);
    }
/* allocate */	    
    size = sizeof(table) + (np-1)*sizeof(TBL_REAL);
    T = (table *) malloc(size);
    T->M = M;
    T->P = P;
    T->np = np;
    T->mysize = size;
    T->min = (TBL_REAL) min;
    T->max = (TBL_REAL) max;
    T->range = (TBL_REAL) max - min;
    T->step = (TBL_REAL) (max - min)/(np-1);
    for (i=0; i<np; i++) {
	if (! fgets(line, MAXLIN, F)) {
	    fprintf(stderr, "ran out of data at %d\n", i);
	    free((char *) T);
	    return(NULL);
	}
#if USE_FLOAT
	if (sscanf(line, "%f", &val) != 1) {
#else
	if (sscanf(line, "%lf", &val) != 1) {
#endif
	    fprintf(stderr, "expected: value on\n\t\"%s\"\n", line);
	    free((char *) T);
	    return(NULL);
	}
	T->YV[i] = (TBL_REAL) val;
    }
    fclose(F);
#if TBL_DEBUG
    printf("read table from %s with %d points over %g to %g\n", fname, np, min, max);
#endif
    return(T);
}


void write_table_to_file(table *T, char *fname) {
    FILE    *F;
    int	    i;
    F = fopen(fname, "w");
    if (! F) {
	fprintf(stderr, "cannot write to \"%s\"\n", fname);
	return;
    }
    fprintf(F, "%d %g %g\n", T->np, (REAL) T->min, (REAL) T->max);
    for (i=0; i<T->np; i++) {
	fprintf(F, "%g\n", (REAL) T->YV[i]);
    }
}

/* this function needs a careful look for the case when x is not real */
TBL_REAL interpolate_table(table *T, TBL_REAL x) {
    int		i;
    TBL_REAL	y, z;
    if (x <= T->min) {
	return(T->YV[0]);
    }
    if (x >= T->max) {
	return(T->YV[T->np-1]);
    }
    i = (int)((x - T->min)/T->step);
#ifdef PARANOID
    assert(i == (int) floor((x - T->min)/T->step));
#endif
    assert(i >= 0 && i <= T->np-1);
    //what is the meaning of step here??
    z = (x - i*T->step)/T->step;
    y = T->YV[i] + z * (T->YV[i+1] - T->YV[i]);
#if TBL_DEBUG
    printf("interpolate %g -> %d + %g -> %g %g -> %g\n", x, i, z, T->YV[i], T->YV[i+1], y);
#endif
    return(y);
}
