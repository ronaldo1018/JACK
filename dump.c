#include    <time.h>
#include    "jack.h"

void dump_good_bad_dose(collector *C) {
#if UNIFORM
    return;
#else
    REAL	good, bad, maxgood, maxbad;
    int		ngood, nbad, idx;
/* if a location file was specified, compute the good/bad dose */
    maxbad = maxgood = AS_REAL(0.0);
    bad = good = AS_REAL(0.0);
    nbad = ngood = 0;
    for (idx=0; idx<C->NTOTAL; idx++) {
	if (ISTUMOR(ARENA->matG[idx])) {
	    good += C->G[idx];
	    if (C->G[idx] > maxgood) maxgood = good;
	    ++ngood;
	} else if (ISEXCLU(ARENA->matG[idx])) {
	    bad += C->G[idx];
	    if (C->G[idx] > maxbad) maxbad = bad;
	    ++nbad;
	}
    }
    if (VERBOSE) {
	printf("Good dose (to tumor %d points)_ = %.2f (max=%.2f)  Bad dose (outside %d points) = %.2f (max=%.2f).\n",
	       ngood, good, maxgood, nbad, bad, maxbad);
    } else {
	printf("%12.5e %12.5e %d %12.5e %12.5e %d\n", good, maxgood, ngood, bad, maxbad, nbad);
    }
#endif
}
	
void dump_dose_text_file(collector *C, char *OUTFIL) {
    FILE	*F;
    REAL	max;
    int		i, ix, iy, iz, ntot, nout;
/* total number of points */
    ntot = C->N[CX] * C->N[CY] * C->N[CZ];
/* dump the dose file */
    F = fopen(OUTFIL, "w");
    if (! F) {
	fprintf(stderr, "unable to open %s\n", OUTFIL);
	return;
    }
/* print the dimensions, only needed if we dump indices */
/*  fprintf(F, "%d %d %d\n", C->N[CX], C->N[CY], C->N[CZ]); */
/* find the maximum */
    max = AS_REAL(0.0);
    for (i=0; i<ntot; i++) {
	if (C->G[i] > max) max = C->G[i];
    }
/* if maximum is zero there is no dose, so return */
    if (max == AS_REAL(0.0)) {
	fclose(F);
	return;
    }
/* dump all locations more than THR % of the maximum  */
    max *= THR;
    nout = 0;
    for (i=0; i<ntot; i++) {
	if (C->G[i] > max) {
	    split(i,ix,iy,iz,C->N[CX],C->N[CY],C->N[CZ]);
	    assert( getidx(ix,iy,iz,C->N[CX],C->N[CY],C->N[CZ]) == i );
	    fprintf(F, "%d %d %d %12.5e\n", ix, iy, iz, C->G[i]);
	    ++nout;
	}
    }
    fclose(F);
/* let me know */
    if (VERBOSE) {
	printf("Wrote dose to text file %s with %d entries.\n", OUTFIL, nout);
    }
}

void dump_dose_binary_file(collector *C, char *OUTFIL) {
    FILE	*F;
    int		ntot;
/* total number of points */
    ntot = C->N[CX] * C->N[CY] * C->N[CZ];
/* open the file */
    F = fopen(OUTFIL, "wb");
    if (! F) {
	fprintf(stderr, "unable to open %s\n", OUTFIL);
	return;
    }
/* dump it out */
    fwrite(C->G, sizeof(ENERGY), ntot, F);
    fclose(F);
/* let me know */
    if (VERBOSE) {
	printf("Wrote dose to binary float file %s with %d floats.\n", OUTFIL, ntot);
    }
}
