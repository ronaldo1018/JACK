#include "jack.h"
/* ********************************************************************** */
int		NE;	/* number of Energy values */
TBL_REAL	*ENRG;	/* energies */
int		NQ;	/* number of quantiles */
TBL_REAL	*PERC;	/* energies */
int		NM;	/* number of materials */
TBL_REAL	*FLUC;	/* table values */
/* ********************************************************************** */
#define	    INDEX(m,e,p)    ((m)*NE*NQ + (e)*NQ + (p))
/* #define	    STANDALONE */
/* ********************************************************************** */
void initialize_fluctuations() {
    int	    i, j, k, idx;
    FILE    *F;
    char    fname[MAXLIN];
    sprintf(fname, "%s/fluctuations.tbl", TPATH);
    F = fopen(fname, "r");
    if (! F) {
	fprintf(stderr, "Unable to open %s!", fname);
	return;
    }
    fscanf(F, "%d", &NE);
    ENRG = (TBL_REAL *) malloc(NE*sizeof(TBL_REAL));
    for (i=0; i<NE; i++) {
#if USE_FLOAT
	fscanf(F, "%f", &(ENRG[i]));
#else
	fscanf(F, "%lf", &(ENRG[i]));
#endif
    }
    fscanf(F, "%d", &NQ);
    PERC = (TBL_REAL *) malloc(NQ*sizeof(TBL_REAL));
    for (i=0; i<NQ; i++) {
#if USE_FLOAT
	fscanf(F, "%f", &(PERC[i]));
#else
	fscanf(F, "%lf", &(PERC[i]));
#endif
    }
    fscanf(F, "%d", &NM);
    FLUC = (TBL_REAL *) malloc(NE*NQ*NM*sizeof(TBL_REAL));
    for (i=0; i<NM; i++) {		/* material */
	for (j=0; j<NE; j++) {		/* energy */
/* these should go away by simply not printing them into the file... */
	    fscanf(F, "%d", &k);
	    fscanf(F, "%d", &k);
	    for (k=0; k<NQ; k++) {	/* percentile */
		idx = INDEX(i,j,k);
#if USE_FLOAT
		fscanf(F, "%f", &(FLUC[idx]));
#else
		fscanf(F, "%lf", &(FLUC[idx]));
#endif
	    }
	}
    }
}
/* ********************************************************************** */
static TBL_REAL interpol8(TBL_REAL xlo, TBL_REAL xhi, TBL_REAL x,
			  TBL_REAL ylo, TBL_REAL yhi,
			  TBL_REAL zlo, TBL_REAL zhi) {
    TBL_REAL	alpha, y, z;
    alpha = (x - xlo) / (xhi - xlo);
    y = ylo + alpha*(yhi - ylo);
    z = zlo + alpha*(zhi - zlo);
/* average */
/*    return(0.5*(y+z)); */
/* uniform distribution */
/*    return(y + drand48()*(z-y)); */
/* equal area */
    return((y+z) - COPYSIGN(SQRT(0.5*(y*y+z*z)), y));
}
/* ********************************************************************** */
TBL_REAL get_fluct(TBL_REAL E, enum material_type M, unsigned int *seed) {
    int		i, j, k;
    TBL_REAL    z, s;
/* We want to find the right energy level */
    if (E < ENRG[0]) return(0.0);
    k = -1;
    for (i=0; i<NE; i++) {
	if (E < ENRG[i]) {
	    k = i;
	    break;
	}
    }
    assert(k > 0);
/* generate a random number from zero to 1 */
    z = (TBL_REAL) fastrand01(seed);
    j = NQ-1;
    for (i=0; i<NQ; i++) {
	if (z < PERC[i]) {
	    j = i;
	    break;
	}
    }
/* so now we want to interpolate */
    s = interpol8(ENRG[k-1], ENRG[k], E,
		  FLUC[INDEX(M,k-1,j-1)], FLUC[INDEX(M,k,j-1)],
		  FLUC[INDEX(M,k-1,j)], FLUC[INDEX(M,k,j)]);
/* the value is hyperbolic sin of 1000 times the value in MeV so... */
    return(0.001*SINH(s));
}
/* ********************************************************************** */
#ifdef	STANDALONE
void PrintFluctuations() {
    int i, j, k, idx;
    printf("%d\n", NE);
    for (i=0; i<NE; i++) {
	printf("%f ", ENRG[i]);
    }
    putchar('\n');

    printf("%d\n", NQ);
    for (i=0; i<NQ; i++) {
	printf("%f ", PERC[i]);
    }
    putchar('\n');

    printf("%d\n", NM);
    for (i=0; i<NM; i++) {		/* material */
	for (j=0; j<NE; j++) {		/* energy */
	    printf("%d %d ", i, j);
	    for (k=0; k<NQ; k++) {	/* percentile */
		idx = INDEX(i,j,k);
		printf("%.5f ", FLUC[idx]);
	    }
	    putchar('\n');
	}
    }
}
/* ********************************************************************** */
int main(int argc, char **argv) {
    int		i;
    TBL_REAL	e;
    InitFluctuations("mluct4.out");
    e = 100.0;
    for (i=0; i<100; i++) {
	printf("%f\n", get_fluct(e, WATER));
    }
//    PrintFluctuations();
}
#endif
