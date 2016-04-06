#include    <time.h>
#include    "jack.h"


//#define DEBUG 1
int			nx, ny, nz;

collector *allocate_collector() {
    int		i, size;
    collector	*C;
    assert(ARENA);
    assert(ARENA->NTOTAL);
    size = sizeof(collector) + (ARENA->NTOTAL - 1)*sizeof(ENERGY);
    C = (collector *) malloc(size);
    C->mysize = size;
    C->NTOTAL = ARENA->NTOTAL;
    C->NSTEP = C->NTRACE = C->NMAX = 0;
    C->N[CX] = ARENA->N[CX];
    C->N[CY] = ARENA->N[CY];
    C->N[CZ] = ARENA->N[CZ];
    for (i=0; i<C->NTOTAL; i++) {
	C->G[i] = (ENERGY) AS_REAL(0.0);
    }
    return(C);
}

void reset_collect(collector *C) {
    int	    i;
    C->NSTEP = C->NTRACE = C->NMAX = 0;
    for (i=0; i<C->NTOTAL; i++) {
	C->G[i] = (ENERGY) AS_REAL(0.0);
    }
}

void accumulate_collector(collector *BASE, collector *C) {
    int	    i;
    assert(BASE->NTOTAL == C->NTOTAL);
    for (i=0; i<BASE->NTOTAL; i++) {
	BASE->G[i] += C->G[i];
    }
}

void collect(collector *C, particle *P, trace *T) {
    int		i, idx;
    ++C->NTRACE;
#if FULL_TRACE
    for (i=1; i<T->np; i++) {
	idx = T->data[i].where;
	if (idx == -1) continue;
	C->G[idx] += (T->data[i-1].e - T->data[i].e);
	++C->NSTEP;
    }
#else
    for (i=0; i<T->np; i++) {
	idx = T->data[i].where;
	if (idx == -1) continue;
	C->G[idx] += T->data[i].e;
	++C->NSTEP;
    }
#endif
    if (T->np > C->NMAX) C->NMAX = T->np;
}

static int get_index(int ix, int iy, int iz) {
    int	    idx;
    assert(ix >= 0);
    assert(iy >= 0);
    assert(iz >= 0);
    assert(ix < ARENA->N[CX]);
    assert(iy < ARENA->N[CY]);
    assert(iz < ARENA->N[CZ]);
    idx = ARENA->N[CX]*ARENA->N[CY]*iz + ARENA->N[CX]*iy + ix;
#if DEBUG
    fprintf(stderr, "%d %d %d -> %d\n", ix, iy, iz, idx);
#endif
    return(idx);
}

arena *initialize_collect(char *MRIFIL, char *LOCFIL) {
    char		*temp, *copy;
    int			size, nrow, ncol, npln;
    POSITION		xlo, xhi, ylo, yhi, zlo, zhi;
    arena		*A;
#if ! UNIFORM
    int			i, j, k, l, idx;
    enum material_type	*B;
    char		line[MAXLIN];
    FILE		*F;
#endif
    if (MRIFIL) {
	copy = strdup(MRIFIL);
/* if an MRI file was given, use it to define the collection region */
	temp = strtok(copy,"_.");
	temp = strtok(0,"_."); nrow = atoi(temp);   /* row = Y */
	temp = strtok(0,"_."); ncol = atoi(temp);   /* col = X */
	temp = strtok(0,"_."); npln = atoi(temp);   /* plane = Z */
	if (VERBOSE) {
	    printf("assuming file %s is X(%d) x Y(%d) x Z(%d)\n", MRIFIL, ncol, nrow, npln);
	}
	nx = ncol;
	xlo = -(ncol/2) * PIXEL;
	xhi = (ncol/2) * PIXEL;
	ny = nrow;
	ylo = -(nrow/2) * PIXEL;
	yhi = (nrow/2) * PIXEL;
	nz = npln;
	zlo = AS_REAL(0.0);
	zhi = npln * PIXEL;
    } else {
/* otherwise, make a default collector */
//	nx = ny = nz = 100;
	xlo = AS_REAL(0.0);	    /* mm */
	xhi = AS_REAL(200.0);	    /* mm */
	ylo = AS_REAL(-10.0);	    /* mm */
	yhi = AS_REAL(10.0);	    /* mm */
	zlo = AS_REAL(-10.0);	    /* mm */
	zhi = AS_REAL(10.0);	    /* mm */
    }
    size = sizeof(arena) + (nx*ny*nz - 1)*sizeof(short);
    A = (arena *) malloc(size);
    A->mysize = size;
/* store the box description */
    A->N[CX] = nx;
    A->MIN[CX] = xlo;
    A->MAX[CX] = xhi;
    A->WID[CX] = xhi - xlo;
    A->INVWID[CX] = A->N[CX]/A->WID[CX];
    A->CEN[CX] = (xlo + xhi)/AS_REAL(2.0);

    A->N[CY] = ny;
    A->MIN[CY] = ylo;
    A->MAX[CY] = yhi;
    A->WID[CY] = yhi - ylo;
    A->INVWID[CY] = A->N[CY]/A->WID[CY];
    A->CEN[CY] = (ylo + yhi)/AS_REAL(2.0);

    A->N[CZ] = nz;
    A->MIN[CZ] = zlo;
    A->MAX[CZ] = zhi;
    A->WID[CZ] = zhi - zlo;
    A->INVWID[CZ] = A->N[CZ]/A->WID[CZ];
    A->CEN[CZ] = (zlo + zhi)/AS_REAL(2.0);

    A->RADIUS = SQRT(A->WID[CX]*A->WID[CX] + A->WID[CY]*A->WID[CY] + A->WID[CZ]*A->WID[CZ]) * 0.5;

    A->NTOTAL = nx*ny*nz; 
    //總共切幾等份
    if (VERBOSE) {
	printf("bounding box: x(%.1f %.1f) y(%.1f %.1f) z(%.1f %.1f) mm\n", xlo, xhi, ylo, yhi, zlo, zhi);
	printf("center: x(%.1f) y(%.1f) z(%.1f) radius=%.2f\n", A->CEN[CX], A->CEN[CY], A->CEN[CZ], A->RADIUS);
	printf("widths: x(%.1f) y(%.1f) z(%.1f)\n", A->WID[CX], A->WID[CY], A->WID[CZ]);
	printf("created is_tumor vector to be %d long\n", A->NTOTAL);
    }
#if UNIFORM
    A->mat = WATER;
    if (VERBOSE && (MRIFIL || LOCFIL)) {
	printf("warning: ignoring MRI and Location files for uniform case!\n");
    }
#else
    if (MRIFIL) {
	if (VERBOSE) printf("attempting to open %s\n", MRIFIL);
	F = fopen(MRIFIL, "rb");
	B = (enum material_type *) malloc(A->NTOTAL * sizeof(enum material_type));
	i = (int) fread(B, sizeof(enum material_type), A->NTOTAL, F);
	if (i != A->NTOTAL) {
	    fprintf(stderr, "found only %d material in file %s\n", i, MRIFIL);
	    exit(-1);
	}
	for (i=0; i<A->NTOTAL; i++) {
	    A->matG[i] = (short) B[i];
	}
	free((char *) B);
	fclose(F);
    } else {
	for (i=0; i<A->NTOTAL; i++) {
	    A->matG[i] = (short) WATER;
	}
    }
    if (LOCFIL) {
	if (VERBOSE) printf("attempting to open %s\n", LOCFIL);
	F = fopen(LOCFIL, "r");
	while (fgets(line, MAXLIN, F)) {
	    sscanf(line, "%d %d %d %d", &i, &j, &k, &l);
	    idx = getidx(i, j, k, A->N[CX], A->N[CY], A->N[CZ]);
	    if (l > 0) A->matG[idx] |= TUMOR_DELTA;    /* is a tumor */
	    if (l < 0) A->matG[idx] |= EXCLU_DELTA;    /* is exclusion region */
	}
	fclose(F);
    }
#endif
/* and return the arena */
    return(A);
}

int WhereAmI(POSITION x, POSITION y, POSITION z) {
    int		idx, m[3];
    m[CX] = (int)( (x - ARENA->MIN[CX])*ARENA->INVWID[CX] );
    m[CY] = (int)( (y - ARENA->MIN[CY])*ARENA->INVWID[CY] );
    m[CZ] = (int)( (z - ARENA->MIN[CZ])*ARENA->INVWID[CZ] );
    if (m[CX] < 0 || m[CX] > ARENA->N[CX] ||
	m[CY] < 0 || m[CY] > ARENA->N[CY] ||
	m[CZ] < 0 || m[CZ] > ARENA->N[CZ]) {
	    idx = -1;
    } else {
	idx = getidx(m[CX],m[CY],m[CZ],ARENA->N[CX],ARENA->N[CY],ARENA->N[CZ]);
    }
#if DEBUG
    printf("WhereAmI %7.3f,%7.3f,%7.3f -> %d,%d,%d -> %d\n", x, y, z, m[CX], m[CY], m[CZ], idx);
#endif
    return(idx);
}
	
void summarize_collect(collector *C, int DOPLOT) {
    FILE	*F;
    REAL	max, *SX, *SY, *SZ;
    int		i, j, k, idx, peak[3];
    char	line[MAXLIN];
    time_t	now;
/* otherwise, we compute some simple statistics */
    if (VERBOSE) {
	max = 0;
	peak[CX] = peak[CY] = peak[CZ] = 0;
/* these arrays store the plane sums */
	SX = (REAL *) malloc(C->N[CX]*sizeof(REAL));
	for (i=0; i<C->N[CX]; i++) SX[i] = AS_REAL(0.0);
	SY = (REAL *) malloc(C->N[CY]*sizeof(REAL));
	for (i=0; i<C->N[CY]; i++) SY[i] = AS_REAL(0.0);
	SZ = (REAL *) malloc(C->N[CZ]*sizeof(REAL));
	for (i=0; i<C->N[CZ]; i++) SZ[i] = AS_REAL(0.0);
/* process the points */
  for (i=0; i<C->N[CX]; i++) {
      for (j=0; j<C->N[CY]; j++) {
          for (k=0; k<C->N[CZ]; k++) {
              idx = get_index(i, j, k);
              SX[i] += C->G[idx];
              SY[j] += C->G[idx];
              SZ[k] += C->G[idx];
              if (C->G[idx] > max) {
                  max = C->G[idx];
                  peak[CX] = i;
                  peak[CY] = j;
                  peak[CZ] = k;
              }
          }
	    }
	}
  
  printf("Total number of traces: %d\n", C->NTRACE);
	printf("Total number of steps: %d\n", C->NSTEP);
	printf("Maximum trace length: %d\n", C->NMAX);
	printf("Depth (mm) vs. Energy (MeV) curve\n");
	for (i=0; i<C->N[CX]; i++) {
	    printf("%6.1f %8.2e ", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]); ++i;
	    printf("%6.1f %8.2e ", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]); ++i;
	    printf("%6.1f %8.2e ", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]); ++i;
	    printf("%6.1f %8.2e ", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]); ++i;
	    printf("%6.1f %8.2e ", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]); ++i;
	    printf("%6.1f %8.2e ", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]); ++i;
	    printf("%6.1f %8.2e ", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]); ++i;
	    printf("%6.1f %8.2e\n", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]);
	}
	if (DOPLOT) {
	    now = time(NULL);
	    ctime_r(&now, &(line[0]));
	    F = fopen("jack.xgr", "w");
	    fprintf(F, "TitleText: Energy (MeV) vs. Distance (mm) curve from Jack run on %s", line);
	    fprintf(F, "\n\"Along_X\"\n\n");
	    for (i=0; i<ARENA->N[CX]; i++) {
		fprintf(F, " %7.4f %12.5e\n", (ARENA->MIN[CX]+i*ARENA->WID[CX]/ARENA->N[CX]), SX[i]);
	    }
	    fprintf(F, "\n\"Along_Y\"\n\n");
	    for (i=0; i<ARENA->N[CY]; i++) {
		fprintf(F, " %7.4f %12.5e\n", (ARENA->MIN[CY]+i*ARENA->WID[CY]/ARENA->N[CY]), SY[i]);
	    }
	    fprintf(F, "\n\"Along_Z\"\n\n");
	    for (i=0; i<ARENA->N[CZ]; i++) {
		fprintf(F, " %7.4f %12.5e\n", (ARENA->MIN[CZ]+i*ARENA->WID[CZ]/ARENA->N[CZ]), SZ[i]);
	    }
	    fclose(F);
	}
	printf("Peak energy = %g\n", max);
	printf("Peak location = %g,%g,%g\n",
	       (ARENA->MIN[CX] + peak[CX]*ARENA->WID[CX]/ARENA->N[CX]), 
	       (ARENA->MIN[CY] + peak[CY]*ARENA->WID[CY]/ARENA->N[CY]), 
	       (ARENA->MIN[CZ] + peak[CZ]*ARENA->WID[CZ]/ARENA->N[CZ]));
    }
}
