#include    "jack.h"

#define	    NITER   10

enum plane_type { NORTH, EAST, SOUTH, WEST, TOP, BOTTOM };
#if GOM_DEBUG
static char *plane_name[6] = { "NORTH", "EAST", "SOUTH", "WEST", "TOP", "BOTTOM" };
#endif

void print_particle(char *str, particle *P) {
    printf("%12s %8s(%7.4f,%7.4f,%7.4f) -> (%6.3f,%6.3f,%6.3f) E=%g (%g,%g)\n", str,
	   particle_name[P->type],
	   P->pos[CX],P->pos[CY],P->pos[CZ],
	   P->dir[CX],P->dir[CY],P->dir[CZ],
	   P->energy, P->to_elastic, P->to_inelastic);
}

static int inside(POSITION *X) {
    if (X[CX] < ARENA->MIN[CX] || X[CX] > ARENA->MAX[CX] ||
	X[CY] < ARENA->MIN[CY] || X[CY] > ARENA->MAX[CY] ||
	X[CZ] < ARENA->MIN[CZ] || X[CZ] > ARENA->MAX[CZ]) return(0);
    return(1);
}

void position_particle(ANGLE phi, ANGLE theta, particle *P) {
    int			i;
    POSITION		A[3], B[3], C[3], xy, d1, d2, dx, dy, dz;
    enum plane_type	S[3];
/* A is the center */
    A[CX] = ARENA->CEN[CX];
    A[CY] = ARENA->CEN[CY];
    A[CZ] = ARENA->CEN[CZ];
/* B is a point with phi/theta at Radius */
    xy = ARENA->RADIUS * SIN(phi);
    B[CX] = A[CX] + xy * COS(theta);
    B[CY] = A[CY] + xy * SIN(theta);
    B[CZ] = A[CZ] + ARENA->RADIUS * COS(phi);
/* check */
    assert(inside(A));
    assert(! inside(B));
#if GOM_DEBUG
    printf("entering position particle: phi=%7.3f theta=%7.3f\n\tA=(%7.3f,%7.3f,%7.3f) B=(%7.3f,%7.3f,%7.3f)\n",
	   phi, theta, A[CX], A[CY], A[CZ], B[CX], B[CY], B[CZ]);
#endif
/* loop for NITER times */
    for (i=0; i<NITER; i++) {
	C[CX] = AS_REAL(0.5)*(A[CX] + B[CX]);
	C[CY] = AS_REAL(0.5)*(A[CY] + B[CY]);
	C[CZ] = AS_REAL(0.5)*(A[CZ] + B[CZ]);
	if (inside(C)) {
	    A[CX] = C[CX]; A[CY] = C[CY]; A[CZ] = C[CZ];
	} else {
	    B[CX] = C[CX]; B[CY] = C[CY]; B[CZ] = C[CZ];
	}
    }
#if GOM_DEBUG
    printf("after binary search C: %7.3f %7.3f %7.3f\n", C[CX], C[CY], C[CZ]);
#endif
/* calculate distances to nearest plane in X/Y/Z */
    d1 = ABS(A[CX]-ARENA->MIN[CX]); d2 = ABS(A[CX]-ARENA->MAX[CX]);
    if (d1 < d2) { dx = d1; S[CX] = WEST; } else { dx = d2; S[CX] = EAST; }
    d1 = ABS(A[CY]-ARENA->MIN[CY]); d2 = ABS(A[CY]-ARENA->MAX[CY]);
    if (d1 < d2) { dy = d1; S[CY] = NORTH; } else { dy = d2; S[CY] = SOUTH; }
    d1 = ABS(A[CZ]-ARENA->MIN[CZ]); d2 = ABS(A[CZ]-ARENA->MAX[CZ]);
    if (d1 < d2) { dz = d1; S[CZ] = BOTTOM; } else { dz = d2; S[CZ] = TOP; }
#if GOM_DEBUG
    printf("distances: %7.3f %s %7.3f %s %7.3f %s\n",
	    dx, plane_name[(int) S[CX]], dy, plane_name[(int) S[CY]], dz, plane_name[(int) S[CZ]]);
#endif
/* find the closest one and snap the appropriate coordinate */
    if (dx <= dy) {
	if (dx <= dz) {
	    C[CX] = (S[CX] == WEST ? ARENA->MIN[CX] : ARENA->MAX[CX]);
	} else {
	    C[CZ] = (S[CZ] == BOTTOM ? ARENA->MIN[CZ] : ARENA->MAX[CZ]);
	}
    } else {
	if (dy <= dz) {
	    C[CY] = (S[CY] == NORTH ? ARENA->MIN[CY] : ARENA->MAX[CY]);
	} else {
	    C[CZ] = (S[CZ] == BOTTOM ? ARENA->MIN[CZ] : ARENA->MAX[CZ]);
	}
    }
    P->pos[CX] = C[CX]; P->pos[CY] = C[CY]; P->pos[CZ] = C[CZ];
    C[CX] -= ARENA->CEN[CX];
    C[CY] -= ARENA->CEN[CY];
    C[CZ] -= ARENA->CEN[CZ];
    d1 = SQRT(C[CX]*C[CX] + C[CY]*C[CY] + C[CZ]*C[CZ]);
    P->dir[CX] = -C[CX]/d1; P->dir[CY] = -C[CY]/d1; P->dir[CZ] = -C[CZ]/d1;
#if GOM_DEBUG
    printf("final C: %7.3f %7.3f %7.3f -> %7.3f %7.3f %7.3f\n", P->pos[CX], P->pos[CY], P->pos[CZ], P->dir[CX], P->dir[CY], P->dir[CZ]);
#endif
}

particle *MakeParticle(particle *P, REAL ENE, REAL ANG, REAL AZI, unsigned int *seed) {
    REAL	angle, azimuth;
    P->type = PROTON;
    P->energy = ENE;
    angle = ANG * M_PI / AS_REAL(180.0);
    azimuth = AZI * M_PI / AS_REAL(180.0);
    position_particle(azimuth, angle, P);
#if DO_ELASTIC
    GenerateElasticDistances(P, seed);
#endif
#if MIN_TRACE
    P->eloss = AS_REAL(0.0);
    P->where = WhereAmI(P->pos[CX], P->pos[CY], P->pos[CZ]);
#endif
#if GOM_DEBUG
    if (VERBOSE) {
	print_particle("Made particle:", P);
    }
#endif
    return(P);
}
