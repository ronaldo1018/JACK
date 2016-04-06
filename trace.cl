#include    "cl_set.h"

//#pragma OPENCL EXTENSION all : enable   
#pragma OPENCL EXTENSION cl_khr_fp64  : enable   
#pragma OPENCL EXTENSION cl_amd_printf  : enable   
//#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

/* turn on for printouts */
#define	DEBUG	    0
/* turn on for table printouts */
#define	TBL_DEBUG   0
/* turn on for queu printouts */
#define	QUE_DEBUG   0
/* turn on for parameter handling */
#define	PAR_DEBUG   0
/* turn on for angle printouts */
#define	DEL_DEBUG   0
/* turn on for elastic printouts */
#define	ELS_DEBUG   0
/* turn on for angle calculation printouts */
#define ANG_DEBUG   0
/* turn on for trace printouts */
#define TRC_DEBUG   0
/* keep off until we have elastic events working */
#define DO_ELASTIC  0
/* if true then real numbers are modelled as floats, else double */
#define USE_FLOAT   1//1
/* if true then use fast transcendentals */
#define USE_FAST    1//1
/* if true then store full traces */
#define FULL_TRACE  0//0
/* if true then only update trace when index changes */
#define MIN_TRACE   1
/* if true then use simple water phantom */
#define UNIFORM	    1
/* is we turn on fast functions, then turn on floats as well */
#if USE_FAST
#define	USE_FLOAT   1
#endif

#define CL_MATH_NATIVE	    1

#if USE_FLOAT
    /* what is a real number */
    #define	REAL	    float
    /* what is a real number for the tables */
    #define	TBL_REAL    float
    /* macro to format a constant */
    #define	AS_REAL(x)  (x##f)
    /* intrinsics */
//    #define	LOG(x)		    logf(x)
//    #define	EXP(x)		    expf(x)
//    #define	SQRT(x)		    sqrtf(x)
//    #define	SIN(x)		    sinf(x)
//    #define	COS(x)		    cosf(x)
//    #define	ABS(x)		    fabsf(x)
//    #define	ATAN2(x,y)	    atan2f(x,y)
//    #define	SINH(x)		    sinhf(x)
//    #define	COPYSIGN(x,y)	    copysignf(x,y)
#else
    /* what is a real number */
    #define	REAL	    double
    /* what is a real number for the tables */
    #define	TBL_REAL    double
    /* macro to format a constant */
    #define	AS_REAL(x)  (x)
    /* intrinsics */
//    #define	LOG(x)		    log(x)
//    #define	EXP(x)		    exp(x)
//    #define	SQRT(x)		    sqrt(x)
//    #define	SIN(x)		    sin(x)
//    #define	COS(x)		    cos(x)
//    #define	ABS(x)		    fabs(x)
//    #define	ATAN2(x,y)	    atan2(x,y)
//    #define	SINH(x)		    sinh(x)
//    #define	COPYSIGN(x,y)	    copysign(x,y)
#endif

    #define	LOG(x)		    log(x)
    #define	EXP(x)		    exp(x)
    #define	SQRT(x)		    sqrt(x)
    #define	SIN(x)		    sin(x)
    #define	COS(x)		    cos(x)
    #define	ABS(x)		    fabs(x)
    #define	ATAN2(x,y)	    atan2(x,y)
    #define	SINH(x)		    sinh(x)
    #define	COPYSIGN(x,y)	    copysign(x,y)

#define TPI_FLOAT    AS_REAL(6.28318531)
#define PI_FLOAT     AS_REAL(3.14159265)
#define PIBY2_FLOAT  AS_REAL(1.57079633)
#define	FOUR_PI	     AS_REAL(1.27323954)
#define	FOUR_PISQ    AS_REAL(0.40528473)

#define	    NITER   10

#if	USE_FAST
    #undef	LOG
    #undef	SIN
    #undef	COS
    #undef	ATAN2
    #undef	SINH
    #undef	EXP
    #define	ATAN2(x,y)  fast_atan2(x,y)
    #define	SINH(x)	    fast_sinh(x)
    
    #if CL_MATH_NATIVE
        #define	SIN(x)		  native_sin(x)
        #define	COS(x)		  native_cos(x)
        #define	EXP(x)		  native_exp(x)
        #define	LOG(x)		  native_log(x)
    #else
        #define	SIN(x)	    fast_sin(x)
        #define	COS(x)	    fast_cos(x)
        #define	EXP(x)	    fast_exp(x)
        #define	LOG(x)	    fast_log(x)
    #endif

#endif


REAL fast_sin(REAL theta) {
    REAL    y;
    while (theta < -PI_FLOAT) { theta += TPI_FLOAT; }
    while (theta > PI_FLOAT) { theta -= TPI_FLOAT; }
    y = FOUR_PI * theta - FOUR_PISQ * theta * ABS(theta);
    y = AS_REAL(0.225) * (y * ABS(y) - y) + y;
    return(y);
}

REAL fast_cos(REAL theta) {
    return(fast_sin(PIBY2_FLOAT - theta));
}

REAL fast_atan2(REAL y, REAL x) {
    REAL    atan, z;
    if (x == AS_REAL(0.0)) {
	if (y > AS_REAL(0.0)) return(PIBY2_FLOAT);
	if (y == AS_REAL(0.0)) return(AS_REAL(0.0));
	return(-PIBY2_FLOAT);
    }
    z = y/x;
    if (ABS(z) < AS_REAL(1.0)) {
	atan = z/(AS_REAL(1.0) + AS_REAL(0.28)*z*z);
	if (x < AS_REAL(0.0)) {
	    if (y < AS_REAL(0.0)) return(atan - PI_FLOAT);
	    return(atan + PI_FLOAT);
	}
    } else {
	atan = PIBY2_FLOAT - z/(z*z + AS_REAL(0.28));
	if (y < AS_REAL(0.0)) return(atan - PI_FLOAT);
    }
    return(atan);
}

REAL fast_log(REAL x) {
#if USE_FLOAT
/* see http://fastapprox.googlecode.com/svn/trunk/fastapprox/src/fastonebigheader.h */
    union { float f; unsigned int i; } vx = { x };
    union { unsigned int i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
    float y = vx.i;
    y *= 1.1920928955078125e-7f;
    return(0.69314718f * (y - 124.22551499f - 1.498030302f * mx.f - 1.72587999f / (0.3520887068f + mx.f)));
#else
    return(log(x));
#endif
}

REAL fastpow2 (REAL p) {
    REAL offset = (p < AS_REAL(0.0)) ? AS_REAL(1.0) : AS_REAL(0.0);
    REAL clipp = (p < AS_REAL(-126.0)) ? AS_REAL(-126.0) : p;
    int w = clipp;
    REAL z = clipp - w + offset;
    union { unsigned int i; REAL f; } v = { (unsigned int) ( (1 << 23) * (clipp + AS_REAL(121.2740575) + AS_REAL(27.7280233) / (AS_REAL(4.84252568) - z) - AS_REAL(1.49012907) * z) ) };
    return v.f;
}

REAL fast_exp(REAL p) {
    return(fastpow2(p * AS_REAL(1.442695040)));
}

REAL fast_sinh(REAL p) {
    REAL   z;
    z = fastpow2(p * AS_REAL(1.442695040));
    return( AS_REAL(0.5) * (z - AS_REAL(1.0)/z) );
}




/* what is a position, units: micrometers */
#define	POSITION    REAL
/* what is an energy, units: kilo-electron-volts */
#define	ENERGY	    REAL
/* what is a direction, units: none */
#define	DIRECTION   REAL
/* what is an angle, unit: radians */
#define	ANGLE	    REAL
/* maximum (constant) distance moved per step, unit is micrometers */
//#define	MAXSTEP	    0.01
/* maximum length of a trace */
#if FULL_TRACE
#define	MAXTRACE    8192
#else
#define	MAXTRACE    128
#endif	/* FULL_TRACE */
/* maximum length of an input line */
#define	MAXLIN	    1024
/* maximum length of particle queu */
#define	NQUEU	    64

enum particle_type { PROTON, ELECTRON, NEUTRON, ALPHA, UNKNOWN_P=-1 };
/* number of material types */
#define	NMATERIAL   7
/* material types and names */
enum material_type { WATER, AIR, ADIPOSETISSUEIRCP, A150TISSUE, MUSCLEWITHSUCROSE, B100BONE, OUTSIDE, UNKNOWN_M=-1 };
/* the three types of events */
enum event_type { IONIZATION, ELASTIC, INELASTIC };
/* types of parameters passed on command line */
enum ParamType	{ P_INTEGER, P_BOOLEAN, P_REAL, P_STRING, P_LIST, P_FLAG, P_PATH, P_READ, P_WRIT };

enum plane_type { NORTH, EAST, SOUTH, WEST, TOP, BOTTOM };


/* to help make the code clearer */
#define	CX   0	    /* columns */
#define	CY   1	    /* rows */
#define	CZ   2	    /* planes */
#define getidx(x,y,z,nx,ny,nz)    ((nx)*(ny)*(z) + (nx)*(y) + (x))
#define	INDEX(m,e,p)    ((m)*NE*NQ + (e)*NQ + (p))

#define	OLD_SCATTER 1

typedef struct trace_point {
#if FULL_TRACE
    POSITION		    x, y, z;	/* position (not really needed) */
    ENERGY		    e;		/* energy at this point */
    enum material_type	    m;		/* type of material (not really needed) */
    int			    where;	/* index in voxel space */
#else
    ENERGY		    e;		/* energy lost in this voxel */
    int			    where;	/* index in voxel space */
#endif
} trace_point;

typedef struct trace {
    int		    np;		    /* number of points in the trace */
#if ! FULL_TRACE
    ENERGY	    last;	    /* previous energy, electron-volt (=1.60217653E-19 Joules) */
#endif
    trace_point	    data[MAXTRACE]; /* trace data */
} trace;

typedef struct particle {
    POSITION	pos[3];	    /* starting position, mm */
    DIRECTION	dir[3];	    /* normalized direction vector */
    ENERGY	energy;	    /* starting energy, electron-volt (=1.60217653E-19 Joules) */
    POSITION	to_elastic; /* distance to next elastic event */
    POSITION	to_inelastic;/* distance to next inelastic event */
    enum particle_type type;/* type */
#if MIN_TRACE
    ENERGY	eloss;
    int		where;
#endif
} particle;

typedef struct arena {
    POSITION		MIN[3], WID[3], INVWID[3], MAX[3], CEN[3], RADIUS;
    int			N[3], NTOTAL;
    int			mysize;	    /* total size in bytes of this arena */
#if UNIFORM
    enum material_type	mat;
#else
    short		matG[1000000];
#endif
} arena;

typedef struct table {
    enum particle_type	P;	    /* type of particle for this table */
    enum material_type	M;	    /* type of particle for this table */
    int			np;	    /* number of points in the table */
    int			mysize;	    /* total size in bytes of this table */
    TBL_REAL		min, max;   /* limits of x variable */
    TBL_REAL		range;	    /* range of x variable */
    TBL_REAL		step;	    /* step size in the x variable */
    TBL_REAL		YV[201];	    /* array of positions, ith position is x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2]  */
} table;

typedef struct collector {
    int			N[3], NTOTAL, NTRACE, NSTEP, NMAX;
    int			mysize;	    /* total size in bytes of this collector */
    ENERGY		G[1000000];
} collector;

//__constant int VERBOSE = 0;
__constant REAL g_scale1 = AS_REAL(1.0)/0x7FFF;
__constant REAL g_scale2 = AS_REAL(2.0)/0x7FFF;

REAL fastrand01(unsigned int *g_seed) { 
    REAL a;
    *g_seed = (214013*(*g_seed)+2531011); 
    a = ((*g_seed>>16)&0x7FFF) * g_scale1;
    return(a);
} 

REAL fastrand(unsigned int *g_seed) { 
    REAL a;
    *g_seed = (214013*(*g_seed)+2531011); 
    a = ((*g_seed>>16)&0x7FFF) * g_scale2 - AS_REAL(1.0);
    return(a);    
} 

//float normal(float mean, float sigma, unsigned int *inseed) {
//    float x1, x2, w, y1;
//    float y2;
//    int use_last = 0;
//    if (use_last) {                
//        y1 = y2;
//        use_last = 0;
//    } else {
//        do {
//            x1 = 2.0 * fastrand(inseed) - 1.0;
//            x2 = 2.0 * fastrand(inseed) - 1.0;
//            w = x1 * x1 + x2 * x2;
//	} while ( w >= 1.0 );
//
//        w = sqrt( (-2.0 * log( w ) ) / w );
//        y1 = x1 * w;
//        y2 = x2 * w;
//        use_last = 1;
//    }
//    return(mean + y1*sigma);
//}

REAL normal(REAL mean, REAL sigma, unsigned int *inseed) {
    REAL x1, x2, w;
    do {
        x1 = fastrand(inseed);
        x2 = fastrand(inseed);
        w = x1 * x1 + x2 * x2;
    } while ( w >= AS_REAL(1.0) );
    w = SQRT( (-AS_REAL(2.0) * LOG( w ) ) / w );
    return(mean + x1*w*sigma);
}

//REAL normal(REAL mean, REAL sigma, unsigned int *inseed) {
//    REAL x1, x2, r, d, u1;
//        x1 = fastrand(inseed);
//        x2 = fastrand(inseed);
//        r = SQRT( AS_REAL(2.0) * LOG( x1 ));
//        d = TPI_FLOAT * x2;
//        u1 = r * COS(d);
////    u2 = r * SIN(d);
//    
//    return(mean + u1*sigma);
//}

//int WhereAmI(POSITION x, POSITION y, POSITION z, global arena *ARENA) {
//    int		idx, j, m[3];
//    m[CX] = (int)( (x - ARENA->MIN[CX])*ARENA->INVWID[CX] );
//    m[CY] = (int)( (y - ARENA->MIN[CY])*ARENA->INVWID[CY] );
//    m[CZ] = (int)( (z - ARENA->MIN[CZ])*ARENA->INVWID[CZ] );
//    for (j=0; j<3; j++) {
//        if (m[j] < 0 || m[j] >= ARENA->N[j]) return(-1);
//    }
//    idx = getidx(m[CX],m[CY],m[CZ],ARENA->N[CX],ARENA->N[CY],ARENA->N[CZ]);
//    return(idx);
//}

int WhereAmI(POSITION x, POSITION y, POSITION z, global arena *ARENA) {
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
	
int inside(POSITION *X, global arena *ARENA) {
    if (X[CX] < ARENA->MIN[CX] || X[CX] > ARENA->MAX[CX] ||
        X[CY] < ARENA->MIN[CY] || X[CY] > ARENA->MAX[CY] ||
        X[CZ] < ARENA->MIN[CZ] || X[CZ] > ARENA->MAX[CZ]) return(0);
    return(1);
}

void position_particle(ANGLE phi, ANGLE theta, 
                       particle *P,
                       global arena *ARENA) {
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
///* check */
//    assert(inside(A));
//    assert(! inside(B));
//#if GOM_DEBUG
//    printf("entering position particle: phi=%7.3f theta=%7.3f\n\tA=(%7.3f,%7.3f,%7.3f) B=(%7.3f,%7.3f,%7.3f)\n",
//	   phi, theta, A[CX], A[CY], A[CZ], B[CX], B[CY], B[CZ]);
//#endif
/* loop for NITER times */
    for (i=0; i<NITER; i++) {
        C[CX] = AS_REAL(0.5)*(A[CX] + B[CX]);
        C[CY] = AS_REAL(0.5)*(A[CY] + B[CY]);
        C[CZ] = AS_REAL(0.5)*(A[CZ] + B[CZ]);
        if (inside(C, ARENA)) {
            A[CX] = C[CX]; A[CY] = C[CY]; A[CZ] = C[CZ];
        } else {
            B[CX] = C[CX]; B[CY] = C[CY]; B[CZ] = C[CZ];
        }
    }
//#if GOM_DEBUG
//    printf("after binary search C: %7.3f %7.3f %7.3f\n", C[CX], C[CY], C[CZ]);
//#endif
/* calculate distances to nearest plane in X/Y/Z */
    d1 = ABS(A[CX]-ARENA->MIN[CX]); d2 = ABS(A[CX]-ARENA->MAX[CX]);
    if (d1 < d2) { dx = d1; S[CX] = WEST; } else { dx = d2; S[CX] = EAST; }
    d1 = ABS(A[CY]-ARENA->MIN[CY]); d2 = ABS(A[CY]-ARENA->MAX[CY]);
    if (d1 < d2) { dy = d1; S[CY] = NORTH; } else { dy = d2; S[CY] = SOUTH; }
    d1 = ABS(A[CZ]-ARENA->MIN[CZ]); d2 = ABS(A[CZ]-ARENA->MAX[CZ]);
    if (d1 < d2) { dz = d1; S[CZ] = BOTTOM; } else { dz = d2; S[CZ] = TOP; }
//#if GOM_DEBUG
//    printf("distances: %7.3f %s %7.3f %s %7.3f %s\n",
//	    dx, plane_name[(int) S[CX]], dy, plane_name[(int) S[CY]], dz, plane_name[(int) S[CZ]]);
//#endif
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
//#if GOM_DEBUG
//    printf("final C: %7.3f %7.3f %7.3f -> %7.3f %7.3f %7.3f\n", P->pos[CX], P->pos[CY], P->pos[CZ], P->dir[CX], P->dir[CY], P->dir[CZ]);
//#endif
}

/* this function needs a careful look for the case when x is not real */
#if L_TBL
TBL_REAL interpolate_table(local table *T, TBL_REAL x) {
#else	
TBL_REAL interpolate_table(constant table *T, TBL_REAL x) {
#endif	
    int		i;
    TBL_REAL	y, z;
    if (x <= T->min) {
        return(T->YV[0]);
    }
    if (x >= T->max) {
        return(T->YV[T->np-1]);
    }
    i = (int)((x - T->min)/T->step);
//    assert(i >= 0 && i <= T->np-1);
    z = (x - i*T->step)/T->step;
    y = T->YV[i] + z * (T->YV[i+1] - T->YV[i]);
#if TBL_DEBUG
    printf("interpolate %g -> %d + %g -> %g %g -> %g\n", x, i, z, T->YV[i], T->YV[i+1], y);
#endif
    return(y);
}

void rotateUz(REAL *oldDir, REAL *newDir) {
    REAL    u1, u2, u3, up, dx, dy, dz, px, py, pz;
    u1 = oldDir[0];
    u2 = oldDir[1];
    u3 = oldDir[2];
    up = u1*u1 + u2*u2;
    dx = newDir[0];
    dy = newDir[1];
    dz = newDir[2];

    if (up>0) {
        up = sqrt(up);
        px = dx;
        py = dy;
        pz = dz;
        dx = (u1*u3*px - u2*py)/up + u1*pz;
        dy = (u2*u3*px + u1*py)/up + u2*pz;
        dz =    -up*px +             u3*pz;
    } else if (u3 < 0.) {
        dx = -dx; dz = -dz; 
    }
    
    newDir[0] = dx;
    newDir[1] = dy;
    newDir[2] = dz;
}

void ProcessDeltaDirection(
                           particle *P,
                           ANGLE Phi, ANGLE Theta) {
    REAL	cth, sth, oldDir[3], newDir[3];
#if DEL_DEBUG
    printf("BEFORE: %7.4f %7.4f %7.4f (theta=%7.4f phi=%7.4f\n",
	  P->dir[CX], P->dir[CY], P->dir[CZ], Theta, Phi);
#endif
/* calculate sin/cos of theta */
    cth = COS(Theta);
    sth = SQRT(1.0 - cth*cth);
/* old direction */
    oldDir[0] = P->dir[0];	oldDir[1] = P->dir[1];	    oldDir[2] = P->dir[2];
/* new direction */
    newDir[0] = sth*COS(Phi);	newDir[1] = sth*SIN(Phi);   newDir[2] = cth;
/* do the rotation */
    rotateUz(oldDir, newDir);
/* place into particle */
    P->dir[0] = newDir[0];  P->dir[1] = newDir[1];      P->dir[2] = newDir[2];
#if DEL_DEBUG
    printf("AFTER:  %7.4f %7.4f %7.4f\n",
	  P->dir[CX], P->dir[CY], P->dir[CZ]);
#endif

#ifdef NOWAY
/* calculate Theta (in X/Y plane) */
    Theta = ATAN2(P->dir[CY], P->dir[CX]);
/* calculate the magnitude of the vector in the X/Y plane */
    XY = SQRT(P->dir[CX]*P->dir[CX] + P->dir[CY]*P->dir[CY]);
/* calculate Phi */
    Phi = ATAN2(XY, P->dir[CZ]); 
#if DEL_DEBUG
    printf("BEFORE: %7.4f %7.4f %7.4f (theta=%7.4f (%7.4f) phi=%7.4f (%7.4f)\n",
	  P->dir[CX], P->dir[CY], P->dir[CZ], Theta, dTheta, Phi, dPhi);
#endif
/* add perturbation */
    Theta += dTheta;
/* add perturbation */
    Phi += dPhi;
/* now determine new direction */
    P->dir[CZ] = COS(Phi);
    tmp = SIN(Phi);
    P->dir[CX] = tmp * COS(Theta);
    P->dir[CY] = tmp * SIN(Theta);
#if DEL_DEBUG
    printf("AFTER:  %7.4f %7.4f %7.4f (theta=%7.4f phi=%7.4f\n",
	   P->dir[CX], P->dir[CY], P->dir[CZ], Theta, Phi);
#endif
#if DEL_DEBUG
    if (isnan(P->dir[CX]) || isnan(P->dir[CY]) || isnan(P->dir[CZ])) {
        fprintf(stderr, "error: generated a NaN with dTheta = %g dPhi = %g \n", dTheta, dPhi);
        fprintf(stderr, "old dir: %7.4f,%7.4f,%7.4f new dir: %7.4f,%7.4f,%7.4f\n",
	      P->dir[CX], P->dir[CY], P->dir[CZ], P->dir[CX], P->dir[CY], P->dir[CZ]);
        exit(-1);
    }
#endif
#endif	/* NOWAY */
}

TBL_REAL interpol8(TBL_REAL xlo, TBL_REAL xhi, TBL_REAL x,
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
    return((y+z) - COPYSIGN((REAL)SQRT(0.5*(y*y+z*z)), y));
}

TBL_REAL get_fluct(TBL_REAL E, enum material_type M, unsigned int *seed,
                   int NE, int NQ, int NM,
                   constant TBL_REAL *ENRG,
                   constant TBL_REAL *PERC,
                   constant TBL_REAL *FLUC) 
{
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
//    assert(k > 0);
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

//REAL	PHI_SCALE[NMATERIAL] = {
__constant REAL	PHI_SCALE[7] = {
/* WATER */		22.634e-3,
/* AIR */		0.7143e-3,
/* ADIPOSETISSUEIRCP */	23.134e-3,
/* A150TISSUE */	19.877e-3,
/* MUSCLEWITHSUCROSE */	22.237e-3,
/* B100BONE */		23.611e-3,
/* OUTSIDE */		0.0		};

ANGLE	SamplePhi(
                  particle *P,
                  enum material_type M, unsigned int *seed) {
    ENERGY  m = PHI_SCALE[(int) M];
    if (P->energy < 1.0) {
	return(0.0);
    } else {
	return( EXP(normal(AS_REAL(0.0), AS_REAL(0.693), seed))*m/P->energy );
    }
}

ANGLE	SampleTheta(
                  particle *P,
                  enum material_type M, unsigned int *seed) {
    return( AS_REAL(2.0) * M_PI * fastrand01(seed) );
}

void add_trace_point(global trace *T, 
                     particle *P,
                     int V, enum material_type M) {
    int		m;
/* if a trace grows too large, kill the particle */
    if (T->np >= MAXTRACE) {
	P->energy = AS_REAL(0.0);
	return;
    }
/* if we are currently outside then... */
    if (M == OUTSIDE) {
        if (T->np) {	/* if particle has been in, kill it */
	          P->energy = AS_REAL(0.0);
        }
        /* in either case, ignore it from the trace point of view */
        return;
    }

    #if FULL_TRACE
        T->data[T->np].x = P->pos[0]; T->data[T->np].y = P->pos[1]; T->data[T->np].z = P->pos[2];
        T->data[T->np].e = P->energy; T->data[T->np].m = M;         T->data[T->np].where = V;
        ++T->np;
    #elif MIN_TRACE
        T->data[T->np].e = P->eloss;  T->data[T->np].where = V;
        ++T->np;
    #else
/*    printf("T->np = %d e = %f last = %f where = %d\n", T->np, P->energy, T->last, V); */
    if (T->np == 0) {
        T->data[T->np].e = 0.0;
        T->data[T->np].where = V;
        ++T->np;
    } else {
        if (V == T->data[T->np-1].where) {
            T->data[T->np-1].e += T->last - P->energy;
        } else {
            T->data[T->np].e = T->last - P->energy;
            T->data[T->np].where = V;
            ++T->np;
        }
    }
    T->last = P->energy;
    #endif	/* FULL_TRACE */
}

void proton_ionization_event(
                             REAL MAXSTEP,                             
                             global arena *ARENA,
                             particle *P,
                             global trace *TRACE,
                             enum material_type M, unsigned int *seed, 
                             int NE, int NQ, int NM,
                             constant TBL_REAL *ENRG,
                             constant TBL_REAL *PERC,
                             constant TBL_REAL *FLUC,
                             #if L_TBL
                             local table *ProtonWater_Energy,
                             local table *ProtonAir_Energy,
                             local table *ProtonAdiposeTissue_Energy,
                             local table *ProtonATissue_Energy,
                             local table *ProtonMuscleWithSucrose_Energy,
                             local table *ProtonBBone_Energy
                             #else                     
                             constant table *ProtonWater_Energy,
                             constant table *ProtonAir_Energy,
                             constant table *ProtonAdiposeTissue_Energy,
                             constant table *ProtonATissue_Energy,
                             constant table *ProtonMuscleWithSucrose_Energy,
                             constant table *ProtonBBone_Energy
                             #endif
                             ) {
/* which table do I use for the delta-energy */
    ENERGY	dE;
    ANGLE	dTheta, dPhi;
    int		new = 0;
/* are we outside? */
    if (M == OUTSIDE) {
/* just update the position for the next event */
        P->pos[0] += P->dir[0] * MAXSTEP;
        P->pos[1] += P->dir[1] * MAXSTEP;
        P->pos[2] += P->dir[2] * MAXSTEP;
	      add_trace_point(TRACE, P, P->where, M);
        return;
    }
 /* look up dE/dx in table */

  switch(M) {
	    case WATER:
          dE = (ENERGY) interpolate_table(ProtonWater_Energy, (TBL_REAL) P->energy);
          break;
	    case AIR:
          dE = (ENERGY) interpolate_table(ProtonAir_Energy, (TBL_REAL) P->energy);
          break;
	    case ADIPOSETISSUEIRCP:
          dE = (ENERGY) interpolate_table(ProtonAdiposeTissue_Energy, (TBL_REAL) P->energy);
          break;
	    case A150TISSUE:
          dE = (ENERGY) interpolate_table(ProtonATissue_Energy, (TBL_REAL) P->energy);
          break;
	    case MUSCLEWITHSUCROSE:
          dE = (ENERGY) interpolate_table(ProtonMuscleWithSucrose_Energy, (TBL_REAL) P->energy);
          break;
	    case B100BONE:
          dE = (ENERGY) interpolate_table(ProtonBBone_Energy, (TBL_REAL) P->energy);
          break;
	    case OUTSIDE:
	    case UNKNOWN_M:
          break;
  }        
/*
   #if L_TBL
   local table* TBL_ptr[6];
   #else
   constant table* TBL_ptr[6];
   #endif
   TBL_ptr[0] = ProtonWater_Energy;
   TBL_ptr[1] = ProtonAir_Energy;
   TBL_ptr[2] = ProtonAdiposeTissue_Energy;
   TBL_ptr[3] = ProtonATissue_Energy;
   TBL_ptr[4] = ProtonMuscleWithSucrose_Energy;
   TBL_ptr[5] = ProtonBBone_Energy;

   dE = (ENERGY) interpolate_table(TBL_ptr[M], (TBL_REAL) P->energy);
*/

//  dE = (ENERGY) interpolate_table(ProtonWater_Energy, (TBL_REAL) P->energy);

//    printf(" --- ProtonWater_Energy->np = %d\n", ProtonWater_Energy->np);
     
//     dEdx = (ENERGY) interpolate_table(&T, (TBL_REAL) P->energy);
//    printf(" ---------- P->energy = %f, dEdx = %f\n", P->energy, dEdx);
 /* look up the statistics of the dE/dx fluctuations in table */
  dE += get_fluct(P->energy, M, seed, NE, NQ, NM, ENRG, PERC, FLUC);
/* and multiply by the step in mm to get the actual dE in MeV */
  dE *= MAXSTEP;

//     assert(dE >= AS_REAL(0.0));
  if (dE > P->energy) {
	    dE = P->energy;
      P->energy = AS_REAL(0.0);
  } else {
      P->energy -= dE;
  }
/* generate the angles */
    dPhi = SamplePhi(P, M, seed);
//    dPhi = 0.0;
    dTheta = SampleTheta(P, M, seed);
 
   /* update particle direction based on the angles */
     ProcessDeltaDirection(P, dTheta, dPhi);
 
 /* update the position for the next event */
     P->pos[0] += P->dir[0] * MAXSTEP;
     P->pos[1] += P->dir[1] * MAXSTEP;
     P->pos[2] += P->dir[2] * MAXSTEP;
#if MIN_TRACE
/* update the loss */
    P->eloss += dE;
/* get the current position */
    new = WhereAmI(P->pos[0], P->pos[1], P->pos[2], ARENA);
/* if not the same position as before, or the particle has died, update the trace */
    if (new != P->where || P->energy <= AS_REAL(0.0)) {
        add_trace_point(TRACE, P, P->where, M);
        P->eloss = AS_REAL(0.0);
        P->where = new;
    }
#endif

}


void proton_event(
                  REAL MAXSTEP,
                  particle *P,
                  global trace *TRACE, 
                  global arena *ARENA,
                  unsigned int *seed, 
                  int NE, int NQ, int NM,
                  constant TBL_REAL *ENRG,
                  constant TBL_REAL *PERC,
                  constant TBL_REAL *FLUC,
                  #if L_TBL
                  local table *ProtonWater_Energy,
                  local table *ProtonAir_Energy,
                  local table *ProtonAdiposeTissue_Energy,
                  local table *ProtonATissue_Energy,
                  local table *ProtonMuscleWithSucrose_Energy,
                  local table *ProtonBBone_Energy
                  #else                     
                  constant table *ProtonWater_Energy,
                  constant table *ProtonAir_Energy,
                  constant table *ProtonAdiposeTissue_Energy,
                  constant table *ProtonATissue_Energy,
                  constant table *ProtonMuscleWithSucrose_Energy,
                  constant table *ProtonBBone_Energy
                  #endif
                 ) 
{
/* the voxel we are in */
    int			    V;
/* what material am I in */
    enum material_type	    M;
#if MIN_TRACE
    V = P->where;
    #if UNIFORM
    M = (V == -1 ? OUTSIDE : ARENA->mat);
    #else
    M = (V == -1 ? OUTSIDE : ARENA->matG[V]);
    #endif
    proton_ionization_event(
                            MAXSTEP, 
                            ARENA, P, TRACE, M, seed, 
                            NE, NQ, NM,
                            ENRG, PERC, FLUC,
                            ProtonWater_Energy,
                            ProtonAir_Energy,
                            ProtonAdiposeTissue_Energy,
                            ProtonATissue_Energy,
                            ProtonMuscleWithSucrose_Energy,
                            ProtonBBone_Energy
                           );
#else
    V = WhereAmI(P->pos[0], P->pos[1], P->pos[2], ARENA);
    #if UNIFORM
    M = (V == -1 ? OUTSIDE : ARENA->mat);
    #else
    M = (V == -1 ? OUTSIDE : ARENA->matG[V]);
    #endif
    proton_ionization_event(
                            MAXSTEP, 
                            ARENA, P, TRACE, M, seed, 
                            NE, NQ, NM,
                            ENRG, PERC, FLUC,
                            ProtonWater_Energy,
                            ProtonAir_Energy,
                            ProtonAdiposeTissue_Energy,
                            ProtonATissue_Energy,
                            ProtonMuscleWithSucrose_Energy,
                            ProtonBBone_Energy
                           );
    add_trace_point(TRACE, P, V, M);
#endif
//    printf("gid = %d V = %d; M = %d == living....\n", gid, V, M);
                    	
}

void  tableG2L(local table *local_table, constant table *global_table, const int EntryNum)
{
    int i;
    local_table->P = global_table->P;
    local_table->M = global_table->M;
    local_table->np = global_table->np;
    local_table->mysize = global_table->mysize;
    local_table->min = global_table->min;
    local_table->max = global_table->max;
    local_table->range = global_table->range;
    local_table->step = global_table->step;

    for(i = 0; i < EntryNum; i++)
        local_table->YV[i] = global_table->YV[i];
        
    barrier(CLK_LOCAL_MEM_FENCE);
}

void collect(global trace *T, global ENERGY* collectG, 
             global int* NTRACE, global int* NSTEP, global int* NMAX) {
    int		i, idx;
    ++*NTRACE;
#if FULL_TRACE
    for (i=1; i<T->np; i++) {
        idx = T->data[i].where;
        if (idx == -1) continue;
        collectG[idx] += (T->data[i-1].e - T->data[i].e);
        ++*NSTEP;
    }
#else
    for (i=0; i<T->np; i++) {
        idx = T->data[i].where;
        if (idx == -1) continue;
        collectG[idx] += T->data[i].e;
        ++*NSTEP;
    }
#endif
    if (T->np > *NMAX) *NMAX = T->np;
}

void reset_trace(global trace *TRACE) {
    TRACE->np = 0;
#if ! (FULL_TRACE || MIN_TRACE)
    TRACE->last = AS_REAL(0.0);
#endif
}

__kernel void trace_kernel(
//                           __local int* weird_factor0,
                           __global int *ctrlflag,
                           const REAL MAXSTEP,
                           const REAL mean,
                           const REAL stdv,
                           const REAL sigma,
                           const REAL angle,
                           const REAL azimuth,
                           const int gridsize,
                           const int BatchItems,
                           __global POSITION *ppos,
                           __global DIRECTION *pdir,
                           __global unsigned int *seed,
                           __global trace* T_A,
                           __global trace* T_B,
                           __global arena* ARENA,
                           const int NE,
                           const int NQ,
                           const int NM,
                           __constant TBL_REAL *ENRG,
                           __constant TBL_REAL *PERC,
                           __constant TBL_REAL *FLUC,
                           const int Num_ProtonWater_Energy,
                           __local table *Loc_ProtonWater_Energy,
                           const int Num_ProtonAir_Energy,
                           __local table *Loc_ProtonAir_Energy,                           
                           const int Num_ProtonAdiposeTissue_Energy,
                           __local table *Loc_ProtonAdiposeTissue_Energy,
                           const int Num_ProtonATissue_Energy,
                           __local table *Loc_ProtonATissue_Energy,
                           const int Num_ProtonMuscleWithSucrose_Energy,
                           __local table *Loc_ProtonMuscleWithSucrose_Energy,
                           const int Num_ProtonBBone_Energy,
                           __local table *Loc_ProtonBBone_Energy,                          
                           __constant table *ProtonWater_Energy,
                           __constant table *ProtonAir_Energy,
                           __constant table *ProtonAdiposeTissue_Energy,
                           __constant table *ProtonATissue_Energy,
                           __constant table *ProtonMuscleWithSucrose_Energy,
                           __constant table *ProtonBBone_Energy
                           
                           )
{         
		                                
//#if L_TBL
  if (*ctrlflag == 11)
  {
    tableG2L(Loc_ProtonWater_Energy, ProtonWater_Energy, Num_ProtonWater_Energy);
    tableG2L(Loc_ProtonAir_Energy, ProtonAir_Energy, Num_ProtonAir_Energy);          
    tableG2L(Loc_ProtonAdiposeTissue_Energy, ProtonAdiposeTissue_Energy, Num_ProtonAdiposeTissue_Energy);
    tableG2L(Loc_ProtonATissue_Energy, ProtonATissue_Energy, Num_ProtonATissue_Energy);
    tableG2L(Loc_ProtonMuscleWithSucrose_Energy, ProtonMuscleWithSucrose_Energy, Num_ProtonMuscleWithSucrose_Energy);
    tableG2L(Loc_ProtonBBone_Energy, ProtonBBone_Energy, Num_ProtonBBone_Energy);
  }   
//#endif
        
  if (*ctrlflag == 67 || *ctrlflag == 66)
  {
    int gid = get_global_id(0);
    int lid = get_local_id(0);
    int gpid = get_group_id(0);
    int gsz = get_global_size(0);
    int lsz = get_local_size(0);
    
    POSITION		A[3], B[3], C[3], xy, d1, d2, dx, dy, dz;
    enum plane_type	S[3];
    REAL ang_nom, azi_nom;
    ANGLE ang, azi;
    unsigned int priv_seed;
    particle priv_p;
    
    priv_seed = seed[gid];    
//    printf(" gid = %d, seed[%d] = %u, prive_seed = %u\n", gid, gid, seed[gid], priv_seed);
    if (sigma > 0.0)
    {	
        ang_nom = angle + normal(0.0, sigma, &priv_seed);
        azi_nom = azimuth + normal(0.0, sigma, &priv_seed);
        ang = ang_nom * M_PI / AS_REAL(180.0);
        azi = azi_nom * M_PI / AS_REAL(180.0);
        position_particle(azi, ang, &priv_p, ARENA);
    }
    else
    {
        priv_p.pos[0] = ppos[0];
        priv_p.pos[1] = ppos[1];
        priv_p.pos[2] = ppos[2];
        priv_p.dir[0] = pdir[0];
        priv_p.dir[1] = pdir[1];
        priv_p.dir[2] = pdir[2];
    }

//    printf("gid = %d ==== living....\n", gid);
//    printf("lid = %d, gid = %d, addr p = %lu\n", lid, gid, &p[lid]);
 
    
#if MIN_TRACE
    priv_p.eloss = AS_REAL(0.0);
    priv_p.where = WhereAmI(priv_p.pos[CX], priv_p.pos[CY], priv_p.pos[CZ], ARENA);
#endif

    priv_p.type = PROTON;
    priv_p.energy = normal(mean, stdv, &priv_seed);    
    //Gen elastic distances
    #if DO_ELASTIC
    priv_p.to_elastic = 1e5*normal(4.0, 1.0, &priv_seed);
    priv_p.to_inelastic = 2e5*normal(4.0, 1.0, &priv_seed);
    #endif
   
    
    int flagidx = *ctrlflag & 1;
  	global trace* T_ptr[2];
    T_ptr[0] = T_A;
    T_ptr[1] = T_B;
    reset_trace(&T_ptr[flagidx][gid]);
//    REAL MAXSTEP = C_MAXSTEP;
//    printf("flag index = %d\n", flagidx);
    
    while (priv_p.energy > AS_REAL(0.0))
    {    
        #if L_TBL
                 proton_event(
                     MAXSTEP, 
                     &priv_p, 
                     &T_ptr[flagidx][gid], 
                     ARENA, &priv_seed, 
                     NE, NQ, NM,
                     ENRG, PERC, FLUC,
                     Loc_ProtonWater_Energy,
                     Loc_ProtonAir_Energy,
                     Loc_ProtonAdiposeTissue_Energy,
                     Loc_ProtonATissue_Energy,
                     Loc_ProtonMuscleWithSucrose_Energy,
                     Loc_ProtonBBone_Energy
                    );
        #else                     
                 proton_event(
                     MAXSTEP, 
                     &priv_p, 
                     &T_ptr[flagidx][gid], 
                     ARENA, &priv_seed, 
                     NE, NQ, NM,
                     ENRG, PERC, FLUC,
                     ProtonWater_Energy,
                     ProtonAir_Energy,
                     ProtonAdiposeTissue_Energy,
                     ProtonATissue_Energy,
                     ProtonMuscleWithSucrose_Energy,
                     ProtonBBone_Energy
                    );
        #endif
    }    
    seed[gid] = priv_seed;
    
  }    
  barrier(CLK_LOCAL_MEM_FENCE);
 
} 

