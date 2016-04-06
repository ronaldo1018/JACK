#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <assert.h>
#include    <math.h>
#include    "coord.h"

/* turn on for printouts */
#define	DEBUG	    0
/* turn on for table printouts */
#define	TBL_DEBUG   0
/* turn on for queu printouts */
#define	QUE_DEBUG   0
/* turn on for parameter handling */
#define	PAR_DEBUG   1
/* turn on for angle printouts */
#define	DEL_DEBUG   0
/* turn on for elastic printouts */
#define	ELS_DEBUG   0
/* turn on for angle calculation printouts */
#define ANG_DEBUG   0
/* turn on for trace printouts */
#define TRC_DEBUG   0
/* turn on for geometry printouts */
#define GOM_DEBUG   0

/* keep off until we have elastic events working */
#define DO_ELASTIC  0
/* if true then real numbers are modelled as floats, else double */
#define USE_FLOAT   1
/* if true then use fast transcendentals */
#define USE_FAST    1
/* if true then store full traces */
#define FULL_TRACE  0
/* if true then only update trace when index changes */
#define MIN_TRACE   1
/* if true then use simple water phantom */
#define UNIFORM	    1
/* is we turn on fast functions, then turn on floats as well */
#if USE_FAST
#define	USE_FLOAT   1
#endif

#if USE_FLOAT
/* what is a real number */
#define	REAL	    float
/* what is a real number for the tables */
#define	TBL_REAL    float
/* macro to format a constant */
#define	AS_REAL(x)  (x##f)
/* intrinsics */
#define	LOG(x)		    logf(x)
#define	EXP(x)		    expf(x)
#define	SQRT(x)		    sqrtf(x)
#define	SIN(x)		    sinf(x)
#define	COS(x)		    cosf(x)
#define	ABS(x)		    fabsf(x)
#define	ATAN2(x,y)	    atan2f(x,y)
#define	SINH(x)		    sinhf(x)
#define	COPYSIGN(x,y)	    copysignf(x,y)
#else
/* what is a real number */
#define	REAL	    double
/* what is a real number for the tables */
#define	TBL_REAL    double
/* macro to format a constant */
#define	AS_REAL(x)  (x)
/* intrinsics */
#define	LOG(x)		    log(x)
#define	EXP(x)		    exp(x)
#define	SQRT(x)		    sqrt(x)
#define	SIN(x)		    sin(x)
#define	COS(x)		    cos(x)
#define	ABS(x)		    fabs(x)
#define	ATAN2(x,y)	    atan2(x,y)
#define	SINH(x)		    sinh(x)
#define	COPYSIGN(x,y)	    copysign(x,y)
#endif

#if	USE_FAST
#undef	SIN
#undef	COS
#undef	ATAN2
#undef	SINH
#undef	EXP
#define	SIN(x)	    fast_sin(x)
#define	COS(x)	    fast_cos(x)
#define	ATAN2(x,y)  fast_atan2(x,y)
#define	SINH(x)	    fast_sinh(x)
#define	EXP(x)	    fast_exp(x)
#endif

#ifndef	M_PI
#define	M_PI	AS_REAL(3.14159265)
#endif

/* what is a position, units: micrometers */
#define	POSITION    REAL
/* what is an energy, units: kilo-electron-volts */
#define	ENERGY	    REAL
/* what is a direction, units: none */
#define	DIRECTION   REAL
/* what is an angle, unit: radians */
#define	ANGLE	    REAL
/* maximum (constant) distance moved per step, unit is mm */
// #define	MAXSTEP	    0.01

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

/* number of particle types */
#define	NPARTICLE   4
/* particle types and names */
enum particle_type { PROTON, ELECTRON, NEUTRON, ALPHA, UNKNOWN_P=-1 };
extern char *particle_name[NPARTICLE];

/* number of material types */
#define	NMATERIAL   7
/* we add these numbers to the material stored in order to mark tumor/exclusion regions */
#define	TUMOR_DELTA 0x0100
#define	EXCLU_DELTA 0x0200
/* this macro takes a material and returns it as the enum */
#define	GETMATER(n) (enum meterial_type ((n) & 0x00ff))
/* this macro checks if a material is tumurous or not */
#define	ISTUMOR(n)  ((n) & TUMOR_DELTA)
#define	ISEXCLU(n)  ((n) & EXCLU_DELTA)
/* material types and names */
enum material_type { WATER, AIR, ADIPOSETISSUEIRCP, A150TISSUE, MUSCLEWITHSUCROSE, B100BONE, OUTSIDE, UNKNOWN_M=-1 };
extern char *material_name[NMATERIAL];

/* the three types of events */
enum event_type { IONIZATION, ELASTIC, INELASTIC };

/* types of parameters passed on command line */
enum ParamType	{ P_INTEGER, P_BOOLEAN, P_REAL, P_STRING, P_LIST, P_FLAG, P_PATH, P_READ, P_WRIT };

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

typedef struct table {
    enum particle_type	P;	    /* type of particle for this table */
    enum material_type	M;	    /* type of particle for this table */
    int			np;	    /* number of points in the table */
    int			mysize;	    /* total size in bytes of this table */
    TBL_REAL		min, max;   /* limits of x variable */
    TBL_REAL		range;	    /* range of x variable */
    TBL_REAL		step;	    /* step size in the x variable */
    TBL_REAL		YV[1];	    /* array of positions, ith position is x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2]  */
} table;

typedef struct ParamDescriptor {
    const char		*key;
    enum ParamType	type;
    const char		*deflt;
    int			required;
} ParamDescriptor;

typedef struct arena {
    POSITION		MIN[3], WID[3], INVWID[3], MAX[3], CEN[3], RADIUS;
    int			N[3], NTOTAL;
    int			mysize;	    /* total size in bytes of this arena */
#if UNIFORM
    enum material_type	mat;
#else
    short		matG[1];
#endif
} arena;

typedef struct collector {
    int			N[3], NTOTAL, NTRACE, NSTEP, NMAX;
    int			mysize;	    /* total size in bytes of this collector */
    ENERGY		G[1];
} collector;

/* function in table.c */
extern table	*read_table_from_file(char *fname);
extern void	write_table_to_file(table *T, char *fname);
extern TBL_REAL	interpolate_table(table *T, TBL_REAL x);

/* functions in tabmgr.c */
extern void initialize_tables();
extern table *GetMaterialEnergyTable(enum particle_type T, enum material_type M);

/* in jack.c */

/* in particle.c */
extern void print_particle(char *str, particle *P);
extern particle *MakeParticle(particle *P, REAL ENE, REAL ANG, REAL AZI, unsigned int *seed);
extern void position_particle(ANGLE phi, ANGLE theta, particle *P);

/* functions in normal.c */
extern void fast_srand(unsigned int *seed);
extern REAL fastrand01(unsigned  int *seed);
extern REAL normal(REAL mean, REAL sigma, unsigned int *seed);

/* functions in dump.c */
extern void dump_good_bad_dose(collector *);
extern void dump_dose_text_file(collector *, char *);
extern void dump_dose_binary_file(collector *, char *);

/* functions in trace.c */
extern trace *alloc_trace();
extern void reset_trace();
extern void add_trace_point(trace *T, particle *P, int V, enum material_type M);
extern void dump_trace(FILE *f, trace *T, particle *P);
extern void free_trace(trace *T);

/* functions in proton.c */
extern void proton_event(particle *P, trace *T, unsigned int *seed);

/* functions in queu.c */
extern void initialize_queu(int n);
extern void reset_queu();
extern void queu_particle(particle *P);
extern particle *next_particle();
extern int empty();
extern int pull(particle *R, int n);

/* functions in collect.c */
extern arena *initialize_collect(char *, char *);
extern collector *allocate_collector();
extern void reset_collect(collector *);
extern void accumulate_collector(collector *, collector *);
extern void collect(collector *C, particle *P, trace *T);
extern void summarize_collect(collector *C, int do_plot);
extern int WhereAmI(POSITION x, POSITION y, POSITION z);

/* functions in param.c */
extern int ParseParams(int argc, char **argv, ParamDescriptor *table, void **destinations);

/* functions in delta.c */
extern void ProcessDeltaDirection(particle *P, ANGLE Theta, ANGLE Phi);
extern void CloneParticle(particle *Q, particle *P);
extern void GenerateElasticDistances(particle *P, unsigned int *seed);

/* functions in fluctuations.c */
extern void initialize_fluctuations();
extern TBL_REAL get_fluct(TBL_REAL E, enum material_type M, unsigned int *seed);

/* functions in fast.c */
extern REAL fast_sin(REAL);
extern REAL fast_cos(REAL);
extern REAL fast_atan2(REAL, REAL);
extern REAL fast_exp(REAL);
extern REAL fast_sinh(REAL);

/* functions in opencl.c */
extern void JackCL(int NMC, int bat_pltsel, int bat_devsel, int bat_wgropitems, int bat_batchitems,
            REAL EMEAN, REAL ESTDV, REAL SIGMA, REAL XYANGLE, REAL AZIMUTH, collector *C,
            int *g_mem_size, int *c_mem_size, int *l_mem_size, 
            int *c_ProtonWater_Energy_size, int *g_traceA_size, int *g_arena_size);


/* path to where the tables are */
extern char *TPATH;

/* control printouts */
extern int  VERBOSE;

/* the geometric arena */
extern arena	*ARENA;

/* size of an MRI pixel (should be in X/Y/Z)  */
extern REAL	PIXEL;

/* threshold for dumping voxels (proportion of maximum)  */
extern REAL	THR;

/* dimensions of the arena (why not using ARENA?) */
extern int	nx, ny, nz;
/* tables for dE/dx */
extern table  *ProtonWater_Energy;             
extern table  *ProtonAir_Energy;               
extern table  *ProtonAdiposeTissue_Energy;     
extern table  *ProtonATissue_Energy;           
extern table  *ProtonMuscleWithSucrose_Energy; 
extern table  *ProtonBBone_Energy;             
/* table for dE/dx fluctuations */
extern int	NE;	/* number of Energy values */
extern TBL_REAL	*ENRG;	/* energies */
extern int	NQ;	/* number of quantiles */
extern TBL_REAL	*PERC;	/* energies */
extern int	NM;	/* number of materials */
extern TBL_REAL	*FLUC;	/* table values */

/* maximum (constant) distance moved per step, unit is mm */
extern REAL MAXSTEP;
