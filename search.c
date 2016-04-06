#include    "jack.h"
#include    <pthread.h>

#define	    MAXTHREAD	32

char *particle_name[NPARTICLE] = { "PROTON", "ELECTRON", "NEUTRON", "ALPHA" };
char *material_name[NMATERIAL] = { "WATER", "AIR", "ADIPOSETISSUEIRCP", "A150TISSUE",
				   "MUSCLEWITHSUCROSE", "B100BONE", "OUTSIDE" };

/* intialize path to tables */
char *TPATH = NULL;
/* intialize path to MRI file */
char *MRIFIL = NULL;
/* intialize path to LOCATION (0/1) file */
char *LOCFIL = NULL;
/* the geometry */
arena *ARENA = NULL;
/* verbosity switch */
int VERBOSE = 0;
/* trace switch */
int SHOW_TRACE = 1;
/* MRI voxel size */
REAL PIXEL = AS_REAL(0.0);
/* threshold for dumping voxels (proportion of maximum)  */
REAL THR = AS_REAL(0.01);

/* this describes the options */
ParamDescriptor _params[] = {
/*  name    type        dflt		mandatory   */
{   "path", P_PATH,	"./tables",	0,	    },	    /* path to where the tables are */
{   "mri",  P_READ,	0,		1,	    },	    /* path to an MRI cube file */
{   "loc",  P_READ,	0,		1,	    },	    /* path to a tumor location file */
{   "pixel",P_REAL,	"600.0",	0,	    },	    /* MRI voxel size in microns */
{   "trace",P_FLAG,	"0",		0,	    },	    /* generate traces? */
{   "thread",P_INTEGER,	"0",		0,	    },	    /* generate traces? */
{   "ratio",P_REAL,	"0.5",		0,	    },	    /* print traces within this much of the maximum */
{   NULL,   P_REAL,	NULL,		0,	    }
};

typedef struct trial {
    ANGLE	    angle, azimuth;
    ENERGY	    energy, good, bad;
    int		    ngood, nbad;
} trial;

trial	    *DATA;
int	    ndata, NTHREAD, ITHR[MAXTHREAD];
pthread_t   PTHR[MAXTHREAD];

void Simulate(int begin, int end) {
    int		    i, j, idx;
    particle	    *P;
    trace	    *T;
/* initialize particle and trace */
    T = alloc_trace();
    P = MakeParticle(NULL, AS_REAL(0.0), AS_REAL(0.0), AS_REAL(0.0));
/* begin */
    for (j=begin; j<end; j++) {
/* create the particle */
	MakeParticle(P, DATA[j].energy*AS_REAL(1000.0), DATA[j].angle, DATA[j].azimuth);
	if (SHOW_TRACE) {
	    printf("particle %lx E=%f pos=(%f,%f,%f) dir=(%f,%f,%f)\n",
		   (long unsigned int) P, P->energy,
		   P->pos[CX], P->pos[CY], P->pos[CZ],
		   P->dir[CX], P->dir[CY], P->dir[CZ]);
	}
/* and trace it */
	reset_trace(T);
	while (P->energy > AS_REAL(0.0)) {
	    proton_event(P, T);
	}
/* print? */
	if (SHOW_TRACE) {
	    dump_trace(stdout, T, P);
	    getchar();
	}
/* then calculate the good and bad energy */
	for (i=1; i<T->np; i++) {
	    idx = T->data[i].where;
	    if (ARENA->is_tumor[idx] > 0) {
		DATA[j].good += (T->data[i-1].e - T->data[i].e);
		++DATA[j].ngood;
	    } else if (ARENA->is_tumor[idx] < 0) {
		DATA[j].bad += (T->data[i-1].e - T->data[i].e);
		++DATA[j].nbad;
	    }
	}
    }
}

void ThreadSimulate(int *n) {
    int	    j, begin, end;
/* determine which region to do */
    if (NTHREAD > 0) {
	j = ndata / NTHREAD;
	begin = *n * j;
	end = begin + j;
	if (*n == NTHREAD-1) end = ndata;
    } else {
	begin = 0;
	end = ndata;
    }
    Simulate(begin, end);
    pthread_exit(0);
}

int main(int argc, char **argv) {
    int		    i, j, k, tmp;
    double	    ratio, pixel;
    ANGLE	    angle, azimuth;
    ENERGY	    energy, avg_good, avg_bad, max_good, max_bad;
    void	    *values[16];
/* grab the parameters from the command line */
    tmp = 0;
    values[tmp++] = (void *) &TPATH;
    values[tmp++] = (void *) &MRIFIL;
    values[tmp++] = (void *) &LOCFIL;
    values[tmp++] = (void *) &pixel;
    values[tmp++] = (void *) &SHOW_TRACE;
    values[tmp++] = (void *) &NTHREAD;
    values[tmp++] = (void *) &ratio;
    if (ParseParams(argc-1, argv+1, _params, values) == -1) {
	fprintf(stderr, "flags were not parsed correctly!\n");
	exit(-1);
    }
    PIXEL = (REAL) pixel;
    if (VERBOSE) {
	printf("MRIFILE: %s\n", MRIFIL);
	printf("LOCFILE: %s\n", LOCFIL);
    }
/* any table or physics initializations go here */
    initialize_tables();
    if (VERBOSE) {
	printf("Tables from: %s\n", TPATH);
    }
/* initialize the particle collection domain */
    ARENA = initialize_collect(MRIFIL, LOCFIL);
    if (VERBOSE) {
	j = k = 0;
	for (i=0; i<ARENA->NTOTAL; i++) {
	    if (ARENA->is_tumor[i] > 0) ++j;
	    if (ARENA->is_tumor[i] < 0) ++k;
	}
	printf("Allocated Arena: %d x %d x %d (good=%d bad=%d)\n",
	       ARENA->N[CX], ARENA->N[CY], ARENA->N[CZ], j, k);
    }
/* and begin */
    DATA = (trial *) malloc(61*61*61*sizeof(trial));
    ndata = 0;
    for (angle=AS_REAL(0.0); angle<AS_REAL(360.0); angle+=AS_REAL(6.0)) {		    /* 60 */
	for (azimuth=AS_REAL(0.0); azimuth<=AS_REAL(180.0); azimuth+=AS_REAL(3.0)) {	    /* 61 */
	    for (energy=AS_REAL(50.0); energy<=AS_REAL(200.0); energy+=AS_REAL(2.5)) {	    /* 61 */
		DATA[ndata].angle = angle;
		DATA[ndata].azimuth = azimuth;
		DATA[ndata].energy = energy;
		DATA[ndata].good = AS_REAL(0.0);
		DATA[ndata].bad = AS_REAL(0.0);
		DATA[ndata].ngood = 0;
		DATA[ndata].nbad = 0;
		++ndata;
	    }
	}
    }
/* do the simulations */
    if (NTHREAD > 1) {
	assert(NTHREAD < MAXTHREAD);
	for (i=0; i<NTHREAD; i++) {
	    ITHR[i] = i;
	    pthread_create(&(PTHR[i]), NULL, (void *) &ThreadSimulate, (void *) &(ITHR[i]));
	}
	for (i=0; i<NTHREAD; i++) {
	    pthread_join(PTHR[i], NULL);
	}
    } else {
	Simulate(0, ndata);
    }
/* compute the average and maximum good dose */
    i = 0;
    avg_good = AS_REAL(0.0);
    max_good = AS_REAL(0.0);
    for (j=0; j<ndata; j++) {
	if (DATA[j].ngood) {
	    avg_good += DATA[j].good;
	    if (DATA[j].good > max_good) max_good = DATA[j].good;
	    ++i;
	}
    }
    if (i) avg_good /= i;
/* compute the average bad dose */
    i = 0;
    avg_bad = AS_REAL(0.0);
    max_bad = AS_REAL(0.0);
    for (j=0; j<ndata; j++) {
	if (DATA[j].nbad) {
	    avg_bad += DATA[j].bad;
	    if (DATA[j].bad > max_bad) max_bad = DATA[j].bad;
	    ++i;
	}
    }
    if (i) avg_bad /= i;
/* print out treatments above the good average and below the bad */
    for (j=0; j<ndata; j++) {
	if (DATA[j].good <= ratio*max_good) continue;
	if (DATA[j].bad  >= (AS_REAL(1.0)-ratio)*max_bad) continue;
	printf("%8.3f %8.3f %8.3f %11.4e %11.4e %5d %5d\n",
	       DATA[j].angle, DATA[j].azimuth, DATA[j].energy,
	       DATA[j].good/max_good, DATA[j].bad/max_bad,
	       DATA[j].ngood, DATA[j].nbad);
    }
    return(0);
}
