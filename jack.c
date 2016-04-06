//---------------------------------------------------------------------
// File:    jack.c
// Date:    2012.10.01
// Create:  
// Modify:
//
// Description:
//  1. Provide command mode CMD_MODE for quickly changes parameters.
//  2. When CMD_MODE set '1'
//       ./jack c [Num]   => run Num traces using cpu single thread
//   eg. ./jack c 1000000 => run 1000000 traces
//  3. When CMD_MODE set '0' : run OpenCL code
//         a. Select Target Platform [platform]:  select which Platform number you need
//         b. Select Target Device [device]  :  select which Device number you want
//         c. [itmes/work group] : IW
//         d. [items/batch] : IB
//            Total batches = total traces / IB,
//            Workloads are shared by (IB/IW) Work Groups 
//     ./jack o [total traces] [platform] [device] [local itmes/work group] [items/batch]
//   eg.
//     ./jack o 1000000 0 0 250 20000
//          1M traces run on device 0 of platfom 0
//          20000 traces per batch, there would run 1M/20k = 50 batches
//          each work group run 250 traces, so there would be 
//           20k/250 = 80 work groups run on compute units.
//  4. ./jack q : query OpenCL hardware info.  
//---------------------------------------------------------------------

#include    "jack.h"
#include    <pthread.h>
#include    <sys/time.h>

#define	    CMD_MODE   1
#define     DEBUG   1
/* maximum (constant) distance moved per step, unit is mm */
REAL MAXSTEP=0.01;

struct timeval gstart, gend; 
void srand48(long int seedval);

char *particle_name[NPARTICLE] = { "PROTON", "ELECTRON", "NEUTRON", "ALPHA" };
char *material_name[NMATERIAL] = { "WATER", "AIR", "ADIPOSETISSUEIRCP", "A150TISSUE",
				   "MUSCLEWITHSUCROSE", "B100BONE", "OUTSIDE" };
/* intialize path to tables */
char *TPATH = NULL;
/* intialize path to MRI file */
char *MRIFIL = NULL;
/* intialize path to LOCATION (0/1) file */
char *LOCFIL = NULL;
/* intialize path to binary output file */
char *BOUTFIL = NULL;
/* intialize path to text output file */
char *TOUTFIL = NULL;
/* verbosity switch */
int VERBOSE = 1;
/* trace switch */
int TRACE = 0;
/* the geometry */
arena* ARENA = NULL;
/* MRI voxel size */
REAL PIXEL = AS_REAL(0.0);
/* threshold for dumping voxels (proportion of maximum)  */
REAL THR = AS_REAL(0.01);
/* this describes the options */
ParamDescriptor _params[] = {
/*	name		type		dflt		mandatory   */
    {   "nmc",		P_INTEGER,	"10000",		0,	    },
    {   "mean",		P_REAL,		"100",		0,	    },	    /* 100MeV */
    {   "stdv",		P_REAL,		"2",		0,	    },	    /* 2MeV */
    {   "sigma",	P_REAL,		"0.01",		0,	    },	    /* initial postion sigma */
    {   "path",		P_PATH,		"./tables",	0,	    },	    /* path to where the tables are */
    {   "mri",		P_READ,		0,		0,	    },	    /* path to an MRI cube file */
    {   "trace",	P_FLAG,		"0",		0,	    },	    /* generate traces? */
    {   "plot",		P_FLAG,		"0",		0,	    },	    /* make plot file? */
    {   "angle",	P_REAL,		"0.0",		0,	    },	    /* angle in XY plane */
    {   "azimuth",	P_REAL,		"90.0",		0,	    },	    /* angle to Z axis */
    {   "dumpt",	P_WRIT,		0,		0,	    },	    /* path to text output dump file */
    {   "dumpb",	P_WRIT,		0,		0,	    },	    /* path to binary output dump file */
    {   "loc",		P_READ,		0,		0,	    },	    /* path to a tumor location file */
    {   "quiet",	P_FLAG,		0,		0,	    },	    /* be verbose */
    {   "pixel",	P_REAL,		"600.0",	0,	    },	    /* MRI voxel size in microns */
    {   "thr",		P_REAL,		"0.01",		0,	    },	    /* threshold for dumping text voxels (as proportion of maximum) */
    {   "thread",	P_INTEGER,	"0",		0,	    },	    /* generate traces? */
    {   "platform",	P_INTEGER,	"-1",		0,	    },	    /* opencl platform */
    {   "device",	P_INTEGER,	"0",		0,	    },	    /* opencl device */
    {   "numworkgrp",	P_INTEGER,	"0",		0,	    },	    /* opencl number of workgroups */
    {   "numbatch",	P_INTEGER,	"1000",		0,	    },	    /* opencl number of batcher */
    {   "collect",	P_INTEGER,	"100",		0,	    },	    /* collect grid points */
    {   NULL,		P_REAL,		NULL,		0,	    }
};

void MonteCarlo(collector *C, int NMC, REAL emean, REAL estdv, REAL sigma, REAL angle, REAL azimuth) {
    particle	    *P;
    trace	    *T;
    int		    NF;
    unsigned int seed;
    NF = 0;
#if DEBUG
    printf("entering MonteCarlo for %d samples and storing into %lx (%d)\n",
	   NMC, (long unsigned int) C, C->NTOTAL);
#endif
/* create the trace and particle */
    fast_srand(&seed);
    T = alloc_trace();
    P = (particle *) malloc(sizeof(particle));
/* the Monte Carlo loop */
//what is sigma??
    if (sigma > AS_REAL(0.0)) {
	while (NF++ < NMC) {
	    MakeParticle(P, normal(emean, estdv,&seed), angle+normal(AS_REAL(0.0), sigma,&seed), azimuth+normal(AS_REAL(0.0), sigma,&seed), &seed);
	    reset_trace(T);
	    while (P->energy > AS_REAL(0.0)) {
		proton_event(P, T, &seed); //stepping()
	    }
	    if (TRACE) dump_trace(stdout, T, P);
	    collect(C, P, T); //Hits collection
	}
    } else {
	while (NF++ < NMC) {
	    MakeParticle(P, normal(emean, estdv,&seed), angle, azimuth, &seed);
	    reset_trace(T);
	    while (P->energy > AS_REAL(0.0)) {
		proton_event(P, T, &seed);
	    }
	    if (TRACE) dump_trace(stdout, T, P);
	    collect(C, P, T);
	}
    }
}

int	    NTHREAD, NMC, *ITHR;
collector   *C, **CTHR;
pthread_t   *PTHR;
double	    EMEAN, ESTDV, SIGMA, XYANGLE, AZIMUTH;

void ThreadSimulate(int *n) {
    MonteCarlo(CTHR[*n], NMC/NTHREAD, (REAL) EMEAN, (REAL) ESTDV, (REAL) SIGMA, (REAL) XYANGLE, (REAL) AZIMUTH);
    pthread_exit(0);
}

int main(int argc, char **argv) {
    //printf("%d %c\n",argc,argv[0][0]);
    int		    i, DOPLOT, quiet, tmp;
    double	    pixel, thr;
    void	    *values[32];
/* opencl options */
    int		    bat_pltsel=0, bat_devsel=0, bat_wgropitems=0, bat_batchitems=0;
    int script_mode = 0;
    int g_mem_size = -1, c_mem_size = -1, l_mem_size = -1, c_ProtonWater_Energy_size = -1, g_traceA_size = -1, g_arena_size = -1;
/* grab the parameters from the command line */
    tmp = 0;
    values[tmp++] = (void *) &NMC;
    values[tmp++] = (void *) &EMEAN;
    values[tmp++] = (void *) &ESTDV;
    values[tmp++] = (void *) &SIGMA;
    values[tmp++] = (void *) &TPATH;
    values[tmp++] = (void *) &MRIFIL;
    values[tmp++] = (void *) &TRACE;
    values[tmp++] = (void *) &DOPLOT;
    values[tmp++] = (void *) &XYANGLE;
    values[tmp++] = (void *) &AZIMUTH;
    values[tmp++] = (void *) &TOUTFIL;
    values[tmp++] = (void *) &BOUTFIL;
    values[tmp++] = (void *) &LOCFIL;
    values[tmp++] = (void *) &quiet;
    values[tmp++] = (void *) &pixel;
    values[tmp++] = (void *) &thr;
    values[tmp++] = (void *) &NTHREAD;
    values[tmp++] = (void *) &bat_pltsel;
    values[tmp++] = (void *) &bat_devsel;
    values[tmp++] = (void *) &bat_wgropitems;
    values[tmp++] = (void *) &bat_batchitems;
    values[tmp++] = (void *) &nx;
    if (ParseParams(argc-1, argv+1, _params, values) == -1) {
	fprintf(stderr, "flags were not parsed correctly!\n");
	exit(-1);
    }    
#if CMD_MODE        
    if (argc >= 2)
    {
        char *intrace = argv[1];
        if (intrace[0] == 'c')
        {
            NMC = atoi(argv[2]);           
            bat_pltsel = -1;
        }        
        if (intrace[0] == 'q')
        {
            quiet = 1;
            bat_pltsel = 99;
        }        
        if (intrace[0] == 'o')
        {
            NMC = atoi(argv[2]);
            bat_pltsel = atoi(argv[3]);
            bat_devsel = atoi(argv[4]);
            bat_wgropitems = atoi(argv[5]);
            bat_batchitems = atoi(argv[6]);
        }
        if (intrace[0] == 's' && intrace[1] == 'o')  //script opencl
        {
            script_mode = 2;
            quiet = 0;
            NMC = atoi(argv[2]);
            NTHREAD = -1;
            bat_pltsel = atoi(argv[3]);
            bat_devsel = atoi(argv[4]);
            bat_wgropitems = atoi(argv[5]);
            bat_batchitems = atoi(argv[6]);
            nx = atoi(argv[7]);
            MAXSTEP = atof(argv[8]);
            EMEAN = atof(argv[9]);
        }
        
        if (intrace[0] == 's' && intrace[1] == 'c')  //script cpu
        {
            script_mode = 1;
            quiet = 0;
            NMC = atoi(argv[2]);
            NTHREAD = atoi(argv[3]);
            bat_pltsel = -1;
            bat_devsel = -1;
            bat_wgropitems = -1;
            bat_batchitems = -1;
            nx = atoi(argv[4]);
            MAXSTEP = atof(argv[5]);
            EMEAN = atof(argv[6]);
        }
    }
#endif

    FILE *csv_file = NULL;

    if (script_mode)
    { 	
        csv_file = fopen("jack_script.csv", "r");
        if (csv_file == NULL)
        {
            csv_file = fopen("jack_script.csv", "a");
            fprintf(csv_file, "NTHREAD, PLATFORM, DEVICE, Items_WG, Items_batch, NMC, dC, dS(mm), E(MeV), Time(sec), Global_MEM(B), Const_MEM(B), Local_MEM(B), Trace Size(B), Arena Size(B), WaterTable Size(B)\n");  
        }
        else
        {
            csv_file = fopen("jack_script.csv", "a");
        }
    }
	ny = nz = nx;
    PIXEL = (REAL) pixel;
    THR = (REAL) thr;
    VERBOSE = ! quiet;  
    if (VERBOSE) {
        printf("Number of Monte-Carlo trials (-nmc): %d\n", NMC);
        printf("Incident Particle Energy (-mean,-stdv,-sigma): N(%g,%g,%g) direction (-angle,-azimuth) (%g,%g)\n",
               EMEAN, ESTDV, SIGMA, XYANGLE, AZIMUTH);
        printf("Tables Read From (-path): %s\n", TPATH);
        if (MRIFIL) printf("MRI Read From (-mri): %s with pixel size %.2f\n", MRIFIL, PIXEL);
        if (LOCFIL) printf("Location Read From (-loc): %s\n", LOCFIL);
        if (BOUTFIL) printf("Output binary dose will be written to (-dump): %s\n", BOUTFIL);
        if (TOUTFIL) printf("Output text dose will be written to (-dump): %s\n", TOUTFIL);
        printf("Dump Traces (-trace): %s\n", TRACE ? "yes" : "no");
        printf("Make Plot File (-plot): %s\n", DOPLOT ? "yes" : "no");
        printf("OpenCL Device: %d Platform: %d WorkGroupItems: %d BatchItems: %d\n", bat_pltsel, bat_devsel, bat_wgropitems, bat_batchitems);
    }
/* initialize */
    initialize_queu(NQUEU);
/* initialize the particle collection domain */
    ARENA = initialize_collect(MRIFIL, LOCFIL);
/* any table or physics initializations go here */
    initialize_tables();
    initialize_fluctuations();
/* create the global collector */
    C = allocate_collector();

    clock_t looper = clock();
    gettimeofday(&gstart, NULL); 
//    printf("bat_pltsel = %d\n", bat_pltsel);
        
/* if opencl device was specified, then run the opencl version */
    if (bat_pltsel != -1) {
        JackCL(NMC, bat_pltsel, bat_devsel, bat_wgropitems, bat_batchitems, EMEAN, ESTDV, SIGMA, XYANGLE, AZIMUTH, C,
               &g_mem_size, &c_mem_size, &l_mem_size, &c_ProtonWater_Energy_size, &g_traceA_size, &g_arena_size);
    } else if (NTHREAD < 2) {
        MonteCarlo(C, NMC, (REAL) EMEAN, (REAL) ESTDV, (REAL) SIGMA, (REAL) XYANGLE, (REAL) AZIMUTH);
    } else {
        ITHR = (int *) malloc(NTHREAD*sizeof(int));
        CTHR = (collector **) malloc(NTHREAD*sizeof(collector *));
        PTHR = (pthread_t *) malloc(NTHREAD*sizeof(pthread_t));
        for (i=0; i<NTHREAD; i++) {
            CTHR[i] = allocate_collector();
        }
        for (i=0; i<NTHREAD; i++) {
            ITHR[i] = i;
            pthread_create(&(PTHR[i]), NULL, (void *) &ThreadSimulate, (void *) &(ITHR[i]));
        }
        for (i=0; i<NTHREAD; i++) {
            pthread_join(PTHR[i], NULL);
        }
        for (i=0; i<NTHREAD; i++) accumulate_collector(C, CTHR[i]);
    }
    if (bat_pltsel == 99) return 0;
    clock_t end = clock();
    gettimeofday(&gend, NULL); 
    float delta = ((gend.tv_sec  - gstart.tv_sec) * 1000000u +           
		   gend.tv_usec - gstart.tv_usec) / 1.e6;   
    float second_time = (float) (end - looper) / CLOCKS_PER_SEC;

/* dump files */
    if (TOUTFIL) dump_dose_text_file(C, TOUTFIL);
    if (BOUTFIL) dump_dose_binary_file(C, BOUTFIL);
    if (LOCFIL) dump_good_bad_dose(C);
/* print the results */
    summarize_collect(C, DOPLOT);
    
    if (script_mode)
    { 	
        fprintf(csv_file, "%d, %d, %d, %d, %d, %d, %d, %.2f, %f, %f, %d, %d, %d, %d, %d, %d\n", 
            NTHREAD, bat_pltsel, bat_devsel, bat_wgropitems, bat_batchitems, NMC, nx, MAXSTEP, EMEAN, delta, g_mem_size, c_mem_size, l_mem_size, g_traceA_size, g_arena_size, c_ProtonWater_Energy_size);  
    }
     
    printf("\nSimulation Loop Time: %f seconds \n", second_time);
    printf("Simulation Loop Time (gettimeofday): %f seconds \n", delta);
    return(0);
}
