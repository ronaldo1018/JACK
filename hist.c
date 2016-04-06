#include "jack.h"

int main(int argc, char **argv) {
    int			i, j, nrow, ncol, npln, ntot;
    char		*fc, *temp;
    FILE		*F;
    short		*v;
    double		valuedensity[8], valueCT[8], density, deltaCT, deltaDensity;
/* these arrays hold the outputs that will eventually be dumped to binary files */
    float		*D;
    enum material_type	*M;
    
    if (argc < 3) {
	fprintf(stderr, "expect: <input file> <material file>\n");
	exit(-1);
    }
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
    printf("assuming file is %d x %d x %d\n", nrow, ncol, npln);

    ntot = nrow*ncol*npln;
    v = (short *) malloc(ntot*sizeof(short));
    M = (enum material_type *) malloc(ntot*sizeof(enum material_type));
    D = (float *) malloc(ntot*sizeof(float));
    if (! v || ! M || ! D) {
	fprintf(stderr, "failed to malloc %d shorts/enums/floats\n", ntot);
	exit(-1);
    }
    i = (int) fread(v, sizeof(short), ntot, F);
    if (i != ntot) {
      fprintf(stderr, "found only %d shorts while expecting %d in file\n", i, ntot);
	exit(-1);
    }
    fclose(F);

    printf("ntot=%d\n", ntot);
/* these numbers are taken from the MRI example in the Geant directory */
    valueCT[0] = -5000.0;
    valueCT[1] = -1000.0;
    valueCT[2] = -400.0;
    valueCT[3] = -150.0;
    valueCT[4] = 100.0;
    valueCT[5] = 300.0;
    valueCT[6] = 2000.0;
    valueCT[7] = 4927.0;
/* ditto */
    valuedensity[0] = 0.0;
    valuedensity[1] = 0.0;
    valuedensity[2] = 0.602;
    valuedensity[3] = 0.924;
    valuedensity[4] = 1.075;
    valuedensity[5] = 1.145;
    valuedensity[6] = 1.856;
    valuedensity[7] = 3.379;
/* go through the file and translate each density to the material enum */
    for (i=0; i<ntot; i++) {
	density = 0.0;
	for (j=1; j<8; j++) {
	    if (v[i] >= valueCT[j-1] && v[i] < valueCT[j]) {
		deltaCT = valueCT[j] - valueCT[j-1];
		deltaDensity = valuedensity[j] - valuedensity[j-1];
/* interpolating linearly */
		D[i] = density = valuedensity[j] - ((valueCT[j] - v[i])*deltaDensity/deltaCT );
		break;
	    }
	}
/* these numbers are also taken from the MRI example in the Geant directory */
	if (density >= 0.0 && density < 0.207) {	    /* Air */
	    M[i] = AIR;
	} else if (density >= 0.207 && density < 0.481) {   /* iLung -> Air */
	    M[i] = AIR;
	} else if (density >= 0.481 && density < 0.919) {    /* eLung -> Air */
	    M[i] = AIR;
	} else if (density >= 0.919 && density < 0.979) {    /* Adipose */
	    M[i] = ADIPOSETISSUEIRCP;
	} else if (density >= 0.979 && density < 1.004) {    /* Breast */
	    M[i] = A150TISSUE;
	} else if (density >= 1.004 && density < 1.043) {    /* Phantom */
	    M[i] = WATER;
	} else if (density >= 1.043 && density < 1.109) {    /* Liver */
	    M[i] = A150TISSUE;
	} else if (density >= 1.109 && density < 1.113) {    /* Muscle */
	    M[i] = MUSCLEWITHSUCROSE;
	} else if (density >= 1.113 && density < 1.496) {    /* TrabecularBone */
	    M[i] = B100BONE;
	} else if (density >= 1.496) {			    /* DenseBone */
	    M[i] = B100BONE;
	} else {
	    printf("WOW density=%f pixel=%i\n", density, v[i]);
	}
    }
/* write material file */
    fc = argv[2];
    F = fopen(fc,"wb");
    if (! F) {
	fprintf(stderr, "unable to open: %s\n", fc);
	exit(-1);
    }
    fwrite(M, sizeof(enum material_type), ntot, F);
    fclose(F);
/* done */
    return(0);
}
