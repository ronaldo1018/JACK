#include    "jack.h"

/* this should be done much better, hacked to get things going now */
table  *ProtonWater_Energy;
table  *ProtonAir_Energy;
table  *ProtonAdiposeTissue_Energy;
table  *ProtonATissue_Energy;
table  *ProtonMuscleWithSucrose_Energy;
table  *ProtonBBone_Energy;

void initialize_tables() {
    char    fname[MAXLIN];
/* read proton-energy-water file */
    sprintf(fname, "%s/%s_%s_energy.tbl", TPATH, particle_name[PROTON], material_name[WATER]);
    ProtonWater_Energy = read_table_from_file(fname);
    if (! ProtonWater_Energy) {
	fprintf(stderr, "cannot find table in \"%s\"!\n", fname);
	exit(-1);
    }
    assert(ProtonWater_Energy->P == PROTON);
    assert(ProtonWater_Energy->M == WATER);
/* read proton-energy-air file */
    sprintf(fname, "%s/%s_%s_energy.tbl", TPATH, particle_name[PROTON], material_name[AIR]);
    ProtonAir_Energy = read_table_from_file(fname);
    if (! ProtonAir_Energy) {
      printf("Cannot find table in %s\n", fname);
	fprintf(stderr, "cannot find table in \"%s\"!\n", fname);
	exit(-1);
    }
    assert(ProtonWater_Energy->P == PROTON);
    assert(ProtonWater_Energy->M == WATER);
/* read proton-energy-adipos file */
    sprintf(fname, "%s/%s_%s_energy.tbl", TPATH, particle_name[PROTON], material_name[ADIPOSETISSUEIRCP]);
    ProtonAdiposeTissue_Energy = read_table_from_file(fname);
    if (! ProtonAdiposeTissue_Energy) {
      printf("Cannot find table in %s\n", fname);
	fprintf(stderr, "cannot find table in \"%s\"!\n", fname);
	exit(-1);
    }
    assert(ProtonWater_Energy->P == PROTON);
    assert(ProtonWater_Energy->M == WATER);
/* read proton-energy-a150 file */
    sprintf(fname, "%s/%s_%s_energy.tbl", TPATH, particle_name[PROTON], material_name[A150TISSUE]);
    ProtonATissue_Energy = read_table_from_file(fname);
    if (! ProtonATissue_Energy) {
      printf("Cannot find table in %s\n", fname);
	fprintf(stderr, "cannot find table in \"%s\"!\n", fname);
	exit(-1);
    }
    assert(ProtonWater_Energy->P == PROTON);
    assert(ProtonWater_Energy->M == WATER);
/* read proton-energy-muscle file */
    sprintf(fname, "%s/%s_%s_energy.tbl", TPATH, particle_name[PROTON], material_name[MUSCLEWITHSUCROSE]);
    ProtonMuscleWithSucrose_Energy = read_table_from_file(fname);
    if (! ProtonMuscleWithSucrose_Energy) {
      printf("Cannot find table in %s\n", fname);
	fprintf(stderr, "cannot find table in \"%s\"!\n", fname);
	exit(-1);
    }
    assert(ProtonWater_Energy->P == PROTON);
    assert(ProtonWater_Energy->M == WATER);
/* read proton-energy-bone file */
    sprintf(fname, "%s/%s_%s_energy.tbl", TPATH, particle_name[PROTON], material_name[B100BONE]);
    ProtonBBone_Energy = read_table_from_file(fname);
    if (! ProtonBBone_Energy) {
      printf("Cannot find table in %s\n", fname);
	fprintf(stderr, "cannot find table in \"%s\"!\n", fname);
	exit(-1);
    }
    assert(ProtonWater_Energy->P == PROTON);
    assert(ProtonWater_Energy->M == WATER);

}

table *GetMaterialEnergyTable(enum particle_type T, enum material_type M) {
    switch(M) {
	case WATER:
	    return ProtonWater_Energy;
	    break;
	case AIR:
	    return ProtonAir_Energy;
	    break;
	case ADIPOSETISSUEIRCP:
	    return ProtonAdiposeTissue_Energy;
	    break;
	case A150TISSUE:
	    return ProtonATissue_Energy;
	    break;
	case MUSCLEWITHSUCROSE:
	    return ProtonMuscleWithSucrose_Energy;
	    break;
	case B100BONE:
	    return ProtonBBone_Energy;
	    break;
	case OUTSIDE:
	case UNKNOWN_M:
	    assert(1);
	    break;
    }
    return(ProtonWater_Energy);
}
