#include    "jack.h"

static REAL	PHI_SCALE[NMATERIAL] = {
/* WATER */		22.634e-3,
/* AIR */		0.7143e-3,
/* ADIPOSETISSUEIRCP */	23.134e-3,
/* A150TISSUE */	19.877e-3,
/* MUSCLEWITHSUCROSE */	22.237e-3,
/* B100BONE */		23.611e-3,
/* OUTSIDE */		0.0		};

static	ANGLE	SamplePhi(particle *P, enum material_type M, unsigned int *seed) {
    ENERGY  m = PHI_SCALE[(int) M];
    if (P->energy < 1.0) {
	return(0.0);
    } else {
	return( EXP(normal(AS_REAL(0.0), AS_REAL(0.693),seed))*m/P->energy );
    }
}

static	ANGLE	SampleTheta(particle *P, enum material_type M, unsigned int *seed) {
    return( AS_REAL(2.0) * M_PI * fastrand01(seed) );
}

static void proton_ionization_event(arena *A, particle *P, enum material_type M, trace *TRACE, unsigned int *seed) {
/* which table do I use for the delta-energy */
    table	*T;
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
    T = GetMaterialEnergyTable(P->type, M);
    dE = (ENERGY) interpolate_table(T, (TBL_REAL) P->energy); //search the table by particle energy
/* look up the statistics of the dE/dx fluctuations in table */
    dE += get_fluct(P->energy, M, seed);
/* and multiply by the step in mm to get the actual dE in MeV */
    dE *= MAXSTEP;
/* update particle energy */
    assert(dE >= AS_REAL(0.0));
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
    new = WhereAmI(P->pos[0], P->pos[1], P->pos[2]);
/* if not the same position as before, or the particle has died, update the trace */
    if (new != P->where || P->energy <= AS_REAL(0.0)) {
        add_trace_point(TRACE, P, P->where, M);
        P->eloss = AS_REAL(0.0);
        P->where = new;
    }
#endif
}

void proton_event(particle *P, trace *TRACE, unsigned int *seed) {
/* the voxel we are in */
    int			    V;
/* what material am I in */
    enum material_type	    M;
#if DEBUG
    if (VERBOSE) {
	print_particle("event", P);
    }
#endif
#if MIN_TRACE                              
    V = P->where; //index of the porticle position
#if UNIFORM
    M = (V == -1 ? OUTSIDE : ARENA->mat);
#else
    M = (V == -1 ? OUTSIDE : ARENA->matG[V]);
#endif
    proton_ionization_event(ARENA, P, M, TRACE, seed);
#else
    V = WhereAmI(P->pos[0], P->pos[1], P->pos[2]);
#if UNIFORM
    M = (V == -1 ? OUTSIDE : ARENA->mat);
#else
    M = (V == -1 ? OUTSIDE : ARENA->matG[V]);
#endif
    proton_ionization_event(ARENA, P, M, TRACE, seed);
    add_trace_point(TRACE, P, V, M);
#endif
}
