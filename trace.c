#include    "jack.h"

trace *alloc_trace() {
    trace   *T;
    T = (trace *) malloc(sizeof(trace));
    T->np = 0;
#if ! (FULL_TRACE || MIN_TRACE)
    T->last = AS_REAL(0.0);
#endif
    return(T);
}

void free_trace(trace *T) {
    assert(T);
    free((char *) T);
}

void reset_trace(trace *T) {
    T->np = 0;
#if ! (FULL_TRACE || MIN_TRACE)
    T->last = AS_REAL(0.0);
#endif
}

void dump_trace(FILE *f, trace *T, particle *P) {
    int	    i, ix, iy, iz;
    fprintf(f, "Trace %lx for particle %lx (%s) (np=%d)\n",
	    (long unsigned int) T, (long unsigned int) P, particle_name[P->type], T->np);
    for (i=0; i<T->np; i++) {
#if FULL_TRACE
	fprintf(f, "%4d @ %9.1f,%9.1f,%9.1f E=%9.3e dE=%9.3e %d/%s\n",
		i, (double) T->data[i].x, (double) T->data[i].y, (double) T->data[i].z,
		(double) T->data[i].e, i>0 ? (double) T->data[i-1].e - T->data[i].e : 0.0,
		T->data[i].where, material_name[T->data[i].m]);
#else
	split(T->data[i].where,ix,iy,iz,ARENA->N[CX],ARENA->N[CY],ARENA->N[CZ]);
	fprintf(f, "%4d @ %3d,%3d,%3d dE=%9.3e\n", i, ix, iy, iz, (double) T->data[i].e);
#endif	/* FULL_TRACE */
    }
}

void dump_trace_header(FILE *f, trace *T, particle *P) {
    int	    i;
    fprintf(f, "Trace %lx for particle %lx (%s) (np=%d)\n",
	    (long unsigned int) T, (long unsigned int) P, particle_name[P->type], T->np);
    for (i=0; i<5; i++) {
#if FULL_TRACE
	fprintf(f, "%4d @ %9.1f,%9.1f,%9.1f E=%9.3e dE=%9.3e %d/%s\n",
		i, (double) T->data[i].x, (double) T->data[i].y, (double) T->data[i].z,
		(double) T->data[i].e, i>0 ? (double) T->data[i-1].e - T->data[i].e : 0.0,
		T->data[i].where, material_name[T->data[i].m]);
#else
	fprintf(f, "%4d @ %07d E=%9.3e dE=%9.3e\n",
		i, T->data[i].where, (double) T->data[i].e,
		i>0 ? (double) T->data[i-1].e - T->data[i].e : 0.0);
#endif	/* FULL_TRACE */
    }
    for (i=T->np-5; i<T->np; i++) {
#if FULL_TRACE
	fprintf(f, "%4d @ %9.1f,%9.1f,%9.1f E=%9.3e dE=%9.3e %d/%s\n",
		i, (double) T->data[i].x, (double) T->data[i].y, (double) T->data[i].z,
		(double) T->data[i].e, i>0 ? (double) T->data[i-1].e - T->data[i].e : 0.0,
		T->data[i].where, material_name[T->data[i].m]);
#else
	fprintf(f, "%4d @ %07d E=%9.3e dE=%9.3e\n",
		i, T->data[i].where, (double) T->data[i].e,
		i>0 ? (double) T->data[i-1].e - T->data[i].e : 0.0);
#endif	/* FULL_TRACE */
    }
}

void add_trace_point(trace *T, particle *P, int V, enum material_type M) {
/* if a trace grows too large, kill the particle */
    if (T->np >= MAXTRACE) {
#if TRC_DEBUG
	if (VERBOSE) { fprintf(stderr, "terminating overflowing trace\n"); dump_trace_header(stderr, T, P); }
#endif
	P->energy = AS_REAL(0.0); return;
    }
/* if we are currently outside then... */
    if (M == OUTSIDE) {
	if (T->np) {	/* if particle has been in, kill it */
#if TRC_DEBUG
	    if (VERBOSE) { fprintf(stderr, "terminating trace leaving region\n"); dump_trace_header(stderr, T, P); }
#endif
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
