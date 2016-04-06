#include    "jack.h"

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

void ProcessDeltaDirection(particle *P, ANGLE Phi, ANGLE Theta) {
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

void CloneParticle(particle *Q, particle *P) {
    P->type = Q->type;
    P->energy = Q->energy;
    P->pos[CX] = Q->pos[CX];
    P->pos[CY] = Q->pos[CY];
    P->pos[CZ] = Q->pos[CZ];
    P->dir[CX] = Q->dir[CX];
    P->dir[CY] = Q->dir[CY];
    P->dir[CZ] = Q->dir[CZ];
    P->to_elastic = Q->to_elastic;
    P->to_inelastic = Q->to_inelastic;
}

void GenerateElasticDistances(particle *P, unsigned int *seed) {
/* this is completely bogus... */
    P->to_elastic = 1e5*normal(4.0, 1.0, seed);
    P->to_inelastic = 2e5*normal(4.0, 1.0, seed);
}
