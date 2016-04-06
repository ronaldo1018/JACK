#include    "jack.h"

#define TPI_FLOAT    AS_REAL(6.28318531)
#define PI_FLOAT     AS_REAL(3.14159265)
#define PIBY2_FLOAT  AS_REAL(1.57079633)
#define	FOUR_PI	     AS_REAL(1.27323954)
#define	FOUR_PISQ    AS_REAL(0.40528473)
/* ********************************************************************** */
REAL fast_sin(REAL theta) {
    REAL    y;
    while (theta < -PI_FLOAT) { theta += TPI_FLOAT; }
    while (theta > PI_FLOAT) { theta -= TPI_FLOAT; }
    y = FOUR_PI * theta - FOUR_PISQ * theta * ABS(theta);
    y = AS_REAL(0.225) * (y * ABS(y) - y) + y;
    return(y);
}
/* ********************************************************************** */
REAL fast_cos(REAL theta) {
    return(fast_sin(PIBY2_FLOAT - theta));
}
/* ********************************************************************** */
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
/* ********************************************************************** */
static inline REAL fastpow2 (REAL p) {
    REAL offset = (p < AS_REAL(0.0)) ? AS_REAL(1.0) : AS_REAL(0.0);
    REAL clipp = (p < AS_REAL(-126.0)) ? AS_REAL(-126.0) : p;
    int w = clipp;
    REAL z = clipp - w + offset;
    union { unsigned int i; REAL f; } v = { (unsigned int) ( (1 << 23) * (clipp + AS_REAL(121.2740575) + AS_REAL(27.7280233) / (AS_REAL(4.84252568) - z) - AS_REAL(1.49012907) * z) ) };
    return v.f;
}
/* ********************************************************************** */
REAL fast_exp(REAL p) {
    return(fastpow2(p * AS_REAL(1.442695040)));
}
/* ********************************************************************** */
REAL fast_sinh(REAL p) {
    REAL   z;
    z = fastpow2(p * AS_REAL(1.442695040));
    return( AS_REAL(0.5) * (z - AS_REAL(1.0)/z) );
}
