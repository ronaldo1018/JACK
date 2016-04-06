#include    "jack.h"

#define  _g_seed  90862346
const REAL _g_scale1 = AS_REAL(1.0)/0x7FFF;
const REAL _g_scale2 = AS_REAL(2.0)/0x7FFF;

inline void fast_srand(unsigned int *seed) { 
    *seed = _g_seed;
} 
 
inline REAL fastrand01(unsigned int *seed) { 
    REAL a;
    *seed = (214013*(*seed)+2531011); 
    a = (((*seed)>>16)&0x7FFF) * _g_scale1;
    return(a);
} 
 
static inline REAL fastrand(unsigned int *seed) { 
    REAL a;
    *seed = (214013*(*seed)+2531011); 
    a = (((*seed)>>16)&0x7FFF) * _g_scale2 - AS_REAL(1.0);
    return(a);
} 

static inline REAL fast_log(REAL x) {
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

#if	USE_FAST
#undef	LOG(x)
#define	LOG(x)	    fast_log(x)
#endif

/* this is an alternative but slower version of the normal random number generator */
#if 0
REAL normal(REAL mean, REAL sigma) {
    REAL q, u, v, x, y;
/*  Generate P = (u,v) uniform in rectangle enclosing acceptance region */
    while (1) {
	u = fastrand01();
	v = fastrand01();
	v = (v - AS_REAL(0.5)) * AS_REAL(1.7156);
	x = u - AS_REAL(0.449871);
	y = ABS(v) + AS_REAL(0.386595);
	q = x*x + y*(AS_REAL(0.196)*y - AS_REAL(0.25472)*x);
/*  Accept P if inside inner ellipse */
	if (q < AS_REAL(0.27597)) {
	    x = v/u;
	    break;
	}
/*  Reject P if outside outer ellipse */
	if (q > AS_REAL(0.27846)) {
	    continue;
	}
/*  Reject P if outside acceptance region */
	if (v*v > AS_REAL(-4.0)*LOG(u)*u*u) {
	    continue;
	}
	x = v/u;
	break;
    }
    return(mean + x*sigma);
}
#else
REAL normal(REAL mean, REAL sigma, unsigned int *seed) {
    REAL x1, x2, w;
    do {
	x1 = fastrand(seed);
	x2 = fastrand(seed);
	w = x1 * x1 + x2 * x2;
    } while ( w >= AS_REAL(1.0) );
    w = SQRT( (-AS_REAL(2.0) * LOG( w ) ) / w );
    return(mean + x1*w*sigma);
}
#endif

