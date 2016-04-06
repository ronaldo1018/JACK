/* to help make the code clearer */
#define	CX   0	    /* columns */
#define	CY   1	    /* rows */
#define	CZ   2	    /* planes */

#define	    getidx(x,y,z,nx,ny,nz)    ((nx)*(ny)*(z) + (nx)*(y) + (x))

#define	    split(idx,x,y,z,nx,ny,nz)    {\
		int tmp = idx % ( (nx)*(ny) ); \
		(x) = tmp % (nx); \
		(y) = tmp / (nx); \
		(z) = idx / ( (nx)*(ny) ); \
	    }
