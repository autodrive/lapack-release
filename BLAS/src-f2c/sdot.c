#line 1 "sdot.f"
/* sdot.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "sdot.f"
/* > \brief \b SDOT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SDOT(N,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*),SY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SDOT forms the dot product of two vectors. */
/* >    uses unrolled loops for increments equal to one. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup single_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal sdot_(integer *n, doublereal *sx, integer *incx, doublereal *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal stemp;


/*  -- Reference BLAS level1 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 75 "sdot.f"
    /* Parameter adjustments */
#line 75 "sdot.f"
    --sy;
#line 75 "sdot.f"
    --sx;
#line 75 "sdot.f"

#line 75 "sdot.f"
    /* Function Body */
#line 75 "sdot.f"
    stemp = 0.;
#line 76 "sdot.f"
    ret_val = 0.;
#line 77 "sdot.f"
    if (*n <= 0) {
#line 77 "sdot.f"
	return ret_val;
#line 77 "sdot.f"
    }
#line 78 "sdot.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 85 "sdot.f"
	m = *n % 5;
#line 86 "sdot.f"
	if (m != 0) {
#line 87 "sdot.f"
	    i__1 = m;
#line 87 "sdot.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 88 "sdot.f"
		stemp += sx[i__] * sy[i__];
#line 89 "sdot.f"
	    }
#line 90 "sdot.f"
	    if (*n < 5) {
#line 91 "sdot.f"
		ret_val = stemp;
#line 92 "sdot.f"
		return ret_val;
#line 93 "sdot.f"
	    }
#line 94 "sdot.f"
	}
#line 95 "sdot.f"
	mp1 = m + 1;
#line 96 "sdot.f"
	i__1 = *n;
#line 96 "sdot.f"
	for (i__ = mp1; i__ <= i__1; i__ += 5) {
#line 97 "sdot.f"
	    stemp = stemp + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] + 
		    sx[i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] + 
		    sx[i__ + 4] * sy[i__ + 4];
#line 99 "sdot.f"
	}
#line 100 "sdot.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 105 "sdot.f"
	ix = 1;
#line 106 "sdot.f"
	iy = 1;
#line 107 "sdot.f"
	if (*incx < 0) {
#line 107 "sdot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 107 "sdot.f"
	}
#line 108 "sdot.f"
	if (*incy < 0) {
#line 108 "sdot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 108 "sdot.f"
	}
#line 109 "sdot.f"
	i__1 = *n;
#line 109 "sdot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 110 "sdot.f"
	    stemp += sx[ix] * sy[iy];
#line 111 "sdot.f"
	    ix += *incx;
#line 112 "sdot.f"
	    iy += *incy;
#line 113 "sdot.f"
	}
#line 114 "sdot.f"
    }
#line 115 "sdot.f"
    ret_val = stemp;
#line 116 "sdot.f"
    return ret_val;
} /* sdot_ */

