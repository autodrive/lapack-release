#line 1 "sswap.f"
/* sswap.f -- translated by f2c (version 20100827).
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

#line 1 "sswap.f"
/* > \brief \b SSWAP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSWAP(N,SX,INCX,SY,INCY) */

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
/* >    interchanges two vectors. */
/* >    uses unrolled loops for increments equal to 1. */
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
/* Subroutine */ int sswap_(integer *n, doublereal *sx, integer *incx, 
	doublereal *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;

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
#line 75 "sswap.f"
    /* Parameter adjustments */
#line 75 "sswap.f"
    --sy;
#line 75 "sswap.f"
    --sx;
#line 75 "sswap.f"

#line 75 "sswap.f"
    /* Function Body */
#line 75 "sswap.f"
    if (*n <= 0) {
#line 75 "sswap.f"
	return 0;
#line 75 "sswap.f"
    }
#line 76 "sswap.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */


/*       clean-up loop */

#line 83 "sswap.f"
	m = *n % 3;
#line 84 "sswap.f"
	if (m != 0) {
#line 85 "sswap.f"
	    i__1 = m;
#line 85 "sswap.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 86 "sswap.f"
		stemp = sx[i__];
#line 87 "sswap.f"
		sx[i__] = sy[i__];
#line 88 "sswap.f"
		sy[i__] = stemp;
#line 89 "sswap.f"
	    }
#line 90 "sswap.f"
	    if (*n < 3) {
#line 90 "sswap.f"
		return 0;
#line 90 "sswap.f"
	    }
#line 91 "sswap.f"
	}
#line 92 "sswap.f"
	mp1 = m + 1;
#line 93 "sswap.f"
	i__1 = *n;
#line 93 "sswap.f"
	for (i__ = mp1; i__ <= i__1; i__ += 3) {
#line 94 "sswap.f"
	    stemp = sx[i__];
#line 95 "sswap.f"
	    sx[i__] = sy[i__];
#line 96 "sswap.f"
	    sy[i__] = stemp;
#line 97 "sswap.f"
	    stemp = sx[i__ + 1];
#line 98 "sswap.f"
	    sx[i__ + 1] = sy[i__ + 1];
#line 99 "sswap.f"
	    sy[i__ + 1] = stemp;
#line 100 "sswap.f"
	    stemp = sx[i__ + 2];
#line 101 "sswap.f"
	    sx[i__ + 2] = sy[i__ + 2];
#line 102 "sswap.f"
	    sy[i__ + 2] = stemp;
#line 103 "sswap.f"
	}
#line 104 "sswap.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 109 "sswap.f"
	ix = 1;
#line 110 "sswap.f"
	iy = 1;
#line 111 "sswap.f"
	if (*incx < 0) {
#line 111 "sswap.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 111 "sswap.f"
	}
#line 112 "sswap.f"
	if (*incy < 0) {
#line 112 "sswap.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 112 "sswap.f"
	}
#line 113 "sswap.f"
	i__1 = *n;
#line 113 "sswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 114 "sswap.f"
	    stemp = sx[ix];
#line 115 "sswap.f"
	    sx[ix] = sy[iy];
#line 116 "sswap.f"
	    sy[iy] = stemp;
#line 117 "sswap.f"
	    ix += *incx;
#line 118 "sswap.f"
	    iy += *incy;
#line 119 "sswap.f"
	}
#line 120 "sswap.f"
    }
#line 121 "sswap.f"
    return 0;
} /* sswap_ */

