#line 1 "scopy.f"
/* scopy.f -- translated by f2c (version 20100827).
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

#line 1 "scopy.f"
/* > \brief \b SCOPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SCOPY(N,SX,INCX,SY,INCY) */

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
/* >    SCOPY copies a vector, x, to a vector, y. */
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
/* Subroutine */ int scopy_(integer *n, doublereal *sx, integer *incx, 
	doublereal *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


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
#line 74 "scopy.f"
    /* Parameter adjustments */
#line 74 "scopy.f"
    --sy;
#line 74 "scopy.f"
    --sx;
#line 74 "scopy.f"

#line 74 "scopy.f"
    /* Function Body */
#line 74 "scopy.f"
    if (*n <= 0) {
#line 74 "scopy.f"
	return 0;
#line 74 "scopy.f"
    }
#line 75 "scopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 82 "scopy.f"
	m = *n % 7;
#line 83 "scopy.f"
	if (m != 0) {
#line 84 "scopy.f"
	    i__1 = m;
#line 84 "scopy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 85 "scopy.f"
		sy[i__] = sx[i__];
#line 86 "scopy.f"
	    }
#line 87 "scopy.f"
	    if (*n < 7) {
#line 87 "scopy.f"
		return 0;
#line 87 "scopy.f"
	    }
#line 88 "scopy.f"
	}
#line 89 "scopy.f"
	mp1 = m + 1;
#line 90 "scopy.f"
	i__1 = *n;
#line 90 "scopy.f"
	for (i__ = mp1; i__ <= i__1; i__ += 7) {
#line 91 "scopy.f"
	    sy[i__] = sx[i__];
#line 92 "scopy.f"
	    sy[i__ + 1] = sx[i__ + 1];
#line 93 "scopy.f"
	    sy[i__ + 2] = sx[i__ + 2];
#line 94 "scopy.f"
	    sy[i__ + 3] = sx[i__ + 3];
#line 95 "scopy.f"
	    sy[i__ + 4] = sx[i__ + 4];
#line 96 "scopy.f"
	    sy[i__ + 5] = sx[i__ + 5];
#line 97 "scopy.f"
	    sy[i__ + 6] = sx[i__ + 6];
#line 98 "scopy.f"
	}
#line 99 "scopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 104 "scopy.f"
	ix = 1;
#line 105 "scopy.f"
	iy = 1;
#line 106 "scopy.f"
	if (*incx < 0) {
#line 106 "scopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 106 "scopy.f"
	}
#line 107 "scopy.f"
	if (*incy < 0) {
#line 107 "scopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 107 "scopy.f"
	}
#line 108 "scopy.f"
	i__1 = *n;
#line 108 "scopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 109 "scopy.f"
	    sy[iy] = sx[ix];
#line 110 "scopy.f"
	    ix += *incx;
#line 111 "scopy.f"
	    iy += *incy;
#line 112 "scopy.f"
	}
#line 113 "scopy.f"
    }
#line 114 "scopy.f"
    return 0;
} /* scopy_ */

