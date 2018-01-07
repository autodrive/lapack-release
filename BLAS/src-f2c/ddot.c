#line 1 "ddot.f"
/* ddot.f -- translated by f2c (version 20100827).
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

#line 1 "ddot.f"
/* > \brief \b DDOT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DX(*),DY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DDOT forms the dot product of two vectors. */
/* >    uses unrolled loops for increments equal to one. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup double_blas_level1 */

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
doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal dtemp;


/*  -- Reference BLAS level1 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 75 "ddot.f"
    /* Parameter adjustments */
#line 75 "ddot.f"
    --dy;
#line 75 "ddot.f"
    --dx;
#line 75 "ddot.f"

#line 75 "ddot.f"
    /* Function Body */
#line 75 "ddot.f"
    ret_val = 0.;
#line 76 "ddot.f"
    dtemp = 0.;
#line 77 "ddot.f"
    if (*n <= 0) {
#line 77 "ddot.f"
	return ret_val;
#line 77 "ddot.f"
    }
#line 78 "ddot.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 85 "ddot.f"
	m = *n % 5;
#line 86 "ddot.f"
	if (m != 0) {
#line 87 "ddot.f"
	    i__1 = m;
#line 87 "ddot.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 88 "ddot.f"
		dtemp += dx[i__] * dy[i__];
#line 89 "ddot.f"
	    }
#line 90 "ddot.f"
	    if (*n < 5) {
#line 91 "ddot.f"
		ret_val = dtemp;
#line 92 "ddot.f"
		return ret_val;
#line 93 "ddot.f"
	    }
#line 94 "ddot.f"
	}
#line 95 "ddot.f"
	mp1 = m + 1;
#line 96 "ddot.f"
	i__1 = *n;
#line 96 "ddot.f"
	for (i__ = mp1; i__ <= i__1; i__ += 5) {
#line 97 "ddot.f"
	    dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + 
		    dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + 
		    dx[i__ + 4] * dy[i__ + 4];
#line 99 "ddot.f"
	}
#line 100 "ddot.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 105 "ddot.f"
	ix = 1;
#line 106 "ddot.f"
	iy = 1;
#line 107 "ddot.f"
	if (*incx < 0) {
#line 107 "ddot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 107 "ddot.f"
	}
#line 108 "ddot.f"
	if (*incy < 0) {
#line 108 "ddot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 108 "ddot.f"
	}
#line 109 "ddot.f"
	i__1 = *n;
#line 109 "ddot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 110 "ddot.f"
	    dtemp += dx[ix] * dy[iy];
#line 111 "ddot.f"
	    ix += *incx;
#line 112 "ddot.f"
	    iy += *incy;
#line 113 "ddot.f"
	}
#line 114 "ddot.f"
    }
#line 115 "ddot.f"
    ret_val = dtemp;
#line 116 "ddot.f"
    return ret_val;
} /* ddot_ */

