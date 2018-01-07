#line 1 "dcopy.f"
/* dcopy.f -- translated by f2c (version 20100827).
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

#line 1 "dcopy.f"
/* > \brief \b DCOPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY) */

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
/* >    DCOPY copies a vector, x, to a vector, y. */
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
/* Subroutine */ int dcopy_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


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
#line 74 "dcopy.f"
    /* Parameter adjustments */
#line 74 "dcopy.f"
    --dy;
#line 74 "dcopy.f"
    --dx;
#line 74 "dcopy.f"

#line 74 "dcopy.f"
    /* Function Body */
#line 74 "dcopy.f"
    if (*n <= 0) {
#line 74 "dcopy.f"
	return 0;
#line 74 "dcopy.f"
    }
#line 75 "dcopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 82 "dcopy.f"
	m = *n % 7;
#line 83 "dcopy.f"
	if (m != 0) {
#line 84 "dcopy.f"
	    i__1 = m;
#line 84 "dcopy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 85 "dcopy.f"
		dy[i__] = dx[i__];
#line 86 "dcopy.f"
	    }
#line 87 "dcopy.f"
	    if (*n < 7) {
#line 87 "dcopy.f"
		return 0;
#line 87 "dcopy.f"
	    }
#line 88 "dcopy.f"
	}
#line 89 "dcopy.f"
	mp1 = m + 1;
#line 90 "dcopy.f"
	i__1 = *n;
#line 90 "dcopy.f"
	for (i__ = mp1; i__ <= i__1; i__ += 7) {
#line 91 "dcopy.f"
	    dy[i__] = dx[i__];
#line 92 "dcopy.f"
	    dy[i__ + 1] = dx[i__ + 1];
#line 93 "dcopy.f"
	    dy[i__ + 2] = dx[i__ + 2];
#line 94 "dcopy.f"
	    dy[i__ + 3] = dx[i__ + 3];
#line 95 "dcopy.f"
	    dy[i__ + 4] = dx[i__ + 4];
#line 96 "dcopy.f"
	    dy[i__ + 5] = dx[i__ + 5];
#line 97 "dcopy.f"
	    dy[i__ + 6] = dx[i__ + 6];
#line 98 "dcopy.f"
	}
#line 99 "dcopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 104 "dcopy.f"
	ix = 1;
#line 105 "dcopy.f"
	iy = 1;
#line 106 "dcopy.f"
	if (*incx < 0) {
#line 106 "dcopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 106 "dcopy.f"
	}
#line 107 "dcopy.f"
	if (*incy < 0) {
#line 107 "dcopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 107 "dcopy.f"
	}
#line 108 "dcopy.f"
	i__1 = *n;
#line 108 "dcopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 109 "dcopy.f"
	    dy[iy] = dx[ix];
#line 110 "dcopy.f"
	    ix += *incx;
#line 111 "dcopy.f"
	    iy += *incy;
#line 112 "dcopy.f"
	}
#line 113 "dcopy.f"
    }
#line 114 "dcopy.f"
    return 0;
} /* dcopy_ */

