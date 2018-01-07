#line 1 "dswap.f"
/* dswap.f -- translated by f2c (version 20100827).
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

#line 1 "dswap.f"
/* > \brief \b DSWAP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY) */

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
/* >    interchanges two vectors. */
/* >    uses unrolled loops for increments equal one. */
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
/* Subroutine */ int dswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

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
#line 75 "dswap.f"
    /* Parameter adjustments */
#line 75 "dswap.f"
    --dy;
#line 75 "dswap.f"
    --dx;
#line 75 "dswap.f"

#line 75 "dswap.f"
    /* Function Body */
#line 75 "dswap.f"
    if (*n <= 0) {
#line 75 "dswap.f"
	return 0;
#line 75 "dswap.f"
    }
#line 76 "dswap.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */


/*       clean-up loop */

#line 83 "dswap.f"
	m = *n % 3;
#line 84 "dswap.f"
	if (m != 0) {
#line 85 "dswap.f"
	    i__1 = m;
#line 85 "dswap.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 86 "dswap.f"
		dtemp = dx[i__];
#line 87 "dswap.f"
		dx[i__] = dy[i__];
#line 88 "dswap.f"
		dy[i__] = dtemp;
#line 89 "dswap.f"
	    }
#line 90 "dswap.f"
	    if (*n < 3) {
#line 90 "dswap.f"
		return 0;
#line 90 "dswap.f"
	    }
#line 91 "dswap.f"
	}
#line 92 "dswap.f"
	mp1 = m + 1;
#line 93 "dswap.f"
	i__1 = *n;
#line 93 "dswap.f"
	for (i__ = mp1; i__ <= i__1; i__ += 3) {
#line 94 "dswap.f"
	    dtemp = dx[i__];
#line 95 "dswap.f"
	    dx[i__] = dy[i__];
#line 96 "dswap.f"
	    dy[i__] = dtemp;
#line 97 "dswap.f"
	    dtemp = dx[i__ + 1];
#line 98 "dswap.f"
	    dx[i__ + 1] = dy[i__ + 1];
#line 99 "dswap.f"
	    dy[i__ + 1] = dtemp;
#line 100 "dswap.f"
	    dtemp = dx[i__ + 2];
#line 101 "dswap.f"
	    dx[i__ + 2] = dy[i__ + 2];
#line 102 "dswap.f"
	    dy[i__ + 2] = dtemp;
#line 103 "dswap.f"
	}
#line 104 "dswap.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 109 "dswap.f"
	ix = 1;
#line 110 "dswap.f"
	iy = 1;
#line 111 "dswap.f"
	if (*incx < 0) {
#line 111 "dswap.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 111 "dswap.f"
	}
#line 112 "dswap.f"
	if (*incy < 0) {
#line 112 "dswap.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 112 "dswap.f"
	}
#line 113 "dswap.f"
	i__1 = *n;
#line 113 "dswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 114 "dswap.f"
	    dtemp = dx[ix];
#line 115 "dswap.f"
	    dx[ix] = dy[iy];
#line 116 "dswap.f"
	    dy[iy] = dtemp;
#line 117 "dswap.f"
	    ix += *incx;
#line 118 "dswap.f"
	    iy += *incy;
#line 119 "dswap.f"
	}
#line 120 "dswap.f"
    }
#line 121 "dswap.f"
    return 0;
} /* dswap_ */

