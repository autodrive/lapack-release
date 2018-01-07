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
/* >    uses unrolled loops for increments equal to 1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] DX */
/* > \verbatim */
/* >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of DX */
/* > \endverbatim */
/* > */
/* > \param[out] DY */
/* > \verbatim */
/* >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of DY */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

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


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 105 "dcopy.f"
    /* Parameter adjustments */
#line 105 "dcopy.f"
    --dy;
#line 105 "dcopy.f"
    --dx;
#line 105 "dcopy.f"

#line 105 "dcopy.f"
    /* Function Body */
#line 105 "dcopy.f"
    if (*n <= 0) {
#line 105 "dcopy.f"
	return 0;
#line 105 "dcopy.f"
    }
#line 106 "dcopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 113 "dcopy.f"
	m = *n % 7;
#line 114 "dcopy.f"
	if (m != 0) {
#line 115 "dcopy.f"
	    i__1 = m;
#line 115 "dcopy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 116 "dcopy.f"
		dy[i__] = dx[i__];
#line 117 "dcopy.f"
	    }
#line 118 "dcopy.f"
	    if (*n < 7) {
#line 118 "dcopy.f"
		return 0;
#line 118 "dcopy.f"
	    }
#line 119 "dcopy.f"
	}
#line 120 "dcopy.f"
	mp1 = m + 1;
#line 121 "dcopy.f"
	i__1 = *n;
#line 121 "dcopy.f"
	for (i__ = mp1; i__ <= i__1; i__ += 7) {
#line 122 "dcopy.f"
	    dy[i__] = dx[i__];
#line 123 "dcopy.f"
	    dy[i__ + 1] = dx[i__ + 1];
#line 124 "dcopy.f"
	    dy[i__ + 2] = dx[i__ + 2];
#line 125 "dcopy.f"
	    dy[i__ + 3] = dx[i__ + 3];
#line 126 "dcopy.f"
	    dy[i__ + 4] = dx[i__ + 4];
#line 127 "dcopy.f"
	    dy[i__ + 5] = dx[i__ + 5];
#line 128 "dcopy.f"
	    dy[i__ + 6] = dx[i__ + 6];
#line 129 "dcopy.f"
	}
#line 130 "dcopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 135 "dcopy.f"
	ix = 1;
#line 136 "dcopy.f"
	iy = 1;
#line 137 "dcopy.f"
	if (*incx < 0) {
#line 137 "dcopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 137 "dcopy.f"
	}
#line 138 "dcopy.f"
	if (*incy < 0) {
#line 138 "dcopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 138 "dcopy.f"
	}
#line 139 "dcopy.f"
	i__1 = *n;
#line 139 "dcopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 140 "dcopy.f"
	    dy[iy] = dx[ix];
#line 141 "dcopy.f"
	    ix += *incx;
#line 142 "dcopy.f"
	    iy += *incy;
#line 143 "dcopy.f"
	}
#line 144 "dcopy.f"
    }
#line 145 "dcopy.f"
    return 0;
} /* dcopy_ */

