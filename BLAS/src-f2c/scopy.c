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

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] SX */
/* > \verbatim */
/* >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of SX */
/* > \endverbatim */
/* > */
/* > \param[out] SY */
/* > \verbatim */
/* >          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of SY */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

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
#line 105 "scopy.f"
    /* Parameter adjustments */
#line 105 "scopy.f"
    --sy;
#line 105 "scopy.f"
    --sx;
#line 105 "scopy.f"

#line 105 "scopy.f"
    /* Function Body */
#line 105 "scopy.f"
    if (*n <= 0) {
#line 105 "scopy.f"
	return 0;
#line 105 "scopy.f"
    }
#line 106 "scopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 113 "scopy.f"
	m = *n % 7;
#line 114 "scopy.f"
	if (m != 0) {
#line 115 "scopy.f"
	    i__1 = m;
#line 115 "scopy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 116 "scopy.f"
		sy[i__] = sx[i__];
#line 117 "scopy.f"
	    }
#line 118 "scopy.f"
	    if (*n < 7) {
#line 118 "scopy.f"
		return 0;
#line 118 "scopy.f"
	    }
#line 119 "scopy.f"
	}
#line 120 "scopy.f"
	mp1 = m + 1;
#line 121 "scopy.f"
	i__1 = *n;
#line 121 "scopy.f"
	for (i__ = mp1; i__ <= i__1; i__ += 7) {
#line 122 "scopy.f"
	    sy[i__] = sx[i__];
#line 123 "scopy.f"
	    sy[i__ + 1] = sx[i__ + 1];
#line 124 "scopy.f"
	    sy[i__ + 2] = sx[i__ + 2];
#line 125 "scopy.f"
	    sy[i__ + 3] = sx[i__ + 3];
#line 126 "scopy.f"
	    sy[i__ + 4] = sx[i__ + 4];
#line 127 "scopy.f"
	    sy[i__ + 5] = sx[i__ + 5];
#line 128 "scopy.f"
	    sy[i__ + 6] = sx[i__ + 6];
#line 129 "scopy.f"
	}
#line 130 "scopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 135 "scopy.f"
	ix = 1;
#line 136 "scopy.f"
	iy = 1;
#line 137 "scopy.f"
	if (*incx < 0) {
#line 137 "scopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 137 "scopy.f"
	}
#line 138 "scopy.f"
	if (*incy < 0) {
#line 138 "scopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 138 "scopy.f"
	}
#line 139 "scopy.f"
	i__1 = *n;
#line 139 "scopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 140 "scopy.f"
	    sy[iy] = sx[ix];
#line 141 "scopy.f"
	    ix += *incx;
#line 142 "scopy.f"
	    iy += *incy;
#line 143 "scopy.f"
	}
#line 144 "scopy.f"
    }
#line 145 "scopy.f"
    return 0;
} /* scopy_ */

