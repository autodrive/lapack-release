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
/* > \param[in] DY */
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
doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal dtemp;


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
#line 106 "ddot.f"
    /* Parameter adjustments */
#line 106 "ddot.f"
    --dy;
#line 106 "ddot.f"
    --dx;
#line 106 "ddot.f"

#line 106 "ddot.f"
    /* Function Body */
#line 106 "ddot.f"
    ret_val = 0.;
#line 107 "ddot.f"
    dtemp = 0.;
#line 108 "ddot.f"
    if (*n <= 0) {
#line 108 "ddot.f"
	return ret_val;
#line 108 "ddot.f"
    }
#line 109 "ddot.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 116 "ddot.f"
	m = *n % 5;
#line 117 "ddot.f"
	if (m != 0) {
#line 118 "ddot.f"
	    i__1 = m;
#line 118 "ddot.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 119 "ddot.f"
		dtemp += dx[i__] * dy[i__];
#line 120 "ddot.f"
	    }
#line 121 "ddot.f"
	    if (*n < 5) {
#line 122 "ddot.f"
		ret_val = dtemp;
#line 123 "ddot.f"
		return ret_val;
#line 124 "ddot.f"
	    }
#line 125 "ddot.f"
	}
#line 126 "ddot.f"
	mp1 = m + 1;
#line 127 "ddot.f"
	i__1 = *n;
#line 127 "ddot.f"
	for (i__ = mp1; i__ <= i__1; i__ += 5) {
#line 128 "ddot.f"
	    dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + 
		    dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + 
		    dx[i__ + 4] * dy[i__ + 4];
#line 130 "ddot.f"
	}
#line 131 "ddot.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 136 "ddot.f"
	ix = 1;
#line 137 "ddot.f"
	iy = 1;
#line 138 "ddot.f"
	if (*incx < 0) {
#line 138 "ddot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 138 "ddot.f"
	}
#line 139 "ddot.f"
	if (*incy < 0) {
#line 139 "ddot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 139 "ddot.f"
	}
#line 140 "ddot.f"
	i__1 = *n;
#line 140 "ddot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 141 "ddot.f"
	    dtemp += dx[ix] * dy[iy];
#line 142 "ddot.f"
	    ix += *incx;
#line 143 "ddot.f"
	    iy += *incy;
#line 144 "ddot.f"
	}
#line 145 "ddot.f"
    }
#line 146 "ddot.f"
    ret_val = dtemp;
#line 147 "ddot.f"
    return ret_val;
} /* ddot_ */

