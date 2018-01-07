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
/* >    SSWAP interchanges two vectors. */
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
/* > \param[in,out] SX */
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
/* > \param[in,out] SY */
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
/* Subroutine */ int sswap_(integer *n, doublereal *sx, integer *incx, 
	doublereal *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal stemp;


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
#line 106 "sswap.f"
    /* Parameter adjustments */
#line 106 "sswap.f"
    --sy;
#line 106 "sswap.f"
    --sx;
#line 106 "sswap.f"

#line 106 "sswap.f"
    /* Function Body */
#line 106 "sswap.f"
    if (*n <= 0) {
#line 106 "sswap.f"
	return 0;
#line 106 "sswap.f"
    }
#line 107 "sswap.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */


/*       clean-up loop */

#line 114 "sswap.f"
	m = *n % 3;
#line 115 "sswap.f"
	if (m != 0) {
#line 116 "sswap.f"
	    i__1 = m;
#line 116 "sswap.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 117 "sswap.f"
		stemp = sx[i__];
#line 118 "sswap.f"
		sx[i__] = sy[i__];
#line 119 "sswap.f"
		sy[i__] = stemp;
#line 120 "sswap.f"
	    }
#line 121 "sswap.f"
	    if (*n < 3) {
#line 121 "sswap.f"
		return 0;
#line 121 "sswap.f"
	    }
#line 122 "sswap.f"
	}
#line 123 "sswap.f"
	mp1 = m + 1;
#line 124 "sswap.f"
	i__1 = *n;
#line 124 "sswap.f"
	for (i__ = mp1; i__ <= i__1; i__ += 3) {
#line 125 "sswap.f"
	    stemp = sx[i__];
#line 126 "sswap.f"
	    sx[i__] = sy[i__];
#line 127 "sswap.f"
	    sy[i__] = stemp;
#line 128 "sswap.f"
	    stemp = sx[i__ + 1];
#line 129 "sswap.f"
	    sx[i__ + 1] = sy[i__ + 1];
#line 130 "sswap.f"
	    sy[i__ + 1] = stemp;
#line 131 "sswap.f"
	    stemp = sx[i__ + 2];
#line 132 "sswap.f"
	    sx[i__ + 2] = sy[i__ + 2];
#line 133 "sswap.f"
	    sy[i__ + 2] = stemp;
#line 134 "sswap.f"
	}
#line 135 "sswap.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 140 "sswap.f"
	ix = 1;
#line 141 "sswap.f"
	iy = 1;
#line 142 "sswap.f"
	if (*incx < 0) {
#line 142 "sswap.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 142 "sswap.f"
	}
#line 143 "sswap.f"
	if (*incy < 0) {
#line 143 "sswap.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 143 "sswap.f"
	}
#line 144 "sswap.f"
	i__1 = *n;
#line 144 "sswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 145 "sswap.f"
	    stemp = sx[ix];
#line 146 "sswap.f"
	    sx[ix] = sy[iy];
#line 147 "sswap.f"
	    sy[iy] = stemp;
#line 148 "sswap.f"
	    ix += *incx;
#line 149 "sswap.f"
	    iy += *incy;
#line 150 "sswap.f"
	}
#line 151 "sswap.f"
    }
#line 152 "sswap.f"
    return 0;
} /* sswap_ */

