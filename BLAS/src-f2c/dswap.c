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
/* >    DSWAP interchanges two vectors. */
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
/* > \param[in,out] DX */
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
/* > \param[in,out] DY */
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
/* Subroutine */ int dswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

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
#line 106 "dswap.f"
    /* Parameter adjustments */
#line 106 "dswap.f"
    --dy;
#line 106 "dswap.f"
    --dx;
#line 106 "dswap.f"

#line 106 "dswap.f"
    /* Function Body */
#line 106 "dswap.f"
    if (*n <= 0) {
#line 106 "dswap.f"
	return 0;
#line 106 "dswap.f"
    }
#line 107 "dswap.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */


/*       clean-up loop */

#line 114 "dswap.f"
	m = *n % 3;
#line 115 "dswap.f"
	if (m != 0) {
#line 116 "dswap.f"
	    i__1 = m;
#line 116 "dswap.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 117 "dswap.f"
		dtemp = dx[i__];
#line 118 "dswap.f"
		dx[i__] = dy[i__];
#line 119 "dswap.f"
		dy[i__] = dtemp;
#line 120 "dswap.f"
	    }
#line 121 "dswap.f"
	    if (*n < 3) {
#line 121 "dswap.f"
		return 0;
#line 121 "dswap.f"
	    }
#line 122 "dswap.f"
	}
#line 123 "dswap.f"
	mp1 = m + 1;
#line 124 "dswap.f"
	i__1 = *n;
#line 124 "dswap.f"
	for (i__ = mp1; i__ <= i__1; i__ += 3) {
#line 125 "dswap.f"
	    dtemp = dx[i__];
#line 126 "dswap.f"
	    dx[i__] = dy[i__];
#line 127 "dswap.f"
	    dy[i__] = dtemp;
#line 128 "dswap.f"
	    dtemp = dx[i__ + 1];
#line 129 "dswap.f"
	    dx[i__ + 1] = dy[i__ + 1];
#line 130 "dswap.f"
	    dy[i__ + 1] = dtemp;
#line 131 "dswap.f"
	    dtemp = dx[i__ + 2];
#line 132 "dswap.f"
	    dx[i__ + 2] = dy[i__ + 2];
#line 133 "dswap.f"
	    dy[i__ + 2] = dtemp;
#line 134 "dswap.f"
	}
#line 135 "dswap.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 140 "dswap.f"
	ix = 1;
#line 141 "dswap.f"
	iy = 1;
#line 142 "dswap.f"
	if (*incx < 0) {
#line 142 "dswap.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 142 "dswap.f"
	}
#line 143 "dswap.f"
	if (*incy < 0) {
#line 143 "dswap.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 143 "dswap.f"
	}
#line 144 "dswap.f"
	i__1 = *n;
#line 144 "dswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 145 "dswap.f"
	    dtemp = dx[ix];
#line 146 "dswap.f"
	    dx[ix] = dy[iy];
#line 147 "dswap.f"
	    dy[iy] = dtemp;
#line 148 "dswap.f"
	    ix += *incx;
#line 149 "dswap.f"
	    iy += *incy;
#line 150 "dswap.f"
	}
#line 151 "dswap.f"
    }
#line 152 "dswap.f"
    return 0;
} /* dswap_ */

