#line 1 "daxpy.f"
/* daxpy.f -- translated by f2c (version 20100827).
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

#line 1 "daxpy.f"
/* > \brief \b DAXPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION DA */
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
/* >    DAXPY constant times a vector plus a vector. */
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
/* > \param[in] DA */
/* > \verbatim */
/* >          DA is DOUBLE PRECISION */
/* >           On entry, DA specifies the scalar alpha. */
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
/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
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
#line 113 "daxpy.f"
    /* Parameter adjustments */
#line 113 "daxpy.f"
    --dy;
#line 113 "daxpy.f"
    --dx;
#line 113 "daxpy.f"

#line 113 "daxpy.f"
    /* Function Body */
#line 113 "daxpy.f"
    if (*n <= 0) {
#line 113 "daxpy.f"
	return 0;
#line 113 "daxpy.f"
    }
#line 114 "daxpy.f"
    if (*da == 0.) {
#line 114 "daxpy.f"
	return 0;
#line 114 "daxpy.f"
    }
#line 115 "daxpy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 122 "daxpy.f"
	m = *n % 4;
#line 123 "daxpy.f"
	if (m != 0) {
#line 124 "daxpy.f"
	    i__1 = m;
#line 124 "daxpy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 125 "daxpy.f"
		dy[i__] += *da * dx[i__];
#line 126 "daxpy.f"
	    }
#line 127 "daxpy.f"
	}
#line 128 "daxpy.f"
	if (*n < 4) {
#line 128 "daxpy.f"
	    return 0;
#line 128 "daxpy.f"
	}
#line 129 "daxpy.f"
	mp1 = m + 1;
#line 130 "daxpy.f"
	i__1 = *n;
#line 130 "daxpy.f"
	for (i__ = mp1; i__ <= i__1; i__ += 4) {
#line 131 "daxpy.f"
	    dy[i__] += *da * dx[i__];
#line 132 "daxpy.f"
	    dy[i__ + 1] += *da * dx[i__ + 1];
#line 133 "daxpy.f"
	    dy[i__ + 2] += *da * dx[i__ + 2];
#line 134 "daxpy.f"
	    dy[i__ + 3] += *da * dx[i__ + 3];
#line 135 "daxpy.f"
	}
#line 136 "daxpy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 141 "daxpy.f"
	ix = 1;
#line 142 "daxpy.f"
	iy = 1;
#line 143 "daxpy.f"
	if (*incx < 0) {
#line 143 "daxpy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 143 "daxpy.f"
	}
#line 144 "daxpy.f"
	if (*incy < 0) {
#line 144 "daxpy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 144 "daxpy.f"
	}
#line 145 "daxpy.f"
	i__1 = *n;
#line 145 "daxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 146 "daxpy.f"
	    dy[iy] += *da * dx[ix];
#line 147 "daxpy.f"
	    ix += *incx;
#line 148 "daxpy.f"
	    iy += *incy;
#line 149 "daxpy.f"
	}
#line 150 "daxpy.f"
    }
#line 151 "daxpy.f"
    return 0;
} /* daxpy_ */

