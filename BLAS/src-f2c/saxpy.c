#line 1 "saxpy.f"
/* saxpy.f -- translated by f2c (version 20100827).
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

#line 1 "saxpy.f"
/* > \brief \b SAXPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL SA */
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
/* >    SAXPY constant times a vector plus a vector. */
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
/* > \param[in] SA */
/* > \verbatim */
/* >          SA is REAL */
/* >           On entry, SA specifies the scalar alpha. */
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
/* Subroutine */ int saxpy_(integer *n, doublereal *sa, doublereal *sx, 
	integer *incx, doublereal *sy, integer *incy)
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
#line 113 "saxpy.f"
    /* Parameter adjustments */
#line 113 "saxpy.f"
    --sy;
#line 113 "saxpy.f"
    --sx;
#line 113 "saxpy.f"

#line 113 "saxpy.f"
    /* Function Body */
#line 113 "saxpy.f"
    if (*n <= 0) {
#line 113 "saxpy.f"
	return 0;
#line 113 "saxpy.f"
    }
#line 114 "saxpy.f"
    if (*sa == 0.) {
#line 114 "saxpy.f"
	return 0;
#line 114 "saxpy.f"
    }
#line 115 "saxpy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 122 "saxpy.f"
	m = *n % 4;
#line 123 "saxpy.f"
	if (m != 0) {
#line 124 "saxpy.f"
	    i__1 = m;
#line 124 "saxpy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 125 "saxpy.f"
		sy[i__] += *sa * sx[i__];
#line 126 "saxpy.f"
	    }
#line 127 "saxpy.f"
	}
#line 128 "saxpy.f"
	if (*n < 4) {
#line 128 "saxpy.f"
	    return 0;
#line 128 "saxpy.f"
	}
#line 129 "saxpy.f"
	mp1 = m + 1;
#line 130 "saxpy.f"
	i__1 = *n;
#line 130 "saxpy.f"
	for (i__ = mp1; i__ <= i__1; i__ += 4) {
#line 131 "saxpy.f"
	    sy[i__] += *sa * sx[i__];
#line 132 "saxpy.f"
	    sy[i__ + 1] += *sa * sx[i__ + 1];
#line 133 "saxpy.f"
	    sy[i__ + 2] += *sa * sx[i__ + 2];
#line 134 "saxpy.f"
	    sy[i__ + 3] += *sa * sx[i__ + 3];
#line 135 "saxpy.f"
	}
#line 136 "saxpy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 141 "saxpy.f"
	ix = 1;
#line 142 "saxpy.f"
	iy = 1;
#line 143 "saxpy.f"
	if (*incx < 0) {
#line 143 "saxpy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 143 "saxpy.f"
	}
#line 144 "saxpy.f"
	if (*incy < 0) {
#line 144 "saxpy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 144 "saxpy.f"
	}
#line 145 "saxpy.f"
	i__1 = *n;
#line 145 "saxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 146 "saxpy.f"
	    sy[iy] += *sa * sx[ix];
#line 147 "saxpy.f"
	    ix += *incx;
#line 148 "saxpy.f"
	    iy += *incy;
#line 149 "saxpy.f"
	}
#line 150 "saxpy.f"
    }
#line 151 "saxpy.f"
    return 0;
} /* saxpy_ */

