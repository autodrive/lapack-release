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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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


/*  -- Reference BLAS level1 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 76 "saxpy.f"
    /* Parameter adjustments */
#line 76 "saxpy.f"
    --sy;
#line 76 "saxpy.f"
    --sx;
#line 76 "saxpy.f"

#line 76 "saxpy.f"
    /* Function Body */
#line 76 "saxpy.f"
    if (*n <= 0) {
#line 76 "saxpy.f"
	return 0;
#line 76 "saxpy.f"
    }
#line 77 "saxpy.f"
    if (*sa == 0.) {
#line 77 "saxpy.f"
	return 0;
#line 77 "saxpy.f"
    }
#line 78 "saxpy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

#line 85 "saxpy.f"
	m = *n % 4;
#line 86 "saxpy.f"
	if (m != 0) {
#line 87 "saxpy.f"
	    i__1 = m;
#line 87 "saxpy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 88 "saxpy.f"
		sy[i__] += *sa * sx[i__];
#line 89 "saxpy.f"
	    }
#line 90 "saxpy.f"
	}
#line 91 "saxpy.f"
	if (*n < 4) {
#line 91 "saxpy.f"
	    return 0;
#line 91 "saxpy.f"
	}
#line 92 "saxpy.f"
	mp1 = m + 1;
#line 93 "saxpy.f"
	i__1 = *n;
#line 93 "saxpy.f"
	for (i__ = mp1; i__ <= i__1; i__ += 4) {
#line 94 "saxpy.f"
	    sy[i__] += *sa * sx[i__];
#line 95 "saxpy.f"
	    sy[i__ + 1] += *sa * sx[i__ + 1];
#line 96 "saxpy.f"
	    sy[i__ + 2] += *sa * sx[i__ + 2];
#line 97 "saxpy.f"
	    sy[i__ + 3] += *sa * sx[i__ + 3];
#line 98 "saxpy.f"
	}
#line 99 "saxpy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 104 "saxpy.f"
	ix = 1;
#line 105 "saxpy.f"
	iy = 1;
#line 106 "saxpy.f"
	if (*incx < 0) {
#line 106 "saxpy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 106 "saxpy.f"
	}
#line 107 "saxpy.f"
	if (*incy < 0) {
#line 107 "saxpy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 107 "saxpy.f"
	}
#line 108 "saxpy.f"
	i__1 = *n;
#line 108 "saxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 109 "saxpy.f"
	    sy[iy] += *sa * sx[ix];
#line 110 "saxpy.f"
	    ix += *incx;
#line 111 "saxpy.f"
	    iy += *incy;
#line 112 "saxpy.f"
	}
#line 113 "saxpy.f"
    }
#line 114 "saxpy.f"
    return 0;
} /* saxpy_ */

