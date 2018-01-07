#line 1 "dscal.f"
/* dscal.f -- translated by f2c (version 20100827).
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

#line 1 "dscal.f"
/* > \brief \b DSCAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSCAL(N,DA,DX,INCX) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION DA */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DSCAL scales a vector by a constant. */
/* >    uses unrolled loops for increment equal to one. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup double_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, mp1, nincx;


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
#line 77 "dscal.f"
    /* Parameter adjustments */
#line 77 "dscal.f"
    --dx;
#line 77 "dscal.f"

#line 77 "dscal.f"
    /* Function Body */
#line 77 "dscal.f"
    if (*n <= 0 || *incx <= 0) {
#line 77 "dscal.f"
	return 0;
#line 77 "dscal.f"
    }
#line 78 "dscal.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */


/*        clean-up loop */

#line 85 "dscal.f"
	m = *n % 5;
#line 86 "dscal.f"
	if (m != 0) {
#line 87 "dscal.f"
	    i__1 = m;
#line 87 "dscal.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 88 "dscal.f"
		dx[i__] = *da * dx[i__];
#line 89 "dscal.f"
	    }
#line 90 "dscal.f"
	    if (*n < 5) {
#line 90 "dscal.f"
		return 0;
#line 90 "dscal.f"
	    }
#line 91 "dscal.f"
	}
#line 92 "dscal.f"
	mp1 = m + 1;
#line 93 "dscal.f"
	i__1 = *n;
#line 93 "dscal.f"
	for (i__ = mp1; i__ <= i__1; i__ += 5) {
#line 94 "dscal.f"
	    dx[i__] = *da * dx[i__];
#line 95 "dscal.f"
	    dx[i__ + 1] = *da * dx[i__ + 1];
#line 96 "dscal.f"
	    dx[i__ + 2] = *da * dx[i__ + 2];
#line 97 "dscal.f"
	    dx[i__ + 3] = *da * dx[i__ + 3];
#line 98 "dscal.f"
	    dx[i__ + 4] = *da * dx[i__ + 4];
#line 99 "dscal.f"
	}
#line 100 "dscal.f"
    } else {

/*        code for increment not equal to 1 */

#line 104 "dscal.f"
	nincx = *n * *incx;
#line 105 "dscal.f"
	i__1 = nincx;
#line 105 "dscal.f"
	i__2 = *incx;
#line 105 "dscal.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 106 "dscal.f"
	    dx[i__] = *da * dx[i__];
#line 107 "dscal.f"
	}
#line 108 "dscal.f"
    }
#line 109 "dscal.f"
    return 0;
} /* dscal_ */

