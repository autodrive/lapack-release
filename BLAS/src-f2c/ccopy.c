#line 1 "ccopy.f"
/* ccopy.f -- translated by f2c (version 20100827).
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

#line 1 "ccopy.f"
/* > \brief \b CCOPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CCOPY(N,CX,INCX,CY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX CX(*),CY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CCOPY copies a vector x to a vector y. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex_blas_level1 */

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
/* Subroutine */ int ccopy_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ix, iy;


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
#line 70 "ccopy.f"
    /* Parameter adjustments */
#line 70 "ccopy.f"
    --cy;
#line 70 "ccopy.f"
    --cx;
#line 70 "ccopy.f"

#line 70 "ccopy.f"
    /* Function Body */
#line 70 "ccopy.f"
    if (*n <= 0) {
#line 70 "ccopy.f"
	return 0;
#line 70 "ccopy.f"
    }
#line 71 "ccopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 75 "ccopy.f"
	i__1 = *n;
#line 75 "ccopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 76 "ccopy.f"
	    i__2 = i__;
#line 76 "ccopy.f"
	    i__3 = i__;
#line 76 "ccopy.f"
	    cy[i__2].r = cx[i__3].r, cy[i__2].i = cx[i__3].i;
#line 77 "ccopy.f"
	}
#line 78 "ccopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 83 "ccopy.f"
	ix = 1;
#line 84 "ccopy.f"
	iy = 1;
#line 85 "ccopy.f"
	if (*incx < 0) {
#line 85 "ccopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 85 "ccopy.f"
	}
#line 86 "ccopy.f"
	if (*incy < 0) {
#line 86 "ccopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 86 "ccopy.f"
	}
#line 87 "ccopy.f"
	i__1 = *n;
#line 87 "ccopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 88 "ccopy.f"
	    i__2 = iy;
#line 88 "ccopy.f"
	    i__3 = ix;
#line 88 "ccopy.f"
	    cy[i__2].r = cx[i__3].r, cy[i__2].i = cx[i__3].i;
#line 89 "ccopy.f"
	    ix += *incx;
#line 90 "ccopy.f"
	    iy += *incy;
#line 91 "ccopy.f"
	}
#line 92 "ccopy.f"
    }
#line 93 "ccopy.f"
    return 0;
} /* ccopy_ */

