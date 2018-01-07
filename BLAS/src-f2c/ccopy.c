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

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] CX */
/* > \verbatim */
/* >          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of CX */
/* > \endverbatim */
/* > */
/* > \param[out] CY */
/* > \verbatim */
/* >          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of CY */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

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
#line 101 "ccopy.f"
    /* Parameter adjustments */
#line 101 "ccopy.f"
    --cy;
#line 101 "ccopy.f"
    --cx;
#line 101 "ccopy.f"

#line 101 "ccopy.f"
    /* Function Body */
#line 101 "ccopy.f"
    if (*n <= 0) {
#line 101 "ccopy.f"
	return 0;
#line 101 "ccopy.f"
    }
#line 102 "ccopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 106 "ccopy.f"
	i__1 = *n;
#line 106 "ccopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 107 "ccopy.f"
	    i__2 = i__;
#line 107 "ccopy.f"
	    i__3 = i__;
#line 107 "ccopy.f"
	    cy[i__2].r = cx[i__3].r, cy[i__2].i = cx[i__3].i;
#line 108 "ccopy.f"
	}
#line 109 "ccopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 114 "ccopy.f"
	ix = 1;
#line 115 "ccopy.f"
	iy = 1;
#line 116 "ccopy.f"
	if (*incx < 0) {
#line 116 "ccopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 116 "ccopy.f"
	}
#line 117 "ccopy.f"
	if (*incy < 0) {
#line 117 "ccopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 117 "ccopy.f"
	}
#line 118 "ccopy.f"
	i__1 = *n;
#line 118 "ccopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 119 "ccopy.f"
	    i__2 = iy;
#line 119 "ccopy.f"
	    i__3 = ix;
#line 119 "ccopy.f"
	    cy[i__2].r = cx[i__3].r, cy[i__2].i = cx[i__3].i;
#line 120 "ccopy.f"
	    ix += *incx;
#line 121 "ccopy.f"
	    iy += *incy;
#line 122 "ccopy.f"
	}
#line 123 "ccopy.f"
    }
#line 124 "ccopy.f"
    return 0;
} /* ccopy_ */

