#line 1 "cswap.f"
/* cswap.f -- translated by f2c (version 20100827).
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

#line 1 "cswap.f"
/* > \brief \b CSWAP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSWAP(N,CX,INCX,CY,INCY) */

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
/* >   CSWAP interchanges two vectors. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
/* Subroutine */ int cswap_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


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
#line 71 "cswap.f"
    /* Parameter adjustments */
#line 71 "cswap.f"
    --cy;
#line 71 "cswap.f"
    --cx;
#line 71 "cswap.f"

#line 71 "cswap.f"
    /* Function Body */
#line 71 "cswap.f"
    if (*n <= 0) {
#line 71 "cswap.f"
	return 0;
#line 71 "cswap.f"
    }
#line 72 "cswap.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */
#line 75 "cswap.f"
	i__1 = *n;
#line 75 "cswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 76 "cswap.f"
	    i__2 = i__;
#line 76 "cswap.f"
	    ctemp.r = cx[i__2].r, ctemp.i = cx[i__2].i;
#line 77 "cswap.f"
	    i__2 = i__;
#line 77 "cswap.f"
	    i__3 = i__;
#line 77 "cswap.f"
	    cx[i__2].r = cy[i__3].r, cx[i__2].i = cy[i__3].i;
#line 78 "cswap.f"
	    i__2 = i__;
#line 78 "cswap.f"
	    cy[i__2].r = ctemp.r, cy[i__2].i = ctemp.i;
#line 79 "cswap.f"
	}
#line 80 "cswap.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 85 "cswap.f"
	ix = 1;
#line 86 "cswap.f"
	iy = 1;
#line 87 "cswap.f"
	if (*incx < 0) {
#line 87 "cswap.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 87 "cswap.f"
	}
#line 88 "cswap.f"
	if (*incy < 0) {
#line 88 "cswap.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 88 "cswap.f"
	}
#line 89 "cswap.f"
	i__1 = *n;
#line 89 "cswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 90 "cswap.f"
	    i__2 = ix;
#line 90 "cswap.f"
	    ctemp.r = cx[i__2].r, ctemp.i = cx[i__2].i;
#line 91 "cswap.f"
	    i__2 = ix;
#line 91 "cswap.f"
	    i__3 = iy;
#line 91 "cswap.f"
	    cx[i__2].r = cy[i__3].r, cx[i__2].i = cy[i__3].i;
#line 92 "cswap.f"
	    i__2 = iy;
#line 92 "cswap.f"
	    cy[i__2].r = ctemp.r, cy[i__2].i = ctemp.i;
#line 93 "cswap.f"
	    ix += *incx;
#line 94 "cswap.f"
	    iy += *incy;
#line 95 "cswap.f"
	}
#line 96 "cswap.f"
    }
#line 97 "cswap.f"
    return 0;
} /* cswap_ */

