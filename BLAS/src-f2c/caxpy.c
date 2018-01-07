#line 1 "caxpy.f"
/* caxpy.f -- translated by f2c (version 20100827).
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

#line 1 "caxpy.f"
/* > \brief \b CAXPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX CA */
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
/* >    CAXPY constant times a vector plus a vector. */
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
/* Subroutine */ int caxpy_(integer *n, doublecomplex *ca, doublecomplex *cx, 
	integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, ix, iy;
    extern doublereal scabs1_(doublecomplex *);


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
/*     .. External Functions .. */
/*     .. */
#line 76 "caxpy.f"
    /* Parameter adjustments */
#line 76 "caxpy.f"
    --cy;
#line 76 "caxpy.f"
    --cx;
#line 76 "caxpy.f"

#line 76 "caxpy.f"
    /* Function Body */
#line 76 "caxpy.f"
    if (*n <= 0) {
#line 76 "caxpy.f"
	return 0;
#line 76 "caxpy.f"
    }
#line 77 "caxpy.f"
    if (scabs1_(ca) == 0.) {
#line 77 "caxpy.f"
	return 0;
#line 77 "caxpy.f"
    }
#line 78 "caxpy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 82 "caxpy.f"
	i__1 = *n;
#line 82 "caxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 83 "caxpy.f"
	    i__2 = i__;
#line 83 "caxpy.f"
	    i__3 = i__;
#line 83 "caxpy.f"
	    i__4 = i__;
#line 83 "caxpy.f"
	    z__2.r = ca->r * cx[i__4].r - ca->i * cx[i__4].i, z__2.i = ca->r *
		     cx[i__4].i + ca->i * cx[i__4].r;
#line 83 "caxpy.f"
	    z__1.r = cy[i__3].r + z__2.r, z__1.i = cy[i__3].i + z__2.i;
#line 83 "caxpy.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 84 "caxpy.f"
	}
#line 85 "caxpy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 90 "caxpy.f"
	ix = 1;
#line 91 "caxpy.f"
	iy = 1;
#line 92 "caxpy.f"
	if (*incx < 0) {
#line 92 "caxpy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 92 "caxpy.f"
	}
#line 93 "caxpy.f"
	if (*incy < 0) {
#line 93 "caxpy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 93 "caxpy.f"
	}
#line 94 "caxpy.f"
	i__1 = *n;
#line 94 "caxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 95 "caxpy.f"
	    i__2 = iy;
#line 95 "caxpy.f"
	    i__3 = iy;
#line 95 "caxpy.f"
	    i__4 = ix;
#line 95 "caxpy.f"
	    z__2.r = ca->r * cx[i__4].r - ca->i * cx[i__4].i, z__2.i = ca->r *
		     cx[i__4].i + ca->i * cx[i__4].r;
#line 95 "caxpy.f"
	    z__1.r = cy[i__3].r + z__2.r, z__1.i = cy[i__3].i + z__2.i;
#line 95 "caxpy.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 96 "caxpy.f"
	    ix += *incx;
#line 97 "caxpy.f"
	    iy += *incy;
#line 98 "caxpy.f"
	}
#line 99 "caxpy.f"
    }

#line 101 "caxpy.f"
    return 0;
} /* caxpy_ */

