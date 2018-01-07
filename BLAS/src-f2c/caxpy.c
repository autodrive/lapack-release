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

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] CA */
/* > \verbatim */
/* >          CA is COMPLEX */
/* >           On entry, CA specifies the scalar alpha. */
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
/* > \param[in,out] CY */
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
/* Subroutine */ int caxpy_(integer *n, doublecomplex *ca, doublecomplex *cx, 
	integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, ix, iy;
    extern doublereal scabs1_(doublecomplex *);


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
/*     .. External Functions .. */
/*     .. */
#line 113 "caxpy.f"
    /* Parameter adjustments */
#line 113 "caxpy.f"
    --cy;
#line 113 "caxpy.f"
    --cx;
#line 113 "caxpy.f"

#line 113 "caxpy.f"
    /* Function Body */
#line 113 "caxpy.f"
    if (*n <= 0) {
#line 113 "caxpy.f"
	return 0;
#line 113 "caxpy.f"
    }
#line 114 "caxpy.f"
    if (scabs1_(ca) == 0.) {
#line 114 "caxpy.f"
	return 0;
#line 114 "caxpy.f"
    }
#line 115 "caxpy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 119 "caxpy.f"
	i__1 = *n;
#line 119 "caxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 120 "caxpy.f"
	    i__2 = i__;
#line 120 "caxpy.f"
	    i__3 = i__;
#line 120 "caxpy.f"
	    i__4 = i__;
#line 120 "caxpy.f"
	    z__2.r = ca->r * cx[i__4].r - ca->i * cx[i__4].i, z__2.i = ca->r *
		     cx[i__4].i + ca->i * cx[i__4].r;
#line 120 "caxpy.f"
	    z__1.r = cy[i__3].r + z__2.r, z__1.i = cy[i__3].i + z__2.i;
#line 120 "caxpy.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 121 "caxpy.f"
	}
#line 122 "caxpy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 127 "caxpy.f"
	ix = 1;
#line 128 "caxpy.f"
	iy = 1;
#line 129 "caxpy.f"
	if (*incx < 0) {
#line 129 "caxpy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 129 "caxpy.f"
	}
#line 130 "caxpy.f"
	if (*incy < 0) {
#line 130 "caxpy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 130 "caxpy.f"
	}
#line 131 "caxpy.f"
	i__1 = *n;
#line 131 "caxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 132 "caxpy.f"
	    i__2 = iy;
#line 132 "caxpy.f"
	    i__3 = iy;
#line 132 "caxpy.f"
	    i__4 = ix;
#line 132 "caxpy.f"
	    z__2.r = ca->r * cx[i__4].r - ca->i * cx[i__4].i, z__2.i = ca->r *
		     cx[i__4].i + ca->i * cx[i__4].r;
#line 132 "caxpy.f"
	    z__1.r = cy[i__3].r + z__2.r, z__1.i = cy[i__3].i + z__2.i;
#line 132 "caxpy.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 133 "caxpy.f"
	    ix += *incx;
#line 134 "caxpy.f"
	    iy += *incy;
#line 135 "caxpy.f"
	}
#line 136 "caxpy.f"
    }

#line 138 "caxpy.f"
    return 0;
} /* caxpy_ */

