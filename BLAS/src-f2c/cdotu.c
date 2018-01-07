#line 1 "cdotu.f"
/* cdotu.f -- translated by f2c (version 20100827).
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

#line 1 "cdotu.f"
/* > \brief \b CDOTU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY) */

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
/* > CDOTU forms the dot product of two complex vectors */
/* >      CDOTU = X^T * Y */
/* > */
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
/* >          CX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of CX */
/* > \endverbatim */
/* > */
/* > \param[in] CY */
/* > \verbatim */
/* >          CY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
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
/* Double Complex */ VOID cdotu_(doublecomplex * ret_val, integer *n, 
	doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


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
#line 104 "cdotu.f"
    /* Parameter adjustments */
#line 104 "cdotu.f"
    --cy;
#line 104 "cdotu.f"
    --cx;
#line 104 "cdotu.f"

#line 104 "cdotu.f"
    /* Function Body */
#line 104 "cdotu.f"
    ctemp.r = 0., ctemp.i = 0.;
#line 105 "cdotu.f"
     ret_val->r = 0.,  ret_val->i = 0.;
#line 106 "cdotu.f"
    if (*n <= 0) {
#line 106 "cdotu.f"
	return ;
#line 106 "cdotu.f"
    }
#line 107 "cdotu.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 111 "cdotu.f"
	i__1 = *n;
#line 111 "cdotu.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 112 "cdotu.f"
	    i__2 = i__;
#line 112 "cdotu.f"
	    i__3 = i__;
#line 112 "cdotu.f"
	    z__2.r = cx[i__2].r * cy[i__3].r - cx[i__2].i * cy[i__3].i, 
		    z__2.i = cx[i__2].r * cy[i__3].i + cx[i__2].i * cy[i__3]
		    .r;
#line 112 "cdotu.f"
	    z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
#line 112 "cdotu.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 113 "cdotu.f"
	}
#line 114 "cdotu.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 119 "cdotu.f"
	ix = 1;
#line 120 "cdotu.f"
	iy = 1;
#line 121 "cdotu.f"
	if (*incx < 0) {
#line 121 "cdotu.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 121 "cdotu.f"
	}
#line 122 "cdotu.f"
	if (*incy < 0) {
#line 122 "cdotu.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 122 "cdotu.f"
	}
#line 123 "cdotu.f"
	i__1 = *n;
#line 123 "cdotu.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 124 "cdotu.f"
	    i__2 = ix;
#line 124 "cdotu.f"
	    i__3 = iy;
#line 124 "cdotu.f"
	    z__2.r = cx[i__2].r * cy[i__3].r - cx[i__2].i * cy[i__3].i, 
		    z__2.i = cx[i__2].r * cy[i__3].i + cx[i__2].i * cy[i__3]
		    .r;
#line 124 "cdotu.f"
	    z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
#line 124 "cdotu.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 125 "cdotu.f"
	    ix += *incx;
#line 126 "cdotu.f"
	    iy += *incy;
#line 127 "cdotu.f"
	}
#line 128 "cdotu.f"
    }
#line 129 "cdotu.f"
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
#line 130 "cdotu.f"
    return ;
} /* cdotu_ */

