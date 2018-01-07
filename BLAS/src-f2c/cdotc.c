#line 1 "cdotc.f"
/* cdotc.f -- translated by f2c (version 20100827).
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

#line 1 "cdotc.f"
/* > \brief \b CDOTC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY) */

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
/* > CDOTC forms the dot product of two complex vectors */
/* >      CDOTC = X^H * Y */
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
/* >     jack dongarra, linpack,  3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Double Complex */ VOID cdotc_(doublecomplex * ret_val, integer *n, 
	doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

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
/*     .. Intrinsic Functions .. */
/*     .. */
#line 107 "cdotc.f"
    /* Parameter adjustments */
#line 107 "cdotc.f"
    --cy;
#line 107 "cdotc.f"
    --cx;
#line 107 "cdotc.f"

#line 107 "cdotc.f"
    /* Function Body */
#line 107 "cdotc.f"
    ctemp.r = 0., ctemp.i = 0.;
#line 108 "cdotc.f"
     ret_val->r = 0.,  ret_val->i = 0.;
#line 109 "cdotc.f"
    if (*n <= 0) {
#line 109 "cdotc.f"
	return ;
#line 109 "cdotc.f"
    }
#line 110 "cdotc.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 114 "cdotc.f"
	i__1 = *n;
#line 114 "cdotc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 115 "cdotc.f"
	    d_cnjg(&z__3, &cx[i__]);
#line 115 "cdotc.f"
	    i__2 = i__;
#line 115 "cdotc.f"
	    z__2.r = z__3.r * cy[i__2].r - z__3.i * cy[i__2].i, z__2.i = 
		    z__3.r * cy[i__2].i + z__3.i * cy[i__2].r;
#line 115 "cdotc.f"
	    z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
#line 115 "cdotc.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 116 "cdotc.f"
	}
#line 117 "cdotc.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 122 "cdotc.f"
	ix = 1;
#line 123 "cdotc.f"
	iy = 1;
#line 124 "cdotc.f"
	if (*incx < 0) {
#line 124 "cdotc.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 124 "cdotc.f"
	}
#line 125 "cdotc.f"
	if (*incy < 0) {
#line 125 "cdotc.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 125 "cdotc.f"
	}
#line 126 "cdotc.f"
	i__1 = *n;
#line 126 "cdotc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 127 "cdotc.f"
	    d_cnjg(&z__3, &cx[ix]);
#line 127 "cdotc.f"
	    i__2 = iy;
#line 127 "cdotc.f"
	    z__2.r = z__3.r * cy[i__2].r - z__3.i * cy[i__2].i, z__2.i = 
		    z__3.r * cy[i__2].i + z__3.i * cy[i__2].r;
#line 127 "cdotc.f"
	    z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
#line 127 "cdotc.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 128 "cdotc.f"
	    ix += *incx;
#line 129 "cdotc.f"
	    iy += *incy;
#line 130 "cdotc.f"
	}
#line 131 "cdotc.f"
    }
#line 132 "cdotc.f"
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
#line 133 "cdotc.f"
    return ;
} /* cdotc_ */

