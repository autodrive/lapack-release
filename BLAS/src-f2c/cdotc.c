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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

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


/*  -- Reference BLAS level1 routine (version 3.6.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 76 "cdotc.f"
    /* Parameter adjustments */
#line 76 "cdotc.f"
    --cy;
#line 76 "cdotc.f"
    --cx;
#line 76 "cdotc.f"

#line 76 "cdotc.f"
    /* Function Body */
#line 76 "cdotc.f"
    ctemp.r = 0., ctemp.i = 0.;
#line 77 "cdotc.f"
     ret_val->r = 0.,  ret_val->i = 0.;
#line 78 "cdotc.f"
    if (*n <= 0) {
#line 78 "cdotc.f"
	return ;
#line 78 "cdotc.f"
    }
#line 79 "cdotc.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 83 "cdotc.f"
	i__1 = *n;
#line 83 "cdotc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 84 "cdotc.f"
	    d_cnjg(&z__3, &cx[i__]);
#line 84 "cdotc.f"
	    i__2 = i__;
#line 84 "cdotc.f"
	    z__2.r = z__3.r * cy[i__2].r - z__3.i * cy[i__2].i, z__2.i = 
		    z__3.r * cy[i__2].i + z__3.i * cy[i__2].r;
#line 84 "cdotc.f"
	    z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
#line 84 "cdotc.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 85 "cdotc.f"
	}
#line 86 "cdotc.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 91 "cdotc.f"
	ix = 1;
#line 92 "cdotc.f"
	iy = 1;
#line 93 "cdotc.f"
	if (*incx < 0) {
#line 93 "cdotc.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 93 "cdotc.f"
	}
#line 94 "cdotc.f"
	if (*incy < 0) {
#line 94 "cdotc.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 94 "cdotc.f"
	}
#line 95 "cdotc.f"
	i__1 = *n;
#line 95 "cdotc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 96 "cdotc.f"
	    d_cnjg(&z__3, &cx[ix]);
#line 96 "cdotc.f"
	    i__2 = iy;
#line 96 "cdotc.f"
	    z__2.r = z__3.r * cy[i__2].r - z__3.i * cy[i__2].i, z__2.i = 
		    z__3.r * cy[i__2].i + z__3.i * cy[i__2].r;
#line 96 "cdotc.f"
	    z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
#line 96 "cdotc.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 97 "cdotc.f"
	    ix += *incx;
#line 98 "cdotc.f"
	    iy += *incy;
#line 99 "cdotc.f"
	}
#line 100 "cdotc.f"
    }
#line 101 "cdotc.f"
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
#line 102 "cdotc.f"
    return ;
} /* cdotc_ */

