#line 1 "zdotc.f"
/* zdotc.f -- translated by f2c (version 20100827).
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

#line 1 "zdotc.f"
/* > \brief \b ZDOTC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 ZX(*),ZY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZDOTC forms the dot product of two complex vectors */
/* >      ZDOTC = X^H * Y */
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
/* > \param[in] ZX */
/* > \verbatim */
/* >          ZX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of ZX */
/* > \endverbatim */
/* > */
/* > \param[in] ZY */
/* > \verbatim */
/* >          ZY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of ZY */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complex16_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, 3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Double Complex */ VOID zdotc_(doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ztemp;


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
#line 107 "zdotc.f"
    /* Parameter adjustments */
#line 107 "zdotc.f"
    --zy;
#line 107 "zdotc.f"
    --zx;
#line 107 "zdotc.f"

#line 107 "zdotc.f"
    /* Function Body */
#line 107 "zdotc.f"
    ztemp.r = 0., ztemp.i = 0.;
#line 108 "zdotc.f"
     ret_val->r = 0.,  ret_val->i = 0.;
#line 109 "zdotc.f"
    if (*n <= 0) {
#line 109 "zdotc.f"
	return ;
#line 109 "zdotc.f"
    }
#line 110 "zdotc.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 114 "zdotc.f"
	i__1 = *n;
#line 114 "zdotc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 115 "zdotc.f"
	    d_cnjg(&z__3, &zx[i__]);
#line 115 "zdotc.f"
	    i__2 = i__;
#line 115 "zdotc.f"
	    z__2.r = z__3.r * zy[i__2].r - z__3.i * zy[i__2].i, z__2.i = 
		    z__3.r * zy[i__2].i + z__3.i * zy[i__2].r;
#line 115 "zdotc.f"
	    z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
#line 115 "zdotc.f"
	    ztemp.r = z__1.r, ztemp.i = z__1.i;
#line 116 "zdotc.f"
	}
#line 117 "zdotc.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 122 "zdotc.f"
	ix = 1;
#line 123 "zdotc.f"
	iy = 1;
#line 124 "zdotc.f"
	if (*incx < 0) {
#line 124 "zdotc.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 124 "zdotc.f"
	}
#line 125 "zdotc.f"
	if (*incy < 0) {
#line 125 "zdotc.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 125 "zdotc.f"
	}
#line 126 "zdotc.f"
	i__1 = *n;
#line 126 "zdotc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 127 "zdotc.f"
	    d_cnjg(&z__3, &zx[ix]);
#line 127 "zdotc.f"
	    i__2 = iy;
#line 127 "zdotc.f"
	    z__2.r = z__3.r * zy[i__2].r - z__3.i * zy[i__2].i, z__2.i = 
		    z__3.r * zy[i__2].i + z__3.i * zy[i__2].r;
#line 127 "zdotc.f"
	    z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
#line 127 "zdotc.f"
	    ztemp.r = z__1.r, ztemp.i = z__1.i;
#line 128 "zdotc.f"
	    ix += *incx;
#line 129 "zdotc.f"
	    iy += *incy;
#line 130 "zdotc.f"
	}
#line 131 "zdotc.f"
    }
#line 132 "zdotc.f"
     ret_val->r = ztemp.r,  ret_val->i = ztemp.i;
#line 133 "zdotc.f"
    return ;
} /* zdotc_ */

