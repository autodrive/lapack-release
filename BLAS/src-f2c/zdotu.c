#line 1 "zdotu.f"
/* zdotu.f -- translated by f2c (version 20100827).
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

#line 1 "zdotu.f"
/* > \brief \b ZDOTU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY) */

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
/* >    ZDOTU forms the dot product of two vectors. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
/* Double Complex */ VOID zdotu_(doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ztemp;


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
#line 71 "zdotu.f"
    /* Parameter adjustments */
#line 71 "zdotu.f"
    --zy;
#line 71 "zdotu.f"
    --zx;
#line 71 "zdotu.f"

#line 71 "zdotu.f"
    /* Function Body */
#line 71 "zdotu.f"
    ztemp.r = 0., ztemp.i = 0.;
#line 72 "zdotu.f"
     ret_val->r = 0.,  ret_val->i = 0.;
#line 73 "zdotu.f"
    if (*n <= 0) {
#line 73 "zdotu.f"
	return ;
#line 73 "zdotu.f"
    }
#line 74 "zdotu.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 78 "zdotu.f"
	i__1 = *n;
#line 78 "zdotu.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 79 "zdotu.f"
	    i__2 = i__;
#line 79 "zdotu.f"
	    i__3 = i__;
#line 79 "zdotu.f"
	    z__2.r = zx[i__2].r * zy[i__3].r - zx[i__2].i * zy[i__3].i, 
		    z__2.i = zx[i__2].r * zy[i__3].i + zx[i__2].i * zy[i__3]
		    .r;
#line 79 "zdotu.f"
	    z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
#line 79 "zdotu.f"
	    ztemp.r = z__1.r, ztemp.i = z__1.i;
#line 80 "zdotu.f"
	}
#line 81 "zdotu.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 86 "zdotu.f"
	ix = 1;
#line 87 "zdotu.f"
	iy = 1;
#line 88 "zdotu.f"
	if (*incx < 0) {
#line 88 "zdotu.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 88 "zdotu.f"
	}
#line 89 "zdotu.f"
	if (*incy < 0) {
#line 89 "zdotu.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 89 "zdotu.f"
	}
#line 90 "zdotu.f"
	i__1 = *n;
#line 90 "zdotu.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 91 "zdotu.f"
	    i__2 = ix;
#line 91 "zdotu.f"
	    i__3 = iy;
#line 91 "zdotu.f"
	    z__2.r = zx[i__2].r * zy[i__3].r - zx[i__2].i * zy[i__3].i, 
		    z__2.i = zx[i__2].r * zy[i__3].i + zx[i__2].i * zy[i__3]
		    .r;
#line 91 "zdotu.f"
	    z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
#line 91 "zdotu.f"
	    ztemp.r = z__1.r, ztemp.i = z__1.i;
#line 92 "zdotu.f"
	    ix += *incx;
#line 93 "zdotu.f"
	    iy += *incy;
#line 94 "zdotu.f"
	}
#line 95 "zdotu.f"
    }
#line 96 "zdotu.f"
     ret_val->r = ztemp.r,  ret_val->i = ztemp.i;
#line 97 "zdotu.f"
    return ;
} /* zdotu_ */

