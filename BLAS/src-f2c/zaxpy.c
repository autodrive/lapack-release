#line 1 "zaxpy.f"
/* zaxpy.f -- translated by f2c (version 20100827).
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

#line 1 "zaxpy.f"
/* > \brief \b ZAXPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ZA */
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
/* >    ZAXPY constant times a vector plus a vector. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

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
/* Subroutine */ int zaxpy_(integer *n, doublecomplex *za, doublecomplex *zx, 
	integer *incx, doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, ix, iy;
    extern doublereal dcabs1_(doublecomplex *);


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
/*     .. External Functions .. */
/*     .. */
#line 76 "zaxpy.f"
    /* Parameter adjustments */
#line 76 "zaxpy.f"
    --zy;
#line 76 "zaxpy.f"
    --zx;
#line 76 "zaxpy.f"

#line 76 "zaxpy.f"
    /* Function Body */
#line 76 "zaxpy.f"
    if (*n <= 0) {
#line 76 "zaxpy.f"
	return 0;
#line 76 "zaxpy.f"
    }
#line 77 "zaxpy.f"
    if (dcabs1_(za) == 0.) {
#line 77 "zaxpy.f"
	return 0;
#line 77 "zaxpy.f"
    }
#line 78 "zaxpy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 82 "zaxpy.f"
	i__1 = *n;
#line 82 "zaxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 83 "zaxpy.f"
	    i__2 = i__;
#line 83 "zaxpy.f"
	    i__3 = i__;
#line 83 "zaxpy.f"
	    i__4 = i__;
#line 83 "zaxpy.f"
	    z__2.r = za->r * zx[i__4].r - za->i * zx[i__4].i, z__2.i = za->r *
		     zx[i__4].i + za->i * zx[i__4].r;
#line 83 "zaxpy.f"
	    z__1.r = zy[i__3].r + z__2.r, z__1.i = zy[i__3].i + z__2.i;
#line 83 "zaxpy.f"
	    zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
#line 84 "zaxpy.f"
	}
#line 85 "zaxpy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 90 "zaxpy.f"
	ix = 1;
#line 91 "zaxpy.f"
	iy = 1;
#line 92 "zaxpy.f"
	if (*incx < 0) {
#line 92 "zaxpy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 92 "zaxpy.f"
	}
#line 93 "zaxpy.f"
	if (*incy < 0) {
#line 93 "zaxpy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 93 "zaxpy.f"
	}
#line 94 "zaxpy.f"
	i__1 = *n;
#line 94 "zaxpy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 95 "zaxpy.f"
	    i__2 = iy;
#line 95 "zaxpy.f"
	    i__3 = iy;
#line 95 "zaxpy.f"
	    i__4 = ix;
#line 95 "zaxpy.f"
	    z__2.r = za->r * zx[i__4].r - za->i * zx[i__4].i, z__2.i = za->r *
		     zx[i__4].i + za->i * zx[i__4].r;
#line 95 "zaxpy.f"
	    z__1.r = zy[i__3].r + z__2.r, z__1.i = zy[i__3].i + z__2.i;
#line 95 "zaxpy.f"
	    zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
#line 96 "zaxpy.f"
	    ix += *incx;
#line 97 "zaxpy.f"
	    iy += *incy;
#line 98 "zaxpy.f"
	}
#line 99 "zaxpy.f"
    }

#line 101 "zaxpy.f"
    return 0;
} /* zaxpy_ */

