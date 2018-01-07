#line 1 "zcopy.f"
/* zcopy.f -- translated by f2c (version 20100827).
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

#line 1 "zcopy.f"
/* > \brief \b ZCOPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY) */

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
/* >    ZCOPY copies a vector, x, to a vector, y. */
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
/* >     jack dongarra, linpack, 4/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zcopy_(integer *n, doublecomplex *zx, integer *incx, 
	doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ix, iy;


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
#line 70 "zcopy.f"
    /* Parameter adjustments */
#line 70 "zcopy.f"
    --zy;
#line 70 "zcopy.f"
    --zx;
#line 70 "zcopy.f"

#line 70 "zcopy.f"
    /* Function Body */
#line 70 "zcopy.f"
    if (*n <= 0) {
#line 70 "zcopy.f"
	return 0;
#line 70 "zcopy.f"
    }
#line 71 "zcopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 75 "zcopy.f"
	i__1 = *n;
#line 75 "zcopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 76 "zcopy.f"
	    i__2 = i__;
#line 76 "zcopy.f"
	    i__3 = i__;
#line 76 "zcopy.f"
	    zy[i__2].r = zx[i__3].r, zy[i__2].i = zx[i__3].i;
#line 77 "zcopy.f"
	}
#line 78 "zcopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 83 "zcopy.f"
	ix = 1;
#line 84 "zcopy.f"
	iy = 1;
#line 85 "zcopy.f"
	if (*incx < 0) {
#line 85 "zcopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 85 "zcopy.f"
	}
#line 86 "zcopy.f"
	if (*incy < 0) {
#line 86 "zcopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 86 "zcopy.f"
	}
#line 87 "zcopy.f"
	i__1 = *n;
#line 87 "zcopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 88 "zcopy.f"
	    i__2 = iy;
#line 88 "zcopy.f"
	    i__3 = ix;
#line 88 "zcopy.f"
	    zy[i__2].r = zx[i__3].r, zy[i__2].i = zx[i__3].i;
#line 89 "zcopy.f"
	    ix += *incx;
#line 90 "zcopy.f"
	    iy += *incy;
#line 91 "zcopy.f"
	}
#line 92 "zcopy.f"
    }
#line 93 "zcopy.f"
    return 0;
} /* zcopy_ */

