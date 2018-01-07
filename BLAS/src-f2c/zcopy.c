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
/* >          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of ZX */
/* > \endverbatim */
/* > */
/* > \param[out] ZY */
/* > \verbatim */
/* >          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
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
#line 101 "zcopy.f"
    /* Parameter adjustments */
#line 101 "zcopy.f"
    --zy;
#line 101 "zcopy.f"
    --zx;
#line 101 "zcopy.f"

#line 101 "zcopy.f"
    /* Function Body */
#line 101 "zcopy.f"
    if (*n <= 0) {
#line 101 "zcopy.f"
	return 0;
#line 101 "zcopy.f"
    }
#line 102 "zcopy.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 106 "zcopy.f"
	i__1 = *n;
#line 106 "zcopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 107 "zcopy.f"
	    i__2 = i__;
#line 107 "zcopy.f"
	    i__3 = i__;
#line 107 "zcopy.f"
	    zy[i__2].r = zx[i__3].r, zy[i__2].i = zx[i__3].i;
#line 108 "zcopy.f"
	}
#line 109 "zcopy.f"
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

#line 114 "zcopy.f"
	ix = 1;
#line 115 "zcopy.f"
	iy = 1;
#line 116 "zcopy.f"
	if (*incx < 0) {
#line 116 "zcopy.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 116 "zcopy.f"
	}
#line 117 "zcopy.f"
	if (*incy < 0) {
#line 117 "zcopy.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 117 "zcopy.f"
	}
#line 118 "zcopy.f"
	i__1 = *n;
#line 118 "zcopy.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 119 "zcopy.f"
	    i__2 = iy;
#line 119 "zcopy.f"
	    i__3 = ix;
#line 119 "zcopy.f"
	    zy[i__2].r = zx[i__3].r, zy[i__2].i = zx[i__3].i;
#line 120 "zcopy.f"
	    ix += *incx;
#line 121 "zcopy.f"
	    iy += *incy;
#line 122 "zcopy.f"
	}
#line 123 "zcopy.f"
    }
#line 124 "zcopy.f"
    return 0;
} /* zcopy_ */

