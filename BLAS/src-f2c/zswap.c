#line 1 "zswap.f"
/* zswap.f -- translated by f2c (version 20100827).
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

#line 1 "zswap.f"
/* > \brief \b ZSWAP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY) */

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
/* >    ZSWAP interchanges two vectors. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZX */
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
/* > \param[in,out] ZY */
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
/* >     jack dongarra, 3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zswap_(integer *n, doublecomplex *zx, integer *incx, 
	doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

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
#line 102 "zswap.f"
    /* Parameter adjustments */
#line 102 "zswap.f"
    --zy;
#line 102 "zswap.f"
    --zx;
#line 102 "zswap.f"

#line 102 "zswap.f"
    /* Function Body */
#line 102 "zswap.f"
    if (*n <= 0) {
#line 102 "zswap.f"
	return 0;
#line 102 "zswap.f"
    }
#line 103 "zswap.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */
#line 106 "zswap.f"
	i__1 = *n;
#line 106 "zswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 107 "zswap.f"
	    i__2 = i__;
#line 107 "zswap.f"
	    ztemp.r = zx[i__2].r, ztemp.i = zx[i__2].i;
#line 108 "zswap.f"
	    i__2 = i__;
#line 108 "zswap.f"
	    i__3 = i__;
#line 108 "zswap.f"
	    zx[i__2].r = zy[i__3].r, zx[i__2].i = zy[i__3].i;
#line 109 "zswap.f"
	    i__2 = i__;
#line 109 "zswap.f"
	    zy[i__2].r = ztemp.r, zy[i__2].i = ztemp.i;
#line 110 "zswap.f"
	}
#line 111 "zswap.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 116 "zswap.f"
	ix = 1;
#line 117 "zswap.f"
	iy = 1;
#line 118 "zswap.f"
	if (*incx < 0) {
#line 118 "zswap.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 118 "zswap.f"
	}
#line 119 "zswap.f"
	if (*incy < 0) {
#line 119 "zswap.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 119 "zswap.f"
	}
#line 120 "zswap.f"
	i__1 = *n;
#line 120 "zswap.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 121 "zswap.f"
	    i__2 = ix;
#line 121 "zswap.f"
	    ztemp.r = zx[i__2].r, ztemp.i = zx[i__2].i;
#line 122 "zswap.f"
	    i__2 = ix;
#line 122 "zswap.f"
	    i__3 = iy;
#line 122 "zswap.f"
	    zx[i__2].r = zy[i__3].r, zx[i__2].i = zy[i__3].i;
#line 123 "zswap.f"
	    i__2 = iy;
#line 123 "zswap.f"
	    zy[i__2].r = ztemp.r, zy[i__2].i = ztemp.i;
#line 124 "zswap.f"
	    ix += *incx;
#line 125 "zswap.f"
	    iy += *incy;
#line 126 "zswap.f"
	}
#line 127 "zswap.f"
    }
#line 128 "zswap.f"
    return 0;
} /* zswap_ */

