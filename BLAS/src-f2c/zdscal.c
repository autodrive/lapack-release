#line 1 "zdscal.f"
/* zdscal.f -- translated by f2c (version 20100827).
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

#line 1 "zdscal.f"
/* > \brief \b ZDSCAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZDSCAL(N,DA,ZX,INCX) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION DA */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 ZX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZDSCAL scales a vector by a constant. */
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
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zdscal_(integer *n, doublereal *da, doublecomplex *zx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, nincx;


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
/*     .. Intrinsic Functions .. */
/*     .. */
#line 76 "zdscal.f"
    /* Parameter adjustments */
#line 76 "zdscal.f"
    --zx;
#line 76 "zdscal.f"

#line 76 "zdscal.f"
    /* Function Body */
#line 76 "zdscal.f"
    if (*n <= 0 || *incx <= 0) {
#line 76 "zdscal.f"
	return 0;
#line 76 "zdscal.f"
    }
#line 77 "zdscal.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 81 "zdscal.f"
	i__1 = *n;
#line 81 "zdscal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 82 "zdscal.f"
	    i__2 = i__;
#line 82 "zdscal.f"
	    z__2.r = *da, z__2.i = 0.;
#line 82 "zdscal.f"
	    i__3 = i__;
#line 82 "zdscal.f"
	    z__1.r = z__2.r * zx[i__3].r - z__2.i * zx[i__3].i, z__1.i = 
		    z__2.r * zx[i__3].i + z__2.i * zx[i__3].r;
#line 82 "zdscal.f"
	    zx[i__2].r = z__1.r, zx[i__2].i = z__1.i;
#line 83 "zdscal.f"
	}
#line 84 "zdscal.f"
    } else {

/*        code for increment not equal to 1 */

#line 88 "zdscal.f"
	nincx = *n * *incx;
#line 89 "zdscal.f"
	i__1 = nincx;
#line 89 "zdscal.f"
	i__2 = *incx;
#line 89 "zdscal.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 90 "zdscal.f"
	    i__3 = i__;
#line 90 "zdscal.f"
	    z__2.r = *da, z__2.i = 0.;
#line 90 "zdscal.f"
	    i__4 = i__;
#line 90 "zdscal.f"
	    z__1.r = z__2.r * zx[i__4].r - z__2.i * zx[i__4].i, z__1.i = 
		    z__2.r * zx[i__4].i + z__2.i * zx[i__4].r;
#line 90 "zdscal.f"
	    zx[i__3].r = z__1.r, zx[i__3].i = z__1.i;
#line 91 "zdscal.f"
	}
#line 92 "zdscal.f"
    }
#line 93 "zdscal.f"
    return 0;
} /* zdscal_ */

