#line 1 "csscal.f"
/* csscal.f -- translated by f2c (version 20100827).
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

#line 1 "csscal.f"
/* > \brief \b CSSCAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSSCAL(N,SA,CX,INCX) */

/*       .. Scalar Arguments .. */
/*       REAL SA */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX CX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CSSCAL scales a complex vector by a real constant. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int csscal_(integer *n, doublereal *sa, doublecomplex *cx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

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
#line 76 "csscal.f"
    /* Parameter adjustments */
#line 76 "csscal.f"
    --cx;
#line 76 "csscal.f"

#line 76 "csscal.f"
    /* Function Body */
#line 76 "csscal.f"
    if (*n <= 0 || *incx <= 0) {
#line 76 "csscal.f"
	return 0;
#line 76 "csscal.f"
    }
#line 77 "csscal.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 81 "csscal.f"
	i__1 = *n;
#line 81 "csscal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 82 "csscal.f"
	    i__2 = i__;
#line 82 "csscal.f"
	    i__3 = i__;
#line 82 "csscal.f"
	    d__1 = *sa * cx[i__3].r;
#line 82 "csscal.f"
	    d__2 = *sa * d_imag(&cx[i__]);
#line 82 "csscal.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 82 "csscal.f"
	    cx[i__2].r = z__1.r, cx[i__2].i = z__1.i;
#line 83 "csscal.f"
	}
#line 84 "csscal.f"
    } else {

/*        code for increment not equal to 1 */

#line 88 "csscal.f"
	nincx = *n * *incx;
#line 89 "csscal.f"
	i__1 = nincx;
#line 89 "csscal.f"
	i__2 = *incx;
#line 89 "csscal.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 90 "csscal.f"
	    i__3 = i__;
#line 90 "csscal.f"
	    i__4 = i__;
#line 90 "csscal.f"
	    d__1 = *sa * cx[i__4].r;
#line 90 "csscal.f"
	    d__2 = *sa * d_imag(&cx[i__]);
#line 90 "csscal.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 90 "csscal.f"
	    cx[i__3].r = z__1.r, cx[i__3].i = z__1.i;
#line 91 "csscal.f"
	}
#line 92 "csscal.f"
    }
#line 93 "csscal.f"
    return 0;
} /* csscal_ */

