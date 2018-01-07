#line 1 "cscal.f"
/* cscal.f -- translated by f2c (version 20100827).
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

#line 1 "cscal.f"
/* > \brief \b CSCAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSCAL(N,CA,CX,INCX) */

/*       .. Scalar Arguments .. */
/*       COMPLEX CA */
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
/* >    CSCAL scales a vector by a constant. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] CA */
/* > \verbatim */
/* >          CA is COMPLEX */
/* >           On entry, CA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* >          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of CX */
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
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cscal_(integer *n, doublecomplex *ca, doublecomplex *cx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, nincx;


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
#line 99 "cscal.f"
    /* Parameter adjustments */
#line 99 "cscal.f"
    --cx;
#line 99 "cscal.f"

#line 99 "cscal.f"
    /* Function Body */
#line 99 "cscal.f"
    if (*n <= 0 || *incx <= 0) {
#line 99 "cscal.f"
	return 0;
#line 99 "cscal.f"
    }
#line 100 "cscal.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 104 "cscal.f"
	i__1 = *n;
#line 104 "cscal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 105 "cscal.f"
	    i__2 = i__;
#line 105 "cscal.f"
	    i__3 = i__;
#line 105 "cscal.f"
	    z__1.r = ca->r * cx[i__3].r - ca->i * cx[i__3].i, z__1.i = ca->r *
		     cx[i__3].i + ca->i * cx[i__3].r;
#line 105 "cscal.f"
	    cx[i__2].r = z__1.r, cx[i__2].i = z__1.i;
#line 106 "cscal.f"
	}
#line 107 "cscal.f"
    } else {

/*        code for increment not equal to 1 */

#line 111 "cscal.f"
	nincx = *n * *incx;
#line 112 "cscal.f"
	i__1 = nincx;
#line 112 "cscal.f"
	i__2 = *incx;
#line 112 "cscal.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 113 "cscal.f"
	    i__3 = i__;
#line 113 "cscal.f"
	    i__4 = i__;
#line 113 "cscal.f"
	    z__1.r = ca->r * cx[i__4].r - ca->i * cx[i__4].i, z__1.i = ca->r *
		     cx[i__4].i + ca->i * cx[i__4].r;
#line 113 "cscal.f"
	    cx[i__3].r = z__1.r, cx[i__3].i = z__1.i;
#line 114 "cscal.f"
	}
#line 115 "cscal.f"
    }
#line 116 "cscal.f"
    return 0;
} /* cscal_ */

