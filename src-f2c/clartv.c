#line 1 "clartv.f"
/* clartv.f -- translated by f2c (version 20100827).
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

#line 1 "clartv.f"
/* > \brief \b CLARTV applies a vector of plane rotations with real cosines and complex sines to the elements 
of a pair of vectors. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARTV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clartv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clartv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clartv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARTV( N, X, INCX, Y, INCY, C, S, INCC ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCC, INCX, INCY, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( * ) */
/*       COMPLEX            S( * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARTV applies a vector of complex plane rotations with real cosines */
/* > to elements of the complex vectors x and y. For i = 1,2,...,n */
/* > */
/* >    ( x(i) ) := (        c(i)   s(i) ) ( x(i) ) */
/* >    ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) ) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of plane rotations to be applied. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (1+(N-1)*INCX) */
/* >          The vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension (1+(N-1)*INCY) */
/* >          The vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >          The increment between elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL array, dimension (1+(N-1)*INCC) */
/* >          The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is COMPLEX array, dimension (1+(N-1)*INCC) */
/* >          The sines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCC */
/* > \verbatim */
/* >          INCC is INTEGER */
/* >          The increment between elements of C and S. INCC > 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clartv_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, doublecomplex *s, 
	integer *incc)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, ic, ix, iy;
    static doublecomplex xi, yi;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
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
/*     .. Executable Statements .. */

#line 134 "clartv.f"
    /* Parameter adjustments */
#line 134 "clartv.f"
    --s;
#line 134 "clartv.f"
    --c__;
#line 134 "clartv.f"
    --y;
#line 134 "clartv.f"
    --x;
#line 134 "clartv.f"

#line 134 "clartv.f"
    /* Function Body */
#line 134 "clartv.f"
    ix = 1;
#line 135 "clartv.f"
    iy = 1;
#line 136 "clartv.f"
    ic = 1;
#line 137 "clartv.f"
    i__1 = *n;
#line 137 "clartv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 138 "clartv.f"
	i__2 = ix;
#line 138 "clartv.f"
	xi.r = x[i__2].r, xi.i = x[i__2].i;
#line 139 "clartv.f"
	i__2 = iy;
#line 139 "clartv.f"
	yi.r = y[i__2].r, yi.i = y[i__2].i;
#line 140 "clartv.f"
	i__2 = ix;
#line 140 "clartv.f"
	i__3 = ic;
#line 140 "clartv.f"
	z__2.r = c__[i__3] * xi.r, z__2.i = c__[i__3] * xi.i;
#line 140 "clartv.f"
	i__4 = ic;
#line 140 "clartv.f"
	z__3.r = s[i__4].r * yi.r - s[i__4].i * yi.i, z__3.i = s[i__4].r * 
		yi.i + s[i__4].i * yi.r;
#line 140 "clartv.f"
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 140 "clartv.f"
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 141 "clartv.f"
	i__2 = iy;
#line 141 "clartv.f"
	i__3 = ic;
#line 141 "clartv.f"
	z__2.r = c__[i__3] * yi.r, z__2.i = c__[i__3] * yi.i;
#line 141 "clartv.f"
	d_cnjg(&z__4, &s[ic]);
#line 141 "clartv.f"
	z__3.r = z__4.r * xi.r - z__4.i * xi.i, z__3.i = z__4.r * xi.i + 
		z__4.i * xi.r;
#line 141 "clartv.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 141 "clartv.f"
	y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 142 "clartv.f"
	ix += *incx;
#line 143 "clartv.f"
	iy += *incy;
#line 144 "clartv.f"
	ic += *incc;
#line 145 "clartv.f"
/* L10: */
#line 145 "clartv.f"
    }
#line 146 "clartv.f"
    return 0;

/*     End of CLARTV */

} /* clartv_ */

