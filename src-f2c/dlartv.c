#line 1 "dlartv.f"
/* dlartv.f -- translated by f2c (version 20100827).
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

#line 1 "dlartv.f"
/* > \brief \b DLARTV applies a vector of plane rotations with real cosines and real sines to the elements of 
a pair of vectors. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARTV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARTV( N, X, INCX, Y, INCY, C, S, INCC ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCC, INCX, INCY, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ), S( * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARTV applies a vector of real plane rotations to elements of the */
/* > real vectors x and y. For i = 1,2,...,n */
/* > */
/* >    ( x(i) ) := (  c(i)  s(i) ) ( x(i) ) */
/* >    ( y(i) )    ( -s(i)  c(i) ) ( y(i) ) */
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
/* >          X is DOUBLE PRECISION array, */
/* >                         dimension (1+(N-1)*INCX) */
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
/* >          Y is DOUBLE PRECISION array, */
/* >                         dimension (1+(N-1)*INCY) */
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
/* >          C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
/* >          The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
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

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlartv_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *c__, doublereal *s, integer 
	*incc)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ic, ix, iy;
    static doublereal xi, yi;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 131 "dlartv.f"
    /* Parameter adjustments */
#line 131 "dlartv.f"
    --s;
#line 131 "dlartv.f"
    --c__;
#line 131 "dlartv.f"
    --y;
#line 131 "dlartv.f"
    --x;
#line 131 "dlartv.f"

#line 131 "dlartv.f"
    /* Function Body */
#line 131 "dlartv.f"
    ix = 1;
#line 132 "dlartv.f"
    iy = 1;
#line 133 "dlartv.f"
    ic = 1;
#line 134 "dlartv.f"
    i__1 = *n;
#line 134 "dlartv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 135 "dlartv.f"
	xi = x[ix];
#line 136 "dlartv.f"
	yi = y[iy];
#line 137 "dlartv.f"
	x[ix] = c__[ic] * xi + s[ic] * yi;
#line 138 "dlartv.f"
	y[iy] = c__[ic] * yi - s[ic] * xi;
#line 139 "dlartv.f"
	ix += *incx;
#line 140 "dlartv.f"
	iy += *incy;
#line 141 "dlartv.f"
	ic += *incc;
#line 142 "dlartv.f"
/* L10: */
#line 142 "dlartv.f"
    }
#line 143 "dlartv.f"
    return 0;

/*     End of DLARTV */

} /* dlartv_ */

