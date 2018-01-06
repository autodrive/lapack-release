#line 1 "dlargv.f"
/* dlargv.f -- translated by f2c (version 20100827).
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

#line 1 "dlargv.f"
/* > \brief \b DLARGV generates a vector of plane rotations with real cosines and real sines. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlargv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlargv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlargv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCC, INCX, INCY, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARGV generates a vector of real plane rotations, determined by */
/* > elements of the real vectors x and y. For i = 1,2,...,n */
/* > */
/* >    (  c(i)  s(i) ) ( x(i) ) = ( a(i) ) */
/* >    ( -s(i)  c(i) ) ( y(i) ) = (   0  ) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of plane rotations to be generated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, */
/* >                         dimension (1+(N-1)*INCX) */
/* >          On entry, the vector x. */
/* >          On exit, x(i) is overwritten by a(i), for i = 1,...,n. */
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
/* >          On entry, the vector y. */
/* >          On exit, the sines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >          The increment between elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
/* >          The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCC */
/* > \verbatim */
/* >          INCC is INTEGER */
/* >          The increment between elements of C. INCC > 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlargv_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *c__, integer *incc)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal f, g;
    static integer i__;
    static doublereal t;
    static integer ic, ix, iy;
    static doublereal tt;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 134 "dlargv.f"
    /* Parameter adjustments */
#line 134 "dlargv.f"
    --c__;
#line 134 "dlargv.f"
    --y;
#line 134 "dlargv.f"
    --x;
#line 134 "dlargv.f"

#line 134 "dlargv.f"
    /* Function Body */
#line 134 "dlargv.f"
    ix = 1;
#line 135 "dlargv.f"
    iy = 1;
#line 136 "dlargv.f"
    ic = 1;
#line 137 "dlargv.f"
    i__1 = *n;
#line 137 "dlargv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 138 "dlargv.f"
	f = x[ix];
#line 139 "dlargv.f"
	g = y[iy];
#line 140 "dlargv.f"
	if (g == 0.) {
#line 141 "dlargv.f"
	    c__[ic] = 1.;
#line 142 "dlargv.f"
	} else if (f == 0.) {
#line 143 "dlargv.f"
	    c__[ic] = 0.;
#line 144 "dlargv.f"
	    y[iy] = 1.;
#line 145 "dlargv.f"
	    x[ix] = g;
#line 146 "dlargv.f"
	} else if (abs(f) > abs(g)) {
#line 147 "dlargv.f"
	    t = g / f;
#line 148 "dlargv.f"
	    tt = sqrt(t * t + 1.);
#line 149 "dlargv.f"
	    c__[ic] = 1. / tt;
#line 150 "dlargv.f"
	    y[iy] = t * c__[ic];
#line 151 "dlargv.f"
	    x[ix] = f * tt;
#line 152 "dlargv.f"
	} else {
#line 153 "dlargv.f"
	    t = f / g;
#line 154 "dlargv.f"
	    tt = sqrt(t * t + 1.);
#line 155 "dlargv.f"
	    y[iy] = 1. / tt;
#line 156 "dlargv.f"
	    c__[ic] = t * y[iy];
#line 157 "dlargv.f"
	    x[ix] = g * tt;
#line 158 "dlargv.f"
	}
#line 159 "dlargv.f"
	ic += *incc;
#line 160 "dlargv.f"
	iy += *incy;
#line 161 "dlargv.f"
	ix += *incx;
#line 162 "dlargv.f"
/* L10: */
#line 162 "dlargv.f"
    }
#line 163 "dlargv.f"
    return 0;

/*     End of DLARGV */

} /* dlargv_ */

