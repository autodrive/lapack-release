#line 1 "slargv.f"
/* slargv.f -- translated by f2c (version 20100827).
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

#line 1 "slargv.f"
/* > \brief \b SLARGV generates a vector of plane rotations with real cosines and real sines. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slargv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slargv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slargv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARGV( N, X, INCX, Y, INCY, C, INCC ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCC, INCX, INCY, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARGV generates a vector of real plane rotations, determined by */
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
/* >          X is REAL array, */
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
/* >          Y is REAL array, */
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
/* >          C is REAL array, dimension (1+(N-1)*INCC) */
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

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slargv_(integer *n, doublereal *x, integer *incx, 
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

#line 134 "slargv.f"
    /* Parameter adjustments */
#line 134 "slargv.f"
    --c__;
#line 134 "slargv.f"
    --y;
#line 134 "slargv.f"
    --x;
#line 134 "slargv.f"

#line 134 "slargv.f"
    /* Function Body */
#line 134 "slargv.f"
    ix = 1;
#line 135 "slargv.f"
    iy = 1;
#line 136 "slargv.f"
    ic = 1;
#line 137 "slargv.f"
    i__1 = *n;
#line 137 "slargv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 138 "slargv.f"
	f = x[ix];
#line 139 "slargv.f"
	g = y[iy];
#line 140 "slargv.f"
	if (g == 0.) {
#line 141 "slargv.f"
	    c__[ic] = 1.;
#line 142 "slargv.f"
	} else if (f == 0.) {
#line 143 "slargv.f"
	    c__[ic] = 0.;
#line 144 "slargv.f"
	    y[iy] = 1.;
#line 145 "slargv.f"
	    x[ix] = g;
#line 146 "slargv.f"
	} else if (abs(f) > abs(g)) {
#line 147 "slargv.f"
	    t = g / f;
#line 148 "slargv.f"
	    tt = sqrt(t * t + 1.);
#line 149 "slargv.f"
	    c__[ic] = 1. / tt;
#line 150 "slargv.f"
	    y[iy] = t * c__[ic];
#line 151 "slargv.f"
	    x[ix] = f * tt;
#line 152 "slargv.f"
	} else {
#line 153 "slargv.f"
	    t = f / g;
#line 154 "slargv.f"
	    tt = sqrt(t * t + 1.);
#line 155 "slargv.f"
	    y[iy] = 1. / tt;
#line 156 "slargv.f"
	    c__[ic] = t * y[iy];
#line 157 "slargv.f"
	    x[ix] = g * tt;
#line 158 "slargv.f"
	}
#line 159 "slargv.f"
	ic += *incc;
#line 160 "slargv.f"
	iy += *incy;
#line 161 "slargv.f"
	ix += *incx;
#line 162 "slargv.f"
/* L10: */
#line 162 "slargv.f"
    }
#line 163 "slargv.f"
    return 0;

/*     End of SLARGV */

} /* slargv_ */

