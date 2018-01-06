#line 1 "dlarfg.f"
/* dlarfg.f -- translated by f2c (version 20100827).
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

#line 1 "dlarfg.f"
/* > \brief \b DLARFG generates an elementary reflector (Householder matrix). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARFG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       DOUBLE PRECISION   ALPHA, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARFG generates a real elementary reflector H of order n, such */
/* > that */
/* > */
/* >       H * ( alpha ) = ( beta ),   H**T * H = I. */
/* >           (   x   )   (   0  ) */
/* > */
/* > where alpha and beta are scalars, and x is an (n-1)-element real */
/* > vector. H is represented in the form */
/* > */
/* >       H = I - tau * ( 1 ) * ( 1 v**T ) , */
/* >                     ( v ) */
/* > */
/* > where tau is a real scalar and v is a real (n-1)-element */
/* > vector. */
/* > */
/* > If the elements of x are all zero, then tau = 0 and H is taken to be */
/* > the unit matrix. */
/* > */
/* > Otherwise  1 <= tau <= 2. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the elementary reflector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION */
/* >          On entry, the value alpha. */
/* >          On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension */
/* >                         (1+(N-2)*abs(INCX)) */
/* >          On entry, the vector x. */
/* >          On exit, it is overwritten with the vector v. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
/* >          The value tau. */
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
/* Subroutine */ int dlarfg_(integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *tau)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer j, knt;
    static doublereal beta;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal xnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    static doublereal safmin, rsafmn;


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 144 "dlarfg.f"
    /* Parameter adjustments */
#line 144 "dlarfg.f"
    --x;
#line 144 "dlarfg.f"

#line 144 "dlarfg.f"
    /* Function Body */
#line 144 "dlarfg.f"
    if (*n <= 1) {
#line 145 "dlarfg.f"
	*tau = 0.;
#line 146 "dlarfg.f"
	return 0;
#line 147 "dlarfg.f"
    }

#line 149 "dlarfg.f"
    i__1 = *n - 1;
#line 149 "dlarfg.f"
    xnorm = dnrm2_(&i__1, &x[1], incx);

#line 151 "dlarfg.f"
    if (xnorm == 0.) {

/*        H  =  I */

#line 155 "dlarfg.f"
	*tau = 0.;
#line 156 "dlarfg.f"
    } else {

/*        general case */

#line 160 "dlarfg.f"
	d__1 = dlapy2_(alpha, &xnorm);
#line 160 "dlarfg.f"
	beta = -d_sign(&d__1, alpha);
#line 161 "dlarfg.f"
	safmin = dlamch_("S", (ftnlen)1) / dlamch_("E", (ftnlen)1);
#line 162 "dlarfg.f"
	knt = 0;
#line 163 "dlarfg.f"
	if (abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

#line 167 "dlarfg.f"
	    rsafmn = 1. / safmin;
#line 168 "dlarfg.f"
L10:
#line 169 "dlarfg.f"
	    ++knt;
#line 170 "dlarfg.f"
	    i__1 = *n - 1;
#line 170 "dlarfg.f"
	    dscal_(&i__1, &rsafmn, &x[1], incx);
#line 171 "dlarfg.f"
	    beta *= rsafmn;
#line 172 "dlarfg.f"
	    *alpha *= rsafmn;
#line 173 "dlarfg.f"
	    if (abs(beta) < safmin) {
#line 173 "dlarfg.f"
		goto L10;
#line 173 "dlarfg.f"
	    }

/*           New BETA is at most 1, at least SAFMIN */

#line 178 "dlarfg.f"
	    i__1 = *n - 1;
#line 178 "dlarfg.f"
	    xnorm = dnrm2_(&i__1, &x[1], incx);
#line 179 "dlarfg.f"
	    d__1 = dlapy2_(alpha, &xnorm);
#line 179 "dlarfg.f"
	    beta = -d_sign(&d__1, alpha);
#line 180 "dlarfg.f"
	}
#line 181 "dlarfg.f"
	*tau = (beta - *alpha) / beta;
#line 182 "dlarfg.f"
	i__1 = *n - 1;
#line 182 "dlarfg.f"
	d__1 = 1. / (*alpha - beta);
#line 182 "dlarfg.f"
	dscal_(&i__1, &d__1, &x[1], incx);

/*        If ALPHA is subnormal, it may lose relative accuracy */

#line 186 "dlarfg.f"
	i__1 = knt;
#line 186 "dlarfg.f"
	for (j = 1; j <= i__1; ++j) {
#line 187 "dlarfg.f"
	    beta *= safmin;
#line 188 "dlarfg.f"
/* L20: */
#line 188 "dlarfg.f"
	}
#line 189 "dlarfg.f"
	*alpha = beta;
#line 190 "dlarfg.f"
    }

#line 192 "dlarfg.f"
    return 0;

/*     End of DLARFG */

} /* dlarfg_ */

