#line 1 "slarfg.f"
/* slarfg.f -- translated by f2c (version 20100827).
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

#line 1 "slarfg.f"
/* > \brief \b SLARFG generates an elementary reflector (Householder matrix). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARFG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       REAL               ALPHA, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARFG generates a real elementary reflector H of order n, such */
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
/* >          ALPHA is REAL */
/* >          On entry, the value alpha. */
/* >          On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, dimension */
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
/* >          TAU is REAL */
/* >          The value tau. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slarfg_(integer *n, doublereal *alpha, doublereal *x, 
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
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal xnorm;
    extern doublereal slapy2_(doublereal *, doublereal *), slamch_(char *, 
	    ftnlen);
    static doublereal safmin, rsafmn;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 144 "slarfg.f"
    /* Parameter adjustments */
#line 144 "slarfg.f"
    --x;
#line 144 "slarfg.f"

#line 144 "slarfg.f"
    /* Function Body */
#line 144 "slarfg.f"
    if (*n <= 1) {
#line 145 "slarfg.f"
	*tau = 0.;
#line 146 "slarfg.f"
	return 0;
#line 147 "slarfg.f"
    }

#line 149 "slarfg.f"
    i__1 = *n - 1;
#line 149 "slarfg.f"
    xnorm = snrm2_(&i__1, &x[1], incx);

#line 151 "slarfg.f"
    if (xnorm == 0.) {

/*        H  =  I */

#line 155 "slarfg.f"
	*tau = 0.;
#line 156 "slarfg.f"
    } else {

/*        general case */

#line 160 "slarfg.f"
	d__1 = slapy2_(alpha, &xnorm);
#line 160 "slarfg.f"
	beta = -d_sign(&d__1, alpha);
#line 161 "slarfg.f"
	safmin = slamch_("S", (ftnlen)1) / slamch_("E", (ftnlen)1);
#line 162 "slarfg.f"
	knt = 0;
#line 163 "slarfg.f"
	if (abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

#line 167 "slarfg.f"
	    rsafmn = 1. / safmin;
#line 168 "slarfg.f"
L10:
#line 169 "slarfg.f"
	    ++knt;
#line 170 "slarfg.f"
	    i__1 = *n - 1;
#line 170 "slarfg.f"
	    sscal_(&i__1, &rsafmn, &x[1], incx);
#line 171 "slarfg.f"
	    beta *= rsafmn;
#line 172 "slarfg.f"
	    *alpha *= rsafmn;
#line 173 "slarfg.f"
	    if (abs(beta) < safmin) {
#line 173 "slarfg.f"
		goto L10;
#line 173 "slarfg.f"
	    }

/*           New BETA is at most 1, at least SAFMIN */

#line 178 "slarfg.f"
	    i__1 = *n - 1;
#line 178 "slarfg.f"
	    xnorm = snrm2_(&i__1, &x[1], incx);
#line 179 "slarfg.f"
	    d__1 = slapy2_(alpha, &xnorm);
#line 179 "slarfg.f"
	    beta = -d_sign(&d__1, alpha);
#line 180 "slarfg.f"
	}
#line 181 "slarfg.f"
	*tau = (beta - *alpha) / beta;
#line 182 "slarfg.f"
	i__1 = *n - 1;
#line 182 "slarfg.f"
	d__1 = 1. / (*alpha - beta);
#line 182 "slarfg.f"
	sscal_(&i__1, &d__1, &x[1], incx);

/*        If ALPHA is subnormal, it may lose relative accuracy */

#line 186 "slarfg.f"
	i__1 = knt;
#line 186 "slarfg.f"
	for (j = 1; j <= i__1; ++j) {
#line 187 "slarfg.f"
	    beta *= safmin;
#line 188 "slarfg.f"
/* L20: */
#line 188 "slarfg.f"
	}
#line 189 "slarfg.f"
	*alpha = beta;
#line 190 "slarfg.f"
    }

#line 192 "slarfg.f"
    return 0;

/*     End of SLARFG */

} /* slarfg_ */

