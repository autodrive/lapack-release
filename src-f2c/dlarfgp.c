#line 1 "dlarfgp.f"
/* dlarfgp.f -- translated by f2c (version 20100827).
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

#line 1 "dlarfgp.f"
/* > \brief \b DLARFGP generates an elementary reflector (Householder matrix) with non-negative beta. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARFGP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfgp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfgp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfgp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARFGP( N, ALPHA, X, INCX, TAU ) */

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
/* > DLARFGP generates a real elementary reflector H of order n, such */
/* > that */
/* > */
/* >       H * ( alpha ) = ( beta ),   H**T * H = I. */
/* >           (   x   )   (   0  ) */
/* > */
/* > where alpha and beta are scalars, beta is non-negative, and x is */
/* > an (n-1)-element real vector.  H is represented in the form */
/* > */
/* >       H = I - tau * ( 1 ) * ( 1 v**T ) , */
/* >                     ( v ) */
/* > */
/* > where tau is a real scalar and v is a real (n-1)-element */
/* > vector. */
/* > */
/* > If the elements of x are all zero, then tau = 0 and H is taken to be */
/* > the unit matrix. */
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

/* > \date November 2015 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlarfgp_(integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *tau)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer j;
    static doublereal savealpha;
    static integer knt;
    static doublereal beta;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal xnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 142 "dlarfgp.f"
    /* Parameter adjustments */
#line 142 "dlarfgp.f"
    --x;
#line 142 "dlarfgp.f"

#line 142 "dlarfgp.f"
    /* Function Body */
#line 142 "dlarfgp.f"
    if (*n <= 0) {
#line 143 "dlarfgp.f"
	*tau = 0.;
#line 144 "dlarfgp.f"
	return 0;
#line 145 "dlarfgp.f"
    }

#line 147 "dlarfgp.f"
    i__1 = *n - 1;
#line 147 "dlarfgp.f"
    xnorm = dnrm2_(&i__1, &x[1], incx);

#line 149 "dlarfgp.f"
    if (xnorm == 0.) {

/*        H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0 */

#line 153 "dlarfgp.f"
	if (*alpha >= 0.) {
/*           When TAU.eq.ZERO, the vector is special-cased to be */
/*           all zeros in the application routines.  We do not need */
/*           to clear it. */
#line 157 "dlarfgp.f"
	    *tau = 0.;
#line 158 "dlarfgp.f"
	} else {
/*           However, the application routines rely on explicit */
/*           zero checks when TAU.ne.ZERO, and we must clear X. */
#line 161 "dlarfgp.f"
	    *tau = 2.;
#line 162 "dlarfgp.f"
	    i__1 = *n - 1;
#line 162 "dlarfgp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 163 "dlarfgp.f"
		x[(j - 1) * *incx + 1] = 0.;
#line 164 "dlarfgp.f"
	    }
#line 165 "dlarfgp.f"
	    *alpha = -(*alpha);
#line 166 "dlarfgp.f"
	}
#line 167 "dlarfgp.f"
    } else {

/*        general case */

#line 171 "dlarfgp.f"
	d__1 = dlapy2_(alpha, &xnorm);
#line 171 "dlarfgp.f"
	beta = d_sign(&d__1, alpha);
#line 172 "dlarfgp.f"
	smlnum = dlamch_("S", (ftnlen)1) / dlamch_("E", (ftnlen)1);
#line 173 "dlarfgp.f"
	knt = 0;
#line 174 "dlarfgp.f"
	if (abs(beta) < smlnum) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

#line 178 "dlarfgp.f"
	    bignum = 1. / smlnum;
#line 179 "dlarfgp.f"
L10:
#line 180 "dlarfgp.f"
	    ++knt;
#line 181 "dlarfgp.f"
	    i__1 = *n - 1;
#line 181 "dlarfgp.f"
	    dscal_(&i__1, &bignum, &x[1], incx);
#line 182 "dlarfgp.f"
	    beta *= bignum;
#line 183 "dlarfgp.f"
	    *alpha *= bignum;
#line 184 "dlarfgp.f"
	    if (abs(beta) < smlnum) {
#line 184 "dlarfgp.f"
		goto L10;
#line 184 "dlarfgp.f"
	    }

/*           New BETA is at most 1, at least SMLNUM */

#line 189 "dlarfgp.f"
	    i__1 = *n - 1;
#line 189 "dlarfgp.f"
	    xnorm = dnrm2_(&i__1, &x[1], incx);
#line 190 "dlarfgp.f"
	    d__1 = dlapy2_(alpha, &xnorm);
#line 190 "dlarfgp.f"
	    beta = d_sign(&d__1, alpha);
#line 191 "dlarfgp.f"
	}
#line 192 "dlarfgp.f"
	savealpha = *alpha;
#line 193 "dlarfgp.f"
	*alpha += beta;
#line 194 "dlarfgp.f"
	if (beta < 0.) {
#line 195 "dlarfgp.f"
	    beta = -beta;
#line 196 "dlarfgp.f"
	    *tau = -(*alpha) / beta;
#line 197 "dlarfgp.f"
	} else {
#line 198 "dlarfgp.f"
	    *alpha = xnorm * (xnorm / *alpha);
#line 199 "dlarfgp.f"
	    *tau = *alpha / beta;
#line 200 "dlarfgp.f"
	    *alpha = -(*alpha);
#line 201 "dlarfgp.f"
	}

#line 203 "dlarfgp.f"
	if (abs(*tau) <= smlnum) {

/*           In the case where the computed TAU ends up being a denormalized number, */
/*           it loses relative accuracy. This is a BIG problem. Solution: flush TAU */
/*           to ZERO. This explains the next IF statement. */

/*           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.) */
/*           (Thanks Pat. Thanks MathWorks.) */

#line 212 "dlarfgp.f"
	    if (savealpha >= 0.) {
#line 213 "dlarfgp.f"
		*tau = 0.;
#line 214 "dlarfgp.f"
	    } else {
#line 215 "dlarfgp.f"
		*tau = 2.;
#line 216 "dlarfgp.f"
		i__1 = *n - 1;
#line 216 "dlarfgp.f"
		for (j = 1; j <= i__1; ++j) {
#line 217 "dlarfgp.f"
		    x[(j - 1) * *incx + 1] = 0.;
#line 218 "dlarfgp.f"
		}
#line 219 "dlarfgp.f"
		beta = -savealpha;
#line 220 "dlarfgp.f"
	    }

#line 222 "dlarfgp.f"
	} else {

/*           This is the general case. */

#line 226 "dlarfgp.f"
	    i__1 = *n - 1;
#line 226 "dlarfgp.f"
	    d__1 = 1. / *alpha;
#line 226 "dlarfgp.f"
	    dscal_(&i__1, &d__1, &x[1], incx);

#line 228 "dlarfgp.f"
	}

/*        If BETA is subnormal, it may lose relative accuracy */

#line 232 "dlarfgp.f"
	i__1 = knt;
#line 232 "dlarfgp.f"
	for (j = 1; j <= i__1; ++j) {
#line 233 "dlarfgp.f"
	    beta *= smlnum;
#line 234 "dlarfgp.f"
/* L20: */
#line 234 "dlarfgp.f"
	}
#line 235 "dlarfgp.f"
	*alpha = beta;
#line 236 "dlarfgp.f"
    }

#line 238 "dlarfgp.f"
    return 0;

/*     End of DLARFGP */

} /* dlarfgp_ */

