#line 1 "slarfgp.f"
/* slarfgp.f -- translated by f2c (version 20100827).
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

#line 1 "slarfgp.f"
/* > \brief \b SLARFGP generates an elementary reflector (Householder matrix) with non-negative beta. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARFGP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfgp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfgp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfgp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARFGP( N, ALPHA, X, INCX, TAU ) */

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
/* > SLARFGP generates a real elementary reflector H of order n, such */
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

/* > \date November 2017 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slarfgp_(integer *n, doublereal *alpha, doublereal *x, 
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
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal xnorm;
    extern doublereal slapy2_(doublereal *, doublereal *), slamch_(char *, 
	    ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 142 "slarfgp.f"
    /* Parameter adjustments */
#line 142 "slarfgp.f"
    --x;
#line 142 "slarfgp.f"

#line 142 "slarfgp.f"
    /* Function Body */
#line 142 "slarfgp.f"
    if (*n <= 0) {
#line 143 "slarfgp.f"
	*tau = 0.;
#line 144 "slarfgp.f"
	return 0;
#line 145 "slarfgp.f"
    }

#line 147 "slarfgp.f"
    i__1 = *n - 1;
#line 147 "slarfgp.f"
    xnorm = snrm2_(&i__1, &x[1], incx);

#line 149 "slarfgp.f"
    if (xnorm == 0.) {

/*        H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0. */

#line 153 "slarfgp.f"
	if (*alpha >= 0.) {
/*           When TAU.eq.ZERO, the vector is special-cased to be */
/*           all zeros in the application routines.  We do not need */
/*           to clear it. */
#line 157 "slarfgp.f"
	    *tau = 0.;
#line 158 "slarfgp.f"
	} else {
/*           However, the application routines rely on explicit */
/*           zero checks when TAU.ne.ZERO, and we must clear X. */
#line 161 "slarfgp.f"
	    *tau = 2.;
#line 162 "slarfgp.f"
	    i__1 = *n - 1;
#line 162 "slarfgp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 163 "slarfgp.f"
		x[(j - 1) * *incx + 1] = 0.;
#line 164 "slarfgp.f"
	    }
#line 165 "slarfgp.f"
	    *alpha = -(*alpha);
#line 166 "slarfgp.f"
	}
#line 167 "slarfgp.f"
    } else {

/*        general case */

#line 171 "slarfgp.f"
	d__1 = slapy2_(alpha, &xnorm);
#line 171 "slarfgp.f"
	beta = d_sign(&d__1, alpha);
#line 172 "slarfgp.f"
	smlnum = slamch_("S", (ftnlen)1) / slamch_("E", (ftnlen)1);
#line 173 "slarfgp.f"
	knt = 0;
#line 174 "slarfgp.f"
	if (abs(beta) < smlnum) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

#line 178 "slarfgp.f"
	    bignum = 1. / smlnum;
#line 179 "slarfgp.f"
L10:
#line 180 "slarfgp.f"
	    ++knt;
#line 181 "slarfgp.f"
	    i__1 = *n - 1;
#line 181 "slarfgp.f"
	    sscal_(&i__1, &bignum, &x[1], incx);
#line 182 "slarfgp.f"
	    beta *= bignum;
#line 183 "slarfgp.f"
	    *alpha *= bignum;
#line 184 "slarfgp.f"
	    if (abs(beta) < smlnum && knt < 20) {
#line 184 "slarfgp.f"
		goto L10;
#line 184 "slarfgp.f"
	    }

/*           New BETA is at most 1, at least SMLNUM */

#line 189 "slarfgp.f"
	    i__1 = *n - 1;
#line 189 "slarfgp.f"
	    xnorm = snrm2_(&i__1, &x[1], incx);
#line 190 "slarfgp.f"
	    d__1 = slapy2_(alpha, &xnorm);
#line 190 "slarfgp.f"
	    beta = d_sign(&d__1, alpha);
#line 191 "slarfgp.f"
	}
#line 192 "slarfgp.f"
	savealpha = *alpha;
#line 193 "slarfgp.f"
	*alpha += beta;
#line 194 "slarfgp.f"
	if (beta < 0.) {
#line 195 "slarfgp.f"
	    beta = -beta;
#line 196 "slarfgp.f"
	    *tau = -(*alpha) / beta;
#line 197 "slarfgp.f"
	} else {
#line 198 "slarfgp.f"
	    *alpha = xnorm * (xnorm / *alpha);
#line 199 "slarfgp.f"
	    *tau = *alpha / beta;
#line 200 "slarfgp.f"
	    *alpha = -(*alpha);
#line 201 "slarfgp.f"
	}

#line 203 "slarfgp.f"
	if (abs(*tau) <= smlnum) {

/*           In the case where the computed TAU ends up being a denormalized number, */
/*           it loses relative accuracy. This is a BIG problem. Solution: flush TAU */
/*           to ZERO. This explains the next IF statement. */

/*           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.) */
/*           (Thanks Pat. Thanks MathWorks.) */

#line 212 "slarfgp.f"
	    if (savealpha >= 0.) {
#line 213 "slarfgp.f"
		*tau = 0.;
#line 214 "slarfgp.f"
	    } else {
#line 215 "slarfgp.f"
		*tau = 2.;
#line 216 "slarfgp.f"
		i__1 = *n - 1;
#line 216 "slarfgp.f"
		for (j = 1; j <= i__1; ++j) {
#line 217 "slarfgp.f"
		    x[(j - 1) * *incx + 1] = 0.;
#line 218 "slarfgp.f"
		}
#line 219 "slarfgp.f"
		beta = -savealpha;
#line 220 "slarfgp.f"
	    }

#line 222 "slarfgp.f"
	} else {

/*           This is the general case. */

#line 226 "slarfgp.f"
	    i__1 = *n - 1;
#line 226 "slarfgp.f"
	    d__1 = 1. / *alpha;
#line 226 "slarfgp.f"
	    sscal_(&i__1, &d__1, &x[1], incx);

#line 228 "slarfgp.f"
	}

/*        If BETA is subnormal, it may lose relative accuracy */

#line 232 "slarfgp.f"
	i__1 = knt;
#line 232 "slarfgp.f"
	for (j = 1; j <= i__1; ++j) {
#line 233 "slarfgp.f"
	    beta *= smlnum;
#line 234 "slarfgp.f"
/* L20: */
#line 234 "slarfgp.f"
	}
#line 235 "slarfgp.f"
	*alpha = beta;
#line 236 "slarfgp.f"
    }

#line 238 "slarfgp.f"
    return 0;

/*     End of SLARFGP */

} /* slarfgp_ */

