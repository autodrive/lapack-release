#line 1 "zlarfgp.f"
/* zlarfgp.f -- translated by f2c (version 20100827).
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

#line 1 "zlarfgp.f"
/* Table of constant values */

static doublecomplex c_b5 = {1.,0.};

/* > \brief \b ZLARFGP generates an elementary reflector (Householder matrix) with non-negatibe beta. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARFGP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfgp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfgp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfgp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARFGP( N, ALPHA, X, INCX, TAU ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       COMPLEX*16         ALPHA, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARFGP generates a complex elementary reflector H of order n, such */
/* > that */
/* > */
/* >       H**H * ( alpha ) = ( beta ),   H**H * H = I. */
/* >              (   x   )   (   0  ) */
/* > */
/* > where alpha and beta are scalars, beta is real and non-negative, and */
/* > x is an (n-1)-element complex vector.  H is represented in the form */
/* > */
/* >       H = I - tau * ( 1 ) * ( 1 v**H ) , */
/* >                     ( v ) */
/* > */
/* > where tau is a complex scalar and v is a complex (n-1)-element */
/* > vector. Note that H is not hermitian. */
/* > */
/* > If the elements of x are all zero and alpha is real, then tau = 0 */
/* > and H is taken to be the unit matrix. */
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
/* >          ALPHA is COMPLEX*16 */
/* >          On entry, the value alpha. */
/* >          On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension */
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
/* >          TAU is COMPLEX*16 */
/* >          The value tau. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlarfgp_(integer *n, doublecomplex *alpha, doublecomplex 
	*x, integer *incx, doublecomplex *tau)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), d_sign(doublereal *, doublereal *), z_abs(
	    doublecomplex *);

    /* Local variables */
    static integer j;
    static doublecomplex savealpha;
    static integer knt;
    static doublereal beta, alphi, alphr;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal xnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlapy3_(doublereal 
	    *, doublereal *, doublereal *), dznrm2_(integer *, doublecomplex *
	    , integer *), dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal bignum;
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static doublereal smlnum;


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

#line 144 "zlarfgp.f"
    /* Parameter adjustments */
#line 144 "zlarfgp.f"
    --x;
#line 144 "zlarfgp.f"

#line 144 "zlarfgp.f"
    /* Function Body */
#line 144 "zlarfgp.f"
    if (*n <= 0) {
#line 145 "zlarfgp.f"
	tau->r = 0., tau->i = 0.;
#line 146 "zlarfgp.f"
	return 0;
#line 147 "zlarfgp.f"
    }

#line 149 "zlarfgp.f"
    i__1 = *n - 1;
#line 149 "zlarfgp.f"
    xnorm = dznrm2_(&i__1, &x[1], incx);
#line 150 "zlarfgp.f"
    alphr = alpha->r;
#line 151 "zlarfgp.f"
    alphi = d_imag(alpha);

#line 153 "zlarfgp.f"
    if (xnorm == 0.) {

/*        H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0. */

#line 157 "zlarfgp.f"
	if (alphi == 0.) {
#line 158 "zlarfgp.f"
	    if (alphr >= 0.) {
/*              When TAU.eq.ZERO, the vector is special-cased to be */
/*              all zeros in the application routines.  We do not need */
/*              to clear it. */
#line 162 "zlarfgp.f"
		tau->r = 0., tau->i = 0.;
#line 163 "zlarfgp.f"
	    } else {
/*              However, the application routines rely on explicit */
/*              zero checks when TAU.ne.ZERO, and we must clear X. */
#line 166 "zlarfgp.f"
		tau->r = 2., tau->i = 0.;
#line 167 "zlarfgp.f"
		i__1 = *n - 1;
#line 167 "zlarfgp.f"
		for (j = 1; j <= i__1; ++j) {
#line 168 "zlarfgp.f"
		    i__2 = (j - 1) * *incx + 1;
#line 168 "zlarfgp.f"
		    x[i__2].r = 0., x[i__2].i = 0.;
#line 169 "zlarfgp.f"
		}
#line 170 "zlarfgp.f"
		z__1.r = -alpha->r, z__1.i = -alpha->i;
#line 170 "zlarfgp.f"
		alpha->r = z__1.r, alpha->i = z__1.i;
#line 171 "zlarfgp.f"
	    }
#line 172 "zlarfgp.f"
	} else {
/*           Only "reflecting" the diagonal entry to be real and non-negative. */
#line 174 "zlarfgp.f"
	    xnorm = dlapy2_(&alphr, &alphi);
#line 175 "zlarfgp.f"
	    d__1 = 1. - alphr / xnorm;
#line 175 "zlarfgp.f"
	    d__2 = -alphi / xnorm;
#line 175 "zlarfgp.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 175 "zlarfgp.f"
	    tau->r = z__1.r, tau->i = z__1.i;
#line 176 "zlarfgp.f"
	    i__1 = *n - 1;
#line 176 "zlarfgp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 177 "zlarfgp.f"
		i__2 = (j - 1) * *incx + 1;
#line 177 "zlarfgp.f"
		x[i__2].r = 0., x[i__2].i = 0.;
#line 178 "zlarfgp.f"
	    }
#line 179 "zlarfgp.f"
	    alpha->r = xnorm, alpha->i = 0.;
#line 180 "zlarfgp.f"
	}
#line 181 "zlarfgp.f"
    } else {

/*        general case */

#line 185 "zlarfgp.f"
	d__1 = dlapy3_(&alphr, &alphi, &xnorm);
#line 185 "zlarfgp.f"
	beta = d_sign(&d__1, &alphr);
#line 186 "zlarfgp.f"
	smlnum = dlamch_("S", (ftnlen)1) / dlamch_("E", (ftnlen)1);
#line 187 "zlarfgp.f"
	bignum = 1. / smlnum;

#line 189 "zlarfgp.f"
	knt = 0;
#line 190 "zlarfgp.f"
	if (abs(beta) < smlnum) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

#line 194 "zlarfgp.f"
L10:
#line 195 "zlarfgp.f"
	    ++knt;
#line 196 "zlarfgp.f"
	    i__1 = *n - 1;
#line 196 "zlarfgp.f"
	    zdscal_(&i__1, &bignum, &x[1], incx);
#line 197 "zlarfgp.f"
	    beta *= bignum;
#line 198 "zlarfgp.f"
	    alphi *= bignum;
#line 199 "zlarfgp.f"
	    alphr *= bignum;
#line 200 "zlarfgp.f"
	    if (abs(beta) < smlnum) {
#line 200 "zlarfgp.f"
		goto L10;
#line 200 "zlarfgp.f"
	    }

/*           New BETA is at most 1, at least SMLNUM */

#line 205 "zlarfgp.f"
	    i__1 = *n - 1;
#line 205 "zlarfgp.f"
	    xnorm = dznrm2_(&i__1, &x[1], incx);
#line 206 "zlarfgp.f"
	    z__1.r = alphr, z__1.i = alphi;
#line 206 "zlarfgp.f"
	    alpha->r = z__1.r, alpha->i = z__1.i;
#line 207 "zlarfgp.f"
	    d__1 = dlapy3_(&alphr, &alphi, &xnorm);
#line 207 "zlarfgp.f"
	    beta = d_sign(&d__1, &alphr);
#line 208 "zlarfgp.f"
	}
#line 209 "zlarfgp.f"
	savealpha.r = alpha->r, savealpha.i = alpha->i;
#line 210 "zlarfgp.f"
	z__1.r = alpha->r + beta, z__1.i = alpha->i;
#line 210 "zlarfgp.f"
	alpha->r = z__1.r, alpha->i = z__1.i;
#line 211 "zlarfgp.f"
	if (beta < 0.) {
#line 212 "zlarfgp.f"
	    beta = -beta;
#line 213 "zlarfgp.f"
	    z__2.r = -alpha->r, z__2.i = -alpha->i;
#line 213 "zlarfgp.f"
	    z__1.r = z__2.r / beta, z__1.i = z__2.i / beta;
#line 213 "zlarfgp.f"
	    tau->r = z__1.r, tau->i = z__1.i;
#line 214 "zlarfgp.f"
	} else {
#line 215 "zlarfgp.f"
	    alphr = alphi * (alphi / alpha->r);
#line 216 "zlarfgp.f"
	    alphr += xnorm * (xnorm / alpha->r);
#line 217 "zlarfgp.f"
	    d__1 = alphr / beta;
#line 217 "zlarfgp.f"
	    d__2 = -alphi / beta;
#line 217 "zlarfgp.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 217 "zlarfgp.f"
	    tau->r = z__1.r, tau->i = z__1.i;
#line 218 "zlarfgp.f"
	    d__1 = -alphr;
#line 218 "zlarfgp.f"
	    z__1.r = d__1, z__1.i = alphi;
#line 218 "zlarfgp.f"
	    alpha->r = z__1.r, alpha->i = z__1.i;
#line 219 "zlarfgp.f"
	}
#line 220 "zlarfgp.f"
	zladiv_(&z__1, &c_b5, alpha);
#line 220 "zlarfgp.f"
	alpha->r = z__1.r, alpha->i = z__1.i;

#line 222 "zlarfgp.f"
	if (z_abs(tau) <= smlnum) {

/*           In the case where the computed TAU ends up being a denormalized number, */
/*           it loses relative accuracy. This is a BIG problem. Solution: flush TAU */
/*           to ZERO (or TWO or whatever makes a nonnegative real number for BETA). */

/*           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.) */
/*           (Thanks Pat. Thanks MathWorks.) */

#line 231 "zlarfgp.f"
	    alphr = savealpha.r;
#line 232 "zlarfgp.f"
	    alphi = d_imag(&savealpha);
#line 233 "zlarfgp.f"
	    if (alphi == 0.) {
#line 234 "zlarfgp.f"
		if (alphr >= 0.) {
#line 235 "zlarfgp.f"
		    tau->r = 0., tau->i = 0.;
#line 236 "zlarfgp.f"
		} else {
#line 237 "zlarfgp.f"
		    tau->r = 2., tau->i = 0.;
#line 238 "zlarfgp.f"
		    i__1 = *n - 1;
#line 238 "zlarfgp.f"
		    for (j = 1; j <= i__1; ++j) {
#line 239 "zlarfgp.f"
			i__2 = (j - 1) * *incx + 1;
#line 239 "zlarfgp.f"
			x[i__2].r = 0., x[i__2].i = 0.;
#line 240 "zlarfgp.f"
		    }
#line 241 "zlarfgp.f"
		    z__1.r = -savealpha.r, z__1.i = -savealpha.i;
#line 241 "zlarfgp.f"
		    beta = z__1.r;
#line 242 "zlarfgp.f"
		}
#line 243 "zlarfgp.f"
	    } else {
#line 244 "zlarfgp.f"
		xnorm = dlapy2_(&alphr, &alphi);
#line 245 "zlarfgp.f"
		d__1 = 1. - alphr / xnorm;
#line 245 "zlarfgp.f"
		d__2 = -alphi / xnorm;
#line 245 "zlarfgp.f"
		z__1.r = d__1, z__1.i = d__2;
#line 245 "zlarfgp.f"
		tau->r = z__1.r, tau->i = z__1.i;
#line 246 "zlarfgp.f"
		i__1 = *n - 1;
#line 246 "zlarfgp.f"
		for (j = 1; j <= i__1; ++j) {
#line 247 "zlarfgp.f"
		    i__2 = (j - 1) * *incx + 1;
#line 247 "zlarfgp.f"
		    x[i__2].r = 0., x[i__2].i = 0.;
#line 248 "zlarfgp.f"
		}
#line 249 "zlarfgp.f"
		beta = xnorm;
#line 250 "zlarfgp.f"
	    }

#line 252 "zlarfgp.f"
	} else {

/*           This is the general case. */

#line 256 "zlarfgp.f"
	    i__1 = *n - 1;
#line 256 "zlarfgp.f"
	    zscal_(&i__1, alpha, &x[1], incx);

#line 258 "zlarfgp.f"
	}

/*        If BETA is subnormal, it may lose relative accuracy */

#line 262 "zlarfgp.f"
	i__1 = knt;
#line 262 "zlarfgp.f"
	for (j = 1; j <= i__1; ++j) {
#line 263 "zlarfgp.f"
	    beta *= smlnum;
#line 264 "zlarfgp.f"
/* L20: */
#line 264 "zlarfgp.f"
	}
#line 265 "zlarfgp.f"
	alpha->r = beta, alpha->i = 0.;
#line 266 "zlarfgp.f"
    }

#line 268 "zlarfgp.f"
    return 0;

/*     End of ZLARFGP */

} /* zlarfgp_ */

