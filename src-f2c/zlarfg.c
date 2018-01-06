#line 1 "zlarfg.f"
/* zlarfg.f -- translated by f2c (version 20100827).
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

#line 1 "zlarfg.f"
/* Table of constant values */

static doublecomplex c_b5 = {1.,0.};

/* > \brief \b ZLARFG generates an elementary reflector (Householder matrix). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARFG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU ) */

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
/* > ZLARFG generates a complex elementary reflector H of order n, such */
/* > that */
/* > */
/* >       H**H * ( alpha ) = ( beta ),   H**H * H = I. */
/* >              (   x   )   (   0  ) */
/* > */
/* > where alpha and beta are scalars, with beta real, and x is an */
/* > (n-1)-element complex vector. H is represented in the form */
/* > */
/* >       H = I - tau * ( 1 ) * ( 1 v**H ) , */
/* >                     ( v ) */
/* > */
/* > where tau is a complex scalar and v is a complex (n-1)-element */
/* > vector. Note that H is not hermitian. */
/* > */
/* > If the elements of x are all zero and alpha is real, then tau = 0 */
/* > and H is taken to be the unit matrix. */
/* > */
/* > Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 . */
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
/* Subroutine */ int zlarfg_(integer *n, doublecomplex *alpha, doublecomplex *
	x, integer *incx, doublecomplex *tau)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer j, knt;
    static doublereal beta, alphi, alphr;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal xnorm;
    extern doublereal dlapy3_(doublereal *, doublereal *, doublereal *), 
	    dznrm2_(integer *, doublecomplex *, integer *), dlamch_(char *, 
	    ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal rsafmn;
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);


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

#line 145 "zlarfg.f"
    /* Parameter adjustments */
#line 145 "zlarfg.f"
    --x;
#line 145 "zlarfg.f"

#line 145 "zlarfg.f"
    /* Function Body */
#line 145 "zlarfg.f"
    if (*n <= 0) {
#line 146 "zlarfg.f"
	tau->r = 0., tau->i = 0.;
#line 147 "zlarfg.f"
	return 0;
#line 148 "zlarfg.f"
    }

#line 150 "zlarfg.f"
    i__1 = *n - 1;
#line 150 "zlarfg.f"
    xnorm = dznrm2_(&i__1, &x[1], incx);
#line 151 "zlarfg.f"
    alphr = alpha->r;
#line 152 "zlarfg.f"
    alphi = d_imag(alpha);

#line 154 "zlarfg.f"
    if (xnorm == 0. && alphi == 0.) {

/*        H  =  I */

#line 158 "zlarfg.f"
	tau->r = 0., tau->i = 0.;
#line 159 "zlarfg.f"
    } else {

/*        general case */

#line 163 "zlarfg.f"
	d__1 = dlapy3_(&alphr, &alphi, &xnorm);
#line 163 "zlarfg.f"
	beta = -d_sign(&d__1, &alphr);
#line 164 "zlarfg.f"
	safmin = dlamch_("S", (ftnlen)1) / dlamch_("E", (ftnlen)1);
#line 165 "zlarfg.f"
	rsafmn = 1. / safmin;

#line 167 "zlarfg.f"
	knt = 0;
#line 168 "zlarfg.f"
	if (abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

#line 172 "zlarfg.f"
L10:
#line 173 "zlarfg.f"
	    ++knt;
#line 174 "zlarfg.f"
	    i__1 = *n - 1;
#line 174 "zlarfg.f"
	    zdscal_(&i__1, &rsafmn, &x[1], incx);
#line 175 "zlarfg.f"
	    beta *= rsafmn;
#line 176 "zlarfg.f"
	    alphi *= rsafmn;
#line 177 "zlarfg.f"
	    alphr *= rsafmn;
#line 178 "zlarfg.f"
	    if (abs(beta) < safmin) {
#line 178 "zlarfg.f"
		goto L10;
#line 178 "zlarfg.f"
	    }

/*           New BETA is at most 1, at least SAFMIN */

#line 183 "zlarfg.f"
	    i__1 = *n - 1;
#line 183 "zlarfg.f"
	    xnorm = dznrm2_(&i__1, &x[1], incx);
#line 184 "zlarfg.f"
	    z__1.r = alphr, z__1.i = alphi;
#line 184 "zlarfg.f"
	    alpha->r = z__1.r, alpha->i = z__1.i;
#line 185 "zlarfg.f"
	    d__1 = dlapy3_(&alphr, &alphi, &xnorm);
#line 185 "zlarfg.f"
	    beta = -d_sign(&d__1, &alphr);
#line 186 "zlarfg.f"
	}
#line 187 "zlarfg.f"
	d__1 = (beta - alphr) / beta;
#line 187 "zlarfg.f"
	d__2 = -alphi / beta;
#line 187 "zlarfg.f"
	z__1.r = d__1, z__1.i = d__2;
#line 187 "zlarfg.f"
	tau->r = z__1.r, tau->i = z__1.i;
#line 188 "zlarfg.f"
	z__2.r = alpha->r - beta, z__2.i = alpha->i;
#line 188 "zlarfg.f"
	zladiv_(&z__1, &c_b5, &z__2);
#line 188 "zlarfg.f"
	alpha->r = z__1.r, alpha->i = z__1.i;
#line 189 "zlarfg.f"
	i__1 = *n - 1;
#line 189 "zlarfg.f"
	zscal_(&i__1, alpha, &x[1], incx);

/*        If ALPHA is subnormal, it may lose relative accuracy */

#line 193 "zlarfg.f"
	i__1 = knt;
#line 193 "zlarfg.f"
	for (j = 1; j <= i__1; ++j) {
#line 194 "zlarfg.f"
	    beta *= safmin;
#line 195 "zlarfg.f"
/* L20: */
#line 195 "zlarfg.f"
	}
#line 196 "zlarfg.f"
	alpha->r = beta, alpha->i = 0.;
#line 197 "zlarfg.f"
    }

#line 199 "zlarfg.f"
    return 0;

/*     End of ZLARFG */

} /* zlarfg_ */

