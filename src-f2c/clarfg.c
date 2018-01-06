#line 1 "clarfg.f"
/* clarfg.f -- translated by f2c (version 20100827).
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

#line 1 "clarfg.f"
/* Table of constant values */

static doublecomplex c_b5 = {1.,0.};

/* > \brief \b CLARFG generates an elementary reflector (Householder matrix). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARFG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       COMPLEX            ALPHA, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFG generates a complex elementary reflector H of order n, such */
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
/* >          ALPHA is COMPLEX */
/* >          On entry, the value alpha. */
/* >          On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension */
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
/* >          TAU is COMPLEX */
/* >          The value tau. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clarfg_(integer *n, doublecomplex *alpha, doublecomplex *
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
    static doublereal beta;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal alphi, alphr, xnorm;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *), slapy3_(
	    doublereal *, doublereal *, doublereal *);
    extern /* Double Complex */ VOID cladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal safmin, rsafmn;


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

#line 145 "clarfg.f"
    /* Parameter adjustments */
#line 145 "clarfg.f"
    --x;
#line 145 "clarfg.f"

#line 145 "clarfg.f"
    /* Function Body */
#line 145 "clarfg.f"
    if (*n <= 0) {
#line 146 "clarfg.f"
	tau->r = 0., tau->i = 0.;
#line 147 "clarfg.f"
	return 0;
#line 148 "clarfg.f"
    }

#line 150 "clarfg.f"
    i__1 = *n - 1;
#line 150 "clarfg.f"
    xnorm = scnrm2_(&i__1, &x[1], incx);
#line 151 "clarfg.f"
    alphr = alpha->r;
#line 152 "clarfg.f"
    alphi = d_imag(alpha);

#line 154 "clarfg.f"
    if (xnorm == 0. && alphi == 0.) {

/*        H  =  I */

#line 158 "clarfg.f"
	tau->r = 0., tau->i = 0.;
#line 159 "clarfg.f"
    } else {

/*        general case */

#line 163 "clarfg.f"
	d__1 = slapy3_(&alphr, &alphi, &xnorm);
#line 163 "clarfg.f"
	beta = -d_sign(&d__1, &alphr);
#line 164 "clarfg.f"
	safmin = slamch_("S", (ftnlen)1) / slamch_("E", (ftnlen)1);
#line 165 "clarfg.f"
	rsafmn = 1. / safmin;

#line 167 "clarfg.f"
	knt = 0;
#line 168 "clarfg.f"
	if (abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

#line 172 "clarfg.f"
L10:
#line 173 "clarfg.f"
	    ++knt;
#line 174 "clarfg.f"
	    i__1 = *n - 1;
#line 174 "clarfg.f"
	    csscal_(&i__1, &rsafmn, &x[1], incx);
#line 175 "clarfg.f"
	    beta *= rsafmn;
#line 176 "clarfg.f"
	    alphi *= rsafmn;
#line 177 "clarfg.f"
	    alphr *= rsafmn;
#line 178 "clarfg.f"
	    if (abs(beta) < safmin && knt < 20) {
#line 178 "clarfg.f"
		goto L10;
#line 178 "clarfg.f"
	    }

/*           New BETA is at most 1, at least SAFMIN */

#line 183 "clarfg.f"
	    i__1 = *n - 1;
#line 183 "clarfg.f"
	    xnorm = scnrm2_(&i__1, &x[1], incx);
#line 184 "clarfg.f"
	    z__1.r = alphr, z__1.i = alphi;
#line 184 "clarfg.f"
	    alpha->r = z__1.r, alpha->i = z__1.i;
#line 185 "clarfg.f"
	    d__1 = slapy3_(&alphr, &alphi, &xnorm);
#line 185 "clarfg.f"
	    beta = -d_sign(&d__1, &alphr);
#line 186 "clarfg.f"
	}
#line 187 "clarfg.f"
	d__1 = (beta - alphr) / beta;
#line 187 "clarfg.f"
	d__2 = -alphi / beta;
#line 187 "clarfg.f"
	z__1.r = d__1, z__1.i = d__2;
#line 187 "clarfg.f"
	tau->r = z__1.r, tau->i = z__1.i;
#line 188 "clarfg.f"
	z__2.r = alpha->r - beta, z__2.i = alpha->i;
#line 188 "clarfg.f"
	cladiv_(&z__1, &c_b5, &z__2);
#line 188 "clarfg.f"
	alpha->r = z__1.r, alpha->i = z__1.i;
#line 189 "clarfg.f"
	i__1 = *n - 1;
#line 189 "clarfg.f"
	cscal_(&i__1, alpha, &x[1], incx);

/*        If ALPHA is subnormal, it may lose relative accuracy */

#line 193 "clarfg.f"
	i__1 = knt;
#line 193 "clarfg.f"
	for (j = 1; j <= i__1; ++j) {
#line 194 "clarfg.f"
	    beta *= safmin;
#line 195 "clarfg.f"
/* L20: */
#line 195 "clarfg.f"
	}
#line 196 "clarfg.f"
	alpha->r = beta, alpha->i = 0.;
#line 197 "clarfg.f"
    }

#line 199 "clarfg.f"
    return 0;

/*     End of CLARFG */

} /* clarfg_ */

