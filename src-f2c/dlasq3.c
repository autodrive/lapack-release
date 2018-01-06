#line 1 "dlasq3.f"
/* dlasq3.f -- translated by f2c (version 20100827).
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

#line 1 "dlasq3.f"
/* > \brief \b DLASQ3 checks for deflation, computes a shift and calls dqds. Used by sbdsqr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, */
/*                          ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, */
/*                          DN2, G, TAU ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            IEEE */
/*       INTEGER            I0, ITER, N0, NDIV, NFAIL, PP */
/*       DOUBLE PRECISION   DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, */
/*      $                   QMAX, SIGMA, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds. */
/* > In case of failure it changes shifts, and tries again until output */
/* > is positive. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] I0 */
/* > \verbatim */
/* >          I0 is INTEGER */
/* >         First index. */
/* > \endverbatim */
/* > */
/* > \param[in,out] N0 */
/* > \verbatim */
/* >          N0 is INTEGER */
/* >         Last index. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( 4*N0 ) */
/* >         Z holds the qd array. */
/* > \endverbatim */
/* > */
/* > \param[in,out] PP */
/* > \verbatim */
/* >          PP is INTEGER */
/* >         PP=0 for ping, PP=1 for pong. */
/* >         PP=2 indicates that flipping was applied to the Z array */
/* >         and that the initial tests for deflation should not be */
/* >         performed. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN */
/* > \verbatim */
/* >          DMIN is DOUBLE PRECISION */
/* >         Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >         Sum of shifts used in current segment. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DESIG */
/* > \verbatim */
/* >          DESIG is DOUBLE PRECISION */
/* >         Lower order part of SIGMA */
/* > \endverbatim */
/* > */
/* > \param[in] QMAX */
/* > \verbatim */
/* >          QMAX is DOUBLE PRECISION */
/* >         Maximum value of q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] NFAIL */
/* > \verbatim */
/* >          NFAIL is INTEGER */
/* >         Increment NFAIL by 1 each time the shift was too big. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ITER */
/* > \verbatim */
/* >          ITER is INTEGER */
/* >         Increment ITER by 1 for each iteration. */
/* > \endverbatim */
/* > */
/* > \param[in,out] NDIV */
/* > \verbatim */
/* >          NDIV is INTEGER */
/* >         Increment NDIV by 1 for each division. */
/* > \endverbatim */
/* > */
/* > \param[in] IEEE */
/* > \verbatim */
/* >          IEEE is LOGICAL */
/* >         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5). */
/* > \endverbatim */
/* > */
/* > \param[in,out] TTYPE */
/* > \verbatim */
/* >          TTYPE is INTEGER */
/* >         Shift type. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DN */
/* > \verbatim */
/* >          DN is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DN1 */
/* > \verbatim */
/* >          DN1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DN2 */
/* > \verbatim */
/* >          DN2 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] G */
/* > \verbatim */
/* >          G is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
/* > */
/* >         These are passed as arguments in order to save their values */
/* >         between calls to DLASQ3. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlasq3_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *dmin__, doublereal *sigma, doublereal *desig,
	 doublereal *qmax, integer *nfail, integer *iter, integer *ndiv, 
	logical *ieee, integer *ttype, doublereal *dmin1, doublereal *dmin2, 
	doublereal *dn, doublereal *dn1, doublereal *dn2, doublereal *g, 
	doublereal *tau)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal s, t;
    static integer j4, nn;
    static doublereal eps, tol;
    static integer n0in, ipn4;
    static doublereal tol2, temp;
    extern /* Subroutine */ int dlasq4_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *), dlasq5_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , doublereal *), dlasq6_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern logical disnan_(doublereal *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Function .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 226 "dlasq3.f"
    /* Parameter adjustments */
#line 226 "dlasq3.f"
    --z__;
#line 226 "dlasq3.f"

#line 226 "dlasq3.f"
    /* Function Body */
#line 226 "dlasq3.f"
    n0in = *n0;
#line 227 "dlasq3.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 228 "dlasq3.f"
    tol = eps * 100.;
/* Computing 2nd power */
#line 229 "dlasq3.f"
    d__1 = tol;
#line 229 "dlasq3.f"
    tol2 = d__1 * d__1;

/*     Check for deflation. */

#line 233 "dlasq3.f"
L10:

#line 235 "dlasq3.f"
    if (*n0 < *i0) {
#line 235 "dlasq3.f"
	return 0;
#line 235 "dlasq3.f"
    }
#line 237 "dlasq3.f"
    if (*n0 == *i0) {
#line 237 "dlasq3.f"
	goto L20;
#line 237 "dlasq3.f"
    }
#line 239 "dlasq3.f"
    nn = (*n0 << 2) + *pp;
#line 240 "dlasq3.f"
    if (*n0 == *i0 + 1) {
#line 240 "dlasq3.f"
	goto L40;
#line 240 "dlasq3.f"
    }

/*     Check whether E(N0-1) is negligible, 1 eigenvalue. */

#line 245 "dlasq3.f"
    if (z__[nn - 5] > tol2 * (*sigma + z__[nn - 3]) && z__[nn - (*pp << 1) - 
	    4] > tol2 * z__[nn - 7]) {
#line 245 "dlasq3.f"
	goto L30;
#line 245 "dlasq3.f"
    }

#line 249 "dlasq3.f"
L20:

#line 251 "dlasq3.f"
    z__[(*n0 << 2) - 3] = z__[(*n0 << 2) + *pp - 3] + *sigma;
#line 252 "dlasq3.f"
    --(*n0);
#line 253 "dlasq3.f"
    goto L10;

/*     Check  whether E(N0-2) is negligible, 2 eigenvalues. */

#line 257 "dlasq3.f"
L30:

#line 259 "dlasq3.f"
    if (z__[nn - 9] > tol2 * *sigma && z__[nn - (*pp << 1) - 8] > tol2 * z__[
	    nn - 11]) {
#line 259 "dlasq3.f"
	goto L50;
#line 259 "dlasq3.f"
    }

#line 263 "dlasq3.f"
L40:

#line 265 "dlasq3.f"
    if (z__[nn - 3] > z__[nn - 7]) {
#line 266 "dlasq3.f"
	s = z__[nn - 3];
#line 267 "dlasq3.f"
	z__[nn - 3] = z__[nn - 7];
#line 268 "dlasq3.f"
	z__[nn - 7] = s;
#line 269 "dlasq3.f"
    }
#line 270 "dlasq3.f"
    t = (z__[nn - 7] - z__[nn - 3] + z__[nn - 5]) * .5;
#line 271 "dlasq3.f"
    if (z__[nn - 5] > z__[nn - 3] * tol2 && t != 0.) {
#line 272 "dlasq3.f"
	s = z__[nn - 3] * (z__[nn - 5] / t);
#line 273 "dlasq3.f"
	if (s <= t) {
#line 274 "dlasq3.f"
	    s = z__[nn - 3] * (z__[nn - 5] / (t * (sqrt(s / t + 1.) + 1.)));
#line 276 "dlasq3.f"
	} else {
#line 277 "dlasq3.f"
	    s = z__[nn - 3] * (z__[nn - 5] / (t + sqrt(t) * sqrt(t + s)));
#line 278 "dlasq3.f"
	}
#line 279 "dlasq3.f"
	t = z__[nn - 7] + (s + z__[nn - 5]);
#line 280 "dlasq3.f"
	z__[nn - 3] *= z__[nn - 7] / t;
#line 281 "dlasq3.f"
	z__[nn - 7] = t;
#line 282 "dlasq3.f"
    }
#line 283 "dlasq3.f"
    z__[(*n0 << 2) - 7] = z__[nn - 7] + *sigma;
#line 284 "dlasq3.f"
    z__[(*n0 << 2) - 3] = z__[nn - 3] + *sigma;
#line 285 "dlasq3.f"
    *n0 += -2;
#line 286 "dlasq3.f"
    goto L10;

#line 288 "dlasq3.f"
L50:
#line 289 "dlasq3.f"
    if (*pp == 2) {
#line 289 "dlasq3.f"
	*pp = 0;
#line 289 "dlasq3.f"
    }

/*     Reverse the qd-array, if warranted. */

#line 294 "dlasq3.f"
    if (*dmin__ <= 0. || *n0 < n0in) {
#line 295 "dlasq3.f"
	if (z__[(*i0 << 2) + *pp - 3] * 1.5 < z__[(*n0 << 2) + *pp - 3]) {
#line 296 "dlasq3.f"
	    ipn4 = *i0 + *n0 << 2;
#line 297 "dlasq3.f"
	    i__1 = *i0 + *n0 - 1 << 1;
#line 297 "dlasq3.f"
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 298 "dlasq3.f"
		temp = z__[j4 - 3];
#line 299 "dlasq3.f"
		z__[j4 - 3] = z__[ipn4 - j4 - 3];
#line 300 "dlasq3.f"
		z__[ipn4 - j4 - 3] = temp;
#line 301 "dlasq3.f"
		temp = z__[j4 - 2];
#line 302 "dlasq3.f"
		z__[j4 - 2] = z__[ipn4 - j4 - 2];
#line 303 "dlasq3.f"
		z__[ipn4 - j4 - 2] = temp;
#line 304 "dlasq3.f"
		temp = z__[j4 - 1];
#line 305 "dlasq3.f"
		z__[j4 - 1] = z__[ipn4 - j4 - 5];
#line 306 "dlasq3.f"
		z__[ipn4 - j4 - 5] = temp;
#line 307 "dlasq3.f"
		temp = z__[j4];
#line 308 "dlasq3.f"
		z__[j4] = z__[ipn4 - j4 - 4];
#line 309 "dlasq3.f"
		z__[ipn4 - j4 - 4] = temp;
#line 310 "dlasq3.f"
/* L60: */
#line 310 "dlasq3.f"
	    }
#line 311 "dlasq3.f"
	    if (*n0 - *i0 <= 4) {
#line 312 "dlasq3.f"
		z__[(*n0 << 2) + *pp - 1] = z__[(*i0 << 2) + *pp - 1];
#line 313 "dlasq3.f"
		z__[(*n0 << 2) - *pp] = z__[(*i0 << 2) - *pp];
#line 314 "dlasq3.f"
	    }
/* Computing MIN */
#line 315 "dlasq3.f"
	    d__1 = *dmin2, d__2 = z__[(*n0 << 2) + *pp - 1];
#line 315 "dlasq3.f"
	    *dmin2 = min(d__1,d__2);
/* Computing MIN */
#line 316 "dlasq3.f"
	    d__1 = z__[(*n0 << 2) + *pp - 1], d__2 = z__[(*i0 << 2) + *pp - 1]
		    , d__1 = min(d__1,d__2), d__2 = z__[(*i0 << 2) + *pp + 3];
#line 316 "dlasq3.f"
	    z__[(*n0 << 2) + *pp - 1] = min(d__1,d__2);
/* Computing MIN */
#line 318 "dlasq3.f"
	    d__1 = z__[(*n0 << 2) - *pp], d__2 = z__[(*i0 << 2) - *pp], d__1 =
		     min(d__1,d__2), d__2 = z__[(*i0 << 2) - *pp + 4];
#line 318 "dlasq3.f"
	    z__[(*n0 << 2) - *pp] = min(d__1,d__2);
/* Computing MAX */
#line 320 "dlasq3.f"
	    d__1 = *qmax, d__2 = z__[(*i0 << 2) + *pp - 3], d__1 = max(d__1,
		    d__2), d__2 = z__[(*i0 << 2) + *pp + 1];
#line 320 "dlasq3.f"
	    *qmax = max(d__1,d__2);
#line 321 "dlasq3.f"
	    *dmin__ = -0.;
#line 322 "dlasq3.f"
	}
#line 323 "dlasq3.f"
    }

/*     Choose a shift. */

#line 327 "dlasq3.f"
    dlasq4_(i0, n0, &z__[1], pp, &n0in, dmin__, dmin1, dmin2, dn, dn1, dn2, 
	    tau, ttype, g);

/*     Call dqds until DMIN > 0. */

#line 332 "dlasq3.f"
L70:

#line 334 "dlasq3.f"
    dlasq5_(i0, n0, &z__[1], pp, tau, sigma, dmin__, dmin1, dmin2, dn, dn1, 
	    dn2, ieee, &eps);

#line 337 "dlasq3.f"
    *ndiv += *n0 - *i0 + 2;
#line 338 "dlasq3.f"
    ++(*iter);

/*     Check status. */

#line 342 "dlasq3.f"
    if (*dmin__ >= 0. && *dmin1 >= 0.) {

/*        Success. */

#line 346 "dlasq3.f"
	goto L90;

#line 348 "dlasq3.f"
    } else if (*dmin__ < 0. && *dmin1 > 0. && z__[(*n0 - 1 << 2) - *pp] < tol 
	    * (*sigma + *dn1) && abs(*dn) < tol * *sigma) {

/*        Convergence hidden by negative DN. */

#line 354 "dlasq3.f"
	z__[(*n0 - 1 << 2) - *pp + 2] = 0.;
#line 355 "dlasq3.f"
	*dmin__ = 0.;
#line 356 "dlasq3.f"
	goto L90;
#line 357 "dlasq3.f"
    } else if (*dmin__ < 0.) {

/*        TAU too big. Select new TAU and try again. */

#line 361 "dlasq3.f"
	++(*nfail);
#line 362 "dlasq3.f"
	if (*ttype < -22) {

/*           Failed twice. Play it safe. */

#line 366 "dlasq3.f"
	    *tau = 0.;
#line 367 "dlasq3.f"
	} else if (*dmin1 > 0.) {

/*           Late failure. Gives excellent shift. */

#line 371 "dlasq3.f"
	    *tau = (*tau + *dmin__) * (1. - eps * 2.);
#line 372 "dlasq3.f"
	    *ttype += -11;
#line 373 "dlasq3.f"
	} else {

/*           Early failure. Divide by 4. */

#line 377 "dlasq3.f"
	    *tau *= .25;
#line 378 "dlasq3.f"
	    *ttype += -12;
#line 379 "dlasq3.f"
	}
#line 380 "dlasq3.f"
	goto L70;
#line 381 "dlasq3.f"
    } else if (disnan_(dmin__)) {

/*        NaN. */

#line 385 "dlasq3.f"
	if (*tau == 0.) {
#line 386 "dlasq3.f"
	    goto L80;
#line 387 "dlasq3.f"
	} else {
#line 388 "dlasq3.f"
	    *tau = 0.;
#line 389 "dlasq3.f"
	    goto L70;
#line 390 "dlasq3.f"
	}
#line 391 "dlasq3.f"
    } else {

/*        Possible underflow. Play it safe. */

#line 395 "dlasq3.f"
	goto L80;
#line 396 "dlasq3.f"
    }

/*     Risk of underflow. */

#line 400 "dlasq3.f"
L80:
#line 401 "dlasq3.f"
    dlasq6_(i0, n0, &z__[1], pp, dmin__, dmin1, dmin2, dn, dn1, dn2);
#line 402 "dlasq3.f"
    *ndiv += *n0 - *i0 + 2;
#line 403 "dlasq3.f"
    ++(*iter);
#line 404 "dlasq3.f"
    *tau = 0.;

#line 406 "dlasq3.f"
L90:
#line 407 "dlasq3.f"
    if (*tau < *sigma) {
#line 408 "dlasq3.f"
	*desig += *tau;
#line 409 "dlasq3.f"
	t = *sigma + *desig;
#line 410 "dlasq3.f"
	*desig -= t - *sigma;
#line 411 "dlasq3.f"
    } else {
#line 412 "dlasq3.f"
	t = *sigma + *tau;
#line 413 "dlasq3.f"
	*desig = *sigma - (t - *tau) + *desig;
#line 414 "dlasq3.f"
    }
#line 415 "dlasq3.f"
    *sigma = t;

#line 417 "dlasq3.f"
    return 0;

/*     End of DLASQ3 */

} /* dlasq3_ */

