#line 1 "dlaed4.f"
/* dlaed4.f -- translated by f2c (version 20100827).
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

#line 1 "dlaed4.f"
/* > \brief \b DLAED4 used by sstedc. Finds a single root of the secular equation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAED4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I, INFO, N */
/*       DOUBLE PRECISION   DLAM, RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), DELTA( * ), Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine computes the I-th updated eigenvalue of a symmetric */
/* > rank-one modification to a diagonal matrix whose elements are */
/* > given in the array d, and that */
/* > */
/* >            D(i) < D(j)  for  i < j */
/* > */
/* > and that RHO > 0.  This is arranged by the calling routine, and is */
/* > no loss in generality.  The rank-one modified system is thus */
/* > */
/* >            diag( D )  +  RHO * Z * Z_transpose. */
/* > */
/* > where we assume the Euclidean norm of Z is 1. */
/* > */
/* > The method consists of approximating the rational functions in the */
/* > secular equation by simpler interpolating rational functions. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The length of all arrays. */
/* > \endverbatim */
/* > */
/* > \param[in] I */
/* > \verbatim */
/* >          I is INTEGER */
/* >         The index of the eigenvalue to be computed.  1 <= I <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         The original eigenvalues.  It is assumed that they are in */
/* >         order, D(I) < D(J)  for I < J. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (N) */
/* >         The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* >          DELTA is DOUBLE PRECISION array, dimension (N) */
/* >         If N .GT. 2, DELTA contains (D(j) - lambda_I) in its  j-th */
/* >         component.  If N = 1, then DELTA(1) = 1. If N = 2, see DLAED5 */
/* >         for detail. The vector DELTA contains the information necessary */
/* >         to construct the eigenvectors by DLAED3 and DLAED9. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >         The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAM */
/* > \verbatim */
/* >          DLAM is DOUBLE PRECISION */
/* >         The computed lambda_I, the I-th updated eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >         = 0:  successful exit */
/* >         > 0:  if INFO = 1, the updating process failed. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  Logical variable ORGATI (origin-at-i?) is used for distinguishing */
/* >  whether D(i) or D(i+1) is treated as the origin. */
/* > */
/* >            ORGATI = .true.    origin at i */
/* >            ORGATI = .false.   origin at i+1 */
/* > */
/* >   Logical variable SWTCH3 (switch-for-3-poles?) is for noting */
/* >   if we are working with THREE poles! */
/* > */
/* >   MAXIT is the maximum number of iterations allowed for each */
/* >   eigenvalue. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ren-Cang Li, Computer Science Division, University of California */
/* >     at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaed4_(integer *n, integer *i__, doublereal *d__, 
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam,
	 integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b, c__;
    static integer j;
    static doublereal w;
    static integer ii;
    static doublereal dw, zz[3];
    static integer ip1;
    static doublereal del, eta, phi, eps, tau, psi;
    static integer iim1, iip1;
    static doublereal dphi, dpsi;
    static integer iter;
    static doublereal temp, prew, temp1, dltlb, dltub, midpt;
    static integer niter;
    static logical swtch;
    extern /* Subroutine */ int dlaed5_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *), dlaed6_(integer *, 
	    logical *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *);
    static logical swtch3;
    extern doublereal dlamch_(char *, ftnlen);
    static logical orgati;
    static doublereal erretm, rhoinv;


/*  -- LAPACK computational routine (version 3.4.2) -- */
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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Since this routine is called in an inner loop, we do no argument */
/*     checking. */

/*     Quick return for N=1 and 2. */

#line 198 "dlaed4.f"
    /* Parameter adjustments */
#line 198 "dlaed4.f"
    --delta;
#line 198 "dlaed4.f"
    --z__;
#line 198 "dlaed4.f"
    --d__;
#line 198 "dlaed4.f"

#line 198 "dlaed4.f"
    /* Function Body */
#line 198 "dlaed4.f"
    *info = 0;
#line 199 "dlaed4.f"
    if (*n == 1) {

/*         Presumably, I=1 upon entry */

#line 203 "dlaed4.f"
	*dlam = d__[1] + *rho * z__[1] * z__[1];
#line 204 "dlaed4.f"
	delta[1] = 1.;
#line 205 "dlaed4.f"
	return 0;
#line 206 "dlaed4.f"
    }
#line 207 "dlaed4.f"
    if (*n == 2) {
#line 208 "dlaed4.f"
	dlaed5_(i__, &d__[1], &z__[1], &delta[1], rho, dlam);
#line 209 "dlaed4.f"
	return 0;
#line 210 "dlaed4.f"
    }

/*     Compute machine epsilon */

#line 214 "dlaed4.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 215 "dlaed4.f"
    rhoinv = 1. / *rho;

/*     The case I = N */

#line 219 "dlaed4.f"
    if (*i__ == *n) {

/*        Initialize some basic variables */

#line 223 "dlaed4.f"
	ii = *n - 1;
#line 224 "dlaed4.f"
	niter = 1;

/*        Calculate initial guess */

#line 228 "dlaed4.f"
	midpt = *rho / 2.;

/*        If ||Z||_2 is not one, then TEMP should be set to */
/*        RHO * ||Z||_2^2 / TWO */

#line 233 "dlaed4.f"
	i__1 = *n;
#line 233 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 234 "dlaed4.f"
	    delta[j] = d__[j] - d__[*i__] - midpt;
#line 235 "dlaed4.f"
/* L10: */
#line 235 "dlaed4.f"
	}

#line 237 "dlaed4.f"
	psi = 0.;
#line 238 "dlaed4.f"
	i__1 = *n - 2;
#line 238 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 239 "dlaed4.f"
	    psi += z__[j] * z__[j] / delta[j];
#line 240 "dlaed4.f"
/* L20: */
#line 240 "dlaed4.f"
	}

#line 242 "dlaed4.f"
	c__ = rhoinv + psi;
#line 243 "dlaed4.f"
	w = c__ + z__[ii] * z__[ii] / delta[ii] + z__[*n] * z__[*n] / delta[*
		n];

#line 246 "dlaed4.f"
	if (w <= 0.) {
#line 247 "dlaed4.f"
	    temp = z__[*n - 1] * z__[*n - 1] / (d__[*n] - d__[*n - 1] + *rho) 
		    + z__[*n] * z__[*n] / *rho;
#line 249 "dlaed4.f"
	    if (c__ <= temp) {
#line 250 "dlaed4.f"
		tau = *rho;
#line 251 "dlaed4.f"
	    } else {
#line 252 "dlaed4.f"
		del = d__[*n] - d__[*n - 1];
#line 253 "dlaed4.f"
		a = -c__ * del + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n]
			;
#line 254 "dlaed4.f"
		b = z__[*n] * z__[*n] * del;
#line 255 "dlaed4.f"
		if (a < 0.) {
#line 256 "dlaed4.f"
		    tau = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 257 "dlaed4.f"
		} else {
#line 258 "dlaed4.f"
		    tau = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 259 "dlaed4.f"
		}
#line 260 "dlaed4.f"
	    }

/*           It can be proved that */
/*               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO */

#line 265 "dlaed4.f"
	    dltlb = midpt;
#line 266 "dlaed4.f"
	    dltub = *rho;
#line 267 "dlaed4.f"
	} else {
#line 268 "dlaed4.f"
	    del = d__[*n] - d__[*n - 1];
#line 269 "dlaed4.f"
	    a = -c__ * del + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
#line 270 "dlaed4.f"
	    b = z__[*n] * z__[*n] * del;
#line 271 "dlaed4.f"
	    if (a < 0.) {
#line 272 "dlaed4.f"
		tau = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 273 "dlaed4.f"
	    } else {
#line 274 "dlaed4.f"
		tau = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 275 "dlaed4.f"
	    }

/*           It can be proved that */
/*               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2 */

#line 280 "dlaed4.f"
	    dltlb = 0.;
#line 281 "dlaed4.f"
	    dltub = midpt;
#line 282 "dlaed4.f"
	}

#line 284 "dlaed4.f"
	i__1 = *n;
#line 284 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 285 "dlaed4.f"
	    delta[j] = d__[j] - d__[*i__] - tau;
#line 286 "dlaed4.f"
/* L30: */
#line 286 "dlaed4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 290 "dlaed4.f"
	dpsi = 0.;
#line 291 "dlaed4.f"
	psi = 0.;
#line 292 "dlaed4.f"
	erretm = 0.;
#line 293 "dlaed4.f"
	i__1 = ii;
#line 293 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 294 "dlaed4.f"
	    temp = z__[j] / delta[j];
#line 295 "dlaed4.f"
	    psi += z__[j] * temp;
#line 296 "dlaed4.f"
	    dpsi += temp * temp;
#line 297 "dlaed4.f"
	    erretm += psi;
#line 298 "dlaed4.f"
/* L40: */
#line 298 "dlaed4.f"
	}
#line 299 "dlaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 303 "dlaed4.f"
	temp = z__[*n] / delta[*n];
#line 304 "dlaed4.f"
	phi = z__[*n] * temp;
#line 305 "dlaed4.f"
	dphi = temp * temp;
#line 306 "dlaed4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi 
		+ dphi);

#line 309 "dlaed4.f"
	w = rhoinv + phi + psi;

/*        Test for convergence */

#line 313 "dlaed4.f"
	if (abs(w) <= eps * erretm) {
#line 314 "dlaed4.f"
	    *dlam = d__[*i__] + tau;
#line 315 "dlaed4.f"
	    goto L250;
#line 316 "dlaed4.f"
	}

#line 318 "dlaed4.f"
	if (w <= 0.) {
#line 319 "dlaed4.f"
	    dltlb = max(dltlb,tau);
#line 320 "dlaed4.f"
	} else {
#line 321 "dlaed4.f"
	    dltub = min(dltub,tau);
#line 322 "dlaed4.f"
	}

/*        Calculate the new step */

#line 326 "dlaed4.f"
	++niter;
#line 327 "dlaed4.f"
	c__ = w - delta[*n - 1] * dpsi - delta[*n] * dphi;
#line 328 "dlaed4.f"
	a = (delta[*n - 1] + delta[*n]) * w - delta[*n - 1] * delta[*n] * (
		dpsi + dphi);
#line 330 "dlaed4.f"
	b = delta[*n - 1] * delta[*n] * w;
#line 331 "dlaed4.f"
	if (c__ < 0.) {
#line 331 "dlaed4.f"
	    c__ = abs(c__);
#line 331 "dlaed4.f"
	}
#line 333 "dlaed4.f"
	if (c__ == 0.) {
/*          ETA = B/A */
/*           ETA = RHO - TAU */
#line 336 "dlaed4.f"
	    eta = dltub - tau;
#line 337 "dlaed4.f"
	} else if (a >= 0.) {
#line 338 "dlaed4.f"
	    eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ 
		    * 2.);
#line 339 "dlaed4.f"
	} else {
#line 340 "dlaed4.f"
	    eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))
		    );
#line 341 "dlaed4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 349 "dlaed4.f"
	if (w * eta > 0.) {
#line 349 "dlaed4.f"
	    eta = -w / (dpsi + dphi);
#line 349 "dlaed4.f"
	}
#line 351 "dlaed4.f"
	temp = tau + eta;
#line 352 "dlaed4.f"
	if (temp > dltub || temp < dltlb) {
#line 353 "dlaed4.f"
	    if (w < 0.) {
#line 354 "dlaed4.f"
		eta = (dltub - tau) / 2.;
#line 355 "dlaed4.f"
	    } else {
#line 356 "dlaed4.f"
		eta = (dltlb - tau) / 2.;
#line 357 "dlaed4.f"
	    }
#line 358 "dlaed4.f"
	}
#line 359 "dlaed4.f"
	i__1 = *n;
#line 359 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 360 "dlaed4.f"
	    delta[j] -= eta;
#line 361 "dlaed4.f"
/* L50: */
#line 361 "dlaed4.f"
	}

#line 363 "dlaed4.f"
	tau += eta;

/*        Evaluate PSI and the derivative DPSI */

#line 367 "dlaed4.f"
	dpsi = 0.;
#line 368 "dlaed4.f"
	psi = 0.;
#line 369 "dlaed4.f"
	erretm = 0.;
#line 370 "dlaed4.f"
	i__1 = ii;
#line 370 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 371 "dlaed4.f"
	    temp = z__[j] / delta[j];
#line 372 "dlaed4.f"
	    psi += z__[j] * temp;
#line 373 "dlaed4.f"
	    dpsi += temp * temp;
#line 374 "dlaed4.f"
	    erretm += psi;
#line 375 "dlaed4.f"
/* L60: */
#line 375 "dlaed4.f"
	}
#line 376 "dlaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 380 "dlaed4.f"
	temp = z__[*n] / delta[*n];
#line 381 "dlaed4.f"
	phi = z__[*n] * temp;
#line 382 "dlaed4.f"
	dphi = temp * temp;
#line 383 "dlaed4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi 
		+ dphi);

#line 386 "dlaed4.f"
	w = rhoinv + phi + psi;

/*        Main loop to update the values of the array   DELTA */

#line 390 "dlaed4.f"
	iter = niter + 1;

#line 392 "dlaed4.f"
	for (niter = iter; niter <= 30; ++niter) {

/*           Test for convergence */

#line 396 "dlaed4.f"
	    if (abs(w) <= eps * erretm) {
#line 397 "dlaed4.f"
		*dlam = d__[*i__] + tau;
#line 398 "dlaed4.f"
		goto L250;
#line 399 "dlaed4.f"
	    }

#line 401 "dlaed4.f"
	    if (w <= 0.) {
#line 402 "dlaed4.f"
		dltlb = max(dltlb,tau);
#line 403 "dlaed4.f"
	    } else {
#line 404 "dlaed4.f"
		dltub = min(dltub,tau);
#line 405 "dlaed4.f"
	    }

/*           Calculate the new step */

#line 409 "dlaed4.f"
	    c__ = w - delta[*n - 1] * dpsi - delta[*n] * dphi;
#line 410 "dlaed4.f"
	    a = (delta[*n - 1] + delta[*n]) * w - delta[*n - 1] * delta[*n] * 
		    (dpsi + dphi);
#line 412 "dlaed4.f"
	    b = delta[*n - 1] * delta[*n] * w;
#line 413 "dlaed4.f"
	    if (a >= 0.) {
#line 414 "dlaed4.f"
		eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 415 "dlaed4.f"
	    } else {
#line 416 "dlaed4.f"
		eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 417 "dlaed4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 425 "dlaed4.f"
	    if (w * eta > 0.) {
#line 425 "dlaed4.f"
		eta = -w / (dpsi + dphi);
#line 425 "dlaed4.f"
	    }
#line 427 "dlaed4.f"
	    temp = tau + eta;
#line 428 "dlaed4.f"
	    if (temp > dltub || temp < dltlb) {
#line 429 "dlaed4.f"
		if (w < 0.) {
#line 430 "dlaed4.f"
		    eta = (dltub - tau) / 2.;
#line 431 "dlaed4.f"
		} else {
#line 432 "dlaed4.f"
		    eta = (dltlb - tau) / 2.;
#line 433 "dlaed4.f"
		}
#line 434 "dlaed4.f"
	    }
#line 435 "dlaed4.f"
	    i__1 = *n;
#line 435 "dlaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 436 "dlaed4.f"
		delta[j] -= eta;
#line 437 "dlaed4.f"
/* L70: */
#line 437 "dlaed4.f"
	    }

#line 439 "dlaed4.f"
	    tau += eta;

/*           Evaluate PSI and the derivative DPSI */

#line 443 "dlaed4.f"
	    dpsi = 0.;
#line 444 "dlaed4.f"
	    psi = 0.;
#line 445 "dlaed4.f"
	    erretm = 0.;
#line 446 "dlaed4.f"
	    i__1 = ii;
#line 446 "dlaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 447 "dlaed4.f"
		temp = z__[j] / delta[j];
#line 448 "dlaed4.f"
		psi += z__[j] * temp;
#line 449 "dlaed4.f"
		dpsi += temp * temp;
#line 450 "dlaed4.f"
		erretm += psi;
#line 451 "dlaed4.f"
/* L80: */
#line 451 "dlaed4.f"
	    }
#line 452 "dlaed4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 456 "dlaed4.f"
	    temp = z__[*n] / delta[*n];
#line 457 "dlaed4.f"
	    phi = z__[*n] * temp;
#line 458 "dlaed4.f"
	    dphi = temp * temp;
#line 459 "dlaed4.f"
	    erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (
		    dpsi + dphi);

#line 462 "dlaed4.f"
	    w = rhoinv + phi + psi;
#line 463 "dlaed4.f"
/* L90: */
#line 463 "dlaed4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 467 "dlaed4.f"
	*info = 1;
#line 468 "dlaed4.f"
	*dlam = d__[*i__] + tau;
#line 469 "dlaed4.f"
	goto L250;

/*        End for the case I = N */

#line 473 "dlaed4.f"
    } else {

/*        The case for I < N */

#line 477 "dlaed4.f"
	niter = 1;
#line 478 "dlaed4.f"
	ip1 = *i__ + 1;

/*        Calculate initial guess */

#line 482 "dlaed4.f"
	del = d__[ip1] - d__[*i__];
#line 483 "dlaed4.f"
	midpt = del / 2.;
#line 484 "dlaed4.f"
	i__1 = *n;
#line 484 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 485 "dlaed4.f"
	    delta[j] = d__[j] - d__[*i__] - midpt;
#line 486 "dlaed4.f"
/* L100: */
#line 486 "dlaed4.f"
	}

#line 488 "dlaed4.f"
	psi = 0.;
#line 489 "dlaed4.f"
	i__1 = *i__ - 1;
#line 489 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 490 "dlaed4.f"
	    psi += z__[j] * z__[j] / delta[j];
#line 491 "dlaed4.f"
/* L110: */
#line 491 "dlaed4.f"
	}

#line 493 "dlaed4.f"
	phi = 0.;
#line 494 "dlaed4.f"
	i__1 = *i__ + 2;
#line 494 "dlaed4.f"
	for (j = *n; j >= i__1; --j) {
#line 495 "dlaed4.f"
	    phi += z__[j] * z__[j] / delta[j];
#line 496 "dlaed4.f"
/* L120: */
#line 496 "dlaed4.f"
	}
#line 497 "dlaed4.f"
	c__ = rhoinv + psi + phi;
#line 498 "dlaed4.f"
	w = c__ + z__[*i__] * z__[*i__] / delta[*i__] + z__[ip1] * z__[ip1] / 
		delta[ip1];

#line 501 "dlaed4.f"
	if (w > 0.) {

/*           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2 */

/*           We choose d(i) as origin. */

#line 507 "dlaed4.f"
	    orgati = TRUE_;
#line 508 "dlaed4.f"
	    a = c__ * del + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
#line 509 "dlaed4.f"
	    b = z__[*i__] * z__[*i__] * del;
#line 510 "dlaed4.f"
	    if (a > 0.) {
#line 511 "dlaed4.f"
		tau = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 512 "dlaed4.f"
	    } else {
#line 513 "dlaed4.f"
		tau = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 514 "dlaed4.f"
	    }
#line 515 "dlaed4.f"
	    dltlb = 0.;
#line 516 "dlaed4.f"
	    dltub = midpt;
#line 517 "dlaed4.f"
	} else {

/*           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1) */

/*           We choose d(i+1) as origin. */

#line 523 "dlaed4.f"
	    orgati = FALSE_;
#line 524 "dlaed4.f"
	    a = c__ * del - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
#line 525 "dlaed4.f"
	    b = z__[ip1] * z__[ip1] * del;
#line 526 "dlaed4.f"
	    if (a < 0.) {
#line 527 "dlaed4.f"
		tau = b * 2. / (a - sqrt((d__1 = a * a + b * 4. * c__, abs(
			d__1))));
#line 528 "dlaed4.f"
	    } else {
#line 529 "dlaed4.f"
		tau = -(a + sqrt((d__1 = a * a + b * 4. * c__, abs(d__1)))) / 
			(c__ * 2.);
#line 530 "dlaed4.f"
	    }
#line 531 "dlaed4.f"
	    dltlb = -midpt;
#line 532 "dlaed4.f"
	    dltub = 0.;
#line 533 "dlaed4.f"
	}

#line 535 "dlaed4.f"
	if (orgati) {
#line 536 "dlaed4.f"
	    i__1 = *n;
#line 536 "dlaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 537 "dlaed4.f"
		delta[j] = d__[j] - d__[*i__] - tau;
#line 538 "dlaed4.f"
/* L130: */
#line 538 "dlaed4.f"
	    }
#line 539 "dlaed4.f"
	} else {
#line 540 "dlaed4.f"
	    i__1 = *n;
#line 540 "dlaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 541 "dlaed4.f"
		delta[j] = d__[j] - d__[ip1] - tau;
#line 542 "dlaed4.f"
/* L140: */
#line 542 "dlaed4.f"
	    }
#line 543 "dlaed4.f"
	}
#line 544 "dlaed4.f"
	if (orgati) {
#line 545 "dlaed4.f"
	    ii = *i__;
#line 546 "dlaed4.f"
	} else {
#line 547 "dlaed4.f"
	    ii = *i__ + 1;
#line 548 "dlaed4.f"
	}
#line 549 "dlaed4.f"
	iim1 = ii - 1;
#line 550 "dlaed4.f"
	iip1 = ii + 1;

/*        Evaluate PSI and the derivative DPSI */

#line 554 "dlaed4.f"
	dpsi = 0.;
#line 555 "dlaed4.f"
	psi = 0.;
#line 556 "dlaed4.f"
	erretm = 0.;
#line 557 "dlaed4.f"
	i__1 = iim1;
#line 557 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 558 "dlaed4.f"
	    temp = z__[j] / delta[j];
#line 559 "dlaed4.f"
	    psi += z__[j] * temp;
#line 560 "dlaed4.f"
	    dpsi += temp * temp;
#line 561 "dlaed4.f"
	    erretm += psi;
#line 562 "dlaed4.f"
/* L150: */
#line 562 "dlaed4.f"
	}
#line 563 "dlaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 567 "dlaed4.f"
	dphi = 0.;
#line 568 "dlaed4.f"
	phi = 0.;
#line 569 "dlaed4.f"
	i__1 = iip1;
#line 569 "dlaed4.f"
	for (j = *n; j >= i__1; --j) {
#line 570 "dlaed4.f"
	    temp = z__[j] / delta[j];
#line 571 "dlaed4.f"
	    phi += z__[j] * temp;
#line 572 "dlaed4.f"
	    dphi += temp * temp;
#line 573 "dlaed4.f"
	    erretm += phi;
#line 574 "dlaed4.f"
/* L160: */
#line 574 "dlaed4.f"
	}

#line 576 "dlaed4.f"
	w = rhoinv + phi + psi;

/*        W is the value of the secular function with */
/*        its ii-th element removed. */

#line 581 "dlaed4.f"
	swtch3 = FALSE_;
#line 582 "dlaed4.f"
	if (orgati) {
#line 583 "dlaed4.f"
	    if (w < 0.) {
#line 583 "dlaed4.f"
		swtch3 = TRUE_;
#line 583 "dlaed4.f"
	    }
#line 585 "dlaed4.f"
	} else {
#line 586 "dlaed4.f"
	    if (w > 0.) {
#line 586 "dlaed4.f"
		swtch3 = TRUE_;
#line 586 "dlaed4.f"
	    }
#line 588 "dlaed4.f"
	}
#line 589 "dlaed4.f"
	if (ii == 1 || ii == *n) {
#line 589 "dlaed4.f"
	    swtch3 = FALSE_;
#line 589 "dlaed4.f"
	}

#line 592 "dlaed4.f"
	temp = z__[ii] / delta[ii];
#line 593 "dlaed4.f"
	dw = dpsi + dphi + temp * temp;
#line 594 "dlaed4.f"
	temp = z__[ii] * temp;
#line 595 "dlaed4.f"
	w += temp;
#line 596 "dlaed4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. + 
		abs(tau) * dw;

/*        Test for convergence */

#line 601 "dlaed4.f"
	if (abs(w) <= eps * erretm) {
#line 602 "dlaed4.f"
	    if (orgati) {
#line 603 "dlaed4.f"
		*dlam = d__[*i__] + tau;
#line 604 "dlaed4.f"
	    } else {
#line 605 "dlaed4.f"
		*dlam = d__[ip1] + tau;
#line 606 "dlaed4.f"
	    }
#line 607 "dlaed4.f"
	    goto L250;
#line 608 "dlaed4.f"
	}

#line 610 "dlaed4.f"
	if (w <= 0.) {
#line 611 "dlaed4.f"
	    dltlb = max(dltlb,tau);
#line 612 "dlaed4.f"
	} else {
#line 613 "dlaed4.f"
	    dltub = min(dltub,tau);
#line 614 "dlaed4.f"
	}

/*        Calculate the new step */

#line 618 "dlaed4.f"
	++niter;
#line 619 "dlaed4.f"
	if (! swtch3) {
#line 620 "dlaed4.f"
	    if (orgati) {
/* Computing 2nd power */
#line 621 "dlaed4.f"
		d__1 = z__[*i__] / delta[*i__];
#line 621 "dlaed4.f"
		c__ = w - delta[ip1] * dw - (d__[*i__] - d__[ip1]) * (d__1 * 
			d__1);
#line 623 "dlaed4.f"
	    } else {
/* Computing 2nd power */
#line 624 "dlaed4.f"
		d__1 = z__[ip1] / delta[ip1];
#line 624 "dlaed4.f"
		c__ = w - delta[*i__] * dw - (d__[ip1] - d__[*i__]) * (d__1 * 
			d__1);
#line 626 "dlaed4.f"
	    }
#line 627 "dlaed4.f"
	    a = (delta[*i__] + delta[ip1]) * w - delta[*i__] * delta[ip1] * 
		    dw;
#line 629 "dlaed4.f"
	    b = delta[*i__] * delta[ip1] * w;
#line 630 "dlaed4.f"
	    if (c__ == 0.) {
#line 631 "dlaed4.f"
		if (a == 0.) {
#line 632 "dlaed4.f"
		    if (orgati) {
#line 633 "dlaed4.f"
			a = z__[*i__] * z__[*i__] + delta[ip1] * delta[ip1] * 
				(dpsi + dphi);
#line 635 "dlaed4.f"
		    } else {
#line 636 "dlaed4.f"
			a = z__[ip1] * z__[ip1] + delta[*i__] * delta[*i__] * 
				(dpsi + dphi);
#line 638 "dlaed4.f"
		    }
#line 639 "dlaed4.f"
		}
#line 640 "dlaed4.f"
		eta = b / a;
#line 641 "dlaed4.f"
	    } else if (a <= 0.) {
#line 642 "dlaed4.f"
		eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 643 "dlaed4.f"
	    } else {
#line 644 "dlaed4.f"
		eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 645 "dlaed4.f"
	    }
#line 646 "dlaed4.f"
	} else {

/*           Interpolation using THREE most relevant poles */

#line 650 "dlaed4.f"
	    temp = rhoinv + psi + phi;
#line 651 "dlaed4.f"
	    if (orgati) {
#line 652 "dlaed4.f"
		temp1 = z__[iim1] / delta[iim1];
#line 653 "dlaed4.f"
		temp1 *= temp1;
#line 654 "dlaed4.f"
		c__ = temp - delta[iip1] * (dpsi + dphi) - (d__[iim1] - d__[
			iip1]) * temp1;
#line 656 "dlaed4.f"
		zz[0] = z__[iim1] * z__[iim1];
#line 657 "dlaed4.f"
		zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + dphi);
#line 659 "dlaed4.f"
	    } else {
#line 660 "dlaed4.f"
		temp1 = z__[iip1] / delta[iip1];
#line 661 "dlaed4.f"
		temp1 *= temp1;
#line 662 "dlaed4.f"
		c__ = temp - delta[iim1] * (dpsi + dphi) - (d__[iip1] - d__[
			iim1]) * temp1;
#line 664 "dlaed4.f"
		zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - temp1));
#line 666 "dlaed4.f"
		zz[2] = z__[iip1] * z__[iip1];
#line 667 "dlaed4.f"
	    }
#line 668 "dlaed4.f"
	    zz[1] = z__[ii] * z__[ii];
#line 669 "dlaed4.f"
	    dlaed6_(&niter, &orgati, &c__, &delta[iim1], zz, &w, &eta, info);
#line 671 "dlaed4.f"
	    if (*info != 0) {
#line 671 "dlaed4.f"
		goto L250;
#line 671 "dlaed4.f"
	    }
#line 673 "dlaed4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 681 "dlaed4.f"
	if (w * eta >= 0.) {
#line 681 "dlaed4.f"
	    eta = -w / dw;
#line 681 "dlaed4.f"
	}
#line 683 "dlaed4.f"
	temp = tau + eta;
#line 684 "dlaed4.f"
	if (temp > dltub || temp < dltlb) {
#line 685 "dlaed4.f"
	    if (w < 0.) {
#line 686 "dlaed4.f"
		eta = (dltub - tau) / 2.;
#line 687 "dlaed4.f"
	    } else {
#line 688 "dlaed4.f"
		eta = (dltlb - tau) / 2.;
#line 689 "dlaed4.f"
	    }
#line 690 "dlaed4.f"
	}

#line 692 "dlaed4.f"
	prew = w;

#line 694 "dlaed4.f"
	i__1 = *n;
#line 694 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 695 "dlaed4.f"
	    delta[j] -= eta;
#line 696 "dlaed4.f"
/* L180: */
#line 696 "dlaed4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 700 "dlaed4.f"
	dpsi = 0.;
#line 701 "dlaed4.f"
	psi = 0.;
#line 702 "dlaed4.f"
	erretm = 0.;
#line 703 "dlaed4.f"
	i__1 = iim1;
#line 703 "dlaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 704 "dlaed4.f"
	    temp = z__[j] / delta[j];
#line 705 "dlaed4.f"
	    psi += z__[j] * temp;
#line 706 "dlaed4.f"
	    dpsi += temp * temp;
#line 707 "dlaed4.f"
	    erretm += psi;
#line 708 "dlaed4.f"
/* L190: */
#line 708 "dlaed4.f"
	}
#line 709 "dlaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 713 "dlaed4.f"
	dphi = 0.;
#line 714 "dlaed4.f"
	phi = 0.;
#line 715 "dlaed4.f"
	i__1 = iip1;
#line 715 "dlaed4.f"
	for (j = *n; j >= i__1; --j) {
#line 716 "dlaed4.f"
	    temp = z__[j] / delta[j];
#line 717 "dlaed4.f"
	    phi += z__[j] * temp;
#line 718 "dlaed4.f"
	    dphi += temp * temp;
#line 719 "dlaed4.f"
	    erretm += phi;
#line 720 "dlaed4.f"
/* L200: */
#line 720 "dlaed4.f"
	}

#line 722 "dlaed4.f"
	temp = z__[ii] / delta[ii];
#line 723 "dlaed4.f"
	dw = dpsi + dphi + temp * temp;
#line 724 "dlaed4.f"
	temp = z__[ii] * temp;
#line 725 "dlaed4.f"
	w = rhoinv + phi + psi + temp;
#line 726 "dlaed4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. + (
		d__1 = tau + eta, abs(d__1)) * dw;

#line 729 "dlaed4.f"
	swtch = FALSE_;
#line 730 "dlaed4.f"
	if (orgati) {
#line 731 "dlaed4.f"
	    if (-w > abs(prew) / 10.) {
#line 731 "dlaed4.f"
		swtch = TRUE_;
#line 731 "dlaed4.f"
	    }
#line 733 "dlaed4.f"
	} else {
#line 734 "dlaed4.f"
	    if (w > abs(prew) / 10.) {
#line 734 "dlaed4.f"
		swtch = TRUE_;
#line 734 "dlaed4.f"
	    }
#line 736 "dlaed4.f"
	}

#line 738 "dlaed4.f"
	tau += eta;

/*        Main loop to update the values of the array   DELTA */

#line 742 "dlaed4.f"
	iter = niter + 1;

#line 744 "dlaed4.f"
	for (niter = iter; niter <= 30; ++niter) {

/*           Test for convergence */

#line 748 "dlaed4.f"
	    if (abs(w) <= eps * erretm) {
#line 749 "dlaed4.f"
		if (orgati) {
#line 750 "dlaed4.f"
		    *dlam = d__[*i__] + tau;
#line 751 "dlaed4.f"
		} else {
#line 752 "dlaed4.f"
		    *dlam = d__[ip1] + tau;
#line 753 "dlaed4.f"
		}
#line 754 "dlaed4.f"
		goto L250;
#line 755 "dlaed4.f"
	    }

#line 757 "dlaed4.f"
	    if (w <= 0.) {
#line 758 "dlaed4.f"
		dltlb = max(dltlb,tau);
#line 759 "dlaed4.f"
	    } else {
#line 760 "dlaed4.f"
		dltub = min(dltub,tau);
#line 761 "dlaed4.f"
	    }

/*           Calculate the new step */

#line 765 "dlaed4.f"
	    if (! swtch3) {
#line 766 "dlaed4.f"
		if (! swtch) {
#line 767 "dlaed4.f"
		    if (orgati) {
/* Computing 2nd power */
#line 768 "dlaed4.f"
			d__1 = z__[*i__] / delta[*i__];
#line 768 "dlaed4.f"
			c__ = w - delta[ip1] * dw - (d__[*i__] - d__[ip1]) * (
				d__1 * d__1);
#line 770 "dlaed4.f"
		    } else {
/* Computing 2nd power */
#line 771 "dlaed4.f"
			d__1 = z__[ip1] / delta[ip1];
#line 771 "dlaed4.f"
			c__ = w - delta[*i__] * dw - (d__[ip1] - d__[*i__]) * 
				(d__1 * d__1);
#line 773 "dlaed4.f"
		    }
#line 774 "dlaed4.f"
		} else {
#line 775 "dlaed4.f"
		    temp = z__[ii] / delta[ii];
#line 776 "dlaed4.f"
		    if (orgati) {
#line 777 "dlaed4.f"
			dpsi += temp * temp;
#line 778 "dlaed4.f"
		    } else {
#line 779 "dlaed4.f"
			dphi += temp * temp;
#line 780 "dlaed4.f"
		    }
#line 781 "dlaed4.f"
		    c__ = w - delta[*i__] * dpsi - delta[ip1] * dphi;
#line 782 "dlaed4.f"
		}
#line 783 "dlaed4.f"
		a = (delta[*i__] + delta[ip1]) * w - delta[*i__] * delta[ip1] 
			* dw;
#line 785 "dlaed4.f"
		b = delta[*i__] * delta[ip1] * w;
#line 786 "dlaed4.f"
		if (c__ == 0.) {
#line 787 "dlaed4.f"
		    if (a == 0.) {
#line 788 "dlaed4.f"
			if (! swtch) {
#line 789 "dlaed4.f"
			    if (orgati) {
#line 790 "dlaed4.f"
				a = z__[*i__] * z__[*i__] + delta[ip1] * 
					delta[ip1] * (dpsi + dphi);
#line 792 "dlaed4.f"
			    } else {
#line 793 "dlaed4.f"
				a = z__[ip1] * z__[ip1] + delta[*i__] * delta[
					*i__] * (dpsi + dphi);
#line 795 "dlaed4.f"
			    }
#line 796 "dlaed4.f"
			} else {
#line 797 "dlaed4.f"
			    a = delta[*i__] * delta[*i__] * dpsi + delta[ip1] 
				    * delta[ip1] * dphi;
#line 799 "dlaed4.f"
			}
#line 800 "dlaed4.f"
		    }
#line 801 "dlaed4.f"
		    eta = b / a;
#line 802 "dlaed4.f"
		} else if (a <= 0.) {
#line 803 "dlaed4.f"
		    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
			     / (c__ * 2.);
#line 804 "dlaed4.f"
		} else {
#line 805 "dlaed4.f"
		    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, 
			    abs(d__1))));
#line 806 "dlaed4.f"
		}
#line 807 "dlaed4.f"
	    } else {

/*              Interpolation using THREE most relevant poles */

#line 811 "dlaed4.f"
		temp = rhoinv + psi + phi;
#line 812 "dlaed4.f"
		if (swtch) {
#line 813 "dlaed4.f"
		    c__ = temp - delta[iim1] * dpsi - delta[iip1] * dphi;
#line 814 "dlaed4.f"
		    zz[0] = delta[iim1] * delta[iim1] * dpsi;
#line 815 "dlaed4.f"
		    zz[2] = delta[iip1] * delta[iip1] * dphi;
#line 816 "dlaed4.f"
		} else {
#line 817 "dlaed4.f"
		    if (orgati) {
#line 818 "dlaed4.f"
			temp1 = z__[iim1] / delta[iim1];
#line 819 "dlaed4.f"
			temp1 *= temp1;
#line 820 "dlaed4.f"
			c__ = temp - delta[iip1] * (dpsi + dphi) - (d__[iim1] 
				- d__[iip1]) * temp1;
#line 822 "dlaed4.f"
			zz[0] = z__[iim1] * z__[iim1];
#line 823 "dlaed4.f"
			zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + 
				dphi);
#line 825 "dlaed4.f"
		    } else {
#line 826 "dlaed4.f"
			temp1 = z__[iip1] / delta[iip1];
#line 827 "dlaed4.f"
			temp1 *= temp1;
#line 828 "dlaed4.f"
			c__ = temp - delta[iim1] * (dpsi + dphi) - (d__[iip1] 
				- d__[iim1]) * temp1;
#line 830 "dlaed4.f"
			zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - 
				temp1));
#line 832 "dlaed4.f"
			zz[2] = z__[iip1] * z__[iip1];
#line 833 "dlaed4.f"
		    }
#line 834 "dlaed4.f"
		}
#line 835 "dlaed4.f"
		dlaed6_(&niter, &orgati, &c__, &delta[iim1], zz, &w, &eta, 
			info);
#line 837 "dlaed4.f"
		if (*info != 0) {
#line 837 "dlaed4.f"
		    goto L250;
#line 837 "dlaed4.f"
		}
#line 839 "dlaed4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 847 "dlaed4.f"
	    if (w * eta >= 0.) {
#line 847 "dlaed4.f"
		eta = -w / dw;
#line 847 "dlaed4.f"
	    }
#line 849 "dlaed4.f"
	    temp = tau + eta;
#line 850 "dlaed4.f"
	    if (temp > dltub || temp < dltlb) {
#line 851 "dlaed4.f"
		if (w < 0.) {
#line 852 "dlaed4.f"
		    eta = (dltub - tau) / 2.;
#line 853 "dlaed4.f"
		} else {
#line 854 "dlaed4.f"
		    eta = (dltlb - tau) / 2.;
#line 855 "dlaed4.f"
		}
#line 856 "dlaed4.f"
	    }

#line 858 "dlaed4.f"
	    i__1 = *n;
#line 858 "dlaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 859 "dlaed4.f"
		delta[j] -= eta;
#line 860 "dlaed4.f"
/* L210: */
#line 860 "dlaed4.f"
	    }

#line 862 "dlaed4.f"
	    tau += eta;
#line 863 "dlaed4.f"
	    prew = w;

/*           Evaluate PSI and the derivative DPSI */

#line 867 "dlaed4.f"
	    dpsi = 0.;
#line 868 "dlaed4.f"
	    psi = 0.;
#line 869 "dlaed4.f"
	    erretm = 0.;
#line 870 "dlaed4.f"
	    i__1 = iim1;
#line 870 "dlaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 871 "dlaed4.f"
		temp = z__[j] / delta[j];
#line 872 "dlaed4.f"
		psi += z__[j] * temp;
#line 873 "dlaed4.f"
		dpsi += temp * temp;
#line 874 "dlaed4.f"
		erretm += psi;
#line 875 "dlaed4.f"
/* L220: */
#line 875 "dlaed4.f"
	    }
#line 876 "dlaed4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 880 "dlaed4.f"
	    dphi = 0.;
#line 881 "dlaed4.f"
	    phi = 0.;
#line 882 "dlaed4.f"
	    i__1 = iip1;
#line 882 "dlaed4.f"
	    for (j = *n; j >= i__1; --j) {
#line 883 "dlaed4.f"
		temp = z__[j] / delta[j];
#line 884 "dlaed4.f"
		phi += z__[j] * temp;
#line 885 "dlaed4.f"
		dphi += temp * temp;
#line 886 "dlaed4.f"
		erretm += phi;
#line 887 "dlaed4.f"
/* L230: */
#line 887 "dlaed4.f"
	    }

#line 889 "dlaed4.f"
	    temp = z__[ii] / delta[ii];
#line 890 "dlaed4.f"
	    dw = dpsi + dphi + temp * temp;
#line 891 "dlaed4.f"
	    temp = z__[ii] * temp;
#line 892 "dlaed4.f"
	    w = rhoinv + phi + psi + temp;
#line 893 "dlaed4.f"
	    erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. 
		    + abs(tau) * dw;
#line 895 "dlaed4.f"
	    if (w * prew > 0. && abs(w) > abs(prew) / 10.) {
#line 895 "dlaed4.f"
		swtch = ! swtch;
#line 895 "dlaed4.f"
	    }

#line 898 "dlaed4.f"
/* L240: */
#line 898 "dlaed4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 902 "dlaed4.f"
	*info = 1;
#line 903 "dlaed4.f"
	if (orgati) {
#line 904 "dlaed4.f"
	    *dlam = d__[*i__] + tau;
#line 905 "dlaed4.f"
	} else {
#line 906 "dlaed4.f"
	    *dlam = d__[ip1] + tau;
#line 907 "dlaed4.f"
	}

#line 909 "dlaed4.f"
    }

#line 911 "dlaed4.f"
L250:

#line 913 "dlaed4.f"
    return 0;

/*     End of DLAED4 */

} /* dlaed4_ */

