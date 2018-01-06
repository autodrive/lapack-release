#line 1 "slaed4.f"
/* slaed4.f -- translated by f2c (version 20100827).
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

#line 1 "slaed4.f"
/* > \brief \b SLAED4 used by sstedc. Finds a single root of the secular equation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAED4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I, INFO, N */
/*       REAL               DLAM, RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), DELTA( * ), Z( * ) */
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
/* >          D is REAL array, dimension (N) */
/* >         The original eigenvalues.  It is assumed that they are in */
/* >         order, D(I) < D(J)  for I < J. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (N) */
/* >         The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* >          DELTA is REAL array, dimension (N) */
/* >         If N .GT. 2, DELTA contains (D(j) - lambda_I) in its  j-th */
/* >         component.  If N = 1, then DELTA(1) = 1. If N = 2, see SLAED5 */
/* >         for detail. The vector DELTA contains the information necessary */
/* >         to construct the eigenvectors by SLAED3 and SLAED9. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is REAL */
/* >         The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAM */
/* > \verbatim */
/* >          DLAM is REAL */
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

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ren-Cang Li, Computer Science Division, University of California */
/* >     at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slaed4_(integer *n, integer *i__, doublereal *d__, 
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
    extern /* Subroutine */ int slaed5_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *), slaed6_(integer *, 
	    logical *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *);
    static logical swtch3;
    extern doublereal slamch_(char *, ftnlen);
    static logical orgati;
    static doublereal erretm, rhoinv;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

#line 198 "slaed4.f"
    /* Parameter adjustments */
#line 198 "slaed4.f"
    --delta;
#line 198 "slaed4.f"
    --z__;
#line 198 "slaed4.f"
    --d__;
#line 198 "slaed4.f"

#line 198 "slaed4.f"
    /* Function Body */
#line 198 "slaed4.f"
    *info = 0;
#line 199 "slaed4.f"
    if (*n == 1) {

/*         Presumably, I=1 upon entry */

#line 203 "slaed4.f"
	*dlam = d__[1] + *rho * z__[1] * z__[1];
#line 204 "slaed4.f"
	delta[1] = 1.;
#line 205 "slaed4.f"
	return 0;
#line 206 "slaed4.f"
    }
#line 207 "slaed4.f"
    if (*n == 2) {
#line 208 "slaed4.f"
	slaed5_(i__, &d__[1], &z__[1], &delta[1], rho, dlam);
#line 209 "slaed4.f"
	return 0;
#line 210 "slaed4.f"
    }

/*     Compute machine epsilon */

#line 214 "slaed4.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 215 "slaed4.f"
    rhoinv = 1. / *rho;

/*     The case I = N */

#line 219 "slaed4.f"
    if (*i__ == *n) {

/*        Initialize some basic variables */

#line 223 "slaed4.f"
	ii = *n - 1;
#line 224 "slaed4.f"
	niter = 1;

/*        Calculate initial guess */

#line 228 "slaed4.f"
	midpt = *rho / 2.;

/*        If ||Z||_2 is not one, then TEMP should be set to */
/*        RHO * ||Z||_2^2 / TWO */

#line 233 "slaed4.f"
	i__1 = *n;
#line 233 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 234 "slaed4.f"
	    delta[j] = d__[j] - d__[*i__] - midpt;
#line 235 "slaed4.f"
/* L10: */
#line 235 "slaed4.f"
	}

#line 237 "slaed4.f"
	psi = 0.;
#line 238 "slaed4.f"
	i__1 = *n - 2;
#line 238 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 239 "slaed4.f"
	    psi += z__[j] * z__[j] / delta[j];
#line 240 "slaed4.f"
/* L20: */
#line 240 "slaed4.f"
	}

#line 242 "slaed4.f"
	c__ = rhoinv + psi;
#line 243 "slaed4.f"
	w = c__ + z__[ii] * z__[ii] / delta[ii] + z__[*n] * z__[*n] / delta[*
		n];

#line 246 "slaed4.f"
	if (w <= 0.) {
#line 247 "slaed4.f"
	    temp = z__[*n - 1] * z__[*n - 1] / (d__[*n] - d__[*n - 1] + *rho) 
		    + z__[*n] * z__[*n] / *rho;
#line 249 "slaed4.f"
	    if (c__ <= temp) {
#line 250 "slaed4.f"
		tau = *rho;
#line 251 "slaed4.f"
	    } else {
#line 252 "slaed4.f"
		del = d__[*n] - d__[*n - 1];
#line 253 "slaed4.f"
		a = -c__ * del + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n]
			;
#line 254 "slaed4.f"
		b = z__[*n] * z__[*n] * del;
#line 255 "slaed4.f"
		if (a < 0.) {
#line 256 "slaed4.f"
		    tau = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 257 "slaed4.f"
		} else {
#line 258 "slaed4.f"
		    tau = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 259 "slaed4.f"
		}
#line 260 "slaed4.f"
	    }

/*           It can be proved that */
/*               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO */

#line 265 "slaed4.f"
	    dltlb = midpt;
#line 266 "slaed4.f"
	    dltub = *rho;
#line 267 "slaed4.f"
	} else {
#line 268 "slaed4.f"
	    del = d__[*n] - d__[*n - 1];
#line 269 "slaed4.f"
	    a = -c__ * del + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
#line 270 "slaed4.f"
	    b = z__[*n] * z__[*n] * del;
#line 271 "slaed4.f"
	    if (a < 0.) {
#line 272 "slaed4.f"
		tau = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 273 "slaed4.f"
	    } else {
#line 274 "slaed4.f"
		tau = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 275 "slaed4.f"
	    }

/*           It can be proved that */
/*               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2 */

#line 280 "slaed4.f"
	    dltlb = 0.;
#line 281 "slaed4.f"
	    dltub = midpt;
#line 282 "slaed4.f"
	}

#line 284 "slaed4.f"
	i__1 = *n;
#line 284 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 285 "slaed4.f"
	    delta[j] = d__[j] - d__[*i__] - tau;
#line 286 "slaed4.f"
/* L30: */
#line 286 "slaed4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 290 "slaed4.f"
	dpsi = 0.;
#line 291 "slaed4.f"
	psi = 0.;
#line 292 "slaed4.f"
	erretm = 0.;
#line 293 "slaed4.f"
	i__1 = ii;
#line 293 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 294 "slaed4.f"
	    temp = z__[j] / delta[j];
#line 295 "slaed4.f"
	    psi += z__[j] * temp;
#line 296 "slaed4.f"
	    dpsi += temp * temp;
#line 297 "slaed4.f"
	    erretm += psi;
#line 298 "slaed4.f"
/* L40: */
#line 298 "slaed4.f"
	}
#line 299 "slaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 303 "slaed4.f"
	temp = z__[*n] / delta[*n];
#line 304 "slaed4.f"
	phi = z__[*n] * temp;
#line 305 "slaed4.f"
	dphi = temp * temp;
#line 306 "slaed4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi 
		+ dphi);

#line 309 "slaed4.f"
	w = rhoinv + phi + psi;

/*        Test for convergence */

#line 313 "slaed4.f"
	if (abs(w) <= eps * erretm) {
#line 314 "slaed4.f"
	    *dlam = d__[*i__] + tau;
#line 315 "slaed4.f"
	    goto L250;
#line 316 "slaed4.f"
	}

#line 318 "slaed4.f"
	if (w <= 0.) {
#line 319 "slaed4.f"
	    dltlb = max(dltlb,tau);
#line 320 "slaed4.f"
	} else {
#line 321 "slaed4.f"
	    dltub = min(dltub,tau);
#line 322 "slaed4.f"
	}

/*        Calculate the new step */

#line 326 "slaed4.f"
	++niter;
#line 327 "slaed4.f"
	c__ = w - delta[*n - 1] * dpsi - delta[*n] * dphi;
#line 328 "slaed4.f"
	a = (delta[*n - 1] + delta[*n]) * w - delta[*n - 1] * delta[*n] * (
		dpsi + dphi);
#line 330 "slaed4.f"
	b = delta[*n - 1] * delta[*n] * w;
#line 331 "slaed4.f"
	if (c__ < 0.) {
#line 331 "slaed4.f"
	    c__ = abs(c__);
#line 331 "slaed4.f"
	}
#line 333 "slaed4.f"
	if (c__ == 0.) {
/*          ETA = B/A */
/*           ETA = RHO - TAU */
#line 336 "slaed4.f"
	    eta = dltub - tau;
#line 337 "slaed4.f"
	} else if (a >= 0.) {
#line 338 "slaed4.f"
	    eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ 
		    * 2.);
#line 339 "slaed4.f"
	} else {
#line 340 "slaed4.f"
	    eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))
		    );
#line 341 "slaed4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 349 "slaed4.f"
	if (w * eta > 0.) {
#line 349 "slaed4.f"
	    eta = -w / (dpsi + dphi);
#line 349 "slaed4.f"
	}
#line 351 "slaed4.f"
	temp = tau + eta;
#line 352 "slaed4.f"
	if (temp > dltub || temp < dltlb) {
#line 353 "slaed4.f"
	    if (w < 0.) {
#line 354 "slaed4.f"
		eta = (dltub - tau) / 2.;
#line 355 "slaed4.f"
	    } else {
#line 356 "slaed4.f"
		eta = (dltlb - tau) / 2.;
#line 357 "slaed4.f"
	    }
#line 358 "slaed4.f"
	}
#line 359 "slaed4.f"
	i__1 = *n;
#line 359 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 360 "slaed4.f"
	    delta[j] -= eta;
#line 361 "slaed4.f"
/* L50: */
#line 361 "slaed4.f"
	}

#line 363 "slaed4.f"
	tau += eta;

/*        Evaluate PSI and the derivative DPSI */

#line 367 "slaed4.f"
	dpsi = 0.;
#line 368 "slaed4.f"
	psi = 0.;
#line 369 "slaed4.f"
	erretm = 0.;
#line 370 "slaed4.f"
	i__1 = ii;
#line 370 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 371 "slaed4.f"
	    temp = z__[j] / delta[j];
#line 372 "slaed4.f"
	    psi += z__[j] * temp;
#line 373 "slaed4.f"
	    dpsi += temp * temp;
#line 374 "slaed4.f"
	    erretm += psi;
#line 375 "slaed4.f"
/* L60: */
#line 375 "slaed4.f"
	}
#line 376 "slaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 380 "slaed4.f"
	temp = z__[*n] / delta[*n];
#line 381 "slaed4.f"
	phi = z__[*n] * temp;
#line 382 "slaed4.f"
	dphi = temp * temp;
#line 383 "slaed4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi 
		+ dphi);

#line 386 "slaed4.f"
	w = rhoinv + phi + psi;

/*        Main loop to update the values of the array   DELTA */

#line 390 "slaed4.f"
	iter = niter + 1;

#line 392 "slaed4.f"
	for (niter = iter; niter <= 30; ++niter) {

/*           Test for convergence */

#line 396 "slaed4.f"
	    if (abs(w) <= eps * erretm) {
#line 397 "slaed4.f"
		*dlam = d__[*i__] + tau;
#line 398 "slaed4.f"
		goto L250;
#line 399 "slaed4.f"
	    }

#line 401 "slaed4.f"
	    if (w <= 0.) {
#line 402 "slaed4.f"
		dltlb = max(dltlb,tau);
#line 403 "slaed4.f"
	    } else {
#line 404 "slaed4.f"
		dltub = min(dltub,tau);
#line 405 "slaed4.f"
	    }

/*           Calculate the new step */

#line 409 "slaed4.f"
	    c__ = w - delta[*n - 1] * dpsi - delta[*n] * dphi;
#line 410 "slaed4.f"
	    a = (delta[*n - 1] + delta[*n]) * w - delta[*n - 1] * delta[*n] * 
		    (dpsi + dphi);
#line 412 "slaed4.f"
	    b = delta[*n - 1] * delta[*n] * w;
#line 413 "slaed4.f"
	    if (a >= 0.) {
#line 414 "slaed4.f"
		eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 415 "slaed4.f"
	    } else {
#line 416 "slaed4.f"
		eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 417 "slaed4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 425 "slaed4.f"
	    if (w * eta > 0.) {
#line 425 "slaed4.f"
		eta = -w / (dpsi + dphi);
#line 425 "slaed4.f"
	    }
#line 427 "slaed4.f"
	    temp = tau + eta;
#line 428 "slaed4.f"
	    if (temp > dltub || temp < dltlb) {
#line 429 "slaed4.f"
		if (w < 0.) {
#line 430 "slaed4.f"
		    eta = (dltub - tau) / 2.;
#line 431 "slaed4.f"
		} else {
#line 432 "slaed4.f"
		    eta = (dltlb - tau) / 2.;
#line 433 "slaed4.f"
		}
#line 434 "slaed4.f"
	    }
#line 435 "slaed4.f"
	    i__1 = *n;
#line 435 "slaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 436 "slaed4.f"
		delta[j] -= eta;
#line 437 "slaed4.f"
/* L70: */
#line 437 "slaed4.f"
	    }

#line 439 "slaed4.f"
	    tau += eta;

/*           Evaluate PSI and the derivative DPSI */

#line 443 "slaed4.f"
	    dpsi = 0.;
#line 444 "slaed4.f"
	    psi = 0.;
#line 445 "slaed4.f"
	    erretm = 0.;
#line 446 "slaed4.f"
	    i__1 = ii;
#line 446 "slaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 447 "slaed4.f"
		temp = z__[j] / delta[j];
#line 448 "slaed4.f"
		psi += z__[j] * temp;
#line 449 "slaed4.f"
		dpsi += temp * temp;
#line 450 "slaed4.f"
		erretm += psi;
#line 451 "slaed4.f"
/* L80: */
#line 451 "slaed4.f"
	    }
#line 452 "slaed4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 456 "slaed4.f"
	    temp = z__[*n] / delta[*n];
#line 457 "slaed4.f"
	    phi = z__[*n] * temp;
#line 458 "slaed4.f"
	    dphi = temp * temp;
#line 459 "slaed4.f"
	    erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (
		    dpsi + dphi);

#line 462 "slaed4.f"
	    w = rhoinv + phi + psi;
#line 463 "slaed4.f"
/* L90: */
#line 463 "slaed4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 467 "slaed4.f"
	*info = 1;
#line 468 "slaed4.f"
	*dlam = d__[*i__] + tau;
#line 469 "slaed4.f"
	goto L250;

/*        End for the case I = N */

#line 473 "slaed4.f"
    } else {

/*        The case for I < N */

#line 477 "slaed4.f"
	niter = 1;
#line 478 "slaed4.f"
	ip1 = *i__ + 1;

/*        Calculate initial guess */

#line 482 "slaed4.f"
	del = d__[ip1] - d__[*i__];
#line 483 "slaed4.f"
	midpt = del / 2.;
#line 484 "slaed4.f"
	i__1 = *n;
#line 484 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 485 "slaed4.f"
	    delta[j] = d__[j] - d__[*i__] - midpt;
#line 486 "slaed4.f"
/* L100: */
#line 486 "slaed4.f"
	}

#line 488 "slaed4.f"
	psi = 0.;
#line 489 "slaed4.f"
	i__1 = *i__ - 1;
#line 489 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 490 "slaed4.f"
	    psi += z__[j] * z__[j] / delta[j];
#line 491 "slaed4.f"
/* L110: */
#line 491 "slaed4.f"
	}

#line 493 "slaed4.f"
	phi = 0.;
#line 494 "slaed4.f"
	i__1 = *i__ + 2;
#line 494 "slaed4.f"
	for (j = *n; j >= i__1; --j) {
#line 495 "slaed4.f"
	    phi += z__[j] * z__[j] / delta[j];
#line 496 "slaed4.f"
/* L120: */
#line 496 "slaed4.f"
	}
#line 497 "slaed4.f"
	c__ = rhoinv + psi + phi;
#line 498 "slaed4.f"
	w = c__ + z__[*i__] * z__[*i__] / delta[*i__] + z__[ip1] * z__[ip1] / 
		delta[ip1];

#line 501 "slaed4.f"
	if (w > 0.) {

/*           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2 */

/*           We choose d(i) as origin. */

#line 507 "slaed4.f"
	    orgati = TRUE_;
#line 508 "slaed4.f"
	    a = c__ * del + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
#line 509 "slaed4.f"
	    b = z__[*i__] * z__[*i__] * del;
#line 510 "slaed4.f"
	    if (a > 0.) {
#line 511 "slaed4.f"
		tau = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 512 "slaed4.f"
	    } else {
#line 513 "slaed4.f"
		tau = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 514 "slaed4.f"
	    }
#line 515 "slaed4.f"
	    dltlb = 0.;
#line 516 "slaed4.f"
	    dltub = midpt;
#line 517 "slaed4.f"
	} else {

/*           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1) */

/*           We choose d(i+1) as origin. */

#line 523 "slaed4.f"
	    orgati = FALSE_;
#line 524 "slaed4.f"
	    a = c__ * del - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
#line 525 "slaed4.f"
	    b = z__[ip1] * z__[ip1] * del;
#line 526 "slaed4.f"
	    if (a < 0.) {
#line 527 "slaed4.f"
		tau = b * 2. / (a - sqrt((d__1 = a * a + b * 4. * c__, abs(
			d__1))));
#line 528 "slaed4.f"
	    } else {
#line 529 "slaed4.f"
		tau = -(a + sqrt((d__1 = a * a + b * 4. * c__, abs(d__1)))) / 
			(c__ * 2.);
#line 530 "slaed4.f"
	    }
#line 531 "slaed4.f"
	    dltlb = -midpt;
#line 532 "slaed4.f"
	    dltub = 0.;
#line 533 "slaed4.f"
	}

#line 535 "slaed4.f"
	if (orgati) {
#line 536 "slaed4.f"
	    i__1 = *n;
#line 536 "slaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 537 "slaed4.f"
		delta[j] = d__[j] - d__[*i__] - tau;
#line 538 "slaed4.f"
/* L130: */
#line 538 "slaed4.f"
	    }
#line 539 "slaed4.f"
	} else {
#line 540 "slaed4.f"
	    i__1 = *n;
#line 540 "slaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 541 "slaed4.f"
		delta[j] = d__[j] - d__[ip1] - tau;
#line 542 "slaed4.f"
/* L140: */
#line 542 "slaed4.f"
	    }
#line 543 "slaed4.f"
	}
#line 544 "slaed4.f"
	if (orgati) {
#line 545 "slaed4.f"
	    ii = *i__;
#line 546 "slaed4.f"
	} else {
#line 547 "slaed4.f"
	    ii = *i__ + 1;
#line 548 "slaed4.f"
	}
#line 549 "slaed4.f"
	iim1 = ii - 1;
#line 550 "slaed4.f"
	iip1 = ii + 1;

/*        Evaluate PSI and the derivative DPSI */

#line 554 "slaed4.f"
	dpsi = 0.;
#line 555 "slaed4.f"
	psi = 0.;
#line 556 "slaed4.f"
	erretm = 0.;
#line 557 "slaed4.f"
	i__1 = iim1;
#line 557 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 558 "slaed4.f"
	    temp = z__[j] / delta[j];
#line 559 "slaed4.f"
	    psi += z__[j] * temp;
#line 560 "slaed4.f"
	    dpsi += temp * temp;
#line 561 "slaed4.f"
	    erretm += psi;
#line 562 "slaed4.f"
/* L150: */
#line 562 "slaed4.f"
	}
#line 563 "slaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 567 "slaed4.f"
	dphi = 0.;
#line 568 "slaed4.f"
	phi = 0.;
#line 569 "slaed4.f"
	i__1 = iip1;
#line 569 "slaed4.f"
	for (j = *n; j >= i__1; --j) {
#line 570 "slaed4.f"
	    temp = z__[j] / delta[j];
#line 571 "slaed4.f"
	    phi += z__[j] * temp;
#line 572 "slaed4.f"
	    dphi += temp * temp;
#line 573 "slaed4.f"
	    erretm += phi;
#line 574 "slaed4.f"
/* L160: */
#line 574 "slaed4.f"
	}

#line 576 "slaed4.f"
	w = rhoinv + phi + psi;

/*        W is the value of the secular function with */
/*        its ii-th element removed. */

#line 581 "slaed4.f"
	swtch3 = FALSE_;
#line 582 "slaed4.f"
	if (orgati) {
#line 583 "slaed4.f"
	    if (w < 0.) {
#line 583 "slaed4.f"
		swtch3 = TRUE_;
#line 583 "slaed4.f"
	    }
#line 585 "slaed4.f"
	} else {
#line 586 "slaed4.f"
	    if (w > 0.) {
#line 586 "slaed4.f"
		swtch3 = TRUE_;
#line 586 "slaed4.f"
	    }
#line 588 "slaed4.f"
	}
#line 589 "slaed4.f"
	if (ii == 1 || ii == *n) {
#line 589 "slaed4.f"
	    swtch3 = FALSE_;
#line 589 "slaed4.f"
	}

#line 592 "slaed4.f"
	temp = z__[ii] / delta[ii];
#line 593 "slaed4.f"
	dw = dpsi + dphi + temp * temp;
#line 594 "slaed4.f"
	temp = z__[ii] * temp;
#line 595 "slaed4.f"
	w += temp;
#line 596 "slaed4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. + 
		abs(tau) * dw;

/*        Test for convergence */

#line 601 "slaed4.f"
	if (abs(w) <= eps * erretm) {
#line 602 "slaed4.f"
	    if (orgati) {
#line 603 "slaed4.f"
		*dlam = d__[*i__] + tau;
#line 604 "slaed4.f"
	    } else {
#line 605 "slaed4.f"
		*dlam = d__[ip1] + tau;
#line 606 "slaed4.f"
	    }
#line 607 "slaed4.f"
	    goto L250;
#line 608 "slaed4.f"
	}

#line 610 "slaed4.f"
	if (w <= 0.) {
#line 611 "slaed4.f"
	    dltlb = max(dltlb,tau);
#line 612 "slaed4.f"
	} else {
#line 613 "slaed4.f"
	    dltub = min(dltub,tau);
#line 614 "slaed4.f"
	}

/*        Calculate the new step */

#line 618 "slaed4.f"
	++niter;
#line 619 "slaed4.f"
	if (! swtch3) {
#line 620 "slaed4.f"
	    if (orgati) {
/* Computing 2nd power */
#line 621 "slaed4.f"
		d__1 = z__[*i__] / delta[*i__];
#line 621 "slaed4.f"
		c__ = w - delta[ip1] * dw - (d__[*i__] - d__[ip1]) * (d__1 * 
			d__1);
#line 623 "slaed4.f"
	    } else {
/* Computing 2nd power */
#line 624 "slaed4.f"
		d__1 = z__[ip1] / delta[ip1];
#line 624 "slaed4.f"
		c__ = w - delta[*i__] * dw - (d__[ip1] - d__[*i__]) * (d__1 * 
			d__1);
#line 626 "slaed4.f"
	    }
#line 627 "slaed4.f"
	    a = (delta[*i__] + delta[ip1]) * w - delta[*i__] * delta[ip1] * 
		    dw;
#line 629 "slaed4.f"
	    b = delta[*i__] * delta[ip1] * w;
#line 630 "slaed4.f"
	    if (c__ == 0.) {
#line 631 "slaed4.f"
		if (a == 0.) {
#line 632 "slaed4.f"
		    if (orgati) {
#line 633 "slaed4.f"
			a = z__[*i__] * z__[*i__] + delta[ip1] * delta[ip1] * 
				(dpsi + dphi);
#line 635 "slaed4.f"
		    } else {
#line 636 "slaed4.f"
			a = z__[ip1] * z__[ip1] + delta[*i__] * delta[*i__] * 
				(dpsi + dphi);
#line 638 "slaed4.f"
		    }
#line 639 "slaed4.f"
		}
#line 640 "slaed4.f"
		eta = b / a;
#line 641 "slaed4.f"
	    } else if (a <= 0.) {
#line 642 "slaed4.f"
		eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 643 "slaed4.f"
	    } else {
#line 644 "slaed4.f"
		eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 645 "slaed4.f"
	    }
#line 646 "slaed4.f"
	} else {

/*           Interpolation using THREE most relevant poles */

#line 650 "slaed4.f"
	    temp = rhoinv + psi + phi;
#line 651 "slaed4.f"
	    if (orgati) {
#line 652 "slaed4.f"
		temp1 = z__[iim1] / delta[iim1];
#line 653 "slaed4.f"
		temp1 *= temp1;
#line 654 "slaed4.f"
		c__ = temp - delta[iip1] * (dpsi + dphi) - (d__[iim1] - d__[
			iip1]) * temp1;
#line 656 "slaed4.f"
		zz[0] = z__[iim1] * z__[iim1];
#line 657 "slaed4.f"
		zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + dphi);
#line 659 "slaed4.f"
	    } else {
#line 660 "slaed4.f"
		temp1 = z__[iip1] / delta[iip1];
#line 661 "slaed4.f"
		temp1 *= temp1;
#line 662 "slaed4.f"
		c__ = temp - delta[iim1] * (dpsi + dphi) - (d__[iip1] - d__[
			iim1]) * temp1;
#line 664 "slaed4.f"
		zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - temp1));
#line 666 "slaed4.f"
		zz[2] = z__[iip1] * z__[iip1];
#line 667 "slaed4.f"
	    }
#line 668 "slaed4.f"
	    zz[1] = z__[ii] * z__[ii];
#line 669 "slaed4.f"
	    slaed6_(&niter, &orgati, &c__, &delta[iim1], zz, &w, &eta, info);
#line 671 "slaed4.f"
	    if (*info != 0) {
#line 671 "slaed4.f"
		goto L250;
#line 671 "slaed4.f"
	    }
#line 673 "slaed4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 681 "slaed4.f"
	if (w * eta >= 0.) {
#line 681 "slaed4.f"
	    eta = -w / dw;
#line 681 "slaed4.f"
	}
#line 683 "slaed4.f"
	temp = tau + eta;
#line 684 "slaed4.f"
	if (temp > dltub || temp < dltlb) {
#line 685 "slaed4.f"
	    if (w < 0.) {
#line 686 "slaed4.f"
		eta = (dltub - tau) / 2.;
#line 687 "slaed4.f"
	    } else {
#line 688 "slaed4.f"
		eta = (dltlb - tau) / 2.;
#line 689 "slaed4.f"
	    }
#line 690 "slaed4.f"
	}

#line 692 "slaed4.f"
	prew = w;

#line 694 "slaed4.f"
	i__1 = *n;
#line 694 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 695 "slaed4.f"
	    delta[j] -= eta;
#line 696 "slaed4.f"
/* L180: */
#line 696 "slaed4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 700 "slaed4.f"
	dpsi = 0.;
#line 701 "slaed4.f"
	psi = 0.;
#line 702 "slaed4.f"
	erretm = 0.;
#line 703 "slaed4.f"
	i__1 = iim1;
#line 703 "slaed4.f"
	for (j = 1; j <= i__1; ++j) {
#line 704 "slaed4.f"
	    temp = z__[j] / delta[j];
#line 705 "slaed4.f"
	    psi += z__[j] * temp;
#line 706 "slaed4.f"
	    dpsi += temp * temp;
#line 707 "slaed4.f"
	    erretm += psi;
#line 708 "slaed4.f"
/* L190: */
#line 708 "slaed4.f"
	}
#line 709 "slaed4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 713 "slaed4.f"
	dphi = 0.;
#line 714 "slaed4.f"
	phi = 0.;
#line 715 "slaed4.f"
	i__1 = iip1;
#line 715 "slaed4.f"
	for (j = *n; j >= i__1; --j) {
#line 716 "slaed4.f"
	    temp = z__[j] / delta[j];
#line 717 "slaed4.f"
	    phi += z__[j] * temp;
#line 718 "slaed4.f"
	    dphi += temp * temp;
#line 719 "slaed4.f"
	    erretm += phi;
#line 720 "slaed4.f"
/* L200: */
#line 720 "slaed4.f"
	}

#line 722 "slaed4.f"
	temp = z__[ii] / delta[ii];
#line 723 "slaed4.f"
	dw = dpsi + dphi + temp * temp;
#line 724 "slaed4.f"
	temp = z__[ii] * temp;
#line 725 "slaed4.f"
	w = rhoinv + phi + psi + temp;
#line 726 "slaed4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. + (
		d__1 = tau + eta, abs(d__1)) * dw;

#line 729 "slaed4.f"
	swtch = FALSE_;
#line 730 "slaed4.f"
	if (orgati) {
#line 731 "slaed4.f"
	    if (-w > abs(prew) / 10.) {
#line 731 "slaed4.f"
		swtch = TRUE_;
#line 731 "slaed4.f"
	    }
#line 733 "slaed4.f"
	} else {
#line 734 "slaed4.f"
	    if (w > abs(prew) / 10.) {
#line 734 "slaed4.f"
		swtch = TRUE_;
#line 734 "slaed4.f"
	    }
#line 736 "slaed4.f"
	}

#line 738 "slaed4.f"
	tau += eta;

/*        Main loop to update the values of the array   DELTA */

#line 742 "slaed4.f"
	iter = niter + 1;

#line 744 "slaed4.f"
	for (niter = iter; niter <= 30; ++niter) {

/*           Test for convergence */

#line 748 "slaed4.f"
	    if (abs(w) <= eps * erretm) {
#line 749 "slaed4.f"
		if (orgati) {
#line 750 "slaed4.f"
		    *dlam = d__[*i__] + tau;
#line 751 "slaed4.f"
		} else {
#line 752 "slaed4.f"
		    *dlam = d__[ip1] + tau;
#line 753 "slaed4.f"
		}
#line 754 "slaed4.f"
		goto L250;
#line 755 "slaed4.f"
	    }

#line 757 "slaed4.f"
	    if (w <= 0.) {
#line 758 "slaed4.f"
		dltlb = max(dltlb,tau);
#line 759 "slaed4.f"
	    } else {
#line 760 "slaed4.f"
		dltub = min(dltub,tau);
#line 761 "slaed4.f"
	    }

/*           Calculate the new step */

#line 765 "slaed4.f"
	    if (! swtch3) {
#line 766 "slaed4.f"
		if (! swtch) {
#line 767 "slaed4.f"
		    if (orgati) {
/* Computing 2nd power */
#line 768 "slaed4.f"
			d__1 = z__[*i__] / delta[*i__];
#line 768 "slaed4.f"
			c__ = w - delta[ip1] * dw - (d__[*i__] - d__[ip1]) * (
				d__1 * d__1);
#line 770 "slaed4.f"
		    } else {
/* Computing 2nd power */
#line 771 "slaed4.f"
			d__1 = z__[ip1] / delta[ip1];
#line 771 "slaed4.f"
			c__ = w - delta[*i__] * dw - (d__[ip1] - d__[*i__]) * 
				(d__1 * d__1);
#line 773 "slaed4.f"
		    }
#line 774 "slaed4.f"
		} else {
#line 775 "slaed4.f"
		    temp = z__[ii] / delta[ii];
#line 776 "slaed4.f"
		    if (orgati) {
#line 777 "slaed4.f"
			dpsi += temp * temp;
#line 778 "slaed4.f"
		    } else {
#line 779 "slaed4.f"
			dphi += temp * temp;
#line 780 "slaed4.f"
		    }
#line 781 "slaed4.f"
		    c__ = w - delta[*i__] * dpsi - delta[ip1] * dphi;
#line 782 "slaed4.f"
		}
#line 783 "slaed4.f"
		a = (delta[*i__] + delta[ip1]) * w - delta[*i__] * delta[ip1] 
			* dw;
#line 785 "slaed4.f"
		b = delta[*i__] * delta[ip1] * w;
#line 786 "slaed4.f"
		if (c__ == 0.) {
#line 787 "slaed4.f"
		    if (a == 0.) {
#line 788 "slaed4.f"
			if (! swtch) {
#line 789 "slaed4.f"
			    if (orgati) {
#line 790 "slaed4.f"
				a = z__[*i__] * z__[*i__] + delta[ip1] * 
					delta[ip1] * (dpsi + dphi);
#line 792 "slaed4.f"
			    } else {
#line 793 "slaed4.f"
				a = z__[ip1] * z__[ip1] + delta[*i__] * delta[
					*i__] * (dpsi + dphi);
#line 795 "slaed4.f"
			    }
#line 796 "slaed4.f"
			} else {
#line 797 "slaed4.f"
			    a = delta[*i__] * delta[*i__] * dpsi + delta[ip1] 
				    * delta[ip1] * dphi;
#line 799 "slaed4.f"
			}
#line 800 "slaed4.f"
		    }
#line 801 "slaed4.f"
		    eta = b / a;
#line 802 "slaed4.f"
		} else if (a <= 0.) {
#line 803 "slaed4.f"
		    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
			     / (c__ * 2.);
#line 804 "slaed4.f"
		} else {
#line 805 "slaed4.f"
		    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, 
			    abs(d__1))));
#line 806 "slaed4.f"
		}
#line 807 "slaed4.f"
	    } else {

/*              Interpolation using THREE most relevant poles */

#line 811 "slaed4.f"
		temp = rhoinv + psi + phi;
#line 812 "slaed4.f"
		if (swtch) {
#line 813 "slaed4.f"
		    c__ = temp - delta[iim1] * dpsi - delta[iip1] * dphi;
#line 814 "slaed4.f"
		    zz[0] = delta[iim1] * delta[iim1] * dpsi;
#line 815 "slaed4.f"
		    zz[2] = delta[iip1] * delta[iip1] * dphi;
#line 816 "slaed4.f"
		} else {
#line 817 "slaed4.f"
		    if (orgati) {
#line 818 "slaed4.f"
			temp1 = z__[iim1] / delta[iim1];
#line 819 "slaed4.f"
			temp1 *= temp1;
#line 820 "slaed4.f"
			c__ = temp - delta[iip1] * (dpsi + dphi) - (d__[iim1] 
				- d__[iip1]) * temp1;
#line 822 "slaed4.f"
			zz[0] = z__[iim1] * z__[iim1];
#line 823 "slaed4.f"
			zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + 
				dphi);
#line 825 "slaed4.f"
		    } else {
#line 826 "slaed4.f"
			temp1 = z__[iip1] / delta[iip1];
#line 827 "slaed4.f"
			temp1 *= temp1;
#line 828 "slaed4.f"
			c__ = temp - delta[iim1] * (dpsi + dphi) - (d__[iip1] 
				- d__[iim1]) * temp1;
#line 830 "slaed4.f"
			zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - 
				temp1));
#line 832 "slaed4.f"
			zz[2] = z__[iip1] * z__[iip1];
#line 833 "slaed4.f"
		    }
#line 834 "slaed4.f"
		}
#line 835 "slaed4.f"
		slaed6_(&niter, &orgati, &c__, &delta[iim1], zz, &w, &eta, 
			info);
#line 837 "slaed4.f"
		if (*info != 0) {
#line 837 "slaed4.f"
		    goto L250;
#line 837 "slaed4.f"
		}
#line 839 "slaed4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 847 "slaed4.f"
	    if (w * eta >= 0.) {
#line 847 "slaed4.f"
		eta = -w / dw;
#line 847 "slaed4.f"
	    }
#line 849 "slaed4.f"
	    temp = tau + eta;
#line 850 "slaed4.f"
	    if (temp > dltub || temp < dltlb) {
#line 851 "slaed4.f"
		if (w < 0.) {
#line 852 "slaed4.f"
		    eta = (dltub - tau) / 2.;
#line 853 "slaed4.f"
		} else {
#line 854 "slaed4.f"
		    eta = (dltlb - tau) / 2.;
#line 855 "slaed4.f"
		}
#line 856 "slaed4.f"
	    }

#line 858 "slaed4.f"
	    i__1 = *n;
#line 858 "slaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 859 "slaed4.f"
		delta[j] -= eta;
#line 860 "slaed4.f"
/* L210: */
#line 860 "slaed4.f"
	    }

#line 862 "slaed4.f"
	    tau += eta;
#line 863 "slaed4.f"
	    prew = w;

/*           Evaluate PSI and the derivative DPSI */

#line 867 "slaed4.f"
	    dpsi = 0.;
#line 868 "slaed4.f"
	    psi = 0.;
#line 869 "slaed4.f"
	    erretm = 0.;
#line 870 "slaed4.f"
	    i__1 = iim1;
#line 870 "slaed4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 871 "slaed4.f"
		temp = z__[j] / delta[j];
#line 872 "slaed4.f"
		psi += z__[j] * temp;
#line 873 "slaed4.f"
		dpsi += temp * temp;
#line 874 "slaed4.f"
		erretm += psi;
#line 875 "slaed4.f"
/* L220: */
#line 875 "slaed4.f"
	    }
#line 876 "slaed4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 880 "slaed4.f"
	    dphi = 0.;
#line 881 "slaed4.f"
	    phi = 0.;
#line 882 "slaed4.f"
	    i__1 = iip1;
#line 882 "slaed4.f"
	    for (j = *n; j >= i__1; --j) {
#line 883 "slaed4.f"
		temp = z__[j] / delta[j];
#line 884 "slaed4.f"
		phi += z__[j] * temp;
#line 885 "slaed4.f"
		dphi += temp * temp;
#line 886 "slaed4.f"
		erretm += phi;
#line 887 "slaed4.f"
/* L230: */
#line 887 "slaed4.f"
	    }

#line 889 "slaed4.f"
	    temp = z__[ii] / delta[ii];
#line 890 "slaed4.f"
	    dw = dpsi + dphi + temp * temp;
#line 891 "slaed4.f"
	    temp = z__[ii] * temp;
#line 892 "slaed4.f"
	    w = rhoinv + phi + psi + temp;
#line 893 "slaed4.f"
	    erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. 
		    + abs(tau) * dw;
#line 895 "slaed4.f"
	    if (w * prew > 0. && abs(w) > abs(prew) / 10.) {
#line 895 "slaed4.f"
		swtch = ! swtch;
#line 895 "slaed4.f"
	    }

#line 898 "slaed4.f"
/* L240: */
#line 898 "slaed4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 902 "slaed4.f"
	*info = 1;
#line 903 "slaed4.f"
	if (orgati) {
#line 904 "slaed4.f"
	    *dlam = d__[*i__] + tau;
#line 905 "slaed4.f"
	} else {
#line 906 "slaed4.f"
	    *dlam = d__[ip1] + tau;
#line 907 "slaed4.f"
	}

#line 909 "slaed4.f"
    }

#line 911 "slaed4.f"
L250:

#line 913 "slaed4.f"
    return 0;

/*     End of SLAED4 */

} /* slaed4_ */

