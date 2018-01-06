#line 1 "slasd4.f"
/* slasd4.f -- translated by f2c (version 20100827).
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

#line 1 "slasd4.f"
/* > \brief \b SLASD4 computes the square root of the i-th updated eigenvalue of a positive symmetric rank-one
 modification to a positive diagonal matrix. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASD4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I, INFO, N */
/*       REAL               RHO, SIGMA */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), DELTA( * ), WORK( * ), Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine computes the square root of the I-th updated */
/* > eigenvalue of a positive symmetric rank-one modification to */
/* > a positive diagonal matrix whose entries are given as the squares */
/* > of the corresponding entries in the array d, and that */
/* > */
/* >        0 <= D(i) < D(j)  for  i < j */
/* > */
/* > and that RHO > 0. This is arranged by the calling routine, and is */
/* > no loss in generality.  The rank-one modified system is thus */
/* > */
/* >        diag( D ) * diag( D ) +  RHO * Z * Z_transpose. */
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
/* >          D is REAL array, dimension ( N ) */
/* >         The original eigenvalues.  It is assumed that they are in */
/* >         order, 0 <= D(I) < D(J)  for I < J. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension ( N ) */
/* >         The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* >          DELTA is REAL array, dimension ( N ) */
/* >         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th */
/* >         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA */
/* >         contains the information necessary to construct the */
/* >         (singular) eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is REAL */
/* >         The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* >          SIGMA is REAL */
/* >         The computed sigma_I, the I-th updated eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension ( N ) */
/* >         If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th */
/* >         component.  If N = 1, then WORK( 1 ) = 1. */
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
/* >  Logical variable SWTCH3 (switch-for-3-poles?) is for noting */
/* >  if we are working with THREE poles! */
/* > */
/* >  MAXIT is the maximum number of iterations allowed for each */
/* >  eigenvalue. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ren-Cang Li, Computer Science Division, University of California */
/* >     at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slasd4_(integer *n, integer *i__, doublereal *d__, 
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *
	sigma, doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b, c__;
    static integer j;
    static doublereal w, dd[3];
    static integer ii;
    static doublereal dw, zz[3];
    static integer ip1;
    static doublereal sq2, eta, phi, eps, tau, psi;
    static integer iim1, iip1;
    static doublereal tau2, dphi, sglb, dpsi, sgub;
    static integer iter;
    static doublereal temp, prew, temp1, temp2, dtiim, delsq, dtiip;
    static integer niter;
    static doublereal dtisq;
    static logical swtch;
    static doublereal dtnsq;
    extern /* Subroutine */ int slaed6_(integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal delsq2;
    extern /* Subroutine */ int slasd5_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal dtnsq1;
    static logical swtch3;
    extern doublereal slamch_(char *, ftnlen);
    static logical orgati;
    static doublereal erretm, dtipsq, rhoinv;
    static logical geomavg;


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Since this routine is called in an inner loop, we do no argument */
/*     checking. */

/*     Quick return for N=1 and 2. */

#line 207 "slasd4.f"
    /* Parameter adjustments */
#line 207 "slasd4.f"
    --work;
#line 207 "slasd4.f"
    --delta;
#line 207 "slasd4.f"
    --z__;
#line 207 "slasd4.f"
    --d__;
#line 207 "slasd4.f"

#line 207 "slasd4.f"
    /* Function Body */
#line 207 "slasd4.f"
    *info = 0;
#line 208 "slasd4.f"
    if (*n == 1) {

/*        Presumably, I=1 upon entry */

#line 212 "slasd4.f"
	*sigma = sqrt(d__[1] * d__[1] + *rho * z__[1] * z__[1]);
#line 213 "slasd4.f"
	delta[1] = 1.;
#line 214 "slasd4.f"
	work[1] = 1.;
#line 215 "slasd4.f"
	return 0;
#line 216 "slasd4.f"
    }
#line 217 "slasd4.f"
    if (*n == 2) {
#line 218 "slasd4.f"
	slasd5_(i__, &d__[1], &z__[1], &delta[1], rho, sigma, &work[1]);
#line 219 "slasd4.f"
	return 0;
#line 220 "slasd4.f"
    }

/*     Compute machine epsilon */

#line 224 "slasd4.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 225 "slasd4.f"
    rhoinv = 1. / *rho;
#line 226 "slasd4.f"
    tau2 = 0.;

/*     The case I = N */

#line 230 "slasd4.f"
    if (*i__ == *n) {

/*        Initialize some basic variables */

#line 234 "slasd4.f"
	ii = *n - 1;
#line 235 "slasd4.f"
	niter = 1;

/*        Calculate initial guess */

#line 239 "slasd4.f"
	temp = *rho / 2.;

/*        If ||Z||_2 is not one, then TEMP should be set to */
/*        RHO * ||Z||_2^2 / TWO */

#line 244 "slasd4.f"
	temp1 = temp / (d__[*n] + sqrt(d__[*n] * d__[*n] + temp));
#line 245 "slasd4.f"
	i__1 = *n;
#line 245 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 246 "slasd4.f"
	    work[j] = d__[j] + d__[*n] + temp1;
#line 247 "slasd4.f"
	    delta[j] = d__[j] - d__[*n] - temp1;
#line 248 "slasd4.f"
/* L10: */
#line 248 "slasd4.f"
	}

#line 250 "slasd4.f"
	psi = 0.;
#line 251 "slasd4.f"
	i__1 = *n - 2;
#line 251 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 252 "slasd4.f"
	    psi += z__[j] * z__[j] / (delta[j] * work[j]);
#line 253 "slasd4.f"
/* L20: */
#line 253 "slasd4.f"
	}

#line 255 "slasd4.f"
	c__ = rhoinv + psi;
#line 256 "slasd4.f"
	w = c__ + z__[ii] * z__[ii] / (delta[ii] * work[ii]) + z__[*n] * z__[*
		n] / (delta[*n] * work[*n]);

#line 259 "slasd4.f"
	if (w <= 0.) {
#line 260 "slasd4.f"
	    temp1 = sqrt(d__[*n] * d__[*n] + *rho);
#line 261 "slasd4.f"
	    temp = z__[*n - 1] * z__[*n - 1] / ((d__[*n - 1] + temp1) * (d__[*
		    n] - d__[*n - 1] + *rho / (d__[*n] + temp1))) + z__[*n] * 
		    z__[*n] / *rho;

/*           The following TAU2 is to approximate */
/*           SIGMA_n^2 - D( N )*D( N ) */

#line 268 "slasd4.f"
	    if (c__ <= temp) {
#line 269 "slasd4.f"
		tau = *rho;
#line 270 "slasd4.f"
	    } else {
#line 271 "slasd4.f"
		delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
#line 272 "slasd4.f"
		a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*
			n];
#line 273 "slasd4.f"
		b = z__[*n] * z__[*n] * delsq;
#line 274 "slasd4.f"
		if (a < 0.) {
#line 275 "slasd4.f"
		    tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 276 "slasd4.f"
		} else {
#line 277 "slasd4.f"
		    tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 278 "slasd4.f"
		}
#line 279 "slasd4.f"
		tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));
#line 280 "slasd4.f"
	    }

/*           It can be proved that */
/*               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO */

#line 285 "slasd4.f"
	} else {
#line 286 "slasd4.f"
	    delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
#line 287 "slasd4.f"
	    a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
#line 288 "slasd4.f"
	    b = z__[*n] * z__[*n] * delsq;

/*           The following TAU2 is to approximate */
/*           SIGMA_n^2 - D( N )*D( N ) */

#line 293 "slasd4.f"
	    if (a < 0.) {
#line 294 "slasd4.f"
		tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 295 "slasd4.f"
	    } else {
#line 296 "slasd4.f"
		tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 297 "slasd4.f"
	    }
#line 298 "slasd4.f"
	    tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));

/*           It can be proved that */
/*           D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2 */

#line 304 "slasd4.f"
	}

/*        The following TAU is to approximate SIGMA_n - D( N ) */

/*         TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) ) */

#line 310 "slasd4.f"
	*sigma = d__[*n] + tau;
#line 311 "slasd4.f"
	i__1 = *n;
#line 311 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 312 "slasd4.f"
	    delta[j] = d__[j] - d__[*n] - tau;
#line 313 "slasd4.f"
	    work[j] = d__[j] + d__[*n] + tau;
#line 314 "slasd4.f"
/* L30: */
#line 314 "slasd4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 318 "slasd4.f"
	dpsi = 0.;
#line 319 "slasd4.f"
	psi = 0.;
#line 320 "slasd4.f"
	erretm = 0.;
#line 321 "slasd4.f"
	i__1 = ii;
#line 321 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 322 "slasd4.f"
	    temp = z__[j] / (delta[j] * work[j]);
#line 323 "slasd4.f"
	    psi += z__[j] * temp;
#line 324 "slasd4.f"
	    dpsi += temp * temp;
#line 325 "slasd4.f"
	    erretm += psi;
#line 326 "slasd4.f"
/* L40: */
#line 326 "slasd4.f"
	}
#line 327 "slasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 331 "slasd4.f"
	temp = z__[*n] / (delta[*n] * work[*n]);
#line 332 "slasd4.f"
	phi = z__[*n] * temp;
#line 333 "slasd4.f"
	dphi = temp * temp;
#line 334 "slasd4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $          + ABS( TAU2 )*( DPSI+DPHI ) */

#line 337 "slasd4.f"
	w = rhoinv + phi + psi;

/*        Test for convergence */

#line 341 "slasd4.f"
	if (abs(w) <= eps * erretm) {
#line 342 "slasd4.f"
	    goto L240;
#line 343 "slasd4.f"
	}

/*        Calculate the new step */

#line 347 "slasd4.f"
	++niter;
#line 348 "slasd4.f"
	dtnsq1 = work[*n - 1] * delta[*n - 1];
#line 349 "slasd4.f"
	dtnsq = work[*n] * delta[*n];
#line 350 "slasd4.f"
	c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
#line 351 "slasd4.f"
	a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
#line 352 "slasd4.f"
	b = dtnsq * dtnsq1 * w;
#line 353 "slasd4.f"
	if (c__ < 0.) {
#line 353 "slasd4.f"
	    c__ = abs(c__);
#line 353 "slasd4.f"
	}
#line 355 "slasd4.f"
	if (c__ == 0.) {
#line 356 "slasd4.f"
	    eta = *rho - *sigma * *sigma;
#line 357 "slasd4.f"
	} else if (a >= 0.) {
#line 358 "slasd4.f"
	    eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ 
		    * 2.);
#line 359 "slasd4.f"
	} else {
#line 360 "slasd4.f"
	    eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))
		    );
#line 361 "slasd4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 369 "slasd4.f"
	if (w * eta > 0.) {
#line 369 "slasd4.f"
	    eta = -w / (dpsi + dphi);
#line 369 "slasd4.f"
	}
#line 371 "slasd4.f"
	temp = eta - dtnsq;
#line 372 "slasd4.f"
	if (temp > *rho) {
#line 372 "slasd4.f"
	    eta = *rho + dtnsq;
#line 372 "slasd4.f"
	}

#line 375 "slasd4.f"
	eta /= *sigma + sqrt(eta + *sigma * *sigma);
#line 376 "slasd4.f"
	tau += eta;
#line 377 "slasd4.f"
	*sigma += eta;

#line 379 "slasd4.f"
	i__1 = *n;
#line 379 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 380 "slasd4.f"
	    delta[j] -= eta;
#line 381 "slasd4.f"
	    work[j] += eta;
#line 382 "slasd4.f"
/* L50: */
#line 382 "slasd4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 386 "slasd4.f"
	dpsi = 0.;
#line 387 "slasd4.f"
	psi = 0.;
#line 388 "slasd4.f"
	erretm = 0.;
#line 389 "slasd4.f"
	i__1 = ii;
#line 389 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 390 "slasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 391 "slasd4.f"
	    psi += z__[j] * temp;
#line 392 "slasd4.f"
	    dpsi += temp * temp;
#line 393 "slasd4.f"
	    erretm += psi;
#line 394 "slasd4.f"
/* L60: */
#line 394 "slasd4.f"
	}
#line 395 "slasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 399 "slasd4.f"
	tau2 = work[*n] * delta[*n];
#line 400 "slasd4.f"
	temp = z__[*n] / tau2;
#line 401 "slasd4.f"
	phi = z__[*n] * temp;
#line 402 "slasd4.f"
	dphi = temp * temp;
#line 403 "slasd4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $          + ABS( TAU2 )*( DPSI+DPHI ) */

#line 406 "slasd4.f"
	w = rhoinv + phi + psi;

/*        Main loop to update the values of the array   DELTA */

#line 410 "slasd4.f"
	iter = niter + 1;

#line 412 "slasd4.f"
	for (niter = iter; niter <= 400; ++niter) {

/*           Test for convergence */

#line 416 "slasd4.f"
	    if (abs(w) <= eps * erretm) {
#line 417 "slasd4.f"
		goto L240;
#line 418 "slasd4.f"
	    }

/*           Calculate the new step */

#line 422 "slasd4.f"
	    dtnsq1 = work[*n - 1] * delta[*n - 1];
#line 423 "slasd4.f"
	    dtnsq = work[*n] * delta[*n];
#line 424 "slasd4.f"
	    c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
#line 425 "slasd4.f"
	    a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
#line 426 "slasd4.f"
	    b = dtnsq1 * dtnsq * w;
#line 427 "slasd4.f"
	    if (a >= 0.) {
#line 428 "slasd4.f"
		eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 429 "slasd4.f"
	    } else {
#line 430 "slasd4.f"
		eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 431 "slasd4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 439 "slasd4.f"
	    if (w * eta > 0.) {
#line 439 "slasd4.f"
		eta = -w / (dpsi + dphi);
#line 439 "slasd4.f"
	    }
#line 441 "slasd4.f"
	    temp = eta - dtnsq;
#line 442 "slasd4.f"
	    if (temp <= 0.) {
#line 442 "slasd4.f"
		eta /= 2.;
#line 442 "slasd4.f"
	    }

#line 445 "slasd4.f"
	    eta /= *sigma + sqrt(eta + *sigma * *sigma);
#line 446 "slasd4.f"
	    tau += eta;
#line 447 "slasd4.f"
	    *sigma += eta;

#line 449 "slasd4.f"
	    i__1 = *n;
#line 449 "slasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 450 "slasd4.f"
		delta[j] -= eta;
#line 451 "slasd4.f"
		work[j] += eta;
#line 452 "slasd4.f"
/* L70: */
#line 452 "slasd4.f"
	    }

/*           Evaluate PSI and the derivative DPSI */

#line 456 "slasd4.f"
	    dpsi = 0.;
#line 457 "slasd4.f"
	    psi = 0.;
#line 458 "slasd4.f"
	    erretm = 0.;
#line 459 "slasd4.f"
	    i__1 = ii;
#line 459 "slasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 460 "slasd4.f"
		temp = z__[j] / (work[j] * delta[j]);
#line 461 "slasd4.f"
		psi += z__[j] * temp;
#line 462 "slasd4.f"
		dpsi += temp * temp;
#line 463 "slasd4.f"
		erretm += psi;
#line 464 "slasd4.f"
/* L80: */
#line 464 "slasd4.f"
	    }
#line 465 "slasd4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 469 "slasd4.f"
	    tau2 = work[*n] * delta[*n];
#line 470 "slasd4.f"
	    temp = z__[*n] / tau2;
#line 471 "slasd4.f"
	    phi = z__[*n] * temp;
#line 472 "slasd4.f"
	    dphi = temp * temp;
#line 473 "slasd4.f"
	    erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $             + ABS( TAU2 )*( DPSI+DPHI ) */

#line 476 "slasd4.f"
	    w = rhoinv + phi + psi;
#line 477 "slasd4.f"
/* L90: */
#line 477 "slasd4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 481 "slasd4.f"
	*info = 1;
#line 482 "slasd4.f"
	goto L240;

/*        End for the case I = N */

#line 486 "slasd4.f"
    } else {

/*        The case for I < N */

#line 490 "slasd4.f"
	niter = 1;
#line 491 "slasd4.f"
	ip1 = *i__ + 1;

/*        Calculate initial guess */

#line 495 "slasd4.f"
	delsq = (d__[ip1] - d__[*i__]) * (d__[ip1] + d__[*i__]);
#line 496 "slasd4.f"
	delsq2 = delsq / 2.;
#line 497 "slasd4.f"
	sq2 = sqrt((d__[*i__] * d__[*i__] + d__[ip1] * d__[ip1]) / 2.);
#line 498 "slasd4.f"
	temp = delsq2 / (d__[*i__] + sq2);
#line 499 "slasd4.f"
	i__1 = *n;
#line 499 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 500 "slasd4.f"
	    work[j] = d__[j] + d__[*i__] + temp;
#line 501 "slasd4.f"
	    delta[j] = d__[j] - d__[*i__] - temp;
#line 502 "slasd4.f"
/* L100: */
#line 502 "slasd4.f"
	}

#line 504 "slasd4.f"
	psi = 0.;
#line 505 "slasd4.f"
	i__1 = *i__ - 1;
#line 505 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 506 "slasd4.f"
	    psi += z__[j] * z__[j] / (work[j] * delta[j]);
#line 507 "slasd4.f"
/* L110: */
#line 507 "slasd4.f"
	}

#line 509 "slasd4.f"
	phi = 0.;
#line 510 "slasd4.f"
	i__1 = *i__ + 2;
#line 510 "slasd4.f"
	for (j = *n; j >= i__1; --j) {
#line 511 "slasd4.f"
	    phi += z__[j] * z__[j] / (work[j] * delta[j]);
#line 512 "slasd4.f"
/* L120: */
#line 512 "slasd4.f"
	}
#line 513 "slasd4.f"
	c__ = rhoinv + psi + phi;
#line 514 "slasd4.f"
	w = c__ + z__[*i__] * z__[*i__] / (work[*i__] * delta[*i__]) + z__[
		ip1] * z__[ip1] / (work[ip1] * delta[ip1]);

#line 517 "slasd4.f"
	geomavg = FALSE_;
#line 518 "slasd4.f"
	if (w > 0.) {

/*           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2 */

/*           We choose d(i) as origin. */

#line 524 "slasd4.f"
	    orgati = TRUE_;
#line 525 "slasd4.f"
	    ii = *i__;
#line 526 "slasd4.f"
	    sglb = 0.;
#line 527 "slasd4.f"
	    sgub = delsq2 / (d__[*i__] + sq2);
#line 528 "slasd4.f"
	    a = c__ * delsq + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
#line 529 "slasd4.f"
	    b = z__[*i__] * z__[*i__] * delsq;
#line 530 "slasd4.f"
	    if (a > 0.) {
#line 531 "slasd4.f"
		tau2 = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 532 "slasd4.f"
	    } else {
#line 533 "slasd4.f"
		tau2 = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / 
			(c__ * 2.);
#line 534 "slasd4.f"
	    }

/*           TAU2 now is an estimation of SIGMA^2 - D( I )^2. The */
/*           following, however, is the corresponding estimation of */
/*           SIGMA - D( I ). */

#line 540 "slasd4.f"
	    tau = tau2 / (d__[*i__] + sqrt(d__[*i__] * d__[*i__] + tau2));
#line 541 "slasd4.f"
	    temp = sqrt(eps);
#line 542 "slasd4.f"
	    if (d__[*i__] <= temp * d__[ip1] && (d__1 = z__[*i__], abs(d__1)) 
		    <= temp && d__[*i__] > 0.) {
/* Computing MIN */
#line 544 "slasd4.f"
		d__1 = d__[*i__] * 10.;
#line 544 "slasd4.f"
		tau = min(d__1,sgub);
#line 545 "slasd4.f"
		geomavg = TRUE_;
#line 546 "slasd4.f"
	    }
#line 547 "slasd4.f"
	} else {

/*           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2 */

/*           We choose d(i+1) as origin. */

#line 553 "slasd4.f"
	    orgati = FALSE_;
#line 554 "slasd4.f"
	    ii = ip1;
#line 555 "slasd4.f"
	    sglb = -delsq2 / (d__[ii] + sq2);
#line 556 "slasd4.f"
	    sgub = 0.;
#line 557 "slasd4.f"
	    a = c__ * delsq - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
#line 558 "slasd4.f"
	    b = z__[ip1] * z__[ip1] * delsq;
#line 559 "slasd4.f"
	    if (a < 0.) {
#line 560 "slasd4.f"
		tau2 = b * 2. / (a - sqrt((d__1 = a * a + b * 4. * c__, abs(
			d__1))));
#line 561 "slasd4.f"
	    } else {
#line 562 "slasd4.f"
		tau2 = -(a + sqrt((d__1 = a * a + b * 4. * c__, abs(d__1)))) /
			 (c__ * 2.);
#line 563 "slasd4.f"
	    }

/*           TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The */
/*           following, however, is the corresponding estimation of */
/*           SIGMA - D( IP1 ). */

#line 569 "slasd4.f"
	    tau = tau2 / (d__[ip1] + sqrt((d__1 = d__[ip1] * d__[ip1] + tau2, 
		    abs(d__1))));
#line 571 "slasd4.f"
	}

#line 573 "slasd4.f"
	*sigma = d__[ii] + tau;
#line 574 "slasd4.f"
	i__1 = *n;
#line 574 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 575 "slasd4.f"
	    work[j] = d__[j] + d__[ii] + tau;
#line 576 "slasd4.f"
	    delta[j] = d__[j] - d__[ii] - tau;
#line 577 "slasd4.f"
/* L130: */
#line 577 "slasd4.f"
	}
#line 578 "slasd4.f"
	iim1 = ii - 1;
#line 579 "slasd4.f"
	iip1 = ii + 1;

/*        Evaluate PSI and the derivative DPSI */

#line 583 "slasd4.f"
	dpsi = 0.;
#line 584 "slasd4.f"
	psi = 0.;
#line 585 "slasd4.f"
	erretm = 0.;
#line 586 "slasd4.f"
	i__1 = iim1;
#line 586 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 587 "slasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 588 "slasd4.f"
	    psi += z__[j] * temp;
#line 589 "slasd4.f"
	    dpsi += temp * temp;
#line 590 "slasd4.f"
	    erretm += psi;
#line 591 "slasd4.f"
/* L150: */
#line 591 "slasd4.f"
	}
#line 592 "slasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 596 "slasd4.f"
	dphi = 0.;
#line 597 "slasd4.f"
	phi = 0.;
#line 598 "slasd4.f"
	i__1 = iip1;
#line 598 "slasd4.f"
	for (j = *n; j >= i__1; --j) {
#line 599 "slasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 600 "slasd4.f"
	    phi += z__[j] * temp;
#line 601 "slasd4.f"
	    dphi += temp * temp;
#line 602 "slasd4.f"
	    erretm += phi;
#line 603 "slasd4.f"
/* L160: */
#line 603 "slasd4.f"
	}

#line 605 "slasd4.f"
	w = rhoinv + phi + psi;

/*        W is the value of the secular function with */
/*        its ii-th element removed. */

#line 610 "slasd4.f"
	swtch3 = FALSE_;
#line 611 "slasd4.f"
	if (orgati) {
#line 612 "slasd4.f"
	    if (w < 0.) {
#line 612 "slasd4.f"
		swtch3 = TRUE_;
#line 612 "slasd4.f"
	    }
#line 614 "slasd4.f"
	} else {
#line 615 "slasd4.f"
	    if (w > 0.) {
#line 615 "slasd4.f"
		swtch3 = TRUE_;
#line 615 "slasd4.f"
	    }
#line 617 "slasd4.f"
	}
#line 618 "slasd4.f"
	if (ii == 1 || ii == *n) {
#line 618 "slasd4.f"
	    swtch3 = FALSE_;
#line 618 "slasd4.f"
	}

#line 621 "slasd4.f"
	temp = z__[ii] / (work[ii] * delta[ii]);
#line 622 "slasd4.f"
	dw = dpsi + dphi + temp * temp;
#line 623 "slasd4.f"
	temp = z__[ii] * temp;
#line 624 "slasd4.f"
	w += temp;
#line 625 "slasd4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $          + ABS( TAU2 )*DW */

/*        Test for convergence */

#line 631 "slasd4.f"
	if (abs(w) <= eps * erretm) {
#line 632 "slasd4.f"
	    goto L240;
#line 633 "slasd4.f"
	}

#line 635 "slasd4.f"
	if (w <= 0.) {
#line 636 "slasd4.f"
	    sglb = max(sglb,tau);
#line 637 "slasd4.f"
	} else {
#line 638 "slasd4.f"
	    sgub = min(sgub,tau);
#line 639 "slasd4.f"
	}

/*        Calculate the new step */

#line 643 "slasd4.f"
	++niter;
#line 644 "slasd4.f"
	if (! swtch3) {
#line 645 "slasd4.f"
	    dtipsq = work[ip1] * delta[ip1];
#line 646 "slasd4.f"
	    dtisq = work[*i__] * delta[*i__];
#line 647 "slasd4.f"
	    if (orgati) {
/* Computing 2nd power */
#line 648 "slasd4.f"
		d__1 = z__[*i__] / dtisq;
#line 648 "slasd4.f"
		c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 649 "slasd4.f"
	    } else {
/* Computing 2nd power */
#line 650 "slasd4.f"
		d__1 = z__[ip1] / dtipsq;
#line 650 "slasd4.f"
		c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 651 "slasd4.f"
	    }
#line 652 "slasd4.f"
	    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 653 "slasd4.f"
	    b = dtipsq * dtisq * w;
#line 654 "slasd4.f"
	    if (c__ == 0.) {
#line 655 "slasd4.f"
		if (a == 0.) {
#line 656 "slasd4.f"
		    if (orgati) {
#line 657 "slasd4.f"
			a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + 
				dphi);
#line 658 "slasd4.f"
		    } else {
#line 659 "slasd4.f"
			a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + 
				dphi);
#line 660 "slasd4.f"
		    }
#line 661 "slasd4.f"
		}
#line 662 "slasd4.f"
		eta = b / a;
#line 663 "slasd4.f"
	    } else if (a <= 0.) {
#line 664 "slasd4.f"
		eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 665 "slasd4.f"
	    } else {
#line 666 "slasd4.f"
		eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 667 "slasd4.f"
	    }
#line 668 "slasd4.f"
	} else {

/*           Interpolation using THREE most relevant poles */

#line 672 "slasd4.f"
	    dtiim = work[iim1] * delta[iim1];
#line 673 "slasd4.f"
	    dtiip = work[iip1] * delta[iip1];
#line 674 "slasd4.f"
	    temp = rhoinv + psi + phi;
#line 675 "slasd4.f"
	    if (orgati) {
#line 676 "slasd4.f"
		temp1 = z__[iim1] / dtiim;
#line 677 "slasd4.f"
		temp1 *= temp1;
#line 678 "slasd4.f"
		c__ = temp - dtiip * (dpsi + dphi) - (d__[iim1] - d__[iip1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
#line 680 "slasd4.f"
		zz[0] = z__[iim1] * z__[iim1];
#line 681 "slasd4.f"
		if (dpsi < temp1) {
#line 682 "slasd4.f"
		    zz[2] = dtiip * dtiip * dphi;
#line 683 "slasd4.f"
		} else {
#line 684 "slasd4.f"
		    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
#line 685 "slasd4.f"
		}
#line 686 "slasd4.f"
	    } else {
#line 687 "slasd4.f"
		temp1 = z__[iip1] / dtiip;
#line 688 "slasd4.f"
		temp1 *= temp1;
#line 689 "slasd4.f"
		c__ = temp - dtiim * (dpsi + dphi) - (d__[iip1] - d__[iim1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
#line 691 "slasd4.f"
		if (dphi < temp1) {
#line 692 "slasd4.f"
		    zz[0] = dtiim * dtiim * dpsi;
#line 693 "slasd4.f"
		} else {
#line 694 "slasd4.f"
		    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
#line 695 "slasd4.f"
		}
#line 696 "slasd4.f"
		zz[2] = z__[iip1] * z__[iip1];
#line 697 "slasd4.f"
	    }
#line 698 "slasd4.f"
	    zz[1] = z__[ii] * z__[ii];
#line 699 "slasd4.f"
	    dd[0] = dtiim;
#line 700 "slasd4.f"
	    dd[1] = delta[ii] * work[ii];
#line 701 "slasd4.f"
	    dd[2] = dtiip;
#line 702 "slasd4.f"
	    slaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);

#line 704 "slasd4.f"
	    if (*info != 0) {

/*              If INFO is not 0, i.e., SLAED6 failed, switch back */
/*              to 2 pole interpolation. */

#line 709 "slasd4.f"
		swtch3 = FALSE_;
#line 710 "slasd4.f"
		*info = 0;
#line 711 "slasd4.f"
		dtipsq = work[ip1] * delta[ip1];
#line 712 "slasd4.f"
		dtisq = work[*i__] * delta[*i__];
#line 713 "slasd4.f"
		if (orgati) {
/* Computing 2nd power */
#line 714 "slasd4.f"
		    d__1 = z__[*i__] / dtisq;
#line 714 "slasd4.f"
		    c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 715 "slasd4.f"
		} else {
/* Computing 2nd power */
#line 716 "slasd4.f"
		    d__1 = z__[ip1] / dtipsq;
#line 716 "slasd4.f"
		    c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 717 "slasd4.f"
		}
#line 718 "slasd4.f"
		a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 719 "slasd4.f"
		b = dtipsq * dtisq * w;
#line 720 "slasd4.f"
		if (c__ == 0.) {
#line 721 "slasd4.f"
		    if (a == 0.) {
#line 722 "slasd4.f"
			if (orgati) {
#line 723 "slasd4.f"
			    a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (
				    dpsi + dphi);
#line 724 "slasd4.f"
			} else {
#line 725 "slasd4.f"
			    a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + 
				    dphi);
#line 726 "slasd4.f"
			}
#line 727 "slasd4.f"
		    }
#line 728 "slasd4.f"
		    eta = b / a;
#line 729 "slasd4.f"
		} else if (a <= 0.) {
#line 730 "slasd4.f"
		    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
			     / (c__ * 2.);
#line 731 "slasd4.f"
		} else {
#line 732 "slasd4.f"
		    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, 
			    abs(d__1))));
#line 733 "slasd4.f"
		}
#line 734 "slasd4.f"
	    }
#line 735 "slasd4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 743 "slasd4.f"
	if (w * eta >= 0.) {
#line 743 "slasd4.f"
	    eta = -w / dw;
#line 743 "slasd4.f"
	}

#line 746 "slasd4.f"
	eta /= *sigma + sqrt(*sigma * *sigma + eta);
#line 747 "slasd4.f"
	temp = tau + eta;
#line 748 "slasd4.f"
	if (temp > sgub || temp < sglb) {
#line 749 "slasd4.f"
	    if (w < 0.) {
#line 750 "slasd4.f"
		eta = (sgub - tau) / 2.;
#line 751 "slasd4.f"
	    } else {
#line 752 "slasd4.f"
		eta = (sglb - tau) / 2.;
#line 753 "slasd4.f"
	    }
#line 754 "slasd4.f"
	    if (geomavg) {
#line 755 "slasd4.f"
		if (w < 0.) {
#line 756 "slasd4.f"
		    if (tau > 0.) {
#line 757 "slasd4.f"
			eta = sqrt(sgub * tau) - tau;
#line 758 "slasd4.f"
		    }
#line 759 "slasd4.f"
		} else {
#line 760 "slasd4.f"
		    if (sglb > 0.) {
#line 761 "slasd4.f"
			eta = sqrt(sglb * tau) - tau;
#line 762 "slasd4.f"
		    }
#line 763 "slasd4.f"
		}
#line 764 "slasd4.f"
	    }
#line 765 "slasd4.f"
	}

#line 767 "slasd4.f"
	prew = w;

#line 769 "slasd4.f"
	tau += eta;
#line 770 "slasd4.f"
	*sigma += eta;

#line 772 "slasd4.f"
	i__1 = *n;
#line 772 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 773 "slasd4.f"
	    work[j] += eta;
#line 774 "slasd4.f"
	    delta[j] -= eta;
#line 775 "slasd4.f"
/* L170: */
#line 775 "slasd4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 779 "slasd4.f"
	dpsi = 0.;
#line 780 "slasd4.f"
	psi = 0.;
#line 781 "slasd4.f"
	erretm = 0.;
#line 782 "slasd4.f"
	i__1 = iim1;
#line 782 "slasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 783 "slasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 784 "slasd4.f"
	    psi += z__[j] * temp;
#line 785 "slasd4.f"
	    dpsi += temp * temp;
#line 786 "slasd4.f"
	    erretm += psi;
#line 787 "slasd4.f"
/* L180: */
#line 787 "slasd4.f"
	}
#line 788 "slasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 792 "slasd4.f"
	dphi = 0.;
#line 793 "slasd4.f"
	phi = 0.;
#line 794 "slasd4.f"
	i__1 = iip1;
#line 794 "slasd4.f"
	for (j = *n; j >= i__1; --j) {
#line 795 "slasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 796 "slasd4.f"
	    phi += z__[j] * temp;
#line 797 "slasd4.f"
	    dphi += temp * temp;
#line 798 "slasd4.f"
	    erretm += phi;
#line 799 "slasd4.f"
/* L190: */
#line 799 "slasd4.f"
	}

#line 801 "slasd4.f"
	tau2 = work[ii] * delta[ii];
#line 802 "slasd4.f"
	temp = z__[ii] / tau2;
#line 803 "slasd4.f"
	dw = dpsi + dphi + temp * temp;
#line 804 "slasd4.f"
	temp = z__[ii] * temp;
#line 805 "slasd4.f"
	w = rhoinv + phi + psi + temp;
#line 806 "slasd4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $          + ABS( TAU2 )*DW */

#line 810 "slasd4.f"
	swtch = FALSE_;
#line 811 "slasd4.f"
	if (orgati) {
#line 812 "slasd4.f"
	    if (-w > abs(prew) / 10.) {
#line 812 "slasd4.f"
		swtch = TRUE_;
#line 812 "slasd4.f"
	    }
#line 814 "slasd4.f"
	} else {
#line 815 "slasd4.f"
	    if (w > abs(prew) / 10.) {
#line 815 "slasd4.f"
		swtch = TRUE_;
#line 815 "slasd4.f"
	    }
#line 817 "slasd4.f"
	}

/*        Main loop to update the values of the array   DELTA and WORK */

#line 821 "slasd4.f"
	iter = niter + 1;

#line 823 "slasd4.f"
	for (niter = iter; niter <= 400; ++niter) {

/*           Test for convergence */

#line 827 "slasd4.f"
	    if (abs(w) <= eps * erretm) {
/*     $          .OR. (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN */
#line 829 "slasd4.f"
		goto L240;
#line 830 "slasd4.f"
	    }

#line 832 "slasd4.f"
	    if (w <= 0.) {
#line 833 "slasd4.f"
		sglb = max(sglb,tau);
#line 834 "slasd4.f"
	    } else {
#line 835 "slasd4.f"
		sgub = min(sgub,tau);
#line 836 "slasd4.f"
	    }

/*           Calculate the new step */

#line 840 "slasd4.f"
	    if (! swtch3) {
#line 841 "slasd4.f"
		dtipsq = work[ip1] * delta[ip1];
#line 842 "slasd4.f"
		dtisq = work[*i__] * delta[*i__];
#line 843 "slasd4.f"
		if (! swtch) {
#line 844 "slasd4.f"
		    if (orgati) {
/* Computing 2nd power */
#line 845 "slasd4.f"
			d__1 = z__[*i__] / dtisq;
#line 845 "slasd4.f"
			c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 846 "slasd4.f"
		    } else {
/* Computing 2nd power */
#line 847 "slasd4.f"
			d__1 = z__[ip1] / dtipsq;
#line 847 "slasd4.f"
			c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 848 "slasd4.f"
		    }
#line 849 "slasd4.f"
		} else {
#line 850 "slasd4.f"
		    temp = z__[ii] / (work[ii] * delta[ii]);
#line 851 "slasd4.f"
		    if (orgati) {
#line 852 "slasd4.f"
			dpsi += temp * temp;
#line 853 "slasd4.f"
		    } else {
#line 854 "slasd4.f"
			dphi += temp * temp;
#line 855 "slasd4.f"
		    }
#line 856 "slasd4.f"
		    c__ = w - dtisq * dpsi - dtipsq * dphi;
#line 857 "slasd4.f"
		}
#line 858 "slasd4.f"
		a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 859 "slasd4.f"
		b = dtipsq * dtisq * w;
#line 860 "slasd4.f"
		if (c__ == 0.) {
#line 861 "slasd4.f"
		    if (a == 0.) {
#line 862 "slasd4.f"
			if (! swtch) {
#line 863 "slasd4.f"
			    if (orgati) {
#line 864 "slasd4.f"
				a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * 
					(dpsi + dphi);
#line 866 "slasd4.f"
			    } else {
#line 867 "slasd4.f"
				a = z__[ip1] * z__[ip1] + dtisq * dtisq * (
					dpsi + dphi);
#line 869 "slasd4.f"
			    }
#line 870 "slasd4.f"
			} else {
#line 871 "slasd4.f"
			    a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
#line 872 "slasd4.f"
			}
#line 873 "slasd4.f"
		    }
#line 874 "slasd4.f"
		    eta = b / a;
#line 875 "slasd4.f"
		} else if (a <= 0.) {
#line 876 "slasd4.f"
		    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
			     / (c__ * 2.);
#line 877 "slasd4.f"
		} else {
#line 878 "slasd4.f"
		    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, 
			    abs(d__1))));
#line 879 "slasd4.f"
		}
#line 880 "slasd4.f"
	    } else {

/*              Interpolation using THREE most relevant poles */

#line 884 "slasd4.f"
		dtiim = work[iim1] * delta[iim1];
#line 885 "slasd4.f"
		dtiip = work[iip1] * delta[iip1];
#line 886 "slasd4.f"
		temp = rhoinv + psi + phi;
#line 887 "slasd4.f"
		if (swtch) {
#line 888 "slasd4.f"
		    c__ = temp - dtiim * dpsi - dtiip * dphi;
#line 889 "slasd4.f"
		    zz[0] = dtiim * dtiim * dpsi;
#line 890 "slasd4.f"
		    zz[2] = dtiip * dtiip * dphi;
#line 891 "slasd4.f"
		} else {
#line 892 "slasd4.f"
		    if (orgati) {
#line 893 "slasd4.f"
			temp1 = z__[iim1] / dtiim;
#line 894 "slasd4.f"
			temp1 *= temp1;
#line 895 "slasd4.f"
			temp2 = (d__[iim1] - d__[iip1]) * (d__[iim1] + d__[
				iip1]) * temp1;
#line 897 "slasd4.f"
			c__ = temp - dtiip * (dpsi + dphi) - temp2;
#line 898 "slasd4.f"
			zz[0] = z__[iim1] * z__[iim1];
#line 899 "slasd4.f"
			if (dpsi < temp1) {
#line 900 "slasd4.f"
			    zz[2] = dtiip * dtiip * dphi;
#line 901 "slasd4.f"
			} else {
#line 902 "slasd4.f"
			    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
#line 903 "slasd4.f"
			}
#line 904 "slasd4.f"
		    } else {
#line 905 "slasd4.f"
			temp1 = z__[iip1] / dtiip;
#line 906 "slasd4.f"
			temp1 *= temp1;
#line 907 "slasd4.f"
			temp2 = (d__[iip1] - d__[iim1]) * (d__[iim1] + d__[
				iip1]) * temp1;
#line 909 "slasd4.f"
			c__ = temp - dtiim * (dpsi + dphi) - temp2;
#line 910 "slasd4.f"
			if (dphi < temp1) {
#line 911 "slasd4.f"
			    zz[0] = dtiim * dtiim * dpsi;
#line 912 "slasd4.f"
			} else {
#line 913 "slasd4.f"
			    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
#line 914 "slasd4.f"
			}
#line 915 "slasd4.f"
			zz[2] = z__[iip1] * z__[iip1];
#line 916 "slasd4.f"
		    }
#line 917 "slasd4.f"
		}
#line 918 "slasd4.f"
		dd[0] = dtiim;
#line 919 "slasd4.f"
		dd[1] = delta[ii] * work[ii];
#line 920 "slasd4.f"
		dd[2] = dtiip;
#line 921 "slasd4.f"
		slaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);

#line 923 "slasd4.f"
		if (*info != 0) {

/*                 If INFO is not 0, i.e., SLAED6 failed, switch */
/*                 back to two pole interpolation */

#line 928 "slasd4.f"
		    swtch3 = FALSE_;
#line 929 "slasd4.f"
		    *info = 0;
#line 930 "slasd4.f"
		    dtipsq = work[ip1] * delta[ip1];
#line 931 "slasd4.f"
		    dtisq = work[*i__] * delta[*i__];
#line 932 "slasd4.f"
		    if (! swtch) {
#line 933 "slasd4.f"
			if (orgati) {
/* Computing 2nd power */
#line 934 "slasd4.f"
			    d__1 = z__[*i__] / dtisq;
#line 934 "slasd4.f"
			    c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 935 "slasd4.f"
			} else {
/* Computing 2nd power */
#line 936 "slasd4.f"
			    d__1 = z__[ip1] / dtipsq;
#line 936 "slasd4.f"
			    c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 937 "slasd4.f"
			}
#line 938 "slasd4.f"
		    } else {
#line 939 "slasd4.f"
			temp = z__[ii] / (work[ii] * delta[ii]);
#line 940 "slasd4.f"
			if (orgati) {
#line 941 "slasd4.f"
			    dpsi += temp * temp;
#line 942 "slasd4.f"
			} else {
#line 943 "slasd4.f"
			    dphi += temp * temp;
#line 944 "slasd4.f"
			}
#line 945 "slasd4.f"
			c__ = w - dtisq * dpsi - dtipsq * dphi;
#line 946 "slasd4.f"
		    }
#line 947 "slasd4.f"
		    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 948 "slasd4.f"
		    b = dtipsq * dtisq * w;
#line 949 "slasd4.f"
		    if (c__ == 0.) {
#line 950 "slasd4.f"
			if (a == 0.) {
#line 951 "slasd4.f"
			    if (! swtch) {
#line 952 "slasd4.f"
				if (orgati) {
#line 953 "slasd4.f"
				    a = z__[*i__] * z__[*i__] + dtipsq * 
					    dtipsq * (dpsi + dphi);
#line 955 "slasd4.f"
				} else {
#line 956 "slasd4.f"
				    a = z__[ip1] * z__[ip1] + dtisq * dtisq * 
					    (dpsi + dphi);
#line 958 "slasd4.f"
				}
#line 959 "slasd4.f"
			    } else {
#line 960 "slasd4.f"
				a = dtisq * dtisq * dpsi + dtipsq * dtipsq * 
					dphi;
#line 961 "slasd4.f"
			    }
#line 962 "slasd4.f"
			}
#line 963 "slasd4.f"
			eta = b / a;
#line 964 "slasd4.f"
		    } else if (a <= 0.) {
#line 965 "slasd4.f"
			eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
				d__1)))) / (c__ * 2.);
#line 966 "slasd4.f"
		    } else {
#line 967 "slasd4.f"
			eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__,
				 abs(d__1))));
#line 968 "slasd4.f"
		    }
#line 969 "slasd4.f"
		}
#line 970 "slasd4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 978 "slasd4.f"
	    if (w * eta >= 0.) {
#line 978 "slasd4.f"
		eta = -w / dw;
#line 978 "slasd4.f"
	    }

#line 981 "slasd4.f"
	    eta /= *sigma + sqrt(*sigma * *sigma + eta);
#line 982 "slasd4.f"
	    temp = tau + eta;
#line 983 "slasd4.f"
	    if (temp > sgub || temp < sglb) {
#line 984 "slasd4.f"
		if (w < 0.) {
#line 985 "slasd4.f"
		    eta = (sgub - tau) / 2.;
#line 986 "slasd4.f"
		} else {
#line 987 "slasd4.f"
		    eta = (sglb - tau) / 2.;
#line 988 "slasd4.f"
		}
#line 989 "slasd4.f"
		if (geomavg) {
#line 990 "slasd4.f"
		    if (w < 0.) {
#line 991 "slasd4.f"
			if (tau > 0.) {
#line 992 "slasd4.f"
			    eta = sqrt(sgub * tau) - tau;
#line 993 "slasd4.f"
			}
#line 994 "slasd4.f"
		    } else {
#line 995 "slasd4.f"
			if (sglb > 0.) {
#line 996 "slasd4.f"
			    eta = sqrt(sglb * tau) - tau;
#line 997 "slasd4.f"
			}
#line 998 "slasd4.f"
		    }
#line 999 "slasd4.f"
		}
#line 1000 "slasd4.f"
	    }

#line 1002 "slasd4.f"
	    prew = w;

#line 1004 "slasd4.f"
	    tau += eta;
#line 1005 "slasd4.f"
	    *sigma += eta;

#line 1007 "slasd4.f"
	    i__1 = *n;
#line 1007 "slasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1008 "slasd4.f"
		work[j] += eta;
#line 1009 "slasd4.f"
		delta[j] -= eta;
#line 1010 "slasd4.f"
/* L200: */
#line 1010 "slasd4.f"
	    }

/*           Evaluate PSI and the derivative DPSI */

#line 1014 "slasd4.f"
	    dpsi = 0.;
#line 1015 "slasd4.f"
	    psi = 0.;
#line 1016 "slasd4.f"
	    erretm = 0.;
#line 1017 "slasd4.f"
	    i__1 = iim1;
#line 1017 "slasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1018 "slasd4.f"
		temp = z__[j] / (work[j] * delta[j]);
#line 1019 "slasd4.f"
		psi += z__[j] * temp;
#line 1020 "slasd4.f"
		dpsi += temp * temp;
#line 1021 "slasd4.f"
		erretm += psi;
#line 1022 "slasd4.f"
/* L210: */
#line 1022 "slasd4.f"
	    }
#line 1023 "slasd4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 1027 "slasd4.f"
	    dphi = 0.;
#line 1028 "slasd4.f"
	    phi = 0.;
#line 1029 "slasd4.f"
	    i__1 = iip1;
#line 1029 "slasd4.f"
	    for (j = *n; j >= i__1; --j) {
#line 1030 "slasd4.f"
		temp = z__[j] / (work[j] * delta[j]);
#line 1031 "slasd4.f"
		phi += z__[j] * temp;
#line 1032 "slasd4.f"
		dphi += temp * temp;
#line 1033 "slasd4.f"
		erretm += phi;
#line 1034 "slasd4.f"
/* L220: */
#line 1034 "slasd4.f"
	    }

#line 1036 "slasd4.f"
	    tau2 = work[ii] * delta[ii];
#line 1037 "slasd4.f"
	    temp = z__[ii] / tau2;
#line 1038 "slasd4.f"
	    dw = dpsi + dphi + temp * temp;
#line 1039 "slasd4.f"
	    temp = z__[ii] * temp;
#line 1040 "slasd4.f"
	    w = rhoinv + phi + psi + temp;
#line 1041 "slasd4.f"
	    erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $             + ABS( TAU2 )*DW */

#line 1045 "slasd4.f"
	    if (w * prew > 0. && abs(w) > abs(prew) / 10.) {
#line 1045 "slasd4.f"
		swtch = ! swtch;
#line 1045 "slasd4.f"
	    }

#line 1048 "slasd4.f"
/* L230: */
#line 1048 "slasd4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 1052 "slasd4.f"
	*info = 1;

#line 1054 "slasd4.f"
    }

#line 1056 "slasd4.f"
L240:
#line 1057 "slasd4.f"
    return 0;

/*     End of SLASD4 */

} /* slasd4_ */

