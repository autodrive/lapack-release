#line 1 "dlasd4.f"
/* dlasd4.f -- translated by f2c (version 20100827).
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

#line 1 "dlasd4.f"
/* > \brief \b DLASD4 computes the square root of the i-th updated eigenvalue of a positive symmetric rank-one
 modification to a positive diagonal matrix. Used by dbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASD4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I, INFO, N */
/*       DOUBLE PRECISION   RHO, SIGMA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), DELTA( * ), WORK( * ), Z( * ) */
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
/* >          D is DOUBLE PRECISION array, dimension ( N ) */
/* >         The original eigenvalues.  It is assumed that they are in */
/* >         order, 0 <= D(I) < D(J)  for I < J. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( N ) */
/* >         The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* >          DELTA is DOUBLE PRECISION array, dimension ( N ) */
/* >         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th */
/* >         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA */
/* >         contains the information necessary to construct the */
/* >         (singular) eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >         The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >         The computed sigma_I, the I-th updated eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension ( N ) */
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
/* Subroutine */ int dlasd4_(integer *n, integer *i__, doublereal *d__, 
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
    extern /* Subroutine */ int dlaed6_(integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , dlasd5_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal delsq2, dtnsq1;
    static logical swtch3;
    extern doublereal dlamch_(char *, ftnlen);
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

#line 207 "dlasd4.f"
    /* Parameter adjustments */
#line 207 "dlasd4.f"
    --work;
#line 207 "dlasd4.f"
    --delta;
#line 207 "dlasd4.f"
    --z__;
#line 207 "dlasd4.f"
    --d__;
#line 207 "dlasd4.f"

#line 207 "dlasd4.f"
    /* Function Body */
#line 207 "dlasd4.f"
    *info = 0;
#line 208 "dlasd4.f"
    if (*n == 1) {

/*        Presumably, I=1 upon entry */

#line 212 "dlasd4.f"
	*sigma = sqrt(d__[1] * d__[1] + *rho * z__[1] * z__[1]);
#line 213 "dlasd4.f"
	delta[1] = 1.;
#line 214 "dlasd4.f"
	work[1] = 1.;
#line 215 "dlasd4.f"
	return 0;
#line 216 "dlasd4.f"
    }
#line 217 "dlasd4.f"
    if (*n == 2) {
#line 218 "dlasd4.f"
	dlasd5_(i__, &d__[1], &z__[1], &delta[1], rho, sigma, &work[1]);
#line 219 "dlasd4.f"
	return 0;
#line 220 "dlasd4.f"
    }

/*     Compute machine epsilon */

#line 224 "dlasd4.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 225 "dlasd4.f"
    rhoinv = 1. / *rho;
#line 226 "dlasd4.f"
    tau2 = 0.;

/*     The case I = N */

#line 230 "dlasd4.f"
    if (*i__ == *n) {

/*        Initialize some basic variables */

#line 234 "dlasd4.f"
	ii = *n - 1;
#line 235 "dlasd4.f"
	niter = 1;

/*        Calculate initial guess */

#line 239 "dlasd4.f"
	temp = *rho / 2.;

/*        If ||Z||_2 is not one, then TEMP should be set to */
/*        RHO * ||Z||_2^2 / TWO */

#line 244 "dlasd4.f"
	temp1 = temp / (d__[*n] + sqrt(d__[*n] * d__[*n] + temp));
#line 245 "dlasd4.f"
	i__1 = *n;
#line 245 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 246 "dlasd4.f"
	    work[j] = d__[j] + d__[*n] + temp1;
#line 247 "dlasd4.f"
	    delta[j] = d__[j] - d__[*n] - temp1;
#line 248 "dlasd4.f"
/* L10: */
#line 248 "dlasd4.f"
	}

#line 250 "dlasd4.f"
	psi = 0.;
#line 251 "dlasd4.f"
	i__1 = *n - 2;
#line 251 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 252 "dlasd4.f"
	    psi += z__[j] * z__[j] / (delta[j] * work[j]);
#line 253 "dlasd4.f"
/* L20: */
#line 253 "dlasd4.f"
	}

#line 255 "dlasd4.f"
	c__ = rhoinv + psi;
#line 256 "dlasd4.f"
	w = c__ + z__[ii] * z__[ii] / (delta[ii] * work[ii]) + z__[*n] * z__[*
		n] / (delta[*n] * work[*n]);

#line 259 "dlasd4.f"
	if (w <= 0.) {
#line 260 "dlasd4.f"
	    temp1 = sqrt(d__[*n] * d__[*n] + *rho);
#line 261 "dlasd4.f"
	    temp = z__[*n - 1] * z__[*n - 1] / ((d__[*n - 1] + temp1) * (d__[*
		    n] - d__[*n - 1] + *rho / (d__[*n] + temp1))) + z__[*n] * 
		    z__[*n] / *rho;

/*           The following TAU2 is to approximate */
/*           SIGMA_n^2 - D( N )*D( N ) */

#line 268 "dlasd4.f"
	    if (c__ <= temp) {
#line 269 "dlasd4.f"
		tau = *rho;
#line 270 "dlasd4.f"
	    } else {
#line 271 "dlasd4.f"
		delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
#line 272 "dlasd4.f"
		a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*
			n];
#line 273 "dlasd4.f"
		b = z__[*n] * z__[*n] * delsq;
#line 274 "dlasd4.f"
		if (a < 0.) {
#line 275 "dlasd4.f"
		    tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 276 "dlasd4.f"
		} else {
#line 277 "dlasd4.f"
		    tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 278 "dlasd4.f"
		}
#line 279 "dlasd4.f"
		tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));
#line 280 "dlasd4.f"
	    }

/*           It can be proved that */
/*               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO */

#line 285 "dlasd4.f"
	} else {
#line 286 "dlasd4.f"
	    delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
#line 287 "dlasd4.f"
	    a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
#line 288 "dlasd4.f"
	    b = z__[*n] * z__[*n] * delsq;

/*           The following TAU2 is to approximate */
/*           SIGMA_n^2 - D( N )*D( N ) */

#line 293 "dlasd4.f"
	    if (a < 0.) {
#line 294 "dlasd4.f"
		tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
#line 295 "dlasd4.f"
	    } else {
#line 296 "dlasd4.f"
		tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
#line 297 "dlasd4.f"
	    }
#line 298 "dlasd4.f"
	    tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));

/*           It can be proved that */
/*           D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2 */

#line 304 "dlasd4.f"
	}

/*        The following TAU is to approximate SIGMA_n - D( N ) */

/*         TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) ) */

#line 310 "dlasd4.f"
	*sigma = d__[*n] + tau;
#line 311 "dlasd4.f"
	i__1 = *n;
#line 311 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 312 "dlasd4.f"
	    delta[j] = d__[j] - d__[*n] - tau;
#line 313 "dlasd4.f"
	    work[j] = d__[j] + d__[*n] + tau;
#line 314 "dlasd4.f"
/* L30: */
#line 314 "dlasd4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 318 "dlasd4.f"
	dpsi = 0.;
#line 319 "dlasd4.f"
	psi = 0.;
#line 320 "dlasd4.f"
	erretm = 0.;
#line 321 "dlasd4.f"
	i__1 = ii;
#line 321 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 322 "dlasd4.f"
	    temp = z__[j] / (delta[j] * work[j]);
#line 323 "dlasd4.f"
	    psi += z__[j] * temp;
#line 324 "dlasd4.f"
	    dpsi += temp * temp;
#line 325 "dlasd4.f"
	    erretm += psi;
#line 326 "dlasd4.f"
/* L40: */
#line 326 "dlasd4.f"
	}
#line 327 "dlasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 331 "dlasd4.f"
	temp = z__[*n] / (delta[*n] * work[*n]);
#line 332 "dlasd4.f"
	phi = z__[*n] * temp;
#line 333 "dlasd4.f"
	dphi = temp * temp;
#line 334 "dlasd4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $          + ABS( TAU2 )*( DPSI+DPHI ) */

#line 337 "dlasd4.f"
	w = rhoinv + phi + psi;

/*        Test for convergence */

#line 341 "dlasd4.f"
	if (abs(w) <= eps * erretm) {
#line 342 "dlasd4.f"
	    goto L240;
#line 343 "dlasd4.f"
	}

/*        Calculate the new step */

#line 347 "dlasd4.f"
	++niter;
#line 348 "dlasd4.f"
	dtnsq1 = work[*n - 1] * delta[*n - 1];
#line 349 "dlasd4.f"
	dtnsq = work[*n] * delta[*n];
#line 350 "dlasd4.f"
	c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
#line 351 "dlasd4.f"
	a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
#line 352 "dlasd4.f"
	b = dtnsq * dtnsq1 * w;
#line 353 "dlasd4.f"
	if (c__ < 0.) {
#line 353 "dlasd4.f"
	    c__ = abs(c__);
#line 353 "dlasd4.f"
	}
#line 355 "dlasd4.f"
	if (c__ == 0.) {
#line 356 "dlasd4.f"
	    eta = *rho - *sigma * *sigma;
#line 357 "dlasd4.f"
	} else if (a >= 0.) {
#line 358 "dlasd4.f"
	    eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ 
		    * 2.);
#line 359 "dlasd4.f"
	} else {
#line 360 "dlasd4.f"
	    eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))
		    );
#line 361 "dlasd4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 369 "dlasd4.f"
	if (w * eta > 0.) {
#line 369 "dlasd4.f"
	    eta = -w / (dpsi + dphi);
#line 369 "dlasd4.f"
	}
#line 371 "dlasd4.f"
	temp = eta - dtnsq;
#line 372 "dlasd4.f"
	if (temp > *rho) {
#line 372 "dlasd4.f"
	    eta = *rho + dtnsq;
#line 372 "dlasd4.f"
	}

#line 375 "dlasd4.f"
	eta /= *sigma + sqrt(eta + *sigma * *sigma);
#line 376 "dlasd4.f"
	tau += eta;
#line 377 "dlasd4.f"
	*sigma += eta;

#line 379 "dlasd4.f"
	i__1 = *n;
#line 379 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 380 "dlasd4.f"
	    delta[j] -= eta;
#line 381 "dlasd4.f"
	    work[j] += eta;
#line 382 "dlasd4.f"
/* L50: */
#line 382 "dlasd4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 386 "dlasd4.f"
	dpsi = 0.;
#line 387 "dlasd4.f"
	psi = 0.;
#line 388 "dlasd4.f"
	erretm = 0.;
#line 389 "dlasd4.f"
	i__1 = ii;
#line 389 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 390 "dlasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 391 "dlasd4.f"
	    psi += z__[j] * temp;
#line 392 "dlasd4.f"
	    dpsi += temp * temp;
#line 393 "dlasd4.f"
	    erretm += psi;
#line 394 "dlasd4.f"
/* L60: */
#line 394 "dlasd4.f"
	}
#line 395 "dlasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 399 "dlasd4.f"
	tau2 = work[*n] * delta[*n];
#line 400 "dlasd4.f"
	temp = z__[*n] / tau2;
#line 401 "dlasd4.f"
	phi = z__[*n] * temp;
#line 402 "dlasd4.f"
	dphi = temp * temp;
#line 403 "dlasd4.f"
	erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $          + ABS( TAU2 )*( DPSI+DPHI ) */

#line 406 "dlasd4.f"
	w = rhoinv + phi + psi;

/*        Main loop to update the values of the array   DELTA */

#line 410 "dlasd4.f"
	iter = niter + 1;

#line 412 "dlasd4.f"
	for (niter = iter; niter <= 400; ++niter) {

/*           Test for convergence */

#line 416 "dlasd4.f"
	    if (abs(w) <= eps * erretm) {
#line 417 "dlasd4.f"
		goto L240;
#line 418 "dlasd4.f"
	    }

/*           Calculate the new step */

#line 422 "dlasd4.f"
	    dtnsq1 = work[*n - 1] * delta[*n - 1];
#line 423 "dlasd4.f"
	    dtnsq = work[*n] * delta[*n];
#line 424 "dlasd4.f"
	    c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
#line 425 "dlasd4.f"
	    a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
#line 426 "dlasd4.f"
	    b = dtnsq1 * dtnsq * w;
#line 427 "dlasd4.f"
	    if (a >= 0.) {
#line 428 "dlasd4.f"
		eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 429 "dlasd4.f"
	    } else {
#line 430 "dlasd4.f"
		eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 431 "dlasd4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 439 "dlasd4.f"
	    if (w * eta > 0.) {
#line 439 "dlasd4.f"
		eta = -w / (dpsi + dphi);
#line 439 "dlasd4.f"
	    }
#line 441 "dlasd4.f"
	    temp = eta - dtnsq;
#line 442 "dlasd4.f"
	    if (temp <= 0.) {
#line 442 "dlasd4.f"
		eta /= 2.;
#line 442 "dlasd4.f"
	    }

#line 445 "dlasd4.f"
	    eta /= *sigma + sqrt(eta + *sigma * *sigma);
#line 446 "dlasd4.f"
	    tau += eta;
#line 447 "dlasd4.f"
	    *sigma += eta;

#line 449 "dlasd4.f"
	    i__1 = *n;
#line 449 "dlasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 450 "dlasd4.f"
		delta[j] -= eta;
#line 451 "dlasd4.f"
		work[j] += eta;
#line 452 "dlasd4.f"
/* L70: */
#line 452 "dlasd4.f"
	    }

/*           Evaluate PSI and the derivative DPSI */

#line 456 "dlasd4.f"
	    dpsi = 0.;
#line 457 "dlasd4.f"
	    psi = 0.;
#line 458 "dlasd4.f"
	    erretm = 0.;
#line 459 "dlasd4.f"
	    i__1 = ii;
#line 459 "dlasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 460 "dlasd4.f"
		temp = z__[j] / (work[j] * delta[j]);
#line 461 "dlasd4.f"
		psi += z__[j] * temp;
#line 462 "dlasd4.f"
		dpsi += temp * temp;
#line 463 "dlasd4.f"
		erretm += psi;
#line 464 "dlasd4.f"
/* L80: */
#line 464 "dlasd4.f"
	    }
#line 465 "dlasd4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 469 "dlasd4.f"
	    tau2 = work[*n] * delta[*n];
#line 470 "dlasd4.f"
	    temp = z__[*n] / tau2;
#line 471 "dlasd4.f"
	    phi = z__[*n] * temp;
#line 472 "dlasd4.f"
	    dphi = temp * temp;
#line 473 "dlasd4.f"
	    erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $             + ABS( TAU2 )*( DPSI+DPHI ) */

#line 476 "dlasd4.f"
	    w = rhoinv + phi + psi;
#line 477 "dlasd4.f"
/* L90: */
#line 477 "dlasd4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 481 "dlasd4.f"
	*info = 1;
#line 482 "dlasd4.f"
	goto L240;

/*        End for the case I = N */

#line 486 "dlasd4.f"
    } else {

/*        The case for I < N */

#line 490 "dlasd4.f"
	niter = 1;
#line 491 "dlasd4.f"
	ip1 = *i__ + 1;

/*        Calculate initial guess */

#line 495 "dlasd4.f"
	delsq = (d__[ip1] - d__[*i__]) * (d__[ip1] + d__[*i__]);
#line 496 "dlasd4.f"
	delsq2 = delsq / 2.;
#line 497 "dlasd4.f"
	sq2 = sqrt((d__[*i__] * d__[*i__] + d__[ip1] * d__[ip1]) / 2.);
#line 498 "dlasd4.f"
	temp = delsq2 / (d__[*i__] + sq2);
#line 499 "dlasd4.f"
	i__1 = *n;
#line 499 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 500 "dlasd4.f"
	    work[j] = d__[j] + d__[*i__] + temp;
#line 501 "dlasd4.f"
	    delta[j] = d__[j] - d__[*i__] - temp;
#line 502 "dlasd4.f"
/* L100: */
#line 502 "dlasd4.f"
	}

#line 504 "dlasd4.f"
	psi = 0.;
#line 505 "dlasd4.f"
	i__1 = *i__ - 1;
#line 505 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 506 "dlasd4.f"
	    psi += z__[j] * z__[j] / (work[j] * delta[j]);
#line 507 "dlasd4.f"
/* L110: */
#line 507 "dlasd4.f"
	}

#line 509 "dlasd4.f"
	phi = 0.;
#line 510 "dlasd4.f"
	i__1 = *i__ + 2;
#line 510 "dlasd4.f"
	for (j = *n; j >= i__1; --j) {
#line 511 "dlasd4.f"
	    phi += z__[j] * z__[j] / (work[j] * delta[j]);
#line 512 "dlasd4.f"
/* L120: */
#line 512 "dlasd4.f"
	}
#line 513 "dlasd4.f"
	c__ = rhoinv + psi + phi;
#line 514 "dlasd4.f"
	w = c__ + z__[*i__] * z__[*i__] / (work[*i__] * delta[*i__]) + z__[
		ip1] * z__[ip1] / (work[ip1] * delta[ip1]);

#line 517 "dlasd4.f"
	geomavg = FALSE_;
#line 518 "dlasd4.f"
	if (w > 0.) {

/*           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2 */

/*           We choose d(i) as origin. */

#line 524 "dlasd4.f"
	    orgati = TRUE_;
#line 525 "dlasd4.f"
	    ii = *i__;
#line 526 "dlasd4.f"
	    sglb = 0.;
#line 527 "dlasd4.f"
	    sgub = delsq2 / (d__[*i__] + sq2);
#line 528 "dlasd4.f"
	    a = c__ * delsq + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
#line 529 "dlasd4.f"
	    b = z__[*i__] * z__[*i__] * delsq;
#line 530 "dlasd4.f"
	    if (a > 0.) {
#line 531 "dlasd4.f"
		tau2 = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 532 "dlasd4.f"
	    } else {
#line 533 "dlasd4.f"
		tau2 = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / 
			(c__ * 2.);
#line 534 "dlasd4.f"
	    }

/*           TAU2 now is an estimation of SIGMA^2 - D( I )^2. The */
/*           following, however, is the corresponding estimation of */
/*           SIGMA - D( I ). */

#line 540 "dlasd4.f"
	    tau = tau2 / (d__[*i__] + sqrt(d__[*i__] * d__[*i__] + tau2));
#line 541 "dlasd4.f"
	    temp = sqrt(eps);
#line 542 "dlasd4.f"
	    if (d__[*i__] <= temp * d__[ip1] && (d__1 = z__[*i__], abs(d__1)) 
		    <= temp && d__[*i__] > 0.) {
/* Computing MIN */
#line 544 "dlasd4.f"
		d__1 = d__[*i__] * 10.;
#line 544 "dlasd4.f"
		tau = min(d__1,sgub);
#line 545 "dlasd4.f"
		geomavg = TRUE_;
#line 546 "dlasd4.f"
	    }
#line 547 "dlasd4.f"
	} else {

/*           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2 */

/*           We choose d(i+1) as origin. */

#line 553 "dlasd4.f"
	    orgati = FALSE_;
#line 554 "dlasd4.f"
	    ii = ip1;
#line 555 "dlasd4.f"
	    sglb = -delsq2 / (d__[ii] + sq2);
#line 556 "dlasd4.f"
	    sgub = 0.;
#line 557 "dlasd4.f"
	    a = c__ * delsq - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
#line 558 "dlasd4.f"
	    b = z__[ip1] * z__[ip1] * delsq;
#line 559 "dlasd4.f"
	    if (a < 0.) {
#line 560 "dlasd4.f"
		tau2 = b * 2. / (a - sqrt((d__1 = a * a + b * 4. * c__, abs(
			d__1))));
#line 561 "dlasd4.f"
	    } else {
#line 562 "dlasd4.f"
		tau2 = -(a + sqrt((d__1 = a * a + b * 4. * c__, abs(d__1)))) /
			 (c__ * 2.);
#line 563 "dlasd4.f"
	    }

/*           TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The */
/*           following, however, is the corresponding estimation of */
/*           SIGMA - D( IP1 ). */

#line 569 "dlasd4.f"
	    tau = tau2 / (d__[ip1] + sqrt((d__1 = d__[ip1] * d__[ip1] + tau2, 
		    abs(d__1))));
#line 571 "dlasd4.f"
	}

#line 573 "dlasd4.f"
	*sigma = d__[ii] + tau;
#line 574 "dlasd4.f"
	i__1 = *n;
#line 574 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 575 "dlasd4.f"
	    work[j] = d__[j] + d__[ii] + tau;
#line 576 "dlasd4.f"
	    delta[j] = d__[j] - d__[ii] - tau;
#line 577 "dlasd4.f"
/* L130: */
#line 577 "dlasd4.f"
	}
#line 578 "dlasd4.f"
	iim1 = ii - 1;
#line 579 "dlasd4.f"
	iip1 = ii + 1;

/*        Evaluate PSI and the derivative DPSI */

#line 583 "dlasd4.f"
	dpsi = 0.;
#line 584 "dlasd4.f"
	psi = 0.;
#line 585 "dlasd4.f"
	erretm = 0.;
#line 586 "dlasd4.f"
	i__1 = iim1;
#line 586 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 587 "dlasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 588 "dlasd4.f"
	    psi += z__[j] * temp;
#line 589 "dlasd4.f"
	    dpsi += temp * temp;
#line 590 "dlasd4.f"
	    erretm += psi;
#line 591 "dlasd4.f"
/* L150: */
#line 591 "dlasd4.f"
	}
#line 592 "dlasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 596 "dlasd4.f"
	dphi = 0.;
#line 597 "dlasd4.f"
	phi = 0.;
#line 598 "dlasd4.f"
	i__1 = iip1;
#line 598 "dlasd4.f"
	for (j = *n; j >= i__1; --j) {
#line 599 "dlasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 600 "dlasd4.f"
	    phi += z__[j] * temp;
#line 601 "dlasd4.f"
	    dphi += temp * temp;
#line 602 "dlasd4.f"
	    erretm += phi;
#line 603 "dlasd4.f"
/* L160: */
#line 603 "dlasd4.f"
	}

#line 605 "dlasd4.f"
	w = rhoinv + phi + psi;

/*        W is the value of the secular function with */
/*        its ii-th element removed. */

#line 610 "dlasd4.f"
	swtch3 = FALSE_;
#line 611 "dlasd4.f"
	if (orgati) {
#line 612 "dlasd4.f"
	    if (w < 0.) {
#line 612 "dlasd4.f"
		swtch3 = TRUE_;
#line 612 "dlasd4.f"
	    }
#line 614 "dlasd4.f"
	} else {
#line 615 "dlasd4.f"
	    if (w > 0.) {
#line 615 "dlasd4.f"
		swtch3 = TRUE_;
#line 615 "dlasd4.f"
	    }
#line 617 "dlasd4.f"
	}
#line 618 "dlasd4.f"
	if (ii == 1 || ii == *n) {
#line 618 "dlasd4.f"
	    swtch3 = FALSE_;
#line 618 "dlasd4.f"
	}

#line 621 "dlasd4.f"
	temp = z__[ii] / (work[ii] * delta[ii]);
#line 622 "dlasd4.f"
	dw = dpsi + dphi + temp * temp;
#line 623 "dlasd4.f"
	temp = z__[ii] * temp;
#line 624 "dlasd4.f"
	w += temp;
#line 625 "dlasd4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $          + ABS( TAU2 )*DW */

/*        Test for convergence */

#line 631 "dlasd4.f"
	if (abs(w) <= eps * erretm) {
#line 632 "dlasd4.f"
	    goto L240;
#line 633 "dlasd4.f"
	}

#line 635 "dlasd4.f"
	if (w <= 0.) {
#line 636 "dlasd4.f"
	    sglb = max(sglb,tau);
#line 637 "dlasd4.f"
	} else {
#line 638 "dlasd4.f"
	    sgub = min(sgub,tau);
#line 639 "dlasd4.f"
	}

/*        Calculate the new step */

#line 643 "dlasd4.f"
	++niter;
#line 644 "dlasd4.f"
	if (! swtch3) {
#line 645 "dlasd4.f"
	    dtipsq = work[ip1] * delta[ip1];
#line 646 "dlasd4.f"
	    dtisq = work[*i__] * delta[*i__];
#line 647 "dlasd4.f"
	    if (orgati) {
/* Computing 2nd power */
#line 648 "dlasd4.f"
		d__1 = z__[*i__] / dtisq;
#line 648 "dlasd4.f"
		c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 649 "dlasd4.f"
	    } else {
/* Computing 2nd power */
#line 650 "dlasd4.f"
		d__1 = z__[ip1] / dtipsq;
#line 650 "dlasd4.f"
		c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 651 "dlasd4.f"
	    }
#line 652 "dlasd4.f"
	    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 653 "dlasd4.f"
	    b = dtipsq * dtisq * w;
#line 654 "dlasd4.f"
	    if (c__ == 0.) {
#line 655 "dlasd4.f"
		if (a == 0.) {
#line 656 "dlasd4.f"
		    if (orgati) {
#line 657 "dlasd4.f"
			a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + 
				dphi);
#line 658 "dlasd4.f"
		    } else {
#line 659 "dlasd4.f"
			a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + 
				dphi);
#line 660 "dlasd4.f"
		    }
#line 661 "dlasd4.f"
		}
#line 662 "dlasd4.f"
		eta = b / a;
#line 663 "dlasd4.f"
	    } else if (a <= 0.) {
#line 664 "dlasd4.f"
		eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
			c__ * 2.);
#line 665 "dlasd4.f"
	    } else {
#line 666 "dlasd4.f"
		eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
			d__1))));
#line 667 "dlasd4.f"
	    }
#line 668 "dlasd4.f"
	} else {

/*           Interpolation using THREE most relevant poles */

#line 672 "dlasd4.f"
	    dtiim = work[iim1] * delta[iim1];
#line 673 "dlasd4.f"
	    dtiip = work[iip1] * delta[iip1];
#line 674 "dlasd4.f"
	    temp = rhoinv + psi + phi;
#line 675 "dlasd4.f"
	    if (orgati) {
#line 676 "dlasd4.f"
		temp1 = z__[iim1] / dtiim;
#line 677 "dlasd4.f"
		temp1 *= temp1;
#line 678 "dlasd4.f"
		c__ = temp - dtiip * (dpsi + dphi) - (d__[iim1] - d__[iip1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
#line 680 "dlasd4.f"
		zz[0] = z__[iim1] * z__[iim1];
#line 681 "dlasd4.f"
		if (dpsi < temp1) {
#line 682 "dlasd4.f"
		    zz[2] = dtiip * dtiip * dphi;
#line 683 "dlasd4.f"
		} else {
#line 684 "dlasd4.f"
		    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
#line 685 "dlasd4.f"
		}
#line 686 "dlasd4.f"
	    } else {
#line 687 "dlasd4.f"
		temp1 = z__[iip1] / dtiip;
#line 688 "dlasd4.f"
		temp1 *= temp1;
#line 689 "dlasd4.f"
		c__ = temp - dtiim * (dpsi + dphi) - (d__[iip1] - d__[iim1]) *
			 (d__[iim1] + d__[iip1]) * temp1;
#line 691 "dlasd4.f"
		if (dphi < temp1) {
#line 692 "dlasd4.f"
		    zz[0] = dtiim * dtiim * dpsi;
#line 693 "dlasd4.f"
		} else {
#line 694 "dlasd4.f"
		    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
#line 695 "dlasd4.f"
		}
#line 696 "dlasd4.f"
		zz[2] = z__[iip1] * z__[iip1];
#line 697 "dlasd4.f"
	    }
#line 698 "dlasd4.f"
	    zz[1] = z__[ii] * z__[ii];
#line 699 "dlasd4.f"
	    dd[0] = dtiim;
#line 700 "dlasd4.f"
	    dd[1] = delta[ii] * work[ii];
#line 701 "dlasd4.f"
	    dd[2] = dtiip;
#line 702 "dlasd4.f"
	    dlaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);

#line 704 "dlasd4.f"
	    if (*info != 0) {

/*              If INFO is not 0, i.e., DLAED6 failed, switch back */
/*              to 2 pole interpolation. */

#line 709 "dlasd4.f"
		swtch3 = FALSE_;
#line 710 "dlasd4.f"
		*info = 0;
#line 711 "dlasd4.f"
		dtipsq = work[ip1] * delta[ip1];
#line 712 "dlasd4.f"
		dtisq = work[*i__] * delta[*i__];
#line 713 "dlasd4.f"
		if (orgati) {
/* Computing 2nd power */
#line 714 "dlasd4.f"
		    d__1 = z__[*i__] / dtisq;
#line 714 "dlasd4.f"
		    c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 715 "dlasd4.f"
		} else {
/* Computing 2nd power */
#line 716 "dlasd4.f"
		    d__1 = z__[ip1] / dtipsq;
#line 716 "dlasd4.f"
		    c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 717 "dlasd4.f"
		}
#line 718 "dlasd4.f"
		a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 719 "dlasd4.f"
		b = dtipsq * dtisq * w;
#line 720 "dlasd4.f"
		if (c__ == 0.) {
#line 721 "dlasd4.f"
		    if (a == 0.) {
#line 722 "dlasd4.f"
			if (orgati) {
#line 723 "dlasd4.f"
			    a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (
				    dpsi + dphi);
#line 724 "dlasd4.f"
			} else {
#line 725 "dlasd4.f"
			    a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + 
				    dphi);
#line 726 "dlasd4.f"
			}
#line 727 "dlasd4.f"
		    }
#line 728 "dlasd4.f"
		    eta = b / a;
#line 729 "dlasd4.f"
		} else if (a <= 0.) {
#line 730 "dlasd4.f"
		    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
			     / (c__ * 2.);
#line 731 "dlasd4.f"
		} else {
#line 732 "dlasd4.f"
		    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, 
			    abs(d__1))));
#line 733 "dlasd4.f"
		}
#line 734 "dlasd4.f"
	    }
#line 735 "dlasd4.f"
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

#line 743 "dlasd4.f"
	if (w * eta >= 0.) {
#line 743 "dlasd4.f"
	    eta = -w / dw;
#line 743 "dlasd4.f"
	}

#line 746 "dlasd4.f"
	eta /= *sigma + sqrt(*sigma * *sigma + eta);
#line 747 "dlasd4.f"
	temp = tau + eta;
#line 748 "dlasd4.f"
	if (temp > sgub || temp < sglb) {
#line 749 "dlasd4.f"
	    if (w < 0.) {
#line 750 "dlasd4.f"
		eta = (sgub - tau) / 2.;
#line 751 "dlasd4.f"
	    } else {
#line 752 "dlasd4.f"
		eta = (sglb - tau) / 2.;
#line 753 "dlasd4.f"
	    }
#line 754 "dlasd4.f"
	    if (geomavg) {
#line 755 "dlasd4.f"
		if (w < 0.) {
#line 756 "dlasd4.f"
		    if (tau > 0.) {
#line 757 "dlasd4.f"
			eta = sqrt(sgub * tau) - tau;
#line 758 "dlasd4.f"
		    }
#line 759 "dlasd4.f"
		} else {
#line 760 "dlasd4.f"
		    if (sglb > 0.) {
#line 761 "dlasd4.f"
			eta = sqrt(sglb * tau) - tau;
#line 762 "dlasd4.f"
		    }
#line 763 "dlasd4.f"
		}
#line 764 "dlasd4.f"
	    }
#line 765 "dlasd4.f"
	}

#line 767 "dlasd4.f"
	prew = w;

#line 769 "dlasd4.f"
	tau += eta;
#line 770 "dlasd4.f"
	*sigma += eta;

#line 772 "dlasd4.f"
	i__1 = *n;
#line 772 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 773 "dlasd4.f"
	    work[j] += eta;
#line 774 "dlasd4.f"
	    delta[j] -= eta;
#line 775 "dlasd4.f"
/* L170: */
#line 775 "dlasd4.f"
	}

/*        Evaluate PSI and the derivative DPSI */

#line 779 "dlasd4.f"
	dpsi = 0.;
#line 780 "dlasd4.f"
	psi = 0.;
#line 781 "dlasd4.f"
	erretm = 0.;
#line 782 "dlasd4.f"
	i__1 = iim1;
#line 782 "dlasd4.f"
	for (j = 1; j <= i__1; ++j) {
#line 783 "dlasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 784 "dlasd4.f"
	    psi += z__[j] * temp;
#line 785 "dlasd4.f"
	    dpsi += temp * temp;
#line 786 "dlasd4.f"
	    erretm += psi;
#line 787 "dlasd4.f"
/* L180: */
#line 787 "dlasd4.f"
	}
#line 788 "dlasd4.f"
	erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

#line 792 "dlasd4.f"
	dphi = 0.;
#line 793 "dlasd4.f"
	phi = 0.;
#line 794 "dlasd4.f"
	i__1 = iip1;
#line 794 "dlasd4.f"
	for (j = *n; j >= i__1; --j) {
#line 795 "dlasd4.f"
	    temp = z__[j] / (work[j] * delta[j]);
#line 796 "dlasd4.f"
	    phi += z__[j] * temp;
#line 797 "dlasd4.f"
	    dphi += temp * temp;
#line 798 "dlasd4.f"
	    erretm += phi;
#line 799 "dlasd4.f"
/* L190: */
#line 799 "dlasd4.f"
	}

#line 801 "dlasd4.f"
	tau2 = work[ii] * delta[ii];
#line 802 "dlasd4.f"
	temp = z__[ii] / tau2;
#line 803 "dlasd4.f"
	dw = dpsi + dphi + temp * temp;
#line 804 "dlasd4.f"
	temp = z__[ii] * temp;
#line 805 "dlasd4.f"
	w = rhoinv + phi + psi + temp;
#line 806 "dlasd4.f"
	erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $          + ABS( TAU2 )*DW */

#line 810 "dlasd4.f"
	swtch = FALSE_;
#line 811 "dlasd4.f"
	if (orgati) {
#line 812 "dlasd4.f"
	    if (-w > abs(prew) / 10.) {
#line 812 "dlasd4.f"
		swtch = TRUE_;
#line 812 "dlasd4.f"
	    }
#line 814 "dlasd4.f"
	} else {
#line 815 "dlasd4.f"
	    if (w > abs(prew) / 10.) {
#line 815 "dlasd4.f"
		swtch = TRUE_;
#line 815 "dlasd4.f"
	    }
#line 817 "dlasd4.f"
	}

/*        Main loop to update the values of the array   DELTA and WORK */

#line 821 "dlasd4.f"
	iter = niter + 1;

#line 823 "dlasd4.f"
	for (niter = iter; niter <= 400; ++niter) {

/*           Test for convergence */

#line 827 "dlasd4.f"
	    if (abs(w) <= eps * erretm) {
/*     $          .OR. (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN */
#line 829 "dlasd4.f"
		goto L240;
#line 830 "dlasd4.f"
	    }

#line 832 "dlasd4.f"
	    if (w <= 0.) {
#line 833 "dlasd4.f"
		sglb = max(sglb,tau);
#line 834 "dlasd4.f"
	    } else {
#line 835 "dlasd4.f"
		sgub = min(sgub,tau);
#line 836 "dlasd4.f"
	    }

/*           Calculate the new step */

#line 840 "dlasd4.f"
	    if (! swtch3) {
#line 841 "dlasd4.f"
		dtipsq = work[ip1] * delta[ip1];
#line 842 "dlasd4.f"
		dtisq = work[*i__] * delta[*i__];
#line 843 "dlasd4.f"
		if (! swtch) {
#line 844 "dlasd4.f"
		    if (orgati) {
/* Computing 2nd power */
#line 845 "dlasd4.f"
			d__1 = z__[*i__] / dtisq;
#line 845 "dlasd4.f"
			c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 846 "dlasd4.f"
		    } else {
/* Computing 2nd power */
#line 847 "dlasd4.f"
			d__1 = z__[ip1] / dtipsq;
#line 847 "dlasd4.f"
			c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 848 "dlasd4.f"
		    }
#line 849 "dlasd4.f"
		} else {
#line 850 "dlasd4.f"
		    temp = z__[ii] / (work[ii] * delta[ii]);
#line 851 "dlasd4.f"
		    if (orgati) {
#line 852 "dlasd4.f"
			dpsi += temp * temp;
#line 853 "dlasd4.f"
		    } else {
#line 854 "dlasd4.f"
			dphi += temp * temp;
#line 855 "dlasd4.f"
		    }
#line 856 "dlasd4.f"
		    c__ = w - dtisq * dpsi - dtipsq * dphi;
#line 857 "dlasd4.f"
		}
#line 858 "dlasd4.f"
		a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 859 "dlasd4.f"
		b = dtipsq * dtisq * w;
#line 860 "dlasd4.f"
		if (c__ == 0.) {
#line 861 "dlasd4.f"
		    if (a == 0.) {
#line 862 "dlasd4.f"
			if (! swtch) {
#line 863 "dlasd4.f"
			    if (orgati) {
#line 864 "dlasd4.f"
				a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * 
					(dpsi + dphi);
#line 866 "dlasd4.f"
			    } else {
#line 867 "dlasd4.f"
				a = z__[ip1] * z__[ip1] + dtisq * dtisq * (
					dpsi + dphi);
#line 869 "dlasd4.f"
			    }
#line 870 "dlasd4.f"
			} else {
#line 871 "dlasd4.f"
			    a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
#line 872 "dlasd4.f"
			}
#line 873 "dlasd4.f"
		    }
#line 874 "dlasd4.f"
		    eta = b / a;
#line 875 "dlasd4.f"
		} else if (a <= 0.) {
#line 876 "dlasd4.f"
		    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
			     / (c__ * 2.);
#line 877 "dlasd4.f"
		} else {
#line 878 "dlasd4.f"
		    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, 
			    abs(d__1))));
#line 879 "dlasd4.f"
		}
#line 880 "dlasd4.f"
	    } else {

/*              Interpolation using THREE most relevant poles */

#line 884 "dlasd4.f"
		dtiim = work[iim1] * delta[iim1];
#line 885 "dlasd4.f"
		dtiip = work[iip1] * delta[iip1];
#line 886 "dlasd4.f"
		temp = rhoinv + psi + phi;
#line 887 "dlasd4.f"
		if (swtch) {
#line 888 "dlasd4.f"
		    c__ = temp - dtiim * dpsi - dtiip * dphi;
#line 889 "dlasd4.f"
		    zz[0] = dtiim * dtiim * dpsi;
#line 890 "dlasd4.f"
		    zz[2] = dtiip * dtiip * dphi;
#line 891 "dlasd4.f"
		} else {
#line 892 "dlasd4.f"
		    if (orgati) {
#line 893 "dlasd4.f"
			temp1 = z__[iim1] / dtiim;
#line 894 "dlasd4.f"
			temp1 *= temp1;
#line 895 "dlasd4.f"
			temp2 = (d__[iim1] - d__[iip1]) * (d__[iim1] + d__[
				iip1]) * temp1;
#line 897 "dlasd4.f"
			c__ = temp - dtiip * (dpsi + dphi) - temp2;
#line 898 "dlasd4.f"
			zz[0] = z__[iim1] * z__[iim1];
#line 899 "dlasd4.f"
			if (dpsi < temp1) {
#line 900 "dlasd4.f"
			    zz[2] = dtiip * dtiip * dphi;
#line 901 "dlasd4.f"
			} else {
#line 902 "dlasd4.f"
			    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
#line 903 "dlasd4.f"
			}
#line 904 "dlasd4.f"
		    } else {
#line 905 "dlasd4.f"
			temp1 = z__[iip1] / dtiip;
#line 906 "dlasd4.f"
			temp1 *= temp1;
#line 907 "dlasd4.f"
			temp2 = (d__[iip1] - d__[iim1]) * (d__[iim1] + d__[
				iip1]) * temp1;
#line 909 "dlasd4.f"
			c__ = temp - dtiim * (dpsi + dphi) - temp2;
#line 910 "dlasd4.f"
			if (dphi < temp1) {
#line 911 "dlasd4.f"
			    zz[0] = dtiim * dtiim * dpsi;
#line 912 "dlasd4.f"
			} else {
#line 913 "dlasd4.f"
			    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
#line 914 "dlasd4.f"
			}
#line 915 "dlasd4.f"
			zz[2] = z__[iip1] * z__[iip1];
#line 916 "dlasd4.f"
		    }
#line 917 "dlasd4.f"
		}
#line 918 "dlasd4.f"
		dd[0] = dtiim;
#line 919 "dlasd4.f"
		dd[1] = delta[ii] * work[ii];
#line 920 "dlasd4.f"
		dd[2] = dtiip;
#line 921 "dlasd4.f"
		dlaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);

#line 923 "dlasd4.f"
		if (*info != 0) {

/*                 If INFO is not 0, i.e., DLAED6 failed, switch */
/*                 back to two pole interpolation */

#line 928 "dlasd4.f"
		    swtch3 = FALSE_;
#line 929 "dlasd4.f"
		    *info = 0;
#line 930 "dlasd4.f"
		    dtipsq = work[ip1] * delta[ip1];
#line 931 "dlasd4.f"
		    dtisq = work[*i__] * delta[*i__];
#line 932 "dlasd4.f"
		    if (! swtch) {
#line 933 "dlasd4.f"
			if (orgati) {
/* Computing 2nd power */
#line 934 "dlasd4.f"
			    d__1 = z__[*i__] / dtisq;
#line 934 "dlasd4.f"
			    c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
#line 935 "dlasd4.f"
			} else {
/* Computing 2nd power */
#line 936 "dlasd4.f"
			    d__1 = z__[ip1] / dtipsq;
#line 936 "dlasd4.f"
			    c__ = w - dtisq * dw - delsq * (d__1 * d__1);
#line 937 "dlasd4.f"
			}
#line 938 "dlasd4.f"
		    } else {
#line 939 "dlasd4.f"
			temp = z__[ii] / (work[ii] * delta[ii]);
#line 940 "dlasd4.f"
			if (orgati) {
#line 941 "dlasd4.f"
			    dpsi += temp * temp;
#line 942 "dlasd4.f"
			} else {
#line 943 "dlasd4.f"
			    dphi += temp * temp;
#line 944 "dlasd4.f"
			}
#line 945 "dlasd4.f"
			c__ = w - dtisq * dpsi - dtipsq * dphi;
#line 946 "dlasd4.f"
		    }
#line 947 "dlasd4.f"
		    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
#line 948 "dlasd4.f"
		    b = dtipsq * dtisq * w;
#line 949 "dlasd4.f"
		    if (c__ == 0.) {
#line 950 "dlasd4.f"
			if (a == 0.) {
#line 951 "dlasd4.f"
			    if (! swtch) {
#line 952 "dlasd4.f"
				if (orgati) {
#line 953 "dlasd4.f"
				    a = z__[*i__] * z__[*i__] + dtipsq * 
					    dtipsq * (dpsi + dphi);
#line 955 "dlasd4.f"
				} else {
#line 956 "dlasd4.f"
				    a = z__[ip1] * z__[ip1] + dtisq * dtisq * 
					    (dpsi + dphi);
#line 958 "dlasd4.f"
				}
#line 959 "dlasd4.f"
			    } else {
#line 960 "dlasd4.f"
				a = dtisq * dtisq * dpsi + dtipsq * dtipsq * 
					dphi;
#line 961 "dlasd4.f"
			    }
#line 962 "dlasd4.f"
			}
#line 963 "dlasd4.f"
			eta = b / a;
#line 964 "dlasd4.f"
		    } else if (a <= 0.) {
#line 965 "dlasd4.f"
			eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
				d__1)))) / (c__ * 2.);
#line 966 "dlasd4.f"
		    } else {
#line 967 "dlasd4.f"
			eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__,
				 abs(d__1))));
#line 968 "dlasd4.f"
		    }
#line 969 "dlasd4.f"
		}
#line 970 "dlasd4.f"
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

#line 978 "dlasd4.f"
	    if (w * eta >= 0.) {
#line 978 "dlasd4.f"
		eta = -w / dw;
#line 978 "dlasd4.f"
	    }

#line 981 "dlasd4.f"
	    eta /= *sigma + sqrt(*sigma * *sigma + eta);
#line 982 "dlasd4.f"
	    temp = tau + eta;
#line 983 "dlasd4.f"
	    if (temp > sgub || temp < sglb) {
#line 984 "dlasd4.f"
		if (w < 0.) {
#line 985 "dlasd4.f"
		    eta = (sgub - tau) / 2.;
#line 986 "dlasd4.f"
		} else {
#line 987 "dlasd4.f"
		    eta = (sglb - tau) / 2.;
#line 988 "dlasd4.f"
		}
#line 989 "dlasd4.f"
		if (geomavg) {
#line 990 "dlasd4.f"
		    if (w < 0.) {
#line 991 "dlasd4.f"
			if (tau > 0.) {
#line 992 "dlasd4.f"
			    eta = sqrt(sgub * tau) - tau;
#line 993 "dlasd4.f"
			}
#line 994 "dlasd4.f"
		    } else {
#line 995 "dlasd4.f"
			if (sglb > 0.) {
#line 996 "dlasd4.f"
			    eta = sqrt(sglb * tau) - tau;
#line 997 "dlasd4.f"
			}
#line 998 "dlasd4.f"
		    }
#line 999 "dlasd4.f"
		}
#line 1000 "dlasd4.f"
	    }

#line 1002 "dlasd4.f"
	    prew = w;

#line 1004 "dlasd4.f"
	    tau += eta;
#line 1005 "dlasd4.f"
	    *sigma += eta;

#line 1007 "dlasd4.f"
	    i__1 = *n;
#line 1007 "dlasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1008 "dlasd4.f"
		work[j] += eta;
#line 1009 "dlasd4.f"
		delta[j] -= eta;
#line 1010 "dlasd4.f"
/* L200: */
#line 1010 "dlasd4.f"
	    }

/*           Evaluate PSI and the derivative DPSI */

#line 1014 "dlasd4.f"
	    dpsi = 0.;
#line 1015 "dlasd4.f"
	    psi = 0.;
#line 1016 "dlasd4.f"
	    erretm = 0.;
#line 1017 "dlasd4.f"
	    i__1 = iim1;
#line 1017 "dlasd4.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1018 "dlasd4.f"
		temp = z__[j] / (work[j] * delta[j]);
#line 1019 "dlasd4.f"
		psi += z__[j] * temp;
#line 1020 "dlasd4.f"
		dpsi += temp * temp;
#line 1021 "dlasd4.f"
		erretm += psi;
#line 1022 "dlasd4.f"
/* L210: */
#line 1022 "dlasd4.f"
	    }
#line 1023 "dlasd4.f"
	    erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

#line 1027 "dlasd4.f"
	    dphi = 0.;
#line 1028 "dlasd4.f"
	    phi = 0.;
#line 1029 "dlasd4.f"
	    i__1 = iip1;
#line 1029 "dlasd4.f"
	    for (j = *n; j >= i__1; --j) {
#line 1030 "dlasd4.f"
		temp = z__[j] / (work[j] * delta[j]);
#line 1031 "dlasd4.f"
		phi += z__[j] * temp;
#line 1032 "dlasd4.f"
		dphi += temp * temp;
#line 1033 "dlasd4.f"
		erretm += phi;
#line 1034 "dlasd4.f"
/* L220: */
#line 1034 "dlasd4.f"
	    }

#line 1036 "dlasd4.f"
	    tau2 = work[ii] * delta[ii];
#line 1037 "dlasd4.f"
	    temp = z__[ii] / tau2;
#line 1038 "dlasd4.f"
	    dw = dpsi + dphi + temp * temp;
#line 1039 "dlasd4.f"
	    temp = z__[ii] * temp;
#line 1040 "dlasd4.f"
	    w = rhoinv + phi + psi + temp;
#line 1041 "dlasd4.f"
	    erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $             + ABS( TAU2 )*DW */

#line 1045 "dlasd4.f"
	    if (w * prew > 0. && abs(w) > abs(prew) / 10.) {
#line 1045 "dlasd4.f"
		swtch = ! swtch;
#line 1045 "dlasd4.f"
	    }

#line 1048 "dlasd4.f"
/* L230: */
#line 1048 "dlasd4.f"
	}

/*        Return with INFO = 1, NITER = MAXIT and not converged */

#line 1052 "dlasd4.f"
	*info = 1;

#line 1054 "dlasd4.f"
    }

#line 1056 "dlasd4.f"
L240:
#line 1057 "dlasd4.f"
    return 0;

/*     End of DLASD4 */

} /* dlasd4_ */

