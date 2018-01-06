#line 1 "slar1v.f"
/* slar1v.f -- translated by f2c (version 20100827).
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

#line 1 "slar1v.f"
/* > \brief \b SLAR1V computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn 
of the tridiagonal matrix LDLT - λI. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAR1V + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slar1v.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slar1v.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slar1v.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD, */
/*                  PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, */
/*                  R, ISUPPZ, NRMINV, RESID, RQCORR, WORK ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTNC */
/*       INTEGER   B1, BN, N, NEGCNT, R */
/*       REAL               GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID, */
/*      $                   RQCORR, ZTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ) */
/*       REAL               D( * ), L( * ), LD( * ), LLD( * ), */
/*      $                  WORK( * ) */
/*       REAL             Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAR1V computes the (scaled) r-th column of the inverse of */
/* > the sumbmatrix in rows B1 through BN of the tridiagonal matrix */
/* > L D L**T - sigma I. When sigma is close to an eigenvalue, the */
/* > computed vector is an accurate eigenvector. Usually, r corresponds */
/* > to the index where the eigenvector is largest in magnitude. */
/* > The following steps accomplish this computation : */
/* > (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T, */
/* > (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T, */
/* > (c) Computation of the diagonal elements of the inverse of */
/* >     L D L**T - sigma I by combining the above transforms, and choosing */
/* >     r as the index where the diagonal of the inverse is (one of the) */
/* >     largest in magnitude. */
/* > (d) Computation of the (scaled) r-th column of the inverse using the */
/* >     twisted factorization obtained by combining the top part of the */
/* >     the stationary and the bottom part of the progressive transform. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           The order of the matrix L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] B1 */
/* > \verbatim */
/* >          B1 is INTEGER */
/* >           First index of the submatrix of L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] BN */
/* > \verbatim */
/* >          BN is INTEGER */
/* >           Last index of the submatrix of L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LAMBDA */
/* > \verbatim */
/* >          LAMBDA is REAL */
/* >           The shift. In order to compute an accurate eigenvector, */
/* >           LAMBDA should be a good approximation to an eigenvalue */
/* >           of L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is REAL array, dimension (N-1) */
/* >           The (n-1) subdiagonal elements of the unit bidiagonal matrix */
/* >           L, in elements 1 to N-1. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >           The n diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LD */
/* > \verbatim */
/* >          LD is REAL array, dimension (N-1) */
/* >           The n-1 elements L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] LLD */
/* > \verbatim */
/* >          LLD is REAL array, dimension (N-1) */
/* >           The n-1 elements L(i)*L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
/* >           The minimum pivot in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[in] GAPTOL */
/* > \verbatim */
/* >          GAPTOL is REAL */
/* >           Tolerance that indicates when eigenvector entries are negligible */
/* >           w.r.t. their contribution to the residual. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (N) */
/* >           On input, all entries of Z must be set to 0. */
/* >           On output, Z contains the (scaled) r-th column of the */
/* >           inverse. The scaling is such that Z(R) equals 1. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTNC */
/* > \verbatim */
/* >          WANTNC is LOGICAL */
/* >           Specifies whether NEGCNT has to be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] NEGCNT */
/* > \verbatim */
/* >          NEGCNT is INTEGER */
/* >           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin */
/* >           in the  matrix factorization L D L**T, and NEGCNT = -1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] ZTZ */
/* > \verbatim */
/* >          ZTZ is REAL */
/* >           The square of the 2-norm of Z. */
/* > \endverbatim */
/* > */
/* > \param[out] MINGMA */
/* > \verbatim */
/* >          MINGMA is REAL */
/* >           The reciprocal of the largest (in magnitude) diagonal */
/* >           element of the inverse of L D L**T - sigma I. */
/* > \endverbatim */
/* > */
/* > \param[in,out] R */
/* > \verbatim */
/* >          R is INTEGER */
/* >           The twist index for the twisted factorization used to */
/* >           compute Z. */
/* >           On input, 0 <= R <= N. If R is input as 0, R is set to */
/* >           the index where (L D L**T - sigma I)^{-1} is largest */
/* >           in magnitude. If 1 <= R <= N, R is unchanged. */
/* >           On output, R contains the twist index used to compute Z. */
/* >           Ideally, R designates the position of the maximum entry in the */
/* >           eigenvector. */
/* > \endverbatim */
/* > */
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* >          ISUPPZ is INTEGER array, dimension (2) */
/* >           The support of the vector in Z, i.e., the vector Z is */
/* >           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ). */
/* > \endverbatim */
/* > */
/* > \param[out] NRMINV */
/* > \verbatim */
/* >          NRMINV is REAL */
/* >           NRMINV = 1/SQRT( ZTZ ) */
/* > \endverbatim */
/* > */
/* > \param[out] RESID */
/* > \verbatim */
/* >          RESID is REAL */
/* >           The residual of the FP vector. */
/* >           RESID = ABS( MINGMA )/SQRT( ZTZ ) */
/* > \endverbatim */
/* > */
/* > \param[out] RQCORR */
/* > \verbatim */
/* >          RQCORR is REAL */
/* >           The Rayleigh Quotient correction to LAMBDA. */
/* >           RQCORR = MINGMA*TMP */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (4*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int slar1v_(integer *n, integer *b1, integer *bn, doublereal 
	*lambda, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *pivmin, doublereal *gaptol, doublereal *z__, logical 
	*wantnc, integer *negcnt, doublereal *ztz, doublereal *mingma, 
	integer *r__, integer *isuppz, doublereal *nrminv, doublereal *resid, 
	doublereal *rqcorr, doublereal *work)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s;
    static integer r1, r2;
    static doublereal eps, tmp;
    static integer neg1, neg2, indp, inds;
    static doublereal dplus;
    extern doublereal slamch_(char *, ftnlen);
    static integer indlpl, indumn;
    extern logical sisnan_(doublereal *);
    static doublereal dminus;
    static logical sawnan1, sawnan2;


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
/*     .. Executable Statements .. */

#line 274 "slar1v.f"
    /* Parameter adjustments */
#line 274 "slar1v.f"
    --work;
#line 274 "slar1v.f"
    --isuppz;
#line 274 "slar1v.f"
    --z__;
#line 274 "slar1v.f"
    --lld;
#line 274 "slar1v.f"
    --ld;
#line 274 "slar1v.f"
    --l;
#line 274 "slar1v.f"
    --d__;
#line 274 "slar1v.f"

#line 274 "slar1v.f"
    /* Function Body */
#line 274 "slar1v.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 277 "slar1v.f"
    if (*r__ == 0) {
#line 278 "slar1v.f"
	r1 = *b1;
#line 279 "slar1v.f"
	r2 = *bn;
#line 280 "slar1v.f"
    } else {
#line 281 "slar1v.f"
	r1 = *r__;
#line 282 "slar1v.f"
	r2 = *r__;
#line 283 "slar1v.f"
    }
/*     Storage for LPLUS */
#line 286 "slar1v.f"
    indlpl = 0;
/*     Storage for UMINUS */
#line 288 "slar1v.f"
    indumn = *n;
#line 289 "slar1v.f"
    inds = (*n << 1) + 1;
#line 290 "slar1v.f"
    indp = *n * 3 + 1;
#line 292 "slar1v.f"
    if (*b1 == 1) {
#line 293 "slar1v.f"
	work[inds] = 0.;
#line 294 "slar1v.f"
    } else {
#line 295 "slar1v.f"
	work[inds + *b1 - 1] = lld[*b1 - 1];
#line 296 "slar1v.f"
    }

/*     Compute the stationary transform (using the differential form) */
/*     until the index R2. */

#line 302 "slar1v.f"
    sawnan1 = FALSE_;
#line 303 "slar1v.f"
    neg1 = 0;
#line 304 "slar1v.f"
    s = work[inds + *b1 - 1] - *lambda;
#line 305 "slar1v.f"
    i__1 = r1 - 1;
#line 305 "slar1v.f"
    for (i__ = *b1; i__ <= i__1; ++i__) {
#line 306 "slar1v.f"
	dplus = d__[i__] + s;
#line 307 "slar1v.f"
	work[indlpl + i__] = ld[i__] / dplus;
#line 308 "slar1v.f"
	if (dplus < 0.) {
#line 308 "slar1v.f"
	    ++neg1;
#line 308 "slar1v.f"
	}
#line 309 "slar1v.f"
	work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 310 "slar1v.f"
	s = work[inds + i__] - *lambda;
#line 311 "slar1v.f"
/* L50: */
#line 311 "slar1v.f"
    }
#line 312 "slar1v.f"
    sawnan1 = sisnan_(&s);
#line 313 "slar1v.f"
    if (sawnan1) {
#line 313 "slar1v.f"
	goto L60;
#line 313 "slar1v.f"
    }
#line 314 "slar1v.f"
    i__1 = r2 - 1;
#line 314 "slar1v.f"
    for (i__ = r1; i__ <= i__1; ++i__) {
#line 315 "slar1v.f"
	dplus = d__[i__] + s;
#line 316 "slar1v.f"
	work[indlpl + i__] = ld[i__] / dplus;
#line 317 "slar1v.f"
	work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 318 "slar1v.f"
	s = work[inds + i__] - *lambda;
#line 319 "slar1v.f"
/* L51: */
#line 319 "slar1v.f"
    }
#line 320 "slar1v.f"
    sawnan1 = sisnan_(&s);

#line 322 "slar1v.f"
L60:
#line 323 "slar1v.f"
    if (sawnan1) {
/*        Runs a slower version of the above loop if a NaN is detected */
#line 325 "slar1v.f"
	neg1 = 0;
#line 326 "slar1v.f"
	s = work[inds + *b1 - 1] - *lambda;
#line 327 "slar1v.f"
	i__1 = r1 - 1;
#line 327 "slar1v.f"
	for (i__ = *b1; i__ <= i__1; ++i__) {
#line 328 "slar1v.f"
	    dplus = d__[i__] + s;
#line 329 "slar1v.f"
	    if (abs(dplus) < *pivmin) {
#line 329 "slar1v.f"
		dplus = -(*pivmin);
#line 329 "slar1v.f"
	    }
#line 330 "slar1v.f"
	    work[indlpl + i__] = ld[i__] / dplus;
#line 331 "slar1v.f"
	    if (dplus < 0.) {
#line 331 "slar1v.f"
		++neg1;
#line 331 "slar1v.f"
	    }
#line 332 "slar1v.f"
	    work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 333 "slar1v.f"
	    if (work[indlpl + i__] == 0.) {
#line 333 "slar1v.f"
		work[inds + i__] = lld[i__];
#line 333 "slar1v.f"
	    }
#line 335 "slar1v.f"
	    s = work[inds + i__] - *lambda;
#line 336 "slar1v.f"
/* L70: */
#line 336 "slar1v.f"
	}
#line 337 "slar1v.f"
	i__1 = r2 - 1;
#line 337 "slar1v.f"
	for (i__ = r1; i__ <= i__1; ++i__) {
#line 338 "slar1v.f"
	    dplus = d__[i__] + s;
#line 339 "slar1v.f"
	    if (abs(dplus) < *pivmin) {
#line 339 "slar1v.f"
		dplus = -(*pivmin);
#line 339 "slar1v.f"
	    }
#line 340 "slar1v.f"
	    work[indlpl + i__] = ld[i__] / dplus;
#line 341 "slar1v.f"
	    work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 342 "slar1v.f"
	    if (work[indlpl + i__] == 0.) {
#line 342 "slar1v.f"
		work[inds + i__] = lld[i__];
#line 342 "slar1v.f"
	    }
#line 344 "slar1v.f"
	    s = work[inds + i__] - *lambda;
#line 345 "slar1v.f"
/* L71: */
#line 345 "slar1v.f"
	}
#line 346 "slar1v.f"
    }

/*     Compute the progressive transform (using the differential form) */
/*     until the index R1 */

#line 351 "slar1v.f"
    sawnan2 = FALSE_;
#line 352 "slar1v.f"
    neg2 = 0;
#line 353 "slar1v.f"
    work[indp + *bn - 1] = d__[*bn] - *lambda;
#line 354 "slar1v.f"
    i__1 = r1;
#line 354 "slar1v.f"
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
#line 355 "slar1v.f"
	dminus = lld[i__] + work[indp + i__];
#line 356 "slar1v.f"
	tmp = d__[i__] / dminus;
#line 357 "slar1v.f"
	if (dminus < 0.) {
#line 357 "slar1v.f"
	    ++neg2;
#line 357 "slar1v.f"
	}
#line 358 "slar1v.f"
	work[indumn + i__] = l[i__] * tmp;
#line 359 "slar1v.f"
	work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
#line 360 "slar1v.f"
/* L80: */
#line 360 "slar1v.f"
    }
#line 361 "slar1v.f"
    tmp = work[indp + r1 - 1];
#line 362 "slar1v.f"
    sawnan2 = sisnan_(&tmp);
#line 364 "slar1v.f"
    if (sawnan2) {
/*        Runs a slower version of the above loop if a NaN is detected */
#line 366 "slar1v.f"
	neg2 = 0;
#line 367 "slar1v.f"
	i__1 = r1;
#line 367 "slar1v.f"
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
#line 368 "slar1v.f"
	    dminus = lld[i__] + work[indp + i__];
#line 369 "slar1v.f"
	    if (abs(dminus) < *pivmin) {
#line 369 "slar1v.f"
		dminus = -(*pivmin);
#line 369 "slar1v.f"
	    }
#line 370 "slar1v.f"
	    tmp = d__[i__] / dminus;
#line 371 "slar1v.f"
	    if (dminus < 0.) {
#line 371 "slar1v.f"
		++neg2;
#line 371 "slar1v.f"
	    }
#line 372 "slar1v.f"
	    work[indumn + i__] = l[i__] * tmp;
#line 373 "slar1v.f"
	    work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
#line 374 "slar1v.f"
	    if (tmp == 0.) {
#line 374 "slar1v.f"
		work[indp + i__ - 1] = d__[i__] - *lambda;
#line 374 "slar1v.f"
	    }
#line 376 "slar1v.f"
/* L100: */
#line 376 "slar1v.f"
	}
#line 377 "slar1v.f"
    }

/*     Find the index (from R1 to R2) of the largest (in magnitude) */
/*     diagonal element of the inverse */

#line 382 "slar1v.f"
    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
#line 383 "slar1v.f"
    if (*mingma < 0.) {
#line 383 "slar1v.f"
	++neg1;
#line 383 "slar1v.f"
    }
#line 384 "slar1v.f"
    if (*wantnc) {
#line 385 "slar1v.f"
	*negcnt = neg1 + neg2;
#line 386 "slar1v.f"
    } else {
#line 387 "slar1v.f"
	*negcnt = -1;
#line 388 "slar1v.f"
    }
#line 389 "slar1v.f"
    if (abs(*mingma) == 0.) {
#line 389 "slar1v.f"
	*mingma = eps * work[inds + r1 - 1];
#line 389 "slar1v.f"
    }
#line 391 "slar1v.f"
    *r__ = r1;
#line 392 "slar1v.f"
    i__1 = r2 - 1;
#line 392 "slar1v.f"
    for (i__ = r1; i__ <= i__1; ++i__) {
#line 393 "slar1v.f"
	tmp = work[inds + i__] + work[indp + i__];
#line 394 "slar1v.f"
	if (tmp == 0.) {
#line 394 "slar1v.f"
	    tmp = eps * work[inds + i__];
#line 394 "slar1v.f"
	}
#line 396 "slar1v.f"
	if (abs(tmp) <= abs(*mingma)) {
#line 397 "slar1v.f"
	    *mingma = tmp;
#line 398 "slar1v.f"
	    *r__ = i__ + 1;
#line 399 "slar1v.f"
	}
#line 400 "slar1v.f"
/* L110: */
#line 400 "slar1v.f"
    }

/*     Compute the FP vector: solve N^T v = e_r */

#line 404 "slar1v.f"
    isuppz[1] = *b1;
#line 405 "slar1v.f"
    isuppz[2] = *bn;
#line 406 "slar1v.f"
    z__[*r__] = 1.;
#line 407 "slar1v.f"
    *ztz = 1.;

/*     Compute the FP vector upwards from R */

#line 411 "slar1v.f"
    if (! sawnan1 && ! sawnan2) {
#line 412 "slar1v.f"
	i__1 = *b1;
#line 412 "slar1v.f"
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
#line 413 "slar1v.f"
	    z__[i__] = -(work[indlpl + i__] * z__[i__ + 1]);
#line 414 "slar1v.f"
	    if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(
		    d__2))) * (d__3 = ld[i__], abs(d__3)) < *gaptol) {
#line 416 "slar1v.f"
		z__[i__] = 0.;
#line 417 "slar1v.f"
		isuppz[1] = i__ + 1;
#line 418 "slar1v.f"
		goto L220;
#line 419 "slar1v.f"
	    }
#line 420 "slar1v.f"
	    *ztz += z__[i__] * z__[i__];
#line 421 "slar1v.f"
/* L210: */
#line 421 "slar1v.f"
	}
#line 422 "slar1v.f"
L220:
#line 423 "slar1v.f"
	;
#line 423 "slar1v.f"
    } else {
/*        Run slower loop if NaN occurred. */
#line 425 "slar1v.f"
	i__1 = *b1;
#line 425 "slar1v.f"
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
#line 426 "slar1v.f"
	    if (z__[i__ + 1] == 0.) {
#line 427 "slar1v.f"
		z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
#line 428 "slar1v.f"
	    } else {
#line 429 "slar1v.f"
		z__[i__] = -(work[indlpl + i__] * z__[i__ + 1]);
#line 430 "slar1v.f"
	    }
#line 431 "slar1v.f"
	    if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(
		    d__2))) * (d__3 = ld[i__], abs(d__3)) < *gaptol) {
#line 433 "slar1v.f"
		z__[i__] = 0.;
#line 434 "slar1v.f"
		isuppz[1] = i__ + 1;
#line 435 "slar1v.f"
		goto L240;
#line 436 "slar1v.f"
	    }
#line 437 "slar1v.f"
	    *ztz += z__[i__] * z__[i__];
#line 438 "slar1v.f"
/* L230: */
#line 438 "slar1v.f"
	}
#line 439 "slar1v.f"
L240:
#line 440 "slar1v.f"
	;
#line 440 "slar1v.f"
    }
/*     Compute the FP vector downwards from R in blocks of size BLKSIZ */
#line 443 "slar1v.f"
    if (! sawnan1 && ! sawnan2) {
#line 444 "slar1v.f"
	i__1 = *bn - 1;
#line 444 "slar1v.f"
	for (i__ = *r__; i__ <= i__1; ++i__) {
#line 445 "slar1v.f"
	    z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
#line 446 "slar1v.f"
	    if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(
		    d__2))) * (d__3 = ld[i__], abs(d__3)) < *gaptol) {
#line 448 "slar1v.f"
		z__[i__ + 1] = 0.;
#line 449 "slar1v.f"
		isuppz[2] = i__;
#line 450 "slar1v.f"
		goto L260;
#line 451 "slar1v.f"
	    }
#line 452 "slar1v.f"
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
#line 453 "slar1v.f"
/* L250: */
#line 453 "slar1v.f"
	}
#line 454 "slar1v.f"
L260:
#line 455 "slar1v.f"
	;
#line 455 "slar1v.f"
    } else {
/*        Run slower loop if NaN occurred. */
#line 457 "slar1v.f"
	i__1 = *bn - 1;
#line 457 "slar1v.f"
	for (i__ = *r__; i__ <= i__1; ++i__) {
#line 458 "slar1v.f"
	    if (z__[i__] == 0.) {
#line 459 "slar1v.f"
		z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
#line 460 "slar1v.f"
	    } else {
#line 461 "slar1v.f"
		z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
#line 462 "slar1v.f"
	    }
#line 463 "slar1v.f"
	    if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(
		    d__2))) * (d__3 = ld[i__], abs(d__3)) < *gaptol) {
#line 465 "slar1v.f"
		z__[i__ + 1] = 0.;
#line 466 "slar1v.f"
		isuppz[2] = i__;
#line 467 "slar1v.f"
		goto L280;
#line 468 "slar1v.f"
	    }
#line 469 "slar1v.f"
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
#line 470 "slar1v.f"
/* L270: */
#line 470 "slar1v.f"
	}
#line 471 "slar1v.f"
L280:
#line 472 "slar1v.f"
	;
#line 472 "slar1v.f"
    }

/*     Compute quantities for convergence test */

#line 476 "slar1v.f"
    tmp = 1. / *ztz;
#line 477 "slar1v.f"
    *nrminv = sqrt(tmp);
#line 478 "slar1v.f"
    *resid = abs(*mingma) * *nrminv;
#line 479 "slar1v.f"
    *rqcorr = *mingma * tmp;


#line 482 "slar1v.f"
    return 0;

/*     End of SLAR1V */

} /* slar1v_ */

