#line 1 "zlar1v.f"
/* zlar1v.f -- translated by f2c (version 20100827).
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

#line 1 "zlar1v.f"
/* > \brief \b ZLAR1V computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn 
of the tridiagonal matrix LDLT - λI. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAR1V + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlar1v.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlar1v.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlar1v.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD, */
/*                  PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, */
/*                  R, ISUPPZ, NRMINV, RESID, RQCORR, WORK ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTNC */
/*       INTEGER   B1, BN, N, NEGCNT, R */
/*       DOUBLE PRECISION   GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID, */
/*      $                   RQCORR, ZTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ) */
/*       DOUBLE PRECISION   D( * ), L( * ), LD( * ), LLD( * ), */
/*      $                  WORK( * ) */
/*       COMPLEX*16       Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAR1V computes the (scaled) r-th column of the inverse of */
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
/* >          LAMBDA is DOUBLE PRECISION */
/* >           The shift. In order to compute an accurate eigenvector, */
/* >           LAMBDA should be a good approximation to an eigenvalue */
/* >           of L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is DOUBLE PRECISION array, dimension (N-1) */
/* >           The (n-1) subdiagonal elements of the unit bidiagonal matrix */
/* >           L, in elements 1 to N-1. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >           The n diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LD */
/* > \verbatim */
/* >          LD is DOUBLE PRECISION array, dimension (N-1) */
/* >           The n-1 elements L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] LLD */
/* > \verbatim */
/* >          LLD is DOUBLE PRECISION array, dimension (N-1) */
/* >           The n-1 elements L(i)*L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
/* >           The minimum pivot in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[in] GAPTOL */
/* > \verbatim */
/* >          GAPTOL is DOUBLE PRECISION */
/* >           Tolerance that indicates when eigenvector entries are negligible */
/* >           w.r.t. their contribution to the residual. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (N) */
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
/* >          ZTZ is DOUBLE PRECISION */
/* >           The square of the 2-norm of Z. */
/* > \endverbatim */
/* > */
/* > \param[out] MINGMA */
/* > \verbatim */
/* >          MINGMA is DOUBLE PRECISION */
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
/* >          NRMINV is DOUBLE PRECISION */
/* >           NRMINV = 1/SQRT( ZTZ ) */
/* > \endverbatim */
/* > */
/* > \param[out] RESID */
/* > \verbatim */
/* >          RESID is DOUBLE PRECISION */
/* >           The residual of the FP vector. */
/* >           RESID = ABS( MINGMA )/SQRT( ZTZ ) */
/* > \endverbatim */
/* > */
/* > \param[out] RQCORR */
/* > \verbatim */
/* >          RQCORR is DOUBLE PRECISION */
/* >           The Rayleigh Quotient correction to LAMBDA. */
/* >           RQCORR = MINGMA*TMP */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int zlar1v_(integer *n, integer *b1, integer *bn, doublereal 
	*lambda, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *pivmin, doublereal *gaptol, doublecomplex *z__, 
	logical *wantnc, integer *negcnt, doublereal *ztz, doublereal *mingma,
	 integer *r__, integer *isuppz, doublereal *nrminv, doublereal *resid,
	 doublereal *rqcorr, doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s;
    static integer r1, r2;
    static doublereal eps, tmp;
    static integer neg1, neg2, indp, inds;
    static doublereal dplus;
    extern doublereal dlamch_(char *, ftnlen);
    extern logical disnan_(doublereal *);
    static integer indlpl, indumn;
    static doublereal dminus;
    static logical sawnan1, sawnan2;


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 276 "zlar1v.f"
    /* Parameter adjustments */
#line 276 "zlar1v.f"
    --work;
#line 276 "zlar1v.f"
    --isuppz;
#line 276 "zlar1v.f"
    --z__;
#line 276 "zlar1v.f"
    --lld;
#line 276 "zlar1v.f"
    --ld;
#line 276 "zlar1v.f"
    --l;
#line 276 "zlar1v.f"
    --d__;
#line 276 "zlar1v.f"

#line 276 "zlar1v.f"
    /* Function Body */
#line 276 "zlar1v.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 279 "zlar1v.f"
    if (*r__ == 0) {
#line 280 "zlar1v.f"
	r1 = *b1;
#line 281 "zlar1v.f"
	r2 = *bn;
#line 282 "zlar1v.f"
    } else {
#line 283 "zlar1v.f"
	r1 = *r__;
#line 284 "zlar1v.f"
	r2 = *r__;
#line 285 "zlar1v.f"
    }
/*     Storage for LPLUS */
#line 288 "zlar1v.f"
    indlpl = 0;
/*     Storage for UMINUS */
#line 290 "zlar1v.f"
    indumn = *n;
#line 291 "zlar1v.f"
    inds = (*n << 1) + 1;
#line 292 "zlar1v.f"
    indp = *n * 3 + 1;
#line 294 "zlar1v.f"
    if (*b1 == 1) {
#line 295 "zlar1v.f"
	work[inds] = 0.;
#line 296 "zlar1v.f"
    } else {
#line 297 "zlar1v.f"
	work[inds + *b1 - 1] = lld[*b1 - 1];
#line 298 "zlar1v.f"
    }

/*     Compute the stationary transform (using the differential form) */
/*     until the index R2. */

#line 304 "zlar1v.f"
    sawnan1 = FALSE_;
#line 305 "zlar1v.f"
    neg1 = 0;
#line 306 "zlar1v.f"
    s = work[inds + *b1 - 1] - *lambda;
#line 307 "zlar1v.f"
    i__1 = r1 - 1;
#line 307 "zlar1v.f"
    for (i__ = *b1; i__ <= i__1; ++i__) {
#line 308 "zlar1v.f"
	dplus = d__[i__] + s;
#line 309 "zlar1v.f"
	work[indlpl + i__] = ld[i__] / dplus;
#line 310 "zlar1v.f"
	if (dplus < 0.) {
#line 310 "zlar1v.f"
	    ++neg1;
#line 310 "zlar1v.f"
	}
#line 311 "zlar1v.f"
	work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 312 "zlar1v.f"
	s = work[inds + i__] - *lambda;
#line 313 "zlar1v.f"
/* L50: */
#line 313 "zlar1v.f"
    }
#line 314 "zlar1v.f"
    sawnan1 = disnan_(&s);
#line 315 "zlar1v.f"
    if (sawnan1) {
#line 315 "zlar1v.f"
	goto L60;
#line 315 "zlar1v.f"
    }
#line 316 "zlar1v.f"
    i__1 = r2 - 1;
#line 316 "zlar1v.f"
    for (i__ = r1; i__ <= i__1; ++i__) {
#line 317 "zlar1v.f"
	dplus = d__[i__] + s;
#line 318 "zlar1v.f"
	work[indlpl + i__] = ld[i__] / dplus;
#line 319 "zlar1v.f"
	work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 320 "zlar1v.f"
	s = work[inds + i__] - *lambda;
#line 321 "zlar1v.f"
/* L51: */
#line 321 "zlar1v.f"
    }
#line 322 "zlar1v.f"
    sawnan1 = disnan_(&s);

#line 324 "zlar1v.f"
L60:
#line 325 "zlar1v.f"
    if (sawnan1) {
/*        Runs a slower version of the above loop if a NaN is detected */
#line 327 "zlar1v.f"
	neg1 = 0;
#line 328 "zlar1v.f"
	s = work[inds + *b1 - 1] - *lambda;
#line 329 "zlar1v.f"
	i__1 = r1 - 1;
#line 329 "zlar1v.f"
	for (i__ = *b1; i__ <= i__1; ++i__) {
#line 330 "zlar1v.f"
	    dplus = d__[i__] + s;
#line 331 "zlar1v.f"
	    if (abs(dplus) < *pivmin) {
#line 331 "zlar1v.f"
		dplus = -(*pivmin);
#line 331 "zlar1v.f"
	    }
#line 332 "zlar1v.f"
	    work[indlpl + i__] = ld[i__] / dplus;
#line 333 "zlar1v.f"
	    if (dplus < 0.) {
#line 333 "zlar1v.f"
		++neg1;
#line 333 "zlar1v.f"
	    }
#line 334 "zlar1v.f"
	    work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 335 "zlar1v.f"
	    if (work[indlpl + i__] == 0.) {
#line 335 "zlar1v.f"
		work[inds + i__] = lld[i__];
#line 335 "zlar1v.f"
	    }
#line 337 "zlar1v.f"
	    s = work[inds + i__] - *lambda;
#line 338 "zlar1v.f"
/* L70: */
#line 338 "zlar1v.f"
	}
#line 339 "zlar1v.f"
	i__1 = r2 - 1;
#line 339 "zlar1v.f"
	for (i__ = r1; i__ <= i__1; ++i__) {
#line 340 "zlar1v.f"
	    dplus = d__[i__] + s;
#line 341 "zlar1v.f"
	    if (abs(dplus) < *pivmin) {
#line 341 "zlar1v.f"
		dplus = -(*pivmin);
#line 341 "zlar1v.f"
	    }
#line 342 "zlar1v.f"
	    work[indlpl + i__] = ld[i__] / dplus;
#line 343 "zlar1v.f"
	    work[inds + i__] = s * work[indlpl + i__] * l[i__];
#line 344 "zlar1v.f"
	    if (work[indlpl + i__] == 0.) {
#line 344 "zlar1v.f"
		work[inds + i__] = lld[i__];
#line 344 "zlar1v.f"
	    }
#line 346 "zlar1v.f"
	    s = work[inds + i__] - *lambda;
#line 347 "zlar1v.f"
/* L71: */
#line 347 "zlar1v.f"
	}
#line 348 "zlar1v.f"
    }

/*     Compute the progressive transform (using the differential form) */
/*     until the index R1 */

#line 353 "zlar1v.f"
    sawnan2 = FALSE_;
#line 354 "zlar1v.f"
    neg2 = 0;
#line 355 "zlar1v.f"
    work[indp + *bn - 1] = d__[*bn] - *lambda;
#line 356 "zlar1v.f"
    i__1 = r1;
#line 356 "zlar1v.f"
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
#line 357 "zlar1v.f"
	dminus = lld[i__] + work[indp + i__];
#line 358 "zlar1v.f"
	tmp = d__[i__] / dminus;
#line 359 "zlar1v.f"
	if (dminus < 0.) {
#line 359 "zlar1v.f"
	    ++neg2;
#line 359 "zlar1v.f"
	}
#line 360 "zlar1v.f"
	work[indumn + i__] = l[i__] * tmp;
#line 361 "zlar1v.f"
	work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
#line 362 "zlar1v.f"
/* L80: */
#line 362 "zlar1v.f"
    }
#line 363 "zlar1v.f"
    tmp = work[indp + r1 - 1];
#line 364 "zlar1v.f"
    sawnan2 = disnan_(&tmp);
#line 366 "zlar1v.f"
    if (sawnan2) {
/*        Runs a slower version of the above loop if a NaN is detected */
#line 368 "zlar1v.f"
	neg2 = 0;
#line 369 "zlar1v.f"
	i__1 = r1;
#line 369 "zlar1v.f"
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
#line 370 "zlar1v.f"
	    dminus = lld[i__] + work[indp + i__];
#line 371 "zlar1v.f"
	    if (abs(dminus) < *pivmin) {
#line 371 "zlar1v.f"
		dminus = -(*pivmin);
#line 371 "zlar1v.f"
	    }
#line 372 "zlar1v.f"
	    tmp = d__[i__] / dminus;
#line 373 "zlar1v.f"
	    if (dminus < 0.) {
#line 373 "zlar1v.f"
		++neg2;
#line 373 "zlar1v.f"
	    }
#line 374 "zlar1v.f"
	    work[indumn + i__] = l[i__] * tmp;
#line 375 "zlar1v.f"
	    work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
#line 376 "zlar1v.f"
	    if (tmp == 0.) {
#line 376 "zlar1v.f"
		work[indp + i__ - 1] = d__[i__] - *lambda;
#line 376 "zlar1v.f"
	    }
#line 378 "zlar1v.f"
/* L100: */
#line 378 "zlar1v.f"
	}
#line 379 "zlar1v.f"
    }

/*     Find the index (from R1 to R2) of the largest (in magnitude) */
/*     diagonal element of the inverse */

#line 384 "zlar1v.f"
    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
#line 385 "zlar1v.f"
    if (*mingma < 0.) {
#line 385 "zlar1v.f"
	++neg1;
#line 385 "zlar1v.f"
    }
#line 386 "zlar1v.f"
    if (*wantnc) {
#line 387 "zlar1v.f"
	*negcnt = neg1 + neg2;
#line 388 "zlar1v.f"
    } else {
#line 389 "zlar1v.f"
	*negcnt = -1;
#line 390 "zlar1v.f"
    }
#line 391 "zlar1v.f"
    if (abs(*mingma) == 0.) {
#line 391 "zlar1v.f"
	*mingma = eps * work[inds + r1 - 1];
#line 391 "zlar1v.f"
    }
#line 393 "zlar1v.f"
    *r__ = r1;
#line 394 "zlar1v.f"
    i__1 = r2 - 1;
#line 394 "zlar1v.f"
    for (i__ = r1; i__ <= i__1; ++i__) {
#line 395 "zlar1v.f"
	tmp = work[inds + i__] + work[indp + i__];
#line 396 "zlar1v.f"
	if (tmp == 0.) {
#line 396 "zlar1v.f"
	    tmp = eps * work[inds + i__];
#line 396 "zlar1v.f"
	}
#line 398 "zlar1v.f"
	if (abs(tmp) <= abs(*mingma)) {
#line 399 "zlar1v.f"
	    *mingma = tmp;
#line 400 "zlar1v.f"
	    *r__ = i__ + 1;
#line 401 "zlar1v.f"
	}
#line 402 "zlar1v.f"
/* L110: */
#line 402 "zlar1v.f"
    }

/*     Compute the FP vector: solve N^T v = e_r */

#line 406 "zlar1v.f"
    isuppz[1] = *b1;
#line 407 "zlar1v.f"
    isuppz[2] = *bn;
#line 408 "zlar1v.f"
    i__1 = *r__;
#line 408 "zlar1v.f"
    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 409 "zlar1v.f"
    *ztz = 1.;

/*     Compute the FP vector upwards from R */

#line 413 "zlar1v.f"
    if (! sawnan1 && ! sawnan2) {
#line 414 "zlar1v.f"
	i__1 = *b1;
#line 414 "zlar1v.f"
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
#line 415 "zlar1v.f"
	    i__2 = i__;
#line 415 "zlar1v.f"
	    i__3 = indlpl + i__;
#line 415 "zlar1v.f"
	    i__4 = i__ + 1;
#line 415 "zlar1v.f"
	    z__2.r = work[i__3] * z__[i__4].r, z__2.i = work[i__3] * z__[i__4]
		    .i;
#line 415 "zlar1v.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 415 "zlar1v.f"
	    z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
#line 416 "zlar1v.f"
	    if ((z_abs(&z__[i__]) + z_abs(&z__[i__ + 1])) * (d__1 = ld[i__], 
		    abs(d__1)) < *gaptol) {
#line 418 "zlar1v.f"
		i__2 = i__;
#line 418 "zlar1v.f"
		z__[i__2].r = 0., z__[i__2].i = 0.;
#line 419 "zlar1v.f"
		isuppz[1] = i__ + 1;
#line 420 "zlar1v.f"
		goto L220;
#line 421 "zlar1v.f"
	    }
#line 422 "zlar1v.f"
	    i__2 = i__;
#line 422 "zlar1v.f"
	    i__3 = i__;
#line 422 "zlar1v.f"
	    z__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i, 
		    z__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[
		    i__3].r;
#line 422 "zlar1v.f"
	    *ztz += z__1.r;
#line 423 "zlar1v.f"
/* L210: */
#line 423 "zlar1v.f"
	}
#line 424 "zlar1v.f"
L220:
#line 425 "zlar1v.f"
	;
#line 425 "zlar1v.f"
    } else {
/*        Run slower loop if NaN occurred. */
#line 427 "zlar1v.f"
	i__1 = *b1;
#line 427 "zlar1v.f"
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
#line 428 "zlar1v.f"
	    i__2 = i__ + 1;
#line 428 "zlar1v.f"
	    if (z__[i__2].r == 0. && z__[i__2].i == 0.) {
#line 429 "zlar1v.f"
		i__2 = i__;
#line 429 "zlar1v.f"
		d__1 = -(ld[i__ + 1] / ld[i__]);
#line 429 "zlar1v.f"
		i__3 = i__ + 2;
#line 429 "zlar1v.f"
		z__1.r = d__1 * z__[i__3].r, z__1.i = d__1 * z__[i__3].i;
#line 429 "zlar1v.f"
		z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
#line 430 "zlar1v.f"
	    } else {
#line 431 "zlar1v.f"
		i__2 = i__;
#line 431 "zlar1v.f"
		i__3 = indlpl + i__;
#line 431 "zlar1v.f"
		i__4 = i__ + 1;
#line 431 "zlar1v.f"
		z__2.r = work[i__3] * z__[i__4].r, z__2.i = work[i__3] * z__[
			i__4].i;
#line 431 "zlar1v.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 431 "zlar1v.f"
		z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
#line 432 "zlar1v.f"
	    }
#line 433 "zlar1v.f"
	    if ((z_abs(&z__[i__]) + z_abs(&z__[i__ + 1])) * (d__1 = ld[i__], 
		    abs(d__1)) < *gaptol) {
#line 435 "zlar1v.f"
		i__2 = i__;
#line 435 "zlar1v.f"
		z__[i__2].r = 0., z__[i__2].i = 0.;
#line 436 "zlar1v.f"
		isuppz[1] = i__ + 1;
#line 437 "zlar1v.f"
		goto L240;
#line 438 "zlar1v.f"
	    }
#line 439 "zlar1v.f"
	    i__2 = i__;
#line 439 "zlar1v.f"
	    i__3 = i__;
#line 439 "zlar1v.f"
	    z__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i, 
		    z__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[
		    i__3].r;
#line 439 "zlar1v.f"
	    *ztz += z__1.r;
#line 440 "zlar1v.f"
/* L230: */
#line 440 "zlar1v.f"
	}
#line 441 "zlar1v.f"
L240:
#line 442 "zlar1v.f"
	;
#line 442 "zlar1v.f"
    }
/*     Compute the FP vector downwards from R in blocks of size BLKSIZ */
#line 445 "zlar1v.f"
    if (! sawnan1 && ! sawnan2) {
#line 446 "zlar1v.f"
	i__1 = *bn - 1;
#line 446 "zlar1v.f"
	for (i__ = *r__; i__ <= i__1; ++i__) {
#line 447 "zlar1v.f"
	    i__2 = i__ + 1;
#line 447 "zlar1v.f"
	    i__3 = indumn + i__;
#line 447 "zlar1v.f"
	    i__4 = i__;
#line 447 "zlar1v.f"
	    z__2.r = work[i__3] * z__[i__4].r, z__2.i = work[i__3] * z__[i__4]
		    .i;
#line 447 "zlar1v.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 447 "zlar1v.f"
	    z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
#line 448 "zlar1v.f"
	    if ((z_abs(&z__[i__]) + z_abs(&z__[i__ + 1])) * (d__1 = ld[i__], 
		    abs(d__1)) < *gaptol) {
#line 450 "zlar1v.f"
		i__2 = i__ + 1;
#line 450 "zlar1v.f"
		z__[i__2].r = 0., z__[i__2].i = 0.;
#line 451 "zlar1v.f"
		isuppz[2] = i__;
#line 452 "zlar1v.f"
		goto L260;
#line 453 "zlar1v.f"
	    }
#line 454 "zlar1v.f"
	    i__2 = i__ + 1;
#line 454 "zlar1v.f"
	    i__3 = i__ + 1;
#line 454 "zlar1v.f"
	    z__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i, 
		    z__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[
		    i__3].r;
#line 454 "zlar1v.f"
	    *ztz += z__1.r;
#line 455 "zlar1v.f"
/* L250: */
#line 455 "zlar1v.f"
	}
#line 456 "zlar1v.f"
L260:
#line 457 "zlar1v.f"
	;
#line 457 "zlar1v.f"
    } else {
/*        Run slower loop if NaN occurred. */
#line 459 "zlar1v.f"
	i__1 = *bn - 1;
#line 459 "zlar1v.f"
	for (i__ = *r__; i__ <= i__1; ++i__) {
#line 460 "zlar1v.f"
	    i__2 = i__;
#line 460 "zlar1v.f"
	    if (z__[i__2].r == 0. && z__[i__2].i == 0.) {
#line 461 "zlar1v.f"
		i__2 = i__ + 1;
#line 461 "zlar1v.f"
		d__1 = -(ld[i__ - 1] / ld[i__]);
#line 461 "zlar1v.f"
		i__3 = i__ - 1;
#line 461 "zlar1v.f"
		z__1.r = d__1 * z__[i__3].r, z__1.i = d__1 * z__[i__3].i;
#line 461 "zlar1v.f"
		z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
#line 462 "zlar1v.f"
	    } else {
#line 463 "zlar1v.f"
		i__2 = i__ + 1;
#line 463 "zlar1v.f"
		i__3 = indumn + i__;
#line 463 "zlar1v.f"
		i__4 = i__;
#line 463 "zlar1v.f"
		z__2.r = work[i__3] * z__[i__4].r, z__2.i = work[i__3] * z__[
			i__4].i;
#line 463 "zlar1v.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 463 "zlar1v.f"
		z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
#line 464 "zlar1v.f"
	    }
#line 465 "zlar1v.f"
	    if ((z_abs(&z__[i__]) + z_abs(&z__[i__ + 1])) * (d__1 = ld[i__], 
		    abs(d__1)) < *gaptol) {
#line 467 "zlar1v.f"
		i__2 = i__ + 1;
#line 467 "zlar1v.f"
		z__[i__2].r = 0., z__[i__2].i = 0.;
#line 468 "zlar1v.f"
		isuppz[2] = i__;
#line 469 "zlar1v.f"
		goto L280;
#line 470 "zlar1v.f"
	    }
#line 471 "zlar1v.f"
	    i__2 = i__ + 1;
#line 471 "zlar1v.f"
	    i__3 = i__ + 1;
#line 471 "zlar1v.f"
	    z__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i, 
		    z__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[
		    i__3].r;
#line 471 "zlar1v.f"
	    *ztz += z__1.r;
#line 472 "zlar1v.f"
/* L270: */
#line 472 "zlar1v.f"
	}
#line 473 "zlar1v.f"
L280:
#line 474 "zlar1v.f"
	;
#line 474 "zlar1v.f"
    }

/*     Compute quantities for convergence test */

#line 478 "zlar1v.f"
    tmp = 1. / *ztz;
#line 479 "zlar1v.f"
    *nrminv = sqrt(tmp);
#line 480 "zlar1v.f"
    *resid = abs(*mingma) * *nrminv;
#line 481 "zlar1v.f"
    *rqcorr = *mingma * tmp;


#line 484 "zlar1v.f"
    return 0;

/*     End of ZLAR1V */

} /* zlar1v_ */

