#line 1 "slarrf.f"
/* slarrf.f -- translated by f2c (version 20100827).
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

#line 1 "slarrf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLARRF finds a new relatively robust representation such that at least one of the eigenvalues i
s relatively isolated. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARRF( N, D, L, LD, CLSTRT, CLEND, */
/*                          W, WGAP, WERR, */
/*                          SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA, */
/*                          DPLUS, LPLUS, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            CLSTRT, CLEND, INFO, N */
/*       REAL               CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), DPLUS( * ), L( * ), LD( * ), */
/*      $          LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given the initial representation L D L^T and its cluster of close */
/* > eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ... */
/* > W( CLEND ), SLARRF finds a new relatively robust representation */
/* > L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the */
/* > eigenvalues of L(+) D(+) L(+)^T is relatively isolated. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix (subblock, if the matrix splitted). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is REAL array, dimension (N-1) */
/* >          The (N-1) subdiagonal elements of the unit bidiagonal */
/* >          matrix L. */
/* > \endverbatim */
/* > */
/* > \param[in] LD */
/* > \verbatim */
/* >          LD is REAL array, dimension (N-1) */
/* >          The (N-1) elements L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] CLSTRT */
/* > \verbatim */
/* >          CLSTRT is INTEGER */
/* >          The index of the first eigenvalue in the cluster. */
/* > \endverbatim */
/* > */
/* > \param[in] CLEND */
/* > \verbatim */
/* >          CLEND is INTEGER */
/* >          The index of the last eigenvalue in the cluster. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is REAL array, dimension */
/* >          dimension is >=  (CLEND-CLSTRT+1) */
/* >          The eigenvalue APPROXIMATIONS of L D L^T in ascending order. */
/* >          W( CLSTRT ) through W( CLEND ) form the cluster of relatively */
/* >          close eigenalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WGAP */
/* > \verbatim */
/* >          WGAP is REAL array, dimension */
/* >          dimension is >=  (CLEND-CLSTRT+1) */
/* >          The separation from the right neighbor eigenvalue in W. */
/* > \endverbatim */
/* > */
/* > \param[in] WERR */
/* > \verbatim */
/* >          WERR is REAL array, dimension */
/* >          dimension is >=  (CLEND-CLSTRT+1) */
/* >          WERR contain the semiwidth of the uncertainty */
/* >          interval of the corresponding eigenvalue APPROXIMATION in W */
/* > \endverbatim */
/* > */
/* > \param[in] SPDIAM */
/* > \verbatim */
/* >          SPDIAM is REAL */
/* >          estimate of the spectral diameter obtained from the */
/* >          Gerschgorin intervals */
/* > \endverbatim */
/* > */
/* > \param[in] CLGAPL */
/* > \verbatim */
/* >          CLGAPL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] CLGAPR */
/* > \verbatim */
/* >          CLGAPR is REAL */
/* >          absolute gap on each end of the cluster. */
/* >          Set by the calling routine to protect against shifts too close */
/* >          to eigenvalues outside the cluster. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
/* >          The minimum pivot allowed in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* >          SIGMA is REAL */
/* >          The shift used to form L(+) D(+) L(+)^T. */
/* > \endverbatim */
/* > */
/* > \param[out] DPLUS */
/* > \verbatim */
/* >          DPLUS is REAL array, dimension (N) */
/* >          The N diagonal elements of the diagonal matrix D(+). */
/* > \endverbatim */
/* > */
/* > \param[out] LPLUS */
/* > \verbatim */
/* >          LPLUS is REAL array, dimension (N-1) */
/* >          The first (N-1) elements of LPLUS contain the subdiagonal */
/* >          elements of the unit bidiagonal matrix L(+). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (2*N) */
/* >          Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          Signals processing OK (=0) or failure (=1) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int slarrf_(integer *n, doublereal *d__, doublereal *l, 
	doublereal *ld, integer *clstrt, integer *clend, doublereal *w, 
	doublereal *wgap, doublereal *werr, doublereal *spdiam, doublereal *
	clgapl, doublereal *clgapr, doublereal *pivmin, doublereal *sigma, 
	doublereal *dplus, doublereal *lplus, doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s, bestshift, smlgrowth, eps, tmp, max1, max2, rrr1, 
	    rrr2, znm2, growthbound, fail, fact, oldp;
    static integer indx;
    static doublereal prod;
    static integer ktry;
    static doublereal fail2, avgap, ldmax, rdmax;
    static integer shift;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical dorrr1;
    static doublereal ldelta;
    extern doublereal slamch_(char *, ftnlen);
    static logical nofail;
    static doublereal mingap, lsigma, rdelta;
    static logical forcer;
    static doublereal rsigma, clwdth;
    extern logical sisnan_(doublereal *);
    static logical sawnan1, sawnan2, tryrrr1;


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 241 "slarrf.f"
    /* Parameter adjustments */
#line 241 "slarrf.f"
    --work;
#line 241 "slarrf.f"
    --lplus;
#line 241 "slarrf.f"
    --dplus;
#line 241 "slarrf.f"
    --werr;
#line 241 "slarrf.f"
    --wgap;
#line 241 "slarrf.f"
    --w;
#line 241 "slarrf.f"
    --ld;
#line 241 "slarrf.f"
    --l;
#line 241 "slarrf.f"
    --d__;
#line 241 "slarrf.f"

#line 241 "slarrf.f"
    /* Function Body */
#line 241 "slarrf.f"
    *info = 0;
#line 242 "slarrf.f"
    fact = 2.;
#line 243 "slarrf.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 244 "slarrf.f"
    shift = 0;
#line 245 "slarrf.f"
    forcer = FALSE_;
/*     Note that we cannot guarantee that for any of the shifts tried, */
/*     the factorization has a small or even moderate element growth. */
/*     There could be Ritz values at both ends of the cluster and despite */
/*     backing off, there are examples where all factorizations tried */
/*     (in IEEE mode, allowing zero pivots & infinities) have INFINITE */
/*     element growth. */
/*     For this reason, we should use PIVMIN in this subroutine so that at */
/*     least the L D L^T factorization exists. It can be checked afterwards */
/*     whether the element growth caused bad residuals/orthogonality. */
/*     Decide whether the code should accept the best among all */
/*     representations despite large element growth or signal INFO=1 */
/*     Setting NOFAIL to .FALSE. for quick fix for bug 113 */
#line 261 "slarrf.f"
    nofail = FALSE_;

/*     Compute the average gap length of the cluster */
#line 265 "slarrf.f"
    clwdth = (d__1 = w[*clend] - w[*clstrt], abs(d__1)) + werr[*clend] + werr[
	    *clstrt];
#line 266 "slarrf.f"
    avgap = clwdth / (doublereal) (*clend - *clstrt);
#line 267 "slarrf.f"
    mingap = min(*clgapl,*clgapr);
/*     Initial values for shifts to both ends of cluster */
/* Computing MIN */
#line 269 "slarrf.f"
    d__1 = w[*clstrt], d__2 = w[*clend];
#line 269 "slarrf.f"
    lsigma = min(d__1,d__2) - werr[*clstrt];
/* Computing MAX */
#line 270 "slarrf.f"
    d__1 = w[*clstrt], d__2 = w[*clend];
#line 270 "slarrf.f"
    rsigma = max(d__1,d__2) + werr[*clend];
/*     Use a small fudge to make sure that we really shift to the outside */
#line 273 "slarrf.f"
    lsigma -= abs(lsigma) * 2. * eps;
#line 274 "slarrf.f"
    rsigma += abs(rsigma) * 2. * eps;
/*     Compute upper bounds for how much to back off the initial shifts */
#line 277 "slarrf.f"
    ldmax = mingap * .25 + *pivmin * 2.;
#line 278 "slarrf.f"
    rdmax = mingap * .25 + *pivmin * 2.;
/* Computing MAX */
#line 280 "slarrf.f"
    d__1 = avgap, d__2 = wgap[*clstrt];
#line 280 "slarrf.f"
    ldelta = max(d__1,d__2) / fact;
/* Computing MAX */
#line 281 "slarrf.f"
    d__1 = avgap, d__2 = wgap[*clend - 1];
#line 281 "slarrf.f"
    rdelta = max(d__1,d__2) / fact;

/*     Initialize the record of the best representation found */

#line 285 "slarrf.f"
    s = slamch_("S", (ftnlen)1);
#line 286 "slarrf.f"
    smlgrowth = 1. / s;
#line 287 "slarrf.f"
    fail = (doublereal) (*n - 1) * mingap / (*spdiam * eps);
#line 288 "slarrf.f"
    fail2 = (doublereal) (*n - 1) * mingap / (*spdiam * sqrt(eps));
#line 289 "slarrf.f"
    bestshift = lsigma;

/*     while (KTRY <= KTRYMAX) */
#line 292 "slarrf.f"
    ktry = 0;
#line 293 "slarrf.f"
    growthbound = *spdiam * 8.;
#line 295 "slarrf.f"
L5:
#line 296 "slarrf.f"
    sawnan1 = FALSE_;
#line 297 "slarrf.f"
    sawnan2 = FALSE_;
/*     Ensure that we do not back off too much of the initial shifts */
#line 299 "slarrf.f"
    ldelta = min(ldmax,ldelta);
#line 300 "slarrf.f"
    rdelta = min(rdmax,rdelta);
/*     Compute the element growth when shifting to both ends of the cluster */
/*     accept the shift if there is no element growth at one of the two ends */
/*     Left end */
#line 306 "slarrf.f"
    s = -lsigma;
#line 307 "slarrf.f"
    dplus[1] = d__[1] + s;
#line 308 "slarrf.f"
    if (abs(dplus[1]) < *pivmin) {
#line 309 "slarrf.f"
	dplus[1] = -(*pivmin);
/*        Need to set SAWNAN1 because refined RRR test should not be used */
/*        in this case */
#line 312 "slarrf.f"
	sawnan1 = TRUE_;
#line 313 "slarrf.f"
    }
#line 314 "slarrf.f"
    max1 = abs(dplus[1]);
#line 315 "slarrf.f"
    i__1 = *n - 1;
#line 315 "slarrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "slarrf.f"
	lplus[i__] = ld[i__] / dplus[i__];
#line 317 "slarrf.f"
	s = s * lplus[i__] * l[i__] - lsigma;
#line 318 "slarrf.f"
	dplus[i__ + 1] = d__[i__ + 1] + s;
#line 319 "slarrf.f"
	if ((d__1 = dplus[i__ + 1], abs(d__1)) < *pivmin) {
#line 320 "slarrf.f"
	    dplus[i__ + 1] = -(*pivmin);
/*           Need to set SAWNAN1 because refined RRR test should not be used */
/*           in this case */
#line 323 "slarrf.f"
	    sawnan1 = TRUE_;
#line 324 "slarrf.f"
	}
/* Computing MAX */
#line 325 "slarrf.f"
	d__2 = max1, d__3 = (d__1 = dplus[i__ + 1], abs(d__1));
#line 325 "slarrf.f"
	max1 = max(d__2,d__3);
#line 326 "slarrf.f"
/* L6: */
#line 326 "slarrf.f"
    }
#line 327 "slarrf.f"
    sawnan1 = sawnan1 || sisnan_(&max1);
#line 329 "slarrf.f"
    if (forcer || max1 <= growthbound && ! sawnan1) {
#line 331 "slarrf.f"
	*sigma = lsigma;
#line 332 "slarrf.f"
	shift = 1;
#line 333 "slarrf.f"
	goto L100;
#line 334 "slarrf.f"
    }
/*     Right end */
#line 337 "slarrf.f"
    s = -rsigma;
#line 338 "slarrf.f"
    work[1] = d__[1] + s;
#line 339 "slarrf.f"
    if (abs(work[1]) < *pivmin) {
#line 340 "slarrf.f"
	work[1] = -(*pivmin);
/*        Need to set SAWNAN2 because refined RRR test should not be used */
/*        in this case */
#line 343 "slarrf.f"
	sawnan2 = TRUE_;
#line 344 "slarrf.f"
    }
#line 345 "slarrf.f"
    max2 = abs(work[1]);
#line 346 "slarrf.f"
    i__1 = *n - 1;
#line 346 "slarrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 347 "slarrf.f"
	work[*n + i__] = ld[i__] / work[i__];
#line 348 "slarrf.f"
	s = s * work[*n + i__] * l[i__] - rsigma;
#line 349 "slarrf.f"
	work[i__ + 1] = d__[i__ + 1] + s;
#line 350 "slarrf.f"
	if ((d__1 = work[i__ + 1], abs(d__1)) < *pivmin) {
#line 351 "slarrf.f"
	    work[i__ + 1] = -(*pivmin);
/*           Need to set SAWNAN2 because refined RRR test should not be used */
/*           in this case */
#line 354 "slarrf.f"
	    sawnan2 = TRUE_;
#line 355 "slarrf.f"
	}
/* Computing MAX */
#line 356 "slarrf.f"
	d__2 = max2, d__3 = (d__1 = work[i__ + 1], abs(d__1));
#line 356 "slarrf.f"
	max2 = max(d__2,d__3);
#line 357 "slarrf.f"
/* L7: */
#line 357 "slarrf.f"
    }
#line 358 "slarrf.f"
    sawnan2 = sawnan2 || sisnan_(&max2);
#line 360 "slarrf.f"
    if (forcer || max2 <= growthbound && ! sawnan2) {
#line 362 "slarrf.f"
	*sigma = rsigma;
#line 363 "slarrf.f"
	shift = 2;
#line 364 "slarrf.f"
	goto L100;
#line 365 "slarrf.f"
    }
/*     If we are at this point, both shifts led to too much element growth */
/*     Record the better of the two shifts (provided it didn't lead to NaN) */
#line 369 "slarrf.f"
    if (sawnan1 && sawnan2) {
/*        both MAX1 and MAX2 are NaN */
#line 371 "slarrf.f"
	goto L50;
#line 372 "slarrf.f"
    } else {
#line 373 "slarrf.f"
	if (! sawnan1) {
#line 374 "slarrf.f"
	    indx = 1;
#line 375 "slarrf.f"
	    if (max1 <= smlgrowth) {
#line 376 "slarrf.f"
		smlgrowth = max1;
#line 377 "slarrf.f"
		bestshift = lsigma;
#line 378 "slarrf.f"
	    }
#line 379 "slarrf.f"
	}
#line 380 "slarrf.f"
	if (! sawnan2) {
#line 381 "slarrf.f"
	    if (sawnan1 || max2 <= max1) {
#line 381 "slarrf.f"
		indx = 2;
#line 381 "slarrf.f"
	    }
#line 382 "slarrf.f"
	    if (max2 <= smlgrowth) {
#line 383 "slarrf.f"
		smlgrowth = max2;
#line 384 "slarrf.f"
		bestshift = rsigma;
#line 385 "slarrf.f"
	    }
#line 386 "slarrf.f"
	}
#line 387 "slarrf.f"
    }
/*     If we are here, both the left and the right shift led to */
/*     element growth. If the element growth is moderate, then */
/*     we may still accept the representation, if it passes a */
/*     refined test for RRR. This test supposes that no NaN occurred. */
/*     Moreover, we use the refined RRR test only for isolated clusters. */
#line 394 "slarrf.f"
    if (clwdth < mingap / 128. && min(max1,max2) < fail2 && ! sawnan1 && ! 
	    sawnan2) {
#line 397 "slarrf.f"
	dorrr1 = TRUE_;
#line 398 "slarrf.f"
    } else {
#line 399 "slarrf.f"
	dorrr1 = FALSE_;
#line 400 "slarrf.f"
    }
#line 401 "slarrf.f"
    tryrrr1 = TRUE_;
#line 402 "slarrf.f"
    if (tryrrr1 && dorrr1) {
#line 403 "slarrf.f"
	if (indx == 1) {
#line 404 "slarrf.f"
	    tmp = (d__1 = dplus[*n], abs(d__1));
#line 405 "slarrf.f"
	    znm2 = 1.;
#line 406 "slarrf.f"
	    prod = 1.;
#line 407 "slarrf.f"
	    oldp = 1.;
#line 408 "slarrf.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 409 "slarrf.f"
		if (prod <= eps) {
#line 410 "slarrf.f"
		    prod = dplus[i__ + 1] * work[*n + i__ + 1] / (dplus[i__] *
			     work[*n + i__]) * oldp;
#line 412 "slarrf.f"
		} else {
#line 413 "slarrf.f"
		    prod *= (d__1 = work[*n + i__], abs(d__1));
#line 414 "slarrf.f"
		}
#line 415 "slarrf.f"
		oldp = prod;
/* Computing 2nd power */
#line 416 "slarrf.f"
		d__1 = prod;
#line 416 "slarrf.f"
		znm2 += d__1 * d__1;
/* Computing MAX */
#line 417 "slarrf.f"
		d__2 = tmp, d__3 = (d__1 = dplus[i__] * prod, abs(d__1));
#line 417 "slarrf.f"
		tmp = max(d__2,d__3);
#line 418 "slarrf.f"
/* L15: */
#line 418 "slarrf.f"
	    }
#line 419 "slarrf.f"
	    rrr1 = tmp / (*spdiam * sqrt(znm2));
#line 420 "slarrf.f"
	    if (rrr1 <= 8.) {
#line 421 "slarrf.f"
		*sigma = lsigma;
#line 422 "slarrf.f"
		shift = 1;
#line 423 "slarrf.f"
		goto L100;
#line 424 "slarrf.f"
	    }
#line 425 "slarrf.f"
	} else if (indx == 2) {
#line 426 "slarrf.f"
	    tmp = (d__1 = work[*n], abs(d__1));
#line 427 "slarrf.f"
	    znm2 = 1.;
#line 428 "slarrf.f"
	    prod = 1.;
#line 429 "slarrf.f"
	    oldp = 1.;
#line 430 "slarrf.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 431 "slarrf.f"
		if (prod <= eps) {
#line 432 "slarrf.f"
		    prod = work[i__ + 1] * lplus[i__ + 1] / (work[i__] * 
			    lplus[i__]) * oldp;
#line 433 "slarrf.f"
		} else {
#line 434 "slarrf.f"
		    prod *= (d__1 = lplus[i__], abs(d__1));
#line 435 "slarrf.f"
		}
#line 436 "slarrf.f"
		oldp = prod;
/* Computing 2nd power */
#line 437 "slarrf.f"
		d__1 = prod;
#line 437 "slarrf.f"
		znm2 += d__1 * d__1;
/* Computing MAX */
#line 438 "slarrf.f"
		d__2 = tmp, d__3 = (d__1 = work[i__] * prod, abs(d__1));
#line 438 "slarrf.f"
		tmp = max(d__2,d__3);
#line 439 "slarrf.f"
/* L16: */
#line 439 "slarrf.f"
	    }
#line 440 "slarrf.f"
	    rrr2 = tmp / (*spdiam * sqrt(znm2));
#line 441 "slarrf.f"
	    if (rrr2 <= 8.) {
#line 442 "slarrf.f"
		*sigma = rsigma;
#line 443 "slarrf.f"
		shift = 2;
#line 444 "slarrf.f"
		goto L100;
#line 445 "slarrf.f"
	    }
#line 446 "slarrf.f"
	}
#line 447 "slarrf.f"
    }
#line 449 "slarrf.f"
L50:
#line 451 "slarrf.f"
    if (ktry < 1) {
/*        If we are here, both shifts failed also the RRR test. */
/*        Back off to the outside */
/* Computing MAX */
#line 454 "slarrf.f"
	d__1 = lsigma - ldelta, d__2 = lsigma - ldmax;
#line 454 "slarrf.f"
	lsigma = max(d__1,d__2);
/* Computing MIN */
#line 456 "slarrf.f"
	d__1 = rsigma + rdelta, d__2 = rsigma + rdmax;
#line 456 "slarrf.f"
	rsigma = min(d__1,d__2);
#line 458 "slarrf.f"
	ldelta *= 2.;
#line 459 "slarrf.f"
	rdelta *= 2.;
#line 460 "slarrf.f"
	++ktry;
#line 461 "slarrf.f"
	goto L5;
#line 462 "slarrf.f"
    } else {
/*        None of the representations investigated satisfied our */
/*        criteria. Take the best one we found. */
#line 465 "slarrf.f"
	if (smlgrowth < fail || nofail) {
#line 466 "slarrf.f"
	    lsigma = bestshift;
#line 467 "slarrf.f"
	    rsigma = bestshift;
#line 468 "slarrf.f"
	    forcer = TRUE_;
#line 469 "slarrf.f"
	    goto L5;
#line 470 "slarrf.f"
	} else {
#line 471 "slarrf.f"
	    *info = 1;
#line 472 "slarrf.f"
	    return 0;
#line 473 "slarrf.f"
	}
#line 474 "slarrf.f"
    }
#line 476 "slarrf.f"
L100:
#line 477 "slarrf.f"
    if (shift == 1) {
#line 478 "slarrf.f"
    } else if (shift == 2) {
/*        store new L and D back into DPLUS, LPLUS */
#line 480 "slarrf.f"
	scopy_(n, &work[1], &c__1, &dplus[1], &c__1);
#line 481 "slarrf.f"
	i__1 = *n - 1;
#line 481 "slarrf.f"
	scopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
#line 482 "slarrf.f"
    }
#line 484 "slarrf.f"
    return 0;

/*     End of SLARRF */

} /* slarrf_ */

