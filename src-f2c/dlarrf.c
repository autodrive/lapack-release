#line 1 "dlarrf.f"
/* dlarrf.f -- translated by f2c (version 20100827).
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

#line 1 "dlarrf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLARRF finds a new relatively robust representation such that at least one of the eigenvalues i
s relatively isolated. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRF( N, D, L, LD, CLSTRT, CLEND, */
/*                          W, WGAP, WERR, */
/*                          SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA, */
/*                          DPLUS, LPLUS, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            CLSTRT, CLEND, INFO, N */
/*       DOUBLE PRECISION   CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), DPLUS( * ), L( * ), LD( * ), */
/*      $          LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given the initial representation L D L^T and its cluster of close */
/* > eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ... */
/* > W( CLEND ), DLARRF finds a new relatively robust representation */
/* > L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the */
/* > eigenvalues of L(+) D(+) L(+)^T is relatively isolated. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix (subblock, if the matrix split). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (N-1) subdiagonal elements of the unit bidiagonal */
/* >          matrix L. */
/* > \endverbatim */
/* > */
/* > \param[in] LD */
/* > \verbatim */
/* >          LD is DOUBLE PRECISION array, dimension (N-1) */
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
/* >          W is DOUBLE PRECISION array, dimension */
/* >          dimension is >=  (CLEND-CLSTRT+1) */
/* >          The eigenvalue APPROXIMATIONS of L D L^T in ascending order. */
/* >          W( CLSTRT ) through W( CLEND ) form the cluster of relatively */
/* >          close eigenalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WGAP */
/* > \verbatim */
/* >          WGAP is DOUBLE PRECISION array, dimension */
/* >          dimension is >=  (CLEND-CLSTRT+1) */
/* >          The separation from the right neighbor eigenvalue in W. */
/* > \endverbatim */
/* > */
/* > \param[in] WERR */
/* > \verbatim */
/* >          WERR is DOUBLE PRECISION array, dimension */
/* >          dimension is  >=  (CLEND-CLSTRT+1) */
/* >          WERR contain the semiwidth of the uncertainty */
/* >          interval of the corresponding eigenvalue APPROXIMATION in W */
/* > \endverbatim */
/* > */
/* > \param[in] SPDIAM */
/* > \verbatim */
/* >          SPDIAM is DOUBLE PRECISION */
/* >          estimate of the spectral diameter obtained from the */
/* >          Gerschgorin intervals */
/* > \endverbatim */
/* > */
/* > \param[in] CLGAPL */
/* > \verbatim */
/* >          CLGAPL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] CLGAPR */
/* > \verbatim */
/* >          CLGAPR is DOUBLE PRECISION */
/* >          absolute gap on each end of the cluster. */
/* >          Set by the calling routine to protect against shifts too close */
/* >          to eigenvalues outside the cluster. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
/* >          The minimum pivot allowed in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >          The shift used to form L(+) D(+) L(+)^T. */
/* > \endverbatim */
/* > */
/* > \param[out] DPLUS */
/* > \verbatim */
/* >          DPLUS is DOUBLE PRECISION array, dimension (N) */
/* >          The N diagonal elements of the diagonal matrix D(+). */
/* > \endverbatim */
/* > */
/* > \param[out] LPLUS */
/* > \verbatim */
/* >          LPLUS is DOUBLE PRECISION array, dimension (N-1) */
/* >          The first (N-1) elements of LPLUS contain the subdiagonal */
/* >          elements of the unit bidiagonal matrix L(+). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \date June 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dlarrf_(integer *n, doublereal *d__, doublereal *l, 
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
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical dorrr1;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal ldelta;
    static logical nofail;
    static doublereal mingap, lsigma, rdelta;
    extern logical disnan_(doublereal *);
    static logical forcer;
    static doublereal rsigma, clwdth;
    static logical sawnan1, sawnan2, tryrrr1;


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 241 "dlarrf.f"
    /* Parameter adjustments */
#line 241 "dlarrf.f"
    --work;
#line 241 "dlarrf.f"
    --lplus;
#line 241 "dlarrf.f"
    --dplus;
#line 241 "dlarrf.f"
    --werr;
#line 241 "dlarrf.f"
    --wgap;
#line 241 "dlarrf.f"
    --w;
#line 241 "dlarrf.f"
    --ld;
#line 241 "dlarrf.f"
    --l;
#line 241 "dlarrf.f"
    --d__;
#line 241 "dlarrf.f"

#line 241 "dlarrf.f"
    /* Function Body */
#line 241 "dlarrf.f"
    *info = 0;

/*     Quick return if possible */

#line 245 "dlarrf.f"
    if (*n <= 0) {
#line 246 "dlarrf.f"
	return 0;
#line 247 "dlarrf.f"
    }

#line 249 "dlarrf.f"
    fact = 2.;
#line 250 "dlarrf.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 251 "dlarrf.f"
    shift = 0;
#line 252 "dlarrf.f"
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
#line 268 "dlarrf.f"
    nofail = FALSE_;

/*     Compute the average gap length of the cluster */
#line 272 "dlarrf.f"
    clwdth = (d__1 = w[*clend] - w[*clstrt], abs(d__1)) + werr[*clend] + werr[
	    *clstrt];
#line 273 "dlarrf.f"
    avgap = clwdth / (doublereal) (*clend - *clstrt);
#line 274 "dlarrf.f"
    mingap = min(*clgapl,*clgapr);
/*     Initial values for shifts to both ends of cluster */
/* Computing MIN */
#line 276 "dlarrf.f"
    d__1 = w[*clstrt], d__2 = w[*clend];
#line 276 "dlarrf.f"
    lsigma = min(d__1,d__2) - werr[*clstrt];
/* Computing MAX */
#line 277 "dlarrf.f"
    d__1 = w[*clstrt], d__2 = w[*clend];
#line 277 "dlarrf.f"
    rsigma = max(d__1,d__2) + werr[*clend];
/*     Use a small fudge to make sure that we really shift to the outside */
#line 280 "dlarrf.f"
    lsigma -= abs(lsigma) * 4. * eps;
#line 281 "dlarrf.f"
    rsigma += abs(rsigma) * 4. * eps;
/*     Compute upper bounds for how much to back off the initial shifts */
#line 284 "dlarrf.f"
    ldmax = mingap * .25 + *pivmin * 2.;
#line 285 "dlarrf.f"
    rdmax = mingap * .25 + *pivmin * 2.;
/* Computing MAX */
#line 287 "dlarrf.f"
    d__1 = avgap, d__2 = wgap[*clstrt];
#line 287 "dlarrf.f"
    ldelta = max(d__1,d__2) / fact;
/* Computing MAX */
#line 288 "dlarrf.f"
    d__1 = avgap, d__2 = wgap[*clend - 1];
#line 288 "dlarrf.f"
    rdelta = max(d__1,d__2) / fact;

/*     Initialize the record of the best representation found */

#line 292 "dlarrf.f"
    s = dlamch_("S", (ftnlen)1);
#line 293 "dlarrf.f"
    smlgrowth = 1. / s;
#line 294 "dlarrf.f"
    fail = (doublereal) (*n - 1) * mingap / (*spdiam * eps);
#line 295 "dlarrf.f"
    fail2 = (doublereal) (*n - 1) * mingap / (*spdiam * sqrt(eps));
#line 296 "dlarrf.f"
    bestshift = lsigma;

/*     while (KTRY <= KTRYMAX) */
#line 299 "dlarrf.f"
    ktry = 0;
#line 300 "dlarrf.f"
    growthbound = *spdiam * 8.;
#line 302 "dlarrf.f"
L5:
#line 303 "dlarrf.f"
    sawnan1 = FALSE_;
#line 304 "dlarrf.f"
    sawnan2 = FALSE_;
/*     Ensure that we do not back off too much of the initial shifts */
#line 306 "dlarrf.f"
    ldelta = min(ldmax,ldelta);
#line 307 "dlarrf.f"
    rdelta = min(rdmax,rdelta);
/*     Compute the element growth when shifting to both ends of the cluster */
/*     accept the shift if there is no element growth at one of the two ends */
/*     Left end */
#line 313 "dlarrf.f"
    s = -lsigma;
#line 314 "dlarrf.f"
    dplus[1] = d__[1] + s;
#line 315 "dlarrf.f"
    if (abs(dplus[1]) < *pivmin) {
#line 316 "dlarrf.f"
	dplus[1] = -(*pivmin);
/*        Need to set SAWNAN1 because refined RRR test should not be used */
/*        in this case */
#line 319 "dlarrf.f"
	sawnan1 = TRUE_;
#line 320 "dlarrf.f"
    }
#line 321 "dlarrf.f"
    max1 = abs(dplus[1]);
#line 322 "dlarrf.f"
    i__1 = *n - 1;
#line 322 "dlarrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 323 "dlarrf.f"
	lplus[i__] = ld[i__] / dplus[i__];
#line 324 "dlarrf.f"
	s = s * lplus[i__] * l[i__] - lsigma;
#line 325 "dlarrf.f"
	dplus[i__ + 1] = d__[i__ + 1] + s;
#line 326 "dlarrf.f"
	if ((d__1 = dplus[i__ + 1], abs(d__1)) < *pivmin) {
#line 327 "dlarrf.f"
	    dplus[i__ + 1] = -(*pivmin);
/*           Need to set SAWNAN1 because refined RRR test should not be used */
/*           in this case */
#line 330 "dlarrf.f"
	    sawnan1 = TRUE_;
#line 331 "dlarrf.f"
	}
/* Computing MAX */
#line 332 "dlarrf.f"
	d__2 = max1, d__3 = (d__1 = dplus[i__ + 1], abs(d__1));
#line 332 "dlarrf.f"
	max1 = max(d__2,d__3);
#line 333 "dlarrf.f"
/* L6: */
#line 333 "dlarrf.f"
    }
#line 334 "dlarrf.f"
    sawnan1 = sawnan1 || disnan_(&max1);
#line 336 "dlarrf.f"
    if (forcer || max1 <= growthbound && ! sawnan1) {
#line 338 "dlarrf.f"
	*sigma = lsigma;
#line 339 "dlarrf.f"
	shift = 1;
#line 340 "dlarrf.f"
	goto L100;
#line 341 "dlarrf.f"
    }
/*     Right end */
#line 344 "dlarrf.f"
    s = -rsigma;
#line 345 "dlarrf.f"
    work[1] = d__[1] + s;
#line 346 "dlarrf.f"
    if (abs(work[1]) < *pivmin) {
#line 347 "dlarrf.f"
	work[1] = -(*pivmin);
/*        Need to set SAWNAN2 because refined RRR test should not be used */
/*        in this case */
#line 350 "dlarrf.f"
	sawnan2 = TRUE_;
#line 351 "dlarrf.f"
    }
#line 352 "dlarrf.f"
    max2 = abs(work[1]);
#line 353 "dlarrf.f"
    i__1 = *n - 1;
#line 353 "dlarrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 354 "dlarrf.f"
	work[*n + i__] = ld[i__] / work[i__];
#line 355 "dlarrf.f"
	s = s * work[*n + i__] * l[i__] - rsigma;
#line 356 "dlarrf.f"
	work[i__ + 1] = d__[i__ + 1] + s;
#line 357 "dlarrf.f"
	if ((d__1 = work[i__ + 1], abs(d__1)) < *pivmin) {
#line 358 "dlarrf.f"
	    work[i__ + 1] = -(*pivmin);
/*           Need to set SAWNAN2 because refined RRR test should not be used */
/*           in this case */
#line 361 "dlarrf.f"
	    sawnan2 = TRUE_;
#line 362 "dlarrf.f"
	}
/* Computing MAX */
#line 363 "dlarrf.f"
	d__2 = max2, d__3 = (d__1 = work[i__ + 1], abs(d__1));
#line 363 "dlarrf.f"
	max2 = max(d__2,d__3);
#line 364 "dlarrf.f"
/* L7: */
#line 364 "dlarrf.f"
    }
#line 365 "dlarrf.f"
    sawnan2 = sawnan2 || disnan_(&max2);
#line 367 "dlarrf.f"
    if (forcer || max2 <= growthbound && ! sawnan2) {
#line 369 "dlarrf.f"
	*sigma = rsigma;
#line 370 "dlarrf.f"
	shift = 2;
#line 371 "dlarrf.f"
	goto L100;
#line 372 "dlarrf.f"
    }
/*     If we are at this point, both shifts led to too much element growth */
/*     Record the better of the two shifts (provided it didn't lead to NaN) */
#line 376 "dlarrf.f"
    if (sawnan1 && sawnan2) {
/*        both MAX1 and MAX2 are NaN */
#line 378 "dlarrf.f"
	goto L50;
#line 379 "dlarrf.f"
    } else {
#line 380 "dlarrf.f"
	if (! sawnan1) {
#line 381 "dlarrf.f"
	    indx = 1;
#line 382 "dlarrf.f"
	    if (max1 <= smlgrowth) {
#line 383 "dlarrf.f"
		smlgrowth = max1;
#line 384 "dlarrf.f"
		bestshift = lsigma;
#line 385 "dlarrf.f"
	    }
#line 386 "dlarrf.f"
	}
#line 387 "dlarrf.f"
	if (! sawnan2) {
#line 388 "dlarrf.f"
	    if (sawnan1 || max2 <= max1) {
#line 388 "dlarrf.f"
		indx = 2;
#line 388 "dlarrf.f"
	    }
#line 389 "dlarrf.f"
	    if (max2 <= smlgrowth) {
#line 390 "dlarrf.f"
		smlgrowth = max2;
#line 391 "dlarrf.f"
		bestshift = rsigma;
#line 392 "dlarrf.f"
	    }
#line 393 "dlarrf.f"
	}
#line 394 "dlarrf.f"
    }
/*     If we are here, both the left and the right shift led to */
/*     element growth. If the element growth is moderate, then */
/*     we may still accept the representation, if it passes a */
/*     refined test for RRR. This test supposes that no NaN occurred. */
/*     Moreover, we use the refined RRR test only for isolated clusters. */
#line 401 "dlarrf.f"
    if (clwdth < mingap / 128. && min(max1,max2) < fail2 && ! sawnan1 && ! 
	    sawnan2) {
#line 404 "dlarrf.f"
	dorrr1 = TRUE_;
#line 405 "dlarrf.f"
    } else {
#line 406 "dlarrf.f"
	dorrr1 = FALSE_;
#line 407 "dlarrf.f"
    }
#line 408 "dlarrf.f"
    tryrrr1 = TRUE_;
#line 409 "dlarrf.f"
    if (tryrrr1 && dorrr1) {
#line 410 "dlarrf.f"
	if (indx == 1) {
#line 411 "dlarrf.f"
	    tmp = (d__1 = dplus[*n], abs(d__1));
#line 412 "dlarrf.f"
	    znm2 = 1.;
#line 413 "dlarrf.f"
	    prod = 1.;
#line 414 "dlarrf.f"
	    oldp = 1.;
#line 415 "dlarrf.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 416 "dlarrf.f"
		if (prod <= eps) {
#line 417 "dlarrf.f"
		    prod = dplus[i__ + 1] * work[*n + i__ + 1] / (dplus[i__] *
			     work[*n + i__]) * oldp;
#line 419 "dlarrf.f"
		} else {
#line 420 "dlarrf.f"
		    prod *= (d__1 = work[*n + i__], abs(d__1));
#line 421 "dlarrf.f"
		}
#line 422 "dlarrf.f"
		oldp = prod;
/* Computing 2nd power */
#line 423 "dlarrf.f"
		d__1 = prod;
#line 423 "dlarrf.f"
		znm2 += d__1 * d__1;
/* Computing MAX */
#line 424 "dlarrf.f"
		d__2 = tmp, d__3 = (d__1 = dplus[i__] * prod, abs(d__1));
#line 424 "dlarrf.f"
		tmp = max(d__2,d__3);
#line 425 "dlarrf.f"
/* L15: */
#line 425 "dlarrf.f"
	    }
#line 426 "dlarrf.f"
	    rrr1 = tmp / (*spdiam * sqrt(znm2));
#line 427 "dlarrf.f"
	    if (rrr1 <= 8.) {
#line 428 "dlarrf.f"
		*sigma = lsigma;
#line 429 "dlarrf.f"
		shift = 1;
#line 430 "dlarrf.f"
		goto L100;
#line 431 "dlarrf.f"
	    }
#line 432 "dlarrf.f"
	} else if (indx == 2) {
#line 433 "dlarrf.f"
	    tmp = (d__1 = work[*n], abs(d__1));
#line 434 "dlarrf.f"
	    znm2 = 1.;
#line 435 "dlarrf.f"
	    prod = 1.;
#line 436 "dlarrf.f"
	    oldp = 1.;
#line 437 "dlarrf.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 438 "dlarrf.f"
		if (prod <= eps) {
#line 439 "dlarrf.f"
		    prod = work[i__ + 1] * lplus[i__ + 1] / (work[i__] * 
			    lplus[i__]) * oldp;
#line 440 "dlarrf.f"
		} else {
#line 441 "dlarrf.f"
		    prod *= (d__1 = lplus[i__], abs(d__1));
#line 442 "dlarrf.f"
		}
#line 443 "dlarrf.f"
		oldp = prod;
/* Computing 2nd power */
#line 444 "dlarrf.f"
		d__1 = prod;
#line 444 "dlarrf.f"
		znm2 += d__1 * d__1;
/* Computing MAX */
#line 445 "dlarrf.f"
		d__2 = tmp, d__3 = (d__1 = work[i__] * prod, abs(d__1));
#line 445 "dlarrf.f"
		tmp = max(d__2,d__3);
#line 446 "dlarrf.f"
/* L16: */
#line 446 "dlarrf.f"
	    }
#line 447 "dlarrf.f"
	    rrr2 = tmp / (*spdiam * sqrt(znm2));
#line 448 "dlarrf.f"
	    if (rrr2 <= 8.) {
#line 449 "dlarrf.f"
		*sigma = rsigma;
#line 450 "dlarrf.f"
		shift = 2;
#line 451 "dlarrf.f"
		goto L100;
#line 452 "dlarrf.f"
	    }
#line 453 "dlarrf.f"
	}
#line 454 "dlarrf.f"
    }
#line 456 "dlarrf.f"
L50:
#line 458 "dlarrf.f"
    if (ktry < 1) {
/*        If we are here, both shifts failed also the RRR test. */
/*        Back off to the outside */
/* Computing MAX */
#line 461 "dlarrf.f"
	d__1 = lsigma - ldelta, d__2 = lsigma - ldmax;
#line 461 "dlarrf.f"
	lsigma = max(d__1,d__2);
/* Computing MIN */
#line 463 "dlarrf.f"
	d__1 = rsigma + rdelta, d__2 = rsigma + rdmax;
#line 463 "dlarrf.f"
	rsigma = min(d__1,d__2);
#line 465 "dlarrf.f"
	ldelta *= 2.;
#line 466 "dlarrf.f"
	rdelta *= 2.;
#line 467 "dlarrf.f"
	++ktry;
#line 468 "dlarrf.f"
	goto L5;
#line 469 "dlarrf.f"
    } else {
/*        None of the representations investigated satisfied our */
/*        criteria. Take the best one we found. */
#line 472 "dlarrf.f"
	if (smlgrowth < fail || nofail) {
#line 473 "dlarrf.f"
	    lsigma = bestshift;
#line 474 "dlarrf.f"
	    rsigma = bestshift;
#line 475 "dlarrf.f"
	    forcer = TRUE_;
#line 476 "dlarrf.f"
	    goto L5;
#line 477 "dlarrf.f"
	} else {
#line 478 "dlarrf.f"
	    *info = 1;
#line 479 "dlarrf.f"
	    return 0;
#line 480 "dlarrf.f"
	}
#line 481 "dlarrf.f"
    }
#line 483 "dlarrf.f"
L100:
#line 484 "dlarrf.f"
    if (shift == 1) {
#line 485 "dlarrf.f"
    } else if (shift == 2) {
/*        store new L and D back into DPLUS, LPLUS */
#line 487 "dlarrf.f"
	dcopy_(n, &work[1], &c__1, &dplus[1], &c__1);
#line 488 "dlarrf.f"
	i__1 = *n - 1;
#line 488 "dlarrf.f"
	dcopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
#line 489 "dlarrf.f"
    }
#line 491 "dlarrf.f"
    return 0;

/*     End of DLARRF */

} /* dlarrf_ */

