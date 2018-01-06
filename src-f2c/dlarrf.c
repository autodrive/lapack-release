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
/* >          The order of the matrix (subblock, if the matrix splitted). */
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

/* > \date September 2012 */

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
#line 242 "dlarrf.f"
    fact = 2.;
#line 243 "dlarrf.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 244 "dlarrf.f"
    shift = 0;
#line 245 "dlarrf.f"
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
#line 260 "dlarrf.f"
    nofail = TRUE_;

/*     Compute the average gap length of the cluster */
#line 264 "dlarrf.f"
    clwdth = (d__1 = w[*clend] - w[*clstrt], abs(d__1)) + werr[*clend] + werr[
	    *clstrt];
#line 265 "dlarrf.f"
    avgap = clwdth / (doublereal) (*clend - *clstrt);
#line 266 "dlarrf.f"
    mingap = min(*clgapl,*clgapr);
/*     Initial values for shifts to both ends of cluster */
/* Computing MIN */
#line 268 "dlarrf.f"
    d__1 = w[*clstrt], d__2 = w[*clend];
#line 268 "dlarrf.f"
    lsigma = min(d__1,d__2) - werr[*clstrt];
/* Computing MAX */
#line 269 "dlarrf.f"
    d__1 = w[*clstrt], d__2 = w[*clend];
#line 269 "dlarrf.f"
    rsigma = max(d__1,d__2) + werr[*clend];
/*     Use a small fudge to make sure that we really shift to the outside */
#line 272 "dlarrf.f"
    lsigma -= abs(lsigma) * 4. * eps;
#line 273 "dlarrf.f"
    rsigma += abs(rsigma) * 4. * eps;
/*     Compute upper bounds for how much to back off the initial shifts */
#line 276 "dlarrf.f"
    ldmax = mingap * .25 + *pivmin * 2.;
#line 277 "dlarrf.f"
    rdmax = mingap * .25 + *pivmin * 2.;
/* Computing MAX */
#line 279 "dlarrf.f"
    d__1 = avgap, d__2 = wgap[*clstrt];
#line 279 "dlarrf.f"
    ldelta = max(d__1,d__2) / fact;
/* Computing MAX */
#line 280 "dlarrf.f"
    d__1 = avgap, d__2 = wgap[*clend - 1];
#line 280 "dlarrf.f"
    rdelta = max(d__1,d__2) / fact;

/*     Initialize the record of the best representation found */

#line 284 "dlarrf.f"
    s = dlamch_("S", (ftnlen)1);
#line 285 "dlarrf.f"
    smlgrowth = 1. / s;
#line 286 "dlarrf.f"
    fail = (doublereal) (*n - 1) * mingap / (*spdiam * eps);
#line 287 "dlarrf.f"
    fail2 = (doublereal) (*n - 1) * mingap / (*spdiam * sqrt(eps));
#line 288 "dlarrf.f"
    bestshift = lsigma;

/*     while (KTRY <= KTRYMAX) */
#line 291 "dlarrf.f"
    ktry = 0;
#line 292 "dlarrf.f"
    growthbound = *spdiam * 8.;
#line 294 "dlarrf.f"
L5:
#line 295 "dlarrf.f"
    sawnan1 = FALSE_;
#line 296 "dlarrf.f"
    sawnan2 = FALSE_;
/*     Ensure that we do not back off too much of the initial shifts */
#line 298 "dlarrf.f"
    ldelta = min(ldmax,ldelta);
#line 299 "dlarrf.f"
    rdelta = min(rdmax,rdelta);
/*     Compute the element growth when shifting to both ends of the cluster */
/*     accept the shift if there is no element growth at one of the two ends */
/*     Left end */
#line 305 "dlarrf.f"
    s = -lsigma;
#line 306 "dlarrf.f"
    dplus[1] = d__[1] + s;
#line 307 "dlarrf.f"
    if (abs(dplus[1]) < *pivmin) {
#line 308 "dlarrf.f"
	dplus[1] = -(*pivmin);
/*        Need to set SAWNAN1 because refined RRR test should not be used */
/*        in this case */
#line 311 "dlarrf.f"
	sawnan1 = TRUE_;
#line 312 "dlarrf.f"
    }
#line 313 "dlarrf.f"
    max1 = abs(dplus[1]);
#line 314 "dlarrf.f"
    i__1 = *n - 1;
#line 314 "dlarrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 315 "dlarrf.f"
	lplus[i__] = ld[i__] / dplus[i__];
#line 316 "dlarrf.f"
	s = s * lplus[i__] * l[i__] - lsigma;
#line 317 "dlarrf.f"
	dplus[i__ + 1] = d__[i__ + 1] + s;
#line 318 "dlarrf.f"
	if ((d__1 = dplus[i__ + 1], abs(d__1)) < *pivmin) {
#line 319 "dlarrf.f"
	    dplus[i__ + 1] = -(*pivmin);
/*           Need to set SAWNAN1 because refined RRR test should not be used */
/*           in this case */
#line 322 "dlarrf.f"
	    sawnan1 = TRUE_;
#line 323 "dlarrf.f"
	}
/* Computing MAX */
#line 324 "dlarrf.f"
	d__2 = max1, d__3 = (d__1 = dplus[i__ + 1], abs(d__1));
#line 324 "dlarrf.f"
	max1 = max(d__2,d__3);
#line 325 "dlarrf.f"
/* L6: */
#line 325 "dlarrf.f"
    }
#line 326 "dlarrf.f"
    sawnan1 = sawnan1 || disnan_(&max1);
#line 328 "dlarrf.f"
    if (forcer || max1 <= growthbound && ! sawnan1) {
#line 330 "dlarrf.f"
	*sigma = lsigma;
#line 331 "dlarrf.f"
	shift = 1;
#line 332 "dlarrf.f"
	goto L100;
#line 333 "dlarrf.f"
    }
/*     Right end */
#line 336 "dlarrf.f"
    s = -rsigma;
#line 337 "dlarrf.f"
    work[1] = d__[1] + s;
#line 338 "dlarrf.f"
    if (abs(work[1]) < *pivmin) {
#line 339 "dlarrf.f"
	work[1] = -(*pivmin);
/*        Need to set SAWNAN2 because refined RRR test should not be used */
/*        in this case */
#line 342 "dlarrf.f"
	sawnan2 = TRUE_;
#line 343 "dlarrf.f"
    }
#line 344 "dlarrf.f"
    max2 = abs(work[1]);
#line 345 "dlarrf.f"
    i__1 = *n - 1;
#line 345 "dlarrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 346 "dlarrf.f"
	work[*n + i__] = ld[i__] / work[i__];
#line 347 "dlarrf.f"
	s = s * work[*n + i__] * l[i__] - rsigma;
#line 348 "dlarrf.f"
	work[i__ + 1] = d__[i__ + 1] + s;
#line 349 "dlarrf.f"
	if ((d__1 = work[i__ + 1], abs(d__1)) < *pivmin) {
#line 350 "dlarrf.f"
	    work[i__ + 1] = -(*pivmin);
/*           Need to set SAWNAN2 because refined RRR test should not be used */
/*           in this case */
#line 353 "dlarrf.f"
	    sawnan2 = TRUE_;
#line 354 "dlarrf.f"
	}
/* Computing MAX */
#line 355 "dlarrf.f"
	d__2 = max2, d__3 = (d__1 = work[i__ + 1], abs(d__1));
#line 355 "dlarrf.f"
	max2 = max(d__2,d__3);
#line 356 "dlarrf.f"
/* L7: */
#line 356 "dlarrf.f"
    }
#line 357 "dlarrf.f"
    sawnan2 = sawnan2 || disnan_(&max2);
#line 359 "dlarrf.f"
    if (forcer || max2 <= growthbound && ! sawnan2) {
#line 361 "dlarrf.f"
	*sigma = rsigma;
#line 362 "dlarrf.f"
	shift = 2;
#line 363 "dlarrf.f"
	goto L100;
#line 364 "dlarrf.f"
    }
/*     If we are at this point, both shifts led to too much element growth */
/*     Record the better of the two shifts (provided it didn't lead to NaN) */
#line 368 "dlarrf.f"
    if (sawnan1 && sawnan2) {
/*        both MAX1 and MAX2 are NaN */
#line 370 "dlarrf.f"
	goto L50;
#line 371 "dlarrf.f"
    } else {
#line 372 "dlarrf.f"
	if (! sawnan1) {
#line 373 "dlarrf.f"
	    indx = 1;
#line 374 "dlarrf.f"
	    if (max1 <= smlgrowth) {
#line 375 "dlarrf.f"
		smlgrowth = max1;
#line 376 "dlarrf.f"
		bestshift = lsigma;
#line 377 "dlarrf.f"
	    }
#line 378 "dlarrf.f"
	}
#line 379 "dlarrf.f"
	if (! sawnan2) {
#line 380 "dlarrf.f"
	    if (sawnan1 || max2 <= max1) {
#line 380 "dlarrf.f"
		indx = 2;
#line 380 "dlarrf.f"
	    }
#line 381 "dlarrf.f"
	    if (max2 <= smlgrowth) {
#line 382 "dlarrf.f"
		smlgrowth = max2;
#line 383 "dlarrf.f"
		bestshift = rsigma;
#line 384 "dlarrf.f"
	    }
#line 385 "dlarrf.f"
	}
#line 386 "dlarrf.f"
    }
/*     If we are here, both the left and the right shift led to */
/*     element growth. If the element growth is moderate, then */
/*     we may still accept the representation, if it passes a */
/*     refined test for RRR. This test supposes that no NaN occurred. */
/*     Moreover, we use the refined RRR test only for isolated clusters. */
#line 393 "dlarrf.f"
    if (clwdth < mingap / 128. && min(max1,max2) < fail2 && ! sawnan1 && ! 
	    sawnan2) {
#line 396 "dlarrf.f"
	dorrr1 = TRUE_;
#line 397 "dlarrf.f"
    } else {
#line 398 "dlarrf.f"
	dorrr1 = FALSE_;
#line 399 "dlarrf.f"
    }
#line 400 "dlarrf.f"
    tryrrr1 = TRUE_;
#line 401 "dlarrf.f"
    if (tryrrr1 && dorrr1) {
#line 402 "dlarrf.f"
	if (indx == 1) {
#line 403 "dlarrf.f"
	    tmp = (d__1 = dplus[*n], abs(d__1));
#line 404 "dlarrf.f"
	    znm2 = 1.;
#line 405 "dlarrf.f"
	    prod = 1.;
#line 406 "dlarrf.f"
	    oldp = 1.;
#line 407 "dlarrf.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 408 "dlarrf.f"
		if (prod <= eps) {
#line 409 "dlarrf.f"
		    prod = dplus[i__ + 1] * work[*n + i__ + 1] / (dplus[i__] *
			     work[*n + i__]) * oldp;
#line 411 "dlarrf.f"
		} else {
#line 412 "dlarrf.f"
		    prod *= (d__1 = work[*n + i__], abs(d__1));
#line 413 "dlarrf.f"
		}
#line 414 "dlarrf.f"
		oldp = prod;
/* Computing 2nd power */
#line 415 "dlarrf.f"
		d__1 = prod;
#line 415 "dlarrf.f"
		znm2 += d__1 * d__1;
/* Computing MAX */
#line 416 "dlarrf.f"
		d__2 = tmp, d__3 = (d__1 = dplus[i__] * prod, abs(d__1));
#line 416 "dlarrf.f"
		tmp = max(d__2,d__3);
#line 417 "dlarrf.f"
/* L15: */
#line 417 "dlarrf.f"
	    }
#line 418 "dlarrf.f"
	    rrr1 = tmp / (*spdiam * sqrt(znm2));
#line 419 "dlarrf.f"
	    if (rrr1 <= 8.) {
#line 420 "dlarrf.f"
		*sigma = lsigma;
#line 421 "dlarrf.f"
		shift = 1;
#line 422 "dlarrf.f"
		goto L100;
#line 423 "dlarrf.f"
	    }
#line 424 "dlarrf.f"
	} else if (indx == 2) {
#line 425 "dlarrf.f"
	    tmp = (d__1 = work[*n], abs(d__1));
#line 426 "dlarrf.f"
	    znm2 = 1.;
#line 427 "dlarrf.f"
	    prod = 1.;
#line 428 "dlarrf.f"
	    oldp = 1.;
#line 429 "dlarrf.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 430 "dlarrf.f"
		if (prod <= eps) {
#line 431 "dlarrf.f"
		    prod = work[i__ + 1] * lplus[i__ + 1] / (work[i__] * 
			    lplus[i__]) * oldp;
#line 432 "dlarrf.f"
		} else {
#line 433 "dlarrf.f"
		    prod *= (d__1 = lplus[i__], abs(d__1));
#line 434 "dlarrf.f"
		}
#line 435 "dlarrf.f"
		oldp = prod;
/* Computing 2nd power */
#line 436 "dlarrf.f"
		d__1 = prod;
#line 436 "dlarrf.f"
		znm2 += d__1 * d__1;
/* Computing MAX */
#line 437 "dlarrf.f"
		d__2 = tmp, d__3 = (d__1 = work[i__] * prod, abs(d__1));
#line 437 "dlarrf.f"
		tmp = max(d__2,d__3);
#line 438 "dlarrf.f"
/* L16: */
#line 438 "dlarrf.f"
	    }
#line 439 "dlarrf.f"
	    rrr2 = tmp / (*spdiam * sqrt(znm2));
#line 440 "dlarrf.f"
	    if (rrr2 <= 8.) {
#line 441 "dlarrf.f"
		*sigma = rsigma;
#line 442 "dlarrf.f"
		shift = 2;
#line 443 "dlarrf.f"
		goto L100;
#line 444 "dlarrf.f"
	    }
#line 445 "dlarrf.f"
	}
#line 446 "dlarrf.f"
    }
#line 448 "dlarrf.f"
L50:
#line 450 "dlarrf.f"
    if (ktry < 1) {
/*        If we are here, both shifts failed also the RRR test. */
/*        Back off to the outside */
/* Computing MAX */
#line 453 "dlarrf.f"
	d__1 = lsigma - ldelta, d__2 = lsigma - ldmax;
#line 453 "dlarrf.f"
	lsigma = max(d__1,d__2);
/* Computing MIN */
#line 455 "dlarrf.f"
	d__1 = rsigma + rdelta, d__2 = rsigma + rdmax;
#line 455 "dlarrf.f"
	rsigma = min(d__1,d__2);
#line 457 "dlarrf.f"
	ldelta *= 2.;
#line 458 "dlarrf.f"
	rdelta *= 2.;
#line 459 "dlarrf.f"
	++ktry;
#line 460 "dlarrf.f"
	goto L5;
#line 461 "dlarrf.f"
    } else {
/*        None of the representations investigated satisfied our */
/*        criteria. Take the best one we found. */
#line 464 "dlarrf.f"
	if (smlgrowth < fail || nofail) {
#line 465 "dlarrf.f"
	    lsigma = bestshift;
#line 466 "dlarrf.f"
	    rsigma = bestshift;
#line 467 "dlarrf.f"
	    forcer = TRUE_;
#line 468 "dlarrf.f"
	    goto L5;
#line 469 "dlarrf.f"
	} else {
#line 470 "dlarrf.f"
	    *info = 1;
#line 471 "dlarrf.f"
	    return 0;
#line 472 "dlarrf.f"
	}
#line 473 "dlarrf.f"
    }
#line 475 "dlarrf.f"
L100:
#line 476 "dlarrf.f"
    if (shift == 1) {
#line 477 "dlarrf.f"
    } else if (shift == 2) {
/*        store new L and D back into DPLUS, LPLUS */
#line 479 "dlarrf.f"
	dcopy_(n, &work[1], &c__1, &dplus[1], &c__1);
#line 480 "dlarrf.f"
	i__1 = *n - 1;
#line 480 "dlarrf.f"
	dcopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
#line 481 "dlarrf.f"
    }
#line 483 "dlarrf.f"
    return 0;

/*     End of DLARRF */

} /* dlarrf_ */

