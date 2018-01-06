#line 1 "slarrb.f"
/* slarrb.f -- translated by f2c (version 20100827).
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

#line 1 "slarrb.f"
/* > \brief \b SLARRB provides limited bisection to locate eigenvalues for more accuracy. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARRB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARRB( N, D, LLD, IFIRST, ILAST, RTOL1, */
/*                          RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK, */
/*                          PIVMIN, SPDIAM, TWIST, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IFIRST, ILAST, INFO, N, OFFSET, TWIST */
/*       REAL               PIVMIN, RTOL1, RTOL2, SPDIAM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               D( * ), LLD( * ), W( * ), */
/*      $                   WERR( * ), WGAP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given the relatively robust representation(RRR) L D L^T, SLARRB */
/* > does "limited" bisection to refine the eigenvalues of L D L^T, */
/* > W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial */
/* > guesses for these eigenvalues are input in W, the corresponding estimate */
/* > of the error in these guesses and their gaps are input in WERR */
/* > and WGAP, respectively. During bisection, intervals */
/* > [left, right] are maintained by storing their mid-points and */
/* > semi-widths in the arrays W and WERR respectively. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LLD */
/* > \verbatim */
/* >          LLD is REAL array, dimension (N-1) */
/* >          The (N-1) elements L(i)*L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] IFIRST */
/* > \verbatim */
/* >          IFIRST is INTEGER */
/* >          The index of the first eigenvalue to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in] ILAST */
/* > \verbatim */
/* >          ILAST is INTEGER */
/* >          The index of the last eigenvalue to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL1 */
/* > \verbatim */
/* >          RTOL1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL2 */
/* > \verbatim */
/* >          RTOL2 is REAL */
/* >          Tolerance for the convergence of the bisection intervals. */
/* >          An interval [LEFT,RIGHT] has converged if */
/* >          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) ) */
/* >          where GAP is the (estimated) distance to the nearest */
/* >          eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[in] OFFSET */
/* > \verbatim */
/* >          OFFSET is INTEGER */
/* >          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET */
/* >          through ILAST-OFFSET elements of these arrays are to be used. */
/* > \endverbatim */
/* > */
/* > \param[in,out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are */
/* >          estimates of the eigenvalues of L D L^T indexed IFIRST throug */
/* >          ILAST. */
/* >          On output, these estimates are refined. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WGAP */
/* > \verbatim */
/* >          WGAP is REAL array, dimension (N-1) */
/* >          On input, the (estimated) gaps between consecutive */
/* >          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between */
/* >          eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST */
/* >          then WGAP(IFIRST-OFFSET) must be set to ZERO. */
/* >          On output, these gaps are refined. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WERR */
/* > \verbatim */
/* >          WERR is REAL array, dimension (N) */
/* >          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are */
/* >          the errors in the estimates of the corresponding elements in W. */
/* >          On output, these errors are refined. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (2*N) */
/* >          Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (2*N) */
/* >          Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
/* >          The minimum pivot in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[in] SPDIAM */
/* > \verbatim */
/* >          SPDIAM is REAL */
/* >          The spectral diameter of the matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TWIST */
/* > \verbatim */
/* >          TWIST is INTEGER */
/* >          The twist index for the twisted factorization that is used */
/* >          for the negcount. */
/* >          TWIST = N: Compute negcount from L D L^T - LAMBDA I = L+ D+ L+^T */
/* >          TWIST = 1: Compute negcount from L D L^T - LAMBDA I = U- D- U-^T */
/* >          TWIST = R: Compute negcount from L D L^T - LAMBDA I = N(r) D(r) N(r) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          Error flag. */
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
/* Subroutine */ int slarrb_(integer *n, doublereal *d__, doublereal *lld, 
	integer *ifirst, integer *ilast, doublereal *rtol1, doublereal *rtol2,
	 integer *offset, doublereal *w, doublereal *wgap, doublereal *werr, 
	doublereal *work, integer *iwork, doublereal *pivmin, doublereal *
	spdiam, integer *twist, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, k, r__, i1, ii, ip;
    static doublereal gap, mid, tmp, back, lgap, rgap, left;
    static integer iter, nint, prev, next;
    static doublereal cvrgd, right, width;
    extern integer slaneg_(integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *);
    static integer negcnt;
    static doublereal mnwdth;
    static integer olnint, maxitr;


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

#line 238 "slarrb.f"
    /* Parameter adjustments */
#line 238 "slarrb.f"
    --iwork;
#line 238 "slarrb.f"
    --work;
#line 238 "slarrb.f"
    --werr;
#line 238 "slarrb.f"
    --wgap;
#line 238 "slarrb.f"
    --w;
#line 238 "slarrb.f"
    --lld;
#line 238 "slarrb.f"
    --d__;
#line 238 "slarrb.f"

#line 238 "slarrb.f"
    /* Function Body */
#line 238 "slarrb.f"
    *info = 0;

#line 240 "slarrb.f"
    maxitr = (integer) ((log(*spdiam + *pivmin) - log(*pivmin)) / log(2.)) + 
	    2;
#line 242 "slarrb.f"
    mnwdth = *pivmin * 2.;

#line 244 "slarrb.f"
    r__ = *twist;
#line 245 "slarrb.f"
    if (r__ < 1 || r__ > *n) {
#line 245 "slarrb.f"
	r__ = *n;
#line 245 "slarrb.f"
    }

/*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ]. */
/*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while */
/*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 ) */
/*     for an unconverged interval is set to the index of the next unconverged */
/*     interval, and is -1 or 0 for a converged interval. Thus a linked */
/*     list of unconverged intervals is set up. */

#line 254 "slarrb.f"
    i1 = *ifirst;
/*     The number of unconverged intervals */
#line 256 "slarrb.f"
    nint = 0;
/*     The last unconverged interval found */
#line 258 "slarrb.f"
    prev = 0;
#line 260 "slarrb.f"
    rgap = wgap[i1 - *offset];
#line 261 "slarrb.f"
    i__1 = *ilast;
#line 261 "slarrb.f"
    for (i__ = i1; i__ <= i__1; ++i__) {
#line 262 "slarrb.f"
	k = i__ << 1;
#line 263 "slarrb.f"
	ii = i__ - *offset;
#line 264 "slarrb.f"
	left = w[ii] - werr[ii];
#line 265 "slarrb.f"
	right = w[ii] + werr[ii];
#line 266 "slarrb.f"
	lgap = rgap;
#line 267 "slarrb.f"
	rgap = wgap[ii];
#line 268 "slarrb.f"
	gap = min(lgap,rgap);
/*        Make sure that [LEFT,RIGHT] contains the desired eigenvalue */
/*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT */

/*        Do while( NEGCNT(LEFT).GT.I-1 ) */

#line 275 "slarrb.f"
	back = werr[ii];
#line 276 "slarrb.f"
L20:
#line 277 "slarrb.f"
	negcnt = slaneg_(n, &d__[1], &lld[1], &left, pivmin, &r__);
#line 278 "slarrb.f"
	if (negcnt > i__ - 1) {
#line 279 "slarrb.f"
	    left -= back;
#line 280 "slarrb.f"
	    back *= 2.;
#line 281 "slarrb.f"
	    goto L20;
#line 282 "slarrb.f"
	}

/*        Do while( NEGCNT(RIGHT).LT.I ) */
/*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT */

#line 287 "slarrb.f"
	back = werr[ii];
#line 288 "slarrb.f"
L50:
#line 290 "slarrb.f"
	negcnt = slaneg_(n, &d__[1], &lld[1], &right, pivmin, &r__);
#line 291 "slarrb.f"
	if (negcnt < i__) {
#line 292 "slarrb.f"
	    right += back;
#line 293 "slarrb.f"
	    back *= 2.;
#line 294 "slarrb.f"
	    goto L50;
#line 295 "slarrb.f"
	}
#line 296 "slarrb.f"
	width = (d__1 = left - right, abs(d__1)) * .5;
/* Computing MAX */
#line 297 "slarrb.f"
	d__1 = abs(left), d__2 = abs(right);
#line 297 "slarrb.f"
	tmp = max(d__1,d__2);
/* Computing MAX */
#line 298 "slarrb.f"
	d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
#line 298 "slarrb.f"
	cvrgd = max(d__1,d__2);
#line 299 "slarrb.f"
	if (width <= cvrgd || width <= mnwdth) {
/*           This interval has already converged and does not need refinement. */
/*           (Note that the gaps might change through refining the */
/*            eigenvalues, however, they can only get bigger.) */
/*           Remove it from the list. */
#line 304 "slarrb.f"
	    iwork[k - 1] = -1;
/*           Make sure that I1 always points to the first unconverged interval */
#line 306 "slarrb.f"
	    if (i__ == i1 && i__ < *ilast) {
#line 306 "slarrb.f"
		i1 = i__ + 1;
#line 306 "slarrb.f"
	    }
#line 307 "slarrb.f"
	    if (prev >= i1 && i__ <= *ilast) {
#line 307 "slarrb.f"
		iwork[(prev << 1) - 1] = i__ + 1;
#line 307 "slarrb.f"
	    }
#line 308 "slarrb.f"
	} else {
/*           unconverged interval found */
#line 310 "slarrb.f"
	    prev = i__;
#line 311 "slarrb.f"
	    ++nint;
#line 312 "slarrb.f"
	    iwork[k - 1] = i__ + 1;
#line 313 "slarrb.f"
	    iwork[k] = negcnt;
#line 314 "slarrb.f"
	}
#line 315 "slarrb.f"
	work[k - 1] = left;
#line 316 "slarrb.f"
	work[k] = right;
#line 317 "slarrb.f"
/* L75: */
#line 317 "slarrb.f"
    }

/*     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals */
/*     and while (ITER.LT.MAXITR) */

#line 323 "slarrb.f"
    iter = 0;
#line 324 "slarrb.f"
L80:
#line 325 "slarrb.f"
    prev = i1 - 1;
#line 326 "slarrb.f"
    i__ = i1;
#line 327 "slarrb.f"
    olnint = nint;
#line 329 "slarrb.f"
    i__1 = olnint;
#line 329 "slarrb.f"
    for (ip = 1; ip <= i__1; ++ip) {
#line 330 "slarrb.f"
	k = i__ << 1;
#line 331 "slarrb.f"
	ii = i__ - *offset;
#line 332 "slarrb.f"
	rgap = wgap[ii];
#line 333 "slarrb.f"
	lgap = rgap;
#line 334 "slarrb.f"
	if (ii > 1) {
#line 334 "slarrb.f"
	    lgap = wgap[ii - 1];
#line 334 "slarrb.f"
	}
#line 335 "slarrb.f"
	gap = min(lgap,rgap);
#line 336 "slarrb.f"
	next = iwork[k - 1];
#line 337 "slarrb.f"
	left = work[k - 1];
#line 338 "slarrb.f"
	right = work[k];
#line 339 "slarrb.f"
	mid = (left + right) * .5;
/*        semiwidth of interval */
#line 342 "slarrb.f"
	width = right - mid;
/* Computing MAX */
#line 343 "slarrb.f"
	d__1 = abs(left), d__2 = abs(right);
#line 343 "slarrb.f"
	tmp = max(d__1,d__2);
/* Computing MAX */
#line 344 "slarrb.f"
	d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
#line 344 "slarrb.f"
	cvrgd = max(d__1,d__2);
#line 345 "slarrb.f"
	if (width <= cvrgd || width <= mnwdth || iter == maxitr) {
/*           reduce number of unconverged intervals */
#line 348 "slarrb.f"
	    --nint;
/*           Mark interval as converged. */
#line 350 "slarrb.f"
	    iwork[k - 1] = 0;
#line 351 "slarrb.f"
	    if (i1 == i__) {
#line 352 "slarrb.f"
		i1 = next;
#line 353 "slarrb.f"
	    } else {
/*              Prev holds the last unconverged interval previously examined */
#line 355 "slarrb.f"
		if (prev >= i1) {
#line 355 "slarrb.f"
		    iwork[(prev << 1) - 1] = next;
#line 355 "slarrb.f"
		}
#line 356 "slarrb.f"
	    }
#line 357 "slarrb.f"
	    i__ = next;
#line 358 "slarrb.f"
	    goto L100;
#line 359 "slarrb.f"
	}
#line 360 "slarrb.f"
	prev = i__;

/*        Perform one bisection step */

#line 364 "slarrb.f"
	negcnt = slaneg_(n, &d__[1], &lld[1], &mid, pivmin, &r__);
#line 365 "slarrb.f"
	if (negcnt <= i__ - 1) {
#line 366 "slarrb.f"
	    work[k - 1] = mid;
#line 367 "slarrb.f"
	} else {
#line 368 "slarrb.f"
	    work[k] = mid;
#line 369 "slarrb.f"
	}
#line 370 "slarrb.f"
	i__ = next;
#line 371 "slarrb.f"
L100:
#line 371 "slarrb.f"
	;
#line 371 "slarrb.f"
    }
#line 372 "slarrb.f"
    ++iter;
/*     do another loop if there are still unconverged intervals */
/*     However, in the last iteration, all intervals are accepted */
/*     since this is the best we can do. */
#line 376 "slarrb.f"
    if (nint > 0 && iter <= maxitr) {
#line 376 "slarrb.f"
	goto L80;
#line 376 "slarrb.f"
    }


/*     At this point, all the intervals have converged */
#line 380 "slarrb.f"
    i__1 = *ilast;
#line 380 "slarrb.f"
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
#line 381 "slarrb.f"
	k = i__ << 1;
#line 382 "slarrb.f"
	ii = i__ - *offset;
/*        All intervals marked by '0' have been refined. */
#line 384 "slarrb.f"
	if (iwork[k - 1] == 0) {
#line 385 "slarrb.f"
	    w[ii] = (work[k - 1] + work[k]) * .5;
#line 386 "slarrb.f"
	    werr[ii] = work[k] - w[ii];
#line 387 "slarrb.f"
	}
#line 388 "slarrb.f"
/* L110: */
#line 388 "slarrb.f"
    }

#line 390 "slarrb.f"
    i__1 = *ilast;
#line 390 "slarrb.f"
    for (i__ = *ifirst + 1; i__ <= i__1; ++i__) {
#line 391 "slarrb.f"
	k = i__ << 1;
#line 392 "slarrb.f"
	ii = i__ - *offset;
/* Computing MAX */
#line 393 "slarrb.f"
	d__1 = 0., d__2 = w[ii] - werr[ii] - w[ii - 1] - werr[ii - 1];
#line 393 "slarrb.f"
	wgap[ii - 1] = max(d__1,d__2);
#line 395 "slarrb.f"
/* L111: */
#line 395 "slarrb.f"
    }
#line 397 "slarrb.f"
    return 0;

/*     End of SLARRB */

} /* slarrb_ */

