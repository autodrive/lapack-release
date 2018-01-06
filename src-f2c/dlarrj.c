#line 1 "dlarrj.f"
/* dlarrj.f -- translated by f2c (version 20100827).
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

#line 1 "dlarrj.f"
/* > \brief \b DLARRJ performs refinement of the initial estimates of the eigenvalues of the matrix T. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRJ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrj.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrj.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrj.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRJ( N, D, E2, IFIRST, ILAST, */
/*                          RTOL, OFFSET, W, WERR, WORK, IWORK, */
/*                          PIVMIN, SPDIAM, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IFIRST, ILAST, INFO, N, OFFSET */
/*       DOUBLE PRECISION   PIVMIN, RTOL, SPDIAM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E2( * ), W( * ), */
/*      $                   WERR( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given the initial eigenvalue approximations of T, DLARRJ */
/* > does  bisection to refine the eigenvalues of T, */
/* > W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial */
/* > guesses for these eigenvalues are input in W, the corresponding estimate */
/* > of the error in these guesses in WERR. During bisection, intervals */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The N diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] E2 */
/* > \verbatim */
/* >          E2 is DOUBLE PRECISION array, dimension (N-1) */
/* >          The Squares of the (N-1) subdiagonal elements of T. */
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
/* > \param[in] RTOL */
/* > \verbatim */
/* >          RTOL is DOUBLE PRECISION */
/* >          Tolerance for the convergence of the bisection intervals. */
/* >          An interval [LEFT,RIGHT] has converged if */
/* >          RIGHT-LEFT.LT.RTOL*MAX(|LEFT|,|RIGHT|). */
/* > \endverbatim */
/* > */
/* > \param[in] OFFSET */
/* > \verbatim */
/* >          OFFSET is INTEGER */
/* >          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET */
/* >          through ILAST-OFFSET elements of these arrays are to be used. */
/* > \endverbatim */
/* > */
/* > \param[in,out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are */
/* >          estimates of the eigenvalues of L D L^T indexed IFIRST through */
/* >          ILAST. */
/* >          On output, these estimates are refined. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WERR */
/* > \verbatim */
/* >          WERR is DOUBLE PRECISION array, dimension (N) */
/* >          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are */
/* >          the errors in the estimates of the corresponding elements in W. */
/* >          On output, these errors are refined. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
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
/* >          PIVMIN is DOUBLE PRECISION */
/* >          The minimum pivot in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[in] SPDIAM */
/* > \verbatim */
/* >          SPDIAM is DOUBLE PRECISION */
/* >          The spectral diameter of T. */
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

/* > \date June 2017 */

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
/* Subroutine */ int dlarrj_(integer *n, doublereal *d__, doublereal *e2, 
	integer *ifirst, integer *ilast, doublereal *rtol, integer *offset, 
	doublereal *w, doublereal *werr, doublereal *work, integer *iwork, 
	doublereal *pivmin, doublereal *spdiam, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, j, k, p;
    static doublereal s;
    static integer i1, i2, ii;
    static doublereal fac, mid;
    static integer cnt;
    static doublereal tmp, left;
    static integer iter, nint, prev, next, savi1;
    static doublereal right, width, dplus;
    static integer olnint, maxitr;


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */

/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 205 "dlarrj.f"
    /* Parameter adjustments */
#line 205 "dlarrj.f"
    --iwork;
#line 205 "dlarrj.f"
    --work;
#line 205 "dlarrj.f"
    --werr;
#line 205 "dlarrj.f"
    --w;
#line 205 "dlarrj.f"
    --e2;
#line 205 "dlarrj.f"
    --d__;
#line 205 "dlarrj.f"

#line 205 "dlarrj.f"
    /* Function Body */
#line 205 "dlarrj.f"
    *info = 0;

/*     Quick return if possible */

#line 209 "dlarrj.f"
    if (*n <= 0) {
#line 210 "dlarrj.f"
	return 0;
#line 211 "dlarrj.f"
    }

#line 213 "dlarrj.f"
    maxitr = (integer) ((log(*spdiam + *pivmin) - log(*pivmin)) / log(2.)) + 
	    2;

/*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ]. */
/*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while */
/*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 ) */
/*     for an unconverged interval is set to the index of the next unconverged */
/*     interval, and is -1 or 0 for a converged interval. Thus a linked */
/*     list of unconverged intervals is set up. */

#line 224 "dlarrj.f"
    i1 = *ifirst;
#line 225 "dlarrj.f"
    i2 = *ilast;
/*     The number of unconverged intervals */
#line 227 "dlarrj.f"
    nint = 0;
/*     The last unconverged interval found */
#line 229 "dlarrj.f"
    prev = 0;
#line 230 "dlarrj.f"
    i__1 = i2;
#line 230 "dlarrj.f"
    for (i__ = i1; i__ <= i__1; ++i__) {
#line 231 "dlarrj.f"
	k = i__ << 1;
#line 232 "dlarrj.f"
	ii = i__ - *offset;
#line 233 "dlarrj.f"
	left = w[ii] - werr[ii];
#line 234 "dlarrj.f"
	mid = w[ii];
#line 235 "dlarrj.f"
	right = w[ii] + werr[ii];
#line 236 "dlarrj.f"
	width = right - mid;
/* Computing MAX */
#line 237 "dlarrj.f"
	d__1 = abs(left), d__2 = abs(right);
#line 237 "dlarrj.f"
	tmp = max(d__1,d__2);
/*        The following test prevents the test of converged intervals */
#line 240 "dlarrj.f"
	if (width < *rtol * tmp) {
/*           This interval has already converged and does not need refinement. */
/*           (Note that the gaps might change through refining the */
/*            eigenvalues, however, they can only get bigger.) */
/*           Remove it from the list. */
#line 245 "dlarrj.f"
	    iwork[k - 1] = -1;
/*           Make sure that I1 always points to the first unconverged interval */
#line 247 "dlarrj.f"
	    if (i__ == i1 && i__ < i2) {
#line 247 "dlarrj.f"
		i1 = i__ + 1;
#line 247 "dlarrj.f"
	    }
#line 248 "dlarrj.f"
	    if (prev >= i1 && i__ <= i2) {
#line 248 "dlarrj.f"
		iwork[(prev << 1) - 1] = i__ + 1;
#line 248 "dlarrj.f"
	    }
#line 249 "dlarrj.f"
	} else {
/*           unconverged interval found */
#line 251 "dlarrj.f"
	    prev = i__;
/*           Make sure that [LEFT,RIGHT] contains the desired eigenvalue */

/*           Do while( CNT(LEFT).GT.I-1 ) */

#line 256 "dlarrj.f"
	    fac = 1.;
#line 257 "dlarrj.f"
L20:
#line 258 "dlarrj.f"
	    cnt = 0;
#line 259 "dlarrj.f"
	    s = left;
#line 260 "dlarrj.f"
	    dplus = d__[1] - s;
#line 261 "dlarrj.f"
	    if (dplus < 0.) {
#line 261 "dlarrj.f"
		++cnt;
#line 261 "dlarrj.f"
	    }
#line 262 "dlarrj.f"
	    i__2 = *n;
#line 262 "dlarrj.f"
	    for (j = 2; j <= i__2; ++j) {
#line 263 "dlarrj.f"
		dplus = d__[j] - s - e2[j - 1] / dplus;
#line 264 "dlarrj.f"
		if (dplus < 0.) {
#line 264 "dlarrj.f"
		    ++cnt;
#line 264 "dlarrj.f"
		}
#line 265 "dlarrj.f"
/* L30: */
#line 265 "dlarrj.f"
	    }
#line 266 "dlarrj.f"
	    if (cnt > i__ - 1) {
#line 267 "dlarrj.f"
		left -= werr[ii] * fac;
#line 268 "dlarrj.f"
		fac *= 2.;
#line 269 "dlarrj.f"
		goto L20;
#line 270 "dlarrj.f"
	    }

/*           Do while( CNT(RIGHT).LT.I ) */

#line 274 "dlarrj.f"
	    fac = 1.;
#line 275 "dlarrj.f"
L50:
#line 276 "dlarrj.f"
	    cnt = 0;
#line 277 "dlarrj.f"
	    s = right;
#line 278 "dlarrj.f"
	    dplus = d__[1] - s;
#line 279 "dlarrj.f"
	    if (dplus < 0.) {
#line 279 "dlarrj.f"
		++cnt;
#line 279 "dlarrj.f"
	    }
#line 280 "dlarrj.f"
	    i__2 = *n;
#line 280 "dlarrj.f"
	    for (j = 2; j <= i__2; ++j) {
#line 281 "dlarrj.f"
		dplus = d__[j] - s - e2[j - 1] / dplus;
#line 282 "dlarrj.f"
		if (dplus < 0.) {
#line 282 "dlarrj.f"
		    ++cnt;
#line 282 "dlarrj.f"
		}
#line 283 "dlarrj.f"
/* L60: */
#line 283 "dlarrj.f"
	    }
#line 284 "dlarrj.f"
	    if (cnt < i__) {
#line 285 "dlarrj.f"
		right += werr[ii] * fac;
#line 286 "dlarrj.f"
		fac *= 2.;
#line 287 "dlarrj.f"
		goto L50;
#line 288 "dlarrj.f"
	    }
#line 289 "dlarrj.f"
	    ++nint;
#line 290 "dlarrj.f"
	    iwork[k - 1] = i__ + 1;
#line 291 "dlarrj.f"
	    iwork[k] = cnt;
#line 292 "dlarrj.f"
	}
#line 293 "dlarrj.f"
	work[k - 1] = left;
#line 294 "dlarrj.f"
	work[k] = right;
#line 295 "dlarrj.f"
/* L75: */
#line 295 "dlarrj.f"
    }
#line 298 "dlarrj.f"
    savi1 = i1;

/*     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals */
/*     and while (ITER.LT.MAXITR) */

#line 303 "dlarrj.f"
    iter = 0;
#line 304 "dlarrj.f"
L80:
#line 305 "dlarrj.f"
    prev = i1 - 1;
#line 306 "dlarrj.f"
    i__ = i1;
#line 307 "dlarrj.f"
    olnint = nint;
#line 309 "dlarrj.f"
    i__1 = olnint;
#line 309 "dlarrj.f"
    for (p = 1; p <= i__1; ++p) {
#line 310 "dlarrj.f"
	k = i__ << 1;
#line 311 "dlarrj.f"
	ii = i__ - *offset;
#line 312 "dlarrj.f"
	next = iwork[k - 1];
#line 313 "dlarrj.f"
	left = work[k - 1];
#line 314 "dlarrj.f"
	right = work[k];
#line 315 "dlarrj.f"
	mid = (left + right) * .5;
/*        semiwidth of interval */
#line 318 "dlarrj.f"
	width = right - mid;
/* Computing MAX */
#line 319 "dlarrj.f"
	d__1 = abs(left), d__2 = abs(right);
#line 319 "dlarrj.f"
	tmp = max(d__1,d__2);
#line 321 "dlarrj.f"
	if (width < *rtol * tmp || iter == maxitr) {
/*           reduce number of unconverged intervals */
#line 324 "dlarrj.f"
	    --nint;
/*           Mark interval as converged. */
#line 326 "dlarrj.f"
	    iwork[k - 1] = 0;
#line 327 "dlarrj.f"
	    if (i1 == i__) {
#line 328 "dlarrj.f"
		i1 = next;
#line 329 "dlarrj.f"
	    } else {
/*              Prev holds the last unconverged interval previously examined */
#line 331 "dlarrj.f"
		if (prev >= i1) {
#line 331 "dlarrj.f"
		    iwork[(prev << 1) - 1] = next;
#line 331 "dlarrj.f"
		}
#line 332 "dlarrj.f"
	    }
#line 333 "dlarrj.f"
	    i__ = next;
#line 334 "dlarrj.f"
	    goto L100;
#line 335 "dlarrj.f"
	}
#line 336 "dlarrj.f"
	prev = i__;

/*        Perform one bisection step */

#line 340 "dlarrj.f"
	cnt = 0;
#line 341 "dlarrj.f"
	s = mid;
#line 342 "dlarrj.f"
	dplus = d__[1] - s;
#line 343 "dlarrj.f"
	if (dplus < 0.) {
#line 343 "dlarrj.f"
	    ++cnt;
#line 343 "dlarrj.f"
	}
#line 344 "dlarrj.f"
	i__2 = *n;
#line 344 "dlarrj.f"
	for (j = 2; j <= i__2; ++j) {
#line 345 "dlarrj.f"
	    dplus = d__[j] - s - e2[j - 1] / dplus;
#line 346 "dlarrj.f"
	    if (dplus < 0.) {
#line 346 "dlarrj.f"
		++cnt;
#line 346 "dlarrj.f"
	    }
#line 347 "dlarrj.f"
/* L90: */
#line 347 "dlarrj.f"
	}
#line 348 "dlarrj.f"
	if (cnt <= i__ - 1) {
#line 349 "dlarrj.f"
	    work[k - 1] = mid;
#line 350 "dlarrj.f"
	} else {
#line 351 "dlarrj.f"
	    work[k] = mid;
#line 352 "dlarrj.f"
	}
#line 353 "dlarrj.f"
	i__ = next;
#line 355 "dlarrj.f"
L100:
#line 355 "dlarrj.f"
	;
#line 355 "dlarrj.f"
    }
#line 356 "dlarrj.f"
    ++iter;
/*     do another loop if there are still unconverged intervals */
/*     However, in the last iteration, all intervals are accepted */
/*     since this is the best we can do. */
#line 360 "dlarrj.f"
    if (nint > 0 && iter <= maxitr) {
#line 360 "dlarrj.f"
	goto L80;
#line 360 "dlarrj.f"
    }


/*     At this point, all the intervals have converged */
#line 364 "dlarrj.f"
    i__1 = *ilast;
#line 364 "dlarrj.f"
    for (i__ = savi1; i__ <= i__1; ++i__) {
#line 365 "dlarrj.f"
	k = i__ << 1;
#line 366 "dlarrj.f"
	ii = i__ - *offset;
/*        All intervals marked by '0' have been refined. */
#line 368 "dlarrj.f"
	if (iwork[k - 1] == 0) {
#line 369 "dlarrj.f"
	    w[ii] = (work[k - 1] + work[k]) * .5;
#line 370 "dlarrj.f"
	    werr[ii] = work[k] - w[ii];
#line 371 "dlarrj.f"
	}
#line 372 "dlarrj.f"
/* L110: */
#line 372 "dlarrj.f"
    }

#line 375 "dlarrj.f"
    return 0;

/*     End of DLARRJ */

} /* dlarrj_ */

