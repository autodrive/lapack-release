#line 1 "slarrj.f"
/* slarrj.f -- translated by f2c (version 20100827).
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

#line 1 "slarrj.f"
/* > \brief \b SLARRJ performs refinement of the initial estimates of the eigenvalues of the matrix T. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARRJ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrj.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrj.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrj.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARRJ( N, D, E2, IFIRST, ILAST, */
/*                          RTOL, OFFSET, W, WERR, WORK, IWORK, */
/*                          PIVMIN, SPDIAM, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IFIRST, ILAST, INFO, N, OFFSET */
/*       REAL               PIVMIN, RTOL, SPDIAM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               D( * ), E2( * ), W( * ), */
/*      $                   WERR( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given the initial eigenvalue approximations of T, SLARRJ */
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
/* >          D is REAL array, dimension (N) */
/* >          The N diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] E2 */
/* > \verbatim */
/* >          E2 is REAL array, dimension (N-1) */
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
/* >          RTOL is REAL */
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
/* >          W is REAL array, dimension (N) */
/* >          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are */
/* >          estimates of the eigenvalues of L D L^T indexed IFIRST through */
/* >          ILAST. */
/* >          On output, these estimates are refined. */
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
/* >          The minimum pivot in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[in] SPDIAM */
/* > \verbatim */
/* >          SPDIAM is REAL */
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
/* Subroutine */ int slarrj_(integer *n, doublereal *d__, doublereal *e2, 
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 205 "slarrj.f"
    /* Parameter adjustments */
#line 205 "slarrj.f"
    --iwork;
#line 205 "slarrj.f"
    --work;
#line 205 "slarrj.f"
    --werr;
#line 205 "slarrj.f"
    --w;
#line 205 "slarrj.f"
    --e2;
#line 205 "slarrj.f"
    --d__;
#line 205 "slarrj.f"

#line 205 "slarrj.f"
    /* Function Body */
#line 205 "slarrj.f"
    *info = 0;

#line 207 "slarrj.f"
    maxitr = (integer) ((log(*spdiam + *pivmin) - log(*pivmin)) / log(2.)) + 
	    2;

/*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ]. */
/*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while */
/*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 ) */
/*     for an unconverged interval is set to the index of the next unconverged */
/*     interval, and is -1 or 0 for a converged interval. Thus a linked */
/*     list of unconverged intervals is set up. */

#line 218 "slarrj.f"
    i1 = *ifirst;
#line 219 "slarrj.f"
    i2 = *ilast;
/*     The number of unconverged intervals */
#line 221 "slarrj.f"
    nint = 0;
/*     The last unconverged interval found */
#line 223 "slarrj.f"
    prev = 0;
#line 224 "slarrj.f"
    i__1 = i2;
#line 224 "slarrj.f"
    for (i__ = i1; i__ <= i__1; ++i__) {
#line 225 "slarrj.f"
	k = i__ << 1;
#line 226 "slarrj.f"
	ii = i__ - *offset;
#line 227 "slarrj.f"
	left = w[ii] - werr[ii];
#line 228 "slarrj.f"
	mid = w[ii];
#line 229 "slarrj.f"
	right = w[ii] + werr[ii];
#line 230 "slarrj.f"
	width = right - mid;
/* Computing MAX */
#line 231 "slarrj.f"
	d__1 = abs(left), d__2 = abs(right);
#line 231 "slarrj.f"
	tmp = max(d__1,d__2);
/*        The following test prevents the test of converged intervals */
#line 234 "slarrj.f"
	if (width < *rtol * tmp) {
/*           This interval has already converged and does not need refinement. */
/*           (Note that the gaps might change through refining the */
/*            eigenvalues, however, they can only get bigger.) */
/*           Remove it from the list. */
#line 239 "slarrj.f"
	    iwork[k - 1] = -1;
/*           Make sure that I1 always points to the first unconverged interval */
#line 241 "slarrj.f"
	    if (i__ == i1 && i__ < i2) {
#line 241 "slarrj.f"
		i1 = i__ + 1;
#line 241 "slarrj.f"
	    }
#line 242 "slarrj.f"
	    if (prev >= i1 && i__ <= i2) {
#line 242 "slarrj.f"
		iwork[(prev << 1) - 1] = i__ + 1;
#line 242 "slarrj.f"
	    }
#line 243 "slarrj.f"
	} else {
/*           unconverged interval found */
#line 245 "slarrj.f"
	    prev = i__;
/*           Make sure that [LEFT,RIGHT] contains the desired eigenvalue */

/*           Do while( CNT(LEFT).GT.I-1 ) */

#line 250 "slarrj.f"
	    fac = 1.;
#line 251 "slarrj.f"
L20:
#line 252 "slarrj.f"
	    cnt = 0;
#line 253 "slarrj.f"
	    s = left;
#line 254 "slarrj.f"
	    dplus = d__[1] - s;
#line 255 "slarrj.f"
	    if (dplus < 0.) {
#line 255 "slarrj.f"
		++cnt;
#line 255 "slarrj.f"
	    }
#line 256 "slarrj.f"
	    i__2 = *n;
#line 256 "slarrj.f"
	    for (j = 2; j <= i__2; ++j) {
#line 257 "slarrj.f"
		dplus = d__[j] - s - e2[j - 1] / dplus;
#line 258 "slarrj.f"
		if (dplus < 0.) {
#line 258 "slarrj.f"
		    ++cnt;
#line 258 "slarrj.f"
		}
#line 259 "slarrj.f"
/* L30: */
#line 259 "slarrj.f"
	    }
#line 260 "slarrj.f"
	    if (cnt > i__ - 1) {
#line 261 "slarrj.f"
		left -= werr[ii] * fac;
#line 262 "slarrj.f"
		fac *= 2.;
#line 263 "slarrj.f"
		goto L20;
#line 264 "slarrj.f"
	    }

/*           Do while( CNT(RIGHT).LT.I ) */

#line 268 "slarrj.f"
	    fac = 1.;
#line 269 "slarrj.f"
L50:
#line 270 "slarrj.f"
	    cnt = 0;
#line 271 "slarrj.f"
	    s = right;
#line 272 "slarrj.f"
	    dplus = d__[1] - s;
#line 273 "slarrj.f"
	    if (dplus < 0.) {
#line 273 "slarrj.f"
		++cnt;
#line 273 "slarrj.f"
	    }
#line 274 "slarrj.f"
	    i__2 = *n;
#line 274 "slarrj.f"
	    for (j = 2; j <= i__2; ++j) {
#line 275 "slarrj.f"
		dplus = d__[j] - s - e2[j - 1] / dplus;
#line 276 "slarrj.f"
		if (dplus < 0.) {
#line 276 "slarrj.f"
		    ++cnt;
#line 276 "slarrj.f"
		}
#line 277 "slarrj.f"
/* L60: */
#line 277 "slarrj.f"
	    }
#line 278 "slarrj.f"
	    if (cnt < i__) {
#line 279 "slarrj.f"
		right += werr[ii] * fac;
#line 280 "slarrj.f"
		fac *= 2.;
#line 281 "slarrj.f"
		goto L50;
#line 282 "slarrj.f"
	    }
#line 283 "slarrj.f"
	    ++nint;
#line 284 "slarrj.f"
	    iwork[k - 1] = i__ + 1;
#line 285 "slarrj.f"
	    iwork[k] = cnt;
#line 286 "slarrj.f"
	}
#line 287 "slarrj.f"
	work[k - 1] = left;
#line 288 "slarrj.f"
	work[k] = right;
#line 289 "slarrj.f"
/* L75: */
#line 289 "slarrj.f"
    }
#line 292 "slarrj.f"
    savi1 = i1;

/*     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals */
/*     and while (ITER.LT.MAXITR) */

#line 297 "slarrj.f"
    iter = 0;
#line 298 "slarrj.f"
L80:
#line 299 "slarrj.f"
    prev = i1 - 1;
#line 300 "slarrj.f"
    i__ = i1;
#line 301 "slarrj.f"
    olnint = nint;
#line 303 "slarrj.f"
    i__1 = olnint;
#line 303 "slarrj.f"
    for (p = 1; p <= i__1; ++p) {
#line 304 "slarrj.f"
	k = i__ << 1;
#line 305 "slarrj.f"
	ii = i__ - *offset;
#line 306 "slarrj.f"
	next = iwork[k - 1];
#line 307 "slarrj.f"
	left = work[k - 1];
#line 308 "slarrj.f"
	right = work[k];
#line 309 "slarrj.f"
	mid = (left + right) * .5;
/*        semiwidth of interval */
#line 312 "slarrj.f"
	width = right - mid;
/* Computing MAX */
#line 313 "slarrj.f"
	d__1 = abs(left), d__2 = abs(right);
#line 313 "slarrj.f"
	tmp = max(d__1,d__2);
#line 315 "slarrj.f"
	if (width < *rtol * tmp || iter == maxitr) {
/*           reduce number of unconverged intervals */
#line 318 "slarrj.f"
	    --nint;
/*           Mark interval as converged. */
#line 320 "slarrj.f"
	    iwork[k - 1] = 0;
#line 321 "slarrj.f"
	    if (i1 == i__) {
#line 322 "slarrj.f"
		i1 = next;
#line 323 "slarrj.f"
	    } else {
/*              Prev holds the last unconverged interval previously examined */
#line 325 "slarrj.f"
		if (prev >= i1) {
#line 325 "slarrj.f"
		    iwork[(prev << 1) - 1] = next;
#line 325 "slarrj.f"
		}
#line 326 "slarrj.f"
	    }
#line 327 "slarrj.f"
	    i__ = next;
#line 328 "slarrj.f"
	    goto L100;
#line 329 "slarrj.f"
	}
#line 330 "slarrj.f"
	prev = i__;

/*        Perform one bisection step */

#line 334 "slarrj.f"
	cnt = 0;
#line 335 "slarrj.f"
	s = mid;
#line 336 "slarrj.f"
	dplus = d__[1] - s;
#line 337 "slarrj.f"
	if (dplus < 0.) {
#line 337 "slarrj.f"
	    ++cnt;
#line 337 "slarrj.f"
	}
#line 338 "slarrj.f"
	i__2 = *n;
#line 338 "slarrj.f"
	for (j = 2; j <= i__2; ++j) {
#line 339 "slarrj.f"
	    dplus = d__[j] - s - e2[j - 1] / dplus;
#line 340 "slarrj.f"
	    if (dplus < 0.) {
#line 340 "slarrj.f"
		++cnt;
#line 340 "slarrj.f"
	    }
#line 341 "slarrj.f"
/* L90: */
#line 341 "slarrj.f"
	}
#line 342 "slarrj.f"
	if (cnt <= i__ - 1) {
#line 343 "slarrj.f"
	    work[k - 1] = mid;
#line 344 "slarrj.f"
	} else {
#line 345 "slarrj.f"
	    work[k] = mid;
#line 346 "slarrj.f"
	}
#line 347 "slarrj.f"
	i__ = next;
#line 349 "slarrj.f"
L100:
#line 349 "slarrj.f"
	;
#line 349 "slarrj.f"
    }
#line 350 "slarrj.f"
    ++iter;
/*     do another loop if there are still unconverged intervals */
/*     However, in the last iteration, all intervals are accepted */
/*     since this is the best we can do. */
#line 354 "slarrj.f"
    if (nint > 0 && iter <= maxitr) {
#line 354 "slarrj.f"
	goto L80;
#line 354 "slarrj.f"
    }


/*     At this point, all the intervals have converged */
#line 358 "slarrj.f"
    i__1 = *ilast;
#line 358 "slarrj.f"
    for (i__ = savi1; i__ <= i__1; ++i__) {
#line 359 "slarrj.f"
	k = i__ << 1;
#line 360 "slarrj.f"
	ii = i__ - *offset;
/*        All intervals marked by '0' have been refined. */
#line 362 "slarrj.f"
	if (iwork[k - 1] == 0) {
#line 363 "slarrj.f"
	    w[ii] = (work[k - 1] + work[k]) * .5;
#line 364 "slarrj.f"
	    werr[ii] = work[k] - w[ii];
#line 365 "slarrj.f"
	}
#line 366 "slarrj.f"
/* L110: */
#line 366 "slarrj.f"
    }

#line 369 "slarrj.f"
    return 0;

/*     End of SLARRJ */

} /* slarrj_ */

