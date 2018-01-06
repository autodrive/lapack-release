#line 1 "sstebz.f"
/* sstebz.f -- translated by f2c (version 20100827).
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

#line 1 "sstebz.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__0 = 0;

/* > \brief \b SSTEBZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEBZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstebz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstebz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstebz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, */
/*                          M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          ORDER, RANGE */
/*       INTEGER            IL, INFO, IU, M, N, NSPLIT */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ) */
/*       REAL               D( * ), E( * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEBZ computes the eigenvalues of a symmetric tridiagonal */
/* > matrix T.  The user may ask for all eigenvalues, all eigenvalues */
/* > in the half-open interval (VL, VU], or the IL-th through IU-th */
/* > eigenvalues. */
/* > */
/* > To avoid overflow, the matrix must be scaled so that its */
/* > largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
 */
/* > accuracy, it should not be much smaller than that. */
/* > */
/* > See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal */
/* > Matrix", Report CS41, Computer Science Dept., Stanford */
/* > University, July 21, 1966. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] RANGE */
/* > \verbatim */
/* >          RANGE is CHARACTER*1 */
/* >          = 'A': ("All")   all eigenvalues will be found. */
/* >          = 'V': ("Value") all eigenvalues in the half-open interval */
/* >                           (VL, VU] will be found. */
/* >          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the */
/* >                           entire matrix) will be found. */
/* > \endverbatim */
/* > */
/* > \param[in] ORDER */
/* > \verbatim */
/* >          ORDER is CHARACTER*1 */
/* >          = 'B': ("By Block") the eigenvalues will be grouped by */
/* >                              split-off block (see IBLOCK, ISPLIT) and */
/* >                              ordered from smallest to largest within */
/* >                              the block. */
/* >          = 'E': ("Entire matrix") */
/* >                              the eigenvalues for the entire matrix */
/* >                              will be ordered from smallest to */
/* >                              largest. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the tridiagonal matrix T.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is REAL */
/* > */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues.  Eigenvalues less than or equal */
/* >          to VL, or greater than VU, will not be returned.  VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* > */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for eigenvalues.  Eigenvalues less than or equal */
/* >          to VL, or greater than VU, will not be returned.  VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* > */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* > */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* >          ABSTOL is REAL */
/* >          The absolute tolerance for the eigenvalues.  An eigenvalue */
/* >          (or cluster) is considered to be located if it has been */
/* >          determined to lie in an interval whose width is ABSTOL or */
/* >          less.  If ABSTOL is less than or equal to zero, then ULP*|T| */
/* >          will be used, where |T| means the 1-norm of T. */
/* > */
/* >          Eigenvalues will be computed most accurately when ABSTOL is */
/* >          set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          The (n-1) off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The actual number of eigenvalues found. 0 <= M <= N. */
/* >          (See also the description of INFO=2,3.) */
/* > \endverbatim */
/* > */
/* > \param[out] NSPLIT */
/* > \verbatim */
/* >          NSPLIT is INTEGER */
/* >          The number of diagonal blocks in the matrix T. */
/* >          1 <= NSPLIT <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          On exit, the first M elements of W will contain the */
/* >          eigenvalues.  (SSTEBZ may use the remaining N-M elements as */
/* >          workspace.) */
/* > \endverbatim */
/* > */
/* > \param[out] IBLOCK */
/* > \verbatim */
/* >          IBLOCK is INTEGER array, dimension (N) */
/* >          At each row/column j where E(j) is zero or small, the */
/* >          matrix T is considered to split into a block diagonal */
/* >          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which */
/* >          block (from 1 to the number of blocks) the eigenvalue W(i) */
/* >          belongs.  (SSTEBZ may use the remaining N-M elements as */
/* >          workspace.) */
/* > \endverbatim */
/* > */
/* > \param[out] ISPLIT */
/* > \verbatim */
/* >          ISPLIT is INTEGER array, dimension (N) */
/* >          The splitting points, at which T breaks up into submatrices. */
/* >          The first submatrix consists of rows/columns 1 to ISPLIT(1), */
/* >          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2), */
/* >          etc., and the NSPLIT-th consists of rows/columns */
/* >          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. */
/* >          (Only the first NSPLIT elements will actually be used, but */
/* >          since the user cannot know a priori what value NSPLIT will */
/* >          have, N words must be reserved for ISPLIT.) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  some or all of the eigenvalues failed to converge or */
/* >                were not computed: */
/* >                =1 or 3: Bisection failed to converge for some */
/* >                        eigenvalues; these eigenvalues are flagged by a */
/* >                        negative block number.  The effect is that the */
/* >                        eigenvalues may not be as accurate as the */
/* >                        absolute and relative tolerances.  This is */
/* >                        generally caused by unexpectedly inaccurate */
/* >                        arithmetic. */
/* >                =2 or 3: RANGE='I' only: Not all of the eigenvalues */
/* >                        IL:IU were found. */
/* >                        Effect: M < IU+1-IL */
/* >                        Cause:  non-monotonic arithmetic, causing the */
/* >                                Sturm sequence to be non-monotonic. */
/* >                        Cure:   recalculate, using RANGE='A', and pick */
/* >                                out eigenvalues IL:IU.  In some cases, */
/* >                                increasing the PARAMETER "FUDGE" may */
/* >                                make things work. */
/* >                = 4:    RANGE='I', and the Gershgorin interval */
/* >                        initially used was too small.  No eigenvalues */
/* >                        were computed. */
/* >                        Probable cause: your machine has sloppy */
/* >                                        floating-point arithmetic. */
/* >                        Cure: Increase the PARAMETER "FUDGE", */
/* >                              recompile, and try again. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  RELFAC  REAL, default = 2.0e0 */
/* >          The relative tolerance.  An interval (a,b] lies within */
/* >          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|), */
/* >          where "ulp" is the machine precision (distance from 1 to */
/* >          the next larger floating point number.) */
/* > */
/* >  FUDGE   REAL, default = 2 */
/* >          A "fudge factor" to widen the Gershgorin intervals.  Ideally, */
/* >          a value of 1 should work, but on machines with sloppy */
/* >          arithmetic, this needs to be larger.  The default for */
/* >          publicly released versions should be large enough to handle */
/* >          the worst machine around.  Note that this has no effect */
/* >          on accuracy of the solution. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sstebz_(char *range, char *order, integer *n, doublereal 
	*vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	doublereal *d__, doublereal *e, integer *m, integer *nsplit, 
	doublereal *w, integer *iblock, integer *isplit, doublereal *work, 
	integer *iwork, integer *info, ftnlen range_len, ftnlen order_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static integer j, ib, jb, ie, je, nb;
    static doublereal gl;
    static integer im, in;
    static doublereal gu;
    static integer iw;
    static doublereal wl, wu;
    static integer nwl;
    static doublereal ulp, wlu, wul;
    static integer nwu;
    static doublereal tmp1, tmp2;
    static integer iend, ioff, iout, itmp1, jdisc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static doublereal atoli;
    static integer iwoff;
    static doublereal bnorm;
    static integer itmax;
    static doublereal wkill, rtoli, tnorm;
    static integer ibegin, irange, idiscl;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safemn;
    static integer idumma[1];
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer idiscu;
    extern /* Subroutine */ int slaebz_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    static integer iorder;
    static logical ncnvrg;
    static doublereal pivmin;
    static logical toofew;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 326 "sstebz.f"
    /* Parameter adjustments */
#line 326 "sstebz.f"
    --iwork;
#line 326 "sstebz.f"
    --work;
#line 326 "sstebz.f"
    --isplit;
#line 326 "sstebz.f"
    --iblock;
#line 326 "sstebz.f"
    --w;
#line 326 "sstebz.f"
    --e;
#line 326 "sstebz.f"
    --d__;
#line 326 "sstebz.f"

#line 326 "sstebz.f"
    /* Function Body */
#line 326 "sstebz.f"
    *info = 0;

/*     Decode RANGE */

#line 330 "sstebz.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 331 "sstebz.f"
	irange = 1;
#line 332 "sstebz.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 333 "sstebz.f"
	irange = 2;
#line 334 "sstebz.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 335 "sstebz.f"
	irange = 3;
#line 336 "sstebz.f"
    } else {
#line 337 "sstebz.f"
	irange = 0;
#line 338 "sstebz.f"
    }

/*     Decode ORDER */

#line 342 "sstebz.f"
    if (lsame_(order, "B", (ftnlen)1, (ftnlen)1)) {
#line 343 "sstebz.f"
	iorder = 2;
#line 344 "sstebz.f"
    } else if (lsame_(order, "E", (ftnlen)1, (ftnlen)1)) {
#line 345 "sstebz.f"
	iorder = 1;
#line 346 "sstebz.f"
    } else {
#line 347 "sstebz.f"
	iorder = 0;
#line 348 "sstebz.f"
    }

/*     Check for Errors */

#line 352 "sstebz.f"
    if (irange <= 0) {
#line 353 "sstebz.f"
	*info = -1;
#line 354 "sstebz.f"
    } else if (iorder <= 0) {
#line 355 "sstebz.f"
	*info = -2;
#line 356 "sstebz.f"
    } else if (*n < 0) {
#line 357 "sstebz.f"
	*info = -3;
#line 358 "sstebz.f"
    } else if (irange == 2) {
#line 359 "sstebz.f"
	if (*vl >= *vu) {
#line 359 "sstebz.f"
	    *info = -5;
#line 359 "sstebz.f"
	}
#line 360 "sstebz.f"
    } else if (irange == 3 && (*il < 1 || *il > max(1,*n))) {
#line 362 "sstebz.f"
	*info = -6;
#line 363 "sstebz.f"
    } else if (irange == 3 && (*iu < min(*n,*il) || *iu > *n)) {
#line 365 "sstebz.f"
	*info = -7;
#line 366 "sstebz.f"
    }

#line 368 "sstebz.f"
    if (*info != 0) {
#line 369 "sstebz.f"
	i__1 = -(*info);
#line 369 "sstebz.f"
	xerbla_("SSTEBZ", &i__1, (ftnlen)6);
#line 370 "sstebz.f"
	return 0;
#line 371 "sstebz.f"
    }

/*     Initialize error flags */

#line 375 "sstebz.f"
    *info = 0;
#line 376 "sstebz.f"
    ncnvrg = FALSE_;
#line 377 "sstebz.f"
    toofew = FALSE_;

/*     Quick return if possible */

#line 381 "sstebz.f"
    *m = 0;
#line 382 "sstebz.f"
    if (*n == 0) {
#line 382 "sstebz.f"
	return 0;
#line 382 "sstebz.f"
    }

/*     Simplifications: */

#line 387 "sstebz.f"
    if (irange == 3 && *il == 1 && *iu == *n) {
#line 387 "sstebz.f"
	irange = 1;
#line 387 "sstebz.f"
    }

/*     Get machine constants */
/*     NB is the minimum vector length for vector bisection, or 0 */
/*     if only scalar is to be done. */

#line 394 "sstebz.f"
    safemn = slamch_("S", (ftnlen)1);
#line 395 "sstebz.f"
    ulp = slamch_("P", (ftnlen)1);
#line 396 "sstebz.f"
    rtoli = ulp * 2.;
#line 397 "sstebz.f"
    nb = ilaenv_(&c__1, "SSTEBZ", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 398 "sstebz.f"
    if (nb <= 1) {
#line 398 "sstebz.f"
	nb = 0;
#line 398 "sstebz.f"
    }

/*     Special Case when N=1 */

#line 403 "sstebz.f"
    if (*n == 1) {
#line 404 "sstebz.f"
	*nsplit = 1;
#line 405 "sstebz.f"
	isplit[1] = 1;
#line 406 "sstebz.f"
	if (irange == 2 && (*vl >= d__[1] || *vu < d__[1])) {
#line 407 "sstebz.f"
	    *m = 0;
#line 408 "sstebz.f"
	} else {
#line 409 "sstebz.f"
	    w[1] = d__[1];
#line 410 "sstebz.f"
	    iblock[1] = 1;
#line 411 "sstebz.f"
	    *m = 1;
#line 412 "sstebz.f"
	}
#line 413 "sstebz.f"
	return 0;
#line 414 "sstebz.f"
    }

/*     Compute Splitting Points */

#line 418 "sstebz.f"
    *nsplit = 1;
#line 419 "sstebz.f"
    work[*n] = 0.;
#line 420 "sstebz.f"
    pivmin = 1.;

#line 422 "sstebz.f"
    i__1 = *n;
#line 422 "sstebz.f"
    for (j = 2; j <= i__1; ++j) {
/* Computing 2nd power */
#line 423 "sstebz.f"
	d__1 = e[j - 1];
#line 423 "sstebz.f"
	tmp1 = d__1 * d__1;
/* Computing 2nd power */
#line 424 "sstebz.f"
	d__2 = ulp;
#line 424 "sstebz.f"
	if ((d__1 = d__[j] * d__[j - 1], abs(d__1)) * (d__2 * d__2) + safemn 
		> tmp1) {
#line 425 "sstebz.f"
	    isplit[*nsplit] = j - 1;
#line 426 "sstebz.f"
	    ++(*nsplit);
#line 427 "sstebz.f"
	    work[j - 1] = 0.;
#line 428 "sstebz.f"
	} else {
#line 429 "sstebz.f"
	    work[j - 1] = tmp1;
#line 430 "sstebz.f"
	    pivmin = max(pivmin,tmp1);
#line 431 "sstebz.f"
	}
#line 432 "sstebz.f"
/* L10: */
#line 432 "sstebz.f"
    }
#line 433 "sstebz.f"
    isplit[*nsplit] = *n;
#line 434 "sstebz.f"
    pivmin *= safemn;

/*     Compute Interval and ATOLI */

#line 438 "sstebz.f"
    if (irange == 3) {

/*        RANGE='I': Compute the interval containing eigenvalues */
/*                   IL through IU. */

/*        Compute Gershgorin interval for entire (split) matrix */
/*        and use it as the initial interval */

#line 446 "sstebz.f"
	gu = d__[1];
#line 447 "sstebz.f"
	gl = d__[1];
#line 448 "sstebz.f"
	tmp1 = 0.;

#line 450 "sstebz.f"
	i__1 = *n - 1;
#line 450 "sstebz.f"
	for (j = 1; j <= i__1; ++j) {
#line 451 "sstebz.f"
	    tmp2 = sqrt(work[j]);
/* Computing MAX */
#line 452 "sstebz.f"
	    d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
#line 452 "sstebz.f"
	    gu = max(d__1,d__2);
/* Computing MIN */
#line 453 "sstebz.f"
	    d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
#line 453 "sstebz.f"
	    gl = min(d__1,d__2);
#line 454 "sstebz.f"
	    tmp1 = tmp2;
#line 455 "sstebz.f"
/* L20: */
#line 455 "sstebz.f"
	}

/* Computing MAX */
#line 457 "sstebz.f"
	d__1 = gu, d__2 = d__[*n] + tmp1;
#line 457 "sstebz.f"
	gu = max(d__1,d__2);
/* Computing MIN */
#line 458 "sstebz.f"
	d__1 = gl, d__2 = d__[*n] - tmp1;
#line 458 "sstebz.f"
	gl = min(d__1,d__2);
/* Computing MAX */
#line 459 "sstebz.f"
	d__1 = abs(gl), d__2 = abs(gu);
#line 459 "sstebz.f"
	tnorm = max(d__1,d__2);
#line 460 "sstebz.f"
	gl = gl - tnorm * 2.1 * ulp * *n - pivmin * 4.2000000000000002;
#line 461 "sstebz.f"
	gu = gu + tnorm * 2.1 * ulp * *n + pivmin * 2.1;

/*        Compute Iteration parameters */

#line 465 "sstebz.f"
	itmax = (integer) ((log(tnorm + pivmin) - log(pivmin)) / log(2.)) + 2;
#line 467 "sstebz.f"
	if (*abstol <= 0.) {
#line 468 "sstebz.f"
	    atoli = ulp * tnorm;
#line 469 "sstebz.f"
	} else {
#line 470 "sstebz.f"
	    atoli = *abstol;
#line 471 "sstebz.f"
	}

#line 473 "sstebz.f"
	work[*n + 1] = gl;
#line 474 "sstebz.f"
	work[*n + 2] = gl;
#line 475 "sstebz.f"
	work[*n + 3] = gu;
#line 476 "sstebz.f"
	work[*n + 4] = gu;
#line 477 "sstebz.f"
	work[*n + 5] = gl;
#line 478 "sstebz.f"
	work[*n + 6] = gu;
#line 479 "sstebz.f"
	iwork[1] = -1;
#line 480 "sstebz.f"
	iwork[2] = -1;
#line 481 "sstebz.f"
	iwork[3] = *n + 1;
#line 482 "sstebz.f"
	iwork[4] = *n + 1;
#line 483 "sstebz.f"
	iwork[5] = *il - 1;
#line 484 "sstebz.f"
	iwork[6] = *iu;

#line 486 "sstebz.f"
	slaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, &pivmin, 
		&d__[1], &e[1], &work[1], &iwork[5], &work[*n + 1], &work[*n 
		+ 5], &iout, &iwork[1], &w[1], &iblock[1], &iinfo);

#line 490 "sstebz.f"
	if (iwork[6] == *iu) {
#line 491 "sstebz.f"
	    wl = work[*n + 1];
#line 492 "sstebz.f"
	    wlu = work[*n + 3];
#line 493 "sstebz.f"
	    nwl = iwork[1];
#line 494 "sstebz.f"
	    wu = work[*n + 4];
#line 495 "sstebz.f"
	    wul = work[*n + 2];
#line 496 "sstebz.f"
	    nwu = iwork[4];
#line 497 "sstebz.f"
	} else {
#line 498 "sstebz.f"
	    wl = work[*n + 2];
#line 499 "sstebz.f"
	    wlu = work[*n + 4];
#line 500 "sstebz.f"
	    nwl = iwork[2];
#line 501 "sstebz.f"
	    wu = work[*n + 3];
#line 502 "sstebz.f"
	    wul = work[*n + 1];
#line 503 "sstebz.f"
	    nwu = iwork[3];
#line 504 "sstebz.f"
	}

#line 506 "sstebz.f"
	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
#line 507 "sstebz.f"
	    *info = 4;
#line 508 "sstebz.f"
	    return 0;
#line 509 "sstebz.f"
	}
#line 510 "sstebz.f"
    } else {

/*        RANGE='A' or 'V' -- Set ATOLI */

/* Computing MAX */
#line 514 "sstebz.f"
	d__3 = abs(d__[1]) + abs(e[1]), d__4 = (d__1 = d__[*n], abs(d__1)) + (
		d__2 = e[*n - 1], abs(d__2));
#line 514 "sstebz.f"
	tnorm = max(d__3,d__4);

#line 517 "sstebz.f"
	i__1 = *n - 1;
#line 517 "sstebz.f"
	for (j = 2; j <= i__1; ++j) {
/* Computing MAX */
#line 518 "sstebz.f"
	    d__4 = tnorm, d__5 = (d__1 = d__[j], abs(d__1)) + (d__2 = e[j - 1]
		    , abs(d__2)) + (d__3 = e[j], abs(d__3));
#line 518 "sstebz.f"
	    tnorm = max(d__4,d__5);
#line 520 "sstebz.f"
/* L30: */
#line 520 "sstebz.f"
	}

#line 522 "sstebz.f"
	if (*abstol <= 0.) {
#line 523 "sstebz.f"
	    atoli = ulp * tnorm;
#line 524 "sstebz.f"
	} else {
#line 525 "sstebz.f"
	    atoli = *abstol;
#line 526 "sstebz.f"
	}

#line 528 "sstebz.f"
	if (irange == 2) {
#line 529 "sstebz.f"
	    wl = *vl;
#line 530 "sstebz.f"
	    wu = *vu;
#line 531 "sstebz.f"
	} else {
#line 532 "sstebz.f"
	    wl = 0.;
#line 533 "sstebz.f"
	    wu = 0.;
#line 534 "sstebz.f"
	}
#line 535 "sstebz.f"
    }

/*     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU. */
/*     NWL accumulates the number of eigenvalues .le. WL, */
/*     NWU accumulates the number of eigenvalues .le. WU */

#line 541 "sstebz.f"
    *m = 0;
#line 542 "sstebz.f"
    iend = 0;
#line 543 "sstebz.f"
    *info = 0;
#line 544 "sstebz.f"
    nwl = 0;
#line 545 "sstebz.f"
    nwu = 0;

#line 547 "sstebz.f"
    i__1 = *nsplit;
#line 547 "sstebz.f"
    for (jb = 1; jb <= i__1; ++jb) {
#line 548 "sstebz.f"
	ioff = iend;
#line 549 "sstebz.f"
	ibegin = ioff + 1;
#line 550 "sstebz.f"
	iend = isplit[jb];
#line 551 "sstebz.f"
	in = iend - ioff;

#line 553 "sstebz.f"
	if (in == 1) {

/*           Special Case -- IN=1 */

#line 557 "sstebz.f"
	    if (irange == 1 || wl >= d__[ibegin] - pivmin) {
#line 557 "sstebz.f"
		++nwl;
#line 557 "sstebz.f"
	    }
#line 559 "sstebz.f"
	    if (irange == 1 || wu >= d__[ibegin] - pivmin) {
#line 559 "sstebz.f"
		++nwu;
#line 559 "sstebz.f"
	    }
#line 561 "sstebz.f"
	    if (irange == 1 || wl < d__[ibegin] - pivmin && wu >= d__[ibegin] 
		    - pivmin) {
#line 563 "sstebz.f"
		++(*m);
#line 564 "sstebz.f"
		w[*m] = d__[ibegin];
#line 565 "sstebz.f"
		iblock[*m] = jb;
#line 566 "sstebz.f"
	    }
#line 567 "sstebz.f"
	} else {

/*           General Case -- IN > 1 */

/*           Compute Gershgorin Interval */
/*           and use it as the initial interval */

#line 574 "sstebz.f"
	    gu = d__[ibegin];
#line 575 "sstebz.f"
	    gl = d__[ibegin];
#line 576 "sstebz.f"
	    tmp1 = 0.;

#line 578 "sstebz.f"
	    i__2 = iend - 1;
#line 578 "sstebz.f"
	    for (j = ibegin; j <= i__2; ++j) {
#line 579 "sstebz.f"
		tmp2 = (d__1 = e[j], abs(d__1));
/* Computing MAX */
#line 580 "sstebz.f"
		d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
#line 580 "sstebz.f"
		gu = max(d__1,d__2);
/* Computing MIN */
#line 581 "sstebz.f"
		d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
#line 581 "sstebz.f"
		gl = min(d__1,d__2);
#line 582 "sstebz.f"
		tmp1 = tmp2;
#line 583 "sstebz.f"
/* L40: */
#line 583 "sstebz.f"
	    }

/* Computing MAX */
#line 585 "sstebz.f"
	    d__1 = gu, d__2 = d__[iend] + tmp1;
#line 585 "sstebz.f"
	    gu = max(d__1,d__2);
/* Computing MIN */
#line 586 "sstebz.f"
	    d__1 = gl, d__2 = d__[iend] - tmp1;
#line 586 "sstebz.f"
	    gl = min(d__1,d__2);
/* Computing MAX */
#line 587 "sstebz.f"
	    d__1 = abs(gl), d__2 = abs(gu);
#line 587 "sstebz.f"
	    bnorm = max(d__1,d__2);
#line 588 "sstebz.f"
	    gl = gl - bnorm * 2.1 * ulp * in - pivmin * 2.1;
#line 589 "sstebz.f"
	    gu = gu + bnorm * 2.1 * ulp * in + pivmin * 2.1;

/*           Compute ATOLI for the current submatrix */

#line 593 "sstebz.f"
	    if (*abstol <= 0.) {
/* Computing MAX */
#line 594 "sstebz.f"
		d__1 = abs(gl), d__2 = abs(gu);
#line 594 "sstebz.f"
		atoli = ulp * max(d__1,d__2);
#line 595 "sstebz.f"
	    } else {
#line 596 "sstebz.f"
		atoli = *abstol;
#line 597 "sstebz.f"
	    }

#line 599 "sstebz.f"
	    if (irange > 1) {
#line 600 "sstebz.f"
		if (gu < wl) {
#line 601 "sstebz.f"
		    nwl += in;
#line 602 "sstebz.f"
		    nwu += in;
#line 603 "sstebz.f"
		    goto L70;
#line 604 "sstebz.f"
		}
#line 605 "sstebz.f"
		gl = max(gl,wl);
#line 606 "sstebz.f"
		gu = min(gu,wu);
#line 607 "sstebz.f"
		if (gl >= gu) {
#line 607 "sstebz.f"
		    goto L70;
#line 607 "sstebz.f"
		}
#line 609 "sstebz.f"
	    }

/*           Set Up Initial Interval */

#line 613 "sstebz.f"
	    work[*n + 1] = gl;
#line 614 "sstebz.f"
	    work[*n + in + 1] = gu;
#line 615 "sstebz.f"
	    slaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);

#line 620 "sstebz.f"
	    nwl += iwork[1];
#line 621 "sstebz.f"
	    nwu += iwork[in + 1];
#line 622 "sstebz.f"
	    iwoff = *m - iwork[1];

/*           Compute Eigenvalues */

#line 626 "sstebz.f"
	    itmax = (integer) ((log(gu - gl + pivmin) - log(pivmin)) / log(2.)
		    ) + 2;
#line 628 "sstebz.f"
	    slaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);

/*           Copy Eigenvalues Into W and IBLOCK */
/*           Use -JB for block number for unconverged eigenvalues. */

#line 636 "sstebz.f"
	    i__2 = iout;
#line 636 "sstebz.f"
	    for (j = 1; j <= i__2; ++j) {
#line 637 "sstebz.f"
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;

/*              Flag non-convergence. */

#line 641 "sstebz.f"
		if (j > iout - iinfo) {
#line 642 "sstebz.f"
		    ncnvrg = TRUE_;
#line 643 "sstebz.f"
		    ib = -jb;
#line 644 "sstebz.f"
		} else {
#line 645 "sstebz.f"
		    ib = jb;
#line 646 "sstebz.f"
		}
#line 647 "sstebz.f"
		i__3 = iwork[j + in] + iwoff;
#line 647 "sstebz.f"
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
#line 649 "sstebz.f"
		    w[je] = tmp1;
#line 650 "sstebz.f"
		    iblock[je] = ib;
#line 651 "sstebz.f"
/* L50: */
#line 651 "sstebz.f"
		}
#line 652 "sstebz.f"
/* L60: */
#line 652 "sstebz.f"
	    }

#line 654 "sstebz.f"
	    *m += im;
#line 655 "sstebz.f"
	}
#line 656 "sstebz.f"
L70:
#line 656 "sstebz.f"
	;
#line 656 "sstebz.f"
    }

/*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU */
/*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues. */

#line 661 "sstebz.f"
    if (irange == 3) {
#line 662 "sstebz.f"
	im = 0;
#line 663 "sstebz.f"
	idiscl = *il - 1 - nwl;
#line 664 "sstebz.f"
	idiscu = nwu - *iu;

#line 666 "sstebz.f"
	if (idiscl > 0 || idiscu > 0) {
#line 667 "sstebz.f"
	    i__1 = *m;
#line 667 "sstebz.f"
	    for (je = 1; je <= i__1; ++je) {
#line 668 "sstebz.f"
		if (w[je] <= wlu && idiscl > 0) {
#line 669 "sstebz.f"
		    --idiscl;
#line 670 "sstebz.f"
		} else if (w[je] >= wul && idiscu > 0) {
#line 671 "sstebz.f"
		    --idiscu;
#line 672 "sstebz.f"
		} else {
#line 673 "sstebz.f"
		    ++im;
#line 674 "sstebz.f"
		    w[im] = w[je];
#line 675 "sstebz.f"
		    iblock[im] = iblock[je];
#line 676 "sstebz.f"
		}
#line 677 "sstebz.f"
/* L80: */
#line 677 "sstebz.f"
	    }
#line 678 "sstebz.f"
	    *m = im;
#line 679 "sstebz.f"
	}
#line 680 "sstebz.f"
	if (idiscl > 0 || idiscu > 0) {

/*           Code to deal with effects of bad arithmetic: */
/*           Some low eigenvalues to be discarded are not in (WL,WLU], */
/*           or high eigenvalues to be discarded are not in (WUL,WU] */
/*           so just kill off the smallest IDISCL/largest IDISCU */
/*           eigenvalues, by simply finding the smallest/largest */
/*           eigenvalue(s). */

/*           (If N(w) is monotone non-decreasing, this should never */
/*               happen.) */

#line 692 "sstebz.f"
	    if (idiscl > 0) {
#line 693 "sstebz.f"
		wkill = wu;
#line 694 "sstebz.f"
		i__1 = idiscl;
#line 694 "sstebz.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 695 "sstebz.f"
		    iw = 0;
#line 696 "sstebz.f"
		    i__2 = *m;
#line 696 "sstebz.f"
		    for (je = 1; je <= i__2; ++je) {
#line 697 "sstebz.f"
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
#line 699 "sstebz.f"
			    iw = je;
#line 700 "sstebz.f"
			    wkill = w[je];
#line 701 "sstebz.f"
			}
#line 702 "sstebz.f"
/* L90: */
#line 702 "sstebz.f"
		    }
#line 703 "sstebz.f"
		    iblock[iw] = 0;
#line 704 "sstebz.f"
/* L100: */
#line 704 "sstebz.f"
		}
#line 705 "sstebz.f"
	    }
#line 706 "sstebz.f"
	    if (idiscu > 0) {

#line 708 "sstebz.f"
		wkill = wl;
#line 709 "sstebz.f"
		i__1 = idiscu;
#line 709 "sstebz.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 710 "sstebz.f"
		    iw = 0;
#line 711 "sstebz.f"
		    i__2 = *m;
#line 711 "sstebz.f"
		    for (je = 1; je <= i__2; ++je) {
#line 712 "sstebz.f"
			if (iblock[je] != 0 && (w[je] > wkill || iw == 0)) {
#line 714 "sstebz.f"
			    iw = je;
#line 715 "sstebz.f"
			    wkill = w[je];
#line 716 "sstebz.f"
			}
#line 717 "sstebz.f"
/* L110: */
#line 717 "sstebz.f"
		    }
#line 718 "sstebz.f"
		    iblock[iw] = 0;
#line 719 "sstebz.f"
/* L120: */
#line 719 "sstebz.f"
		}
#line 720 "sstebz.f"
	    }
#line 721 "sstebz.f"
	    im = 0;
#line 722 "sstebz.f"
	    i__1 = *m;
#line 722 "sstebz.f"
	    for (je = 1; je <= i__1; ++je) {
#line 723 "sstebz.f"
		if (iblock[je] != 0) {
#line 724 "sstebz.f"
		    ++im;
#line 725 "sstebz.f"
		    w[im] = w[je];
#line 726 "sstebz.f"
		    iblock[im] = iblock[je];
#line 727 "sstebz.f"
		}
#line 728 "sstebz.f"
/* L130: */
#line 728 "sstebz.f"
	    }
#line 729 "sstebz.f"
	    *m = im;
#line 730 "sstebz.f"
	}
#line 731 "sstebz.f"
	if (idiscl < 0 || idiscu < 0) {
#line 732 "sstebz.f"
	    toofew = TRUE_;
#line 733 "sstebz.f"
	}
#line 734 "sstebz.f"
    }

/*     If ORDER='B', do nothing -- the eigenvalues are already sorted */
/*        by block. */
/*     If ORDER='E', sort the eigenvalues from smallest to largest */

#line 740 "sstebz.f"
    if (iorder == 1 && *nsplit > 1) {
#line 741 "sstebz.f"
	i__1 = *m - 1;
#line 741 "sstebz.f"
	for (je = 1; je <= i__1; ++je) {
#line 742 "sstebz.f"
	    ie = 0;
#line 743 "sstebz.f"
	    tmp1 = w[je];
#line 744 "sstebz.f"
	    i__2 = *m;
#line 744 "sstebz.f"
	    for (j = je + 1; j <= i__2; ++j) {
#line 745 "sstebz.f"
		if (w[j] < tmp1) {
#line 746 "sstebz.f"
		    ie = j;
#line 747 "sstebz.f"
		    tmp1 = w[j];
#line 748 "sstebz.f"
		}
#line 749 "sstebz.f"
/* L140: */
#line 749 "sstebz.f"
	    }

#line 751 "sstebz.f"
	    if (ie != 0) {
#line 752 "sstebz.f"
		itmp1 = iblock[ie];
#line 753 "sstebz.f"
		w[ie] = w[je];
#line 754 "sstebz.f"
		iblock[ie] = iblock[je];
#line 755 "sstebz.f"
		w[je] = tmp1;
#line 756 "sstebz.f"
		iblock[je] = itmp1;
#line 757 "sstebz.f"
	    }
#line 758 "sstebz.f"
/* L150: */
#line 758 "sstebz.f"
	}
#line 759 "sstebz.f"
    }

#line 761 "sstebz.f"
    *info = 0;
#line 762 "sstebz.f"
    if (ncnvrg) {
#line 762 "sstebz.f"
	++(*info);
#line 762 "sstebz.f"
    }
#line 764 "sstebz.f"
    if (toofew) {
#line 764 "sstebz.f"
	*info += 2;
#line 764 "sstebz.f"
    }
#line 766 "sstebz.f"
    return 0;

/*     End of SSTEBZ */

} /* sstebz_ */

