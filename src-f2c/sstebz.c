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
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* > */
/* >          If RANGE='V', the lower and upper bounds of the interval to */
/* >          be searched for eigenvalues.  Eigenvalues less than or equal */
/* >          to VL, or greater than VU, will not be returned.  VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* > */
/* >          If RANGE='I', the indices (in ascending order) of the */
/* >          smallest and largest eigenvalues to be returned. */
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

/* > \date November 2011 */

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


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 316 "sstebz.f"
    /* Parameter adjustments */
#line 316 "sstebz.f"
    --iwork;
#line 316 "sstebz.f"
    --work;
#line 316 "sstebz.f"
    --isplit;
#line 316 "sstebz.f"
    --iblock;
#line 316 "sstebz.f"
    --w;
#line 316 "sstebz.f"
    --e;
#line 316 "sstebz.f"
    --d__;
#line 316 "sstebz.f"

#line 316 "sstebz.f"
    /* Function Body */
#line 316 "sstebz.f"
    *info = 0;

/*     Decode RANGE */

#line 320 "sstebz.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 321 "sstebz.f"
	irange = 1;
#line 322 "sstebz.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 323 "sstebz.f"
	irange = 2;
#line 324 "sstebz.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 325 "sstebz.f"
	irange = 3;
#line 326 "sstebz.f"
    } else {
#line 327 "sstebz.f"
	irange = 0;
#line 328 "sstebz.f"
    }

/*     Decode ORDER */

#line 332 "sstebz.f"
    if (lsame_(order, "B", (ftnlen)1, (ftnlen)1)) {
#line 333 "sstebz.f"
	iorder = 2;
#line 334 "sstebz.f"
    } else if (lsame_(order, "E", (ftnlen)1, (ftnlen)1)) {
#line 335 "sstebz.f"
	iorder = 1;
#line 336 "sstebz.f"
    } else {
#line 337 "sstebz.f"
	iorder = 0;
#line 338 "sstebz.f"
    }

/*     Check for Errors */

#line 342 "sstebz.f"
    if (irange <= 0) {
#line 343 "sstebz.f"
	*info = -1;
#line 344 "sstebz.f"
    } else if (iorder <= 0) {
#line 345 "sstebz.f"
	*info = -2;
#line 346 "sstebz.f"
    } else if (*n < 0) {
#line 347 "sstebz.f"
	*info = -3;
#line 348 "sstebz.f"
    } else if (irange == 2) {
#line 349 "sstebz.f"
	if (*vl >= *vu) {
#line 349 "sstebz.f"
	    *info = -5;
#line 349 "sstebz.f"
	}
#line 350 "sstebz.f"
    } else if (irange == 3 && (*il < 1 || *il > max(1,*n))) {
#line 352 "sstebz.f"
	*info = -6;
#line 353 "sstebz.f"
    } else if (irange == 3 && (*iu < min(*n,*il) || *iu > *n)) {
#line 355 "sstebz.f"
	*info = -7;
#line 356 "sstebz.f"
    }

#line 358 "sstebz.f"
    if (*info != 0) {
#line 359 "sstebz.f"
	i__1 = -(*info);
#line 359 "sstebz.f"
	xerbla_("SSTEBZ", &i__1, (ftnlen)6);
#line 360 "sstebz.f"
	return 0;
#line 361 "sstebz.f"
    }

/*     Initialize error flags */

#line 365 "sstebz.f"
    *info = 0;
#line 366 "sstebz.f"
    ncnvrg = FALSE_;
#line 367 "sstebz.f"
    toofew = FALSE_;

/*     Quick return if possible */

#line 371 "sstebz.f"
    *m = 0;
#line 372 "sstebz.f"
    if (*n == 0) {
#line 372 "sstebz.f"
	return 0;
#line 372 "sstebz.f"
    }

/*     Simplifications: */

#line 377 "sstebz.f"
    if (irange == 3 && *il == 1 && *iu == *n) {
#line 377 "sstebz.f"
	irange = 1;
#line 377 "sstebz.f"
    }

/*     Get machine constants */
/*     NB is the minimum vector length for vector bisection, or 0 */
/*     if only scalar is to be done. */

#line 384 "sstebz.f"
    safemn = slamch_("S", (ftnlen)1);
#line 385 "sstebz.f"
    ulp = slamch_("P", (ftnlen)1);
#line 386 "sstebz.f"
    rtoli = ulp * 2.;
#line 387 "sstebz.f"
    nb = ilaenv_(&c__1, "SSTEBZ", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 388 "sstebz.f"
    if (nb <= 1) {
#line 388 "sstebz.f"
	nb = 0;
#line 388 "sstebz.f"
    }

/*     Special Case when N=1 */

#line 393 "sstebz.f"
    if (*n == 1) {
#line 394 "sstebz.f"
	*nsplit = 1;
#line 395 "sstebz.f"
	isplit[1] = 1;
#line 396 "sstebz.f"
	if (irange == 2 && (*vl >= d__[1] || *vu < d__[1])) {
#line 397 "sstebz.f"
	    *m = 0;
#line 398 "sstebz.f"
	} else {
#line 399 "sstebz.f"
	    w[1] = d__[1];
#line 400 "sstebz.f"
	    iblock[1] = 1;
#line 401 "sstebz.f"
	    *m = 1;
#line 402 "sstebz.f"
	}
#line 403 "sstebz.f"
	return 0;
#line 404 "sstebz.f"
    }

/*     Compute Splitting Points */

#line 408 "sstebz.f"
    *nsplit = 1;
#line 409 "sstebz.f"
    work[*n] = 0.;
#line 410 "sstebz.f"
    pivmin = 1.;

#line 412 "sstebz.f"
    i__1 = *n;
#line 412 "sstebz.f"
    for (j = 2; j <= i__1; ++j) {
/* Computing 2nd power */
#line 413 "sstebz.f"
	d__1 = e[j - 1];
#line 413 "sstebz.f"
	tmp1 = d__1 * d__1;
/* Computing 2nd power */
#line 414 "sstebz.f"
	d__2 = ulp;
#line 414 "sstebz.f"
	if ((d__1 = d__[j] * d__[j - 1], abs(d__1)) * (d__2 * d__2) + safemn 
		> tmp1) {
#line 415 "sstebz.f"
	    isplit[*nsplit] = j - 1;
#line 416 "sstebz.f"
	    ++(*nsplit);
#line 417 "sstebz.f"
	    work[j - 1] = 0.;
#line 418 "sstebz.f"
	} else {
#line 419 "sstebz.f"
	    work[j - 1] = tmp1;
#line 420 "sstebz.f"
	    pivmin = max(pivmin,tmp1);
#line 421 "sstebz.f"
	}
#line 422 "sstebz.f"
/* L10: */
#line 422 "sstebz.f"
    }
#line 423 "sstebz.f"
    isplit[*nsplit] = *n;
#line 424 "sstebz.f"
    pivmin *= safemn;

/*     Compute Interval and ATOLI */

#line 428 "sstebz.f"
    if (irange == 3) {

/*        RANGE='I': Compute the interval containing eigenvalues */
/*                   IL through IU. */

/*        Compute Gershgorin interval for entire (split) matrix */
/*        and use it as the initial interval */

#line 436 "sstebz.f"
	gu = d__[1];
#line 437 "sstebz.f"
	gl = d__[1];
#line 438 "sstebz.f"
	tmp1 = 0.;

#line 440 "sstebz.f"
	i__1 = *n - 1;
#line 440 "sstebz.f"
	for (j = 1; j <= i__1; ++j) {
#line 441 "sstebz.f"
	    tmp2 = sqrt(work[j]);
/* Computing MAX */
#line 442 "sstebz.f"
	    d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
#line 442 "sstebz.f"
	    gu = max(d__1,d__2);
/* Computing MIN */
#line 443 "sstebz.f"
	    d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
#line 443 "sstebz.f"
	    gl = min(d__1,d__2);
#line 444 "sstebz.f"
	    tmp1 = tmp2;
#line 445 "sstebz.f"
/* L20: */
#line 445 "sstebz.f"
	}

/* Computing MAX */
#line 447 "sstebz.f"
	d__1 = gu, d__2 = d__[*n] + tmp1;
#line 447 "sstebz.f"
	gu = max(d__1,d__2);
/* Computing MIN */
#line 448 "sstebz.f"
	d__1 = gl, d__2 = d__[*n] - tmp1;
#line 448 "sstebz.f"
	gl = min(d__1,d__2);
/* Computing MAX */
#line 449 "sstebz.f"
	d__1 = abs(gl), d__2 = abs(gu);
#line 449 "sstebz.f"
	tnorm = max(d__1,d__2);
#line 450 "sstebz.f"
	gl = gl - tnorm * 2.1 * ulp * *n - pivmin * 4.2000000000000002;
#line 451 "sstebz.f"
	gu = gu + tnorm * 2.1 * ulp * *n + pivmin * 2.1;

/*        Compute Iteration parameters */

#line 455 "sstebz.f"
	itmax = (integer) ((log(tnorm + pivmin) - log(pivmin)) / log(2.)) + 2;
#line 457 "sstebz.f"
	if (*abstol <= 0.) {
#line 458 "sstebz.f"
	    atoli = ulp * tnorm;
#line 459 "sstebz.f"
	} else {
#line 460 "sstebz.f"
	    atoli = *abstol;
#line 461 "sstebz.f"
	}

#line 463 "sstebz.f"
	work[*n + 1] = gl;
#line 464 "sstebz.f"
	work[*n + 2] = gl;
#line 465 "sstebz.f"
	work[*n + 3] = gu;
#line 466 "sstebz.f"
	work[*n + 4] = gu;
#line 467 "sstebz.f"
	work[*n + 5] = gl;
#line 468 "sstebz.f"
	work[*n + 6] = gu;
#line 469 "sstebz.f"
	iwork[1] = -1;
#line 470 "sstebz.f"
	iwork[2] = -1;
#line 471 "sstebz.f"
	iwork[3] = *n + 1;
#line 472 "sstebz.f"
	iwork[4] = *n + 1;
#line 473 "sstebz.f"
	iwork[5] = *il - 1;
#line 474 "sstebz.f"
	iwork[6] = *iu;

#line 476 "sstebz.f"
	slaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, &pivmin, 
		&d__[1], &e[1], &work[1], &iwork[5], &work[*n + 1], &work[*n 
		+ 5], &iout, &iwork[1], &w[1], &iblock[1], &iinfo);

#line 480 "sstebz.f"
	if (iwork[6] == *iu) {
#line 481 "sstebz.f"
	    wl = work[*n + 1];
#line 482 "sstebz.f"
	    wlu = work[*n + 3];
#line 483 "sstebz.f"
	    nwl = iwork[1];
#line 484 "sstebz.f"
	    wu = work[*n + 4];
#line 485 "sstebz.f"
	    wul = work[*n + 2];
#line 486 "sstebz.f"
	    nwu = iwork[4];
#line 487 "sstebz.f"
	} else {
#line 488 "sstebz.f"
	    wl = work[*n + 2];
#line 489 "sstebz.f"
	    wlu = work[*n + 4];
#line 490 "sstebz.f"
	    nwl = iwork[2];
#line 491 "sstebz.f"
	    wu = work[*n + 3];
#line 492 "sstebz.f"
	    wul = work[*n + 1];
#line 493 "sstebz.f"
	    nwu = iwork[3];
#line 494 "sstebz.f"
	}

#line 496 "sstebz.f"
	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
#line 497 "sstebz.f"
	    *info = 4;
#line 498 "sstebz.f"
	    return 0;
#line 499 "sstebz.f"
	}
#line 500 "sstebz.f"
    } else {

/*        RANGE='A' or 'V' -- Set ATOLI */

/* Computing MAX */
#line 504 "sstebz.f"
	d__3 = abs(d__[1]) + abs(e[1]), d__4 = (d__1 = d__[*n], abs(d__1)) + (
		d__2 = e[*n - 1], abs(d__2));
#line 504 "sstebz.f"
	tnorm = max(d__3,d__4);

#line 507 "sstebz.f"
	i__1 = *n - 1;
#line 507 "sstebz.f"
	for (j = 2; j <= i__1; ++j) {
/* Computing MAX */
#line 508 "sstebz.f"
	    d__4 = tnorm, d__5 = (d__1 = d__[j], abs(d__1)) + (d__2 = e[j - 1]
		    , abs(d__2)) + (d__3 = e[j], abs(d__3));
#line 508 "sstebz.f"
	    tnorm = max(d__4,d__5);
#line 510 "sstebz.f"
/* L30: */
#line 510 "sstebz.f"
	}

#line 512 "sstebz.f"
	if (*abstol <= 0.) {
#line 513 "sstebz.f"
	    atoli = ulp * tnorm;
#line 514 "sstebz.f"
	} else {
#line 515 "sstebz.f"
	    atoli = *abstol;
#line 516 "sstebz.f"
	}

#line 518 "sstebz.f"
	if (irange == 2) {
#line 519 "sstebz.f"
	    wl = *vl;
#line 520 "sstebz.f"
	    wu = *vu;
#line 521 "sstebz.f"
	} else {
#line 522 "sstebz.f"
	    wl = 0.;
#line 523 "sstebz.f"
	    wu = 0.;
#line 524 "sstebz.f"
	}
#line 525 "sstebz.f"
    }

/*     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU. */
/*     NWL accumulates the number of eigenvalues .le. WL, */
/*     NWU accumulates the number of eigenvalues .le. WU */

#line 531 "sstebz.f"
    *m = 0;
#line 532 "sstebz.f"
    iend = 0;
#line 533 "sstebz.f"
    *info = 0;
#line 534 "sstebz.f"
    nwl = 0;
#line 535 "sstebz.f"
    nwu = 0;

#line 537 "sstebz.f"
    i__1 = *nsplit;
#line 537 "sstebz.f"
    for (jb = 1; jb <= i__1; ++jb) {
#line 538 "sstebz.f"
	ioff = iend;
#line 539 "sstebz.f"
	ibegin = ioff + 1;
#line 540 "sstebz.f"
	iend = isplit[jb];
#line 541 "sstebz.f"
	in = iend - ioff;

#line 543 "sstebz.f"
	if (in == 1) {

/*           Special Case -- IN=1 */

#line 547 "sstebz.f"
	    if (irange == 1 || wl >= d__[ibegin] - pivmin) {
#line 547 "sstebz.f"
		++nwl;
#line 547 "sstebz.f"
	    }
#line 549 "sstebz.f"
	    if (irange == 1 || wu >= d__[ibegin] - pivmin) {
#line 549 "sstebz.f"
		++nwu;
#line 549 "sstebz.f"
	    }
#line 551 "sstebz.f"
	    if (irange == 1 || wl < d__[ibegin] - pivmin && wu >= d__[ibegin] 
		    - pivmin) {
#line 553 "sstebz.f"
		++(*m);
#line 554 "sstebz.f"
		w[*m] = d__[ibegin];
#line 555 "sstebz.f"
		iblock[*m] = jb;
#line 556 "sstebz.f"
	    }
#line 557 "sstebz.f"
	} else {

/*           General Case -- IN > 1 */

/*           Compute Gershgorin Interval */
/*           and use it as the initial interval */

#line 564 "sstebz.f"
	    gu = d__[ibegin];
#line 565 "sstebz.f"
	    gl = d__[ibegin];
#line 566 "sstebz.f"
	    tmp1 = 0.;

#line 568 "sstebz.f"
	    i__2 = iend - 1;
#line 568 "sstebz.f"
	    for (j = ibegin; j <= i__2; ++j) {
#line 569 "sstebz.f"
		tmp2 = (d__1 = e[j], abs(d__1));
/* Computing MAX */
#line 570 "sstebz.f"
		d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
#line 570 "sstebz.f"
		gu = max(d__1,d__2);
/* Computing MIN */
#line 571 "sstebz.f"
		d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
#line 571 "sstebz.f"
		gl = min(d__1,d__2);
#line 572 "sstebz.f"
		tmp1 = tmp2;
#line 573 "sstebz.f"
/* L40: */
#line 573 "sstebz.f"
	    }

/* Computing MAX */
#line 575 "sstebz.f"
	    d__1 = gu, d__2 = d__[iend] + tmp1;
#line 575 "sstebz.f"
	    gu = max(d__1,d__2);
/* Computing MIN */
#line 576 "sstebz.f"
	    d__1 = gl, d__2 = d__[iend] - tmp1;
#line 576 "sstebz.f"
	    gl = min(d__1,d__2);
/* Computing MAX */
#line 577 "sstebz.f"
	    d__1 = abs(gl), d__2 = abs(gu);
#line 577 "sstebz.f"
	    bnorm = max(d__1,d__2);
#line 578 "sstebz.f"
	    gl = gl - bnorm * 2.1 * ulp * in - pivmin * 2.1;
#line 579 "sstebz.f"
	    gu = gu + bnorm * 2.1 * ulp * in + pivmin * 2.1;

/*           Compute ATOLI for the current submatrix */

#line 583 "sstebz.f"
	    if (*abstol <= 0.) {
/* Computing MAX */
#line 584 "sstebz.f"
		d__1 = abs(gl), d__2 = abs(gu);
#line 584 "sstebz.f"
		atoli = ulp * max(d__1,d__2);
#line 585 "sstebz.f"
	    } else {
#line 586 "sstebz.f"
		atoli = *abstol;
#line 587 "sstebz.f"
	    }

#line 589 "sstebz.f"
	    if (irange > 1) {
#line 590 "sstebz.f"
		if (gu < wl) {
#line 591 "sstebz.f"
		    nwl += in;
#line 592 "sstebz.f"
		    nwu += in;
#line 593 "sstebz.f"
		    goto L70;
#line 594 "sstebz.f"
		}
#line 595 "sstebz.f"
		gl = max(gl,wl);
#line 596 "sstebz.f"
		gu = min(gu,wu);
#line 597 "sstebz.f"
		if (gl >= gu) {
#line 597 "sstebz.f"
		    goto L70;
#line 597 "sstebz.f"
		}
#line 599 "sstebz.f"
	    }

/*           Set Up Initial Interval */

#line 603 "sstebz.f"
	    work[*n + 1] = gl;
#line 604 "sstebz.f"
	    work[*n + in + 1] = gu;
#line 605 "sstebz.f"
	    slaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);

#line 610 "sstebz.f"
	    nwl += iwork[1];
#line 611 "sstebz.f"
	    nwu += iwork[in + 1];
#line 612 "sstebz.f"
	    iwoff = *m - iwork[1];

/*           Compute Eigenvalues */

#line 616 "sstebz.f"
	    itmax = (integer) ((log(gu - gl + pivmin) - log(pivmin)) / log(2.)
		    ) + 2;
#line 618 "sstebz.f"
	    slaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);

/*           Copy Eigenvalues Into W and IBLOCK */
/*           Use -JB for block number for unconverged eigenvalues. */

#line 626 "sstebz.f"
	    i__2 = iout;
#line 626 "sstebz.f"
	    for (j = 1; j <= i__2; ++j) {
#line 627 "sstebz.f"
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;

/*              Flag non-convergence. */

#line 631 "sstebz.f"
		if (j > iout - iinfo) {
#line 632 "sstebz.f"
		    ncnvrg = TRUE_;
#line 633 "sstebz.f"
		    ib = -jb;
#line 634 "sstebz.f"
		} else {
#line 635 "sstebz.f"
		    ib = jb;
#line 636 "sstebz.f"
		}
#line 637 "sstebz.f"
		i__3 = iwork[j + in] + iwoff;
#line 637 "sstebz.f"
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
#line 639 "sstebz.f"
		    w[je] = tmp1;
#line 640 "sstebz.f"
		    iblock[je] = ib;
#line 641 "sstebz.f"
/* L50: */
#line 641 "sstebz.f"
		}
#line 642 "sstebz.f"
/* L60: */
#line 642 "sstebz.f"
	    }

#line 644 "sstebz.f"
	    *m += im;
#line 645 "sstebz.f"
	}
#line 646 "sstebz.f"
L70:
#line 646 "sstebz.f"
	;
#line 646 "sstebz.f"
    }

/*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU */
/*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues. */

#line 651 "sstebz.f"
    if (irange == 3) {
#line 652 "sstebz.f"
	im = 0;
#line 653 "sstebz.f"
	idiscl = *il - 1 - nwl;
#line 654 "sstebz.f"
	idiscu = nwu - *iu;

#line 656 "sstebz.f"
	if (idiscl > 0 || idiscu > 0) {
#line 657 "sstebz.f"
	    i__1 = *m;
#line 657 "sstebz.f"
	    for (je = 1; je <= i__1; ++je) {
#line 658 "sstebz.f"
		if (w[je] <= wlu && idiscl > 0) {
#line 659 "sstebz.f"
		    --idiscl;
#line 660 "sstebz.f"
		} else if (w[je] >= wul && idiscu > 0) {
#line 661 "sstebz.f"
		    --idiscu;
#line 662 "sstebz.f"
		} else {
#line 663 "sstebz.f"
		    ++im;
#line 664 "sstebz.f"
		    w[im] = w[je];
#line 665 "sstebz.f"
		    iblock[im] = iblock[je];
#line 666 "sstebz.f"
		}
#line 667 "sstebz.f"
/* L80: */
#line 667 "sstebz.f"
	    }
#line 668 "sstebz.f"
	    *m = im;
#line 669 "sstebz.f"
	}
#line 670 "sstebz.f"
	if (idiscl > 0 || idiscu > 0) {

/*           Code to deal with effects of bad arithmetic: */
/*           Some low eigenvalues to be discarded are not in (WL,WLU], */
/*           or high eigenvalues to be discarded are not in (WUL,WU] */
/*           so just kill off the smallest IDISCL/largest IDISCU */
/*           eigenvalues, by simply finding the smallest/largest */
/*           eigenvalue(s). */

/*           (If N(w) is monotone non-decreasing, this should never */
/*               happen.) */

#line 682 "sstebz.f"
	    if (idiscl > 0) {
#line 683 "sstebz.f"
		wkill = wu;
#line 684 "sstebz.f"
		i__1 = idiscl;
#line 684 "sstebz.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 685 "sstebz.f"
		    iw = 0;
#line 686 "sstebz.f"
		    i__2 = *m;
#line 686 "sstebz.f"
		    for (je = 1; je <= i__2; ++je) {
#line 687 "sstebz.f"
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
#line 689 "sstebz.f"
			    iw = je;
#line 690 "sstebz.f"
			    wkill = w[je];
#line 691 "sstebz.f"
			}
#line 692 "sstebz.f"
/* L90: */
#line 692 "sstebz.f"
		    }
#line 693 "sstebz.f"
		    iblock[iw] = 0;
#line 694 "sstebz.f"
/* L100: */
#line 694 "sstebz.f"
		}
#line 695 "sstebz.f"
	    }
#line 696 "sstebz.f"
	    if (idiscu > 0) {

#line 698 "sstebz.f"
		wkill = wl;
#line 699 "sstebz.f"
		i__1 = idiscu;
#line 699 "sstebz.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 700 "sstebz.f"
		    iw = 0;
#line 701 "sstebz.f"
		    i__2 = *m;
#line 701 "sstebz.f"
		    for (je = 1; je <= i__2; ++je) {
#line 702 "sstebz.f"
			if (iblock[je] != 0 && (w[je] > wkill || iw == 0)) {
#line 704 "sstebz.f"
			    iw = je;
#line 705 "sstebz.f"
			    wkill = w[je];
#line 706 "sstebz.f"
			}
#line 707 "sstebz.f"
/* L110: */
#line 707 "sstebz.f"
		    }
#line 708 "sstebz.f"
		    iblock[iw] = 0;
#line 709 "sstebz.f"
/* L120: */
#line 709 "sstebz.f"
		}
#line 710 "sstebz.f"
	    }
#line 711 "sstebz.f"
	    im = 0;
#line 712 "sstebz.f"
	    i__1 = *m;
#line 712 "sstebz.f"
	    for (je = 1; je <= i__1; ++je) {
#line 713 "sstebz.f"
		if (iblock[je] != 0) {
#line 714 "sstebz.f"
		    ++im;
#line 715 "sstebz.f"
		    w[im] = w[je];
#line 716 "sstebz.f"
		    iblock[im] = iblock[je];
#line 717 "sstebz.f"
		}
#line 718 "sstebz.f"
/* L130: */
#line 718 "sstebz.f"
	    }
#line 719 "sstebz.f"
	    *m = im;
#line 720 "sstebz.f"
	}
#line 721 "sstebz.f"
	if (idiscl < 0 || idiscu < 0) {
#line 722 "sstebz.f"
	    toofew = TRUE_;
#line 723 "sstebz.f"
	}
#line 724 "sstebz.f"
    }

/*     If ORDER='B', do nothing -- the eigenvalues are already sorted */
/*        by block. */
/*     If ORDER='E', sort the eigenvalues from smallest to largest */

#line 730 "sstebz.f"
    if (iorder == 1 && *nsplit > 1) {
#line 731 "sstebz.f"
	i__1 = *m - 1;
#line 731 "sstebz.f"
	for (je = 1; je <= i__1; ++je) {
#line 732 "sstebz.f"
	    ie = 0;
#line 733 "sstebz.f"
	    tmp1 = w[je];
#line 734 "sstebz.f"
	    i__2 = *m;
#line 734 "sstebz.f"
	    for (j = je + 1; j <= i__2; ++j) {
#line 735 "sstebz.f"
		if (w[j] < tmp1) {
#line 736 "sstebz.f"
		    ie = j;
#line 737 "sstebz.f"
		    tmp1 = w[j];
#line 738 "sstebz.f"
		}
#line 739 "sstebz.f"
/* L140: */
#line 739 "sstebz.f"
	    }

#line 741 "sstebz.f"
	    if (ie != 0) {
#line 742 "sstebz.f"
		itmp1 = iblock[ie];
#line 743 "sstebz.f"
		w[ie] = w[je];
#line 744 "sstebz.f"
		iblock[ie] = iblock[je];
#line 745 "sstebz.f"
		w[je] = tmp1;
#line 746 "sstebz.f"
		iblock[je] = itmp1;
#line 747 "sstebz.f"
	    }
#line 748 "sstebz.f"
/* L150: */
#line 748 "sstebz.f"
	}
#line 749 "sstebz.f"
    }

#line 751 "sstebz.f"
    *info = 0;
#line 752 "sstebz.f"
    if (ncnvrg) {
#line 752 "sstebz.f"
	++(*info);
#line 752 "sstebz.f"
    }
#line 754 "sstebz.f"
    if (toofew) {
#line 754 "sstebz.f"
	*info += 2;
#line 754 "sstebz.f"
    }
#line 756 "sstebz.f"
    return 0;

/*     End of SSTEBZ */

} /* sstebz_ */

