#line 1 "dstebz.f"
/* dstebz.f -- translated by f2c (version 20100827).
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

#line 1 "dstebz.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__0 = 0;

/* > \brief \b DSTEBZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEBZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstebz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstebz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstebz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, */
/*                          M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          ORDER, RANGE */
/*       INTEGER            IL, INFO, IU, M, N, NSPLIT */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEBZ computes the eigenvalues of a symmetric tridiagonal */
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
/* >          VL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
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
/* >          ABSTOL is DOUBLE PRECISION */
/* >          The absolute tolerance for the eigenvalues.  An eigenvalue */
/* >          (or cluster) is considered to be located if it has been */
/* >          determined to lie in an interval whose width is ABSTOL or */
/* >          less.  If ABSTOL is less than or equal to zero, then ULP*|T| */
/* >          will be used, where |T| means the 1-norm of T. */
/* > */
/* >          Eigenvalues will be computed most accurately when ABSTOL is */
/* >          set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          On exit, the first M elements of W will contain the */
/* >          eigenvalues.  (DSTEBZ may use the remaining N-M elements as */
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
/* >          belongs.  (DSTEBZ may use the remaining N-M elements as */
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
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
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
/* >  RELFAC  DOUBLE PRECISION, default = 2.0e0 */
/* >          The relative tolerance.  An interval (a,b] lies within */
/* >          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|), */
/* >          where "ulp" is the machine precision (distance from 1 to */
/* >          the next larger floating point number.) */
/* > */
/* >  FUDGE   DOUBLE PRECISION, default = 2 */
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
/* Subroutine */ int dstebz_(char *range, char *order, integer *n, doublereal 
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
    extern doublereal dlamch_(char *, ftnlen);
    static integer ibegin;
    extern /* Subroutine */ int dlaebz_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    static integer irange, idiscl;
    static doublereal safemn;
    static integer idumma[1];
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer idiscu, iorder;
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

#line 316 "dstebz.f"
    /* Parameter adjustments */
#line 316 "dstebz.f"
    --iwork;
#line 316 "dstebz.f"
    --work;
#line 316 "dstebz.f"
    --isplit;
#line 316 "dstebz.f"
    --iblock;
#line 316 "dstebz.f"
    --w;
#line 316 "dstebz.f"
    --e;
#line 316 "dstebz.f"
    --d__;
#line 316 "dstebz.f"

#line 316 "dstebz.f"
    /* Function Body */
#line 316 "dstebz.f"
    *info = 0;

/*     Decode RANGE */

#line 320 "dstebz.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 321 "dstebz.f"
	irange = 1;
#line 322 "dstebz.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 323 "dstebz.f"
	irange = 2;
#line 324 "dstebz.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 325 "dstebz.f"
	irange = 3;
#line 326 "dstebz.f"
    } else {
#line 327 "dstebz.f"
	irange = 0;
#line 328 "dstebz.f"
    }

/*     Decode ORDER */

#line 332 "dstebz.f"
    if (lsame_(order, "B", (ftnlen)1, (ftnlen)1)) {
#line 333 "dstebz.f"
	iorder = 2;
#line 334 "dstebz.f"
    } else if (lsame_(order, "E", (ftnlen)1, (ftnlen)1)) {
#line 335 "dstebz.f"
	iorder = 1;
#line 336 "dstebz.f"
    } else {
#line 337 "dstebz.f"
	iorder = 0;
#line 338 "dstebz.f"
    }

/*     Check for Errors */

#line 342 "dstebz.f"
    if (irange <= 0) {
#line 343 "dstebz.f"
	*info = -1;
#line 344 "dstebz.f"
    } else if (iorder <= 0) {
#line 345 "dstebz.f"
	*info = -2;
#line 346 "dstebz.f"
    } else if (*n < 0) {
#line 347 "dstebz.f"
	*info = -3;
#line 348 "dstebz.f"
    } else if (irange == 2) {
#line 349 "dstebz.f"
	if (*vl >= *vu) {
#line 349 "dstebz.f"
	    *info = -5;
#line 349 "dstebz.f"
	}
#line 351 "dstebz.f"
    } else if (irange == 3 && (*il < 1 || *il > max(1,*n))) {
#line 353 "dstebz.f"
	*info = -6;
#line 354 "dstebz.f"
    } else if (irange == 3 && (*iu < min(*n,*il) || *iu > *n)) {
#line 356 "dstebz.f"
	*info = -7;
#line 357 "dstebz.f"
    }

#line 359 "dstebz.f"
    if (*info != 0) {
#line 360 "dstebz.f"
	i__1 = -(*info);
#line 360 "dstebz.f"
	xerbla_("DSTEBZ", &i__1, (ftnlen)6);
#line 361 "dstebz.f"
	return 0;
#line 362 "dstebz.f"
    }

/*     Initialize error flags */

#line 366 "dstebz.f"
    *info = 0;
#line 367 "dstebz.f"
    ncnvrg = FALSE_;
#line 368 "dstebz.f"
    toofew = FALSE_;

/*     Quick return if possible */

#line 372 "dstebz.f"
    *m = 0;
#line 373 "dstebz.f"
    if (*n == 0) {
#line 373 "dstebz.f"
	return 0;
#line 373 "dstebz.f"
    }

/*     Simplifications: */

#line 378 "dstebz.f"
    if (irange == 3 && *il == 1 && *iu == *n) {
#line 378 "dstebz.f"
	irange = 1;
#line 378 "dstebz.f"
    }

/*     Get machine constants */
/*     NB is the minimum vector length for vector bisection, or 0 */
/*     if only scalar is to be done. */

#line 385 "dstebz.f"
    safemn = dlamch_("S", (ftnlen)1);
#line 386 "dstebz.f"
    ulp = dlamch_("P", (ftnlen)1);
#line 387 "dstebz.f"
    rtoli = ulp * 2.;
#line 388 "dstebz.f"
    nb = ilaenv_(&c__1, "DSTEBZ", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 389 "dstebz.f"
    if (nb <= 1) {
#line 389 "dstebz.f"
	nb = 0;
#line 389 "dstebz.f"
    }

/*     Special Case when N=1 */

#line 394 "dstebz.f"
    if (*n == 1) {
#line 395 "dstebz.f"
	*nsplit = 1;
#line 396 "dstebz.f"
	isplit[1] = 1;
#line 397 "dstebz.f"
	if (irange == 2 && (*vl >= d__[1] || *vu < d__[1])) {
#line 398 "dstebz.f"
	    *m = 0;
#line 399 "dstebz.f"
	} else {
#line 400 "dstebz.f"
	    w[1] = d__[1];
#line 401 "dstebz.f"
	    iblock[1] = 1;
#line 402 "dstebz.f"
	    *m = 1;
#line 403 "dstebz.f"
	}
#line 404 "dstebz.f"
	return 0;
#line 405 "dstebz.f"
    }

/*     Compute Splitting Points */

#line 409 "dstebz.f"
    *nsplit = 1;
#line 410 "dstebz.f"
    work[*n] = 0.;
#line 411 "dstebz.f"
    pivmin = 1.;

#line 413 "dstebz.f"
    i__1 = *n;
#line 413 "dstebz.f"
    for (j = 2; j <= i__1; ++j) {
/* Computing 2nd power */
#line 414 "dstebz.f"
	d__1 = e[j - 1];
#line 414 "dstebz.f"
	tmp1 = d__1 * d__1;
/* Computing 2nd power */
#line 415 "dstebz.f"
	d__2 = ulp;
#line 415 "dstebz.f"
	if ((d__1 = d__[j] * d__[j - 1], abs(d__1)) * (d__2 * d__2) + safemn 
		> tmp1) {
#line 416 "dstebz.f"
	    isplit[*nsplit] = j - 1;
#line 417 "dstebz.f"
	    ++(*nsplit);
#line 418 "dstebz.f"
	    work[j - 1] = 0.;
#line 419 "dstebz.f"
	} else {
#line 420 "dstebz.f"
	    work[j - 1] = tmp1;
#line 421 "dstebz.f"
	    pivmin = max(pivmin,tmp1);
#line 422 "dstebz.f"
	}
#line 423 "dstebz.f"
/* L10: */
#line 423 "dstebz.f"
    }
#line 424 "dstebz.f"
    isplit[*nsplit] = *n;
#line 425 "dstebz.f"
    pivmin *= safemn;

/*     Compute Interval and ATOLI */

#line 429 "dstebz.f"
    if (irange == 3) {

/*        RANGE='I': Compute the interval containing eigenvalues */
/*                   IL through IU. */

/*        Compute Gershgorin interval for entire (split) matrix */
/*        and use it as the initial interval */

#line 437 "dstebz.f"
	gu = d__[1];
#line 438 "dstebz.f"
	gl = d__[1];
#line 439 "dstebz.f"
	tmp1 = 0.;

#line 441 "dstebz.f"
	i__1 = *n - 1;
#line 441 "dstebz.f"
	for (j = 1; j <= i__1; ++j) {
#line 442 "dstebz.f"
	    tmp2 = sqrt(work[j]);
/* Computing MAX */
#line 443 "dstebz.f"
	    d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
#line 443 "dstebz.f"
	    gu = max(d__1,d__2);
/* Computing MIN */
#line 444 "dstebz.f"
	    d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
#line 444 "dstebz.f"
	    gl = min(d__1,d__2);
#line 445 "dstebz.f"
	    tmp1 = tmp2;
#line 446 "dstebz.f"
/* L20: */
#line 446 "dstebz.f"
	}

/* Computing MAX */
#line 448 "dstebz.f"
	d__1 = gu, d__2 = d__[*n] + tmp1;
#line 448 "dstebz.f"
	gu = max(d__1,d__2);
/* Computing MIN */
#line 449 "dstebz.f"
	d__1 = gl, d__2 = d__[*n] - tmp1;
#line 449 "dstebz.f"
	gl = min(d__1,d__2);
/* Computing MAX */
#line 450 "dstebz.f"
	d__1 = abs(gl), d__2 = abs(gu);
#line 450 "dstebz.f"
	tnorm = max(d__1,d__2);
#line 451 "dstebz.f"
	gl = gl - tnorm * 2.1 * ulp * *n - pivmin * 4.2000000000000002;
#line 452 "dstebz.f"
	gu = gu + tnorm * 2.1 * ulp * *n + pivmin * 2.1;

/*        Compute Iteration parameters */

#line 456 "dstebz.f"
	itmax = (integer) ((log(tnorm + pivmin) - log(pivmin)) / log(2.)) + 2;
#line 458 "dstebz.f"
	if (*abstol <= 0.) {
#line 459 "dstebz.f"
	    atoli = ulp * tnorm;
#line 460 "dstebz.f"
	} else {
#line 461 "dstebz.f"
	    atoli = *abstol;
#line 462 "dstebz.f"
	}

#line 464 "dstebz.f"
	work[*n + 1] = gl;
#line 465 "dstebz.f"
	work[*n + 2] = gl;
#line 466 "dstebz.f"
	work[*n + 3] = gu;
#line 467 "dstebz.f"
	work[*n + 4] = gu;
#line 468 "dstebz.f"
	work[*n + 5] = gl;
#line 469 "dstebz.f"
	work[*n + 6] = gu;
#line 470 "dstebz.f"
	iwork[1] = -1;
#line 471 "dstebz.f"
	iwork[2] = -1;
#line 472 "dstebz.f"
	iwork[3] = *n + 1;
#line 473 "dstebz.f"
	iwork[4] = *n + 1;
#line 474 "dstebz.f"
	iwork[5] = *il - 1;
#line 475 "dstebz.f"
	iwork[6] = *iu;

#line 477 "dstebz.f"
	dlaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, &pivmin, 
		&d__[1], &e[1], &work[1], &iwork[5], &work[*n + 1], &work[*n 
		+ 5], &iout, &iwork[1], &w[1], &iblock[1], &iinfo);

#line 481 "dstebz.f"
	if (iwork[6] == *iu) {
#line 482 "dstebz.f"
	    wl = work[*n + 1];
#line 483 "dstebz.f"
	    wlu = work[*n + 3];
#line 484 "dstebz.f"
	    nwl = iwork[1];
#line 485 "dstebz.f"
	    wu = work[*n + 4];
#line 486 "dstebz.f"
	    wul = work[*n + 2];
#line 487 "dstebz.f"
	    nwu = iwork[4];
#line 488 "dstebz.f"
	} else {
#line 489 "dstebz.f"
	    wl = work[*n + 2];
#line 490 "dstebz.f"
	    wlu = work[*n + 4];
#line 491 "dstebz.f"
	    nwl = iwork[2];
#line 492 "dstebz.f"
	    wu = work[*n + 3];
#line 493 "dstebz.f"
	    wul = work[*n + 1];
#line 494 "dstebz.f"
	    nwu = iwork[3];
#line 495 "dstebz.f"
	}

#line 497 "dstebz.f"
	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
#line 498 "dstebz.f"
	    *info = 4;
#line 499 "dstebz.f"
	    return 0;
#line 500 "dstebz.f"
	}
#line 501 "dstebz.f"
    } else {

/*        RANGE='A' or 'V' -- Set ATOLI */

/* Computing MAX */
#line 505 "dstebz.f"
	d__3 = abs(d__[1]) + abs(e[1]), d__4 = (d__1 = d__[*n], abs(d__1)) + (
		d__2 = e[*n - 1], abs(d__2));
#line 505 "dstebz.f"
	tnorm = max(d__3,d__4);

#line 508 "dstebz.f"
	i__1 = *n - 1;
#line 508 "dstebz.f"
	for (j = 2; j <= i__1; ++j) {
/* Computing MAX */
#line 509 "dstebz.f"
	    d__4 = tnorm, d__5 = (d__1 = d__[j], abs(d__1)) + (d__2 = e[j - 1]
		    , abs(d__2)) + (d__3 = e[j], abs(d__3));
#line 509 "dstebz.f"
	    tnorm = max(d__4,d__5);
#line 511 "dstebz.f"
/* L30: */
#line 511 "dstebz.f"
	}

#line 513 "dstebz.f"
	if (*abstol <= 0.) {
#line 514 "dstebz.f"
	    atoli = ulp * tnorm;
#line 515 "dstebz.f"
	} else {
#line 516 "dstebz.f"
	    atoli = *abstol;
#line 517 "dstebz.f"
	}

#line 519 "dstebz.f"
	if (irange == 2) {
#line 520 "dstebz.f"
	    wl = *vl;
#line 521 "dstebz.f"
	    wu = *vu;
#line 522 "dstebz.f"
	} else {
#line 523 "dstebz.f"
	    wl = 0.;
#line 524 "dstebz.f"
	    wu = 0.;
#line 525 "dstebz.f"
	}
#line 526 "dstebz.f"
    }

/*     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU. */
/*     NWL accumulates the number of eigenvalues .le. WL, */
/*     NWU accumulates the number of eigenvalues .le. WU */

#line 532 "dstebz.f"
    *m = 0;
#line 533 "dstebz.f"
    iend = 0;
#line 534 "dstebz.f"
    *info = 0;
#line 535 "dstebz.f"
    nwl = 0;
#line 536 "dstebz.f"
    nwu = 0;

#line 538 "dstebz.f"
    i__1 = *nsplit;
#line 538 "dstebz.f"
    for (jb = 1; jb <= i__1; ++jb) {
#line 539 "dstebz.f"
	ioff = iend;
#line 540 "dstebz.f"
	ibegin = ioff + 1;
#line 541 "dstebz.f"
	iend = isplit[jb];
#line 542 "dstebz.f"
	in = iend - ioff;

#line 544 "dstebz.f"
	if (in == 1) {

/*           Special Case -- IN=1 */

#line 548 "dstebz.f"
	    if (irange == 1 || wl >= d__[ibegin] - pivmin) {
#line 548 "dstebz.f"
		++nwl;
#line 548 "dstebz.f"
	    }
#line 550 "dstebz.f"
	    if (irange == 1 || wu >= d__[ibegin] - pivmin) {
#line 550 "dstebz.f"
		++nwu;
#line 550 "dstebz.f"
	    }
#line 552 "dstebz.f"
	    if (irange == 1 || wl < d__[ibegin] - pivmin && wu >= d__[ibegin] 
		    - pivmin) {
#line 554 "dstebz.f"
		++(*m);
#line 555 "dstebz.f"
		w[*m] = d__[ibegin];
#line 556 "dstebz.f"
		iblock[*m] = jb;
#line 557 "dstebz.f"
	    }
#line 558 "dstebz.f"
	} else {

/*           General Case -- IN > 1 */

/*           Compute Gershgorin Interval */
/*           and use it as the initial interval */

#line 565 "dstebz.f"
	    gu = d__[ibegin];
#line 566 "dstebz.f"
	    gl = d__[ibegin];
#line 567 "dstebz.f"
	    tmp1 = 0.;

#line 569 "dstebz.f"
	    i__2 = iend - 1;
#line 569 "dstebz.f"
	    for (j = ibegin; j <= i__2; ++j) {
#line 570 "dstebz.f"
		tmp2 = (d__1 = e[j], abs(d__1));
/* Computing MAX */
#line 571 "dstebz.f"
		d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
#line 571 "dstebz.f"
		gu = max(d__1,d__2);
/* Computing MIN */
#line 572 "dstebz.f"
		d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
#line 572 "dstebz.f"
		gl = min(d__1,d__2);
#line 573 "dstebz.f"
		tmp1 = tmp2;
#line 574 "dstebz.f"
/* L40: */
#line 574 "dstebz.f"
	    }

/* Computing MAX */
#line 576 "dstebz.f"
	    d__1 = gu, d__2 = d__[iend] + tmp1;
#line 576 "dstebz.f"
	    gu = max(d__1,d__2);
/* Computing MIN */
#line 577 "dstebz.f"
	    d__1 = gl, d__2 = d__[iend] - tmp1;
#line 577 "dstebz.f"
	    gl = min(d__1,d__2);
/* Computing MAX */
#line 578 "dstebz.f"
	    d__1 = abs(gl), d__2 = abs(gu);
#line 578 "dstebz.f"
	    bnorm = max(d__1,d__2);
#line 579 "dstebz.f"
	    gl = gl - bnorm * 2.1 * ulp * in - pivmin * 2.1;
#line 580 "dstebz.f"
	    gu = gu + bnorm * 2.1 * ulp * in + pivmin * 2.1;

/*           Compute ATOLI for the current submatrix */

#line 584 "dstebz.f"
	    if (*abstol <= 0.) {
/* Computing MAX */
#line 585 "dstebz.f"
		d__1 = abs(gl), d__2 = abs(gu);
#line 585 "dstebz.f"
		atoli = ulp * max(d__1,d__2);
#line 586 "dstebz.f"
	    } else {
#line 587 "dstebz.f"
		atoli = *abstol;
#line 588 "dstebz.f"
	    }

#line 590 "dstebz.f"
	    if (irange > 1) {
#line 591 "dstebz.f"
		if (gu < wl) {
#line 592 "dstebz.f"
		    nwl += in;
#line 593 "dstebz.f"
		    nwu += in;
#line 594 "dstebz.f"
		    goto L70;
#line 595 "dstebz.f"
		}
#line 596 "dstebz.f"
		gl = max(gl,wl);
#line 597 "dstebz.f"
		gu = min(gu,wu);
#line 598 "dstebz.f"
		if (gl >= gu) {
#line 598 "dstebz.f"
		    goto L70;
#line 598 "dstebz.f"
		}
#line 600 "dstebz.f"
	    }

/*           Set Up Initial Interval */

#line 604 "dstebz.f"
	    work[*n + 1] = gl;
#line 605 "dstebz.f"
	    work[*n + in + 1] = gu;
#line 606 "dstebz.f"
	    dlaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);

#line 611 "dstebz.f"
	    nwl += iwork[1];
#line 612 "dstebz.f"
	    nwu += iwork[in + 1];
#line 613 "dstebz.f"
	    iwoff = *m - iwork[1];

/*           Compute Eigenvalues */

#line 617 "dstebz.f"
	    itmax = (integer) ((log(gu - gl + pivmin) - log(pivmin)) / log(2.)
		    ) + 2;
#line 619 "dstebz.f"
	    dlaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, &
		    pivmin, &d__[ibegin], &e[ibegin], &work[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);

/*           Copy Eigenvalues Into W and IBLOCK */
/*           Use -JB for block number for unconverged eigenvalues. */

#line 627 "dstebz.f"
	    i__2 = iout;
#line 627 "dstebz.f"
	    for (j = 1; j <= i__2; ++j) {
#line 628 "dstebz.f"
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;

/*              Flag non-convergence. */

#line 632 "dstebz.f"
		if (j > iout - iinfo) {
#line 633 "dstebz.f"
		    ncnvrg = TRUE_;
#line 634 "dstebz.f"
		    ib = -jb;
#line 635 "dstebz.f"
		} else {
#line 636 "dstebz.f"
		    ib = jb;
#line 637 "dstebz.f"
		}
#line 638 "dstebz.f"
		i__3 = iwork[j + in] + iwoff;
#line 638 "dstebz.f"
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
#line 640 "dstebz.f"
		    w[je] = tmp1;
#line 641 "dstebz.f"
		    iblock[je] = ib;
#line 642 "dstebz.f"
/* L50: */
#line 642 "dstebz.f"
		}
#line 643 "dstebz.f"
/* L60: */
#line 643 "dstebz.f"
	    }

#line 645 "dstebz.f"
	    *m += im;
#line 646 "dstebz.f"
	}
#line 647 "dstebz.f"
L70:
#line 647 "dstebz.f"
	;
#line 647 "dstebz.f"
    }

/*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU */
/*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues. */

#line 652 "dstebz.f"
    if (irange == 3) {
#line 653 "dstebz.f"
	im = 0;
#line 654 "dstebz.f"
	idiscl = *il - 1 - nwl;
#line 655 "dstebz.f"
	idiscu = nwu - *iu;

#line 657 "dstebz.f"
	if (idiscl > 0 || idiscu > 0) {
#line 658 "dstebz.f"
	    i__1 = *m;
#line 658 "dstebz.f"
	    for (je = 1; je <= i__1; ++je) {
#line 659 "dstebz.f"
		if (w[je] <= wlu && idiscl > 0) {
#line 660 "dstebz.f"
		    --idiscl;
#line 661 "dstebz.f"
		} else if (w[je] >= wul && idiscu > 0) {
#line 662 "dstebz.f"
		    --idiscu;
#line 663 "dstebz.f"
		} else {
#line 664 "dstebz.f"
		    ++im;
#line 665 "dstebz.f"
		    w[im] = w[je];
#line 666 "dstebz.f"
		    iblock[im] = iblock[je];
#line 667 "dstebz.f"
		}
#line 668 "dstebz.f"
/* L80: */
#line 668 "dstebz.f"
	    }
#line 669 "dstebz.f"
	    *m = im;
#line 670 "dstebz.f"
	}
#line 671 "dstebz.f"
	if (idiscl > 0 || idiscu > 0) {

/*           Code to deal with effects of bad arithmetic: */
/*           Some low eigenvalues to be discarded are not in (WL,WLU], */
/*           or high eigenvalues to be discarded are not in (WUL,WU] */
/*           so just kill off the smallest IDISCL/largest IDISCU */
/*           eigenvalues, by simply finding the smallest/largest */
/*           eigenvalue(s). */

/*           (If N(w) is monotone non-decreasing, this should never */
/*               happen.) */

#line 683 "dstebz.f"
	    if (idiscl > 0) {
#line 684 "dstebz.f"
		wkill = wu;
#line 685 "dstebz.f"
		i__1 = idiscl;
#line 685 "dstebz.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 686 "dstebz.f"
		    iw = 0;
#line 687 "dstebz.f"
		    i__2 = *m;
#line 687 "dstebz.f"
		    for (je = 1; je <= i__2; ++je) {
#line 688 "dstebz.f"
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
#line 690 "dstebz.f"
			    iw = je;
#line 691 "dstebz.f"
			    wkill = w[je];
#line 692 "dstebz.f"
			}
#line 693 "dstebz.f"
/* L90: */
#line 693 "dstebz.f"
		    }
#line 694 "dstebz.f"
		    iblock[iw] = 0;
#line 695 "dstebz.f"
/* L100: */
#line 695 "dstebz.f"
		}
#line 696 "dstebz.f"
	    }
#line 697 "dstebz.f"
	    if (idiscu > 0) {

#line 699 "dstebz.f"
		wkill = wl;
#line 700 "dstebz.f"
		i__1 = idiscu;
#line 700 "dstebz.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 701 "dstebz.f"
		    iw = 0;
#line 702 "dstebz.f"
		    i__2 = *m;
#line 702 "dstebz.f"
		    for (je = 1; je <= i__2; ++je) {
#line 703 "dstebz.f"
			if (iblock[je] != 0 && (w[je] > wkill || iw == 0)) {
#line 705 "dstebz.f"
			    iw = je;
#line 706 "dstebz.f"
			    wkill = w[je];
#line 707 "dstebz.f"
			}
#line 708 "dstebz.f"
/* L110: */
#line 708 "dstebz.f"
		    }
#line 709 "dstebz.f"
		    iblock[iw] = 0;
#line 710 "dstebz.f"
/* L120: */
#line 710 "dstebz.f"
		}
#line 711 "dstebz.f"
	    }
#line 712 "dstebz.f"
	    im = 0;
#line 713 "dstebz.f"
	    i__1 = *m;
#line 713 "dstebz.f"
	    for (je = 1; je <= i__1; ++je) {
#line 714 "dstebz.f"
		if (iblock[je] != 0) {
#line 715 "dstebz.f"
		    ++im;
#line 716 "dstebz.f"
		    w[im] = w[je];
#line 717 "dstebz.f"
		    iblock[im] = iblock[je];
#line 718 "dstebz.f"
		}
#line 719 "dstebz.f"
/* L130: */
#line 719 "dstebz.f"
	    }
#line 720 "dstebz.f"
	    *m = im;
#line 721 "dstebz.f"
	}
#line 722 "dstebz.f"
	if (idiscl < 0 || idiscu < 0) {
#line 723 "dstebz.f"
	    toofew = TRUE_;
#line 724 "dstebz.f"
	}
#line 725 "dstebz.f"
    }

/*     If ORDER='B', do nothing -- the eigenvalues are already sorted */
/*        by block. */
/*     If ORDER='E', sort the eigenvalues from smallest to largest */

#line 731 "dstebz.f"
    if (iorder == 1 && *nsplit > 1) {
#line 732 "dstebz.f"
	i__1 = *m - 1;
#line 732 "dstebz.f"
	for (je = 1; je <= i__1; ++je) {
#line 733 "dstebz.f"
	    ie = 0;
#line 734 "dstebz.f"
	    tmp1 = w[je];
#line 735 "dstebz.f"
	    i__2 = *m;
#line 735 "dstebz.f"
	    for (j = je + 1; j <= i__2; ++j) {
#line 736 "dstebz.f"
		if (w[j] < tmp1) {
#line 737 "dstebz.f"
		    ie = j;
#line 738 "dstebz.f"
		    tmp1 = w[j];
#line 739 "dstebz.f"
		}
#line 740 "dstebz.f"
/* L140: */
#line 740 "dstebz.f"
	    }

#line 742 "dstebz.f"
	    if (ie != 0) {
#line 743 "dstebz.f"
		itmp1 = iblock[ie];
#line 744 "dstebz.f"
		w[ie] = w[je];
#line 745 "dstebz.f"
		iblock[ie] = iblock[je];
#line 746 "dstebz.f"
		w[je] = tmp1;
#line 747 "dstebz.f"
		iblock[je] = itmp1;
#line 748 "dstebz.f"
	    }
#line 749 "dstebz.f"
/* L150: */
#line 749 "dstebz.f"
	}
#line 750 "dstebz.f"
    }

#line 752 "dstebz.f"
    *info = 0;
#line 753 "dstebz.f"
    if (ncnvrg) {
#line 753 "dstebz.f"
	++(*info);
#line 753 "dstebz.f"
    }
#line 755 "dstebz.f"
    if (toofew) {
#line 755 "dstebz.f"
	*info += 2;
#line 755 "dstebz.f"
    }
#line 757 "dstebz.f"
    return 0;

/*     End of DSTEBZ */

} /* dstebz_ */

