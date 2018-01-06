#line 1 "slarrd.f"
/* slarrd.f -- translated by f2c (version 20100827).
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

#line 1 "slarrd.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__0 = 0;

/* > \brief \b SLARRD computes the eigenvalues of a symmetric tridiagonal matrix to suitable accuracy. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARRD( RANGE, ORDER, N, VL, VU, IL, IU, GERS, */
/*                           RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT, */
/*                           M, W, WERR, WL, WU, IBLOCK, INDEXW, */
/*                           WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          ORDER, RANGE */
/*       INTEGER            IL, INFO, IU, M, N, NSPLIT */
/*       REAL                PIVMIN, RELTOL, VL, VU, WL, WU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), INDEXW( * ), */
/*      $                   ISPLIT( * ), IWORK( * ) */
/*       REAL               D( * ), E( * ), E2( * ), */
/*      $                   GERS( * ), W( * ), WERR( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARRD computes the eigenvalues of a symmetric tridiagonal */
/* > matrix T to suitable accuracy. This is an auxiliary code to be */
/* > called from SSTEMR. */
/* > The user may ask for all eigenvalues, all eigenvalues */
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
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues.  Eigenvalues less than or equal */
/* >          to VL, or greater than VU, will not be returned.  VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for eigenvalues.  Eigenvalues less than or equal */
/* >          to VL, or greater than VU, will not be returned.  VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] GERS */
/* > \verbatim */
/* >          GERS is REAL array, dimension (2*N) */
/* >          The N Gerschgorin intervals (the i-th Gerschgorin interval */
/* >          is (GERS(2*i-1), GERS(2*i)). */
/* > \endverbatim */
/* > */
/* > \param[in] RELTOL */
/* > \verbatim */
/* >          RELTOL is REAL */
/* >          The minimum relative width of an interval.  When an interval */
/* >          is narrower than RELTOL times the larger (in */
/* >          magnitude) endpoint, then it is considered to be */
/* >          sufficiently small, i.e., converged.  Note: this should */
/* >          always be at least radix*machine epsilon. */
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
/* > \param[in] E2 */
/* > \verbatim */
/* >          E2 is REAL array, dimension (N-1) */
/* >          The (n-1) squared off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
/* >          The minimum pivot allowed in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[in] NSPLIT */
/* > \verbatim */
/* >          NSPLIT is INTEGER */
/* >          The number of diagonal blocks in the matrix T. */
/* >          1 <= NSPLIT <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] ISPLIT */
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
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The actual number of eigenvalues found. 0 <= M <= N. */
/* >          (See also the description of INFO=2,3.) */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          On exit, the first M elements of W will contain the */
/* >          eigenvalue approximations. SLARRD computes an interval */
/* >          I_j = (a_j, b_j] that includes eigenvalue j. The eigenvalue */
/* >          approximation is given as the interval midpoint */
/* >          W(j)= ( a_j + b_j)/2. The corresponding error is bounded by */
/* >          WERR(j) = abs( a_j - b_j)/2 */
/* > \endverbatim */
/* > */
/* > \param[out] WERR */
/* > \verbatim */
/* >          WERR is REAL array, dimension (N) */
/* >          The error bound on the corresponding eigenvalue approximation */
/* >          in W. */
/* > \endverbatim */
/* > */
/* > \param[out] WL */
/* > \verbatim */
/* >          WL is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] WU */
/* > \verbatim */
/* >          WU is REAL */
/* >          The interval (WL, WU] contains all the wanted eigenvalues. */
/* >          If RANGE='V', then WL=VL and WU=VU. */
/* >          If RANGE='A', then WL and WU are the global Gerschgorin bounds */
/* >                        on the spectrum. */
/* >          If RANGE='I', then WL and WU are computed by SLAEBZ from the */
/* >                        index range specified. */
/* > \endverbatim */
/* > */
/* > \param[out] IBLOCK */
/* > \verbatim */
/* >          IBLOCK is INTEGER array, dimension (N) */
/* >          At each row/column j where E(j) is zero or small, the */
/* >          matrix T is considered to split into a block diagonal */
/* >          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which */
/* >          block (from 1 to the number of blocks) the eigenvalue W(i) */
/* >          belongs.  (SLARRD may use the remaining N-M elements as */
/* >          workspace.) */
/* > \endverbatim */
/* > */
/* > \param[out] INDEXW */
/* > \verbatim */
/* >          INDEXW is INTEGER array, dimension (N) */
/* >          The indices of the eigenvalues within each block (submatrix); */
/* >          for example, INDEXW(i)= j and IBLOCK(i)=k imply that the */
/* >          i-th eigenvalue W(i) is the j-th eigenvalue in block k. */
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
/* >  FUDGE   REAL, default = 2 */
/* >          A "fudge factor" to widen the Gershgorin intervals.  Ideally, */
/* >          a value of 1 should work, but on machines with sloppy */
/* >          arithmetic, this needs to be larger.  The default for */
/* >          publicly released versions should be large enough to handle */
/* >          the worst machine around.  Note that this has no effect */
/* >          on accuracy of the solution. */
/* > \endverbatim */
/* > */
/* > \par Contributors: */
/*  ================== */
/* > */
/* >     W. Kahan, University of California, Berkeley, USA \n */
/* >     Beresford Parlett, University of California, Berkeley, USA \n */
/* >     Jim Demmel, University of California, Berkeley, USA \n */
/* >     Inderjit Dhillon, University of Texas, Austin, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */
/* >     Christof Voemel, University of California, Berkeley, USA \n */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slarrd_(char *range, char *order, integer *n, doublereal 
	*vl, doublereal *vu, integer *il, integer *iu, doublereal *gers, 
	doublereal *reltol, doublereal *d__, doublereal *e, doublereal *e2, 
	doublereal *pivmin, integer *nsplit, integer *isplit, integer *m, 
	doublereal *w, doublereal *werr, doublereal *wl, doublereal *wu, 
	integer *iblock, integer *indexw, doublereal *work, integer *iwork, 
	integer *info, ftnlen range_len, ftnlen order_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, j, ib, ie, je, nb;
    static doublereal gl;
    static integer im, in;
    static doublereal gu;
    static integer iw, jee;
    static doublereal eps;
    static integer nwl;
    static doublereal wlu, wul;
    static integer nwu;
    static doublereal tmp1, tmp2;
    static integer iend, jblk, ioff, iout, itmp1, itmp2, jdisc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static doublereal atoli;
    static integer iwoff, itmax;
    static doublereal wkill, rtoli, uflow, tnorm;
    static integer ibegin, irange, idiscl;
    extern doublereal slamch_(char *, ftnlen);
    static integer idumma[1];
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer idiscu;
    extern /* Subroutine */ int slaebz_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    static logical ncnvrg, toofew;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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

#line 386 "slarrd.f"
    /* Parameter adjustments */
#line 386 "slarrd.f"
    --iwork;
#line 386 "slarrd.f"
    --work;
#line 386 "slarrd.f"
    --indexw;
#line 386 "slarrd.f"
    --iblock;
#line 386 "slarrd.f"
    --werr;
#line 386 "slarrd.f"
    --w;
#line 386 "slarrd.f"
    --isplit;
#line 386 "slarrd.f"
    --e2;
#line 386 "slarrd.f"
    --e;
#line 386 "slarrd.f"
    --d__;
#line 386 "slarrd.f"
    --gers;
#line 386 "slarrd.f"

#line 386 "slarrd.f"
    /* Function Body */
#line 386 "slarrd.f"
    *info = 0;

/*     Decode RANGE */

#line 390 "slarrd.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 391 "slarrd.f"
	irange = 1;
#line 392 "slarrd.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 393 "slarrd.f"
	irange = 2;
#line 394 "slarrd.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 395 "slarrd.f"
	irange = 3;
#line 396 "slarrd.f"
    } else {
#line 397 "slarrd.f"
	irange = 0;
#line 398 "slarrd.f"
    }

/*     Check for Errors */

#line 402 "slarrd.f"
    if (irange <= 0) {
#line 403 "slarrd.f"
	*info = -1;
#line 404 "slarrd.f"
    } else if (! (lsame_(order, "B", (ftnlen)1, (ftnlen)1) || lsame_(order, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 405 "slarrd.f"
	*info = -2;
#line 406 "slarrd.f"
    } else if (*n < 0) {
#line 407 "slarrd.f"
	*info = -3;
#line 408 "slarrd.f"
    } else if (irange == 2) {
#line 409 "slarrd.f"
	if (*vl >= *vu) {
#line 409 "slarrd.f"
	    *info = -5;
#line 409 "slarrd.f"
	}
#line 411 "slarrd.f"
    } else if (irange == 3 && (*il < 1 || *il > max(1,*n))) {
#line 413 "slarrd.f"
	*info = -6;
#line 414 "slarrd.f"
    } else if (irange == 3 && (*iu < min(*n,*il) || *iu > *n)) {
#line 416 "slarrd.f"
	*info = -7;
#line 417 "slarrd.f"
    }

#line 419 "slarrd.f"
    if (*info != 0) {
#line 420 "slarrd.f"
	return 0;
#line 421 "slarrd.f"
    }
/*     Initialize error flags */
#line 424 "slarrd.f"
    *info = 0;
#line 425 "slarrd.f"
    ncnvrg = FALSE_;
#line 426 "slarrd.f"
    toofew = FALSE_;
/*     Quick return if possible */
#line 429 "slarrd.f"
    *m = 0;
#line 430 "slarrd.f"
    if (*n == 0) {
#line 430 "slarrd.f"
	return 0;
#line 430 "slarrd.f"
    }
/*     Simplification: */
#line 433 "slarrd.f"
    if (irange == 3 && *il == 1 && *iu == *n) {
#line 433 "slarrd.f"
	irange = 1;
#line 433 "slarrd.f"
    }
/*     Get machine constants */
#line 436 "slarrd.f"
    eps = slamch_("P", (ftnlen)1);
#line 437 "slarrd.f"
    uflow = slamch_("U", (ftnlen)1);
/*     Special Case when N=1 */
/*     Treat case of 1x1 matrix for quick return */
#line 442 "slarrd.f"
    if (*n == 1) {
#line 443 "slarrd.f"
	if (irange == 1 || irange == 2 && d__[1] > *vl && d__[1] <= *vu || 
		irange == 3 && *il == 1 && *iu == 1) {
#line 446 "slarrd.f"
	    *m = 1;
#line 447 "slarrd.f"
	    w[1] = d__[1];
/*           The computation error of the eigenvalue is zero */
#line 449 "slarrd.f"
	    werr[1] = 0.;
#line 450 "slarrd.f"
	    iblock[1] = 1;
#line 451 "slarrd.f"
	    indexw[1] = 1;
#line 452 "slarrd.f"
	}
#line 453 "slarrd.f"
	return 0;
#line 454 "slarrd.f"
    }
/*     NB is the minimum vector length for vector bisection, or 0 */
/*     if only scalar is to be done. */
#line 458 "slarrd.f"
    nb = ilaenv_(&c__1, "SSTEBZ", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 459 "slarrd.f"
    if (nb <= 1) {
#line 459 "slarrd.f"
	nb = 0;
#line 459 "slarrd.f"
    }
/*     Find global spectral radius */
#line 462 "slarrd.f"
    gl = d__[1];
#line 463 "slarrd.f"
    gu = d__[1];
#line 464 "slarrd.f"
    i__1 = *n;
#line 464 "slarrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 465 "slarrd.f"
	d__1 = gl, d__2 = gers[(i__ << 1) - 1];
#line 465 "slarrd.f"
	gl = min(d__1,d__2);
/* Computing MAX */
#line 466 "slarrd.f"
	d__1 = gu, d__2 = gers[i__ * 2];
#line 466 "slarrd.f"
	gu = max(d__1,d__2);
#line 467 "slarrd.f"
/* L5: */
#line 467 "slarrd.f"
    }
/*     Compute global Gerschgorin bounds and spectral diameter */
/* Computing MAX */
#line 469 "slarrd.f"
    d__1 = abs(gl), d__2 = abs(gu);
#line 469 "slarrd.f"
    tnorm = max(d__1,d__2);
#line 470 "slarrd.f"
    gl = gl - tnorm * 2. * eps * *n - *pivmin * 4.;
#line 471 "slarrd.f"
    gu = gu + tnorm * 2. * eps * *n + *pivmin * 4.;
/*     [JAN/28/2009] remove the line below since SPDIAM variable not use */
/*     SPDIAM = GU - GL */
/*     Input arguments for SLAEBZ: */
/*     The relative tolerance.  An interval (a,b] lies within */
/*     "relative tolerance" if  b-a < RELTOL*max(|a|,|b|), */
#line 477 "slarrd.f"
    rtoli = *reltol;
/*     Set the absolute tolerance for interval convergence to zero to force */
/*     interval convergence based on relative size of the interval. */
/*     This is dangerous because intervals might not converge when RELTOL is */
/*     small. But at least a very small number should be selected so that for */
/*     strongly graded matrices, the code can get relatively accurate */
/*     eigenvalues. */
#line 484 "slarrd.f"
    atoli = uflow * 4. + *pivmin * 4.;
#line 486 "slarrd.f"
    if (irange == 3) {
/*        RANGE='I': Compute an interval containing eigenvalues */
/*        IL through IU. The initial interval [GL,GU] from the global */
/*        Gerschgorin bounds GL and GU is refined by SLAEBZ. */
#line 491 "slarrd.f"
	itmax = (integer) ((log(tnorm + *pivmin) - log(*pivmin)) / log(2.)) + 
		2;
#line 493 "slarrd.f"
	work[*n + 1] = gl;
#line 494 "slarrd.f"
	work[*n + 2] = gl;
#line 495 "slarrd.f"
	work[*n + 3] = gu;
#line 496 "slarrd.f"
	work[*n + 4] = gu;
#line 497 "slarrd.f"
	work[*n + 5] = gl;
#line 498 "slarrd.f"
	work[*n + 6] = gu;
#line 499 "slarrd.f"
	iwork[1] = -1;
#line 500 "slarrd.f"
	iwork[2] = -1;
#line 501 "slarrd.f"
	iwork[3] = *n + 1;
#line 502 "slarrd.f"
	iwork[4] = *n + 1;
#line 503 "slarrd.f"
	iwork[5] = *il - 1;
#line 504 "slarrd.f"
	iwork[6] = *iu;

#line 506 "slarrd.f"
	slaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, pivmin, &
		d__[1], &e[1], &e2[1], &iwork[5], &work[*n + 1], &work[*n + 5]
		, &iout, &iwork[1], &w[1], &iblock[1], &iinfo);
#line 509 "slarrd.f"
	if (iinfo != 0) {
#line 510 "slarrd.f"
	    *info = iinfo;
#line 511 "slarrd.f"
	    return 0;
#line 512 "slarrd.f"
	}
/*        On exit, output intervals may not be ordered by ascending negcount */
#line 514 "slarrd.f"
	if (iwork[6] == *iu) {
#line 515 "slarrd.f"
	    *wl = work[*n + 1];
#line 516 "slarrd.f"
	    wlu = work[*n + 3];
#line 517 "slarrd.f"
	    nwl = iwork[1];
#line 518 "slarrd.f"
	    *wu = work[*n + 4];
#line 519 "slarrd.f"
	    wul = work[*n + 2];
#line 520 "slarrd.f"
	    nwu = iwork[4];
#line 521 "slarrd.f"
	} else {
#line 522 "slarrd.f"
	    *wl = work[*n + 2];
#line 523 "slarrd.f"
	    wlu = work[*n + 4];
#line 524 "slarrd.f"
	    nwl = iwork[2];
#line 525 "slarrd.f"
	    *wu = work[*n + 3];
#line 526 "slarrd.f"
	    wul = work[*n + 1];
#line 527 "slarrd.f"
	    nwu = iwork[3];
#line 528 "slarrd.f"
	}
/*        On exit, the interval [WL, WLU] contains a value with negcount NWL, */
/*        and [WUL, WU] contains a value with negcount NWU. */
#line 531 "slarrd.f"
	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
#line 532 "slarrd.f"
	    *info = 4;
#line 533 "slarrd.f"
	    return 0;
#line 534 "slarrd.f"
	}
#line 536 "slarrd.f"
    } else if (irange == 2) {
#line 537 "slarrd.f"
	*wl = *vl;
#line 538 "slarrd.f"
	*wu = *vu;
#line 540 "slarrd.f"
    } else if (irange == 1) {
#line 541 "slarrd.f"
	*wl = gl;
#line 542 "slarrd.f"
	*wu = gu;
#line 543 "slarrd.f"
    }
/*     Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU. */
/*     NWL accumulates the number of eigenvalues .le. WL, */
/*     NWU accumulates the number of eigenvalues .le. WU */
#line 550 "slarrd.f"
    *m = 0;
#line 551 "slarrd.f"
    iend = 0;
#line 552 "slarrd.f"
    *info = 0;
#line 553 "slarrd.f"
    nwl = 0;
#line 554 "slarrd.f"
    nwu = 0;

#line 556 "slarrd.f"
    i__1 = *nsplit;
#line 556 "slarrd.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 557 "slarrd.f"
	ioff = iend;
#line 558 "slarrd.f"
	ibegin = ioff + 1;
#line 559 "slarrd.f"
	iend = isplit[jblk];
#line 560 "slarrd.f"
	in = iend - ioff;

#line 562 "slarrd.f"
	if (in == 1) {
/*           1x1 block */
#line 564 "slarrd.f"
	    if (*wl >= d__[ibegin] - *pivmin) {
#line 564 "slarrd.f"
		++nwl;
#line 564 "slarrd.f"
	    }
#line 566 "slarrd.f"
	    if (*wu >= d__[ibegin] - *pivmin) {
#line 566 "slarrd.f"
		++nwu;
#line 566 "slarrd.f"
	    }
#line 568 "slarrd.f"
	    if (irange == 1 || *wl < d__[ibegin] - *pivmin && *wu >= d__[
		    ibegin] - *pivmin) {
#line 571 "slarrd.f"
		++(*m);
#line 572 "slarrd.f"
		w[*m] = d__[ibegin];
#line 573 "slarrd.f"
		werr[*m] = 0.;
/*              The gap for a single block doesn't matter for the later */
/*              algorithm and is assigned an arbitrary large value */
#line 576 "slarrd.f"
		iblock[*m] = jblk;
#line 577 "slarrd.f"
		indexw[*m] = 1;
#line 578 "slarrd.f"
	    }
/*        Disabled 2x2 case because of a failure on the following matrix */
/*        RANGE = 'I', IL = IU = 4 */
/*          Original Tridiagonal, d = [ */
/*           -0.150102010615740E+00 */
/*           -0.849897989384260E+00 */
/*           -0.128208148052635E-15 */
/*            0.128257718286320E-15 */
/*          ]; */
/*          e = [ */
/*           -0.357171383266986E+00 */
/*           -0.180411241501588E-15 */
/*           -0.175152352710251E-15 */
/*          ]; */

/*         ELSE IF( IN.EQ.2 ) THEN */
/* *           2x2 block */
/*            DISC = SQRT( (HALF*(D(IBEGIN)-D(IEND)))**2 + E(IBEGIN)**2 ) */
/*            TMP1 = HALF*(D(IBEGIN)+D(IEND)) */
/*            L1 = TMP1 - DISC */
/*            IF( WL.GE. L1-PIVMIN ) */
/*     $         NWL = NWL + 1 */
/*            IF( WU.GE. L1-PIVMIN ) */
/*     $         NWU = NWU + 1 */
/*            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L1-PIVMIN .AND. WU.GE. */
/*     $          L1-PIVMIN ) ) THEN */
/*               M = M + 1 */
/*               W( M ) = L1 */
/* *              The uncertainty of eigenvalues of a 2x2 matrix is very small */
/*               WERR( M ) = EPS * ABS( W( M ) ) * TWO */
/*               IBLOCK( M ) = JBLK */
/*               INDEXW( M ) = 1 */
/*            ENDIF */
/*            L2 = TMP1 + DISC */
/*            IF( WL.GE. L2-PIVMIN ) */
/*     $         NWL = NWL + 1 */
/*            IF( WU.GE. L2-PIVMIN ) */
/*     $         NWU = NWU + 1 */
/*            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L2-PIVMIN .AND. WU.GE. */
/*     $          L2-PIVMIN ) ) THEN */
/*               M = M + 1 */
/*               W( M ) = L2 */
/* *              The uncertainty of eigenvalues of a 2x2 matrix is very small */
/*               WERR( M ) = EPS * ABS( W( M ) ) * TWO */
/*               IBLOCK( M ) = JBLK */
/*               INDEXW( M ) = 2 */
/*            ENDIF */
#line 626 "slarrd.f"
	} else {
/*           General Case - block of size IN >= 2 */
/*           Compute local Gerschgorin interval and use it as the initial */
/*           interval for SLAEBZ */
#line 630 "slarrd.f"
	    gu = d__[ibegin];
#line 631 "slarrd.f"
	    gl = d__[ibegin];
#line 632 "slarrd.f"
	    tmp1 = 0.;
#line 634 "slarrd.f"
	    i__2 = iend;
#line 634 "slarrd.f"
	    for (j = ibegin; j <= i__2; ++j) {
/* Computing MIN */
#line 635 "slarrd.f"
		d__1 = gl, d__2 = gers[(j << 1) - 1];
#line 635 "slarrd.f"
		gl = min(d__1,d__2);
/* Computing MAX */
#line 636 "slarrd.f"
		d__1 = gu, d__2 = gers[j * 2];
#line 636 "slarrd.f"
		gu = max(d__1,d__2);
#line 637 "slarrd.f"
/* L40: */
#line 637 "slarrd.f"
	    }
/*           [JAN/28/2009] */
/*           change SPDIAM by TNORM in lines 2 and 3 thereafter */
/*           line 1: remove computation of SPDIAM (not useful anymore) */
/*           SPDIAM = GU - GL */
/*           GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN */
/*           GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN */
#line 644 "slarrd.f"
	    gl = gl - tnorm * 2. * eps * in - *pivmin * 2.;
#line 645 "slarrd.f"
	    gu = gu + tnorm * 2. * eps * in + *pivmin * 2.;

#line 647 "slarrd.f"
	    if (irange > 1) {
#line 648 "slarrd.f"
		if (gu < *wl) {
/*                 the local block contains none of the wanted eigenvalues */
#line 650 "slarrd.f"
		    nwl += in;
#line 651 "slarrd.f"
		    nwu += in;
#line 652 "slarrd.f"
		    goto L70;
#line 653 "slarrd.f"
		}
/*              refine search interval if possible, only range (WL,WU] matters */
#line 655 "slarrd.f"
		gl = max(gl,*wl);
#line 656 "slarrd.f"
		gu = min(gu,*wu);
#line 657 "slarrd.f"
		if (gl >= gu) {
#line 657 "slarrd.f"
		    goto L70;
#line 657 "slarrd.f"
		}
#line 659 "slarrd.f"
	    }
/*           Find negcount of initial interval boundaries GL and GU */
#line 662 "slarrd.f"
	    work[*n + 1] = gl;
#line 663 "slarrd.f"
	    work[*n + in + 1] = gu;
#line 664 "slarrd.f"
	    slaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, 
		    pivmin, &d__[ibegin], &e[ibegin], &e2[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);
#line 668 "slarrd.f"
	    if (iinfo != 0) {
#line 669 "slarrd.f"
		*info = iinfo;
#line 670 "slarrd.f"
		return 0;
#line 671 "slarrd.f"
	    }

#line 673 "slarrd.f"
	    nwl += iwork[1];
#line 674 "slarrd.f"
	    nwu += iwork[in + 1];
#line 675 "slarrd.f"
	    iwoff = *m - iwork[1];
/*           Compute Eigenvalues */
#line 678 "slarrd.f"
	    itmax = (integer) ((log(gu - gl + *pivmin) - log(*pivmin)) / log(
		    2.)) + 2;
#line 680 "slarrd.f"
	    slaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, 
		    pivmin, &d__[ibegin], &e[ibegin], &e2[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);
#line 684 "slarrd.f"
	    if (iinfo != 0) {
#line 685 "slarrd.f"
		*info = iinfo;
#line 686 "slarrd.f"
		return 0;
#line 687 "slarrd.f"
	    }

/*           Copy eigenvalues into W and IBLOCK */
/*           Use -JBLK for block number for unconverged eigenvalues. */
/*           Loop over the number of output intervals from SLAEBZ */
#line 692 "slarrd.f"
	    i__2 = iout;
#line 692 "slarrd.f"
	    for (j = 1; j <= i__2; ++j) {
/*              eigenvalue approximation is middle point of interval */
#line 694 "slarrd.f"
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;
/*              semi length of error interval */
#line 696 "slarrd.f"
		tmp2 = (d__1 = work[j + *n] - work[j + in + *n], abs(d__1)) * 
			.5;
#line 697 "slarrd.f"
		if (j > iout - iinfo) {
/*                 Flag non-convergence. */
#line 699 "slarrd.f"
		    ncnvrg = TRUE_;
#line 700 "slarrd.f"
		    ib = -jblk;
#line 701 "slarrd.f"
		} else {
#line 702 "slarrd.f"
		    ib = jblk;
#line 703 "slarrd.f"
		}
#line 704 "slarrd.f"
		i__3 = iwork[j + in] + iwoff;
#line 704 "slarrd.f"
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
#line 706 "slarrd.f"
		    w[je] = tmp1;
#line 707 "slarrd.f"
		    werr[je] = tmp2;
#line 708 "slarrd.f"
		    indexw[je] = je - iwoff;
#line 709 "slarrd.f"
		    iblock[je] = ib;
#line 710 "slarrd.f"
/* L50: */
#line 710 "slarrd.f"
		}
#line 711 "slarrd.f"
/* L60: */
#line 711 "slarrd.f"
	    }

#line 713 "slarrd.f"
	    *m += im;
#line 714 "slarrd.f"
	}
#line 715 "slarrd.f"
L70:
#line 715 "slarrd.f"
	;
#line 715 "slarrd.f"
    }
/*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU */
/*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues. */
#line 719 "slarrd.f"
    if (irange == 3) {
#line 720 "slarrd.f"
	idiscl = *il - 1 - nwl;
#line 721 "slarrd.f"
	idiscu = nwu - *iu;

#line 723 "slarrd.f"
	if (idiscl > 0) {
#line 724 "slarrd.f"
	    im = 0;
#line 725 "slarrd.f"
	    i__1 = *m;
#line 725 "slarrd.f"
	    for (je = 1; je <= i__1; ++je) {
/*              Remove some of the smallest eigenvalues from the left so that */
/*              at the end IDISCL =0. Move all eigenvalues up to the left. */
#line 728 "slarrd.f"
		if (w[je] <= wlu && idiscl > 0) {
#line 729 "slarrd.f"
		    --idiscl;
#line 730 "slarrd.f"
		} else {
#line 731 "slarrd.f"
		    ++im;
#line 732 "slarrd.f"
		    w[im] = w[je];
#line 733 "slarrd.f"
		    werr[im] = werr[je];
#line 734 "slarrd.f"
		    indexw[im] = indexw[je];
#line 735 "slarrd.f"
		    iblock[im] = iblock[je];
#line 736 "slarrd.f"
		}
#line 737 "slarrd.f"
/* L80: */
#line 737 "slarrd.f"
	    }
#line 738 "slarrd.f"
	    *m = im;
#line 739 "slarrd.f"
	}
#line 740 "slarrd.f"
	if (idiscu > 0) {
/*           Remove some of the largest eigenvalues from the right so that */
/*           at the end IDISCU =0. Move all eigenvalues up to the left. */
#line 743 "slarrd.f"
	    im = *m + 1;
#line 744 "slarrd.f"
	    for (je = *m; je >= 1; --je) {
#line 745 "slarrd.f"
		if (w[je] >= wul && idiscu > 0) {
#line 746 "slarrd.f"
		    --idiscu;
#line 747 "slarrd.f"
		} else {
#line 748 "slarrd.f"
		    --im;
#line 749 "slarrd.f"
		    w[im] = w[je];
#line 750 "slarrd.f"
		    werr[im] = werr[je];
#line 751 "slarrd.f"
		    indexw[im] = indexw[je];
#line 752 "slarrd.f"
		    iblock[im] = iblock[je];
#line 753 "slarrd.f"
		}
#line 754 "slarrd.f"
/* L81: */
#line 754 "slarrd.f"
	    }
#line 755 "slarrd.f"
	    jee = 0;
#line 756 "slarrd.f"
	    i__1 = *m;
#line 756 "slarrd.f"
	    for (je = im; je <= i__1; ++je) {
#line 757 "slarrd.f"
		++jee;
#line 758 "slarrd.f"
		w[jee] = w[je];
#line 759 "slarrd.f"
		werr[jee] = werr[je];
#line 760 "slarrd.f"
		indexw[jee] = indexw[je];
#line 761 "slarrd.f"
		iblock[jee] = iblock[je];
#line 762 "slarrd.f"
/* L82: */
#line 762 "slarrd.f"
	    }
#line 763 "slarrd.f"
	    *m = *m - im + 1;
#line 764 "slarrd.f"
	}
#line 766 "slarrd.f"
	if (idiscl > 0 || idiscu > 0) {
/*           Code to deal with effects of bad arithmetic. (If N(w) is */
/*           monotone non-decreasing, this should never happen.) */
/*           Some low eigenvalues to be discarded are not in (WL,WLU], */
/*           or high eigenvalues to be discarded are not in (WUL,WU] */
/*           so just kill off the smallest IDISCL/largest IDISCU */
/*           eigenvalues, by marking the corresponding IBLOCK = 0 */
#line 773 "slarrd.f"
	    if (idiscl > 0) {
#line 774 "slarrd.f"
		wkill = *wu;
#line 775 "slarrd.f"
		i__1 = idiscl;
#line 775 "slarrd.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 776 "slarrd.f"
		    iw = 0;
#line 777 "slarrd.f"
		    i__2 = *m;
#line 777 "slarrd.f"
		    for (je = 1; je <= i__2; ++je) {
#line 778 "slarrd.f"
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
#line 780 "slarrd.f"
			    iw = je;
#line 781 "slarrd.f"
			    wkill = w[je];
#line 782 "slarrd.f"
			}
#line 783 "slarrd.f"
/* L90: */
#line 783 "slarrd.f"
		    }
#line 784 "slarrd.f"
		    iblock[iw] = 0;
#line 785 "slarrd.f"
/* L100: */
#line 785 "slarrd.f"
		}
#line 786 "slarrd.f"
	    }
#line 787 "slarrd.f"
	    if (idiscu > 0) {
#line 788 "slarrd.f"
		wkill = *wl;
#line 789 "slarrd.f"
		i__1 = idiscu;
#line 789 "slarrd.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 790 "slarrd.f"
		    iw = 0;
#line 791 "slarrd.f"
		    i__2 = *m;
#line 791 "slarrd.f"
		    for (je = 1; je <= i__2; ++je) {
#line 792 "slarrd.f"
			if (iblock[je] != 0 && (w[je] >= wkill || iw == 0)) {
#line 794 "slarrd.f"
			    iw = je;
#line 795 "slarrd.f"
			    wkill = w[je];
#line 796 "slarrd.f"
			}
#line 797 "slarrd.f"
/* L110: */
#line 797 "slarrd.f"
		    }
#line 798 "slarrd.f"
		    iblock[iw] = 0;
#line 799 "slarrd.f"
/* L120: */
#line 799 "slarrd.f"
		}
#line 800 "slarrd.f"
	    }
/*           Now erase all eigenvalues with IBLOCK set to zero */
#line 802 "slarrd.f"
	    im = 0;
#line 803 "slarrd.f"
	    i__1 = *m;
#line 803 "slarrd.f"
	    for (je = 1; je <= i__1; ++je) {
#line 804 "slarrd.f"
		if (iblock[je] != 0) {
#line 805 "slarrd.f"
		    ++im;
#line 806 "slarrd.f"
		    w[im] = w[je];
#line 807 "slarrd.f"
		    werr[im] = werr[je];
#line 808 "slarrd.f"
		    indexw[im] = indexw[je];
#line 809 "slarrd.f"
		    iblock[im] = iblock[je];
#line 810 "slarrd.f"
		}
#line 811 "slarrd.f"
/* L130: */
#line 811 "slarrd.f"
	    }
#line 812 "slarrd.f"
	    *m = im;
#line 813 "slarrd.f"
	}
#line 814 "slarrd.f"
	if (idiscl < 0 || idiscu < 0) {
#line 815 "slarrd.f"
	    toofew = TRUE_;
#line 816 "slarrd.f"
	}
#line 817 "slarrd.f"
    }

#line 819 "slarrd.f"
    if (irange == 1 && *m != *n || irange == 3 && *m != *iu - *il + 1) {
#line 821 "slarrd.f"
	toofew = TRUE_;
#line 822 "slarrd.f"
    }
/*     If ORDER='B', do nothing the eigenvalues are already sorted by */
/*        block. */
/*     If ORDER='E', sort the eigenvalues from smallest to largest */
#line 828 "slarrd.f"
    if (lsame_(order, "E", (ftnlen)1, (ftnlen)1) && *nsplit > 1) {
#line 829 "slarrd.f"
	i__1 = *m - 1;
#line 829 "slarrd.f"
	for (je = 1; je <= i__1; ++je) {
#line 830 "slarrd.f"
	    ie = 0;
#line 831 "slarrd.f"
	    tmp1 = w[je];
#line 832 "slarrd.f"
	    i__2 = *m;
#line 832 "slarrd.f"
	    for (j = je + 1; j <= i__2; ++j) {
#line 833 "slarrd.f"
		if (w[j] < tmp1) {
#line 834 "slarrd.f"
		    ie = j;
#line 835 "slarrd.f"
		    tmp1 = w[j];
#line 836 "slarrd.f"
		}
#line 837 "slarrd.f"
/* L140: */
#line 837 "slarrd.f"
	    }
#line 838 "slarrd.f"
	    if (ie != 0) {
#line 839 "slarrd.f"
		tmp2 = werr[ie];
#line 840 "slarrd.f"
		itmp1 = iblock[ie];
#line 841 "slarrd.f"
		itmp2 = indexw[ie];
#line 842 "slarrd.f"
		w[ie] = w[je];
#line 843 "slarrd.f"
		werr[ie] = werr[je];
#line 844 "slarrd.f"
		iblock[ie] = iblock[je];
#line 845 "slarrd.f"
		indexw[ie] = indexw[je];
#line 846 "slarrd.f"
		w[je] = tmp1;
#line 847 "slarrd.f"
		werr[je] = tmp2;
#line 848 "slarrd.f"
		iblock[je] = itmp1;
#line 849 "slarrd.f"
		indexw[je] = itmp2;
#line 850 "slarrd.f"
	    }
#line 851 "slarrd.f"
/* L150: */
#line 851 "slarrd.f"
	}
#line 852 "slarrd.f"
    }

#line 854 "slarrd.f"
    *info = 0;
#line 855 "slarrd.f"
    if (ncnvrg) {
#line 855 "slarrd.f"
	++(*info);
#line 855 "slarrd.f"
    }
#line 857 "slarrd.f"
    if (toofew) {
#line 857 "slarrd.f"
	*info += 2;
#line 857 "slarrd.f"
    }
#line 859 "slarrd.f"
    return 0;

/*     End of SLARRD */

} /* slarrd_ */

