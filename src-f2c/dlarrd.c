#line 1 "dlarrd.f"
/* dlarrd.f -- translated by f2c (version 20100827).
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

#line 1 "dlarrd.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__0 = 0;

/* > \brief \b DLARRD computes the eigenvalues of a symmetric tridiagonal matrix to suitable accuracy. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRD( RANGE, ORDER, N, VL, VU, IL, IU, GERS, */
/*                           RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT, */
/*                           M, W, WERR, WL, WU, IBLOCK, INDEXW, */
/*                           WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          ORDER, RANGE */
/*       INTEGER            IL, INFO, IU, M, N, NSPLIT */
/*       DOUBLE PRECISION    PIVMIN, RELTOL, VL, VU, WL, WU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), INDEXW( * ), */
/*      $                   ISPLIT( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), E2( * ), */
/*      $                   GERS( * ), W( * ), WERR( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARRD computes the eigenvalues of a symmetric tridiagonal */
/* > matrix T to suitable accuracy. This is an auxiliary code to be */
/* > called from DSTEMR. */
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
/* >          VL is DOUBLE PRECISION */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues.  Eigenvalues less than or equal */
/* >          to VL, or greater than VU, will not be returned.  VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
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
/* >          GERS is DOUBLE PRECISION array, dimension (2*N) */
/* >          The N Gerschgorin intervals (the i-th Gerschgorin interval */
/* >          is (GERS(2*i-1), GERS(2*i)). */
/* > \endverbatim */
/* > */
/* > \param[in] RELTOL */
/* > \verbatim */
/* >          RELTOL is DOUBLE PRECISION */
/* >          The minimum relative width of an interval.  When an interval */
/* >          is narrower than RELTOL times the larger (in */
/* >          magnitude) endpoint, then it is considered to be */
/* >          sufficiently small, i.e., converged.  Note: this should */
/* >          always be at least radix*machine epsilon. */
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
/* > \param[in] E2 */
/* > \verbatim */
/* >          E2 is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) squared off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          On exit, the first M elements of W will contain the */
/* >          eigenvalue approximations. DLARRD computes an interval */
/* >          I_j = (a_j, b_j] that includes eigenvalue j. The eigenvalue */
/* >          approximation is given as the interval midpoint */
/* >          W(j)= ( a_j + b_j)/2. The corresponding error is bounded by */
/* >          WERR(j) = abs( a_j - b_j)/2 */
/* > \endverbatim */
/* > */
/* > \param[out] WERR */
/* > \verbatim */
/* >          WERR is DOUBLE PRECISION array, dimension (N) */
/* >          The error bound on the corresponding eigenvalue approximation */
/* >          in W. */
/* > \endverbatim */
/* > */
/* > \param[out] WL */
/* > \verbatim */
/* >          WL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] WU */
/* > \verbatim */
/* >          WU is DOUBLE PRECISION */
/* >          The interval (WL, WU] contains all the wanted eigenvalues. */
/* >          If RANGE='V', then WL=VL and WU=VU. */
/* >          If RANGE='A', then WL and WU are the global Gerschgorin bounds */
/* >                        on the spectrum. */
/* >          If RANGE='I', then WL and WU are computed by DLAEBZ from the */
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
/* >          belongs.  (DLARRD may use the remaining N-M elements as */
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
/* >  FUDGE   DOUBLE PRECISION, default = 2 */
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
/* Subroutine */ int dlarrd_(char *range, char *order, integer *n, doublereal 
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
    extern doublereal dlamch_(char *, ftnlen);
    static integer ibegin;
    extern /* Subroutine */ int dlaebz_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    static integer irange, idiscl, idumma[1];
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer idiscu;
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

#line 386 "dlarrd.f"
    /* Parameter adjustments */
#line 386 "dlarrd.f"
    --iwork;
#line 386 "dlarrd.f"
    --work;
#line 386 "dlarrd.f"
    --indexw;
#line 386 "dlarrd.f"
    --iblock;
#line 386 "dlarrd.f"
    --werr;
#line 386 "dlarrd.f"
    --w;
#line 386 "dlarrd.f"
    --isplit;
#line 386 "dlarrd.f"
    --e2;
#line 386 "dlarrd.f"
    --e;
#line 386 "dlarrd.f"
    --d__;
#line 386 "dlarrd.f"
    --gers;
#line 386 "dlarrd.f"

#line 386 "dlarrd.f"
    /* Function Body */
#line 386 "dlarrd.f"
    *info = 0;

/*     Decode RANGE */

#line 390 "dlarrd.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 391 "dlarrd.f"
	irange = 1;
#line 392 "dlarrd.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 393 "dlarrd.f"
	irange = 2;
#line 394 "dlarrd.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 395 "dlarrd.f"
	irange = 3;
#line 396 "dlarrd.f"
    } else {
#line 397 "dlarrd.f"
	irange = 0;
#line 398 "dlarrd.f"
    }

/*     Check for Errors */

#line 402 "dlarrd.f"
    if (irange <= 0) {
#line 403 "dlarrd.f"
	*info = -1;
#line 404 "dlarrd.f"
    } else if (! (lsame_(order, "B", (ftnlen)1, (ftnlen)1) || lsame_(order, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 405 "dlarrd.f"
	*info = -2;
#line 406 "dlarrd.f"
    } else if (*n < 0) {
#line 407 "dlarrd.f"
	*info = -3;
#line 408 "dlarrd.f"
    } else if (irange == 2) {
#line 409 "dlarrd.f"
	if (*vl >= *vu) {
#line 409 "dlarrd.f"
	    *info = -5;
#line 409 "dlarrd.f"
	}
#line 411 "dlarrd.f"
    } else if (irange == 3 && (*il < 1 || *il > max(1,*n))) {
#line 413 "dlarrd.f"
	*info = -6;
#line 414 "dlarrd.f"
    } else if (irange == 3 && (*iu < min(*n,*il) || *iu > *n)) {
#line 416 "dlarrd.f"
	*info = -7;
#line 417 "dlarrd.f"
    }

#line 419 "dlarrd.f"
    if (*info != 0) {
#line 420 "dlarrd.f"
	return 0;
#line 421 "dlarrd.f"
    }
/*     Initialize error flags */
#line 424 "dlarrd.f"
    *info = 0;
#line 425 "dlarrd.f"
    ncnvrg = FALSE_;
#line 426 "dlarrd.f"
    toofew = FALSE_;
/*     Quick return if possible */
#line 429 "dlarrd.f"
    *m = 0;
#line 430 "dlarrd.f"
    if (*n == 0) {
#line 430 "dlarrd.f"
	return 0;
#line 430 "dlarrd.f"
    }
/*     Simplification: */
#line 433 "dlarrd.f"
    if (irange == 3 && *il == 1 && *iu == *n) {
#line 433 "dlarrd.f"
	irange = 1;
#line 433 "dlarrd.f"
    }
/*     Get machine constants */
#line 436 "dlarrd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 437 "dlarrd.f"
    uflow = dlamch_("U", (ftnlen)1);
/*     Special Case when N=1 */
/*     Treat case of 1x1 matrix for quick return */
#line 442 "dlarrd.f"
    if (*n == 1) {
#line 443 "dlarrd.f"
	if (irange == 1 || irange == 2 && d__[1] > *vl && d__[1] <= *vu || 
		irange == 3 && *il == 1 && *iu == 1) {
#line 446 "dlarrd.f"
	    *m = 1;
#line 447 "dlarrd.f"
	    w[1] = d__[1];
/*           The computation error of the eigenvalue is zero */
#line 449 "dlarrd.f"
	    werr[1] = 0.;
#line 450 "dlarrd.f"
	    iblock[1] = 1;
#line 451 "dlarrd.f"
	    indexw[1] = 1;
#line 452 "dlarrd.f"
	}
#line 453 "dlarrd.f"
	return 0;
#line 454 "dlarrd.f"
    }
/*     NB is the minimum vector length for vector bisection, or 0 */
/*     if only scalar is to be done. */
#line 458 "dlarrd.f"
    nb = ilaenv_(&c__1, "DSTEBZ", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 459 "dlarrd.f"
    if (nb <= 1) {
#line 459 "dlarrd.f"
	nb = 0;
#line 459 "dlarrd.f"
    }
/*     Find global spectral radius */
#line 462 "dlarrd.f"
    gl = d__[1];
#line 463 "dlarrd.f"
    gu = d__[1];
#line 464 "dlarrd.f"
    i__1 = *n;
#line 464 "dlarrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 465 "dlarrd.f"
	d__1 = gl, d__2 = gers[(i__ << 1) - 1];
#line 465 "dlarrd.f"
	gl = min(d__1,d__2);
/* Computing MAX */
#line 466 "dlarrd.f"
	d__1 = gu, d__2 = gers[i__ * 2];
#line 466 "dlarrd.f"
	gu = max(d__1,d__2);
#line 467 "dlarrd.f"
/* L5: */
#line 467 "dlarrd.f"
    }
/*     Compute global Gerschgorin bounds and spectral diameter */
/* Computing MAX */
#line 469 "dlarrd.f"
    d__1 = abs(gl), d__2 = abs(gu);
#line 469 "dlarrd.f"
    tnorm = max(d__1,d__2);
#line 470 "dlarrd.f"
    gl = gl - tnorm * 2. * eps * *n - *pivmin * 4.;
#line 471 "dlarrd.f"
    gu = gu + tnorm * 2. * eps * *n + *pivmin * 4.;
/*     [JAN/28/2009] remove the line below since SPDIAM variable not use */
/*     SPDIAM = GU - GL */
/*     Input arguments for DLAEBZ: */
/*     The relative tolerance.  An interval (a,b] lies within */
/*     "relative tolerance" if  b-a < RELTOL*max(|a|,|b|), */
#line 477 "dlarrd.f"
    rtoli = *reltol;
/*     Set the absolute tolerance for interval convergence to zero to force */
/*     interval convergence based on relative size of the interval. */
/*     This is dangerous because intervals might not converge when RELTOL is */
/*     small. But at least a very small number should be selected so that for */
/*     strongly graded matrices, the code can get relatively accurate */
/*     eigenvalues. */
#line 484 "dlarrd.f"
    atoli = uflow * 4. + *pivmin * 4.;
#line 486 "dlarrd.f"
    if (irange == 3) {
/*        RANGE='I': Compute an interval containing eigenvalues */
/*        IL through IU. The initial interval [GL,GU] from the global */
/*        Gerschgorin bounds GL and GU is refined by DLAEBZ. */
#line 491 "dlarrd.f"
	itmax = (integer) ((log(tnorm + *pivmin) - log(*pivmin)) / log(2.)) + 
		2;
#line 493 "dlarrd.f"
	work[*n + 1] = gl;
#line 494 "dlarrd.f"
	work[*n + 2] = gl;
#line 495 "dlarrd.f"
	work[*n + 3] = gu;
#line 496 "dlarrd.f"
	work[*n + 4] = gu;
#line 497 "dlarrd.f"
	work[*n + 5] = gl;
#line 498 "dlarrd.f"
	work[*n + 6] = gu;
#line 499 "dlarrd.f"
	iwork[1] = -1;
#line 500 "dlarrd.f"
	iwork[2] = -1;
#line 501 "dlarrd.f"
	iwork[3] = *n + 1;
#line 502 "dlarrd.f"
	iwork[4] = *n + 1;
#line 503 "dlarrd.f"
	iwork[5] = *il - 1;
#line 504 "dlarrd.f"
	iwork[6] = *iu;

#line 506 "dlarrd.f"
	dlaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, pivmin, &
		d__[1], &e[1], &e2[1], &iwork[5], &work[*n + 1], &work[*n + 5]
		, &iout, &iwork[1], &w[1], &iblock[1], &iinfo);
#line 509 "dlarrd.f"
	if (iinfo != 0) {
#line 510 "dlarrd.f"
	    *info = iinfo;
#line 511 "dlarrd.f"
	    return 0;
#line 512 "dlarrd.f"
	}
/*        On exit, output intervals may not be ordered by ascending negcount */
#line 514 "dlarrd.f"
	if (iwork[6] == *iu) {
#line 515 "dlarrd.f"
	    *wl = work[*n + 1];
#line 516 "dlarrd.f"
	    wlu = work[*n + 3];
#line 517 "dlarrd.f"
	    nwl = iwork[1];
#line 518 "dlarrd.f"
	    *wu = work[*n + 4];
#line 519 "dlarrd.f"
	    wul = work[*n + 2];
#line 520 "dlarrd.f"
	    nwu = iwork[4];
#line 521 "dlarrd.f"
	} else {
#line 522 "dlarrd.f"
	    *wl = work[*n + 2];
#line 523 "dlarrd.f"
	    wlu = work[*n + 4];
#line 524 "dlarrd.f"
	    nwl = iwork[2];
#line 525 "dlarrd.f"
	    *wu = work[*n + 3];
#line 526 "dlarrd.f"
	    wul = work[*n + 1];
#line 527 "dlarrd.f"
	    nwu = iwork[3];
#line 528 "dlarrd.f"
	}
/*        On exit, the interval [WL, WLU] contains a value with negcount NWL, */
/*        and [WUL, WU] contains a value with negcount NWU. */
#line 531 "dlarrd.f"
	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
#line 532 "dlarrd.f"
	    *info = 4;
#line 533 "dlarrd.f"
	    return 0;
#line 534 "dlarrd.f"
	}
#line 536 "dlarrd.f"
    } else if (irange == 2) {
#line 537 "dlarrd.f"
	*wl = *vl;
#line 538 "dlarrd.f"
	*wu = *vu;
#line 540 "dlarrd.f"
    } else if (irange == 1) {
#line 541 "dlarrd.f"
	*wl = gl;
#line 542 "dlarrd.f"
	*wu = gu;
#line 543 "dlarrd.f"
    }
/*     Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU. */
/*     NWL accumulates the number of eigenvalues .le. WL, */
/*     NWU accumulates the number of eigenvalues .le. WU */
#line 550 "dlarrd.f"
    *m = 0;
#line 551 "dlarrd.f"
    iend = 0;
#line 552 "dlarrd.f"
    *info = 0;
#line 553 "dlarrd.f"
    nwl = 0;
#line 554 "dlarrd.f"
    nwu = 0;

#line 556 "dlarrd.f"
    i__1 = *nsplit;
#line 556 "dlarrd.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 557 "dlarrd.f"
	ioff = iend;
#line 558 "dlarrd.f"
	ibegin = ioff + 1;
#line 559 "dlarrd.f"
	iend = isplit[jblk];
#line 560 "dlarrd.f"
	in = iend - ioff;

#line 562 "dlarrd.f"
	if (in == 1) {
/*           1x1 block */
#line 564 "dlarrd.f"
	    if (*wl >= d__[ibegin] - *pivmin) {
#line 564 "dlarrd.f"
		++nwl;
#line 564 "dlarrd.f"
	    }
#line 566 "dlarrd.f"
	    if (*wu >= d__[ibegin] - *pivmin) {
#line 566 "dlarrd.f"
		++nwu;
#line 566 "dlarrd.f"
	    }
#line 568 "dlarrd.f"
	    if (irange == 1 || *wl < d__[ibegin] - *pivmin && *wu >= d__[
		    ibegin] - *pivmin) {
#line 571 "dlarrd.f"
		++(*m);
#line 572 "dlarrd.f"
		w[*m] = d__[ibegin];
#line 573 "dlarrd.f"
		werr[*m] = 0.;
/*              The gap for a single block doesn't matter for the later */
/*              algorithm and is assigned an arbitrary large value */
#line 576 "dlarrd.f"
		iblock[*m] = jblk;
#line 577 "dlarrd.f"
		indexw[*m] = 1;
#line 578 "dlarrd.f"
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
#line 626 "dlarrd.f"
	} else {
/*           General Case - block of size IN >= 2 */
/*           Compute local Gerschgorin interval and use it as the initial */
/*           interval for DLAEBZ */
#line 630 "dlarrd.f"
	    gu = d__[ibegin];
#line 631 "dlarrd.f"
	    gl = d__[ibegin];
#line 632 "dlarrd.f"
	    tmp1 = 0.;
#line 634 "dlarrd.f"
	    i__2 = iend;
#line 634 "dlarrd.f"
	    for (j = ibegin; j <= i__2; ++j) {
/* Computing MIN */
#line 635 "dlarrd.f"
		d__1 = gl, d__2 = gers[(j << 1) - 1];
#line 635 "dlarrd.f"
		gl = min(d__1,d__2);
/* Computing MAX */
#line 636 "dlarrd.f"
		d__1 = gu, d__2 = gers[j * 2];
#line 636 "dlarrd.f"
		gu = max(d__1,d__2);
#line 637 "dlarrd.f"
/* L40: */
#line 637 "dlarrd.f"
	    }
/*           [JAN/28/2009] */
/*           change SPDIAM by TNORM in lines 2 and 3 thereafter */
/*           line 1: remove computation of SPDIAM (not useful anymore) */
/*           SPDIAM = GU - GL */
/*           GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN */
/*           GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN */
#line 644 "dlarrd.f"
	    gl = gl - tnorm * 2. * eps * in - *pivmin * 2.;
#line 645 "dlarrd.f"
	    gu = gu + tnorm * 2. * eps * in + *pivmin * 2.;

#line 647 "dlarrd.f"
	    if (irange > 1) {
#line 648 "dlarrd.f"
		if (gu < *wl) {
/*                 the local block contains none of the wanted eigenvalues */
#line 650 "dlarrd.f"
		    nwl += in;
#line 651 "dlarrd.f"
		    nwu += in;
#line 652 "dlarrd.f"
		    goto L70;
#line 653 "dlarrd.f"
		}
/*              refine search interval if possible, only range (WL,WU] matters */
#line 655 "dlarrd.f"
		gl = max(gl,*wl);
#line 656 "dlarrd.f"
		gu = min(gu,*wu);
#line 657 "dlarrd.f"
		if (gl >= gu) {
#line 657 "dlarrd.f"
		    goto L70;
#line 657 "dlarrd.f"
		}
#line 659 "dlarrd.f"
	    }
/*           Find negcount of initial interval boundaries GL and GU */
#line 662 "dlarrd.f"
	    work[*n + 1] = gl;
#line 663 "dlarrd.f"
	    work[*n + in + 1] = gu;
#line 664 "dlarrd.f"
	    dlaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, 
		    pivmin, &d__[ibegin], &e[ibegin], &e2[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);
#line 668 "dlarrd.f"
	    if (iinfo != 0) {
#line 669 "dlarrd.f"
		*info = iinfo;
#line 670 "dlarrd.f"
		return 0;
#line 671 "dlarrd.f"
	    }

#line 673 "dlarrd.f"
	    nwl += iwork[1];
#line 674 "dlarrd.f"
	    nwu += iwork[in + 1];
#line 675 "dlarrd.f"
	    iwoff = *m - iwork[1];
/*           Compute Eigenvalues */
#line 678 "dlarrd.f"
	    itmax = (integer) ((log(gu - gl + *pivmin) - log(*pivmin)) / log(
		    2.)) + 2;
#line 680 "dlarrd.f"
	    dlaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, 
		    pivmin, &d__[ibegin], &e[ibegin], &e2[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);
#line 684 "dlarrd.f"
	    if (iinfo != 0) {
#line 685 "dlarrd.f"
		*info = iinfo;
#line 686 "dlarrd.f"
		return 0;
#line 687 "dlarrd.f"
	    }

/*           Copy eigenvalues into W and IBLOCK */
/*           Use -JBLK for block number for unconverged eigenvalues. */
/*           Loop over the number of output intervals from DLAEBZ */
#line 692 "dlarrd.f"
	    i__2 = iout;
#line 692 "dlarrd.f"
	    for (j = 1; j <= i__2; ++j) {
/*              eigenvalue approximation is middle point of interval */
#line 694 "dlarrd.f"
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;
/*              semi length of error interval */
#line 696 "dlarrd.f"
		tmp2 = (d__1 = work[j + *n] - work[j + in + *n], abs(d__1)) * 
			.5;
#line 697 "dlarrd.f"
		if (j > iout - iinfo) {
/*                 Flag non-convergence. */
#line 699 "dlarrd.f"
		    ncnvrg = TRUE_;
#line 700 "dlarrd.f"
		    ib = -jblk;
#line 701 "dlarrd.f"
		} else {
#line 702 "dlarrd.f"
		    ib = jblk;
#line 703 "dlarrd.f"
		}
#line 704 "dlarrd.f"
		i__3 = iwork[j + in] + iwoff;
#line 704 "dlarrd.f"
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
#line 706 "dlarrd.f"
		    w[je] = tmp1;
#line 707 "dlarrd.f"
		    werr[je] = tmp2;
#line 708 "dlarrd.f"
		    indexw[je] = je - iwoff;
#line 709 "dlarrd.f"
		    iblock[je] = ib;
#line 710 "dlarrd.f"
/* L50: */
#line 710 "dlarrd.f"
		}
#line 711 "dlarrd.f"
/* L60: */
#line 711 "dlarrd.f"
	    }

#line 713 "dlarrd.f"
	    *m += im;
#line 714 "dlarrd.f"
	}
#line 715 "dlarrd.f"
L70:
#line 715 "dlarrd.f"
	;
#line 715 "dlarrd.f"
    }
/*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU */
/*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues. */
#line 719 "dlarrd.f"
    if (irange == 3) {
#line 720 "dlarrd.f"
	idiscl = *il - 1 - nwl;
#line 721 "dlarrd.f"
	idiscu = nwu - *iu;

#line 723 "dlarrd.f"
	if (idiscl > 0) {
#line 724 "dlarrd.f"
	    im = 0;
#line 725 "dlarrd.f"
	    i__1 = *m;
#line 725 "dlarrd.f"
	    for (je = 1; je <= i__1; ++je) {
/*              Remove some of the smallest eigenvalues from the left so that */
/*              at the end IDISCL =0. Move all eigenvalues up to the left. */
#line 728 "dlarrd.f"
		if (w[je] <= wlu && idiscl > 0) {
#line 729 "dlarrd.f"
		    --idiscl;
#line 730 "dlarrd.f"
		} else {
#line 731 "dlarrd.f"
		    ++im;
#line 732 "dlarrd.f"
		    w[im] = w[je];
#line 733 "dlarrd.f"
		    werr[im] = werr[je];
#line 734 "dlarrd.f"
		    indexw[im] = indexw[je];
#line 735 "dlarrd.f"
		    iblock[im] = iblock[je];
#line 736 "dlarrd.f"
		}
#line 737 "dlarrd.f"
/* L80: */
#line 737 "dlarrd.f"
	    }
#line 738 "dlarrd.f"
	    *m = im;
#line 739 "dlarrd.f"
	}
#line 740 "dlarrd.f"
	if (idiscu > 0) {
/*           Remove some of the largest eigenvalues from the right so that */
/*           at the end IDISCU =0. Move all eigenvalues up to the left. */
#line 743 "dlarrd.f"
	    im = *m + 1;
#line 744 "dlarrd.f"
	    for (je = *m; je >= 1; --je) {
#line 745 "dlarrd.f"
		if (w[je] >= wul && idiscu > 0) {
#line 746 "dlarrd.f"
		    --idiscu;
#line 747 "dlarrd.f"
		} else {
#line 748 "dlarrd.f"
		    --im;
#line 749 "dlarrd.f"
		    w[im] = w[je];
#line 750 "dlarrd.f"
		    werr[im] = werr[je];
#line 751 "dlarrd.f"
		    indexw[im] = indexw[je];
#line 752 "dlarrd.f"
		    iblock[im] = iblock[je];
#line 753 "dlarrd.f"
		}
#line 754 "dlarrd.f"
/* L81: */
#line 754 "dlarrd.f"
	    }
#line 755 "dlarrd.f"
	    jee = 0;
#line 756 "dlarrd.f"
	    i__1 = *m;
#line 756 "dlarrd.f"
	    for (je = im; je <= i__1; ++je) {
#line 757 "dlarrd.f"
		++jee;
#line 758 "dlarrd.f"
		w[jee] = w[je];
#line 759 "dlarrd.f"
		werr[jee] = werr[je];
#line 760 "dlarrd.f"
		indexw[jee] = indexw[je];
#line 761 "dlarrd.f"
		iblock[jee] = iblock[je];
#line 762 "dlarrd.f"
/* L82: */
#line 762 "dlarrd.f"
	    }
#line 763 "dlarrd.f"
	    *m = *m - im + 1;
#line 764 "dlarrd.f"
	}
#line 766 "dlarrd.f"
	if (idiscl > 0 || idiscu > 0) {
/*           Code to deal with effects of bad arithmetic. (If N(w) is */
/*           monotone non-decreasing, this should never happen.) */
/*           Some low eigenvalues to be discarded are not in (WL,WLU], */
/*           or high eigenvalues to be discarded are not in (WUL,WU] */
/*           so just kill off the smallest IDISCL/largest IDISCU */
/*           eigenvalues, by marking the corresponding IBLOCK = 0 */
#line 773 "dlarrd.f"
	    if (idiscl > 0) {
#line 774 "dlarrd.f"
		wkill = *wu;
#line 775 "dlarrd.f"
		i__1 = idiscl;
#line 775 "dlarrd.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 776 "dlarrd.f"
		    iw = 0;
#line 777 "dlarrd.f"
		    i__2 = *m;
#line 777 "dlarrd.f"
		    for (je = 1; je <= i__2; ++je) {
#line 778 "dlarrd.f"
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
#line 780 "dlarrd.f"
			    iw = je;
#line 781 "dlarrd.f"
			    wkill = w[je];
#line 782 "dlarrd.f"
			}
#line 783 "dlarrd.f"
/* L90: */
#line 783 "dlarrd.f"
		    }
#line 784 "dlarrd.f"
		    iblock[iw] = 0;
#line 785 "dlarrd.f"
/* L100: */
#line 785 "dlarrd.f"
		}
#line 786 "dlarrd.f"
	    }
#line 787 "dlarrd.f"
	    if (idiscu > 0) {
#line 788 "dlarrd.f"
		wkill = *wl;
#line 789 "dlarrd.f"
		i__1 = idiscu;
#line 789 "dlarrd.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 790 "dlarrd.f"
		    iw = 0;
#line 791 "dlarrd.f"
		    i__2 = *m;
#line 791 "dlarrd.f"
		    for (je = 1; je <= i__2; ++je) {
#line 792 "dlarrd.f"
			if (iblock[je] != 0 && (w[je] >= wkill || iw == 0)) {
#line 794 "dlarrd.f"
			    iw = je;
#line 795 "dlarrd.f"
			    wkill = w[je];
#line 796 "dlarrd.f"
			}
#line 797 "dlarrd.f"
/* L110: */
#line 797 "dlarrd.f"
		    }
#line 798 "dlarrd.f"
		    iblock[iw] = 0;
#line 799 "dlarrd.f"
/* L120: */
#line 799 "dlarrd.f"
		}
#line 800 "dlarrd.f"
	    }
/*           Now erase all eigenvalues with IBLOCK set to zero */
#line 802 "dlarrd.f"
	    im = 0;
#line 803 "dlarrd.f"
	    i__1 = *m;
#line 803 "dlarrd.f"
	    for (je = 1; je <= i__1; ++je) {
#line 804 "dlarrd.f"
		if (iblock[je] != 0) {
#line 805 "dlarrd.f"
		    ++im;
#line 806 "dlarrd.f"
		    w[im] = w[je];
#line 807 "dlarrd.f"
		    werr[im] = werr[je];
#line 808 "dlarrd.f"
		    indexw[im] = indexw[je];
#line 809 "dlarrd.f"
		    iblock[im] = iblock[je];
#line 810 "dlarrd.f"
		}
#line 811 "dlarrd.f"
/* L130: */
#line 811 "dlarrd.f"
	    }
#line 812 "dlarrd.f"
	    *m = im;
#line 813 "dlarrd.f"
	}
#line 814 "dlarrd.f"
	if (idiscl < 0 || idiscu < 0) {
#line 815 "dlarrd.f"
	    toofew = TRUE_;
#line 816 "dlarrd.f"
	}
#line 817 "dlarrd.f"
    }

#line 819 "dlarrd.f"
    if (irange == 1 && *m != *n || irange == 3 && *m != *iu - *il + 1) {
#line 821 "dlarrd.f"
	toofew = TRUE_;
#line 822 "dlarrd.f"
    }
/*     If ORDER='B', do nothing the eigenvalues are already sorted by */
/*        block. */
/*     If ORDER='E', sort the eigenvalues from smallest to largest */
#line 828 "dlarrd.f"
    if (lsame_(order, "E", (ftnlen)1, (ftnlen)1) && *nsplit > 1) {
#line 829 "dlarrd.f"
	i__1 = *m - 1;
#line 829 "dlarrd.f"
	for (je = 1; je <= i__1; ++je) {
#line 830 "dlarrd.f"
	    ie = 0;
#line 831 "dlarrd.f"
	    tmp1 = w[je];
#line 832 "dlarrd.f"
	    i__2 = *m;
#line 832 "dlarrd.f"
	    for (j = je + 1; j <= i__2; ++j) {
#line 833 "dlarrd.f"
		if (w[j] < tmp1) {
#line 834 "dlarrd.f"
		    ie = j;
#line 835 "dlarrd.f"
		    tmp1 = w[j];
#line 836 "dlarrd.f"
		}
#line 837 "dlarrd.f"
/* L140: */
#line 837 "dlarrd.f"
	    }
#line 838 "dlarrd.f"
	    if (ie != 0) {
#line 839 "dlarrd.f"
		tmp2 = werr[ie];
#line 840 "dlarrd.f"
		itmp1 = iblock[ie];
#line 841 "dlarrd.f"
		itmp2 = indexw[ie];
#line 842 "dlarrd.f"
		w[ie] = w[je];
#line 843 "dlarrd.f"
		werr[ie] = werr[je];
#line 844 "dlarrd.f"
		iblock[ie] = iblock[je];
#line 845 "dlarrd.f"
		indexw[ie] = indexw[je];
#line 846 "dlarrd.f"
		w[je] = tmp1;
#line 847 "dlarrd.f"
		werr[je] = tmp2;
#line 848 "dlarrd.f"
		iblock[je] = itmp1;
#line 849 "dlarrd.f"
		indexw[je] = itmp2;
#line 850 "dlarrd.f"
	    }
#line 851 "dlarrd.f"
/* L150: */
#line 851 "dlarrd.f"
	}
#line 852 "dlarrd.f"
    }

#line 854 "dlarrd.f"
    *info = 0;
#line 855 "dlarrd.f"
    if (ncnvrg) {
#line 855 "dlarrd.f"
	++(*info);
#line 855 "dlarrd.f"
    }
#line 857 "dlarrd.f"
    if (toofew) {
#line 857 "dlarrd.f"
	*info += 2;
#line 857 "dlarrd.f"
    }
#line 859 "dlarrd.f"
    return 0;

/*     End of DLARRD */

} /* dlarrd_ */

