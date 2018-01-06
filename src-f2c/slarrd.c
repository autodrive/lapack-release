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

/*     Quick return if possible */

#line 390 "slarrd.f"
    if (*n <= 0) {
#line 391 "slarrd.f"
	return 0;
#line 392 "slarrd.f"
    }

/*     Decode RANGE */

#line 396 "slarrd.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 397 "slarrd.f"
	irange = 1;
#line 398 "slarrd.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 399 "slarrd.f"
	irange = 2;
#line 400 "slarrd.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 401 "slarrd.f"
	irange = 3;
#line 402 "slarrd.f"
    } else {
#line 403 "slarrd.f"
	irange = 0;
#line 404 "slarrd.f"
    }

/*     Check for Errors */

#line 408 "slarrd.f"
    if (irange <= 0) {
#line 409 "slarrd.f"
	*info = -1;
#line 410 "slarrd.f"
    } else if (! (lsame_(order, "B", (ftnlen)1, (ftnlen)1) || lsame_(order, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 411 "slarrd.f"
	*info = -2;
#line 412 "slarrd.f"
    } else if (*n < 0) {
#line 413 "slarrd.f"
	*info = -3;
#line 414 "slarrd.f"
    } else if (irange == 2) {
#line 415 "slarrd.f"
	if (*vl >= *vu) {
#line 415 "slarrd.f"
	    *info = -5;
#line 415 "slarrd.f"
	}
#line 417 "slarrd.f"
    } else if (irange == 3 && (*il < 1 || *il > max(1,*n))) {
#line 419 "slarrd.f"
	*info = -6;
#line 420 "slarrd.f"
    } else if (irange == 3 && (*iu < min(*n,*il) || *iu > *n)) {
#line 422 "slarrd.f"
	*info = -7;
#line 423 "slarrd.f"
    }

#line 425 "slarrd.f"
    if (*info != 0) {
#line 426 "slarrd.f"
	return 0;
#line 427 "slarrd.f"
    }
/*     Initialize error flags */
#line 430 "slarrd.f"
    *info = 0;
#line 431 "slarrd.f"
    ncnvrg = FALSE_;
#line 432 "slarrd.f"
    toofew = FALSE_;
/*     Quick return if possible */
#line 435 "slarrd.f"
    *m = 0;
#line 436 "slarrd.f"
    if (*n == 0) {
#line 436 "slarrd.f"
	return 0;
#line 436 "slarrd.f"
    }
/*     Simplification: */
#line 439 "slarrd.f"
    if (irange == 3 && *il == 1 && *iu == *n) {
#line 439 "slarrd.f"
	irange = 1;
#line 439 "slarrd.f"
    }
/*     Get machine constants */
#line 442 "slarrd.f"
    eps = slamch_("P", (ftnlen)1);
#line 443 "slarrd.f"
    uflow = slamch_("U", (ftnlen)1);
/*     Special Case when N=1 */
/*     Treat case of 1x1 matrix for quick return */
#line 448 "slarrd.f"
    if (*n == 1) {
#line 449 "slarrd.f"
	if (irange == 1 || irange == 2 && d__[1] > *vl && d__[1] <= *vu || 
		irange == 3 && *il == 1 && *iu == 1) {
#line 452 "slarrd.f"
	    *m = 1;
#line 453 "slarrd.f"
	    w[1] = d__[1];
/*           The computation error of the eigenvalue is zero */
#line 455 "slarrd.f"
	    werr[1] = 0.;
#line 456 "slarrd.f"
	    iblock[1] = 1;
#line 457 "slarrd.f"
	    indexw[1] = 1;
#line 458 "slarrd.f"
	}
#line 459 "slarrd.f"
	return 0;
#line 460 "slarrd.f"
    }
/*     NB is the minimum vector length for vector bisection, or 0 */
/*     if only scalar is to be done. */
#line 464 "slarrd.f"
    nb = ilaenv_(&c__1, "SSTEBZ", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 465 "slarrd.f"
    if (nb <= 1) {
#line 465 "slarrd.f"
	nb = 0;
#line 465 "slarrd.f"
    }
/*     Find global spectral radius */
#line 468 "slarrd.f"
    gl = d__[1];
#line 469 "slarrd.f"
    gu = d__[1];
#line 470 "slarrd.f"
    i__1 = *n;
#line 470 "slarrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 471 "slarrd.f"
	d__1 = gl, d__2 = gers[(i__ << 1) - 1];
#line 471 "slarrd.f"
	gl = min(d__1,d__2);
/* Computing MAX */
#line 472 "slarrd.f"
	d__1 = gu, d__2 = gers[i__ * 2];
#line 472 "slarrd.f"
	gu = max(d__1,d__2);
#line 473 "slarrd.f"
/* L5: */
#line 473 "slarrd.f"
    }
/*     Compute global Gerschgorin bounds and spectral diameter */
/* Computing MAX */
#line 475 "slarrd.f"
    d__1 = abs(gl), d__2 = abs(gu);
#line 475 "slarrd.f"
    tnorm = max(d__1,d__2);
#line 476 "slarrd.f"
    gl = gl - tnorm * 2. * eps * *n - *pivmin * 4.;
#line 477 "slarrd.f"
    gu = gu + tnorm * 2. * eps * *n + *pivmin * 4.;
/*     [JAN/28/2009] remove the line below since SPDIAM variable not use */
/*     SPDIAM = GU - GL */
/*     Input arguments for SLAEBZ: */
/*     The relative tolerance.  An interval (a,b] lies within */
/*     "relative tolerance" if  b-a < RELTOL*max(|a|,|b|), */
#line 483 "slarrd.f"
    rtoli = *reltol;
/*     Set the absolute tolerance for interval convergence to zero to force */
/*     interval convergence based on relative size of the interval. */
/*     This is dangerous because intervals might not converge when RELTOL is */
/*     small. But at least a very small number should be selected so that for */
/*     strongly graded matrices, the code can get relatively accurate */
/*     eigenvalues. */
#line 490 "slarrd.f"
    atoli = uflow * 4. + *pivmin * 4.;
#line 492 "slarrd.f"
    if (irange == 3) {
/*        RANGE='I': Compute an interval containing eigenvalues */
/*        IL through IU. The initial interval [GL,GU] from the global */
/*        Gerschgorin bounds GL and GU is refined by SLAEBZ. */
#line 497 "slarrd.f"
	itmax = (integer) ((log(tnorm + *pivmin) - log(*pivmin)) / log(2.)) + 
		2;
#line 499 "slarrd.f"
	work[*n + 1] = gl;
#line 500 "slarrd.f"
	work[*n + 2] = gl;
#line 501 "slarrd.f"
	work[*n + 3] = gu;
#line 502 "slarrd.f"
	work[*n + 4] = gu;
#line 503 "slarrd.f"
	work[*n + 5] = gl;
#line 504 "slarrd.f"
	work[*n + 6] = gu;
#line 505 "slarrd.f"
	iwork[1] = -1;
#line 506 "slarrd.f"
	iwork[2] = -1;
#line 507 "slarrd.f"
	iwork[3] = *n + 1;
#line 508 "slarrd.f"
	iwork[4] = *n + 1;
#line 509 "slarrd.f"
	iwork[5] = *il - 1;
#line 510 "slarrd.f"
	iwork[6] = *iu;

#line 512 "slarrd.f"
	slaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, pivmin, &
		d__[1], &e[1], &e2[1], &iwork[5], &work[*n + 1], &work[*n + 5]
		, &iout, &iwork[1], &w[1], &iblock[1], &iinfo);
#line 515 "slarrd.f"
	if (iinfo != 0) {
#line 516 "slarrd.f"
	    *info = iinfo;
#line 517 "slarrd.f"
	    return 0;
#line 518 "slarrd.f"
	}
/*        On exit, output intervals may not be ordered by ascending negcount */
#line 520 "slarrd.f"
	if (iwork[6] == *iu) {
#line 521 "slarrd.f"
	    *wl = work[*n + 1];
#line 522 "slarrd.f"
	    wlu = work[*n + 3];
#line 523 "slarrd.f"
	    nwl = iwork[1];
#line 524 "slarrd.f"
	    *wu = work[*n + 4];
#line 525 "slarrd.f"
	    wul = work[*n + 2];
#line 526 "slarrd.f"
	    nwu = iwork[4];
#line 527 "slarrd.f"
	} else {
#line 528 "slarrd.f"
	    *wl = work[*n + 2];
#line 529 "slarrd.f"
	    wlu = work[*n + 4];
#line 530 "slarrd.f"
	    nwl = iwork[2];
#line 531 "slarrd.f"
	    *wu = work[*n + 3];
#line 532 "slarrd.f"
	    wul = work[*n + 1];
#line 533 "slarrd.f"
	    nwu = iwork[3];
#line 534 "slarrd.f"
	}
/*        On exit, the interval [WL, WLU] contains a value with negcount NWL, */
/*        and [WUL, WU] contains a value with negcount NWU. */
#line 537 "slarrd.f"
	if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
#line 538 "slarrd.f"
	    *info = 4;
#line 539 "slarrd.f"
	    return 0;
#line 540 "slarrd.f"
	}
#line 542 "slarrd.f"
    } else if (irange == 2) {
#line 543 "slarrd.f"
	*wl = *vl;
#line 544 "slarrd.f"
	*wu = *vu;
#line 546 "slarrd.f"
    } else if (irange == 1) {
#line 547 "slarrd.f"
	*wl = gl;
#line 548 "slarrd.f"
	*wu = gu;
#line 549 "slarrd.f"
    }
/*     Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU. */
/*     NWL accumulates the number of eigenvalues .le. WL, */
/*     NWU accumulates the number of eigenvalues .le. WU */
#line 556 "slarrd.f"
    *m = 0;
#line 557 "slarrd.f"
    iend = 0;
#line 558 "slarrd.f"
    *info = 0;
#line 559 "slarrd.f"
    nwl = 0;
#line 560 "slarrd.f"
    nwu = 0;

#line 562 "slarrd.f"
    i__1 = *nsplit;
#line 562 "slarrd.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 563 "slarrd.f"
	ioff = iend;
#line 564 "slarrd.f"
	ibegin = ioff + 1;
#line 565 "slarrd.f"
	iend = isplit[jblk];
#line 566 "slarrd.f"
	in = iend - ioff;

#line 568 "slarrd.f"
	if (in == 1) {
/*           1x1 block */
#line 570 "slarrd.f"
	    if (*wl >= d__[ibegin] - *pivmin) {
#line 570 "slarrd.f"
		++nwl;
#line 570 "slarrd.f"
	    }
#line 572 "slarrd.f"
	    if (*wu >= d__[ibegin] - *pivmin) {
#line 572 "slarrd.f"
		++nwu;
#line 572 "slarrd.f"
	    }
#line 574 "slarrd.f"
	    if (irange == 1 || *wl < d__[ibegin] - *pivmin && *wu >= d__[
		    ibegin] - *pivmin) {
#line 577 "slarrd.f"
		++(*m);
#line 578 "slarrd.f"
		w[*m] = d__[ibegin];
#line 579 "slarrd.f"
		werr[*m] = 0.;
/*              The gap for a single block doesn't matter for the later */
/*              algorithm and is assigned an arbitrary large value */
#line 582 "slarrd.f"
		iblock[*m] = jblk;
#line 583 "slarrd.f"
		indexw[*m] = 1;
#line 584 "slarrd.f"
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
#line 632 "slarrd.f"
	} else {
/*           General Case - block of size IN >= 2 */
/*           Compute local Gerschgorin interval and use it as the initial */
/*           interval for SLAEBZ */
#line 636 "slarrd.f"
	    gu = d__[ibegin];
#line 637 "slarrd.f"
	    gl = d__[ibegin];
#line 638 "slarrd.f"
	    tmp1 = 0.;
#line 640 "slarrd.f"
	    i__2 = iend;
#line 640 "slarrd.f"
	    for (j = ibegin; j <= i__2; ++j) {
/* Computing MIN */
#line 641 "slarrd.f"
		d__1 = gl, d__2 = gers[(j << 1) - 1];
#line 641 "slarrd.f"
		gl = min(d__1,d__2);
/* Computing MAX */
#line 642 "slarrd.f"
		d__1 = gu, d__2 = gers[j * 2];
#line 642 "slarrd.f"
		gu = max(d__1,d__2);
#line 643 "slarrd.f"
/* L40: */
#line 643 "slarrd.f"
	    }
/*           [JAN/28/2009] */
/*           change SPDIAM by TNORM in lines 2 and 3 thereafter */
/*           line 1: remove computation of SPDIAM (not useful anymore) */
/*           SPDIAM = GU - GL */
/*           GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN */
/*           GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN */
#line 650 "slarrd.f"
	    gl = gl - tnorm * 2. * eps * in - *pivmin * 2.;
#line 651 "slarrd.f"
	    gu = gu + tnorm * 2. * eps * in + *pivmin * 2.;

#line 653 "slarrd.f"
	    if (irange > 1) {
#line 654 "slarrd.f"
		if (gu < *wl) {
/*                 the local block contains none of the wanted eigenvalues */
#line 656 "slarrd.f"
		    nwl += in;
#line 657 "slarrd.f"
		    nwu += in;
#line 658 "slarrd.f"
		    goto L70;
#line 659 "slarrd.f"
		}
/*              refine search interval if possible, only range (WL,WU] matters */
#line 661 "slarrd.f"
		gl = max(gl,*wl);
#line 662 "slarrd.f"
		gu = min(gu,*wu);
#line 663 "slarrd.f"
		if (gl >= gu) {
#line 663 "slarrd.f"
		    goto L70;
#line 663 "slarrd.f"
		}
#line 665 "slarrd.f"
	    }
/*           Find negcount of initial interval boundaries GL and GU */
#line 668 "slarrd.f"
	    work[*n + 1] = gl;
#line 669 "slarrd.f"
	    work[*n + in + 1] = gu;
#line 670 "slarrd.f"
	    slaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, 
		    pivmin, &d__[ibegin], &e[ibegin], &e2[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &im, &iwork[1], &
		    w[*m + 1], &iblock[*m + 1], &iinfo);
#line 674 "slarrd.f"
	    if (iinfo != 0) {
#line 675 "slarrd.f"
		*info = iinfo;
#line 676 "slarrd.f"
		return 0;
#line 677 "slarrd.f"
	    }

#line 679 "slarrd.f"
	    nwl += iwork[1];
#line 680 "slarrd.f"
	    nwu += iwork[in + 1];
#line 681 "slarrd.f"
	    iwoff = *m - iwork[1];
/*           Compute Eigenvalues */
#line 684 "slarrd.f"
	    itmax = (integer) ((log(gu - gl + *pivmin) - log(*pivmin)) / log(
		    2.)) + 2;
#line 686 "slarrd.f"
	    slaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, 
		    pivmin, &d__[ibegin], &e[ibegin], &e2[ibegin], idumma, &
		    work[*n + 1], &work[*n + (in << 1) + 1], &iout, &iwork[1],
		     &w[*m + 1], &iblock[*m + 1], &iinfo);
#line 690 "slarrd.f"
	    if (iinfo != 0) {
#line 691 "slarrd.f"
		*info = iinfo;
#line 692 "slarrd.f"
		return 0;
#line 693 "slarrd.f"
	    }

/*           Copy eigenvalues into W and IBLOCK */
/*           Use -JBLK for block number for unconverged eigenvalues. */
/*           Loop over the number of output intervals from SLAEBZ */
#line 698 "slarrd.f"
	    i__2 = iout;
#line 698 "slarrd.f"
	    for (j = 1; j <= i__2; ++j) {
/*              eigenvalue approximation is middle point of interval */
#line 700 "slarrd.f"
		tmp1 = (work[j + *n] + work[j + in + *n]) * .5;
/*              semi length of error interval */
#line 702 "slarrd.f"
		tmp2 = (d__1 = work[j + *n] - work[j + in + *n], abs(d__1)) * 
			.5;
#line 703 "slarrd.f"
		if (j > iout - iinfo) {
/*                 Flag non-convergence. */
#line 705 "slarrd.f"
		    ncnvrg = TRUE_;
#line 706 "slarrd.f"
		    ib = -jblk;
#line 707 "slarrd.f"
		} else {
#line 708 "slarrd.f"
		    ib = jblk;
#line 709 "slarrd.f"
		}
#line 710 "slarrd.f"
		i__3 = iwork[j + in] + iwoff;
#line 710 "slarrd.f"
		for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
#line 712 "slarrd.f"
		    w[je] = tmp1;
#line 713 "slarrd.f"
		    werr[je] = tmp2;
#line 714 "slarrd.f"
		    indexw[je] = je - iwoff;
#line 715 "slarrd.f"
		    iblock[je] = ib;
#line 716 "slarrd.f"
/* L50: */
#line 716 "slarrd.f"
		}
#line 717 "slarrd.f"
/* L60: */
#line 717 "slarrd.f"
	    }

#line 719 "slarrd.f"
	    *m += im;
#line 720 "slarrd.f"
	}
#line 721 "slarrd.f"
L70:
#line 721 "slarrd.f"
	;
#line 721 "slarrd.f"
    }
/*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU */
/*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues. */
#line 725 "slarrd.f"
    if (irange == 3) {
#line 726 "slarrd.f"
	idiscl = *il - 1 - nwl;
#line 727 "slarrd.f"
	idiscu = nwu - *iu;

#line 729 "slarrd.f"
	if (idiscl > 0) {
#line 730 "slarrd.f"
	    im = 0;
#line 731 "slarrd.f"
	    i__1 = *m;
#line 731 "slarrd.f"
	    for (je = 1; je <= i__1; ++je) {
/*              Remove some of the smallest eigenvalues from the left so that */
/*              at the end IDISCL =0. Move all eigenvalues up to the left. */
#line 734 "slarrd.f"
		if (w[je] <= wlu && idiscl > 0) {
#line 735 "slarrd.f"
		    --idiscl;
#line 736 "slarrd.f"
		} else {
#line 737 "slarrd.f"
		    ++im;
#line 738 "slarrd.f"
		    w[im] = w[je];
#line 739 "slarrd.f"
		    werr[im] = werr[je];
#line 740 "slarrd.f"
		    indexw[im] = indexw[je];
#line 741 "slarrd.f"
		    iblock[im] = iblock[je];
#line 742 "slarrd.f"
		}
#line 743 "slarrd.f"
/* L80: */
#line 743 "slarrd.f"
	    }
#line 744 "slarrd.f"
	    *m = im;
#line 745 "slarrd.f"
	}
#line 746 "slarrd.f"
	if (idiscu > 0) {
/*           Remove some of the largest eigenvalues from the right so that */
/*           at the end IDISCU =0. Move all eigenvalues up to the left. */
#line 749 "slarrd.f"
	    im = *m + 1;
#line 750 "slarrd.f"
	    for (je = *m; je >= 1; --je) {
#line 751 "slarrd.f"
		if (w[je] >= wul && idiscu > 0) {
#line 752 "slarrd.f"
		    --idiscu;
#line 753 "slarrd.f"
		} else {
#line 754 "slarrd.f"
		    --im;
#line 755 "slarrd.f"
		    w[im] = w[je];
#line 756 "slarrd.f"
		    werr[im] = werr[je];
#line 757 "slarrd.f"
		    indexw[im] = indexw[je];
#line 758 "slarrd.f"
		    iblock[im] = iblock[je];
#line 759 "slarrd.f"
		}
#line 760 "slarrd.f"
/* L81: */
#line 760 "slarrd.f"
	    }
#line 761 "slarrd.f"
	    jee = 0;
#line 762 "slarrd.f"
	    i__1 = *m;
#line 762 "slarrd.f"
	    for (je = im; je <= i__1; ++je) {
#line 763 "slarrd.f"
		++jee;
#line 764 "slarrd.f"
		w[jee] = w[je];
#line 765 "slarrd.f"
		werr[jee] = werr[je];
#line 766 "slarrd.f"
		indexw[jee] = indexw[je];
#line 767 "slarrd.f"
		iblock[jee] = iblock[je];
#line 768 "slarrd.f"
/* L82: */
#line 768 "slarrd.f"
	    }
#line 769 "slarrd.f"
	    *m = *m - im + 1;
#line 770 "slarrd.f"
	}
#line 772 "slarrd.f"
	if (idiscl > 0 || idiscu > 0) {
/*           Code to deal with effects of bad arithmetic. (If N(w) is */
/*           monotone non-decreasing, this should never happen.) */
/*           Some low eigenvalues to be discarded are not in (WL,WLU], */
/*           or high eigenvalues to be discarded are not in (WUL,WU] */
/*           so just kill off the smallest IDISCL/largest IDISCU */
/*           eigenvalues, by marking the corresponding IBLOCK = 0 */
#line 779 "slarrd.f"
	    if (idiscl > 0) {
#line 780 "slarrd.f"
		wkill = *wu;
#line 781 "slarrd.f"
		i__1 = idiscl;
#line 781 "slarrd.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 782 "slarrd.f"
		    iw = 0;
#line 783 "slarrd.f"
		    i__2 = *m;
#line 783 "slarrd.f"
		    for (je = 1; je <= i__2; ++je) {
#line 784 "slarrd.f"
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
#line 786 "slarrd.f"
			    iw = je;
#line 787 "slarrd.f"
			    wkill = w[je];
#line 788 "slarrd.f"
			}
#line 789 "slarrd.f"
/* L90: */
#line 789 "slarrd.f"
		    }
#line 790 "slarrd.f"
		    iblock[iw] = 0;
#line 791 "slarrd.f"
/* L100: */
#line 791 "slarrd.f"
		}
#line 792 "slarrd.f"
	    }
#line 793 "slarrd.f"
	    if (idiscu > 0) {
#line 794 "slarrd.f"
		wkill = *wl;
#line 795 "slarrd.f"
		i__1 = idiscu;
#line 795 "slarrd.f"
		for (jdisc = 1; jdisc <= i__1; ++jdisc) {
#line 796 "slarrd.f"
		    iw = 0;
#line 797 "slarrd.f"
		    i__2 = *m;
#line 797 "slarrd.f"
		    for (je = 1; je <= i__2; ++je) {
#line 798 "slarrd.f"
			if (iblock[je] != 0 && (w[je] >= wkill || iw == 0)) {
#line 800 "slarrd.f"
			    iw = je;
#line 801 "slarrd.f"
			    wkill = w[je];
#line 802 "slarrd.f"
			}
#line 803 "slarrd.f"
/* L110: */
#line 803 "slarrd.f"
		    }
#line 804 "slarrd.f"
		    iblock[iw] = 0;
#line 805 "slarrd.f"
/* L120: */
#line 805 "slarrd.f"
		}
#line 806 "slarrd.f"
	    }
/*           Now erase all eigenvalues with IBLOCK set to zero */
#line 808 "slarrd.f"
	    im = 0;
#line 809 "slarrd.f"
	    i__1 = *m;
#line 809 "slarrd.f"
	    for (je = 1; je <= i__1; ++je) {
#line 810 "slarrd.f"
		if (iblock[je] != 0) {
#line 811 "slarrd.f"
		    ++im;
#line 812 "slarrd.f"
		    w[im] = w[je];
#line 813 "slarrd.f"
		    werr[im] = werr[je];
#line 814 "slarrd.f"
		    indexw[im] = indexw[je];
#line 815 "slarrd.f"
		    iblock[im] = iblock[je];
#line 816 "slarrd.f"
		}
#line 817 "slarrd.f"
/* L130: */
#line 817 "slarrd.f"
	    }
#line 818 "slarrd.f"
	    *m = im;
#line 819 "slarrd.f"
	}
#line 820 "slarrd.f"
	if (idiscl < 0 || idiscu < 0) {
#line 821 "slarrd.f"
	    toofew = TRUE_;
#line 822 "slarrd.f"
	}
#line 823 "slarrd.f"
    }

#line 825 "slarrd.f"
    if (irange == 1 && *m != *n || irange == 3 && *m != *iu - *il + 1) {
#line 827 "slarrd.f"
	toofew = TRUE_;
#line 828 "slarrd.f"
    }
/*     If ORDER='B', do nothing the eigenvalues are already sorted by */
/*        block. */
/*     If ORDER='E', sort the eigenvalues from smallest to largest */
#line 834 "slarrd.f"
    if (lsame_(order, "E", (ftnlen)1, (ftnlen)1) && *nsplit > 1) {
#line 835 "slarrd.f"
	i__1 = *m - 1;
#line 835 "slarrd.f"
	for (je = 1; je <= i__1; ++je) {
#line 836 "slarrd.f"
	    ie = 0;
#line 837 "slarrd.f"
	    tmp1 = w[je];
#line 838 "slarrd.f"
	    i__2 = *m;
#line 838 "slarrd.f"
	    for (j = je + 1; j <= i__2; ++j) {
#line 839 "slarrd.f"
		if (w[j] < tmp1) {
#line 840 "slarrd.f"
		    ie = j;
#line 841 "slarrd.f"
		    tmp1 = w[j];
#line 842 "slarrd.f"
		}
#line 843 "slarrd.f"
/* L140: */
#line 843 "slarrd.f"
	    }
#line 844 "slarrd.f"
	    if (ie != 0) {
#line 845 "slarrd.f"
		tmp2 = werr[ie];
#line 846 "slarrd.f"
		itmp1 = iblock[ie];
#line 847 "slarrd.f"
		itmp2 = indexw[ie];
#line 848 "slarrd.f"
		w[ie] = w[je];
#line 849 "slarrd.f"
		werr[ie] = werr[je];
#line 850 "slarrd.f"
		iblock[ie] = iblock[je];
#line 851 "slarrd.f"
		indexw[ie] = indexw[je];
#line 852 "slarrd.f"
		w[je] = tmp1;
#line 853 "slarrd.f"
		werr[je] = tmp2;
#line 854 "slarrd.f"
		iblock[je] = itmp1;
#line 855 "slarrd.f"
		indexw[je] = itmp2;
#line 856 "slarrd.f"
	    }
#line 857 "slarrd.f"
/* L150: */
#line 857 "slarrd.f"
	}
#line 858 "slarrd.f"
    }

#line 860 "slarrd.f"
    *info = 0;
#line 861 "slarrd.f"
    if (ncnvrg) {
#line 861 "slarrd.f"
	++(*info);
#line 861 "slarrd.f"
    }
#line 863 "slarrd.f"
    if (toofew) {
#line 863 "slarrd.f"
	*info += 2;
#line 863 "slarrd.f"
    }
#line 865 "slarrd.f"
    return 0;

/*     End of SLARRD */

} /* slarrd_ */

