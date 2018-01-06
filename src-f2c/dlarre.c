#line 1 "dlarre.f"
/* dlarre.f -- translated by f2c (version 20100827).
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

#line 1 "dlarre.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b DLARRE given the tridiagonal matrix T, sets small off-diagonal elements to zero and for each un
reduced block Ti, finds base representations and eigenvalues. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarre.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarre.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarre.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRE( RANGE, N, VL, VU, IL, IU, D, E, E2, */
/*                           RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, M, */
/*                           W, WERR, WGAP, IBLOCK, INDEXW, GERS, PIVMIN, */
/*                           WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          RANGE */
/*       INTEGER            IL, INFO, IU, M, N, NSPLIT */
/*       DOUBLE PRECISION  PIVMIN, RTOL1, RTOL2, SPLTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ), */
/*      $                   INDEXW( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), E2( * ), GERS( * ), */
/*      $                   W( * ),WERR( * ), WGAP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > To find the desired eigenvalues of a given real symmetric */
/* > tridiagonal matrix T, DLARRE sets any "small" off-diagonal */
/* > elements to zero, and for each unreduced block T_i, it finds */
/* > (a) a suitable shift at one end of the block's spectrum, */
/* > (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and */
/* > (c) eigenvalues of each L_i D_i L_i^T. */
/* > The representations and eigenvalues found are then used by */
/* > DSTEMR to compute the eigenvectors of T. */
/* > The accuracy varies depending on whether bisection is used to */
/* > find a few eigenvalues or the dqds algorithm (subroutine DLASQ2) to */
/* > conpute all and then discard any unwanted one. */
/* > As an added benefit, DLARRE also outputs the n */
/* > Gerschgorin intervals for the matrices L_i D_i L_i^T. */
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
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix. N > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION */
/* >          If RANGE='V', the lower bound for the eigenvalues. */
/* >          Eigenvalues less than or equal to VL, or greater than VU, */
/* >          will not be returned.  VL < VU. */
/* >          If RANGE='I' or ='A', DLARRE computes bounds on the desired */
/* >          part of the spectrum. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          If RANGE='V', the upper bound for the eigenvalues. */
/* >          Eigenvalues less than or equal to VL, or greater than VU, */
/* >          will not be returned.  VL < VU. */
/* >          If RANGE='I' or ='A', DLARRE computes bounds on the desired */
/* >          part of the spectrum. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the N diagonal elements of the tridiagonal */
/* >          matrix T. */
/* >          On exit, the N diagonal elements of the diagonal */
/* >          matrices D_i. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the first (N-1) entries contain the subdiagonal */
/* >          elements of the tridiagonal matrix T; E(N) need not be set. */
/* >          On exit, E contains the subdiagonal elements of the unit */
/* >          bidiagonal matrices L_i. The entries E( ISPLIT( I ) ), */
/* >          1 <= I <= NSPLIT, contain the base points sigma_i on output. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E2 */
/* > \verbatim */
/* >          E2 is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the first (N-1) entries contain the SQUARES of the */
/* >          subdiagonal elements of the tridiagonal matrix T; */
/* >          E2(N) need not be set. */
/* >          On exit, the entries E2( ISPLIT( I ) ), */
/* >          1 <= I <= NSPLIT, have been set to zero */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL1 */
/* > \verbatim */
/* >          RTOL1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL2 */
/* > \verbatim */
/* >          RTOL2 is DOUBLE PRECISION */
/* >           Parameters for bisection. */
/* >           An interval [LEFT,RIGHT] has converged if */
/* >           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) ) */
/* > \endverbatim */
/* > */
/* > \param[in] SPLTOL */
/* > \verbatim */
/* >          SPLTOL is DOUBLE PRECISION */
/* >          The threshold for splitting. */
/* > \endverbatim */
/* > */
/* > \param[out] NSPLIT */
/* > \verbatim */
/* >          NSPLIT is INTEGER */
/* >          The number of blocks T splits into. 1 <= NSPLIT <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] ISPLIT */
/* > \verbatim */
/* >          ISPLIT is INTEGER array, dimension (N) */
/* >          The splitting points, at which T breaks up into blocks. */
/* >          The first block consists of rows/columns 1 to ISPLIT(1), */
/* >          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2), */
/* >          etc., and the NSPLIT-th consists of rows/columns */
/* >          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The total number of eigenvalues (of all L_i D_i L_i^T) */
/* >          found. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          The first M elements contain the eigenvalues. The */
/* >          eigenvalues of each of the blocks, L_i D_i L_i^T, are */
/* >          sorted in ascending order ( DLARRE may use the */
/* >          remaining N-M elements as workspace). */
/* > \endverbatim */
/* > */
/* > \param[out] WERR */
/* > \verbatim */
/* >          WERR is DOUBLE PRECISION array, dimension (N) */
/* >          The error bound on the corresponding eigenvalue in W. */
/* > \endverbatim */
/* > */
/* > \param[out] WGAP */
/* > \verbatim */
/* >          WGAP is DOUBLE PRECISION array, dimension (N) */
/* >          The separation from the right neighbor eigenvalue in W. */
/* >          The gap is only with respect to the eigenvalues of the same block */
/* >          as each block has its own representation tree. */
/* >          Exception: at the right end of a block we store the left gap */
/* > \endverbatim */
/* > */
/* > \param[out] IBLOCK */
/* > \verbatim */
/* >          IBLOCK is INTEGER array, dimension (N) */
/* >          The indices of the blocks (submatrices) associated with the */
/* >          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue */
/* >          W(i) belongs to the first block from the top, =2 if W(i) */
/* >          belongs to the second block, etc. */
/* > \endverbatim */
/* > */
/* > \param[out] INDEXW */
/* > \verbatim */
/* >          INDEXW is INTEGER array, dimension (N) */
/* >          The indices of the eigenvalues within each block (submatrix); */
/* >          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the */
/* >          i-th eigenvalue W(i) is the 10-th eigenvalue in block 2 */
/* > \endverbatim */
/* > */
/* > \param[out] GERS */
/* > \verbatim */
/* >          GERS is DOUBLE PRECISION array, dimension (2*N) */
/* >          The N Gerschgorin intervals (the i-th Gerschgorin interval */
/* >          is (GERS(2*i-1), GERS(2*i)). */
/* > \endverbatim */
/* > */
/* > \param[out] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
/* >          The minimum pivot in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (6*N) */
/* >          Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (5*N) */
/* >          Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          > 0:  A problem occurred in DLARRE. */
/* >          < 0:  One of the called subroutines signaled an internal problem. */
/* >                Needs inspection of the corresponding parameter IINFO */
/* >                for further information. */
/* > */
/* >          =-1:  Problem in DLARRD. */
/* >          = 2:  No base representation could be found in MAXTRY iterations. */
/* >                Increasing MAXTRY and recompilation might be a remedy. */
/* >          =-3:  Problem in DLARRB when computing the refined root */
/* >                representation for DLASQ2. */
/* >          =-4:  Problem in DLARRB when preforming bisection on the */
/* >                desired part of the spectrum. */
/* >          =-5:  Problem in DLASQ2. */
/* >          =-6:  Problem in DLASQ2. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The base representations are required to suffer very little */
/* >  element growth and consequently define all their eigenvalues to */
/* >  high relative accuracy. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Beresford Parlett, University of California, Berkeley, USA \n */
/* >     Jim Demmel, University of California, Berkeley, USA \n */
/* >     Inderjit Dhillon, University of Texas, Austin, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */
/* >     Christof Voemel, University of California, Berkeley, USA \n */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlarre_(char *range, integer *n, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *d__, doublereal 
	*e, doublereal *e2, doublereal *rtol1, doublereal *rtol2, doublereal *
	spltol, integer *nsplit, integer *isplit, integer *m, doublereal *w, 
	doublereal *werr, doublereal *wgap, integer *iblock, integer *indexw, 
	doublereal *gers, doublereal *pivmin, doublereal *work, integer *
	iwork, integer *info, ftnlen range_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal s1, s2;
    static integer mb;
    static doublereal gl;
    static integer in, mm;
    static doublereal gu;
    static integer cnt;
    static doublereal eps, tau, tmp, rtl;
    static integer cnt1, cnt2;
    static doublereal tmp1, eabs;
    static integer iend, jblk;
    static doublereal eold;
    static integer indl;
    static doublereal dmax__, emax;
    static integer wend, idum, indu;
    static doublereal rtol;
    static integer iseed[4];
    static doublereal avgap, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical norep;
    extern /* Subroutine */ int dlasq2_(integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static integer ibegin;
    static logical forceb;
    static integer irange;
    static doublereal sgndef;
    extern /* Subroutine */ int dlarra_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *), dlarrb_(integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *), dlarrc_(char *
	    , integer *, doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *, integer *, integer *, integer *, 
	    ftnlen);
    static integer wbegin;
    extern /* Subroutine */ int dlarrd_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static doublereal safmin, spdiam;
    extern /* Subroutine */ int dlarrk_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static logical usedqd;
    static doublereal clwdth, isleft;
    extern /* Subroutine */ int dlarnv_(integer *, integer *, integer *, 
	    doublereal *);
    static doublereal isrght, bsrtol, dpivot;


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

#line 372 "dlarre.f"
    /* Parameter adjustments */
#line 372 "dlarre.f"
    --iwork;
#line 372 "dlarre.f"
    --work;
#line 372 "dlarre.f"
    --gers;
#line 372 "dlarre.f"
    --indexw;
#line 372 "dlarre.f"
    --iblock;
#line 372 "dlarre.f"
    --wgap;
#line 372 "dlarre.f"
    --werr;
#line 372 "dlarre.f"
    --w;
#line 372 "dlarre.f"
    --isplit;
#line 372 "dlarre.f"
    --e2;
#line 372 "dlarre.f"
    --e;
#line 372 "dlarre.f"
    --d__;
#line 372 "dlarre.f"

#line 372 "dlarre.f"
    /* Function Body */
#line 372 "dlarre.f"
    *info = 0;

/*     Decode RANGE */

#line 377 "dlarre.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 378 "dlarre.f"
	irange = 1;
#line 379 "dlarre.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 380 "dlarre.f"
	irange = 3;
#line 381 "dlarre.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 382 "dlarre.f"
	irange = 2;
#line 383 "dlarre.f"
    }
#line 385 "dlarre.f"
    *m = 0;
/*     Get machine constants */
#line 388 "dlarre.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 389 "dlarre.f"
    eps = dlamch_("P", (ftnlen)1);
/*     Set parameters */
#line 392 "dlarre.f"
    rtl = sqrt(eps);
#line 393 "dlarre.f"
    bsrtol = sqrt(eps);
/*     Treat case of 1x1 matrix for quick return */
#line 396 "dlarre.f"
    if (*n == 1) {
#line 397 "dlarre.f"
	if (irange == 1 || irange == 3 && d__[1] > *vl && d__[1] <= *vu || 
		irange == 2 && *il == 1 && *iu == 1) {
#line 400 "dlarre.f"
	    *m = 1;
#line 401 "dlarre.f"
	    w[1] = d__[1];
/*           The computation error of the eigenvalue is zero */
#line 403 "dlarre.f"
	    werr[1] = 0.;
#line 404 "dlarre.f"
	    wgap[1] = 0.;
#line 405 "dlarre.f"
	    iblock[1] = 1;
#line 406 "dlarre.f"
	    indexw[1] = 1;
#line 407 "dlarre.f"
	    gers[1] = d__[1];
#line 408 "dlarre.f"
	    gers[2] = d__[1];
#line 409 "dlarre.f"
	}
/*        store the shift for the initial RRR, which is zero in this case */
#line 411 "dlarre.f"
	e[1] = 0.;
#line 412 "dlarre.f"
	return 0;
#line 413 "dlarre.f"
    }
/*     General case: tridiagonal matrix of order > 1 */

/*     Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter. */
/*     Compute maximum off-diagonal entry and pivmin. */
#line 419 "dlarre.f"
    gl = d__[1];
#line 420 "dlarre.f"
    gu = d__[1];
#line 421 "dlarre.f"
    eold = 0.;
#line 422 "dlarre.f"
    emax = 0.;
#line 423 "dlarre.f"
    e[*n] = 0.;
#line 424 "dlarre.f"
    i__1 = *n;
#line 424 "dlarre.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 425 "dlarre.f"
	werr[i__] = 0.;
#line 426 "dlarre.f"
	wgap[i__] = 0.;
#line 427 "dlarre.f"
	eabs = (d__1 = e[i__], abs(d__1));
#line 428 "dlarre.f"
	if (eabs >= emax) {
#line 429 "dlarre.f"
	    emax = eabs;
#line 430 "dlarre.f"
	}
#line 431 "dlarre.f"
	tmp1 = eabs + eold;
#line 432 "dlarre.f"
	gers[(i__ << 1) - 1] = d__[i__] - tmp1;
/* Computing MIN */
#line 433 "dlarre.f"
	d__1 = gl, d__2 = gers[(i__ << 1) - 1];
#line 433 "dlarre.f"
	gl = min(d__1,d__2);
#line 434 "dlarre.f"
	gers[i__ * 2] = d__[i__] + tmp1;
/* Computing MAX */
#line 435 "dlarre.f"
	d__1 = gu, d__2 = gers[i__ * 2];
#line 435 "dlarre.f"
	gu = max(d__1,d__2);
#line 436 "dlarre.f"
	eold = eabs;
#line 437 "dlarre.f"
/* L5: */
#line 437 "dlarre.f"
    }
/*     The minimum pivot allowed in the Sturm sequence for T */
/* Computing MAX */
/* Computing 2nd power */
#line 439 "dlarre.f"
    d__3 = emax;
#line 439 "dlarre.f"
    d__1 = 1., d__2 = d__3 * d__3;
#line 439 "dlarre.f"
    *pivmin = safmin * max(d__1,d__2);
/*     Compute spectral diameter. The Gerschgorin bounds give an */
/*     estimate that is wrong by at most a factor of SQRT(2) */
#line 442 "dlarre.f"
    spdiam = gu - gl;
/*     Compute splitting points */
#line 445 "dlarre.f"
    dlarra_(n, &d__[1], &e[1], &e2[1], spltol, &spdiam, nsplit, &isplit[1], &
	    iinfo);
/*     Can force use of bisection instead of faster DQDS. */
/*     Option left in the code for future multisection work. */
#line 450 "dlarre.f"
    forceb = FALSE_;
/*     Initialize USEDQD, DQDS should be used for ALLRNG unless someone */
/*     explicitly wants bisection. */
#line 454 "dlarre.f"
    usedqd = irange == 1 && ! forceb;
#line 456 "dlarre.f"
    if (irange == 1 && ! forceb) {
/*        Set interval [VL,VU] that contains all eigenvalues */
#line 458 "dlarre.f"
	*vl = gl;
#line 459 "dlarre.f"
	*vu = gu;
#line 460 "dlarre.f"
    } else {
/*        We call DLARRD to find crude approximations to the eigenvalues */
/*        in the desired range. In case IRANGE = INDRNG, we also obtain the */
/*        interval (VL,VU] that contains all the wanted eigenvalues. */
/*        An interval [LEFT,RIGHT] has converged if */
/*        RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT)) */
/*        DLARRD needs a WORK of size 4*N, IWORK of size 3*N */
#line 467 "dlarre.f"
	dlarrd_(range, "B", n, vl, vu, il, iu, &gers[1], &bsrtol, &d__[1], &e[
		1], &e2[1], pivmin, nsplit, &isplit[1], &mm, &w[1], &werr[1], 
		vl, vu, &iblock[1], &indexw[1], &work[1], &iwork[1], &iinfo, (
		ftnlen)1, (ftnlen)1);
#line 471 "dlarre.f"
	if (iinfo != 0) {
#line 472 "dlarre.f"
	    *info = -1;
#line 473 "dlarre.f"
	    return 0;
#line 474 "dlarre.f"
	}
/*        Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0 */
#line 476 "dlarre.f"
	i__1 = *n;
#line 476 "dlarre.f"
	for (i__ = mm + 1; i__ <= i__1; ++i__) {
#line 477 "dlarre.f"
	    w[i__] = 0.;
#line 478 "dlarre.f"
	    werr[i__] = 0.;
#line 479 "dlarre.f"
	    iblock[i__] = 0;
#line 480 "dlarre.f"
	    indexw[i__] = 0;
#line 481 "dlarre.f"
/* L14: */
#line 481 "dlarre.f"
	}
#line 482 "dlarre.f"
    }
/* ** */
/*     Loop over unreduced blocks */
#line 487 "dlarre.f"
    ibegin = 1;
#line 488 "dlarre.f"
    wbegin = 1;
#line 489 "dlarre.f"
    i__1 = *nsplit;
#line 489 "dlarre.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 490 "dlarre.f"
	iend = isplit[jblk];
#line 491 "dlarre.f"
	in = iend - ibegin + 1;
/*        1 X 1 block */
#line 494 "dlarre.f"
	if (in == 1) {
#line 495 "dlarre.f"
	    if (irange == 1 || irange == 3 && d__[ibegin] > *vl && d__[ibegin]
		     <= *vu || irange == 2 && iblock[wbegin] == jblk) {
#line 499 "dlarre.f"
		++(*m);
#line 500 "dlarre.f"
		w[*m] = d__[ibegin];
#line 501 "dlarre.f"
		werr[*m] = 0.;
/*              The gap for a single block doesn't matter for the later */
/*              algorithm and is assigned an arbitrary large value */
#line 504 "dlarre.f"
		wgap[*m] = 0.;
#line 505 "dlarre.f"
		iblock[*m] = jblk;
#line 506 "dlarre.f"
		indexw[*m] = 1;
#line 507 "dlarre.f"
		++wbegin;
#line 508 "dlarre.f"
	    }
/*           E( IEND ) holds the shift for the initial RRR */
#line 510 "dlarre.f"
	    e[iend] = 0.;
#line 511 "dlarre.f"
	    ibegin = iend + 1;
#line 512 "dlarre.f"
	    goto L170;
#line 513 "dlarre.f"
	}

/*        Blocks of size larger than 1x1 */

/*        E( IEND ) will hold the shift for the initial RRR, for now set it =0 */
#line 518 "dlarre.f"
	e[iend] = 0.;

/*        Find local outer bounds GL,GU for the block */
#line 521 "dlarre.f"
	gl = d__[ibegin];
#line 522 "dlarre.f"
	gu = d__[ibegin];
#line 523 "dlarre.f"
	i__2 = iend;
#line 523 "dlarre.f"
	for (i__ = ibegin; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 524 "dlarre.f"
	    d__1 = gers[(i__ << 1) - 1];
#line 524 "dlarre.f"
	    gl = min(d__1,gl);
/* Computing MAX */
#line 525 "dlarre.f"
	    d__1 = gers[i__ * 2];
#line 525 "dlarre.f"
	    gu = max(d__1,gu);
#line 526 "dlarre.f"
/* L15: */
#line 526 "dlarre.f"
	}
#line 527 "dlarre.f"
	spdiam = gu - gl;
#line 529 "dlarre.f"
	if (! (irange == 1 && ! forceb)) {
/*           Count the number of eigenvalues in the current block. */
#line 531 "dlarre.f"
	    mb = 0;
#line 532 "dlarre.f"
	    i__2 = mm;
#line 532 "dlarre.f"
	    for (i__ = wbegin; i__ <= i__2; ++i__) {
#line 533 "dlarre.f"
		if (iblock[i__] == jblk) {
#line 534 "dlarre.f"
		    ++mb;
#line 535 "dlarre.f"
		} else {
#line 536 "dlarre.f"
		    goto L21;
#line 537 "dlarre.f"
		}
#line 538 "dlarre.f"
/* L20: */
#line 538 "dlarre.f"
	    }
#line 539 "dlarre.f"
L21:
#line 541 "dlarre.f"
	    if (mb == 0) {
/*              No eigenvalue in the current block lies in the desired range */
/*              E( IEND ) holds the shift for the initial RRR */
#line 544 "dlarre.f"
		e[iend] = 0.;
#line 545 "dlarre.f"
		ibegin = iend + 1;
#line 546 "dlarre.f"
		goto L170;
#line 547 "dlarre.f"
	    } else {
/*              Decide whether dqds or bisection is more efficient */
#line 550 "dlarre.f"
		usedqd = (doublereal) mb > in * .5 && ! forceb;
#line 551 "dlarre.f"
		wend = wbegin + mb - 1;
/*              Calculate gaps for the current block */
/*              In later stages, when representations for individual */
/*              eigenvalues are different, we use SIGMA = E( IEND ). */
#line 555 "dlarre.f"
		sigma = 0.;
#line 556 "dlarre.f"
		i__2 = wend - 1;
#line 556 "dlarre.f"
		for (i__ = wbegin; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 557 "dlarre.f"
		    d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + 
			    werr[i__]);
#line 557 "dlarre.f"
		    wgap[i__] = max(d__1,d__2);
#line 559 "dlarre.f"
/* L30: */
#line 559 "dlarre.f"
		}
/* Computing MAX */
#line 560 "dlarre.f"
		d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
#line 560 "dlarre.f"
		wgap[wend] = max(d__1,d__2);
/*              Find local index of the first and last desired evalue. */
#line 563 "dlarre.f"
		indl = indexw[wbegin];
#line 564 "dlarre.f"
		indu = indexw[wend];
#line 565 "dlarre.f"
	    }
#line 566 "dlarre.f"
	}
#line 567 "dlarre.f"
	if (irange == 1 && ! forceb || usedqd) {
/*           Case of DQDS */
/*           Find approximations to the extremal eigenvalues of the block */
#line 570 "dlarre.f"
	    dlarrk_(&in, &c__1, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &
		    rtl, &tmp, &tmp1, &iinfo);
#line 572 "dlarre.f"
	    if (iinfo != 0) {
#line 573 "dlarre.f"
		*info = -1;
#line 574 "dlarre.f"
		return 0;
#line 575 "dlarre.f"
	    }
/* Computing MAX */
#line 576 "dlarre.f"
	    d__2 = gl, d__3 = tmp - tmp1 - eps * 100. * (d__1 = tmp - tmp1, 
		    abs(d__1));
#line 576 "dlarre.f"
	    isleft = max(d__2,d__3);
#line 579 "dlarre.f"
	    dlarrk_(&in, &in, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &
		    rtl, &tmp, &tmp1, &iinfo);
#line 581 "dlarre.f"
	    if (iinfo != 0) {
#line 582 "dlarre.f"
		*info = -1;
#line 583 "dlarre.f"
		return 0;
#line 584 "dlarre.f"
	    }
/* Computing MIN */
#line 585 "dlarre.f"
	    d__2 = gu, d__3 = tmp + tmp1 + eps * 100. * (d__1 = tmp + tmp1, 
		    abs(d__1));
#line 585 "dlarre.f"
	    isrght = min(d__2,d__3);
/*           Improve the estimate of the spectral diameter */
#line 588 "dlarre.f"
	    spdiam = isrght - isleft;
#line 589 "dlarre.f"
	} else {
/*           Case of bisection */
/*           Find approximations to the wanted extremal eigenvalues */
/* Computing MAX */
#line 592 "dlarre.f"
	    d__2 = gl, d__3 = w[wbegin] - werr[wbegin] - eps * 100. * (d__1 = 
		    w[wbegin] - werr[wbegin], abs(d__1));
#line 592 "dlarre.f"
	    isleft = max(d__2,d__3);
/* Computing MIN */
#line 594 "dlarre.f"
	    d__2 = gu, d__3 = w[wend] + werr[wend] + eps * 100. * (d__1 = w[
		    wend] + werr[wend], abs(d__1));
#line 594 "dlarre.f"
	    isrght = min(d__2,d__3);
#line 596 "dlarre.f"
	}
/*        Decide whether the base representation for the current block */
/*        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I */
/*        should be on the left or the right end of the current block. */
/*        The strategy is to shift to the end which is "more populated" */
/*        Furthermore, decide whether to use DQDS for the computation of */
/*        the eigenvalue approximations at the end of DLARRE or bisection. */
/*        dqds is chosen if all eigenvalues are desired or the number of */
/*        eigenvalues to be computed is large compared to the blocksize. */
#line 607 "dlarre.f"
	if (irange == 1 && ! forceb) {
/*           If all the eigenvalues have to be computed, we use dqd */
#line 609 "dlarre.f"
	    usedqd = TRUE_;
/*           INDL is the local index of the first eigenvalue to compute */
#line 611 "dlarre.f"
	    indl = 1;
#line 612 "dlarre.f"
	    indu = in;
/*           MB =  number of eigenvalues to compute */
#line 614 "dlarre.f"
	    mb = in;
#line 615 "dlarre.f"
	    wend = wbegin + mb - 1;
/*           Define 1/4 and 3/4 points of the spectrum */
#line 617 "dlarre.f"
	    s1 = isleft + spdiam * .25;
#line 618 "dlarre.f"
	    s2 = isrght - spdiam * .25;
#line 619 "dlarre.f"
	} else {
/*           DLARRD has computed IBLOCK and INDEXW for each eigenvalue */
/*           approximation. */
/*           choose sigma */
#line 623 "dlarre.f"
	    if (usedqd) {
#line 624 "dlarre.f"
		s1 = isleft + spdiam * .25;
#line 625 "dlarre.f"
		s2 = isrght - spdiam * .25;
#line 626 "dlarre.f"
	    } else {
#line 627 "dlarre.f"
		tmp = min(isrght,*vu) - max(isleft,*vl);
#line 628 "dlarre.f"
		s1 = max(isleft,*vl) + tmp * .25;
#line 629 "dlarre.f"
		s2 = min(isrght,*vu) - tmp * .25;
#line 630 "dlarre.f"
	    }
#line 631 "dlarre.f"
	}
/*        Compute the negcount at the 1/4 and 3/4 points */
#line 634 "dlarre.f"
	if (mb > 1) {
#line 635 "dlarre.f"
	    dlarrc_("T", &in, &s1, &s2, &d__[ibegin], &e[ibegin], pivmin, &
		    cnt, &cnt1, &cnt2, &iinfo, (ftnlen)1);
#line 637 "dlarre.f"
	}
#line 639 "dlarre.f"
	if (mb == 1) {
#line 640 "dlarre.f"
	    sigma = gl;
#line 641 "dlarre.f"
	    sgndef = 1.;
#line 642 "dlarre.f"
	} else if (cnt1 - indl >= indu - cnt2) {
#line 643 "dlarre.f"
	    if (irange == 1 && ! forceb) {
#line 644 "dlarre.f"
		sigma = max(isleft,gl);
#line 645 "dlarre.f"
	    } else if (usedqd) {
/*              use Gerschgorin bound as shift to get pos def matrix */
/*              for dqds */
#line 648 "dlarre.f"
		sigma = isleft;
#line 649 "dlarre.f"
	    } else {
/*              use approximation of the first desired eigenvalue of the */
/*              block as shift */
#line 652 "dlarre.f"
		sigma = max(isleft,*vl);
#line 653 "dlarre.f"
	    }
#line 654 "dlarre.f"
	    sgndef = 1.;
#line 655 "dlarre.f"
	} else {
#line 656 "dlarre.f"
	    if (irange == 1 && ! forceb) {
#line 657 "dlarre.f"
		sigma = min(isrght,gu);
#line 658 "dlarre.f"
	    } else if (usedqd) {
/*              use Gerschgorin bound as shift to get neg def matrix */
/*              for dqds */
#line 661 "dlarre.f"
		sigma = isrght;
#line 662 "dlarre.f"
	    } else {
/*              use approximation of the first desired eigenvalue of the */
/*              block as shift */
#line 665 "dlarre.f"
		sigma = min(isrght,*vu);
#line 666 "dlarre.f"
	    }
#line 667 "dlarre.f"
	    sgndef = -1.;
#line 668 "dlarre.f"
	}
/*        An initial SIGMA has been chosen that will be used for computing */
/*        T - SIGMA I = L D L^T */
/*        Define the increment TAU of the shift in case the initial shift */
/*        needs to be refined to obtain a factorization with not too much */
/*        element growth. */
#line 676 "dlarre.f"
	if (usedqd) {
/*           The initial SIGMA was to the outer end of the spectrum */
/*           the matrix is definite and we need not retreat. */
#line 679 "dlarre.f"
	    tau = spdiam * eps * *n + *pivmin * 2.;
/* Computing MAX */
#line 680 "dlarre.f"
	    d__1 = tau, d__2 = eps * 2. * abs(sigma);
#line 680 "dlarre.f"
	    tau = max(d__1,d__2);
#line 681 "dlarre.f"
	} else {
#line 682 "dlarre.f"
	    if (mb > 1) {
#line 683 "dlarre.f"
		clwdth = w[wend] + werr[wend] - w[wbegin] - werr[wbegin];
#line 684 "dlarre.f"
		avgap = (d__1 = clwdth / (doublereal) (wend - wbegin), abs(
			d__1));
#line 685 "dlarre.f"
		if (sgndef == 1.) {
/* Computing MAX */
#line 686 "dlarre.f"
		    d__1 = wgap[wbegin];
#line 686 "dlarre.f"
		    tau = max(d__1,avgap) * .5;
/* Computing MAX */
#line 687 "dlarre.f"
		    d__1 = tau, d__2 = werr[wbegin];
#line 687 "dlarre.f"
		    tau = max(d__1,d__2);
#line 688 "dlarre.f"
		} else {
/* Computing MAX */
#line 689 "dlarre.f"
		    d__1 = wgap[wend - 1];
#line 689 "dlarre.f"
		    tau = max(d__1,avgap) * .5;
/* Computing MAX */
#line 690 "dlarre.f"
		    d__1 = tau, d__2 = werr[wend];
#line 690 "dlarre.f"
		    tau = max(d__1,d__2);
#line 691 "dlarre.f"
		}
#line 692 "dlarre.f"
	    } else {
#line 693 "dlarre.f"
		tau = werr[wbegin];
#line 694 "dlarre.f"
	    }
#line 695 "dlarre.f"
	}

#line 697 "dlarre.f"
	for (idum = 1; idum <= 6; ++idum) {
/*           Compute L D L^T factorization of tridiagonal matrix T - sigma I. */
/*           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of */
/*           pivots in WORK(2*IN+1:3*IN) */
#line 701 "dlarre.f"
	    dpivot = d__[ibegin] - sigma;
#line 702 "dlarre.f"
	    work[1] = dpivot;
#line 703 "dlarre.f"
	    dmax__ = abs(work[1]);
#line 704 "dlarre.f"
	    j = ibegin;
#line 705 "dlarre.f"
	    i__2 = in - 1;
#line 705 "dlarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 706 "dlarre.f"
		work[(in << 1) + i__] = 1. / work[i__];
#line 707 "dlarre.f"
		tmp = e[j] * work[(in << 1) + i__];
#line 708 "dlarre.f"
		work[in + i__] = tmp;
#line 709 "dlarre.f"
		dpivot = d__[j + 1] - sigma - tmp * e[j];
#line 710 "dlarre.f"
		work[i__ + 1] = dpivot;
/* Computing MAX */
#line 711 "dlarre.f"
		d__1 = dmax__, d__2 = abs(dpivot);
#line 711 "dlarre.f"
		dmax__ = max(d__1,d__2);
#line 712 "dlarre.f"
		++j;
#line 713 "dlarre.f"
/* L70: */
#line 713 "dlarre.f"
	    }
/*           check for element growth */
#line 715 "dlarre.f"
	    if (dmax__ > spdiam * 64.) {
#line 716 "dlarre.f"
		norep = TRUE_;
#line 717 "dlarre.f"
	    } else {
#line 718 "dlarre.f"
		norep = FALSE_;
#line 719 "dlarre.f"
	    }
#line 720 "dlarre.f"
	    if (usedqd && ! norep) {
/*              Ensure the definiteness of the representation */
/*              All entries of D (of L D L^T) must have the same sign */
#line 723 "dlarre.f"
		i__2 = in;
#line 723 "dlarre.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 724 "dlarre.f"
		    tmp = sgndef * work[i__];
#line 725 "dlarre.f"
		    if (tmp < 0.) {
#line 725 "dlarre.f"
			norep = TRUE_;
#line 725 "dlarre.f"
		    }
#line 726 "dlarre.f"
/* L71: */
#line 726 "dlarre.f"
		}
#line 727 "dlarre.f"
	    }
#line 728 "dlarre.f"
	    if (norep) {
/*              Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin */
/*              shift which makes the matrix definite. So we should end up */
/*              here really only in the case of IRANGE = VALRNG or INDRNG. */
#line 732 "dlarre.f"
		if (idum == 5) {
#line 733 "dlarre.f"
		    if (sgndef == 1.) {
/*                    The fudged Gerschgorin shift should succeed */
#line 735 "dlarre.f"
			sigma = gl - spdiam * 2. * eps * *n - *pivmin * 4.;
#line 737 "dlarre.f"
		    } else {
#line 738 "dlarre.f"
			sigma = gu + spdiam * 2. * eps * *n + *pivmin * 4.;
#line 740 "dlarre.f"
		    }
#line 741 "dlarre.f"
		} else {
#line 742 "dlarre.f"
		    sigma -= sgndef * tau;
#line 743 "dlarre.f"
		    tau *= 2.;
#line 744 "dlarre.f"
		}
#line 745 "dlarre.f"
	    } else {
/*              an initial RRR is found */
#line 747 "dlarre.f"
		goto L83;
#line 748 "dlarre.f"
	    }
#line 749 "dlarre.f"
/* L80: */
#line 749 "dlarre.f"
	}
/*        if the program reaches this point, no base representation could be */
/*        found in MAXTRY iterations. */
#line 752 "dlarre.f"
	*info = 2;
#line 753 "dlarre.f"
	return 0;
#line 755 "dlarre.f"
L83:
/*        At this point, we have found an initial base representation */
/*        T - SIGMA I = L D L^T with not too much element growth. */
/*        Store the shift. */
#line 759 "dlarre.f"
	e[iend] = sigma;
/*        Store D and L. */
#line 761 "dlarre.f"
	dcopy_(&in, &work[1], &c__1, &d__[ibegin], &c__1);
#line 762 "dlarre.f"
	i__2 = in - 1;
#line 762 "dlarre.f"
	dcopy_(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
#line 765 "dlarre.f"
	if (mb > 1) {

/*           Perturb each entry of the base representation by a small */
/*           (but random) relative amount to overcome difficulties with */
/*           glued matrices. */

#line 771 "dlarre.f"
	    for (i__ = 1; i__ <= 4; ++i__) {
#line 772 "dlarre.f"
		iseed[i__ - 1] = 1;
#line 773 "dlarre.f"
/* L122: */
#line 773 "dlarre.f"
	    }
#line 775 "dlarre.f"
	    i__2 = (in << 1) - 1;
#line 775 "dlarre.f"
	    dlarnv_(&c__2, iseed, &i__2, &work[1]);
#line 776 "dlarre.f"
	    i__2 = in - 1;
#line 776 "dlarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 777 "dlarre.f"
		d__[ibegin + i__ - 1] *= eps * 8. * work[i__] + 1.;
#line 778 "dlarre.f"
		e[ibegin + i__ - 1] *= eps * 8. * work[in + i__] + 1.;
#line 779 "dlarre.f"
/* L125: */
#line 779 "dlarre.f"
	    }
#line 780 "dlarre.f"
	    d__[iend] *= eps * 4. * work[in] + 1.;

#line 782 "dlarre.f"
	}

/*        Don't update the Gerschgorin intervals because keeping track */
/*        of the updates would be too much work in DLARRV. */
/*        We update W instead and use it to locate the proper Gerschgorin */
/*        intervals. */
/*        Compute the required eigenvalues of L D L' by bisection or dqds */
#line 790 "dlarre.f"
	if (! usedqd) {
/*           If DLARRD has been used, shift the eigenvalue approximations */
/*           according to their representation. This is necessary for */
/*           a uniform DLARRV since dqds computes eigenvalues of the */
/*           shifted representation. In DLARRV, W will always hold the */
/*           UNshifted eigenvalue approximation. */
#line 796 "dlarre.f"
	    i__2 = wend;
#line 796 "dlarre.f"
	    for (j = wbegin; j <= i__2; ++j) {
#line 797 "dlarre.f"
		w[j] -= sigma;
#line 798 "dlarre.f"
		werr[j] += (d__1 = w[j], abs(d__1)) * eps;
#line 799 "dlarre.f"
/* L134: */
#line 799 "dlarre.f"
	    }
/*           call DLARRB to reduce eigenvalue error of the approximations */
/*           from DLARRD */
#line 802 "dlarre.f"
	    i__2 = iend - 1;
#line 802 "dlarre.f"
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
/* Computing 2nd power */
#line 803 "dlarre.f"
		d__1 = e[i__];
#line 803 "dlarre.f"
		work[i__] = d__[i__] * (d__1 * d__1);
#line 804 "dlarre.f"
/* L135: */
#line 804 "dlarre.f"
	    }
/*           use bisection to find EV from INDL to INDU */
#line 806 "dlarre.f"
	    i__2 = indl - 1;
#line 806 "dlarre.f"
	    dlarrb_(&in, &d__[ibegin], &work[ibegin], &indl, &indu, rtol1, 
		    rtol2, &i__2, &w[wbegin], &wgap[wbegin], &werr[wbegin], &
		    work[(*n << 1) + 1], &iwork[1], pivmin, &spdiam, &in, &
		    iinfo);
#line 811 "dlarre.f"
	    if (iinfo != 0) {
#line 812 "dlarre.f"
		*info = -4;
#line 813 "dlarre.f"
		return 0;
#line 814 "dlarre.f"
	    }
/*           DLARRB computes all gaps correctly except for the last one */
/*           Record distance to VU/GU */
/* Computing MAX */
#line 817 "dlarre.f"
	    d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
#line 817 "dlarre.f"
	    wgap[wend] = max(d__1,d__2);
#line 819 "dlarre.f"
	    i__2 = indu;
#line 819 "dlarre.f"
	    for (i__ = indl; i__ <= i__2; ++i__) {
#line 820 "dlarre.f"
		++(*m);
#line 821 "dlarre.f"
		iblock[*m] = jblk;
#line 822 "dlarre.f"
		indexw[*m] = i__;
#line 823 "dlarre.f"
/* L138: */
#line 823 "dlarre.f"
	    }
#line 824 "dlarre.f"
	} else {
/*           Call dqds to get all eigs (and then possibly delete unwanted */
/*           eigenvalues). */
/*           Note that dqds finds the eigenvalues of the L D L^T representation */
/*           of T to high relative accuracy. High relative accuracy */
/*           might be lost when the shift of the RRR is subtracted to obtain */
/*           the eigenvalues of T. However, T is not guaranteed to define its */
/*           eigenvalues to high relative accuracy anyway. */
/*           Set RTOL to the order of the tolerance used in DLASQ2 */
/*           This is an ESTIMATED error, the worst case bound is 4*N*EPS */
/*           which is usually too large and requires unnecessary work to be */
/*           done by bisection when computing the eigenvectors */
#line 836 "dlarre.f"
	    rtol = log((doublereal) in) * 4. * eps;
#line 837 "dlarre.f"
	    j = ibegin;
#line 838 "dlarre.f"
	    i__2 = in - 1;
#line 838 "dlarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 839 "dlarre.f"
		work[(i__ << 1) - 1] = (d__1 = d__[j], abs(d__1));
#line 840 "dlarre.f"
		work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
#line 841 "dlarre.f"
		++j;
#line 842 "dlarre.f"
/* L140: */
#line 842 "dlarre.f"
	    }
#line 843 "dlarre.f"
	    work[(in << 1) - 1] = (d__1 = d__[iend], abs(d__1));
#line 844 "dlarre.f"
	    work[in * 2] = 0.;
#line 845 "dlarre.f"
	    dlasq2_(&in, &work[1], &iinfo);
#line 846 "dlarre.f"
	    if (iinfo != 0) {
/*              If IINFO = -5 then an index is part of a tight cluster */
/*              and should be changed. The index is in IWORK(1) and the */
/*              gap is in WORK(N+1) */
#line 850 "dlarre.f"
		*info = -5;
#line 851 "dlarre.f"
		return 0;
#line 852 "dlarre.f"
	    } else {
/*              Test that all eigenvalues are positive as expected */
#line 854 "dlarre.f"
		i__2 = in;
#line 854 "dlarre.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 855 "dlarre.f"
		    if (work[i__] < 0.) {
#line 856 "dlarre.f"
			*info = -6;
#line 857 "dlarre.f"
			return 0;
#line 858 "dlarre.f"
		    }
#line 859 "dlarre.f"
/* L149: */
#line 859 "dlarre.f"
		}
#line 860 "dlarre.f"
	    }
#line 861 "dlarre.f"
	    if (sgndef > 0.) {
#line 862 "dlarre.f"
		i__2 = indu;
#line 862 "dlarre.f"
		for (i__ = indl; i__ <= i__2; ++i__) {
#line 863 "dlarre.f"
		    ++(*m);
#line 864 "dlarre.f"
		    w[*m] = work[in - i__ + 1];
#line 865 "dlarre.f"
		    iblock[*m] = jblk;
#line 866 "dlarre.f"
		    indexw[*m] = i__;
#line 867 "dlarre.f"
/* L150: */
#line 867 "dlarre.f"
		}
#line 868 "dlarre.f"
	    } else {
#line 869 "dlarre.f"
		i__2 = indu;
#line 869 "dlarre.f"
		for (i__ = indl; i__ <= i__2; ++i__) {
#line 870 "dlarre.f"
		    ++(*m);
#line 871 "dlarre.f"
		    w[*m] = -work[i__];
#line 872 "dlarre.f"
		    iblock[*m] = jblk;
#line 873 "dlarre.f"
		    indexw[*m] = i__;
#line 874 "dlarre.f"
/* L160: */
#line 874 "dlarre.f"
		}
#line 875 "dlarre.f"
	    }
#line 877 "dlarre.f"
	    i__2 = *m;
#line 877 "dlarre.f"
	    for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
/*              the value of RTOL below should be the tolerance in DLASQ2 */
#line 879 "dlarre.f"
		werr[i__] = rtol * (d__1 = w[i__], abs(d__1));
#line 880 "dlarre.f"
/* L165: */
#line 880 "dlarre.f"
	    }
#line 881 "dlarre.f"
	    i__2 = *m - 1;
#line 881 "dlarre.f"
	    for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
/*              compute the right gap between the intervals */
/* Computing MAX */
#line 883 "dlarre.f"
		d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + werr[
			i__]);
#line 883 "dlarre.f"
		wgap[i__] = max(d__1,d__2);
#line 885 "dlarre.f"
/* L166: */
#line 885 "dlarre.f"
	    }
/* Computing MAX */
#line 886 "dlarre.f"
	    d__1 = 0., d__2 = *vu - sigma - (w[*m] + werr[*m]);
#line 886 "dlarre.f"
	    wgap[*m] = max(d__1,d__2);
#line 888 "dlarre.f"
	}
/*        proceed with next block */
#line 890 "dlarre.f"
	ibegin = iend + 1;
#line 891 "dlarre.f"
	wbegin = wend + 1;
#line 892 "dlarre.f"
L170:
#line 892 "dlarre.f"
	;
#line 892 "dlarre.f"
    }

#line 895 "dlarre.f"
    return 0;

/*     end of DLARRE */

} /* dlarre_ */

