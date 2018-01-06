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
/* > \endverbatim */
/* > */
/* > \param[in,out] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          If RANGE='V', the lower and upper bounds for the eigenvalues. */
/* >          Eigenvalues less than or equal to VL, or greater than VU, */
/* >          will not be returned.  VL < VU. */
/* >          If RANGE='I' or ='A', DLARRE computes bounds on the desired */
/* >          part of the spectrum. */
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
/* >          If RANGE='I', the indices (in ascending order) of the */
/* >          smallest and largest eigenvalues to be returned. */
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
/* >          > 0:  A problem occured in DLARRE. */
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

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 364 "dlarre.f"
    /* Parameter adjustments */
#line 364 "dlarre.f"
    --iwork;
#line 364 "dlarre.f"
    --work;
#line 364 "dlarre.f"
    --gers;
#line 364 "dlarre.f"
    --indexw;
#line 364 "dlarre.f"
    --iblock;
#line 364 "dlarre.f"
    --wgap;
#line 364 "dlarre.f"
    --werr;
#line 364 "dlarre.f"
    --w;
#line 364 "dlarre.f"
    --isplit;
#line 364 "dlarre.f"
    --e2;
#line 364 "dlarre.f"
    --e;
#line 364 "dlarre.f"
    --d__;
#line 364 "dlarre.f"

#line 364 "dlarre.f"
    /* Function Body */
#line 364 "dlarre.f"
    *info = 0;

/*     Decode RANGE */

#line 369 "dlarre.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 370 "dlarre.f"
	irange = 1;
#line 371 "dlarre.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 372 "dlarre.f"
	irange = 3;
#line 373 "dlarre.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 374 "dlarre.f"
	irange = 2;
#line 375 "dlarre.f"
    }
#line 377 "dlarre.f"
    *m = 0;
/*     Get machine constants */
#line 380 "dlarre.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 381 "dlarre.f"
    eps = dlamch_("P", (ftnlen)1);
/*     Set parameters */
#line 384 "dlarre.f"
    rtl = sqrt(eps);
#line 385 "dlarre.f"
    bsrtol = sqrt(eps);
/*     Treat case of 1x1 matrix for quick return */
#line 388 "dlarre.f"
    if (*n == 1) {
#line 389 "dlarre.f"
	if (irange == 1 || irange == 3 && d__[1] > *vl && d__[1] <= *vu || 
		irange == 2 && *il == 1 && *iu == 1) {
#line 392 "dlarre.f"
	    *m = 1;
#line 393 "dlarre.f"
	    w[1] = d__[1];
/*           The computation error of the eigenvalue is zero */
#line 395 "dlarre.f"
	    werr[1] = 0.;
#line 396 "dlarre.f"
	    wgap[1] = 0.;
#line 397 "dlarre.f"
	    iblock[1] = 1;
#line 398 "dlarre.f"
	    indexw[1] = 1;
#line 399 "dlarre.f"
	    gers[1] = d__[1];
#line 400 "dlarre.f"
	    gers[2] = d__[1];
#line 401 "dlarre.f"
	}
/*        store the shift for the initial RRR, which is zero in this case */
#line 403 "dlarre.f"
	e[1] = 0.;
#line 404 "dlarre.f"
	return 0;
#line 405 "dlarre.f"
    }
/*     General case: tridiagonal matrix of order > 1 */

/*     Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter. */
/*     Compute maximum off-diagonal entry and pivmin. */
#line 411 "dlarre.f"
    gl = d__[1];
#line 412 "dlarre.f"
    gu = d__[1];
#line 413 "dlarre.f"
    eold = 0.;
#line 414 "dlarre.f"
    emax = 0.;
#line 415 "dlarre.f"
    e[*n] = 0.;
#line 416 "dlarre.f"
    i__1 = *n;
#line 416 "dlarre.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 417 "dlarre.f"
	werr[i__] = 0.;
#line 418 "dlarre.f"
	wgap[i__] = 0.;
#line 419 "dlarre.f"
	eabs = (d__1 = e[i__], abs(d__1));
#line 420 "dlarre.f"
	if (eabs >= emax) {
#line 421 "dlarre.f"
	    emax = eabs;
#line 422 "dlarre.f"
	}
#line 423 "dlarre.f"
	tmp1 = eabs + eold;
#line 424 "dlarre.f"
	gers[(i__ << 1) - 1] = d__[i__] - tmp1;
/* Computing MIN */
#line 425 "dlarre.f"
	d__1 = gl, d__2 = gers[(i__ << 1) - 1];
#line 425 "dlarre.f"
	gl = min(d__1,d__2);
#line 426 "dlarre.f"
	gers[i__ * 2] = d__[i__] + tmp1;
/* Computing MAX */
#line 427 "dlarre.f"
	d__1 = gu, d__2 = gers[i__ * 2];
#line 427 "dlarre.f"
	gu = max(d__1,d__2);
#line 428 "dlarre.f"
	eold = eabs;
#line 429 "dlarre.f"
/* L5: */
#line 429 "dlarre.f"
    }
/*     The minimum pivot allowed in the Sturm sequence for T */
/* Computing MAX */
/* Computing 2nd power */
#line 431 "dlarre.f"
    d__3 = emax;
#line 431 "dlarre.f"
    d__1 = 1., d__2 = d__3 * d__3;
#line 431 "dlarre.f"
    *pivmin = safmin * max(d__1,d__2);
/*     Compute spectral diameter. The Gerschgorin bounds give an */
/*     estimate that is wrong by at most a factor of SQRT(2) */
#line 434 "dlarre.f"
    spdiam = gu - gl;
/*     Compute splitting points */
#line 437 "dlarre.f"
    dlarra_(n, &d__[1], &e[1], &e2[1], spltol, &spdiam, nsplit, &isplit[1], &
	    iinfo);
/*     Can force use of bisection instead of faster DQDS. */
/*     Option left in the code for future multisection work. */
#line 442 "dlarre.f"
    forceb = FALSE_;
/*     Initialize USEDQD, DQDS should be used for ALLRNG unless someone */
/*     explicitly wants bisection. */
#line 446 "dlarre.f"
    usedqd = irange == 1 && ! forceb;
#line 448 "dlarre.f"
    if (irange == 1 && ! forceb) {
/*        Set interval [VL,VU] that contains all eigenvalues */
#line 450 "dlarre.f"
	*vl = gl;
#line 451 "dlarre.f"
	*vu = gu;
#line 452 "dlarre.f"
    } else {
/*        We call DLARRD to find crude approximations to the eigenvalues */
/*        in the desired range. In case IRANGE = INDRNG, we also obtain the */
/*        interval (VL,VU] that contains all the wanted eigenvalues. */
/*        An interval [LEFT,RIGHT] has converged if */
/*        RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT)) */
/*        DLARRD needs a WORK of size 4*N, IWORK of size 3*N */
#line 459 "dlarre.f"
	dlarrd_(range, "B", n, vl, vu, il, iu, &gers[1], &bsrtol, &d__[1], &e[
		1], &e2[1], pivmin, nsplit, &isplit[1], &mm, &w[1], &werr[1], 
		vl, vu, &iblock[1], &indexw[1], &work[1], &iwork[1], &iinfo, (
		ftnlen)1, (ftnlen)1);
#line 463 "dlarre.f"
	if (iinfo != 0) {
#line 464 "dlarre.f"
	    *info = -1;
#line 465 "dlarre.f"
	    return 0;
#line 466 "dlarre.f"
	}
/*        Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0 */
#line 468 "dlarre.f"
	i__1 = *n;
#line 468 "dlarre.f"
	for (i__ = mm + 1; i__ <= i__1; ++i__) {
#line 469 "dlarre.f"
	    w[i__] = 0.;
#line 470 "dlarre.f"
	    werr[i__] = 0.;
#line 471 "dlarre.f"
	    iblock[i__] = 0;
#line 472 "dlarre.f"
	    indexw[i__] = 0;
#line 473 "dlarre.f"
/* L14: */
#line 473 "dlarre.f"
	}
#line 474 "dlarre.f"
    }
/* ** */
/*     Loop over unreduced blocks */
#line 479 "dlarre.f"
    ibegin = 1;
#line 480 "dlarre.f"
    wbegin = 1;
#line 481 "dlarre.f"
    i__1 = *nsplit;
#line 481 "dlarre.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 482 "dlarre.f"
	iend = isplit[jblk];
#line 483 "dlarre.f"
	in = iend - ibegin + 1;
/*        1 X 1 block */
#line 486 "dlarre.f"
	if (in == 1) {
#line 487 "dlarre.f"
	    if (irange == 1 || irange == 3 && d__[ibegin] > *vl && d__[ibegin]
		     <= *vu || irange == 2 && iblock[wbegin] == jblk) {
#line 491 "dlarre.f"
		++(*m);
#line 492 "dlarre.f"
		w[*m] = d__[ibegin];
#line 493 "dlarre.f"
		werr[*m] = 0.;
/*              The gap for a single block doesn't matter for the later */
/*              algorithm and is assigned an arbitrary large value */
#line 496 "dlarre.f"
		wgap[*m] = 0.;
#line 497 "dlarre.f"
		iblock[*m] = jblk;
#line 498 "dlarre.f"
		indexw[*m] = 1;
#line 499 "dlarre.f"
		++wbegin;
#line 500 "dlarre.f"
	    }
/*           E( IEND ) holds the shift for the initial RRR */
#line 502 "dlarre.f"
	    e[iend] = 0.;
#line 503 "dlarre.f"
	    ibegin = iend + 1;
#line 504 "dlarre.f"
	    goto L170;
#line 505 "dlarre.f"
	}

/*        Blocks of size larger than 1x1 */

/*        E( IEND ) will hold the shift for the initial RRR, for now set it =0 */
#line 510 "dlarre.f"
	e[iend] = 0.;

/*        Find local outer bounds GL,GU for the block */
#line 513 "dlarre.f"
	gl = d__[ibegin];
#line 514 "dlarre.f"
	gu = d__[ibegin];
#line 515 "dlarre.f"
	i__2 = iend;
#line 515 "dlarre.f"
	for (i__ = ibegin; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 516 "dlarre.f"
	    d__1 = gers[(i__ << 1) - 1];
#line 516 "dlarre.f"
	    gl = min(d__1,gl);
/* Computing MAX */
#line 517 "dlarre.f"
	    d__1 = gers[i__ * 2];
#line 517 "dlarre.f"
	    gu = max(d__1,gu);
#line 518 "dlarre.f"
/* L15: */
#line 518 "dlarre.f"
	}
#line 519 "dlarre.f"
	spdiam = gu - gl;
#line 521 "dlarre.f"
	if (! (irange == 1 && ! forceb)) {
/*           Count the number of eigenvalues in the current block. */
#line 523 "dlarre.f"
	    mb = 0;
#line 524 "dlarre.f"
	    i__2 = mm;
#line 524 "dlarre.f"
	    for (i__ = wbegin; i__ <= i__2; ++i__) {
#line 525 "dlarre.f"
		if (iblock[i__] == jblk) {
#line 526 "dlarre.f"
		    ++mb;
#line 527 "dlarre.f"
		} else {
#line 528 "dlarre.f"
		    goto L21;
#line 529 "dlarre.f"
		}
#line 530 "dlarre.f"
/* L20: */
#line 530 "dlarre.f"
	    }
#line 531 "dlarre.f"
L21:
#line 533 "dlarre.f"
	    if (mb == 0) {
/*              No eigenvalue in the current block lies in the desired range */
/*              E( IEND ) holds the shift for the initial RRR */
#line 536 "dlarre.f"
		e[iend] = 0.;
#line 537 "dlarre.f"
		ibegin = iend + 1;
#line 538 "dlarre.f"
		goto L170;
#line 539 "dlarre.f"
	    } else {
/*              Decide whether dqds or bisection is more efficient */
#line 542 "dlarre.f"
		usedqd = (doublereal) mb > in * .5 && ! forceb;
#line 543 "dlarre.f"
		wend = wbegin + mb - 1;
/*              Calculate gaps for the current block */
/*              In later stages, when representations for individual */
/*              eigenvalues are different, we use SIGMA = E( IEND ). */
#line 547 "dlarre.f"
		sigma = 0.;
#line 548 "dlarre.f"
		i__2 = wend - 1;
#line 548 "dlarre.f"
		for (i__ = wbegin; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 549 "dlarre.f"
		    d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + 
			    werr[i__]);
#line 549 "dlarre.f"
		    wgap[i__] = max(d__1,d__2);
#line 551 "dlarre.f"
/* L30: */
#line 551 "dlarre.f"
		}
/* Computing MAX */
#line 552 "dlarre.f"
		d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
#line 552 "dlarre.f"
		wgap[wend] = max(d__1,d__2);
/*              Find local index of the first and last desired evalue. */
#line 555 "dlarre.f"
		indl = indexw[wbegin];
#line 556 "dlarre.f"
		indu = indexw[wend];
#line 557 "dlarre.f"
	    }
#line 558 "dlarre.f"
	}
#line 559 "dlarre.f"
	if (irange == 1 && ! forceb || usedqd) {
/*           Case of DQDS */
/*           Find approximations to the extremal eigenvalues of the block */
#line 562 "dlarre.f"
	    dlarrk_(&in, &c__1, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &
		    rtl, &tmp, &tmp1, &iinfo);
#line 564 "dlarre.f"
	    if (iinfo != 0) {
#line 565 "dlarre.f"
		*info = -1;
#line 566 "dlarre.f"
		return 0;
#line 567 "dlarre.f"
	    }
/* Computing MAX */
#line 568 "dlarre.f"
	    d__2 = gl, d__3 = tmp - tmp1 - eps * 100. * (d__1 = tmp - tmp1, 
		    abs(d__1));
#line 568 "dlarre.f"
	    isleft = max(d__2,d__3);
#line 571 "dlarre.f"
	    dlarrk_(&in, &in, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &
		    rtl, &tmp, &tmp1, &iinfo);
#line 573 "dlarre.f"
	    if (iinfo != 0) {
#line 574 "dlarre.f"
		*info = -1;
#line 575 "dlarre.f"
		return 0;
#line 576 "dlarre.f"
	    }
/* Computing MIN */
#line 577 "dlarre.f"
	    d__2 = gu, d__3 = tmp + tmp1 + eps * 100. * (d__1 = tmp + tmp1, 
		    abs(d__1));
#line 577 "dlarre.f"
	    isrght = min(d__2,d__3);
/*           Improve the estimate of the spectral diameter */
#line 580 "dlarre.f"
	    spdiam = isrght - isleft;
#line 581 "dlarre.f"
	} else {
/*           Case of bisection */
/*           Find approximations to the wanted extremal eigenvalues */
/* Computing MAX */
#line 584 "dlarre.f"
	    d__2 = gl, d__3 = w[wbegin] - werr[wbegin] - eps * 100. * (d__1 = 
		    w[wbegin] - werr[wbegin], abs(d__1));
#line 584 "dlarre.f"
	    isleft = max(d__2,d__3);
/* Computing MIN */
#line 586 "dlarre.f"
	    d__2 = gu, d__3 = w[wend] + werr[wend] + eps * 100. * (d__1 = w[
		    wend] + werr[wend], abs(d__1));
#line 586 "dlarre.f"
	    isrght = min(d__2,d__3);
#line 588 "dlarre.f"
	}
/*        Decide whether the base representation for the current block */
/*        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I */
/*        should be on the left or the right end of the current block. */
/*        The strategy is to shift to the end which is "more populated" */
/*        Furthermore, decide whether to use DQDS for the computation of */
/*        the eigenvalue approximations at the end of DLARRE or bisection. */
/*        dqds is chosen if all eigenvalues are desired or the number of */
/*        eigenvalues to be computed is large compared to the blocksize. */
#line 599 "dlarre.f"
	if (irange == 1 && ! forceb) {
/*           If all the eigenvalues have to be computed, we use dqd */
#line 601 "dlarre.f"
	    usedqd = TRUE_;
/*           INDL is the local index of the first eigenvalue to compute */
#line 603 "dlarre.f"
	    indl = 1;
#line 604 "dlarre.f"
	    indu = in;
/*           MB =  number of eigenvalues to compute */
#line 606 "dlarre.f"
	    mb = in;
#line 607 "dlarre.f"
	    wend = wbegin + mb - 1;
/*           Define 1/4 and 3/4 points of the spectrum */
#line 609 "dlarre.f"
	    s1 = isleft + spdiam * .25;
#line 610 "dlarre.f"
	    s2 = isrght - spdiam * .25;
#line 611 "dlarre.f"
	} else {
/*           DLARRD has computed IBLOCK and INDEXW for each eigenvalue */
/*           approximation. */
/*           choose sigma */
#line 615 "dlarre.f"
	    if (usedqd) {
#line 616 "dlarre.f"
		s1 = isleft + spdiam * .25;
#line 617 "dlarre.f"
		s2 = isrght - spdiam * .25;
#line 618 "dlarre.f"
	    } else {
#line 619 "dlarre.f"
		tmp = min(isrght,*vu) - max(isleft,*vl);
#line 620 "dlarre.f"
		s1 = max(isleft,*vl) + tmp * .25;
#line 621 "dlarre.f"
		s2 = min(isrght,*vu) - tmp * .25;
#line 622 "dlarre.f"
	    }
#line 623 "dlarre.f"
	}
/*        Compute the negcount at the 1/4 and 3/4 points */
#line 626 "dlarre.f"
	if (mb > 1) {
#line 627 "dlarre.f"
	    dlarrc_("T", &in, &s1, &s2, &d__[ibegin], &e[ibegin], pivmin, &
		    cnt, &cnt1, &cnt2, &iinfo, (ftnlen)1);
#line 629 "dlarre.f"
	}
#line 631 "dlarre.f"
	if (mb == 1) {
#line 632 "dlarre.f"
	    sigma = gl;
#line 633 "dlarre.f"
	    sgndef = 1.;
#line 634 "dlarre.f"
	} else if (cnt1 - indl >= indu - cnt2) {
#line 635 "dlarre.f"
	    if (irange == 1 && ! forceb) {
#line 636 "dlarre.f"
		sigma = max(isleft,gl);
#line 637 "dlarre.f"
	    } else if (usedqd) {
/*              use Gerschgorin bound as shift to get pos def matrix */
/*              for dqds */
#line 640 "dlarre.f"
		sigma = isleft;
#line 641 "dlarre.f"
	    } else {
/*              use approximation of the first desired eigenvalue of the */
/*              block as shift */
#line 644 "dlarre.f"
		sigma = max(isleft,*vl);
#line 645 "dlarre.f"
	    }
#line 646 "dlarre.f"
	    sgndef = 1.;
#line 647 "dlarre.f"
	} else {
#line 648 "dlarre.f"
	    if (irange == 1 && ! forceb) {
#line 649 "dlarre.f"
		sigma = min(isrght,gu);
#line 650 "dlarre.f"
	    } else if (usedqd) {
/*              use Gerschgorin bound as shift to get neg def matrix */
/*              for dqds */
#line 653 "dlarre.f"
		sigma = isrght;
#line 654 "dlarre.f"
	    } else {
/*              use approximation of the first desired eigenvalue of the */
/*              block as shift */
#line 657 "dlarre.f"
		sigma = min(isrght,*vu);
#line 658 "dlarre.f"
	    }
#line 659 "dlarre.f"
	    sgndef = -1.;
#line 660 "dlarre.f"
	}
/*        An initial SIGMA has been chosen that will be used for computing */
/*        T - SIGMA I = L D L^T */
/*        Define the increment TAU of the shift in case the initial shift */
/*        needs to be refined to obtain a factorization with not too much */
/*        element growth. */
#line 668 "dlarre.f"
	if (usedqd) {
/*           The initial SIGMA was to the outer end of the spectrum */
/*           the matrix is definite and we need not retreat. */
#line 671 "dlarre.f"
	    tau = spdiam * eps * *n + *pivmin * 2.;
/* Computing MAX */
#line 672 "dlarre.f"
	    d__1 = tau, d__2 = eps * 2. * abs(sigma);
#line 672 "dlarre.f"
	    tau = max(d__1,d__2);
#line 673 "dlarre.f"
	} else {
#line 674 "dlarre.f"
	    if (mb > 1) {
#line 675 "dlarre.f"
		clwdth = w[wend] + werr[wend] - w[wbegin] - werr[wbegin];
#line 676 "dlarre.f"
		avgap = (d__1 = clwdth / (doublereal) (wend - wbegin), abs(
			d__1));
#line 677 "dlarre.f"
		if (sgndef == 1.) {
/* Computing MAX */
#line 678 "dlarre.f"
		    d__1 = wgap[wbegin];
#line 678 "dlarre.f"
		    tau = max(d__1,avgap) * .5;
/* Computing MAX */
#line 679 "dlarre.f"
		    d__1 = tau, d__2 = werr[wbegin];
#line 679 "dlarre.f"
		    tau = max(d__1,d__2);
#line 680 "dlarre.f"
		} else {
/* Computing MAX */
#line 681 "dlarre.f"
		    d__1 = wgap[wend - 1];
#line 681 "dlarre.f"
		    tau = max(d__1,avgap) * .5;
/* Computing MAX */
#line 682 "dlarre.f"
		    d__1 = tau, d__2 = werr[wend];
#line 682 "dlarre.f"
		    tau = max(d__1,d__2);
#line 683 "dlarre.f"
		}
#line 684 "dlarre.f"
	    } else {
#line 685 "dlarre.f"
		tau = werr[wbegin];
#line 686 "dlarre.f"
	    }
#line 687 "dlarre.f"
	}

#line 689 "dlarre.f"
	for (idum = 1; idum <= 6; ++idum) {
/*           Compute L D L^T factorization of tridiagonal matrix T - sigma I. */
/*           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of */
/*           pivots in WORK(2*IN+1:3*IN) */
#line 693 "dlarre.f"
	    dpivot = d__[ibegin] - sigma;
#line 694 "dlarre.f"
	    work[1] = dpivot;
#line 695 "dlarre.f"
	    dmax__ = abs(work[1]);
#line 696 "dlarre.f"
	    j = ibegin;
#line 697 "dlarre.f"
	    i__2 = in - 1;
#line 697 "dlarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 698 "dlarre.f"
		work[(in << 1) + i__] = 1. / work[i__];
#line 699 "dlarre.f"
		tmp = e[j] * work[(in << 1) + i__];
#line 700 "dlarre.f"
		work[in + i__] = tmp;
#line 701 "dlarre.f"
		dpivot = d__[j + 1] - sigma - tmp * e[j];
#line 702 "dlarre.f"
		work[i__ + 1] = dpivot;
/* Computing MAX */
#line 703 "dlarre.f"
		d__1 = dmax__, d__2 = abs(dpivot);
#line 703 "dlarre.f"
		dmax__ = max(d__1,d__2);
#line 704 "dlarre.f"
		++j;
#line 705 "dlarre.f"
/* L70: */
#line 705 "dlarre.f"
	    }
/*           check for element growth */
#line 707 "dlarre.f"
	    if (dmax__ > spdiam * 64.) {
#line 708 "dlarre.f"
		norep = TRUE_;
#line 709 "dlarre.f"
	    } else {
#line 710 "dlarre.f"
		norep = FALSE_;
#line 711 "dlarre.f"
	    }
#line 712 "dlarre.f"
	    if (usedqd && ! norep) {
/*              Ensure the definiteness of the representation */
/*              All entries of D (of L D L^T) must have the same sign */
#line 715 "dlarre.f"
		i__2 = in;
#line 715 "dlarre.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 716 "dlarre.f"
		    tmp = sgndef * work[i__];
#line 717 "dlarre.f"
		    if (tmp < 0.) {
#line 717 "dlarre.f"
			norep = TRUE_;
#line 717 "dlarre.f"
		    }
#line 718 "dlarre.f"
/* L71: */
#line 718 "dlarre.f"
		}
#line 719 "dlarre.f"
	    }
#line 720 "dlarre.f"
	    if (norep) {
/*              Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin */
/*              shift which makes the matrix definite. So we should end up */
/*              here really only in the case of IRANGE = VALRNG or INDRNG. */
#line 724 "dlarre.f"
		if (idum == 5) {
#line 725 "dlarre.f"
		    if (sgndef == 1.) {
/*                    The fudged Gerschgorin shift should succeed */
#line 727 "dlarre.f"
			sigma = gl - spdiam * 2. * eps * *n - *pivmin * 4.;
#line 729 "dlarre.f"
		    } else {
#line 730 "dlarre.f"
			sigma = gu + spdiam * 2. * eps * *n + *pivmin * 4.;
#line 732 "dlarre.f"
		    }
#line 733 "dlarre.f"
		} else {
#line 734 "dlarre.f"
		    sigma -= sgndef * tau;
#line 735 "dlarre.f"
		    tau *= 2.;
#line 736 "dlarre.f"
		}
#line 737 "dlarre.f"
	    } else {
/*              an initial RRR is found */
#line 739 "dlarre.f"
		goto L83;
#line 740 "dlarre.f"
	    }
#line 741 "dlarre.f"
/* L80: */
#line 741 "dlarre.f"
	}
/*        if the program reaches this point, no base representation could be */
/*        found in MAXTRY iterations. */
#line 744 "dlarre.f"
	*info = 2;
#line 745 "dlarre.f"
	return 0;
#line 747 "dlarre.f"
L83:
/*        At this point, we have found an initial base representation */
/*        T - SIGMA I = L D L^T with not too much element growth. */
/*        Store the shift. */
#line 751 "dlarre.f"
	e[iend] = sigma;
/*        Store D and L. */
#line 753 "dlarre.f"
	dcopy_(&in, &work[1], &c__1, &d__[ibegin], &c__1);
#line 754 "dlarre.f"
	i__2 = in - 1;
#line 754 "dlarre.f"
	dcopy_(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
#line 757 "dlarre.f"
	if (mb > 1) {

/*           Perturb each entry of the base representation by a small */
/*           (but random) relative amount to overcome difficulties with */
/*           glued matrices. */

#line 763 "dlarre.f"
	    for (i__ = 1; i__ <= 4; ++i__) {
#line 764 "dlarre.f"
		iseed[i__ - 1] = 1;
#line 765 "dlarre.f"
/* L122: */
#line 765 "dlarre.f"
	    }
#line 767 "dlarre.f"
	    i__2 = (in << 1) - 1;
#line 767 "dlarre.f"
	    dlarnv_(&c__2, iseed, &i__2, &work[1]);
#line 768 "dlarre.f"
	    i__2 = in - 1;
#line 768 "dlarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 769 "dlarre.f"
		d__[ibegin + i__ - 1] *= eps * 8. * work[i__] + 1.;
#line 770 "dlarre.f"
		e[ibegin + i__ - 1] *= eps * 8. * work[in + i__] + 1.;
#line 771 "dlarre.f"
/* L125: */
#line 771 "dlarre.f"
	    }
#line 772 "dlarre.f"
	    d__[iend] *= eps * 4. * work[in] + 1.;

#line 774 "dlarre.f"
	}

/*        Don't update the Gerschgorin intervals because keeping track */
/*        of the updates would be too much work in DLARRV. */
/*        We update W instead and use it to locate the proper Gerschgorin */
/*        intervals. */
/*        Compute the required eigenvalues of L D L' by bisection or dqds */
#line 782 "dlarre.f"
	if (! usedqd) {
/*           If DLARRD has been used, shift the eigenvalue approximations */
/*           according to their representation. This is necessary for */
/*           a uniform DLARRV since dqds computes eigenvalues of the */
/*           shifted representation. In DLARRV, W will always hold the */
/*           UNshifted eigenvalue approximation. */
#line 788 "dlarre.f"
	    i__2 = wend;
#line 788 "dlarre.f"
	    for (j = wbegin; j <= i__2; ++j) {
#line 789 "dlarre.f"
		w[j] -= sigma;
#line 790 "dlarre.f"
		werr[j] += (d__1 = w[j], abs(d__1)) * eps;
#line 791 "dlarre.f"
/* L134: */
#line 791 "dlarre.f"
	    }
/*           call DLARRB to reduce eigenvalue error of the approximations */
/*           from DLARRD */
#line 794 "dlarre.f"
	    i__2 = iend - 1;
#line 794 "dlarre.f"
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
/* Computing 2nd power */
#line 795 "dlarre.f"
		d__1 = e[i__];
#line 795 "dlarre.f"
		work[i__] = d__[i__] * (d__1 * d__1);
#line 796 "dlarre.f"
/* L135: */
#line 796 "dlarre.f"
	    }
/*           use bisection to find EV from INDL to INDU */
#line 798 "dlarre.f"
	    i__2 = indl - 1;
#line 798 "dlarre.f"
	    dlarrb_(&in, &d__[ibegin], &work[ibegin], &indl, &indu, rtol1, 
		    rtol2, &i__2, &w[wbegin], &wgap[wbegin], &werr[wbegin], &
		    work[(*n << 1) + 1], &iwork[1], pivmin, &spdiam, &in, &
		    iinfo);
#line 803 "dlarre.f"
	    if (iinfo != 0) {
#line 804 "dlarre.f"
		*info = -4;
#line 805 "dlarre.f"
		return 0;
#line 806 "dlarre.f"
	    }
/*           DLARRB computes all gaps correctly except for the last one */
/*           Record distance to VU/GU */
/* Computing MAX */
#line 809 "dlarre.f"
	    d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
#line 809 "dlarre.f"
	    wgap[wend] = max(d__1,d__2);
#line 811 "dlarre.f"
	    i__2 = indu;
#line 811 "dlarre.f"
	    for (i__ = indl; i__ <= i__2; ++i__) {
#line 812 "dlarre.f"
		++(*m);
#line 813 "dlarre.f"
		iblock[*m] = jblk;
#line 814 "dlarre.f"
		indexw[*m] = i__;
#line 815 "dlarre.f"
/* L138: */
#line 815 "dlarre.f"
	    }
#line 816 "dlarre.f"
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
#line 828 "dlarre.f"
	    rtol = log((doublereal) in) * 4. * eps;
#line 829 "dlarre.f"
	    j = ibegin;
#line 830 "dlarre.f"
	    i__2 = in - 1;
#line 830 "dlarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 831 "dlarre.f"
		work[(i__ << 1) - 1] = (d__1 = d__[j], abs(d__1));
#line 832 "dlarre.f"
		work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
#line 833 "dlarre.f"
		++j;
#line 834 "dlarre.f"
/* L140: */
#line 834 "dlarre.f"
	    }
#line 835 "dlarre.f"
	    work[(in << 1) - 1] = (d__1 = d__[iend], abs(d__1));
#line 836 "dlarre.f"
	    work[in * 2] = 0.;
#line 837 "dlarre.f"
	    dlasq2_(&in, &work[1], &iinfo);
#line 838 "dlarre.f"
	    if (iinfo != 0) {
/*              If IINFO = -5 then an index is part of a tight cluster */
/*              and should be changed. The index is in IWORK(1) and the */
/*              gap is in WORK(N+1) */
#line 842 "dlarre.f"
		*info = -5;
#line 843 "dlarre.f"
		return 0;
#line 844 "dlarre.f"
	    } else {
/*              Test that all eigenvalues are positive as expected */
#line 846 "dlarre.f"
		i__2 = in;
#line 846 "dlarre.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 847 "dlarre.f"
		    if (work[i__] < 0.) {
#line 848 "dlarre.f"
			*info = -6;
#line 849 "dlarre.f"
			return 0;
#line 850 "dlarre.f"
		    }
#line 851 "dlarre.f"
/* L149: */
#line 851 "dlarre.f"
		}
#line 852 "dlarre.f"
	    }
#line 853 "dlarre.f"
	    if (sgndef > 0.) {
#line 854 "dlarre.f"
		i__2 = indu;
#line 854 "dlarre.f"
		for (i__ = indl; i__ <= i__2; ++i__) {
#line 855 "dlarre.f"
		    ++(*m);
#line 856 "dlarre.f"
		    w[*m] = work[in - i__ + 1];
#line 857 "dlarre.f"
		    iblock[*m] = jblk;
#line 858 "dlarre.f"
		    indexw[*m] = i__;
#line 859 "dlarre.f"
/* L150: */
#line 859 "dlarre.f"
		}
#line 860 "dlarre.f"
	    } else {
#line 861 "dlarre.f"
		i__2 = indu;
#line 861 "dlarre.f"
		for (i__ = indl; i__ <= i__2; ++i__) {
#line 862 "dlarre.f"
		    ++(*m);
#line 863 "dlarre.f"
		    w[*m] = -work[i__];
#line 864 "dlarre.f"
		    iblock[*m] = jblk;
#line 865 "dlarre.f"
		    indexw[*m] = i__;
#line 866 "dlarre.f"
/* L160: */
#line 866 "dlarre.f"
		}
#line 867 "dlarre.f"
	    }
#line 869 "dlarre.f"
	    i__2 = *m;
#line 869 "dlarre.f"
	    for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
/*              the value of RTOL below should be the tolerance in DLASQ2 */
#line 871 "dlarre.f"
		werr[i__] = rtol * (d__1 = w[i__], abs(d__1));
#line 872 "dlarre.f"
/* L165: */
#line 872 "dlarre.f"
	    }
#line 873 "dlarre.f"
	    i__2 = *m - 1;
#line 873 "dlarre.f"
	    for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
/*              compute the right gap between the intervals */
/* Computing MAX */
#line 875 "dlarre.f"
		d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + werr[
			i__]);
#line 875 "dlarre.f"
		wgap[i__] = max(d__1,d__2);
#line 877 "dlarre.f"
/* L166: */
#line 877 "dlarre.f"
	    }
/* Computing MAX */
#line 878 "dlarre.f"
	    d__1 = 0., d__2 = *vu - sigma - (w[*m] + werr[*m]);
#line 878 "dlarre.f"
	    wgap[*m] = max(d__1,d__2);
#line 880 "dlarre.f"
	}
/*        proceed with next block */
#line 882 "dlarre.f"
	ibegin = iend + 1;
#line 883 "dlarre.f"
	wbegin = wend + 1;
#line 884 "dlarre.f"
L170:
#line 884 "dlarre.f"
	;
#line 884 "dlarre.f"
    }

#line 887 "dlarre.f"
    return 0;

/*     end of DLARRE */

} /* dlarre_ */

