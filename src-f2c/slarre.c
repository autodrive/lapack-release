#line 1 "slarre.f"
/* slarre.f -- translated by f2c (version 20100827).
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

#line 1 "slarre.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SLARRE given the tridiagonal matrix T, sets small off-diagonal elements to zero and for each un
reduced block Ti, finds base representations and eigenvalues. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARRE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarre.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarre.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarre.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARRE( RANGE, N, VL, VU, IL, IU, D, E, E2, */
/*                           RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, M, */
/*                           W, WERR, WGAP, IBLOCK, INDEXW, GERS, PIVMIN, */
/*                           WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          RANGE */
/*       INTEGER            IL, INFO, IU, M, N, NSPLIT */
/*       REAL               PIVMIN, RTOL1, RTOL2, SPLTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ), */
/*      $                   INDEXW( * ) */
/*       REAL               D( * ), E( * ), E2( * ), GERS( * ), */
/*      $                   W( * ),WERR( * ), WGAP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > To find the desired eigenvalues of a given real symmetric */
/* > tridiagonal matrix T, SLARRE sets any "small" off-diagonal */
/* > elements to zero, and for each unreduced block T_i, it finds */
/* > (a) a suitable shift at one end of the block's spectrum, */
/* > (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and */
/* > (c) eigenvalues of each L_i D_i L_i^T. */
/* > The representations and eigenvalues found are then used by */
/* > SSTEMR to compute the eigenvectors of T. */
/* > The accuracy varies depending on whether bisection is used to */
/* > find a few eigenvalues or the dqds algorithm (subroutine SLASQ2) to */
/* > conpute all and then discard any unwanted one. */
/* > As an added benefit, SLARRE also outputs the n */
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
/* >          VL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] VU */
/* > \verbatim */
/* >          VU is REAL */
/* >          If RANGE='V', the lower and upper bounds for the eigenvalues. */
/* >          Eigenvalues less than or equal to VL, or greater than VU, */
/* >          will not be returned.  VL < VU. */
/* >          If RANGE='I' or ='A', SLARRE computes bounds on the desired */
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
/* >          D is REAL array, dimension (N) */
/* >          On entry, the N diagonal elements of the tridiagonal */
/* >          matrix T. */
/* >          On exit, the N diagonal elements of the diagonal */
/* >          matrices D_i. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* >          On entry, the first (N-1) entries contain the subdiagonal */
/* >          elements of the tridiagonal matrix T; E(N) need not be set. */
/* >          On exit, E contains the subdiagonal elements of the unit */
/* >          bidiagonal matrices L_i. The entries E( ISPLIT( I ) ), */
/* >          1 <= I <= NSPLIT, contain the base points sigma_i on output. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E2 */
/* > \verbatim */
/* >          E2 is REAL array, dimension (N) */
/* >          On entry, the first (N-1) entries contain the SQUARES of the */
/* >          subdiagonal elements of the tridiagonal matrix T; */
/* >          E2(N) need not be set. */
/* >          On exit, the entries E2( ISPLIT( I ) ), */
/* >          1 <= I <= NSPLIT, have been set to zero */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL1 */
/* > \verbatim */
/* >          RTOL1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL2 */
/* > \verbatim */
/* >          RTOL2 is REAL */
/* >           Parameters for bisection. */
/* >           An interval [LEFT,RIGHT] has converged if */
/* >           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) ) */
/* > \endverbatim */
/* > */
/* > \param[in] SPLTOL */
/* > \verbatim */
/* >          SPLTOL is REAL */
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
/* >          W is REAL array, dimension (N) */
/* >          The first M elements contain the eigenvalues. The */
/* >          eigenvalues of each of the blocks, L_i D_i L_i^T, are */
/* >          sorted in ascending order ( SLARRE may use the */
/* >          remaining N-M elements as workspace). */
/* > \endverbatim */
/* > */
/* > \param[out] WERR */
/* > \verbatim */
/* >          WERR is REAL array, dimension (N) */
/* >          The error bound on the corresponding eigenvalue in W. */
/* > \endverbatim */
/* > */
/* > \param[out] WGAP */
/* > \verbatim */
/* >          WGAP is REAL array, dimension (N) */
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
/* >          GERS is REAL array, dimension (2*N) */
/* >          The N Gerschgorin intervals (the i-th Gerschgorin interval */
/* >          is (GERS(2*i-1), GERS(2*i)). */
/* > \endverbatim */
/* > */
/* > \param[out] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
/* >          The minimum pivot in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (6*N) */
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
/* >          > 0:  A problem occured in SLARRE. */
/* >          < 0:  One of the called subroutines signaled an internal problem. */
/* >                Needs inspection of the corresponding parameter IINFO */
/* >                for further information. */
/* > */
/* >          =-1:  Problem in SLARRD. */
/* >          = 2:  No base representation could be found in MAXTRY iterations. */
/* >                Increasing MAXTRY and recompilation might be a remedy. */
/* >          =-3:  Problem in SLARRB when computing the refined root */
/* >                representation for SLASQ2. */
/* >          =-4:  Problem in SLARRB when preforming bisection on the */
/* >                desired part of the spectrum. */
/* >          =-5:  Problem in SLASQ2. */
/* >          =-6:  Problem in SLASQ2. */
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
/* Subroutine */ int slarre_(char *range, integer *n, doublereal *vl, 
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
    static logical norep;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slasq2_(integer *, doublereal *, 
	    integer *);
    static integer ibegin;
    static logical forceb;
    static integer irange;
    static doublereal sgndef;
    extern doublereal slamch_(char *, ftnlen);
    static integer wbegin;
    static doublereal safmin, spdiam;
    extern /* Subroutine */ int slarra_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *);
    static logical usedqd;
    static doublereal clwdth, isleft;
    extern /* Subroutine */ int slarrb_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *), slarrc_(char *
	    , integer *, doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *, integer *, integer *, integer *, 
	    ftnlen), slarrd_(char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), slarrk_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal isrght, bsrtol, dpivot;
    extern /* Subroutine */ int slarnv_(integer *, integer *, integer *, 
	    doublereal *);


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

#line 364 "slarre.f"
    /* Parameter adjustments */
#line 364 "slarre.f"
    --iwork;
#line 364 "slarre.f"
    --work;
#line 364 "slarre.f"
    --gers;
#line 364 "slarre.f"
    --indexw;
#line 364 "slarre.f"
    --iblock;
#line 364 "slarre.f"
    --wgap;
#line 364 "slarre.f"
    --werr;
#line 364 "slarre.f"
    --w;
#line 364 "slarre.f"
    --isplit;
#line 364 "slarre.f"
    --e2;
#line 364 "slarre.f"
    --e;
#line 364 "slarre.f"
    --d__;
#line 364 "slarre.f"

#line 364 "slarre.f"
    /* Function Body */
#line 364 "slarre.f"
    *info = 0;

/*     Decode RANGE */

#line 369 "slarre.f"
    if (lsame_(range, "A", (ftnlen)1, (ftnlen)1)) {
#line 370 "slarre.f"
	irange = 1;
#line 371 "slarre.f"
    } else if (lsame_(range, "V", (ftnlen)1, (ftnlen)1)) {
#line 372 "slarre.f"
	irange = 3;
#line 373 "slarre.f"
    } else if (lsame_(range, "I", (ftnlen)1, (ftnlen)1)) {
#line 374 "slarre.f"
	irange = 2;
#line 375 "slarre.f"
    }
#line 377 "slarre.f"
    *m = 0;
/*     Get machine constants */
#line 380 "slarre.f"
    safmin = slamch_("S", (ftnlen)1);
#line 381 "slarre.f"
    eps = slamch_("P", (ftnlen)1);
/*     Set parameters */
#line 384 "slarre.f"
    rtl = eps * 100.;
/*     If one were ever to ask for less initial precision in BSRTOL, */
/*     one should keep in mind that for the subset case, the extremal */
/*     eigenvalues must be at least as accurate as the current setting */
/*     (eigenvalues in the middle need not as much accuracy) */
#line 389 "slarre.f"
    bsrtol = sqrt(eps) * 5e-4;
/*     Treat case of 1x1 matrix for quick return */
#line 392 "slarre.f"
    if (*n == 1) {
#line 393 "slarre.f"
	if (irange == 1 || irange == 3 && d__[1] > *vl && d__[1] <= *vu || 
		irange == 2 && *il == 1 && *iu == 1) {
#line 396 "slarre.f"
	    *m = 1;
#line 397 "slarre.f"
	    w[1] = d__[1];
/*           The computation error of the eigenvalue is zero */
#line 399 "slarre.f"
	    werr[1] = 0.;
#line 400 "slarre.f"
	    wgap[1] = 0.;
#line 401 "slarre.f"
	    iblock[1] = 1;
#line 402 "slarre.f"
	    indexw[1] = 1;
#line 403 "slarre.f"
	    gers[1] = d__[1];
#line 404 "slarre.f"
	    gers[2] = d__[1];
#line 405 "slarre.f"
	}
/*        store the shift for the initial RRR, which is zero in this case */
#line 407 "slarre.f"
	e[1] = 0.;
#line 408 "slarre.f"
	return 0;
#line 409 "slarre.f"
    }
/*     General case: tridiagonal matrix of order > 1 */

/*     Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter. */
/*     Compute maximum off-diagonal entry and pivmin. */
#line 415 "slarre.f"
    gl = d__[1];
#line 416 "slarre.f"
    gu = d__[1];
#line 417 "slarre.f"
    eold = 0.;
#line 418 "slarre.f"
    emax = 0.;
#line 419 "slarre.f"
    e[*n] = 0.;
#line 420 "slarre.f"
    i__1 = *n;
#line 420 "slarre.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 421 "slarre.f"
	werr[i__] = 0.;
#line 422 "slarre.f"
	wgap[i__] = 0.;
#line 423 "slarre.f"
	eabs = (d__1 = e[i__], abs(d__1));
#line 424 "slarre.f"
	if (eabs >= emax) {
#line 425 "slarre.f"
	    emax = eabs;
#line 426 "slarre.f"
	}
#line 427 "slarre.f"
	tmp1 = eabs + eold;
#line 428 "slarre.f"
	gers[(i__ << 1) - 1] = d__[i__] - tmp1;
/* Computing MIN */
#line 429 "slarre.f"
	d__1 = gl, d__2 = gers[(i__ << 1) - 1];
#line 429 "slarre.f"
	gl = min(d__1,d__2);
#line 430 "slarre.f"
	gers[i__ * 2] = d__[i__] + tmp1;
/* Computing MAX */
#line 431 "slarre.f"
	d__1 = gu, d__2 = gers[i__ * 2];
#line 431 "slarre.f"
	gu = max(d__1,d__2);
#line 432 "slarre.f"
	eold = eabs;
#line 433 "slarre.f"
/* L5: */
#line 433 "slarre.f"
    }
/*     The minimum pivot allowed in the Sturm sequence for T */
/* Computing MAX */
/* Computing 2nd power */
#line 435 "slarre.f"
    d__3 = emax;
#line 435 "slarre.f"
    d__1 = 1., d__2 = d__3 * d__3;
#line 435 "slarre.f"
    *pivmin = safmin * max(d__1,d__2);
/*     Compute spectral diameter. The Gerschgorin bounds give an */
/*     estimate that is wrong by at most a factor of SQRT(2) */
#line 438 "slarre.f"
    spdiam = gu - gl;
/*     Compute splitting points */
#line 441 "slarre.f"
    slarra_(n, &d__[1], &e[1], &e2[1], spltol, &spdiam, nsplit, &isplit[1], &
	    iinfo);
/*     Can force use of bisection instead of faster DQDS. */
/*     Option left in the code for future multisection work. */
#line 446 "slarre.f"
    forceb = FALSE_;
/*     Initialize USEDQD, DQDS should be used for ALLRNG unless someone */
/*     explicitly wants bisection. */
#line 450 "slarre.f"
    usedqd = irange == 1 && ! forceb;
#line 452 "slarre.f"
    if (irange == 1 && ! forceb) {
/*        Set interval [VL,VU] that contains all eigenvalues */
#line 454 "slarre.f"
	*vl = gl;
#line 455 "slarre.f"
	*vu = gu;
#line 456 "slarre.f"
    } else {
/*        We call SLARRD to find crude approximations to the eigenvalues */
/*        in the desired range. In case IRANGE = INDRNG, we also obtain the */
/*        interval (VL,VU] that contains all the wanted eigenvalues. */
/*        An interval [LEFT,RIGHT] has converged if */
/*        RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT)) */
/*        SLARRD needs a WORK of size 4*N, IWORK of size 3*N */
#line 463 "slarre.f"
	slarrd_(range, "B", n, vl, vu, il, iu, &gers[1], &bsrtol, &d__[1], &e[
		1], &e2[1], pivmin, nsplit, &isplit[1], &mm, &w[1], &werr[1], 
		vl, vu, &iblock[1], &indexw[1], &work[1], &iwork[1], &iinfo, (
		ftnlen)1, (ftnlen)1);
#line 467 "slarre.f"
	if (iinfo != 0) {
#line 468 "slarre.f"
	    *info = -1;
#line 469 "slarre.f"
	    return 0;
#line 470 "slarre.f"
	}
/*        Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0 */
#line 472 "slarre.f"
	i__1 = *n;
#line 472 "slarre.f"
	for (i__ = mm + 1; i__ <= i__1; ++i__) {
#line 473 "slarre.f"
	    w[i__] = 0.;
#line 474 "slarre.f"
	    werr[i__] = 0.;
#line 475 "slarre.f"
	    iblock[i__] = 0;
#line 476 "slarre.f"
	    indexw[i__] = 0;
#line 477 "slarre.f"
/* L14: */
#line 477 "slarre.f"
	}
#line 478 "slarre.f"
    }
/* ** */
/*     Loop over unreduced blocks */
#line 483 "slarre.f"
    ibegin = 1;
#line 484 "slarre.f"
    wbegin = 1;
#line 485 "slarre.f"
    i__1 = *nsplit;
#line 485 "slarre.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 486 "slarre.f"
	iend = isplit[jblk];
#line 487 "slarre.f"
	in = iend - ibegin + 1;
/*        1 X 1 block */
#line 490 "slarre.f"
	if (in == 1) {
#line 491 "slarre.f"
	    if (irange == 1 || irange == 3 && d__[ibegin] > *vl && d__[ibegin]
		     <= *vu || irange == 2 && iblock[wbegin] == jblk) {
#line 495 "slarre.f"
		++(*m);
#line 496 "slarre.f"
		w[*m] = d__[ibegin];
#line 497 "slarre.f"
		werr[*m] = 0.;
/*              The gap for a single block doesn't matter for the later */
/*              algorithm and is assigned an arbitrary large value */
#line 500 "slarre.f"
		wgap[*m] = 0.;
#line 501 "slarre.f"
		iblock[*m] = jblk;
#line 502 "slarre.f"
		indexw[*m] = 1;
#line 503 "slarre.f"
		++wbegin;
#line 504 "slarre.f"
	    }
/*           E( IEND ) holds the shift for the initial RRR */
#line 506 "slarre.f"
	    e[iend] = 0.;
#line 507 "slarre.f"
	    ibegin = iend + 1;
#line 508 "slarre.f"
	    goto L170;
#line 509 "slarre.f"
	}

/*        Blocks of size larger than 1x1 */

/*        E( IEND ) will hold the shift for the initial RRR, for now set it =0 */
#line 514 "slarre.f"
	e[iend] = 0.;

/*        Find local outer bounds GL,GU for the block */
#line 517 "slarre.f"
	gl = d__[ibegin];
#line 518 "slarre.f"
	gu = d__[ibegin];
#line 519 "slarre.f"
	i__2 = iend;
#line 519 "slarre.f"
	for (i__ = ibegin; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 520 "slarre.f"
	    d__1 = gers[(i__ << 1) - 1];
#line 520 "slarre.f"
	    gl = min(d__1,gl);
/* Computing MAX */
#line 521 "slarre.f"
	    d__1 = gers[i__ * 2];
#line 521 "slarre.f"
	    gu = max(d__1,gu);
#line 522 "slarre.f"
/* L15: */
#line 522 "slarre.f"
	}
#line 523 "slarre.f"
	spdiam = gu - gl;
#line 525 "slarre.f"
	if (! (irange == 1 && ! forceb)) {
/*           Count the number of eigenvalues in the current block. */
#line 527 "slarre.f"
	    mb = 0;
#line 528 "slarre.f"
	    i__2 = mm;
#line 528 "slarre.f"
	    for (i__ = wbegin; i__ <= i__2; ++i__) {
#line 529 "slarre.f"
		if (iblock[i__] == jblk) {
#line 530 "slarre.f"
		    ++mb;
#line 531 "slarre.f"
		} else {
#line 532 "slarre.f"
		    goto L21;
#line 533 "slarre.f"
		}
#line 534 "slarre.f"
/* L20: */
#line 534 "slarre.f"
	    }
#line 535 "slarre.f"
L21:
#line 537 "slarre.f"
	    if (mb == 0) {
/*              No eigenvalue in the current block lies in the desired range */
/*              E( IEND ) holds the shift for the initial RRR */
#line 540 "slarre.f"
		e[iend] = 0.;
#line 541 "slarre.f"
		ibegin = iend + 1;
#line 542 "slarre.f"
		goto L170;
#line 543 "slarre.f"
	    } else {
/*              Decide whether dqds or bisection is more efficient */
#line 546 "slarre.f"
		usedqd = (doublereal) mb > in * .5 && ! forceb;
#line 547 "slarre.f"
		wend = wbegin + mb - 1;
/*              Calculate gaps for the current block */
/*              In later stages, when representations for individual */
/*              eigenvalues are different, we use SIGMA = E( IEND ). */
#line 551 "slarre.f"
		sigma = 0.;
#line 552 "slarre.f"
		i__2 = wend - 1;
#line 552 "slarre.f"
		for (i__ = wbegin; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 553 "slarre.f"
		    d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + 
			    werr[i__]);
#line 553 "slarre.f"
		    wgap[i__] = max(d__1,d__2);
#line 555 "slarre.f"
/* L30: */
#line 555 "slarre.f"
		}
/* Computing MAX */
#line 556 "slarre.f"
		d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
#line 556 "slarre.f"
		wgap[wend] = max(d__1,d__2);
/*              Find local index of the first and last desired evalue. */
#line 559 "slarre.f"
		indl = indexw[wbegin];
#line 560 "slarre.f"
		indu = indexw[wend];
#line 561 "slarre.f"
	    }
#line 562 "slarre.f"
	}
#line 563 "slarre.f"
	if (irange == 1 && ! forceb || usedqd) {
/*           Case of DQDS */
/*           Find approximations to the extremal eigenvalues of the block */
#line 566 "slarre.f"
	    slarrk_(&in, &c__1, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &
		    rtl, &tmp, &tmp1, &iinfo);
#line 568 "slarre.f"
	    if (iinfo != 0) {
#line 569 "slarre.f"
		*info = -1;
#line 570 "slarre.f"
		return 0;
#line 571 "slarre.f"
	    }
/* Computing MAX */
#line 572 "slarre.f"
	    d__2 = gl, d__3 = tmp - tmp1 - eps * 100. * (d__1 = tmp - tmp1, 
		    abs(d__1));
#line 572 "slarre.f"
	    isleft = max(d__2,d__3);
#line 575 "slarre.f"
	    slarrk_(&in, &in, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &
		    rtl, &tmp, &tmp1, &iinfo);
#line 577 "slarre.f"
	    if (iinfo != 0) {
#line 578 "slarre.f"
		*info = -1;
#line 579 "slarre.f"
		return 0;
#line 580 "slarre.f"
	    }
/* Computing MIN */
#line 581 "slarre.f"
	    d__2 = gu, d__3 = tmp + tmp1 + eps * 100. * (d__1 = tmp + tmp1, 
		    abs(d__1));
#line 581 "slarre.f"
	    isrght = min(d__2,d__3);
/*           Improve the estimate of the spectral diameter */
#line 584 "slarre.f"
	    spdiam = isrght - isleft;
#line 585 "slarre.f"
	} else {
/*           Case of bisection */
/*           Find approximations to the wanted extremal eigenvalues */
/* Computing MAX */
#line 588 "slarre.f"
	    d__2 = gl, d__3 = w[wbegin] - werr[wbegin] - eps * 100. * (d__1 = 
		    w[wbegin] - werr[wbegin], abs(d__1));
#line 588 "slarre.f"
	    isleft = max(d__2,d__3);
/* Computing MIN */
#line 590 "slarre.f"
	    d__2 = gu, d__3 = w[wend] + werr[wend] + eps * 100. * (d__1 = w[
		    wend] + werr[wend], abs(d__1));
#line 590 "slarre.f"
	    isrght = min(d__2,d__3);
#line 592 "slarre.f"
	}
/*        Decide whether the base representation for the current block */
/*        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I */
/*        should be on the left or the right end of the current block. */
/*        The strategy is to shift to the end which is "more populated" */
/*        Furthermore, decide whether to use DQDS for the computation of */
/*        the eigenvalue approximations at the end of SLARRE or bisection. */
/*        dqds is chosen if all eigenvalues are desired or the number of */
/*        eigenvalues to be computed is large compared to the blocksize. */
#line 603 "slarre.f"
	if (irange == 1 && ! forceb) {
/*           If all the eigenvalues have to be computed, we use dqd */
#line 605 "slarre.f"
	    usedqd = TRUE_;
/*           INDL is the local index of the first eigenvalue to compute */
#line 607 "slarre.f"
	    indl = 1;
#line 608 "slarre.f"
	    indu = in;
/*           MB =  number of eigenvalues to compute */
#line 610 "slarre.f"
	    mb = in;
#line 611 "slarre.f"
	    wend = wbegin + mb - 1;
/*           Define 1/4 and 3/4 points of the spectrum */
#line 613 "slarre.f"
	    s1 = isleft + spdiam * .25;
#line 614 "slarre.f"
	    s2 = isrght - spdiam * .25;
#line 615 "slarre.f"
	} else {
/*           SLARRD has computed IBLOCK and INDEXW for each eigenvalue */
/*           approximation. */
/*           choose sigma */
#line 619 "slarre.f"
	    if (usedqd) {
#line 620 "slarre.f"
		s1 = isleft + spdiam * .25;
#line 621 "slarre.f"
		s2 = isrght - spdiam * .25;
#line 622 "slarre.f"
	    } else {
#line 623 "slarre.f"
		tmp = min(isrght,*vu) - max(isleft,*vl);
#line 624 "slarre.f"
		s1 = max(isleft,*vl) + tmp * .25;
#line 625 "slarre.f"
		s2 = min(isrght,*vu) - tmp * .25;
#line 626 "slarre.f"
	    }
#line 627 "slarre.f"
	}
/*        Compute the negcount at the 1/4 and 3/4 points */
#line 630 "slarre.f"
	if (mb > 1) {
#line 631 "slarre.f"
	    slarrc_("T", &in, &s1, &s2, &d__[ibegin], &e[ibegin], pivmin, &
		    cnt, &cnt1, &cnt2, &iinfo, (ftnlen)1);
#line 633 "slarre.f"
	}
#line 635 "slarre.f"
	if (mb == 1) {
#line 636 "slarre.f"
	    sigma = gl;
#line 637 "slarre.f"
	    sgndef = 1.;
#line 638 "slarre.f"
	} else if (cnt1 - indl >= indu - cnt2) {
#line 639 "slarre.f"
	    if (irange == 1 && ! forceb) {
#line 640 "slarre.f"
		sigma = max(isleft,gl);
#line 641 "slarre.f"
	    } else if (usedqd) {
/*              use Gerschgorin bound as shift to get pos def matrix */
/*              for dqds */
#line 644 "slarre.f"
		sigma = isleft;
#line 645 "slarre.f"
	    } else {
/*              use approximation of the first desired eigenvalue of the */
/*              block as shift */
#line 648 "slarre.f"
		sigma = max(isleft,*vl);
#line 649 "slarre.f"
	    }
#line 650 "slarre.f"
	    sgndef = 1.;
#line 651 "slarre.f"
	} else {
#line 652 "slarre.f"
	    if (irange == 1 && ! forceb) {
#line 653 "slarre.f"
		sigma = min(isrght,gu);
#line 654 "slarre.f"
	    } else if (usedqd) {
/*              use Gerschgorin bound as shift to get neg def matrix */
/*              for dqds */
#line 657 "slarre.f"
		sigma = isrght;
#line 658 "slarre.f"
	    } else {
/*              use approximation of the first desired eigenvalue of the */
/*              block as shift */
#line 661 "slarre.f"
		sigma = min(isrght,*vu);
#line 662 "slarre.f"
	    }
#line 663 "slarre.f"
	    sgndef = -1.;
#line 664 "slarre.f"
	}
/*        An initial SIGMA has been chosen that will be used for computing */
/*        T - SIGMA I = L D L^T */
/*        Define the increment TAU of the shift in case the initial shift */
/*        needs to be refined to obtain a factorization with not too much */
/*        element growth. */
#line 672 "slarre.f"
	if (usedqd) {
/*           The initial SIGMA was to the outer end of the spectrum */
/*           the matrix is definite and we need not retreat. */
#line 675 "slarre.f"
	    tau = spdiam * eps * *n + *pivmin * 2.;
/* Computing MAX */
#line 676 "slarre.f"
	    d__1 = tau, d__2 = eps * 2. * abs(sigma);
#line 676 "slarre.f"
	    tau = max(d__1,d__2);
#line 677 "slarre.f"
	} else {
#line 678 "slarre.f"
	    if (mb > 1) {
#line 679 "slarre.f"
		clwdth = w[wend] + werr[wend] - w[wbegin] - werr[wbegin];
#line 680 "slarre.f"
		avgap = (d__1 = clwdth / (doublereal) (wend - wbegin), abs(
			d__1));
#line 681 "slarre.f"
		if (sgndef == 1.) {
/* Computing MAX */
#line 682 "slarre.f"
		    d__1 = wgap[wbegin];
#line 682 "slarre.f"
		    tau = max(d__1,avgap) * .5;
/* Computing MAX */
#line 683 "slarre.f"
		    d__1 = tau, d__2 = werr[wbegin];
#line 683 "slarre.f"
		    tau = max(d__1,d__2);
#line 684 "slarre.f"
		} else {
/* Computing MAX */
#line 685 "slarre.f"
		    d__1 = wgap[wend - 1];
#line 685 "slarre.f"
		    tau = max(d__1,avgap) * .5;
/* Computing MAX */
#line 686 "slarre.f"
		    d__1 = tau, d__2 = werr[wend];
#line 686 "slarre.f"
		    tau = max(d__1,d__2);
#line 687 "slarre.f"
		}
#line 688 "slarre.f"
	    } else {
#line 689 "slarre.f"
		tau = werr[wbegin];
#line 690 "slarre.f"
	    }
#line 691 "slarre.f"
	}

#line 693 "slarre.f"
	for (idum = 1; idum <= 6; ++idum) {
/*           Compute L D L^T factorization of tridiagonal matrix T - sigma I. */
/*           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of */
/*           pivots in WORK(2*IN+1:3*IN) */
#line 697 "slarre.f"
	    dpivot = d__[ibegin] - sigma;
#line 698 "slarre.f"
	    work[1] = dpivot;
#line 699 "slarre.f"
	    dmax__ = abs(work[1]);
#line 700 "slarre.f"
	    j = ibegin;
#line 701 "slarre.f"
	    i__2 = in - 1;
#line 701 "slarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 702 "slarre.f"
		work[(in << 1) + i__] = 1. / work[i__];
#line 703 "slarre.f"
		tmp = e[j] * work[(in << 1) + i__];
#line 704 "slarre.f"
		work[in + i__] = tmp;
#line 705 "slarre.f"
		dpivot = d__[j + 1] - sigma - tmp * e[j];
#line 706 "slarre.f"
		work[i__ + 1] = dpivot;
/* Computing MAX */
#line 707 "slarre.f"
		d__1 = dmax__, d__2 = abs(dpivot);
#line 707 "slarre.f"
		dmax__ = max(d__1,d__2);
#line 708 "slarre.f"
		++j;
#line 709 "slarre.f"
/* L70: */
#line 709 "slarre.f"
	    }
/*           check for element growth */
#line 711 "slarre.f"
	    if (dmax__ > spdiam * 64.) {
#line 712 "slarre.f"
		norep = TRUE_;
#line 713 "slarre.f"
	    } else {
#line 714 "slarre.f"
		norep = FALSE_;
#line 715 "slarre.f"
	    }
#line 716 "slarre.f"
	    if (usedqd && ! norep) {
/*              Ensure the definiteness of the representation */
/*              All entries of D (of L D L^T) must have the same sign */
#line 719 "slarre.f"
		i__2 = in;
#line 719 "slarre.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 720 "slarre.f"
		    tmp = sgndef * work[i__];
#line 721 "slarre.f"
		    if (tmp < 0.) {
#line 721 "slarre.f"
			norep = TRUE_;
#line 721 "slarre.f"
		    }
#line 722 "slarre.f"
/* L71: */
#line 722 "slarre.f"
		}
#line 723 "slarre.f"
	    }
#line 724 "slarre.f"
	    if (norep) {
/*              Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin */
/*              shift which makes the matrix definite. So we should end up */
/*              here really only in the case of IRANGE = VALRNG or INDRNG. */
#line 728 "slarre.f"
		if (idum == 5) {
#line 729 "slarre.f"
		    if (sgndef == 1.) {
/*                    The fudged Gerschgorin shift should succeed */
#line 731 "slarre.f"
			sigma = gl - spdiam * 2. * eps * *n - *pivmin * 4.;
#line 733 "slarre.f"
		    } else {
#line 734 "slarre.f"
			sigma = gu + spdiam * 2. * eps * *n + *pivmin * 4.;
#line 736 "slarre.f"
		    }
#line 737 "slarre.f"
		} else {
#line 738 "slarre.f"
		    sigma -= sgndef * tau;
#line 739 "slarre.f"
		    tau *= 2.;
#line 740 "slarre.f"
		}
#line 741 "slarre.f"
	    } else {
/*              an initial RRR is found */
#line 743 "slarre.f"
		goto L83;
#line 744 "slarre.f"
	    }
#line 745 "slarre.f"
/* L80: */
#line 745 "slarre.f"
	}
/*        if the program reaches this point, no base representation could be */
/*        found in MAXTRY iterations. */
#line 748 "slarre.f"
	*info = 2;
#line 749 "slarre.f"
	return 0;
#line 751 "slarre.f"
L83:
/*        At this point, we have found an initial base representation */
/*        T - SIGMA I = L D L^T with not too much element growth. */
/*        Store the shift. */
#line 755 "slarre.f"
	e[iend] = sigma;
/*        Store D and L. */
#line 757 "slarre.f"
	scopy_(&in, &work[1], &c__1, &d__[ibegin], &c__1);
#line 758 "slarre.f"
	i__2 = in - 1;
#line 758 "slarre.f"
	scopy_(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
#line 761 "slarre.f"
	if (mb > 1) {

/*           Perturb each entry of the base representation by a small */
/*           (but random) relative amount to overcome difficulties with */
/*           glued matrices. */

#line 767 "slarre.f"
	    for (i__ = 1; i__ <= 4; ++i__) {
#line 768 "slarre.f"
		iseed[i__ - 1] = 1;
#line 769 "slarre.f"
/* L122: */
#line 769 "slarre.f"
	    }
#line 771 "slarre.f"
	    i__2 = (in << 1) - 1;
#line 771 "slarre.f"
	    slarnv_(&c__2, iseed, &i__2, &work[1]);
#line 772 "slarre.f"
	    i__2 = in - 1;
#line 772 "slarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 773 "slarre.f"
		d__[ibegin + i__ - 1] *= eps * 4. * work[i__] + 1.;
#line 774 "slarre.f"
		e[ibegin + i__ - 1] *= eps * 4. * work[in + i__] + 1.;
#line 775 "slarre.f"
/* L125: */
#line 775 "slarre.f"
	    }
#line 776 "slarre.f"
	    d__[iend] *= eps * 4. * work[in] + 1.;

#line 778 "slarre.f"
	}

/*        Don't update the Gerschgorin intervals because keeping track */
/*        of the updates would be too much work in SLARRV. */
/*        We update W instead and use it to locate the proper Gerschgorin */
/*        intervals. */
/*        Compute the required eigenvalues of L D L' by bisection or dqds */
#line 786 "slarre.f"
	if (! usedqd) {
/*           If SLARRD has been used, shift the eigenvalue approximations */
/*           according to their representation. This is necessary for */
/*           a uniform SLARRV since dqds computes eigenvalues of the */
/*           shifted representation. In SLARRV, W will always hold the */
/*           UNshifted eigenvalue approximation. */
#line 792 "slarre.f"
	    i__2 = wend;
#line 792 "slarre.f"
	    for (j = wbegin; j <= i__2; ++j) {
#line 793 "slarre.f"
		w[j] -= sigma;
#line 794 "slarre.f"
		werr[j] += (d__1 = w[j], abs(d__1)) * eps;
#line 795 "slarre.f"
/* L134: */
#line 795 "slarre.f"
	    }
/*           call SLARRB to reduce eigenvalue error of the approximations */
/*           from SLARRD */
#line 798 "slarre.f"
	    i__2 = iend - 1;
#line 798 "slarre.f"
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
/* Computing 2nd power */
#line 799 "slarre.f"
		d__1 = e[i__];
#line 799 "slarre.f"
		work[i__] = d__[i__] * (d__1 * d__1);
#line 800 "slarre.f"
/* L135: */
#line 800 "slarre.f"
	    }
/*           use bisection to find EV from INDL to INDU */
#line 802 "slarre.f"
	    i__2 = indl - 1;
#line 802 "slarre.f"
	    slarrb_(&in, &d__[ibegin], &work[ibegin], &indl, &indu, rtol1, 
		    rtol2, &i__2, &w[wbegin], &wgap[wbegin], &werr[wbegin], &
		    work[(*n << 1) + 1], &iwork[1], pivmin, &spdiam, &in, &
		    iinfo);
#line 807 "slarre.f"
	    if (iinfo != 0) {
#line 808 "slarre.f"
		*info = -4;
#line 809 "slarre.f"
		return 0;
#line 810 "slarre.f"
	    }
/*           SLARRB computes all gaps correctly except for the last one */
/*           Record distance to VU/GU */
/* Computing MAX */
#line 813 "slarre.f"
	    d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
#line 813 "slarre.f"
	    wgap[wend] = max(d__1,d__2);
#line 815 "slarre.f"
	    i__2 = indu;
#line 815 "slarre.f"
	    for (i__ = indl; i__ <= i__2; ++i__) {
#line 816 "slarre.f"
		++(*m);
#line 817 "slarre.f"
		iblock[*m] = jblk;
#line 818 "slarre.f"
		indexw[*m] = i__;
#line 819 "slarre.f"
/* L138: */
#line 819 "slarre.f"
	    }
#line 820 "slarre.f"
	} else {
/*           Call dqds to get all eigs (and then possibly delete unwanted */
/*           eigenvalues). */
/*           Note that dqds finds the eigenvalues of the L D L^T representation */
/*           of T to high relative accuracy. High relative accuracy */
/*           might be lost when the shift of the RRR is subtracted to obtain */
/*           the eigenvalues of T. However, T is not guaranteed to define its */
/*           eigenvalues to high relative accuracy anyway. */
/*           Set RTOL to the order of the tolerance used in SLASQ2 */
/*           This is an ESTIMATED error, the worst case bound is 4*N*EPS */
/*           which is usually too large and requires unnecessary work to be */
/*           done by bisection when computing the eigenvectors */
#line 832 "slarre.f"
	    rtol = log((doublereal) in) * 4. * eps;
#line 833 "slarre.f"
	    j = ibegin;
#line 834 "slarre.f"
	    i__2 = in - 1;
#line 834 "slarre.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 835 "slarre.f"
		work[(i__ << 1) - 1] = (d__1 = d__[j], abs(d__1));
#line 836 "slarre.f"
		work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
#line 837 "slarre.f"
		++j;
#line 838 "slarre.f"
/* L140: */
#line 838 "slarre.f"
	    }
#line 839 "slarre.f"
	    work[(in << 1) - 1] = (d__1 = d__[iend], abs(d__1));
#line 840 "slarre.f"
	    work[in * 2] = 0.;
#line 841 "slarre.f"
	    slasq2_(&in, &work[1], &iinfo);
#line 842 "slarre.f"
	    if (iinfo != 0) {
/*              If IINFO = -5 then an index is part of a tight cluster */
/*              and should be changed. The index is in IWORK(1) and the */
/*              gap is in WORK(N+1) */
#line 846 "slarre.f"
		*info = -5;
#line 847 "slarre.f"
		return 0;
#line 848 "slarre.f"
	    } else {
/*              Test that all eigenvalues are positive as expected */
#line 850 "slarre.f"
		i__2 = in;
#line 850 "slarre.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 851 "slarre.f"
		    if (work[i__] < 0.) {
#line 852 "slarre.f"
			*info = -6;
#line 853 "slarre.f"
			return 0;
#line 854 "slarre.f"
		    }
#line 855 "slarre.f"
/* L149: */
#line 855 "slarre.f"
		}
#line 856 "slarre.f"
	    }
#line 857 "slarre.f"
	    if (sgndef > 0.) {
#line 858 "slarre.f"
		i__2 = indu;
#line 858 "slarre.f"
		for (i__ = indl; i__ <= i__2; ++i__) {
#line 859 "slarre.f"
		    ++(*m);
#line 860 "slarre.f"
		    w[*m] = work[in - i__ + 1];
#line 861 "slarre.f"
		    iblock[*m] = jblk;
#line 862 "slarre.f"
		    indexw[*m] = i__;
#line 863 "slarre.f"
/* L150: */
#line 863 "slarre.f"
		}
#line 864 "slarre.f"
	    } else {
#line 865 "slarre.f"
		i__2 = indu;
#line 865 "slarre.f"
		for (i__ = indl; i__ <= i__2; ++i__) {
#line 866 "slarre.f"
		    ++(*m);
#line 867 "slarre.f"
		    w[*m] = -work[i__];
#line 868 "slarre.f"
		    iblock[*m] = jblk;
#line 869 "slarre.f"
		    indexw[*m] = i__;
#line 870 "slarre.f"
/* L160: */
#line 870 "slarre.f"
		}
#line 871 "slarre.f"
	    }
#line 873 "slarre.f"
	    i__2 = *m;
#line 873 "slarre.f"
	    for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
/*              the value of RTOL below should be the tolerance in SLASQ2 */
#line 875 "slarre.f"
		werr[i__] = rtol * (d__1 = w[i__], abs(d__1));
#line 876 "slarre.f"
/* L165: */
#line 876 "slarre.f"
	    }
#line 877 "slarre.f"
	    i__2 = *m - 1;
#line 877 "slarre.f"
	    for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
/*              compute the right gap between the intervals */
/* Computing MAX */
#line 879 "slarre.f"
		d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + werr[
			i__]);
#line 879 "slarre.f"
		wgap[i__] = max(d__1,d__2);
#line 881 "slarre.f"
/* L166: */
#line 881 "slarre.f"
	    }
/* Computing MAX */
#line 882 "slarre.f"
	    d__1 = 0., d__2 = *vu - sigma - (w[*m] + werr[*m]);
#line 882 "slarre.f"
	    wgap[*m] = max(d__1,d__2);
#line 884 "slarre.f"
	}
/*        proceed with next block */
#line 886 "slarre.f"
	ibegin = iend + 1;
#line 887 "slarre.f"
	wbegin = wend + 1;
#line 888 "slarre.f"
L170:
#line 888 "slarre.f"
	;
#line 888 "slarre.f"
    }

#line 891 "slarre.f"
    return 0;

/*     end of SLARRE */

} /* slarre_ */

