#line 1 "slarrv.f"
/* slarrv.f -- translated by f2c (version 20100827).
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

#line 1 "slarrv.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SLARRV computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenv
alues of L D LT. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARRV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARRV( N, VL, VU, D, L, PIVMIN, */
/*                          ISPLIT, M, DOL, DOU, MINRGP, */
/*                          RTOL1, RTOL2, W, WERR, WGAP, */
/*                          IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            DOL, DOU, INFO, LDZ, M, N */
/*       REAL               MINRGP, PIVMIN, RTOL1, RTOL2, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), INDEXW( * ), ISPLIT( * ), */
/*      $                   ISUPPZ( * ), IWORK( * ) */
/*       REAL               D( * ), GERS( * ), L( * ), W( * ), WERR( * ), */
/*      $                   WGAP( * ), WORK( * ) */
/*       REAL              Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARRV computes the eigenvectors of the tridiagonal matrix */
/* > T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T. */
/* > The input eigenvalues should have been computed by SLARRE. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
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
/* >          Lower and upper bounds of the interval that contains the desired */
/* >          eigenvalues. VL < VU. Needed to compute gaps on the left or right */
/* >          end of the extremal eigenvalues in the desired RANGE. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the N diagonal elements of the diagonal matrix D. */
/* >          On exit, D may be overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in,out] L */
/* > \verbatim */
/* >          L is REAL array, dimension (N) */
/* >          On entry, the (N-1) subdiagonal elements of the unit */
/* >          bidiagonal matrix L are in elements 1 to N-1 of L */
/* >          (if the matrix is not splitted.) At the end of each block */
/* >          is stored the corresponding shift as given by SLARRE. */
/* >          On exit, L is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
/* >          The minimum pivot allowed in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[in] ISPLIT */
/* > \verbatim */
/* >          ISPLIT is INTEGER array, dimension (N) */
/* >          The splitting points, at which T breaks up into blocks. */
/* >          The first block consists of rows/columns 1 to */
/* >          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1 */
/* >          through ISPLIT( 2 ), etc. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The total number of input eigenvalues.  0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] DOL */
/* > \verbatim */
/* >          DOL is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] DOU */
/* > \verbatim */
/* >          DOU is INTEGER */
/* >          If the user wants to compute only selected eigenvectors from all */
/* >          the eigenvalues supplied, he can specify an index range DOL:DOU. */
/* >          Or else the setting DOL=1, DOU=M should be applied. */
/* >          Note that DOL and DOU refer to the order in which the eigenvalues */
/* >          are stored in W. */
/* >          If the user wants to compute only selected eigenpairs, then */
/* >          the columns DOL-1 to DOU+1 of the eigenvector space Z contain the */
/* >          computed eigenvectors. All other columns of Z are set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] MINRGP */
/* > \verbatim */
/* >          MINRGP is REAL */
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
/* > \param[in,out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          The first M elements of W contain the APPROXIMATE eigenvalues for */
/* >          which eigenvectors are to be computed.  The eigenvalues */
/* >          should be grouped by split-off block and ordered from */
/* >          smallest to largest within the block ( The output array */
/* >          W from SLARRE is expected here ). Furthermore, they are with */
/* >          respect to the shift of the corresponding root representation */
/* >          for their block. On exit, W holds the eigenvalues of the */
/* >          UNshifted matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WERR */
/* > \verbatim */
/* >          WERR is REAL array, dimension (N) */
/* >          The first M elements contain the semiwidth of the uncertainty */
/* >          interval of the corresponding eigenvalue in W */
/* > \endverbatim */
/* > */
/* > \param[in,out] WGAP */
/* > \verbatim */
/* >          WGAP is REAL array, dimension (N) */
/* >          The separation from the right neighbor eigenvalue in W. */
/* > \endverbatim */
/* > */
/* > \param[in] IBLOCK */
/* > \verbatim */
/* >          IBLOCK is INTEGER array, dimension (N) */
/* >          The indices of the blocks (submatrices) associated with the */
/* >          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue */
/* >          W(i) belongs to the first block from the top, =2 if W(i) */
/* >          belongs to the second block, etc. */
/* > \endverbatim */
/* > */
/* > \param[in] INDEXW */
/* > \verbatim */
/* >          INDEXW is INTEGER array, dimension (N) */
/* >          The indices of the eigenvalues within each block (submatrix); */
/* >          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the */
/* >          i-th eigenvalue W(i) is the 10-th eigenvalue in the second block. */
/* > \endverbatim */
/* > */
/* > \param[in] GERS */
/* > \verbatim */
/* >          GERS is REAL array, dimension (2*N) */
/* >          The N Gerschgorin intervals (the i-th Gerschgorin interval */
/* >          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should */
/* >          be computed from the original UNshifted matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, max(1,M) ) */
/* >          If INFO = 0, the first M columns of Z contain the */
/* >          orthonormal eigenvectors of the matrix T */
/* >          corresponding to the input eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          Note: the user must ensure that at least max(1,M) columns are */
/* >          supplied in the array Z. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* >          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) ) */
/* >          The support of the eigenvectors in Z, i.e., the indices */
/* >          indicating the nonzero elements in Z. The I-th eigenvector */
/* >          is nonzero only in elements ISUPPZ( 2*I-1 ) through */
/* >          ISUPPZ( 2*I ). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (12*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (7*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* > */
/* >          > 0:  A problem occured in SLARRV. */
/* >          < 0:  One of the called subroutines signaled an internal problem. */
/* >                Needs inspection of the corresponding parameter IINFO */
/* >                for further information. */
/* > */
/* >          =-1:  Problem in SLARRB when refining a child's eigenvalues. */
/* >          =-2:  Problem in SLARRF when computing the RRR of a child. */
/* >                When a child is inside a tight cluster, it can be difficult */
/* >                to find an RRR. A partial remedy from the user's point of */
/* >                view is to make the parameter MINRGP smaller and recompile. */
/* >                However, as the orthogonality of the computed vectors is */
/* >                proportional to 1/MINRGP, the user should be aware that */
/* >                he might be trading in precision when he decreases MINRGP. */
/* >          =-3:  Problem in SLARRB when refining a single eigenvalue */
/* >                after the Rayleigh correction was rejected. */
/* >          = 5:  The Rayleigh Quotient Iteration failed to converge to */
/* >                full accuracy in MAXITR steps. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup realOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int slarrv_(integer *n, doublereal *vl, doublereal *vu, 
	doublereal *d__, doublereal *l, doublereal *pivmin, integer *isplit, 
	integer *m, integer *dol, integer *dou, doublereal *minrgp, 
	doublereal *rtol1, doublereal *rtol2, doublereal *w, doublereal *werr,
	 doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers,
	 doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *iwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    logical L__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer minwsize, i__, j, k, p, q, miniwsize, ii;
    static doublereal gl;
    static integer im, in;
    static doublereal gu, gap, eps, tau, tol, tmp;
    static integer zto;
    static doublereal ztz;
    static integer iend, jblk;
    static doublereal lgap;
    static integer done;
    static doublereal rgap, left;
    static integer wend, iter;
    static doublereal bstw;
    static integer itmp1, indld;
    static doublereal fudge;
    static integer idone;
    static doublereal sigma;
    static integer iinfo, iindr;
    static doublereal resid;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical eskip;
    static doublereal right;
    static integer nclus, zfrom;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal rqtol;
    static integer iindc1, iindc2;
    extern /* Subroutine */ int slar1v_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static logical stp2ii;
    static doublereal lambda;
    static integer ibegin, indeig;
    static logical needbs;
    static integer indlld;
    static doublereal sgndef, mingma;
    extern doublereal slamch_(char *, ftnlen);
    static integer oldien, oldncl, wbegin;
    static doublereal spdiam;
    static integer negcnt, oldcls;
    static doublereal savgap;
    static integer ndepth;
    static doublereal ssigma;
    static logical usedbs;
    static integer iindwk, offset;
    static doublereal gaptol;
    extern /* Subroutine */ int slarrb_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *), slarrf_(
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer newcls, oldfst, indwrk, windex, oldlst;
    static logical usedrq;
    static integer newfst, newftt, parity, windmn, isupmn, newlst, windpl, 
	    zusedl, newsiz, zusedu, zusedw;
    static doublereal bstres, nrminv;
    static logical tryrqc;
    static integer isupmx;
    static doublereal rqcorr;
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
#line 343 "slarrv.f"
    /* Parameter adjustments */
#line 343 "slarrv.f"
    --d__;
#line 343 "slarrv.f"
    --l;
#line 343 "slarrv.f"
    --isplit;
#line 343 "slarrv.f"
    --w;
#line 343 "slarrv.f"
    --werr;
#line 343 "slarrv.f"
    --wgap;
#line 343 "slarrv.f"
    --iblock;
#line 343 "slarrv.f"
    --indexw;
#line 343 "slarrv.f"
    --gers;
#line 343 "slarrv.f"
    z_dim1 = *ldz;
#line 343 "slarrv.f"
    z_offset = 1 + z_dim1;
#line 343 "slarrv.f"
    z__ -= z_offset;
#line 343 "slarrv.f"
    --isuppz;
#line 343 "slarrv.f"
    --work;
#line 343 "slarrv.f"
    --iwork;
#line 343 "slarrv.f"

#line 343 "slarrv.f"
    /* Function Body */
#line 343 "slarrv.f"
    *info = 0;
/*     The first N entries of WORK are reserved for the eigenvalues */
#line 345 "slarrv.f"
    indld = *n + 1;
#line 346 "slarrv.f"
    indlld = (*n << 1) + 1;
#line 347 "slarrv.f"
    indwrk = *n * 3 + 1;
#line 348 "slarrv.f"
    minwsize = *n * 12;
#line 350 "slarrv.f"
    i__1 = minwsize;
#line 350 "slarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 351 "slarrv.f"
	work[i__] = 0.;
#line 352 "slarrv.f"
/* L5: */
#line 352 "slarrv.f"
    }
/*     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the */
/*     factorization used to compute the FP vector */
#line 356 "slarrv.f"
    iindr = 0;
/*     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current */
/*     layer and the one above. */
#line 359 "slarrv.f"
    iindc1 = *n;
#line 360 "slarrv.f"
    iindc2 = *n << 1;
#line 361 "slarrv.f"
    iindwk = *n * 3 + 1;
#line 363 "slarrv.f"
    miniwsize = *n * 7;
#line 364 "slarrv.f"
    i__1 = miniwsize;
#line 364 "slarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 365 "slarrv.f"
	iwork[i__] = 0;
#line 366 "slarrv.f"
/* L10: */
#line 366 "slarrv.f"
    }
#line 368 "slarrv.f"
    zusedl = 1;
#line 369 "slarrv.f"
    if (*dol > 1) {
/*        Set lower bound for use of Z */
#line 371 "slarrv.f"
	zusedl = *dol - 1;
#line 372 "slarrv.f"
    }
#line 373 "slarrv.f"
    zusedu = *m;
#line 374 "slarrv.f"
    if (*dou < *m) {
/*        Set lower bound for use of Z */
#line 376 "slarrv.f"
	zusedu = *dou + 1;
#line 377 "slarrv.f"
    }
/*     The width of the part of Z that is used */
#line 379 "slarrv.f"
    zusedw = zusedu - zusedl + 1;
#line 382 "slarrv.f"
    slaset_("Full", n, &zusedw, &c_b5, &c_b5, &z__[zusedl * z_dim1 + 1], ldz, 
	    (ftnlen)4);
#line 385 "slarrv.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 386 "slarrv.f"
    rqtol = eps * 2.;

/*     Set expert flags for standard code. */
#line 389 "slarrv.f"
    tryrqc = TRUE_;
#line 391 "slarrv.f"
    if (*dol == 1 && *dou == *m) {
#line 392 "slarrv.f"
    } else {
/*        Only selected eigenpairs are computed. Since the other evalues */
/*        are not refined by RQ iteration, bisection has to compute to full */
/*        accuracy. */
#line 396 "slarrv.f"
	*rtol1 = eps * 4.;
#line 397 "slarrv.f"
	*rtol2 = eps * 4.;
#line 398 "slarrv.f"
    }
/*     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the */
/*     desired eigenvalues. The support of the nonzero eigenvector */
/*     entries is contained in the interval IBEGIN:IEND. */
/*     Remark that if k eigenpairs are desired, then the eigenvectors */
/*     are stored in k contiguous columns of Z. */
/*     DONE is the number of eigenvectors already computed */
#line 407 "slarrv.f"
    done = 0;
#line 408 "slarrv.f"
    ibegin = 1;
#line 409 "slarrv.f"
    wbegin = 1;
#line 410 "slarrv.f"
    i__1 = iblock[*m];
#line 410 "slarrv.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 411 "slarrv.f"
	iend = isplit[jblk];
#line 412 "slarrv.f"
	sigma = l[iend];
/*        Find the eigenvectors of the submatrix indexed IBEGIN */
/*        through IEND. */
#line 415 "slarrv.f"
	wend = wbegin - 1;
#line 416 "slarrv.f"
L15:
#line 417 "slarrv.f"
	if (wend < *m) {
#line 418 "slarrv.f"
	    if (iblock[wend + 1] == jblk) {
#line 419 "slarrv.f"
		++wend;
#line 420 "slarrv.f"
		goto L15;
#line 421 "slarrv.f"
	    }
#line 422 "slarrv.f"
	}
#line 423 "slarrv.f"
	if (wend < wbegin) {
#line 424 "slarrv.f"
	    ibegin = iend + 1;
#line 425 "slarrv.f"
	    goto L170;
#line 426 "slarrv.f"
	} else if (wend < *dol || wbegin > *dou) {
#line 427 "slarrv.f"
	    ibegin = iend + 1;
#line 428 "slarrv.f"
	    wbegin = wend + 1;
#line 429 "slarrv.f"
	    goto L170;
#line 430 "slarrv.f"
	}
/*        Find local spectral diameter of the block */
#line 433 "slarrv.f"
	gl = gers[(ibegin << 1) - 1];
#line 434 "slarrv.f"
	gu = gers[ibegin * 2];
#line 435 "slarrv.f"
	i__2 = iend;
#line 435 "slarrv.f"
	for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 436 "slarrv.f"
	    d__1 = gers[(i__ << 1) - 1];
#line 436 "slarrv.f"
	    gl = min(d__1,gl);
/* Computing MAX */
#line 437 "slarrv.f"
	    d__1 = gers[i__ * 2];
#line 437 "slarrv.f"
	    gu = max(d__1,gu);
#line 438 "slarrv.f"
/* L20: */
#line 438 "slarrv.f"
	}
#line 439 "slarrv.f"
	spdiam = gu - gl;
/*        OLDIEN is the last index of the previous block */
#line 442 "slarrv.f"
	oldien = ibegin - 1;
/*        Calculate the size of the current block */
#line 444 "slarrv.f"
	in = iend - ibegin + 1;
/*        The number of eigenvalues in the current block */
#line 446 "slarrv.f"
	im = wend - wbegin + 1;
/*        This is for a 1x1 block */
#line 449 "slarrv.f"
	if (ibegin == iend) {
#line 450 "slarrv.f"
	    ++done;
#line 451 "slarrv.f"
	    z__[ibegin + wbegin * z_dim1] = 1.;
#line 452 "slarrv.f"
	    isuppz[(wbegin << 1) - 1] = ibegin;
#line 453 "slarrv.f"
	    isuppz[wbegin * 2] = ibegin;
#line 454 "slarrv.f"
	    w[wbegin] += sigma;
#line 455 "slarrv.f"
	    work[wbegin] = w[wbegin];
#line 456 "slarrv.f"
	    ibegin = iend + 1;
#line 457 "slarrv.f"
	    ++wbegin;
#line 458 "slarrv.f"
	    goto L170;
#line 459 "slarrv.f"
	}
/*        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND) */
/*        Note that these can be approximations, in this case, the corresp. */
/*        entries of WERR give the size of the uncertainty interval. */
/*        The eigenvalue approximations will be refined when necessary as */
/*        high relative accuracy is required for the computation of the */
/*        corresponding eigenvectors. */
#line 467 "slarrv.f"
	scopy_(&im, &w[wbegin], &c__1, &work[wbegin], &c__1);
/*        We store in W the eigenvalue approximations w.r.t. the original */
/*        matrix T. */
#line 472 "slarrv.f"
	i__2 = im;
#line 472 "slarrv.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 473 "slarrv.f"
	    w[wbegin + i__ - 1] += sigma;
#line 474 "slarrv.f"
/* L30: */
#line 474 "slarrv.f"
	}
/*        NDEPTH is the current depth of the representation tree */
#line 478 "slarrv.f"
	ndepth = 0;
/*        PARITY is either 1 or 0 */
#line 480 "slarrv.f"
	parity = 1;
/*        NCLUS is the number of clusters for the next level of the */
/*        representation tree, we start with NCLUS = 1 for the root */
#line 483 "slarrv.f"
	nclus = 1;
#line 484 "slarrv.f"
	iwork[iindc1 + 1] = 1;
#line 485 "slarrv.f"
	iwork[iindc1 + 2] = im;
/*        IDONE is the number of eigenvectors already computed in the current */
/*        block */
#line 489 "slarrv.f"
	idone = 0;
/*        loop while( IDONE.LT.IM ) */
/*        generate the representation tree for the current block and */
/*        compute the eigenvectors */
#line 493 "slarrv.f"
L40:
#line 494 "slarrv.f"
	if (idone < im) {
/*           This is a crude protection against infinitely deep trees */
#line 496 "slarrv.f"
	    if (ndepth > *m) {
#line 497 "slarrv.f"
		*info = -2;
#line 498 "slarrv.f"
		return 0;
#line 499 "slarrv.f"
	    }
/*           breadth first processing of the current level of the representation */
/*           tree: OLDNCL = number of clusters on current level */
#line 502 "slarrv.f"
	    oldncl = nclus;
/*           reset NCLUS to count the number of child clusters */
#line 504 "slarrv.f"
	    nclus = 0;

#line 506 "slarrv.f"
	    parity = 1 - parity;
#line 507 "slarrv.f"
	    if (parity == 0) {
#line 508 "slarrv.f"
		oldcls = iindc1;
#line 509 "slarrv.f"
		newcls = iindc2;
#line 510 "slarrv.f"
	    } else {
#line 511 "slarrv.f"
		oldcls = iindc2;
#line 512 "slarrv.f"
		newcls = iindc1;
#line 513 "slarrv.f"
	    }
/*           Process the clusters on the current level */
#line 515 "slarrv.f"
	    i__2 = oldncl;
#line 515 "slarrv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 516 "slarrv.f"
		j = oldcls + (i__ << 1);
/*              OLDFST, OLDLST = first, last index of current cluster. */
/*                               cluster indices start with 1 and are relative */
/*                               to WBEGIN when accessing W, WGAP, WERR, Z */
#line 520 "slarrv.f"
		oldfst = iwork[j - 1];
#line 521 "slarrv.f"
		oldlst = iwork[j];
#line 522 "slarrv.f"
		if (ndepth > 0) {
/*                 Retrieve relatively robust representation (RRR) of cluster */
/*                 that has been computed at the previous level */
/*                 The RRR is stored in Z and overwritten once the eigenvectors */
/*                 have been computed or when the cluster is refined */
#line 528 "slarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Get representation from location of the leftmost evalue */
/*                    of the cluster */
#line 531 "slarrv.f"
			j = wbegin + oldfst - 1;
#line 532 "slarrv.f"
		    } else {
#line 533 "slarrv.f"
			if (wbegin + oldfst - 1 < *dol) {
/*                       Get representation from the left end of Z array */
#line 535 "slarrv.f"
			    j = *dol - 1;
#line 536 "slarrv.f"
			} else if (wbegin + oldfst - 1 > *dou) {
/*                       Get representation from the right end of Z array */
#line 538 "slarrv.f"
			    j = *dou;
#line 539 "slarrv.f"
			} else {
#line 540 "slarrv.f"
			    j = wbegin + oldfst - 1;
#line 541 "slarrv.f"
			}
#line 542 "slarrv.f"
		    }
#line 543 "slarrv.f"
		    scopy_(&in, &z__[ibegin + j * z_dim1], &c__1, &d__[ibegin]
			    , &c__1);
#line 544 "slarrv.f"
		    i__3 = in - 1;
#line 544 "slarrv.f"
		    scopy_(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1, &l[
			    ibegin], &c__1);
#line 546 "slarrv.f"
		    sigma = z__[iend + (j + 1) * z_dim1];
/*                 Set the corresponding entries in Z to zero */
#line 549 "slarrv.f"
		    slaset_("Full", &in, &c__2, &c_b5, &c_b5, &z__[ibegin + j 
			    * z_dim1], ldz, (ftnlen)4);
#line 551 "slarrv.f"
		}
/*              Compute DL and DLL of current RRR */
#line 554 "slarrv.f"
		i__3 = iend - 1;
#line 554 "slarrv.f"
		for (j = ibegin; j <= i__3; ++j) {
#line 555 "slarrv.f"
		    tmp = d__[j] * l[j];
#line 556 "slarrv.f"
		    work[indld - 1 + j] = tmp;
#line 557 "slarrv.f"
		    work[indlld - 1 + j] = tmp * l[j];
#line 558 "slarrv.f"
/* L50: */
#line 558 "slarrv.f"
		}
#line 560 "slarrv.f"
		if (ndepth > 0) {
/*                 P and Q are index of the first and last eigenvalue to compute */
/*                 within the current block */
#line 563 "slarrv.f"
		    p = indexw[wbegin - 1 + oldfst];
#line 564 "slarrv.f"
		    q = indexw[wbegin - 1 + oldlst];
/*                 Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET */
/*                 through the Q-OFFSET elements of these arrays are to be used. */
/*                  OFFSET = P-OLDFST */
#line 568 "slarrv.f"
		    offset = indexw[wbegin] - 1;
/*                 perform limited bisection (if necessary) to get approximate */
/*                 eigenvalues to the precision needed. */
#line 571 "slarrv.f"
		    slarrb_(&in, &d__[ibegin], &work[indlld + ibegin - 1], &p,
			     &q, rtol1, rtol2, &offset, &work[wbegin], &wgap[
			    wbegin], &werr[wbegin], &work[indwrk], &iwork[
			    iindwk], pivmin, &spdiam, &in, &iinfo);
#line 577 "slarrv.f"
		    if (iinfo != 0) {
#line 578 "slarrv.f"
			*info = -1;
#line 579 "slarrv.f"
			return 0;
#line 580 "slarrv.f"
		    }
/*                 We also recompute the extremal gaps. W holds all eigenvalues */
/*                 of the unshifted matrix and must be used for computation */
/*                 of WGAP, the entries of WORK might stem from RRRs with */
/*                 different shifts. The gaps from WBEGIN-1+OLDFST to */
/*                 WBEGIN-1+OLDLST are correctly computed in SLARRB. */
/*                 However, we only allow the gaps to become greater since */
/*                 this is what should happen when we decrease WERR */
#line 588 "slarrv.f"
		    if (oldfst > 1) {
/* Computing MAX */
#line 589 "slarrv.f"
			d__1 = wgap[wbegin + oldfst - 2], d__2 = w[wbegin + 
				oldfst - 1] - werr[wbegin + oldfst - 1] - w[
				wbegin + oldfst - 2] - werr[wbegin + oldfst - 
				2];
#line 589 "slarrv.f"
			wgap[wbegin + oldfst - 2] = max(d__1,d__2);
#line 593 "slarrv.f"
		    }
#line 594 "slarrv.f"
		    if (wbegin + oldlst - 1 < wend) {
/* Computing MAX */
#line 595 "slarrv.f"
			d__1 = wgap[wbegin + oldlst - 1], d__2 = w[wbegin + 
				oldlst] - werr[wbegin + oldlst] - w[wbegin + 
				oldlst - 1] - werr[wbegin + oldlst - 1];
#line 595 "slarrv.f"
			wgap[wbegin + oldlst - 1] = max(d__1,d__2);
#line 599 "slarrv.f"
		    }
/*                 Each time the eigenvalues in WORK get refined, we store */
/*                 the newly found approximation with all shifts applied in W */
#line 602 "slarrv.f"
		    i__3 = oldlst;
#line 602 "slarrv.f"
		    for (j = oldfst; j <= i__3; ++j) {
#line 603 "slarrv.f"
			w[wbegin + j - 1] = work[wbegin + j - 1] + sigma;
#line 604 "slarrv.f"
/* L53: */
#line 604 "slarrv.f"
		    }
#line 605 "slarrv.f"
		}
/*              Process the current node. */
#line 608 "slarrv.f"
		newfst = oldfst;
#line 609 "slarrv.f"
		i__3 = oldlst;
#line 609 "slarrv.f"
		for (j = oldfst; j <= i__3; ++j) {
#line 610 "slarrv.f"
		    if (j == oldlst) {
/*                    we are at the right end of the cluster, this is also the */
/*                    boundary of the child cluster */
#line 613 "slarrv.f"
			newlst = j;
#line 614 "slarrv.f"
		    } else if (wgap[wbegin + j - 1] >= *minrgp * (d__1 = work[
			    wbegin + j - 1], abs(d__1))) {
/*                    the right relative gap is big enough, the child cluster */
/*                    (NEWFST,..,NEWLST) is well separated from the following */
#line 618 "slarrv.f"
			newlst = j;
#line 619 "slarrv.f"
		    } else {
/*                    inside a child cluster, the relative gap is not */
/*                    big enough. */
#line 622 "slarrv.f"
			goto L140;
#line 623 "slarrv.f"
		    }
/*                 Compute size of child cluster found */
#line 626 "slarrv.f"
		    newsiz = newlst - newfst + 1;
/*                 NEWFTT is the place in Z where the new RRR or the computed */
/*                 eigenvector is to be stored */
#line 630 "slarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Store representation at location of the leftmost evalue */
/*                    of the cluster */
#line 633 "slarrv.f"
			newftt = wbegin + newfst - 1;
#line 634 "slarrv.f"
		    } else {
#line 635 "slarrv.f"
			if (wbegin + newfst - 1 < *dol) {
/*                       Store representation at the left end of Z array */
#line 637 "slarrv.f"
			    newftt = *dol - 1;
#line 638 "slarrv.f"
			} else if (wbegin + newfst - 1 > *dou) {
/*                       Store representation at the right end of Z array */
#line 640 "slarrv.f"
			    newftt = *dou;
#line 641 "slarrv.f"
			} else {
#line 642 "slarrv.f"
			    newftt = wbegin + newfst - 1;
#line 643 "slarrv.f"
			}
#line 644 "slarrv.f"
		    }
#line 646 "slarrv.f"
		    if (newsiz > 1) {

/*                    Current child is not a singleton but a cluster. */
/*                    Compute and store new representation of child. */


/*                    Compute left and right cluster gap. */

/*                    LGAP and RGAP are not computed from WORK because */
/*                    the eigenvalue approximations may stem from RRRs */
/*                    different shifts. However, W hold all eigenvalues */
/*                    of the unshifted matrix. Still, the entries in WGAP */
/*                    have to be computed from WORK since the entries */
/*                    in W might be of the same order so that gaps are not */
/*                    exhibited correctly for very close eigenvalues. */
#line 661 "slarrv.f"
			if (newfst == 1) {
/* Computing MAX */
#line 662 "slarrv.f"
			    d__1 = 0., d__2 = w[wbegin] - werr[wbegin] - *vl;
#line 662 "slarrv.f"
			    lgap = max(d__1,d__2);
#line 664 "slarrv.f"
			} else {
#line 665 "slarrv.f"
			    lgap = wgap[wbegin + newfst - 2];
#line 666 "slarrv.f"
			}
#line 667 "slarrv.f"
			rgap = wgap[wbegin + newlst - 1];

/*                    Compute left- and rightmost eigenvalue of child */
/*                    to high precision in order to shift as close */
/*                    as possible and obtain as large relative gaps */
/*                    as possible */

#line 674 "slarrv.f"
			for (k = 1; k <= 2; ++k) {
#line 675 "slarrv.f"
			    if (k == 1) {
#line 676 "slarrv.f"
				p = indexw[wbegin - 1 + newfst];
#line 677 "slarrv.f"
			    } else {
#line 678 "slarrv.f"
				p = indexw[wbegin - 1 + newlst];
#line 679 "slarrv.f"
			    }
#line 680 "slarrv.f"
			    offset = indexw[wbegin] - 1;
#line 681 "slarrv.f"
			    slarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &p, &p, &rqtol, &rqtol, &offset, &
				    work[wbegin], &wgap[wbegin], &werr[wbegin]
				    , &work[indwrk], &iwork[iindwk], pivmin, &
				    spdiam, &in, &iinfo);
#line 688 "slarrv.f"
/* L55: */
#line 688 "slarrv.f"
			}

#line 690 "slarrv.f"
			if (wbegin + newlst - 1 < *dol || wbegin + newfst - 1 
				> *dou) {
/*                       if the cluster contains no desired eigenvalues */
/*                       skip the computation of that branch of the rep. tree */

/*                       We could skip before the refinement of the extremal */
/*                       eigenvalues of the child, but then the representation */
/*                       tree could be different from the one when nothing is */
/*                       skipped. For this reason we skip at this place. */
#line 699 "slarrv.f"
			    idone = idone + newlst - newfst + 1;
#line 700 "slarrv.f"
			    goto L139;
#line 701 "slarrv.f"
			}

/*                    Compute RRR of child cluster. */
/*                    Note that the new RRR is stored in Z */

/*                    SLARRF needs LWORK = 2*N */
#line 707 "slarrv.f"
			slarrf_(&in, &d__[ibegin], &l[ibegin], &work[indld + 
				ibegin - 1], &newfst, &newlst, &work[wbegin], 
				&wgap[wbegin], &werr[wbegin], &spdiam, &lgap, 
				&rgap, pivmin, &tau, &z__[ibegin + newftt * 
				z_dim1], &z__[ibegin + (newftt + 1) * z_dim1],
				 &work[indwrk], &iinfo);
#line 714 "slarrv.f"
			if (iinfo == 0) {
/*                       a new RRR for the cluster was found by SLARRF */
/*                       update shift and store it */
#line 717 "slarrv.f"
			    ssigma = sigma + tau;
#line 718 "slarrv.f"
			    z__[iend + (newftt + 1) * z_dim1] = ssigma;
/*                       WORK() are the midpoints and WERR() the semi-width */
/*                       Note that the entries in W are unchanged. */
#line 721 "slarrv.f"
			    i__4 = newlst;
#line 721 "slarrv.f"
			    for (k = newfst; k <= i__4; ++k) {
#line 722 "slarrv.f"
				fudge = eps * 3. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
#line 724 "slarrv.f"
				work[wbegin + k - 1] -= tau;
#line 726 "slarrv.f"
				fudge += eps * 4. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
/*                          Fudge errors */
#line 729 "slarrv.f"
				werr[wbegin + k - 1] += fudge;
/*                          Gaps are not fudged. Provided that WERR is small */
/*                          when eigenvalues are close, a zero gap indicates */
/*                          that a new representation is needed for resolving */
/*                          the cluster. A fudge could lead to a wrong decision */
/*                          of judging eigenvalues 'separated' which in */
/*                          reality are not. This could have a negative impact */
/*                          on the orthogonality of the computed eigenvectors. */
#line 738 "slarrv.f"
/* L116: */
#line 738 "slarrv.f"
			    }
#line 740 "slarrv.f"
			    ++nclus;
#line 741 "slarrv.f"
			    k = newcls + (nclus << 1);
#line 742 "slarrv.f"
			    iwork[k - 1] = newfst;
#line 743 "slarrv.f"
			    iwork[k] = newlst;
#line 744 "slarrv.f"
			} else {
#line 745 "slarrv.f"
			    *info = -2;
#line 746 "slarrv.f"
			    return 0;
#line 747 "slarrv.f"
			}
#line 748 "slarrv.f"
		    } else {

/*                    Compute eigenvector of singleton */

#line 752 "slarrv.f"
			iter = 0;

#line 754 "slarrv.f"
			tol = log((doublereal) in) * 4. * eps;

#line 756 "slarrv.f"
			k = newfst;
#line 757 "slarrv.f"
			windex = wbegin + k - 1;
/* Computing MAX */
#line 758 "slarrv.f"
			i__4 = windex - 1;
#line 758 "slarrv.f"
			windmn = max(i__4,1);
/* Computing MIN */
#line 759 "slarrv.f"
			i__4 = windex + 1;
#line 759 "slarrv.f"
			windpl = min(i__4,*m);
#line 760 "slarrv.f"
			lambda = work[windex];
#line 761 "slarrv.f"
			++done;
/*                    Check if eigenvector computation is to be skipped */
#line 763 "slarrv.f"
			if (windex < *dol || windex > *dou) {
#line 765 "slarrv.f"
			    eskip = TRUE_;
#line 766 "slarrv.f"
			    goto L125;
#line 767 "slarrv.f"
			} else {
#line 768 "slarrv.f"
			    eskip = FALSE_;
#line 769 "slarrv.f"
			}
#line 770 "slarrv.f"
			left = work[windex] - werr[windex];
#line 771 "slarrv.f"
			right = work[windex] + werr[windex];
#line 772 "slarrv.f"
			indeig = indexw[windex];
/*                    Note that since we compute the eigenpairs for a child, */
/*                    all eigenvalue approximations are w.r.t the same shift. */
/*                    In this case, the entries in WORK should be used for */
/*                    computing the gaps since they exhibit even very small */
/*                    differences in the eigenvalues, as opposed to the */
/*                    entries in W which might "look" the same. */
#line 780 "slarrv.f"
			if (k == 1) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VL, the formula */
/*                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA ) */
/*                       can lead to an overestimation of the left gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small left gap. */
/* Computing MAX */
#line 787 "slarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 787 "slarrv.f"
			    lgap = eps * max(d__1,d__2);
#line 788 "slarrv.f"
			} else {
#line 789 "slarrv.f"
			    lgap = wgap[windmn];
#line 790 "slarrv.f"
			}
#line 791 "slarrv.f"
			if (k == im) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VU, the formula */
/*                       can lead to an overestimation of the right gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small right gap. */
/* Computing MAX */
#line 797 "slarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 797 "slarrv.f"
			    rgap = eps * max(d__1,d__2);
#line 798 "slarrv.f"
			} else {
#line 799 "slarrv.f"
			    rgap = wgap[windex];
#line 800 "slarrv.f"
			}
#line 801 "slarrv.f"
			gap = min(lgap,rgap);
#line 802 "slarrv.f"
			if (k == 1 || k == im) {
/*                       The eigenvector support can become wrong */
/*                       because significant entries could be cut off due to a */
/*                       large GAPTOL parameter in LAR1V. Prevent this. */
#line 806 "slarrv.f"
			    gaptol = 0.;
#line 807 "slarrv.f"
			} else {
#line 808 "slarrv.f"
			    gaptol = gap * eps;
#line 809 "slarrv.f"
			}
#line 810 "slarrv.f"
			isupmn = in;
#line 811 "slarrv.f"
			isupmx = 1;
/*                    Update WGAP so that it holds the minimum gap */
/*                    to the left or the right. This is crucial in the */
/*                    case where bisection is used to ensure that the */
/*                    eigenvalue is refined up to the required precision. */
/*                    The correct value is restored afterwards. */
#line 817 "slarrv.f"
			savgap = wgap[windex];
#line 818 "slarrv.f"
			wgap[windex] = gap;
/*                    We want to use the Rayleigh Quotient Correction */
/*                    as often as possible since it converges quadratically */
/*                    when we are close enough to the desired eigenvalue. */
/*                    However, the Rayleigh Quotient can have the wrong sign */
/*                    and lead us away from the desired eigenvalue. In this */
/*                    case, the best we can do is to use bisection. */
#line 825 "slarrv.f"
			usedbs = FALSE_;
#line 826 "slarrv.f"
			usedrq = FALSE_;
/*                    Bisection is initially turned off unless it is forced */
#line 828 "slarrv.f"
			needbs = ! tryrqc;
#line 829 "slarrv.f"
L120:
/*                    Check if bisection should be used to refine eigenvalue */
#line 831 "slarrv.f"
			if (needbs) {
/*                       Take the bisection as new iterate */
#line 833 "slarrv.f"
			    usedbs = TRUE_;
#line 834 "slarrv.f"
			    itmp1 = iwork[iindr + windex];
#line 835 "slarrv.f"
			    offset = indexw[wbegin] - 1;
#line 836 "slarrv.f"
			    d__1 = eps * 2.;
#line 836 "slarrv.f"
			    slarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &indeig, &indeig, &c_b5, &d__1, &
				    offset, &work[wbegin], &wgap[wbegin], &
				    werr[wbegin], &work[indwrk], &iwork[
				    iindwk], pivmin, &spdiam, &itmp1, &iinfo);
#line 843 "slarrv.f"
			    if (iinfo != 0) {
#line 844 "slarrv.f"
				*info = -3;
#line 845 "slarrv.f"
				return 0;
#line 846 "slarrv.f"
			    }
#line 847 "slarrv.f"
			    lambda = work[windex];
/*                       Reset twist index from inaccurate LAMBDA to */
/*                       force computation of true MINGMA */
#line 850 "slarrv.f"
			    iwork[iindr + windex] = 0;
#line 851 "slarrv.f"
			}
/*                    Given LAMBDA, compute the eigenvector. */
#line 853 "slarrv.f"
			L__1 = ! usedbs;
#line 853 "slarrv.f"
			slar1v_(&in, &c__1, &in, &lambda, &d__[ibegin], &l[
				ibegin], &work[indld + ibegin - 1], &work[
				indlld + ibegin - 1], pivmin, &gaptol, &z__[
				ibegin + windex * z_dim1], &L__1, &negcnt, &
				ztz, &mingma, &iwork[iindr + windex], &isuppz[
				(windex << 1) - 1], &nrminv, &resid, &rqcorr, 
				&work[indwrk]);
#line 860 "slarrv.f"
			if (iter == 0) {
#line 861 "slarrv.f"
			    bstres = resid;
#line 862 "slarrv.f"
			    bstw = lambda;
#line 863 "slarrv.f"
			} else if (resid < bstres) {
#line 864 "slarrv.f"
			    bstres = resid;
#line 865 "slarrv.f"
			    bstw = lambda;
#line 866 "slarrv.f"
			}
/* Computing MIN */
#line 867 "slarrv.f"
			i__4 = isupmn, i__5 = isuppz[(windex << 1) - 1];
#line 867 "slarrv.f"
			isupmn = min(i__4,i__5);
/* Computing MAX */
#line 868 "slarrv.f"
			i__4 = isupmx, i__5 = isuppz[windex * 2];
#line 868 "slarrv.f"
			isupmx = max(i__4,i__5);
#line 869 "slarrv.f"
			++iter;
/*                    sin alpha <= |resid|/gap */
/*                    Note that both the residual and the gap are */
/*                    proportional to the matrix, so ||T|| doesn't play */
/*                    a role in the quotient */

/*                    Convergence test for Rayleigh-Quotient iteration */
/*                    (omitted when Bisection has been used) */

#line 880 "slarrv.f"
			if (resid > tol * gap && abs(rqcorr) > rqtol * abs(
				lambda) && ! usedbs) {
/*                       We need to check that the RQCORR update doesn't */
/*                       move the eigenvalue away from the desired one and */
/*                       towards a neighbor. -> protection with bisection */
#line 886 "slarrv.f"
			    if (indeig <= negcnt) {
/*                          The wanted eigenvalue lies to the left */
#line 888 "slarrv.f"
				sgndef = -1.;
#line 889 "slarrv.f"
			    } else {
/*                          The wanted eigenvalue lies to the right */
#line 891 "slarrv.f"
				sgndef = 1.;
#line 892 "slarrv.f"
			    }
/*                       We only use the RQCORR if it improves the */
/*                       the iterate reasonably. */
#line 895 "slarrv.f"
			    if (rqcorr * sgndef >= 0. && lambda + rqcorr <= 
				    right && lambda + rqcorr >= left) {
#line 899 "slarrv.f"
				usedrq = TRUE_;
/*                          Store new midpoint of bisection interval in WORK */
#line 901 "slarrv.f"
				if (sgndef == 1.) {
/*                             The current LAMBDA is on the left of the true */
/*                             eigenvalue */
#line 904 "slarrv.f"
				    left = lambda;
/*                             We prefer to assume that the error estimate */
/*                             is correct. We could make the interval not */
/*                             as a bracket but to be modified if the RQCORR */
/*                             chooses to. In this case, the RIGHT side should */
/*                             be modified as follows: */
/*                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR) */
#line 911 "slarrv.f"
				} else {
/*                             The current LAMBDA is on the right of the true */
/*                             eigenvalue */
#line 914 "slarrv.f"
				    right = lambda;
/*                             See comment about assuming the error estimate is */
/*                             correct above. */
/*                              LEFT = MIN(LEFT, LAMBDA + RQCORR) */
#line 918 "slarrv.f"
				}
#line 919 "slarrv.f"
				work[windex] = (right + left) * .5;
/*                          Take RQCORR since it has the correct sign and */
/*                          improves the iterate reasonably */
#line 923 "slarrv.f"
				lambda += rqcorr;
/*                          Update width of error interval */
#line 925 "slarrv.f"
				werr[windex] = (right - left) * .5;
#line 927 "slarrv.f"
			    } else {
#line 928 "slarrv.f"
				needbs = TRUE_;
#line 929 "slarrv.f"
			    }
#line 930 "slarrv.f"
			    if (right - left < rqtol * abs(lambda)) {
/*                             The eigenvalue is computed to bisection accuracy */
/*                             compute eigenvector and stop */
#line 933 "slarrv.f"
				usedbs = TRUE_;
#line 934 "slarrv.f"
				goto L120;
#line 935 "slarrv.f"
			    } else if (iter < 10) {
#line 936 "slarrv.f"
				goto L120;
#line 937 "slarrv.f"
			    } else if (iter == 10) {
#line 938 "slarrv.f"
				needbs = TRUE_;
#line 939 "slarrv.f"
				goto L120;
#line 940 "slarrv.f"
			    } else {
#line 941 "slarrv.f"
				*info = 5;
#line 942 "slarrv.f"
				return 0;
#line 943 "slarrv.f"
			    }
#line 944 "slarrv.f"
			} else {
#line 945 "slarrv.f"
			    stp2ii = FALSE_;
#line 946 "slarrv.f"
			    if (usedrq && usedbs && bstres <= resid) {
#line 948 "slarrv.f"
				lambda = bstw;
#line 949 "slarrv.f"
				stp2ii = TRUE_;
#line 950 "slarrv.f"
			    }
#line 951 "slarrv.f"
			    if (stp2ii) {
/*                          improve error angle by second step */
#line 953 "slarrv.f"
				L__1 = ! usedbs;
#line 953 "slarrv.f"
				slar1v_(&in, &c__1, &in, &lambda, &d__[ibegin]
					, &l[ibegin], &work[indld + ibegin - 
					1], &work[indlld + ibegin - 1], 
					pivmin, &gaptol, &z__[ibegin + windex 
					* z_dim1], &L__1, &negcnt, &ztz, &
					mingma, &iwork[iindr + windex], &
					isuppz[(windex << 1) - 1], &nrminv, &
					resid, &rqcorr, &work[indwrk]);
#line 962 "slarrv.f"
			    }
#line 963 "slarrv.f"
			    work[windex] = lambda;
#line 964 "slarrv.f"
			}

/*                    Compute FP-vector support w.r.t. whole matrix */

#line 968 "slarrv.f"
			isuppz[(windex << 1) - 1] += oldien;
#line 969 "slarrv.f"
			isuppz[windex * 2] += oldien;
#line 970 "slarrv.f"
			zfrom = isuppz[(windex << 1) - 1];
#line 971 "slarrv.f"
			zto = isuppz[windex * 2];
#line 972 "slarrv.f"
			isupmn += oldien;
#line 973 "slarrv.f"
			isupmx += oldien;
/*                    Ensure vector is ok if support in the RQI has changed */
#line 975 "slarrv.f"
			if (isupmn < zfrom) {
#line 976 "slarrv.f"
			    i__4 = zfrom - 1;
#line 976 "slarrv.f"
			    for (ii = isupmn; ii <= i__4; ++ii) {
#line 977 "slarrv.f"
				z__[ii + windex * z_dim1] = 0.;
#line 978 "slarrv.f"
/* L122: */
#line 978 "slarrv.f"
			    }
#line 979 "slarrv.f"
			}
#line 980 "slarrv.f"
			if (isupmx > zto) {
#line 981 "slarrv.f"
			    i__4 = isupmx;
#line 981 "slarrv.f"
			    for (ii = zto + 1; ii <= i__4; ++ii) {
#line 982 "slarrv.f"
				z__[ii + windex * z_dim1] = 0.;
#line 983 "slarrv.f"
/* L123: */
#line 983 "slarrv.f"
			    }
#line 984 "slarrv.f"
			}
#line 985 "slarrv.f"
			i__4 = zto - zfrom + 1;
#line 985 "slarrv.f"
			sscal_(&i__4, &nrminv, &z__[zfrom + windex * z_dim1], 
				&c__1);
#line 987 "slarrv.f"
L125:
/*                    Update W */
#line 989 "slarrv.f"
			w[windex] = lambda + sigma;
/*                    Recompute the gaps on the left and right */
/*                    But only allow them to become larger and not */
/*                    smaller (which can only happen through "bad" */
/*                    cancellation and doesn't reflect the theory */
/*                    where the initial gaps are underestimated due */
/*                    to WERR being too crude.) */
#line 996 "slarrv.f"
			if (! eskip) {
#line 997 "slarrv.f"
			    if (k > 1) {
/* Computing MAX */
#line 998 "slarrv.f"
				d__1 = wgap[windmn], d__2 = w[windex] - werr[
					windex] - w[windmn] - werr[windmn];
#line 998 "slarrv.f"
				wgap[windmn] = max(d__1,d__2);
#line 1001 "slarrv.f"
			    }
#line 1002 "slarrv.f"
			    if (windex < wend) {
/* Computing MAX */
#line 1003 "slarrv.f"
				d__1 = savgap, d__2 = w[windpl] - werr[windpl]
					 - w[windex] - werr[windex];
#line 1003 "slarrv.f"
				wgap[windex] = max(d__1,d__2);
#line 1006 "slarrv.f"
			    }
#line 1007 "slarrv.f"
			}
#line 1008 "slarrv.f"
			++idone;
#line 1009 "slarrv.f"
		    }
/*                 here ends the code for the current child */

#line 1012 "slarrv.f"
L139:
/*                 Proceed to any remaining child nodes */
#line 1014 "slarrv.f"
		    newfst = j + 1;
#line 1015 "slarrv.f"
L140:
#line 1015 "slarrv.f"
		    ;
#line 1015 "slarrv.f"
		}
#line 1016 "slarrv.f"
/* L150: */
#line 1016 "slarrv.f"
	    }
#line 1017 "slarrv.f"
	    ++ndepth;
#line 1018 "slarrv.f"
	    goto L40;
#line 1019 "slarrv.f"
	}
#line 1020 "slarrv.f"
	ibegin = iend + 1;
#line 1021 "slarrv.f"
	wbegin = wend + 1;
#line 1022 "slarrv.f"
L170:
#line 1022 "slarrv.f"
	;
#line 1022 "slarrv.f"
    }

#line 1025 "slarrv.f"
    return 0;

/*     End of SLARRV */

} /* slarrv_ */

