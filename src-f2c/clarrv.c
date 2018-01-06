#line 1 "clarrv.f"
/* clarrv.f -- translated by f2c (version 20100827).
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

#line 1 "clarrv.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b28 = 0.;

/* > \brief \b CLARRV computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenv
alues of L D LT. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARRV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarrv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarrv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarrv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARRV( N, VL, VU, D, L, PIVMIN, */
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
/*       COMPLEX           Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARRV computes the eigenvectors of the tridiagonal matrix */
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
/* >          Z is array, dimension (LDZ, max(1,M) ) */
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
/* >          > 0:  A problem occured in CLARRV. */
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

/* > \ingroup complexOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int clarrv_(integer *n, doublereal *vl, doublereal *vu, 
	doublereal *d__, doublereal *l, doublereal *pivmin, integer *isplit, 
	integer *m, integer *dol, integer *dou, doublereal *minrgp, 
	doublereal *rtol1, doublereal *rtol2, doublereal *w, doublereal *werr,
	 doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers,
	 doublecomplex *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *iwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;
    doublecomplex z__1;
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
    static logical eskip;
    static doublereal right;
    static integer nclus, zfrom;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal rqtol;
    static integer iindc1, iindc2, indin1, indin2;
    extern /* Subroutine */ int clar1v_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublecomplex *, 
	    logical *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static logical stp2ii;
    static doublereal lambda;
    static integer ibegin, indeig;
    static logical needbs;
    static integer indlld;
    static doublereal sgndef, mingma;
    extern doublereal slamch_(char *, ftnlen);
    static integer oldien, oldncl, wbegin;
    static doublereal spdiam;
    static integer negcnt;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static integer oldcls;
    static doublereal savgap;
    static integer ndepth;
    static doublereal ssigma;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static logical usedbs;
    static integer iindwk, offset;
    static doublereal gaptol;
    extern /* Subroutine */ int slarrb_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *);
    static integer newcls, oldfst, indwrk, windex, oldlst;
    static logical usedrq;
    static integer newfst, newftt, parity, windmn, windpl, isupmn, newlst, 
	    zusedl;
    static doublereal bstres;
    static integer newsiz, zusedu, zusedw;
    static doublereal nrminv, rqcorr;
    static logical tryrqc;
    static integer isupmx;
    extern /* Subroutine */ int slarrf_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);


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
#line 347 "clarrv.f"
    /* Parameter adjustments */
#line 347 "clarrv.f"
    --d__;
#line 347 "clarrv.f"
    --l;
#line 347 "clarrv.f"
    --isplit;
#line 347 "clarrv.f"
    --w;
#line 347 "clarrv.f"
    --werr;
#line 347 "clarrv.f"
    --wgap;
#line 347 "clarrv.f"
    --iblock;
#line 347 "clarrv.f"
    --indexw;
#line 347 "clarrv.f"
    --gers;
#line 347 "clarrv.f"
    z_dim1 = *ldz;
#line 347 "clarrv.f"
    z_offset = 1 + z_dim1;
#line 347 "clarrv.f"
    z__ -= z_offset;
#line 347 "clarrv.f"
    --isuppz;
#line 347 "clarrv.f"
    --work;
#line 347 "clarrv.f"
    --iwork;
#line 347 "clarrv.f"

#line 347 "clarrv.f"
    /* Function Body */
#line 347 "clarrv.f"
    *info = 0;
/*     The first N entries of WORK are reserved for the eigenvalues */
#line 349 "clarrv.f"
    indld = *n + 1;
#line 350 "clarrv.f"
    indlld = (*n << 1) + 1;
#line 351 "clarrv.f"
    indin1 = *n * 3 + 1;
#line 352 "clarrv.f"
    indin2 = (*n << 2) + 1;
#line 353 "clarrv.f"
    indwrk = *n * 5 + 1;
#line 354 "clarrv.f"
    minwsize = *n * 12;
#line 356 "clarrv.f"
    i__1 = minwsize;
#line 356 "clarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 357 "clarrv.f"
	work[i__] = 0.;
#line 358 "clarrv.f"
/* L5: */
#line 358 "clarrv.f"
    }
/*     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the */
/*     factorization used to compute the FP vector */
#line 362 "clarrv.f"
    iindr = 0;
/*     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current */
/*     layer and the one above. */
#line 365 "clarrv.f"
    iindc1 = *n;
#line 366 "clarrv.f"
    iindc2 = *n << 1;
#line 367 "clarrv.f"
    iindwk = *n * 3 + 1;
#line 369 "clarrv.f"
    miniwsize = *n * 7;
#line 370 "clarrv.f"
    i__1 = miniwsize;
#line 370 "clarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 371 "clarrv.f"
	iwork[i__] = 0;
#line 372 "clarrv.f"
/* L10: */
#line 372 "clarrv.f"
    }
#line 374 "clarrv.f"
    zusedl = 1;
#line 375 "clarrv.f"
    if (*dol > 1) {
/*        Set lower bound for use of Z */
#line 377 "clarrv.f"
	zusedl = *dol - 1;
#line 378 "clarrv.f"
    }
#line 379 "clarrv.f"
    zusedu = *m;
#line 380 "clarrv.f"
    if (*dou < *m) {
/*        Set lower bound for use of Z */
#line 382 "clarrv.f"
	zusedu = *dou + 1;
#line 383 "clarrv.f"
    }
/*     The width of the part of Z that is used */
#line 385 "clarrv.f"
    zusedw = zusedu - zusedl + 1;
#line 388 "clarrv.f"
    claset_("Full", n, &zusedw, &c_b1, &c_b1, &z__[zusedl * z_dim1 + 1], ldz, 
	    (ftnlen)4);
#line 391 "clarrv.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 392 "clarrv.f"
    rqtol = eps * 2.;

/*     Set expert flags for standard code. */
#line 395 "clarrv.f"
    tryrqc = TRUE_;
#line 397 "clarrv.f"
    if (*dol == 1 && *dou == *m) {
#line 398 "clarrv.f"
    } else {
/*        Only selected eigenpairs are computed. Since the other evalues */
/*        are not refined by RQ iteration, bisection has to compute to full */
/*        accuracy. */
#line 402 "clarrv.f"
	*rtol1 = eps * 4.;
#line 403 "clarrv.f"
	*rtol2 = eps * 4.;
#line 404 "clarrv.f"
    }
/*     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the */
/*     desired eigenvalues. The support of the nonzero eigenvector */
/*     entries is contained in the interval IBEGIN:IEND. */
/*     Remark that if k eigenpairs are desired, then the eigenvectors */
/*     are stored in k contiguous columns of Z. */
/*     DONE is the number of eigenvectors already computed */
#line 413 "clarrv.f"
    done = 0;
#line 414 "clarrv.f"
    ibegin = 1;
#line 415 "clarrv.f"
    wbegin = 1;
#line 416 "clarrv.f"
    i__1 = iblock[*m];
#line 416 "clarrv.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 417 "clarrv.f"
	iend = isplit[jblk];
#line 418 "clarrv.f"
	sigma = l[iend];
/*        Find the eigenvectors of the submatrix indexed IBEGIN */
/*        through IEND. */
#line 421 "clarrv.f"
	wend = wbegin - 1;
#line 422 "clarrv.f"
L15:
#line 423 "clarrv.f"
	if (wend < *m) {
#line 424 "clarrv.f"
	    if (iblock[wend + 1] == jblk) {
#line 425 "clarrv.f"
		++wend;
#line 426 "clarrv.f"
		goto L15;
#line 427 "clarrv.f"
	    }
#line 428 "clarrv.f"
	}
#line 429 "clarrv.f"
	if (wend < wbegin) {
#line 430 "clarrv.f"
	    ibegin = iend + 1;
#line 431 "clarrv.f"
	    goto L170;
#line 432 "clarrv.f"
	} else if (wend < *dol || wbegin > *dou) {
#line 433 "clarrv.f"
	    ibegin = iend + 1;
#line 434 "clarrv.f"
	    wbegin = wend + 1;
#line 435 "clarrv.f"
	    goto L170;
#line 436 "clarrv.f"
	}
/*        Find local spectral diameter of the block */
#line 439 "clarrv.f"
	gl = gers[(ibegin << 1) - 1];
#line 440 "clarrv.f"
	gu = gers[ibegin * 2];
#line 441 "clarrv.f"
	i__2 = iend;
#line 441 "clarrv.f"
	for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 442 "clarrv.f"
	    d__1 = gers[(i__ << 1) - 1];
#line 442 "clarrv.f"
	    gl = min(d__1,gl);
/* Computing MAX */
#line 443 "clarrv.f"
	    d__1 = gers[i__ * 2];
#line 443 "clarrv.f"
	    gu = max(d__1,gu);
#line 444 "clarrv.f"
/* L20: */
#line 444 "clarrv.f"
	}
#line 445 "clarrv.f"
	spdiam = gu - gl;
/*        OLDIEN is the last index of the previous block */
#line 448 "clarrv.f"
	oldien = ibegin - 1;
/*        Calculate the size of the current block */
#line 450 "clarrv.f"
	in = iend - ibegin + 1;
/*        The number of eigenvalues in the current block */
#line 452 "clarrv.f"
	im = wend - wbegin + 1;
/*        This is for a 1x1 block */
#line 455 "clarrv.f"
	if (ibegin == iend) {
#line 456 "clarrv.f"
	    ++done;
#line 457 "clarrv.f"
	    i__2 = ibegin + wbegin * z_dim1;
#line 457 "clarrv.f"
	    z__[i__2].r = 1., z__[i__2].i = 0.;
#line 458 "clarrv.f"
	    isuppz[(wbegin << 1) - 1] = ibegin;
#line 459 "clarrv.f"
	    isuppz[wbegin * 2] = ibegin;
#line 460 "clarrv.f"
	    w[wbegin] += sigma;
#line 461 "clarrv.f"
	    work[wbegin] = w[wbegin];
#line 462 "clarrv.f"
	    ibegin = iend + 1;
#line 463 "clarrv.f"
	    ++wbegin;
#line 464 "clarrv.f"
	    goto L170;
#line 465 "clarrv.f"
	}
/*        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND) */
/*        Note that these can be approximations, in this case, the corresp. */
/*        entries of WERR give the size of the uncertainty interval. */
/*        The eigenvalue approximations will be refined when necessary as */
/*        high relative accuracy is required for the computation of the */
/*        corresponding eigenvectors. */
#line 473 "clarrv.f"
	scopy_(&im, &w[wbegin], &c__1, &work[wbegin], &c__1);
/*        We store in W the eigenvalue approximations w.r.t. the original */
/*        matrix T. */
#line 478 "clarrv.f"
	i__2 = im;
#line 478 "clarrv.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 479 "clarrv.f"
	    w[wbegin + i__ - 1] += sigma;
#line 480 "clarrv.f"
/* L30: */
#line 480 "clarrv.f"
	}
/*        NDEPTH is the current depth of the representation tree */
#line 484 "clarrv.f"
	ndepth = 0;
/*        PARITY is either 1 or 0 */
#line 486 "clarrv.f"
	parity = 1;
/*        NCLUS is the number of clusters for the next level of the */
/*        representation tree, we start with NCLUS = 1 for the root */
#line 489 "clarrv.f"
	nclus = 1;
#line 490 "clarrv.f"
	iwork[iindc1 + 1] = 1;
#line 491 "clarrv.f"
	iwork[iindc1 + 2] = im;
/*        IDONE is the number of eigenvectors already computed in the current */
/*        block */
#line 495 "clarrv.f"
	idone = 0;
/*        loop while( IDONE.LT.IM ) */
/*        generate the representation tree for the current block and */
/*        compute the eigenvectors */
#line 499 "clarrv.f"
L40:
#line 500 "clarrv.f"
	if (idone < im) {
/*           This is a crude protection against infinitely deep trees */
#line 502 "clarrv.f"
	    if (ndepth > *m) {
#line 503 "clarrv.f"
		*info = -2;
#line 504 "clarrv.f"
		return 0;
#line 505 "clarrv.f"
	    }
/*           breadth first processing of the current level of the representation */
/*           tree: OLDNCL = number of clusters on current level */
#line 508 "clarrv.f"
	    oldncl = nclus;
/*           reset NCLUS to count the number of child clusters */
#line 510 "clarrv.f"
	    nclus = 0;

#line 512 "clarrv.f"
	    parity = 1 - parity;
#line 513 "clarrv.f"
	    if (parity == 0) {
#line 514 "clarrv.f"
		oldcls = iindc1;
#line 515 "clarrv.f"
		newcls = iindc2;
#line 516 "clarrv.f"
	    } else {
#line 517 "clarrv.f"
		oldcls = iindc2;
#line 518 "clarrv.f"
		newcls = iindc1;
#line 519 "clarrv.f"
	    }
/*           Process the clusters on the current level */
#line 521 "clarrv.f"
	    i__2 = oldncl;
#line 521 "clarrv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 522 "clarrv.f"
		j = oldcls + (i__ << 1);
/*              OLDFST, OLDLST = first, last index of current cluster. */
/*                               cluster indices start with 1 and are relative */
/*                               to WBEGIN when accessing W, WGAP, WERR, Z */
#line 526 "clarrv.f"
		oldfst = iwork[j - 1];
#line 527 "clarrv.f"
		oldlst = iwork[j];
#line 528 "clarrv.f"
		if (ndepth > 0) {
/*                 Retrieve relatively robust representation (RRR) of cluster */
/*                 that has been computed at the previous level */
/*                 The RRR is stored in Z and overwritten once the eigenvectors */
/*                 have been computed or when the cluster is refined */
#line 534 "clarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Get representation from location of the leftmost evalue */
/*                    of the cluster */
#line 537 "clarrv.f"
			j = wbegin + oldfst - 1;
#line 538 "clarrv.f"
		    } else {
#line 539 "clarrv.f"
			if (wbegin + oldfst - 1 < *dol) {
/*                       Get representation from the left end of Z array */
#line 541 "clarrv.f"
			    j = *dol - 1;
#line 542 "clarrv.f"
			} else if (wbegin + oldfst - 1 > *dou) {
/*                       Get representation from the right end of Z array */
#line 544 "clarrv.f"
			    j = *dou;
#line 545 "clarrv.f"
			} else {
#line 546 "clarrv.f"
			    j = wbegin + oldfst - 1;
#line 547 "clarrv.f"
			}
#line 548 "clarrv.f"
		    }
#line 549 "clarrv.f"
		    i__3 = in - 1;
#line 549 "clarrv.f"
		    for (k = 1; k <= i__3; ++k) {
#line 550 "clarrv.f"
			i__4 = ibegin + k - 1 + j * z_dim1;
#line 550 "clarrv.f"
			d__[ibegin + k - 1] = z__[i__4].r;
#line 552 "clarrv.f"
			i__4 = ibegin + k - 1 + (j + 1) * z_dim1;
#line 552 "clarrv.f"
			l[ibegin + k - 1] = z__[i__4].r;
#line 554 "clarrv.f"
/* L45: */
#line 554 "clarrv.f"
		    }
#line 555 "clarrv.f"
		    i__3 = iend + j * z_dim1;
#line 555 "clarrv.f"
		    d__[iend] = z__[i__3].r;
#line 556 "clarrv.f"
		    i__3 = iend + (j + 1) * z_dim1;
#line 556 "clarrv.f"
		    sigma = z__[i__3].r;
/*                 Set the corresponding entries in Z to zero */
#line 559 "clarrv.f"
		    claset_("Full", &in, &c__2, &c_b1, &c_b1, &z__[ibegin + j 
			    * z_dim1], ldz, (ftnlen)4);
#line 561 "clarrv.f"
		}
/*              Compute DL and DLL of current RRR */
#line 564 "clarrv.f"
		i__3 = iend - 1;
#line 564 "clarrv.f"
		for (j = ibegin; j <= i__3; ++j) {
#line 565 "clarrv.f"
		    tmp = d__[j] * l[j];
#line 566 "clarrv.f"
		    work[indld - 1 + j] = tmp;
#line 567 "clarrv.f"
		    work[indlld - 1 + j] = tmp * l[j];
#line 568 "clarrv.f"
/* L50: */
#line 568 "clarrv.f"
		}
#line 570 "clarrv.f"
		if (ndepth > 0) {
/*                 P and Q are index of the first and last eigenvalue to compute */
/*                 within the current block */
#line 573 "clarrv.f"
		    p = indexw[wbegin - 1 + oldfst];
#line 574 "clarrv.f"
		    q = indexw[wbegin - 1 + oldlst];
/*                 Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET */
/*                 through the Q-OFFSET elements of these arrays are to be used. */
/*                  OFFSET = P-OLDFST */
#line 578 "clarrv.f"
		    offset = indexw[wbegin] - 1;
/*                 perform limited bisection (if necessary) to get approximate */
/*                 eigenvalues to the precision needed. */
#line 581 "clarrv.f"
		    slarrb_(&in, &d__[ibegin], &work[indlld + ibegin - 1], &p,
			     &q, rtol1, rtol2, &offset, &work[wbegin], &wgap[
			    wbegin], &werr[wbegin], &work[indwrk], &iwork[
			    iindwk], pivmin, &spdiam, &in, &iinfo);
#line 587 "clarrv.f"
		    if (iinfo != 0) {
#line 588 "clarrv.f"
			*info = -1;
#line 589 "clarrv.f"
			return 0;
#line 590 "clarrv.f"
		    }
/*                 We also recompute the extremal gaps. W holds all eigenvalues */
/*                 of the unshifted matrix and must be used for computation */
/*                 of WGAP, the entries of WORK might stem from RRRs with */
/*                 different shifts. The gaps from WBEGIN-1+OLDFST to */
/*                 WBEGIN-1+OLDLST are correctly computed in SLARRB. */
/*                 However, we only allow the gaps to become greater since */
/*                 this is what should happen when we decrease WERR */
#line 598 "clarrv.f"
		    if (oldfst > 1) {
/* Computing MAX */
#line 599 "clarrv.f"
			d__1 = wgap[wbegin + oldfst - 2], d__2 = w[wbegin + 
				oldfst - 1] - werr[wbegin + oldfst - 1] - w[
				wbegin + oldfst - 2] - werr[wbegin + oldfst - 
				2];
#line 599 "clarrv.f"
			wgap[wbegin + oldfst - 2] = max(d__1,d__2);
#line 603 "clarrv.f"
		    }
#line 604 "clarrv.f"
		    if (wbegin + oldlst - 1 < wend) {
/* Computing MAX */
#line 605 "clarrv.f"
			d__1 = wgap[wbegin + oldlst - 1], d__2 = w[wbegin + 
				oldlst] - werr[wbegin + oldlst] - w[wbegin + 
				oldlst - 1] - werr[wbegin + oldlst - 1];
#line 605 "clarrv.f"
			wgap[wbegin + oldlst - 1] = max(d__1,d__2);
#line 609 "clarrv.f"
		    }
/*                 Each time the eigenvalues in WORK get refined, we store */
/*                 the newly found approximation with all shifts applied in W */
#line 612 "clarrv.f"
		    i__3 = oldlst;
#line 612 "clarrv.f"
		    for (j = oldfst; j <= i__3; ++j) {
#line 613 "clarrv.f"
			w[wbegin + j - 1] = work[wbegin + j - 1] + sigma;
#line 614 "clarrv.f"
/* L53: */
#line 614 "clarrv.f"
		    }
#line 615 "clarrv.f"
		}
/*              Process the current node. */
#line 618 "clarrv.f"
		newfst = oldfst;
#line 619 "clarrv.f"
		i__3 = oldlst;
#line 619 "clarrv.f"
		for (j = oldfst; j <= i__3; ++j) {
#line 620 "clarrv.f"
		    if (j == oldlst) {
/*                    we are at the right end of the cluster, this is also the */
/*                    boundary of the child cluster */
#line 623 "clarrv.f"
			newlst = j;
#line 624 "clarrv.f"
		    } else if (wgap[wbegin + j - 1] >= *minrgp * (d__1 = work[
			    wbegin + j - 1], abs(d__1))) {
/*                    the right relative gap is big enough, the child cluster */
/*                    (NEWFST,..,NEWLST) is well separated from the following */
#line 628 "clarrv.f"
			newlst = j;
#line 629 "clarrv.f"
		    } else {
/*                    inside a child cluster, the relative gap is not */
/*                    big enough. */
#line 632 "clarrv.f"
			goto L140;
#line 633 "clarrv.f"
		    }
/*                 Compute size of child cluster found */
#line 636 "clarrv.f"
		    newsiz = newlst - newfst + 1;
/*                 NEWFTT is the place in Z where the new RRR or the computed */
/*                 eigenvector is to be stored */
#line 640 "clarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Store representation at location of the leftmost evalue */
/*                    of the cluster */
#line 643 "clarrv.f"
			newftt = wbegin + newfst - 1;
#line 644 "clarrv.f"
		    } else {
#line 645 "clarrv.f"
			if (wbegin + newfst - 1 < *dol) {
/*                       Store representation at the left end of Z array */
#line 647 "clarrv.f"
			    newftt = *dol - 1;
#line 648 "clarrv.f"
			} else if (wbegin + newfst - 1 > *dou) {
/*                       Store representation at the right end of Z array */
#line 650 "clarrv.f"
			    newftt = *dou;
#line 651 "clarrv.f"
			} else {
#line 652 "clarrv.f"
			    newftt = wbegin + newfst - 1;
#line 653 "clarrv.f"
			}
#line 654 "clarrv.f"
		    }
#line 656 "clarrv.f"
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
#line 671 "clarrv.f"
			if (newfst == 1) {
/* Computing MAX */
#line 672 "clarrv.f"
			    d__1 = 0., d__2 = w[wbegin] - werr[wbegin] - *vl;
#line 672 "clarrv.f"
			    lgap = max(d__1,d__2);
#line 674 "clarrv.f"
			} else {
#line 675 "clarrv.f"
			    lgap = wgap[wbegin + newfst - 2];
#line 676 "clarrv.f"
			}
#line 677 "clarrv.f"
			rgap = wgap[wbegin + newlst - 1];

/*                    Compute left- and rightmost eigenvalue of child */
/*                    to high precision in order to shift as close */
/*                    as possible and obtain as large relative gaps */
/*                    as possible */

#line 684 "clarrv.f"
			for (k = 1; k <= 2; ++k) {
#line 685 "clarrv.f"
			    if (k == 1) {
#line 686 "clarrv.f"
				p = indexw[wbegin - 1 + newfst];
#line 687 "clarrv.f"
			    } else {
#line 688 "clarrv.f"
				p = indexw[wbegin - 1 + newlst];
#line 689 "clarrv.f"
			    }
#line 690 "clarrv.f"
			    offset = indexw[wbegin] - 1;
#line 691 "clarrv.f"
			    slarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &p, &p, &rqtol, &rqtol, &offset, &
				    work[wbegin], &wgap[wbegin], &werr[wbegin]
				    , &work[indwrk], &iwork[iindwk], pivmin, &
				    spdiam, &in, &iinfo);
#line 698 "clarrv.f"
/* L55: */
#line 698 "clarrv.f"
			}

#line 700 "clarrv.f"
			if (wbegin + newlst - 1 < *dol || wbegin + newfst - 1 
				> *dou) {
/*                       if the cluster contains no desired eigenvalues */
/*                       skip the computation of that branch of the rep. tree */

/*                       We could skip before the refinement of the extremal */
/*                       eigenvalues of the child, but then the representation */
/*                       tree could be different from the one when nothing is */
/*                       skipped. For this reason we skip at this place. */
#line 709 "clarrv.f"
			    idone = idone + newlst - newfst + 1;
#line 710 "clarrv.f"
			    goto L139;
#line 711 "clarrv.f"
			}

/*                    Compute RRR of child cluster. */
/*                    Note that the new RRR is stored in Z */

/*                    SLARRF needs LWORK = 2*N */
#line 717 "clarrv.f"
			slarrf_(&in, &d__[ibegin], &l[ibegin], &work[indld + 
				ibegin - 1], &newfst, &newlst, &work[wbegin], 
				&wgap[wbegin], &werr[wbegin], &spdiam, &lgap, 
				&rgap, pivmin, &tau, &work[indin1], &work[
				indin2], &work[indwrk], &iinfo);
/*                    In the complex case, SLARRF cannot write */
/*                    the new RRR directly into Z and needs an intermediate */
/*                    workspace */
#line 727 "clarrv.f"
			i__4 = in - 1;
#line 727 "clarrv.f"
			for (k = 1; k <= i__4; ++k) {
#line 728 "clarrv.f"
			    i__5 = ibegin + k - 1 + newftt * z_dim1;
#line 728 "clarrv.f"
			    i__6 = indin1 + k - 1;
#line 728 "clarrv.f"
			    z__1.r = work[i__6], z__1.i = 0.;
#line 728 "clarrv.f"
			    z__[i__5].r = z__1.r, z__[i__5].i = z__1.i;
#line 730 "clarrv.f"
			    i__5 = ibegin + k - 1 + (newftt + 1) * z_dim1;
#line 730 "clarrv.f"
			    i__6 = indin2 + k - 1;
#line 730 "clarrv.f"
			    z__1.r = work[i__6], z__1.i = 0.;
#line 730 "clarrv.f"
			    z__[i__5].r = z__1.r, z__[i__5].i = z__1.i;
#line 732 "clarrv.f"
/* L56: */
#line 732 "clarrv.f"
			}
#line 733 "clarrv.f"
			i__4 = iend + newftt * z_dim1;
#line 733 "clarrv.f"
			i__5 = indin1 + in - 1;
#line 733 "clarrv.f"
			z__1.r = work[i__5], z__1.i = 0.;
#line 733 "clarrv.f"
			z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 735 "clarrv.f"
			if (iinfo == 0) {
/*                       a new RRR for the cluster was found by SLARRF */
/*                       update shift and store it */
#line 738 "clarrv.f"
			    ssigma = sigma + tau;
#line 739 "clarrv.f"
			    i__4 = iend + (newftt + 1) * z_dim1;
#line 739 "clarrv.f"
			    z__1.r = ssigma, z__1.i = 0.;
#line 739 "clarrv.f"
			    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
/*                       WORK() are the midpoints and WERR() the semi-width */
/*                       Note that the entries in W are unchanged. */
#line 742 "clarrv.f"
			    i__4 = newlst;
#line 742 "clarrv.f"
			    for (k = newfst; k <= i__4; ++k) {
#line 743 "clarrv.f"
				fudge = eps * 3. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
#line 745 "clarrv.f"
				work[wbegin + k - 1] -= tau;
#line 747 "clarrv.f"
				fudge += eps * 4. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
/*                          Fudge errors */
#line 750 "clarrv.f"
				werr[wbegin + k - 1] += fudge;
/*                          Gaps are not fudged. Provided that WERR is small */
/*                          when eigenvalues are close, a zero gap indicates */
/*                          that a new representation is needed for resolving */
/*                          the cluster. A fudge could lead to a wrong decision */
/*                          of judging eigenvalues 'separated' which in */
/*                          reality are not. This could have a negative impact */
/*                          on the orthogonality of the computed eigenvectors. */
#line 759 "clarrv.f"
/* L116: */
#line 759 "clarrv.f"
			    }
#line 761 "clarrv.f"
			    ++nclus;
#line 762 "clarrv.f"
			    k = newcls + (nclus << 1);
#line 763 "clarrv.f"
			    iwork[k - 1] = newfst;
#line 764 "clarrv.f"
			    iwork[k] = newlst;
#line 765 "clarrv.f"
			} else {
#line 766 "clarrv.f"
			    *info = -2;
#line 767 "clarrv.f"
			    return 0;
#line 768 "clarrv.f"
			}
#line 769 "clarrv.f"
		    } else {

/*                    Compute eigenvector of singleton */

#line 773 "clarrv.f"
			iter = 0;

#line 775 "clarrv.f"
			tol = log((doublereal) in) * 4. * eps;

#line 777 "clarrv.f"
			k = newfst;
#line 778 "clarrv.f"
			windex = wbegin + k - 1;
/* Computing MAX */
#line 779 "clarrv.f"
			i__4 = windex - 1;
#line 779 "clarrv.f"
			windmn = max(i__4,1);
/* Computing MIN */
#line 780 "clarrv.f"
			i__4 = windex + 1;
#line 780 "clarrv.f"
			windpl = min(i__4,*m);
#line 781 "clarrv.f"
			lambda = work[windex];
#line 782 "clarrv.f"
			++done;
/*                    Check if eigenvector computation is to be skipped */
#line 784 "clarrv.f"
			if (windex < *dol || windex > *dou) {
#line 786 "clarrv.f"
			    eskip = TRUE_;
#line 787 "clarrv.f"
			    goto L125;
#line 788 "clarrv.f"
			} else {
#line 789 "clarrv.f"
			    eskip = FALSE_;
#line 790 "clarrv.f"
			}
#line 791 "clarrv.f"
			left = work[windex] - werr[windex];
#line 792 "clarrv.f"
			right = work[windex] + werr[windex];
#line 793 "clarrv.f"
			indeig = indexw[windex];
/*                    Note that since we compute the eigenpairs for a child, */
/*                    all eigenvalue approximations are w.r.t the same shift. */
/*                    In this case, the entries in WORK should be used for */
/*                    computing the gaps since they exhibit even very small */
/*                    differences in the eigenvalues, as opposed to the */
/*                    entries in W which might "look" the same. */
#line 801 "clarrv.f"
			if (k == 1) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VL, the formula */
/*                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA ) */
/*                       can lead to an overestimation of the left gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small left gap. */
/* Computing MAX */
#line 808 "clarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 808 "clarrv.f"
			    lgap = eps * max(d__1,d__2);
#line 809 "clarrv.f"
			} else {
#line 810 "clarrv.f"
			    lgap = wgap[windmn];
#line 811 "clarrv.f"
			}
#line 812 "clarrv.f"
			if (k == im) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VU, the formula */
/*                       can lead to an overestimation of the right gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small right gap. */
/* Computing MAX */
#line 818 "clarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 818 "clarrv.f"
			    rgap = eps * max(d__1,d__2);
#line 819 "clarrv.f"
			} else {
#line 820 "clarrv.f"
			    rgap = wgap[windex];
#line 821 "clarrv.f"
			}
#line 822 "clarrv.f"
			gap = min(lgap,rgap);
#line 823 "clarrv.f"
			if (k == 1 || k == im) {
/*                       The eigenvector support can become wrong */
/*                       because significant entries could be cut off due to a */
/*                       large GAPTOL parameter in LAR1V. Prevent this. */
#line 827 "clarrv.f"
			    gaptol = 0.;
#line 828 "clarrv.f"
			} else {
#line 829 "clarrv.f"
			    gaptol = gap * eps;
#line 830 "clarrv.f"
			}
#line 831 "clarrv.f"
			isupmn = in;
#line 832 "clarrv.f"
			isupmx = 1;
/*                    Update WGAP so that it holds the minimum gap */
/*                    to the left or the right. This is crucial in the */
/*                    case where bisection is used to ensure that the */
/*                    eigenvalue is refined up to the required precision. */
/*                    The correct value is restored afterwards. */
#line 838 "clarrv.f"
			savgap = wgap[windex];
#line 839 "clarrv.f"
			wgap[windex] = gap;
/*                    We want to use the Rayleigh Quotient Correction */
/*                    as often as possible since it converges quadratically */
/*                    when we are close enough to the desired eigenvalue. */
/*                    However, the Rayleigh Quotient can have the wrong sign */
/*                    and lead us away from the desired eigenvalue. In this */
/*                    case, the best we can do is to use bisection. */
#line 846 "clarrv.f"
			usedbs = FALSE_;
#line 847 "clarrv.f"
			usedrq = FALSE_;
/*                    Bisection is initially turned off unless it is forced */
#line 849 "clarrv.f"
			needbs = ! tryrqc;
#line 850 "clarrv.f"
L120:
/*                    Check if bisection should be used to refine eigenvalue */
#line 852 "clarrv.f"
			if (needbs) {
/*                       Take the bisection as new iterate */
#line 854 "clarrv.f"
			    usedbs = TRUE_;
#line 855 "clarrv.f"
			    itmp1 = iwork[iindr + windex];
#line 856 "clarrv.f"
			    offset = indexw[wbegin] - 1;
#line 857 "clarrv.f"
			    d__1 = eps * 2.;
#line 857 "clarrv.f"
			    slarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &indeig, &indeig, &c_b28, &d__1, &
				    offset, &work[wbegin], &wgap[wbegin], &
				    werr[wbegin], &work[indwrk], &iwork[
				    iindwk], pivmin, &spdiam, &itmp1, &iinfo);
#line 864 "clarrv.f"
			    if (iinfo != 0) {
#line 865 "clarrv.f"
				*info = -3;
#line 866 "clarrv.f"
				return 0;
#line 867 "clarrv.f"
			    }
#line 868 "clarrv.f"
			    lambda = work[windex];
/*                       Reset twist index from inaccurate LAMBDA to */
/*                       force computation of true MINGMA */
#line 871 "clarrv.f"
			    iwork[iindr + windex] = 0;
#line 872 "clarrv.f"
			}
/*                    Given LAMBDA, compute the eigenvector. */
#line 874 "clarrv.f"
			L__1 = ! usedbs;
#line 874 "clarrv.f"
			clar1v_(&in, &c__1, &in, &lambda, &d__[ibegin], &l[
				ibegin], &work[indld + ibegin - 1], &work[
				indlld + ibegin - 1], pivmin, &gaptol, &z__[
				ibegin + windex * z_dim1], &L__1, &negcnt, &
				ztz, &mingma, &iwork[iindr + windex], &isuppz[
				(windex << 1) - 1], &nrminv, &resid, &rqcorr, 
				&work[indwrk]);
#line 881 "clarrv.f"
			if (iter == 0) {
#line 882 "clarrv.f"
			    bstres = resid;
#line 883 "clarrv.f"
			    bstw = lambda;
#line 884 "clarrv.f"
			} else if (resid < bstres) {
#line 885 "clarrv.f"
			    bstres = resid;
#line 886 "clarrv.f"
			    bstw = lambda;
#line 887 "clarrv.f"
			}
/* Computing MIN */
#line 888 "clarrv.f"
			i__4 = isupmn, i__5 = isuppz[(windex << 1) - 1];
#line 888 "clarrv.f"
			isupmn = min(i__4,i__5);
/* Computing MAX */
#line 889 "clarrv.f"
			i__4 = isupmx, i__5 = isuppz[windex * 2];
#line 889 "clarrv.f"
			isupmx = max(i__4,i__5);
#line 890 "clarrv.f"
			++iter;
/*                    sin alpha <= |resid|/gap */
/*                    Note that both the residual and the gap are */
/*                    proportional to the matrix, so ||T|| doesn't play */
/*                    a role in the quotient */

/*                    Convergence test for Rayleigh-Quotient iteration */
/*                    (omitted when Bisection has been used) */

#line 901 "clarrv.f"
			if (resid > tol * gap && abs(rqcorr) > rqtol * abs(
				lambda) && ! usedbs) {
/*                       We need to check that the RQCORR update doesn't */
/*                       move the eigenvalue away from the desired one and */
/*                       towards a neighbor. -> protection with bisection */
#line 907 "clarrv.f"
			    if (indeig <= negcnt) {
/*                          The wanted eigenvalue lies to the left */
#line 909 "clarrv.f"
				sgndef = -1.;
#line 910 "clarrv.f"
			    } else {
/*                          The wanted eigenvalue lies to the right */
#line 912 "clarrv.f"
				sgndef = 1.;
#line 913 "clarrv.f"
			    }
/*                       We only use the RQCORR if it improves the */
/*                       the iterate reasonably. */
#line 916 "clarrv.f"
			    if (rqcorr * sgndef >= 0. && lambda + rqcorr <= 
				    right && lambda + rqcorr >= left) {
#line 920 "clarrv.f"
				usedrq = TRUE_;
/*                          Store new midpoint of bisection interval in WORK */
#line 922 "clarrv.f"
				if (sgndef == 1.) {
/*                             The current LAMBDA is on the left of the true */
/*                             eigenvalue */
#line 925 "clarrv.f"
				    left = lambda;
/*                             We prefer to assume that the error estimate */
/*                             is correct. We could make the interval not */
/*                             as a bracket but to be modified if the RQCORR */
/*                             chooses to. In this case, the RIGHT side should */
/*                             be modified as follows: */
/*                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR) */
#line 932 "clarrv.f"
				} else {
/*                             The current LAMBDA is on the right of the true */
/*                             eigenvalue */
#line 935 "clarrv.f"
				    right = lambda;
/*                             See comment about assuming the error estimate is */
/*                             correct above. */
/*                              LEFT = MIN(LEFT, LAMBDA + RQCORR) */
#line 939 "clarrv.f"
				}
#line 940 "clarrv.f"
				work[windex] = (right + left) * .5;
/*                          Take RQCORR since it has the correct sign and */
/*                          improves the iterate reasonably */
#line 944 "clarrv.f"
				lambda += rqcorr;
/*                          Update width of error interval */
#line 946 "clarrv.f"
				werr[windex] = (right - left) * .5;
#line 948 "clarrv.f"
			    } else {
#line 949 "clarrv.f"
				needbs = TRUE_;
#line 950 "clarrv.f"
			    }
#line 951 "clarrv.f"
			    if (right - left < rqtol * abs(lambda)) {
/*                             The eigenvalue is computed to bisection accuracy */
/*                             compute eigenvector and stop */
#line 954 "clarrv.f"
				usedbs = TRUE_;
#line 955 "clarrv.f"
				goto L120;
#line 956 "clarrv.f"
			    } else if (iter < 10) {
#line 957 "clarrv.f"
				goto L120;
#line 958 "clarrv.f"
			    } else if (iter == 10) {
#line 959 "clarrv.f"
				needbs = TRUE_;
#line 960 "clarrv.f"
				goto L120;
#line 961 "clarrv.f"
			    } else {
#line 962 "clarrv.f"
				*info = 5;
#line 963 "clarrv.f"
				return 0;
#line 964 "clarrv.f"
			    }
#line 965 "clarrv.f"
			} else {
#line 966 "clarrv.f"
			    stp2ii = FALSE_;
#line 967 "clarrv.f"
			    if (usedrq && usedbs && bstres <= resid) {
#line 969 "clarrv.f"
				lambda = bstw;
#line 970 "clarrv.f"
				stp2ii = TRUE_;
#line 971 "clarrv.f"
			    }
#line 972 "clarrv.f"
			    if (stp2ii) {
/*                          improve error angle by second step */
#line 974 "clarrv.f"
				L__1 = ! usedbs;
#line 974 "clarrv.f"
				clar1v_(&in, &c__1, &in, &lambda, &d__[ibegin]
					, &l[ibegin], &work[indld + ibegin - 
					1], &work[indlld + ibegin - 1], 
					pivmin, &gaptol, &z__[ibegin + windex 
					* z_dim1], &L__1, &negcnt, &ztz, &
					mingma, &iwork[iindr + windex], &
					isuppz[(windex << 1) - 1], &nrminv, &
					resid, &rqcorr, &work[indwrk]);
#line 983 "clarrv.f"
			    }
#line 984 "clarrv.f"
			    work[windex] = lambda;
#line 985 "clarrv.f"
			}

/*                    Compute FP-vector support w.r.t. whole matrix */

#line 989 "clarrv.f"
			isuppz[(windex << 1) - 1] += oldien;
#line 990 "clarrv.f"
			isuppz[windex * 2] += oldien;
#line 991 "clarrv.f"
			zfrom = isuppz[(windex << 1) - 1];
#line 992 "clarrv.f"
			zto = isuppz[windex * 2];
#line 993 "clarrv.f"
			isupmn += oldien;
#line 994 "clarrv.f"
			isupmx += oldien;
/*                    Ensure vector is ok if support in the RQI has changed */
#line 996 "clarrv.f"
			if (isupmn < zfrom) {
#line 997 "clarrv.f"
			    i__4 = zfrom - 1;
#line 997 "clarrv.f"
			    for (ii = isupmn; ii <= i__4; ++ii) {
#line 998 "clarrv.f"
				i__5 = ii + windex * z_dim1;
#line 998 "clarrv.f"
				z__[i__5].r = 0., z__[i__5].i = 0.;
#line 999 "clarrv.f"
/* L122: */
#line 999 "clarrv.f"
			    }
#line 1000 "clarrv.f"
			}
#line 1001 "clarrv.f"
			if (isupmx > zto) {
#line 1002 "clarrv.f"
			    i__4 = isupmx;
#line 1002 "clarrv.f"
			    for (ii = zto + 1; ii <= i__4; ++ii) {
#line 1003 "clarrv.f"
				i__5 = ii + windex * z_dim1;
#line 1003 "clarrv.f"
				z__[i__5].r = 0., z__[i__5].i = 0.;
#line 1004 "clarrv.f"
/* L123: */
#line 1004 "clarrv.f"
			    }
#line 1005 "clarrv.f"
			}
#line 1006 "clarrv.f"
			i__4 = zto - zfrom + 1;
#line 1006 "clarrv.f"
			csscal_(&i__4, &nrminv, &z__[zfrom + windex * z_dim1],
				 &c__1);
#line 1008 "clarrv.f"
L125:
/*                    Update W */
#line 1010 "clarrv.f"
			w[windex] = lambda + sigma;
/*                    Recompute the gaps on the left and right */
/*                    But only allow them to become larger and not */
/*                    smaller (which can only happen through "bad" */
/*                    cancellation and doesn't reflect the theory */
/*                    where the initial gaps are underestimated due */
/*                    to WERR being too crude.) */
#line 1017 "clarrv.f"
			if (! eskip) {
#line 1018 "clarrv.f"
			    if (k > 1) {
/* Computing MAX */
#line 1019 "clarrv.f"
				d__1 = wgap[windmn], d__2 = w[windex] - werr[
					windex] - w[windmn] - werr[windmn];
#line 1019 "clarrv.f"
				wgap[windmn] = max(d__1,d__2);
#line 1022 "clarrv.f"
			    }
#line 1023 "clarrv.f"
			    if (windex < wend) {
/* Computing MAX */
#line 1024 "clarrv.f"
				d__1 = savgap, d__2 = w[windpl] - werr[windpl]
					 - w[windex] - werr[windex];
#line 1024 "clarrv.f"
				wgap[windex] = max(d__1,d__2);
#line 1027 "clarrv.f"
			    }
#line 1028 "clarrv.f"
			}
#line 1029 "clarrv.f"
			++idone;
#line 1030 "clarrv.f"
		    }
/*                 here ends the code for the current child */

#line 1033 "clarrv.f"
L139:
/*                 Proceed to any remaining child nodes */
#line 1035 "clarrv.f"
		    newfst = j + 1;
#line 1036 "clarrv.f"
L140:
#line 1036 "clarrv.f"
		    ;
#line 1036 "clarrv.f"
		}
#line 1037 "clarrv.f"
/* L150: */
#line 1037 "clarrv.f"
	    }
#line 1038 "clarrv.f"
	    ++ndepth;
#line 1039 "clarrv.f"
	    goto L40;
#line 1040 "clarrv.f"
	}
#line 1041 "clarrv.f"
	ibegin = iend + 1;
#line 1042 "clarrv.f"
	wbegin = wend + 1;
#line 1043 "clarrv.f"
L170:
#line 1043 "clarrv.f"
	;
#line 1043 "clarrv.f"
    }

#line 1046 "clarrv.f"
    return 0;

/*     End of CLARRV */

} /* clarrv_ */

