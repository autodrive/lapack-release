#line 1 "dlarrv.f"
/* dlarrv.f -- translated by f2c (version 20100827).
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

#line 1 "dlarrv.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b DLARRV computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenv
alues of L D LT. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRV( N, VL, VU, D, L, PIVMIN, */
/*                          ISPLIT, M, DOL, DOU, MINRGP, */
/*                          RTOL1, RTOL2, W, WERR, WGAP, */
/*                          IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            DOL, DOU, INFO, LDZ, M, N */
/*       DOUBLE PRECISION   MINRGP, PIVMIN, RTOL1, RTOL2, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), INDEXW( * ), ISPLIT( * ), */
/*      $                   ISUPPZ( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), GERS( * ), L( * ), W( * ), WERR( * ), */
/*      $                   WGAP( * ), WORK( * ) */
/*       DOUBLE PRECISION  Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARRV computes the eigenvectors of the tridiagonal matrix */
/* > T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T. */
/* > The input eigenvalues should have been computed by DLARRE. */
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
/* >          VL is DOUBLE PRECISION */
/* >          Lower bound of the interval that contains the desired */
/* >          eigenvalues. VL < VU. Needed to compute gaps on the left or right */
/* >          end of the extremal eigenvalues in the desired RANGE. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          Upper bound of the interval that contains the desired */
/* >          eigenvalues. VL < VU. */
/* >          Note: VU is currently not used by this implementation of DLARRV, VU is */
/* >          passed to DLARRV because it could be used compute gaps on the right end */
/* >          of the extremal eigenvalues. However, with not much initial accuracy in */
/* >          LAMBDA and VU, the formula can lead to an overestimation of the right gap */
/* >          and thus to inadequately early RQI 'convergence'. This is currently */
/* >          prevented this by forcing a small right gap. And so it turns out that VU */
/* >          is currently not used by this implementation of DLARRV. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the N diagonal elements of the diagonal matrix D. */
/* >          On exit, D may be overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in,out] L */
/* > \verbatim */
/* >          L is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the (N-1) subdiagonal elements of the unit */
/* >          bidiagonal matrix L are in elements 1 to N-1 of L */
/* >          (if the matrix is not split.) At the end of each block */
/* >          is stored the corresponding shift as given by DLARRE. */
/* >          On exit, L is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
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
/* >          MINRGP is DOUBLE PRECISION */
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
/* > \param[in,out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          The first M elements of W contain the APPROXIMATE eigenvalues for */
/* >          which eigenvectors are to be computed.  The eigenvalues */
/* >          should be grouped by split-off block and ordered from */
/* >          smallest to largest within the block ( The output array */
/* >          W from DLARRE is expected here ). Furthermore, they are with */
/* >          respect to the shift of the corresponding root representation */
/* >          for their block. On exit, W holds the eigenvalues of the */
/* >          UNshifted matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WERR */
/* > \verbatim */
/* >          WERR is DOUBLE PRECISION array, dimension (N) */
/* >          The first M elements contain the semiwidth of the uncertainty */
/* >          interval of the corresponding eigenvalue in W */
/* > \endverbatim */
/* > */
/* > \param[in,out] WGAP */
/* > \verbatim */
/* >          WGAP is DOUBLE PRECISION array, dimension (N) */
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
/* >          GERS is DOUBLE PRECISION array, dimension (2*N) */
/* >          The N Gerschgorin intervals (the i-th Gerschgorin interval */
/* >          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should */
/* >          be computed from the original UNshifted matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) ) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (12*N) */
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
/* >          > 0:  A problem occurred in DLARRV. */
/* >          < 0:  One of the called subroutines signaled an internal problem. */
/* >                Needs inspection of the corresponding parameter IINFO */
/* >                for further information. */
/* > */
/* >          =-1:  Problem in DLARRB when refining a child's eigenvalues. */
/* >          =-2:  Problem in DLARRF when computing the RRR of a child. */
/* >                When a child is inside a tight cluster, it can be difficult */
/* >                to find an RRR. A partial remedy from the user's point of */
/* >                view is to make the parameter MINRGP smaller and recompile. */
/* >                However, as the orthogonality of the computed vectors is */
/* >                proportional to 1/MINRGP, the user should be aware that */
/* >                he might be trading in precision when he decreases MINRGP. */
/* >          =-3:  Problem in DLARRB when refining a single eigenvalue */
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

/* > \date June 2016 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dlarrv_(integer *n, doublereal *vl, doublereal *vu, 
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
    static integer itmp1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer indld;
    static doublereal fudge;
    static integer idone;
    static doublereal sigma;
    static integer iinfo, iindr;
    static doublereal resid;
    static logical eskip;
    static doublereal right;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nclus, zfrom;
    static doublereal rqtol;
    static integer iindc1, iindc2;
    extern /* Subroutine */ int dlar1v_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static logical stp2ii;
    static doublereal lambda;
    extern doublereal dlamch_(char *, ftnlen);
    static integer ibegin, indeig;
    static logical needbs;
    static integer indlld;
    static doublereal sgndef, mingma;
    extern /* Subroutine */ int dlarrb_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *);
    static integer oldien, oldncl, wbegin;
    static doublereal spdiam;
    static integer negcnt;
    extern /* Subroutine */ int dlarrf_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static integer oldcls;
    static doublereal savgap;
    static integer ndepth;
    static doublereal ssigma;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static logical usedbs;
    static integer iindwk, offset;
    static doublereal gaptol;
    static integer newcls, oldfst, indwrk, windex, oldlst;
    static logical usedrq;
    static integer newfst, newftt, parity, windmn, windpl, isupmn, newlst, 
	    zusedl;
    static doublereal bstres;
    static integer newsiz, zusedu, zusedw;
    static doublereal nrminv, rqcorr;
    static logical tryrqc;
    static integer isupmx;


/*  -- LAPACK auxiliary routine (version 3.8.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
#line 352 "dlarrv.f"
    /* Parameter adjustments */
#line 352 "dlarrv.f"
    --d__;
#line 352 "dlarrv.f"
    --l;
#line 352 "dlarrv.f"
    --isplit;
#line 352 "dlarrv.f"
    --w;
#line 352 "dlarrv.f"
    --werr;
#line 352 "dlarrv.f"
    --wgap;
#line 352 "dlarrv.f"
    --iblock;
#line 352 "dlarrv.f"
    --indexw;
#line 352 "dlarrv.f"
    --gers;
#line 352 "dlarrv.f"
    z_dim1 = *ldz;
#line 352 "dlarrv.f"
    z_offset = 1 + z_dim1;
#line 352 "dlarrv.f"
    z__ -= z_offset;
#line 352 "dlarrv.f"
    --isuppz;
#line 352 "dlarrv.f"
    --work;
#line 352 "dlarrv.f"
    --iwork;
#line 352 "dlarrv.f"

#line 352 "dlarrv.f"
    /* Function Body */
#line 352 "dlarrv.f"
    *info = 0;

/*     Quick return if possible */

#line 356 "dlarrv.f"
    if (*n <= 0) {
#line 357 "dlarrv.f"
	return 0;
#line 358 "dlarrv.f"
    }

/*     The first N entries of WORK are reserved for the eigenvalues */
#line 361 "dlarrv.f"
    indld = *n + 1;
#line 362 "dlarrv.f"
    indlld = (*n << 1) + 1;
#line 363 "dlarrv.f"
    indwrk = *n * 3 + 1;
#line 364 "dlarrv.f"
    minwsize = *n * 12;
#line 366 "dlarrv.f"
    i__1 = minwsize;
#line 366 "dlarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 367 "dlarrv.f"
	work[i__] = 0.;
#line 368 "dlarrv.f"
/* L5: */
#line 368 "dlarrv.f"
    }
/*     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the */
/*     factorization used to compute the FP vector */
#line 372 "dlarrv.f"
    iindr = 0;
/*     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current */
/*     layer and the one above. */
#line 375 "dlarrv.f"
    iindc1 = *n;
#line 376 "dlarrv.f"
    iindc2 = *n << 1;
#line 377 "dlarrv.f"
    iindwk = *n * 3 + 1;
#line 379 "dlarrv.f"
    miniwsize = *n * 7;
#line 380 "dlarrv.f"
    i__1 = miniwsize;
#line 380 "dlarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 381 "dlarrv.f"
	iwork[i__] = 0;
#line 382 "dlarrv.f"
/* L10: */
#line 382 "dlarrv.f"
    }
#line 384 "dlarrv.f"
    zusedl = 1;
#line 385 "dlarrv.f"
    if (*dol > 1) {
/*        Set lower bound for use of Z */
#line 387 "dlarrv.f"
	zusedl = *dol - 1;
#line 388 "dlarrv.f"
    }
#line 389 "dlarrv.f"
    zusedu = *m;
#line 390 "dlarrv.f"
    if (*dou < *m) {
/*        Set lower bound for use of Z */
#line 392 "dlarrv.f"
	zusedu = *dou + 1;
#line 393 "dlarrv.f"
    }
/*     The width of the part of Z that is used */
#line 395 "dlarrv.f"
    zusedw = zusedu - zusedl + 1;
#line 398 "dlarrv.f"
    dlaset_("Full", n, &zusedw, &c_b5, &c_b5, &z__[zusedl * z_dim1 + 1], ldz, 
	    (ftnlen)4);
#line 401 "dlarrv.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 402 "dlarrv.f"
    rqtol = eps * 2.;

/*     Set expert flags for standard code. */
#line 405 "dlarrv.f"
    tryrqc = TRUE_;
#line 407 "dlarrv.f"
    if (*dol == 1 && *dou == *m) {
#line 408 "dlarrv.f"
    } else {
/*        Only selected eigenpairs are computed. Since the other evalues */
/*        are not refined by RQ iteration, bisection has to compute to full */
/*        accuracy. */
#line 412 "dlarrv.f"
	*rtol1 = eps * 4.;
#line 413 "dlarrv.f"
	*rtol2 = eps * 4.;
#line 414 "dlarrv.f"
    }
/*     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the */
/*     desired eigenvalues. The support of the nonzero eigenvector */
/*     entries is contained in the interval IBEGIN:IEND. */
/*     Remark that if k eigenpairs are desired, then the eigenvectors */
/*     are stored in k contiguous columns of Z. */
/*     DONE is the number of eigenvectors already computed */
#line 423 "dlarrv.f"
    done = 0;
#line 424 "dlarrv.f"
    ibegin = 1;
#line 425 "dlarrv.f"
    wbegin = 1;
#line 426 "dlarrv.f"
    i__1 = iblock[*m];
#line 426 "dlarrv.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 427 "dlarrv.f"
	iend = isplit[jblk];
#line 428 "dlarrv.f"
	sigma = l[iend];
/*        Find the eigenvectors of the submatrix indexed IBEGIN */
/*        through IEND. */
#line 431 "dlarrv.f"
	wend = wbegin - 1;
#line 432 "dlarrv.f"
L15:
#line 433 "dlarrv.f"
	if (wend < *m) {
#line 434 "dlarrv.f"
	    if (iblock[wend + 1] == jblk) {
#line 435 "dlarrv.f"
		++wend;
#line 436 "dlarrv.f"
		goto L15;
#line 437 "dlarrv.f"
	    }
#line 438 "dlarrv.f"
	}
#line 439 "dlarrv.f"
	if (wend < wbegin) {
#line 440 "dlarrv.f"
	    ibegin = iend + 1;
#line 441 "dlarrv.f"
	    goto L170;
#line 442 "dlarrv.f"
	} else if (wend < *dol || wbegin > *dou) {
#line 443 "dlarrv.f"
	    ibegin = iend + 1;
#line 444 "dlarrv.f"
	    wbegin = wend + 1;
#line 445 "dlarrv.f"
	    goto L170;
#line 446 "dlarrv.f"
	}
/*        Find local spectral diameter of the block */
#line 449 "dlarrv.f"
	gl = gers[(ibegin << 1) - 1];
#line 450 "dlarrv.f"
	gu = gers[ibegin * 2];
#line 451 "dlarrv.f"
	i__2 = iend;
#line 451 "dlarrv.f"
	for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 452 "dlarrv.f"
	    d__1 = gers[(i__ << 1) - 1];
#line 452 "dlarrv.f"
	    gl = min(d__1,gl);
/* Computing MAX */
#line 453 "dlarrv.f"
	    d__1 = gers[i__ * 2];
#line 453 "dlarrv.f"
	    gu = max(d__1,gu);
#line 454 "dlarrv.f"
/* L20: */
#line 454 "dlarrv.f"
	}
#line 455 "dlarrv.f"
	spdiam = gu - gl;
/*        OLDIEN is the last index of the previous block */
#line 458 "dlarrv.f"
	oldien = ibegin - 1;
/*        Calculate the size of the current block */
#line 460 "dlarrv.f"
	in = iend - ibegin + 1;
/*        The number of eigenvalues in the current block */
#line 462 "dlarrv.f"
	im = wend - wbegin + 1;
/*        This is for a 1x1 block */
#line 465 "dlarrv.f"
	if (ibegin == iend) {
#line 466 "dlarrv.f"
	    ++done;
#line 467 "dlarrv.f"
	    z__[ibegin + wbegin * z_dim1] = 1.;
#line 468 "dlarrv.f"
	    isuppz[(wbegin << 1) - 1] = ibegin;
#line 469 "dlarrv.f"
	    isuppz[wbegin * 2] = ibegin;
#line 470 "dlarrv.f"
	    w[wbegin] += sigma;
#line 471 "dlarrv.f"
	    work[wbegin] = w[wbegin];
#line 472 "dlarrv.f"
	    ibegin = iend + 1;
#line 473 "dlarrv.f"
	    ++wbegin;
#line 474 "dlarrv.f"
	    goto L170;
#line 475 "dlarrv.f"
	}
/*        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND) */
/*        Note that these can be approximations, in this case, the corresp. */
/*        entries of WERR give the size of the uncertainty interval. */
/*        The eigenvalue approximations will be refined when necessary as */
/*        high relative accuracy is required for the computation of the */
/*        corresponding eigenvectors. */
#line 483 "dlarrv.f"
	dcopy_(&im, &w[wbegin], &c__1, &work[wbegin], &c__1);
/*        We store in W the eigenvalue approximations w.r.t. the original */
/*        matrix T. */
#line 488 "dlarrv.f"
	i__2 = im;
#line 488 "dlarrv.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 489 "dlarrv.f"
	    w[wbegin + i__ - 1] += sigma;
#line 490 "dlarrv.f"
/* L30: */
#line 490 "dlarrv.f"
	}
/*        NDEPTH is the current depth of the representation tree */
#line 494 "dlarrv.f"
	ndepth = 0;
/*        PARITY is either 1 or 0 */
#line 496 "dlarrv.f"
	parity = 1;
/*        NCLUS is the number of clusters for the next level of the */
/*        representation tree, we start with NCLUS = 1 for the root */
#line 499 "dlarrv.f"
	nclus = 1;
#line 500 "dlarrv.f"
	iwork[iindc1 + 1] = 1;
#line 501 "dlarrv.f"
	iwork[iindc1 + 2] = im;
/*        IDONE is the number of eigenvectors already computed in the current */
/*        block */
#line 505 "dlarrv.f"
	idone = 0;
/*        loop while( IDONE.LT.IM ) */
/*        generate the representation tree for the current block and */
/*        compute the eigenvectors */
#line 509 "dlarrv.f"
L40:
#line 510 "dlarrv.f"
	if (idone < im) {
/*           This is a crude protection against infinitely deep trees */
#line 512 "dlarrv.f"
	    if (ndepth > *m) {
#line 513 "dlarrv.f"
		*info = -2;
#line 514 "dlarrv.f"
		return 0;
#line 515 "dlarrv.f"
	    }
/*           breadth first processing of the current level of the representation */
/*           tree: OLDNCL = number of clusters on current level */
#line 518 "dlarrv.f"
	    oldncl = nclus;
/*           reset NCLUS to count the number of child clusters */
#line 520 "dlarrv.f"
	    nclus = 0;

#line 522 "dlarrv.f"
	    parity = 1 - parity;
#line 523 "dlarrv.f"
	    if (parity == 0) {
#line 524 "dlarrv.f"
		oldcls = iindc1;
#line 525 "dlarrv.f"
		newcls = iindc2;
#line 526 "dlarrv.f"
	    } else {
#line 527 "dlarrv.f"
		oldcls = iindc2;
#line 528 "dlarrv.f"
		newcls = iindc1;
#line 529 "dlarrv.f"
	    }
/*           Process the clusters on the current level */
#line 531 "dlarrv.f"
	    i__2 = oldncl;
#line 531 "dlarrv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 532 "dlarrv.f"
		j = oldcls + (i__ << 1);
/*              OLDFST, OLDLST = first, last index of current cluster. */
/*                               cluster indices start with 1 and are relative */
/*                               to WBEGIN when accessing W, WGAP, WERR, Z */
#line 536 "dlarrv.f"
		oldfst = iwork[j - 1];
#line 537 "dlarrv.f"
		oldlst = iwork[j];
#line 538 "dlarrv.f"
		if (ndepth > 0) {
/*                 Retrieve relatively robust representation (RRR) of cluster */
/*                 that has been computed at the previous level */
/*                 The RRR is stored in Z and overwritten once the eigenvectors */
/*                 have been computed or when the cluster is refined */
#line 544 "dlarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Get representation from location of the leftmost evalue */
/*                    of the cluster */
#line 547 "dlarrv.f"
			j = wbegin + oldfst - 1;
#line 548 "dlarrv.f"
		    } else {
#line 549 "dlarrv.f"
			if (wbegin + oldfst - 1 < *dol) {
/*                       Get representation from the left end of Z array */
#line 551 "dlarrv.f"
			    j = *dol - 1;
#line 552 "dlarrv.f"
			} else if (wbegin + oldfst - 1 > *dou) {
/*                       Get representation from the right end of Z array */
#line 554 "dlarrv.f"
			    j = *dou;
#line 555 "dlarrv.f"
			} else {
#line 556 "dlarrv.f"
			    j = wbegin + oldfst - 1;
#line 557 "dlarrv.f"
			}
#line 558 "dlarrv.f"
		    }
#line 559 "dlarrv.f"
		    dcopy_(&in, &z__[ibegin + j * z_dim1], &c__1, &d__[ibegin]
			    , &c__1);
#line 560 "dlarrv.f"
		    i__3 = in - 1;
#line 560 "dlarrv.f"
		    dcopy_(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1, &l[
			    ibegin], &c__1);
#line 562 "dlarrv.f"
		    sigma = z__[iend + (j + 1) * z_dim1];
/*                 Set the corresponding entries in Z to zero */
#line 565 "dlarrv.f"
		    dlaset_("Full", &in, &c__2, &c_b5, &c_b5, &z__[ibegin + j 
			    * z_dim1], ldz, (ftnlen)4);
#line 567 "dlarrv.f"
		}
/*              Compute DL and DLL of current RRR */
#line 570 "dlarrv.f"
		i__3 = iend - 1;
#line 570 "dlarrv.f"
		for (j = ibegin; j <= i__3; ++j) {
#line 571 "dlarrv.f"
		    tmp = d__[j] * l[j];
#line 572 "dlarrv.f"
		    work[indld - 1 + j] = tmp;
#line 573 "dlarrv.f"
		    work[indlld - 1 + j] = tmp * l[j];
#line 574 "dlarrv.f"
/* L50: */
#line 574 "dlarrv.f"
		}
#line 576 "dlarrv.f"
		if (ndepth > 0) {
/*                 P and Q are index of the first and last eigenvalue to compute */
/*                 within the current block */
#line 579 "dlarrv.f"
		    p = indexw[wbegin - 1 + oldfst];
#line 580 "dlarrv.f"
		    q = indexw[wbegin - 1 + oldlst];
/*                 Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET */
/*                 through the Q-OFFSET elements of these arrays are to be used. */
/*                  OFFSET = P-OLDFST */
#line 584 "dlarrv.f"
		    offset = indexw[wbegin] - 1;
/*                 perform limited bisection (if necessary) to get approximate */
/*                 eigenvalues to the precision needed. */
#line 587 "dlarrv.f"
		    dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin - 1], &p,
			     &q, rtol1, rtol2, &offset, &work[wbegin], &wgap[
			    wbegin], &werr[wbegin], &work[indwrk], &iwork[
			    iindwk], pivmin, &spdiam, &in, &iinfo);
#line 593 "dlarrv.f"
		    if (iinfo != 0) {
#line 594 "dlarrv.f"
			*info = -1;
#line 595 "dlarrv.f"
			return 0;
#line 596 "dlarrv.f"
		    }
/*                 We also recompute the extremal gaps. W holds all eigenvalues */
/*                 of the unshifted matrix and must be used for computation */
/*                 of WGAP, the entries of WORK might stem from RRRs with */
/*                 different shifts. The gaps from WBEGIN-1+OLDFST to */
/*                 WBEGIN-1+OLDLST are correctly computed in DLARRB. */
/*                 However, we only allow the gaps to become greater since */
/*                 this is what should happen when we decrease WERR */
#line 604 "dlarrv.f"
		    if (oldfst > 1) {
/* Computing MAX */
#line 605 "dlarrv.f"
			d__1 = wgap[wbegin + oldfst - 2], d__2 = w[wbegin + 
				oldfst - 1] - werr[wbegin + oldfst - 1] - w[
				wbegin + oldfst - 2] - werr[wbegin + oldfst - 
				2];
#line 605 "dlarrv.f"
			wgap[wbegin + oldfst - 2] = max(d__1,d__2);
#line 609 "dlarrv.f"
		    }
#line 610 "dlarrv.f"
		    if (wbegin + oldlst - 1 < wend) {
/* Computing MAX */
#line 611 "dlarrv.f"
			d__1 = wgap[wbegin + oldlst - 1], d__2 = w[wbegin + 
				oldlst] - werr[wbegin + oldlst] - w[wbegin + 
				oldlst - 1] - werr[wbegin + oldlst - 1];
#line 611 "dlarrv.f"
			wgap[wbegin + oldlst - 1] = max(d__1,d__2);
#line 615 "dlarrv.f"
		    }
/*                 Each time the eigenvalues in WORK get refined, we store */
/*                 the newly found approximation with all shifts applied in W */
#line 618 "dlarrv.f"
		    i__3 = oldlst;
#line 618 "dlarrv.f"
		    for (j = oldfst; j <= i__3; ++j) {
#line 619 "dlarrv.f"
			w[wbegin + j - 1] = work[wbegin + j - 1] + sigma;
#line 620 "dlarrv.f"
/* L53: */
#line 620 "dlarrv.f"
		    }
#line 621 "dlarrv.f"
		}
/*              Process the current node. */
#line 624 "dlarrv.f"
		newfst = oldfst;
#line 625 "dlarrv.f"
		i__3 = oldlst;
#line 625 "dlarrv.f"
		for (j = oldfst; j <= i__3; ++j) {
#line 626 "dlarrv.f"
		    if (j == oldlst) {
/*                    we are at the right end of the cluster, this is also the */
/*                    boundary of the child cluster */
#line 629 "dlarrv.f"
			newlst = j;
#line 630 "dlarrv.f"
		    } else if (wgap[wbegin + j - 1] >= *minrgp * (d__1 = work[
			    wbegin + j - 1], abs(d__1))) {
/*                    the right relative gap is big enough, the child cluster */
/*                    (NEWFST,..,NEWLST) is well separated from the following */
#line 634 "dlarrv.f"
			newlst = j;
#line 635 "dlarrv.f"
		    } else {
/*                    inside a child cluster, the relative gap is not */
/*                    big enough. */
#line 638 "dlarrv.f"
			goto L140;
#line 639 "dlarrv.f"
		    }
/*                 Compute size of child cluster found */
#line 642 "dlarrv.f"
		    newsiz = newlst - newfst + 1;
/*                 NEWFTT is the place in Z where the new RRR or the computed */
/*                 eigenvector is to be stored */
#line 646 "dlarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Store representation at location of the leftmost evalue */
/*                    of the cluster */
#line 649 "dlarrv.f"
			newftt = wbegin + newfst - 1;
#line 650 "dlarrv.f"
		    } else {
#line 651 "dlarrv.f"
			if (wbegin + newfst - 1 < *dol) {
/*                       Store representation at the left end of Z array */
#line 653 "dlarrv.f"
			    newftt = *dol - 1;
#line 654 "dlarrv.f"
			} else if (wbegin + newfst - 1 > *dou) {
/*                       Store representation at the right end of Z array */
#line 656 "dlarrv.f"
			    newftt = *dou;
#line 657 "dlarrv.f"
			} else {
#line 658 "dlarrv.f"
			    newftt = wbegin + newfst - 1;
#line 659 "dlarrv.f"
			}
#line 660 "dlarrv.f"
		    }
#line 662 "dlarrv.f"
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
#line 677 "dlarrv.f"
			if (newfst == 1) {
/* Computing MAX */
#line 678 "dlarrv.f"
			    d__1 = 0., d__2 = w[wbegin] - werr[wbegin] - *vl;
#line 678 "dlarrv.f"
			    lgap = max(d__1,d__2);
#line 680 "dlarrv.f"
			} else {
#line 681 "dlarrv.f"
			    lgap = wgap[wbegin + newfst - 2];
#line 682 "dlarrv.f"
			}
#line 683 "dlarrv.f"
			rgap = wgap[wbegin + newlst - 1];

/*                    Compute left- and rightmost eigenvalue of child */
/*                    to high precision in order to shift as close */
/*                    as possible and obtain as large relative gaps */
/*                    as possible */

#line 690 "dlarrv.f"
			for (k = 1; k <= 2; ++k) {
#line 691 "dlarrv.f"
			    if (k == 1) {
#line 692 "dlarrv.f"
				p = indexw[wbegin - 1 + newfst];
#line 693 "dlarrv.f"
			    } else {
#line 694 "dlarrv.f"
				p = indexw[wbegin - 1 + newlst];
#line 695 "dlarrv.f"
			    }
#line 696 "dlarrv.f"
			    offset = indexw[wbegin] - 1;
#line 697 "dlarrv.f"
			    dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &p, &p, &rqtol, &rqtol, &offset, &
				    work[wbegin], &wgap[wbegin], &werr[wbegin]
				    , &work[indwrk], &iwork[iindwk], pivmin, &
				    spdiam, &in, &iinfo);
#line 704 "dlarrv.f"
/* L55: */
#line 704 "dlarrv.f"
			}

#line 706 "dlarrv.f"
			if (wbegin + newlst - 1 < *dol || wbegin + newfst - 1 
				> *dou) {
/*                       if the cluster contains no desired eigenvalues */
/*                       skip the computation of that branch of the rep. tree */

/*                       We could skip before the refinement of the extremal */
/*                       eigenvalues of the child, but then the representation */
/*                       tree could be different from the one when nothing is */
/*                       skipped. For this reason we skip at this place. */
#line 715 "dlarrv.f"
			    idone = idone + newlst - newfst + 1;
#line 716 "dlarrv.f"
			    goto L139;
#line 717 "dlarrv.f"
			}

/*                    Compute RRR of child cluster. */
/*                    Note that the new RRR is stored in Z */

/*                    DLARRF needs LWORK = 2*N */
#line 723 "dlarrv.f"
			dlarrf_(&in, &d__[ibegin], &l[ibegin], &work[indld + 
				ibegin - 1], &newfst, &newlst, &work[wbegin], 
				&wgap[wbegin], &werr[wbegin], &spdiam, &lgap, 
				&rgap, pivmin, &tau, &z__[ibegin + newftt * 
				z_dim1], &z__[ibegin + (newftt + 1) * z_dim1],
				 &work[indwrk], &iinfo);
#line 730 "dlarrv.f"
			if (iinfo == 0) {
/*                       a new RRR for the cluster was found by DLARRF */
/*                       update shift and store it */
#line 733 "dlarrv.f"
			    ssigma = sigma + tau;
#line 734 "dlarrv.f"
			    z__[iend + (newftt + 1) * z_dim1] = ssigma;
/*                       WORK() are the midpoints and WERR() the semi-width */
/*                       Note that the entries in W are unchanged. */
#line 737 "dlarrv.f"
			    i__4 = newlst;
#line 737 "dlarrv.f"
			    for (k = newfst; k <= i__4; ++k) {
#line 738 "dlarrv.f"
				fudge = eps * 3. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
#line 740 "dlarrv.f"
				work[wbegin + k - 1] -= tau;
#line 742 "dlarrv.f"
				fudge += eps * 4. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
/*                          Fudge errors */
#line 745 "dlarrv.f"
				werr[wbegin + k - 1] += fudge;
/*                          Gaps are not fudged. Provided that WERR is small */
/*                          when eigenvalues are close, a zero gap indicates */
/*                          that a new representation is needed for resolving */
/*                          the cluster. A fudge could lead to a wrong decision */
/*                          of judging eigenvalues 'separated' which in */
/*                          reality are not. This could have a negative impact */
/*                          on the orthogonality of the computed eigenvectors. */
#line 754 "dlarrv.f"
/* L116: */
#line 754 "dlarrv.f"
			    }
#line 756 "dlarrv.f"
			    ++nclus;
#line 757 "dlarrv.f"
			    k = newcls + (nclus << 1);
#line 758 "dlarrv.f"
			    iwork[k - 1] = newfst;
#line 759 "dlarrv.f"
			    iwork[k] = newlst;
#line 760 "dlarrv.f"
			} else {
#line 761 "dlarrv.f"
			    *info = -2;
#line 762 "dlarrv.f"
			    return 0;
#line 763 "dlarrv.f"
			}
#line 764 "dlarrv.f"
		    } else {

/*                    Compute eigenvector of singleton */

#line 768 "dlarrv.f"
			iter = 0;

#line 770 "dlarrv.f"
			tol = log((doublereal) in) * 4. * eps;

#line 772 "dlarrv.f"
			k = newfst;
#line 773 "dlarrv.f"
			windex = wbegin + k - 1;
/* Computing MAX */
#line 774 "dlarrv.f"
			i__4 = windex - 1;
#line 774 "dlarrv.f"
			windmn = max(i__4,1);
/* Computing MIN */
#line 775 "dlarrv.f"
			i__4 = windex + 1;
#line 775 "dlarrv.f"
			windpl = min(i__4,*m);
#line 776 "dlarrv.f"
			lambda = work[windex];
#line 777 "dlarrv.f"
			++done;
/*                    Check if eigenvector computation is to be skipped */
#line 779 "dlarrv.f"
			if (windex < *dol || windex > *dou) {
#line 781 "dlarrv.f"
			    eskip = TRUE_;
#line 782 "dlarrv.f"
			    goto L125;
#line 783 "dlarrv.f"
			} else {
#line 784 "dlarrv.f"
			    eskip = FALSE_;
#line 785 "dlarrv.f"
			}
#line 786 "dlarrv.f"
			left = work[windex] - werr[windex];
#line 787 "dlarrv.f"
			right = work[windex] + werr[windex];
#line 788 "dlarrv.f"
			indeig = indexw[windex];
/*                    Note that since we compute the eigenpairs for a child, */
/*                    all eigenvalue approximations are w.r.t the same shift. */
/*                    In this case, the entries in WORK should be used for */
/*                    computing the gaps since they exhibit even very small */
/*                    differences in the eigenvalues, as opposed to the */
/*                    entries in W which might "look" the same. */
#line 796 "dlarrv.f"
			if (k == 1) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VL, the formula */
/*                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA ) */
/*                       can lead to an overestimation of the left gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small left gap. */
/* Computing MAX */
#line 803 "dlarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 803 "dlarrv.f"
			    lgap = eps * max(d__1,d__2);
#line 804 "dlarrv.f"
			} else {
#line 805 "dlarrv.f"
			    lgap = wgap[windmn];
#line 806 "dlarrv.f"
			}
#line 807 "dlarrv.f"
			if (k == im) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VU, the formula */
/*                       can lead to an overestimation of the right gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small right gap. */
/* Computing MAX */
#line 813 "dlarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 813 "dlarrv.f"
			    rgap = eps * max(d__1,d__2);
#line 814 "dlarrv.f"
			} else {
#line 815 "dlarrv.f"
			    rgap = wgap[windex];
#line 816 "dlarrv.f"
			}
#line 817 "dlarrv.f"
			gap = min(lgap,rgap);
#line 818 "dlarrv.f"
			if (k == 1 || k == im) {
/*                       The eigenvector support can become wrong */
/*                       because significant entries could be cut off due to a */
/*                       large GAPTOL parameter in LAR1V. Prevent this. */
#line 822 "dlarrv.f"
			    gaptol = 0.;
#line 823 "dlarrv.f"
			} else {
#line 824 "dlarrv.f"
			    gaptol = gap * eps;
#line 825 "dlarrv.f"
			}
#line 826 "dlarrv.f"
			isupmn = in;
#line 827 "dlarrv.f"
			isupmx = 1;
/*                    Update WGAP so that it holds the minimum gap */
/*                    to the left or the right. This is crucial in the */
/*                    case where bisection is used to ensure that the */
/*                    eigenvalue is refined up to the required precision. */
/*                    The correct value is restored afterwards. */
#line 833 "dlarrv.f"
			savgap = wgap[windex];
#line 834 "dlarrv.f"
			wgap[windex] = gap;
/*                    We want to use the Rayleigh Quotient Correction */
/*                    as often as possible since it converges quadratically */
/*                    when we are close enough to the desired eigenvalue. */
/*                    However, the Rayleigh Quotient can have the wrong sign */
/*                    and lead us away from the desired eigenvalue. In this */
/*                    case, the best we can do is to use bisection. */
#line 841 "dlarrv.f"
			usedbs = FALSE_;
#line 842 "dlarrv.f"
			usedrq = FALSE_;
/*                    Bisection is initially turned off unless it is forced */
#line 844 "dlarrv.f"
			needbs = ! tryrqc;
#line 845 "dlarrv.f"
L120:
/*                    Check if bisection should be used to refine eigenvalue */
#line 847 "dlarrv.f"
			if (needbs) {
/*                       Take the bisection as new iterate */
#line 849 "dlarrv.f"
			    usedbs = TRUE_;
#line 850 "dlarrv.f"
			    itmp1 = iwork[iindr + windex];
#line 851 "dlarrv.f"
			    offset = indexw[wbegin] - 1;
#line 852 "dlarrv.f"
			    d__1 = eps * 2.;
#line 852 "dlarrv.f"
			    dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &indeig, &indeig, &c_b5, &d__1, &
				    offset, &work[wbegin], &wgap[wbegin], &
				    werr[wbegin], &work[indwrk], &iwork[
				    iindwk], pivmin, &spdiam, &itmp1, &iinfo);
#line 859 "dlarrv.f"
			    if (iinfo != 0) {
#line 860 "dlarrv.f"
				*info = -3;
#line 861 "dlarrv.f"
				return 0;
#line 862 "dlarrv.f"
			    }
#line 863 "dlarrv.f"
			    lambda = work[windex];
/*                       Reset twist index from inaccurate LAMBDA to */
/*                       force computation of true MINGMA */
#line 866 "dlarrv.f"
			    iwork[iindr + windex] = 0;
#line 867 "dlarrv.f"
			}
/*                    Given LAMBDA, compute the eigenvector. */
#line 869 "dlarrv.f"
			L__1 = ! usedbs;
#line 869 "dlarrv.f"
			dlar1v_(&in, &c__1, &in, &lambda, &d__[ibegin], &l[
				ibegin], &work[indld + ibegin - 1], &work[
				indlld + ibegin - 1], pivmin, &gaptol, &z__[
				ibegin + windex * z_dim1], &L__1, &negcnt, &
				ztz, &mingma, &iwork[iindr + windex], &isuppz[
				(windex << 1) - 1], &nrminv, &resid, &rqcorr, 
				&work[indwrk]);
#line 876 "dlarrv.f"
			if (iter == 0) {
#line 877 "dlarrv.f"
			    bstres = resid;
#line 878 "dlarrv.f"
			    bstw = lambda;
#line 879 "dlarrv.f"
			} else if (resid < bstres) {
#line 880 "dlarrv.f"
			    bstres = resid;
#line 881 "dlarrv.f"
			    bstw = lambda;
#line 882 "dlarrv.f"
			}
/* Computing MIN */
#line 883 "dlarrv.f"
			i__4 = isupmn, i__5 = isuppz[(windex << 1) - 1];
#line 883 "dlarrv.f"
			isupmn = min(i__4,i__5);
/* Computing MAX */
#line 884 "dlarrv.f"
			i__4 = isupmx, i__5 = isuppz[windex * 2];
#line 884 "dlarrv.f"
			isupmx = max(i__4,i__5);
#line 885 "dlarrv.f"
			++iter;
/*                    sin alpha <= |resid|/gap */
/*                    Note that both the residual and the gap are */
/*                    proportional to the matrix, so ||T|| doesn't play */
/*                    a role in the quotient */

/*                    Convergence test for Rayleigh-Quotient iteration */
/*                    (omitted when Bisection has been used) */

#line 896 "dlarrv.f"
			if (resid > tol * gap && abs(rqcorr) > rqtol * abs(
				lambda) && ! usedbs) {
/*                       We need to check that the RQCORR update doesn't */
/*                       move the eigenvalue away from the desired one and */
/*                       towards a neighbor. -> protection with bisection */
#line 902 "dlarrv.f"
			    if (indeig <= negcnt) {
/*                          The wanted eigenvalue lies to the left */
#line 904 "dlarrv.f"
				sgndef = -1.;
#line 905 "dlarrv.f"
			    } else {
/*                          The wanted eigenvalue lies to the right */
#line 907 "dlarrv.f"
				sgndef = 1.;
#line 908 "dlarrv.f"
			    }
/*                       We only use the RQCORR if it improves the */
/*                       the iterate reasonably. */
#line 911 "dlarrv.f"
			    if (rqcorr * sgndef >= 0. && lambda + rqcorr <= 
				    right && lambda + rqcorr >= left) {
#line 915 "dlarrv.f"
				usedrq = TRUE_;
/*                          Store new midpoint of bisection interval in WORK */
#line 917 "dlarrv.f"
				if (sgndef == 1.) {
/*                             The current LAMBDA is on the left of the true */
/*                             eigenvalue */
#line 920 "dlarrv.f"
				    left = lambda;
/*                             We prefer to assume that the error estimate */
/*                             is correct. We could make the interval not */
/*                             as a bracket but to be modified if the RQCORR */
/*                             chooses to. In this case, the RIGHT side should */
/*                             be modified as follows: */
/*                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR) */
#line 927 "dlarrv.f"
				} else {
/*                             The current LAMBDA is on the right of the true */
/*                             eigenvalue */
#line 930 "dlarrv.f"
				    right = lambda;
/*                             See comment about assuming the error estimate is */
/*                             correct above. */
/*                              LEFT = MIN(LEFT, LAMBDA + RQCORR) */
#line 934 "dlarrv.f"
				}
#line 935 "dlarrv.f"
				work[windex] = (right + left) * .5;
/*                          Take RQCORR since it has the correct sign and */
/*                          improves the iterate reasonably */
#line 939 "dlarrv.f"
				lambda += rqcorr;
/*                          Update width of error interval */
#line 941 "dlarrv.f"
				werr[windex] = (right - left) * .5;
#line 943 "dlarrv.f"
			    } else {
#line 944 "dlarrv.f"
				needbs = TRUE_;
#line 945 "dlarrv.f"
			    }
#line 946 "dlarrv.f"
			    if (right - left < rqtol * abs(lambda)) {
/*                             The eigenvalue is computed to bisection accuracy */
/*                             compute eigenvector and stop */
#line 949 "dlarrv.f"
				usedbs = TRUE_;
#line 950 "dlarrv.f"
				goto L120;
#line 951 "dlarrv.f"
			    } else if (iter < 10) {
#line 952 "dlarrv.f"
				goto L120;
#line 953 "dlarrv.f"
			    } else if (iter == 10) {
#line 954 "dlarrv.f"
				needbs = TRUE_;
#line 955 "dlarrv.f"
				goto L120;
#line 956 "dlarrv.f"
			    } else {
#line 957 "dlarrv.f"
				*info = 5;
#line 958 "dlarrv.f"
				return 0;
#line 959 "dlarrv.f"
			    }
#line 960 "dlarrv.f"
			} else {
#line 961 "dlarrv.f"
			    stp2ii = FALSE_;
#line 962 "dlarrv.f"
			    if (usedrq && usedbs && bstres <= resid) {
#line 964 "dlarrv.f"
				lambda = bstw;
#line 965 "dlarrv.f"
				stp2ii = TRUE_;
#line 966 "dlarrv.f"
			    }
#line 967 "dlarrv.f"
			    if (stp2ii) {
/*                          improve error angle by second step */
#line 969 "dlarrv.f"
				L__1 = ! usedbs;
#line 969 "dlarrv.f"
				dlar1v_(&in, &c__1, &in, &lambda, &d__[ibegin]
					, &l[ibegin], &work[indld + ibegin - 
					1], &work[indlld + ibegin - 1], 
					pivmin, &gaptol, &z__[ibegin + windex 
					* z_dim1], &L__1, &negcnt, &ztz, &
					mingma, &iwork[iindr + windex], &
					isuppz[(windex << 1) - 1], &nrminv, &
					resid, &rqcorr, &work[indwrk]);
#line 978 "dlarrv.f"
			    }
#line 979 "dlarrv.f"
			    work[windex] = lambda;
#line 980 "dlarrv.f"
			}

/*                    Compute FP-vector support w.r.t. whole matrix */

#line 984 "dlarrv.f"
			isuppz[(windex << 1) - 1] += oldien;
#line 985 "dlarrv.f"
			isuppz[windex * 2] += oldien;
#line 986 "dlarrv.f"
			zfrom = isuppz[(windex << 1) - 1];
#line 987 "dlarrv.f"
			zto = isuppz[windex * 2];
#line 988 "dlarrv.f"
			isupmn += oldien;
#line 989 "dlarrv.f"
			isupmx += oldien;
/*                    Ensure vector is ok if support in the RQI has changed */
#line 991 "dlarrv.f"
			if (isupmn < zfrom) {
#line 992 "dlarrv.f"
			    i__4 = zfrom - 1;
#line 992 "dlarrv.f"
			    for (ii = isupmn; ii <= i__4; ++ii) {
#line 993 "dlarrv.f"
				z__[ii + windex * z_dim1] = 0.;
#line 994 "dlarrv.f"
/* L122: */
#line 994 "dlarrv.f"
			    }
#line 995 "dlarrv.f"
			}
#line 996 "dlarrv.f"
			if (isupmx > zto) {
#line 997 "dlarrv.f"
			    i__4 = isupmx;
#line 997 "dlarrv.f"
			    for (ii = zto + 1; ii <= i__4; ++ii) {
#line 998 "dlarrv.f"
				z__[ii + windex * z_dim1] = 0.;
#line 999 "dlarrv.f"
/* L123: */
#line 999 "dlarrv.f"
			    }
#line 1000 "dlarrv.f"
			}
#line 1001 "dlarrv.f"
			i__4 = zto - zfrom + 1;
#line 1001 "dlarrv.f"
			dscal_(&i__4, &nrminv, &z__[zfrom + windex * z_dim1], 
				&c__1);
#line 1003 "dlarrv.f"
L125:
/*                    Update W */
#line 1005 "dlarrv.f"
			w[windex] = lambda + sigma;
/*                    Recompute the gaps on the left and right */
/*                    But only allow them to become larger and not */
/*                    smaller (which can only happen through "bad" */
/*                    cancellation and doesn't reflect the theory */
/*                    where the initial gaps are underestimated due */
/*                    to WERR being too crude.) */
#line 1012 "dlarrv.f"
			if (! eskip) {
#line 1013 "dlarrv.f"
			    if (k > 1) {
/* Computing MAX */
#line 1014 "dlarrv.f"
				d__1 = wgap[windmn], d__2 = w[windex] - werr[
					windex] - w[windmn] - werr[windmn];
#line 1014 "dlarrv.f"
				wgap[windmn] = max(d__1,d__2);
#line 1017 "dlarrv.f"
			    }
#line 1018 "dlarrv.f"
			    if (windex < wend) {
/* Computing MAX */
#line 1019 "dlarrv.f"
				d__1 = savgap, d__2 = w[windpl] - werr[windpl]
					 - w[windex] - werr[windex];
#line 1019 "dlarrv.f"
				wgap[windex] = max(d__1,d__2);
#line 1022 "dlarrv.f"
			    }
#line 1023 "dlarrv.f"
			}
#line 1024 "dlarrv.f"
			++idone;
#line 1025 "dlarrv.f"
		    }
/*                 here ends the code for the current child */

#line 1028 "dlarrv.f"
L139:
/*                 Proceed to any remaining child nodes */
#line 1030 "dlarrv.f"
		    newfst = j + 1;
#line 1031 "dlarrv.f"
L140:
#line 1031 "dlarrv.f"
		    ;
#line 1031 "dlarrv.f"
		}
#line 1032 "dlarrv.f"
/* L150: */
#line 1032 "dlarrv.f"
	    }
#line 1033 "dlarrv.f"
	    ++ndepth;
#line 1034 "dlarrv.f"
	    goto L40;
#line 1035 "dlarrv.f"
	}
#line 1036 "dlarrv.f"
	ibegin = iend + 1;
#line 1037 "dlarrv.f"
	wbegin = wend + 1;
#line 1038 "dlarrv.f"
L170:
#line 1038 "dlarrv.f"
	;
#line 1038 "dlarrv.f"
    }

#line 1041 "dlarrv.f"
    return 0;

/*     End of DLARRV */

} /* dlarrv_ */

