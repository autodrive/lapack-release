#line 1 "zlarrv.f"
/* zlarrv.f -- translated by f2c (version 20100827).
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

#line 1 "zlarrv.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b28 = 0.;

/* > \brief \b ZLARRV computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenv
alues of L D LT. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARRV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarrv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarrv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarrv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARRV( N, VL, VU, D, L, PIVMIN, */
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
/*       COMPLEX*16        Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARRV computes the eigenvectors of the tridiagonal matrix */
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
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          Lower and upper bounds of the interval that contains the desired */
/* >          eigenvalues. VL < VU. Needed to compute gaps on the left or right */
/* >          end of the extremal eigenvalues in the desired RANGE. */
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
/* >          (if the matrix is not splitted.) At the end of each block */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ, max(1,M) ) */
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
/* >          > 0:  A problem occured in ZLARRV. */
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

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int zlarrv_(integer *n, doublereal *vl, doublereal *vu, 
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
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nclus, zfrom;
    static doublereal rqtol;
    static integer iindc1, iindc2, indin1, indin2;
    static logical stp2ii;
    extern /* Subroutine */ int zlar1v_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublecomplex *, 
	    logical *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
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
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static logical usedbs;
    static integer iindwk, offset;
    static doublereal gaptol;
    static integer newcls, oldfst, indwrk, windex, oldlst;
    static logical usedrq;
    static integer newfst, newftt, parity, windmn, windpl, isupmn, newlst, 
	    zusedl;
    static doublereal bstres;
    static integer newsiz, zusedu, zusedw;
    static doublereal nrminv;
    static logical tryrqc;
    static integer isupmx;
    static doublereal rqcorr;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
/*     The first N entries of WORK are reserved for the eigenvalues */
#line 348 "zlarrv.f"
    /* Parameter adjustments */
#line 348 "zlarrv.f"
    --d__;
#line 348 "zlarrv.f"
    --l;
#line 348 "zlarrv.f"
    --isplit;
#line 348 "zlarrv.f"
    --w;
#line 348 "zlarrv.f"
    --werr;
#line 348 "zlarrv.f"
    --wgap;
#line 348 "zlarrv.f"
    --iblock;
#line 348 "zlarrv.f"
    --indexw;
#line 348 "zlarrv.f"
    --gers;
#line 348 "zlarrv.f"
    z_dim1 = *ldz;
#line 348 "zlarrv.f"
    z_offset = 1 + z_dim1;
#line 348 "zlarrv.f"
    z__ -= z_offset;
#line 348 "zlarrv.f"
    --isuppz;
#line 348 "zlarrv.f"
    --work;
#line 348 "zlarrv.f"
    --iwork;
#line 348 "zlarrv.f"

#line 348 "zlarrv.f"
    /* Function Body */
#line 348 "zlarrv.f"
    indld = *n + 1;
#line 349 "zlarrv.f"
    indlld = (*n << 1) + 1;
#line 350 "zlarrv.f"
    indin1 = *n * 3 + 1;
#line 351 "zlarrv.f"
    indin2 = (*n << 2) + 1;
#line 352 "zlarrv.f"
    indwrk = *n * 5 + 1;
#line 353 "zlarrv.f"
    minwsize = *n * 12;
#line 355 "zlarrv.f"
    i__1 = minwsize;
#line 355 "zlarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 356 "zlarrv.f"
	work[i__] = 0.;
#line 357 "zlarrv.f"
/* L5: */
#line 357 "zlarrv.f"
    }
/*     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the */
/*     factorization used to compute the FP vector */
#line 361 "zlarrv.f"
    iindr = 0;
/*     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current */
/*     layer and the one above. */
#line 364 "zlarrv.f"
    iindc1 = *n;
#line 365 "zlarrv.f"
    iindc2 = *n << 1;
#line 366 "zlarrv.f"
    iindwk = *n * 3 + 1;
#line 368 "zlarrv.f"
    miniwsize = *n * 7;
#line 369 "zlarrv.f"
    i__1 = miniwsize;
#line 369 "zlarrv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 370 "zlarrv.f"
	iwork[i__] = 0;
#line 371 "zlarrv.f"
/* L10: */
#line 371 "zlarrv.f"
    }
#line 373 "zlarrv.f"
    zusedl = 1;
#line 374 "zlarrv.f"
    if (*dol > 1) {
/*        Set lower bound for use of Z */
#line 376 "zlarrv.f"
	zusedl = *dol - 1;
#line 377 "zlarrv.f"
    }
#line 378 "zlarrv.f"
    zusedu = *m;
#line 379 "zlarrv.f"
    if (*dou < *m) {
/*        Set lower bound for use of Z */
#line 381 "zlarrv.f"
	zusedu = *dou + 1;
#line 382 "zlarrv.f"
    }
/*     The width of the part of Z that is used */
#line 384 "zlarrv.f"
    zusedw = zusedu - zusedl + 1;
#line 387 "zlarrv.f"
    zlaset_("Full", n, &zusedw, &c_b1, &c_b1, &z__[zusedl * z_dim1 + 1], ldz, 
	    (ftnlen)4);
#line 390 "zlarrv.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 391 "zlarrv.f"
    rqtol = eps * 2.;

/*     Set expert flags for standard code. */
#line 394 "zlarrv.f"
    tryrqc = TRUE_;
#line 396 "zlarrv.f"
    if (*dol == 1 && *dou == *m) {
#line 397 "zlarrv.f"
    } else {
/*        Only selected eigenpairs are computed. Since the other evalues */
/*        are not refined by RQ iteration, bisection has to compute to full */
/*        accuracy. */
#line 401 "zlarrv.f"
	*rtol1 = eps * 4.;
#line 402 "zlarrv.f"
	*rtol2 = eps * 4.;
#line 403 "zlarrv.f"
    }
/*     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the */
/*     desired eigenvalues. The support of the nonzero eigenvector */
/*     entries is contained in the interval IBEGIN:IEND. */
/*     Remark that if k eigenpairs are desired, then the eigenvectors */
/*     are stored in k contiguous columns of Z. */
/*     DONE is the number of eigenvectors already computed */
#line 412 "zlarrv.f"
    done = 0;
#line 413 "zlarrv.f"
    ibegin = 1;
#line 414 "zlarrv.f"
    wbegin = 1;
#line 415 "zlarrv.f"
    i__1 = iblock[*m];
#line 415 "zlarrv.f"
    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 416 "zlarrv.f"
	iend = isplit[jblk];
#line 417 "zlarrv.f"
	sigma = l[iend];
/*        Find the eigenvectors of the submatrix indexed IBEGIN */
/*        through IEND. */
#line 420 "zlarrv.f"
	wend = wbegin - 1;
#line 421 "zlarrv.f"
L15:
#line 422 "zlarrv.f"
	if (wend < *m) {
#line 423 "zlarrv.f"
	    if (iblock[wend + 1] == jblk) {
#line 424 "zlarrv.f"
		++wend;
#line 425 "zlarrv.f"
		goto L15;
#line 426 "zlarrv.f"
	    }
#line 427 "zlarrv.f"
	}
#line 428 "zlarrv.f"
	if (wend < wbegin) {
#line 429 "zlarrv.f"
	    ibegin = iend + 1;
#line 430 "zlarrv.f"
	    goto L170;
#line 431 "zlarrv.f"
	} else if (wend < *dol || wbegin > *dou) {
#line 432 "zlarrv.f"
	    ibegin = iend + 1;
#line 433 "zlarrv.f"
	    wbegin = wend + 1;
#line 434 "zlarrv.f"
	    goto L170;
#line 435 "zlarrv.f"
	}
/*        Find local spectral diameter of the block */
#line 438 "zlarrv.f"
	gl = gers[(ibegin << 1) - 1];
#line 439 "zlarrv.f"
	gu = gers[ibegin * 2];
#line 440 "zlarrv.f"
	i__2 = iend;
#line 440 "zlarrv.f"
	for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 441 "zlarrv.f"
	    d__1 = gers[(i__ << 1) - 1];
#line 441 "zlarrv.f"
	    gl = min(d__1,gl);
/* Computing MAX */
#line 442 "zlarrv.f"
	    d__1 = gers[i__ * 2];
#line 442 "zlarrv.f"
	    gu = max(d__1,gu);
#line 443 "zlarrv.f"
/* L20: */
#line 443 "zlarrv.f"
	}
#line 444 "zlarrv.f"
	spdiam = gu - gl;
/*        OLDIEN is the last index of the previous block */
#line 447 "zlarrv.f"
	oldien = ibegin - 1;
/*        Calculate the size of the current block */
#line 449 "zlarrv.f"
	in = iend - ibegin + 1;
/*        The number of eigenvalues in the current block */
#line 451 "zlarrv.f"
	im = wend - wbegin + 1;
/*        This is for a 1x1 block */
#line 454 "zlarrv.f"
	if (ibegin == iend) {
#line 455 "zlarrv.f"
	    ++done;
#line 456 "zlarrv.f"
	    i__2 = ibegin + wbegin * z_dim1;
#line 456 "zlarrv.f"
	    z__[i__2].r = 1., z__[i__2].i = 0.;
#line 457 "zlarrv.f"
	    isuppz[(wbegin << 1) - 1] = ibegin;
#line 458 "zlarrv.f"
	    isuppz[wbegin * 2] = ibegin;
#line 459 "zlarrv.f"
	    w[wbegin] += sigma;
#line 460 "zlarrv.f"
	    work[wbegin] = w[wbegin];
#line 461 "zlarrv.f"
	    ibegin = iend + 1;
#line 462 "zlarrv.f"
	    ++wbegin;
#line 463 "zlarrv.f"
	    goto L170;
#line 464 "zlarrv.f"
	}
/*        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND) */
/*        Note that these can be approximations, in this case, the corresp. */
/*        entries of WERR give the size of the uncertainty interval. */
/*        The eigenvalue approximations will be refined when necessary as */
/*        high relative accuracy is required for the computation of the */
/*        corresponding eigenvectors. */
#line 472 "zlarrv.f"
	dcopy_(&im, &w[wbegin], &c__1, &work[wbegin], &c__1);
/*        We store in W the eigenvalue approximations w.r.t. the original */
/*        matrix T. */
#line 477 "zlarrv.f"
	i__2 = im;
#line 477 "zlarrv.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 478 "zlarrv.f"
	    w[wbegin + i__ - 1] += sigma;
#line 479 "zlarrv.f"
/* L30: */
#line 479 "zlarrv.f"
	}
/*        NDEPTH is the current depth of the representation tree */
#line 483 "zlarrv.f"
	ndepth = 0;
/*        PARITY is either 1 or 0 */
#line 485 "zlarrv.f"
	parity = 1;
/*        NCLUS is the number of clusters for the next level of the */
/*        representation tree, we start with NCLUS = 1 for the root */
#line 488 "zlarrv.f"
	nclus = 1;
#line 489 "zlarrv.f"
	iwork[iindc1 + 1] = 1;
#line 490 "zlarrv.f"
	iwork[iindc1 + 2] = im;
/*        IDONE is the number of eigenvectors already computed in the current */
/*        block */
#line 494 "zlarrv.f"
	idone = 0;
/*        loop while( IDONE.LT.IM ) */
/*        generate the representation tree for the current block and */
/*        compute the eigenvectors */
#line 498 "zlarrv.f"
L40:
#line 499 "zlarrv.f"
	if (idone < im) {
/*           This is a crude protection against infinitely deep trees */
#line 501 "zlarrv.f"
	    if (ndepth > *m) {
#line 502 "zlarrv.f"
		*info = -2;
#line 503 "zlarrv.f"
		return 0;
#line 504 "zlarrv.f"
	    }
/*           breadth first processing of the current level of the representation */
/*           tree: OLDNCL = number of clusters on current level */
#line 507 "zlarrv.f"
	    oldncl = nclus;
/*           reset NCLUS to count the number of child clusters */
#line 509 "zlarrv.f"
	    nclus = 0;

#line 511 "zlarrv.f"
	    parity = 1 - parity;
#line 512 "zlarrv.f"
	    if (parity == 0) {
#line 513 "zlarrv.f"
		oldcls = iindc1;
#line 514 "zlarrv.f"
		newcls = iindc2;
#line 515 "zlarrv.f"
	    } else {
#line 516 "zlarrv.f"
		oldcls = iindc2;
#line 517 "zlarrv.f"
		newcls = iindc1;
#line 518 "zlarrv.f"
	    }
/*           Process the clusters on the current level */
#line 520 "zlarrv.f"
	    i__2 = oldncl;
#line 520 "zlarrv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 521 "zlarrv.f"
		j = oldcls + (i__ << 1);
/*              OLDFST, OLDLST = first, last index of current cluster. */
/*                               cluster indices start with 1 and are relative */
/*                               to WBEGIN when accessing W, WGAP, WERR, Z */
#line 525 "zlarrv.f"
		oldfst = iwork[j - 1];
#line 526 "zlarrv.f"
		oldlst = iwork[j];
#line 527 "zlarrv.f"
		if (ndepth > 0) {
/*                 Retrieve relatively robust representation (RRR) of cluster */
/*                 that has been computed at the previous level */
/*                 The RRR is stored in Z and overwritten once the eigenvectors */
/*                 have been computed or when the cluster is refined */
#line 533 "zlarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Get representation from location of the leftmost evalue */
/*                    of the cluster */
#line 536 "zlarrv.f"
			j = wbegin + oldfst - 1;
#line 537 "zlarrv.f"
		    } else {
#line 538 "zlarrv.f"
			if (wbegin + oldfst - 1 < *dol) {
/*                       Get representation from the left end of Z array */
#line 540 "zlarrv.f"
			    j = *dol - 1;
#line 541 "zlarrv.f"
			} else if (wbegin + oldfst - 1 > *dou) {
/*                       Get representation from the right end of Z array */
#line 543 "zlarrv.f"
			    j = *dou;
#line 544 "zlarrv.f"
			} else {
#line 545 "zlarrv.f"
			    j = wbegin + oldfst - 1;
#line 546 "zlarrv.f"
			}
#line 547 "zlarrv.f"
		    }
#line 548 "zlarrv.f"
		    i__3 = in - 1;
#line 548 "zlarrv.f"
		    for (k = 1; k <= i__3; ++k) {
#line 549 "zlarrv.f"
			i__4 = ibegin + k - 1 + j * z_dim1;
#line 549 "zlarrv.f"
			d__[ibegin + k - 1] = z__[i__4].r;
#line 551 "zlarrv.f"
			i__4 = ibegin + k - 1 + (j + 1) * z_dim1;
#line 551 "zlarrv.f"
			l[ibegin + k - 1] = z__[i__4].r;
#line 553 "zlarrv.f"
/* L45: */
#line 553 "zlarrv.f"
		    }
#line 554 "zlarrv.f"
		    i__3 = iend + j * z_dim1;
#line 554 "zlarrv.f"
		    d__[iend] = z__[i__3].r;
#line 555 "zlarrv.f"
		    i__3 = iend + (j + 1) * z_dim1;
#line 555 "zlarrv.f"
		    sigma = z__[i__3].r;
/*                 Set the corresponding entries in Z to zero */
#line 558 "zlarrv.f"
		    zlaset_("Full", &in, &c__2, &c_b1, &c_b1, &z__[ibegin + j 
			    * z_dim1], ldz, (ftnlen)4);
#line 560 "zlarrv.f"
		}
/*              Compute DL and DLL of current RRR */
#line 563 "zlarrv.f"
		i__3 = iend - 1;
#line 563 "zlarrv.f"
		for (j = ibegin; j <= i__3; ++j) {
#line 564 "zlarrv.f"
		    tmp = d__[j] * l[j];
#line 565 "zlarrv.f"
		    work[indld - 1 + j] = tmp;
#line 566 "zlarrv.f"
		    work[indlld - 1 + j] = tmp * l[j];
#line 567 "zlarrv.f"
/* L50: */
#line 567 "zlarrv.f"
		}
#line 569 "zlarrv.f"
		if (ndepth > 0) {
/*                 P and Q are index of the first and last eigenvalue to compute */
/*                 within the current block */
#line 572 "zlarrv.f"
		    p = indexw[wbegin - 1 + oldfst];
#line 573 "zlarrv.f"
		    q = indexw[wbegin - 1 + oldlst];
/*                 Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET */
/*                 through the Q-OFFSET elements of these arrays are to be used. */
/*                  OFFSET = P-OLDFST */
#line 577 "zlarrv.f"
		    offset = indexw[wbegin] - 1;
/*                 perform limited bisection (if necessary) to get approximate */
/*                 eigenvalues to the precision needed. */
#line 580 "zlarrv.f"
		    dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin - 1], &p,
			     &q, rtol1, rtol2, &offset, &work[wbegin], &wgap[
			    wbegin], &werr[wbegin], &work[indwrk], &iwork[
			    iindwk], pivmin, &spdiam, &in, &iinfo);
#line 586 "zlarrv.f"
		    if (iinfo != 0) {
#line 587 "zlarrv.f"
			*info = -1;
#line 588 "zlarrv.f"
			return 0;
#line 589 "zlarrv.f"
		    }
/*                 We also recompute the extremal gaps. W holds all eigenvalues */
/*                 of the unshifted matrix and must be used for computation */
/*                 of WGAP, the entries of WORK might stem from RRRs with */
/*                 different shifts. The gaps from WBEGIN-1+OLDFST to */
/*                 WBEGIN-1+OLDLST are correctly computed in DLARRB. */
/*                 However, we only allow the gaps to become greater since */
/*                 this is what should happen when we decrease WERR */
#line 597 "zlarrv.f"
		    if (oldfst > 1) {
/* Computing MAX */
#line 598 "zlarrv.f"
			d__1 = wgap[wbegin + oldfst - 2], d__2 = w[wbegin + 
				oldfst - 1] - werr[wbegin + oldfst - 1] - w[
				wbegin + oldfst - 2] - werr[wbegin + oldfst - 
				2];
#line 598 "zlarrv.f"
			wgap[wbegin + oldfst - 2] = max(d__1,d__2);
#line 602 "zlarrv.f"
		    }
#line 603 "zlarrv.f"
		    if (wbegin + oldlst - 1 < wend) {
/* Computing MAX */
#line 604 "zlarrv.f"
			d__1 = wgap[wbegin + oldlst - 1], d__2 = w[wbegin + 
				oldlst] - werr[wbegin + oldlst] - w[wbegin + 
				oldlst - 1] - werr[wbegin + oldlst - 1];
#line 604 "zlarrv.f"
			wgap[wbegin + oldlst - 1] = max(d__1,d__2);
#line 608 "zlarrv.f"
		    }
/*                 Each time the eigenvalues in WORK get refined, we store */
/*                 the newly found approximation with all shifts applied in W */
#line 611 "zlarrv.f"
		    i__3 = oldlst;
#line 611 "zlarrv.f"
		    for (j = oldfst; j <= i__3; ++j) {
#line 612 "zlarrv.f"
			w[wbegin + j - 1] = work[wbegin + j - 1] + sigma;
#line 613 "zlarrv.f"
/* L53: */
#line 613 "zlarrv.f"
		    }
#line 614 "zlarrv.f"
		}
/*              Process the current node. */
#line 617 "zlarrv.f"
		newfst = oldfst;
#line 618 "zlarrv.f"
		i__3 = oldlst;
#line 618 "zlarrv.f"
		for (j = oldfst; j <= i__3; ++j) {
#line 619 "zlarrv.f"
		    if (j == oldlst) {
/*                    we are at the right end of the cluster, this is also the */
/*                    boundary of the child cluster */
#line 622 "zlarrv.f"
			newlst = j;
#line 623 "zlarrv.f"
		    } else if (wgap[wbegin + j - 1] >= *minrgp * (d__1 = work[
			    wbegin + j - 1], abs(d__1))) {
/*                    the right relative gap is big enough, the child cluster */
/*                    (NEWFST,..,NEWLST) is well separated from the following */
#line 627 "zlarrv.f"
			newlst = j;
#line 628 "zlarrv.f"
		    } else {
/*                    inside a child cluster, the relative gap is not */
/*                    big enough. */
#line 631 "zlarrv.f"
			goto L140;
#line 632 "zlarrv.f"
		    }
/*                 Compute size of child cluster found */
#line 635 "zlarrv.f"
		    newsiz = newlst - newfst + 1;
/*                 NEWFTT is the place in Z where the new RRR or the computed */
/*                 eigenvector is to be stored */
#line 639 "zlarrv.f"
		    if (*dol == 1 && *dou == *m) {
/*                    Store representation at location of the leftmost evalue */
/*                    of the cluster */
#line 642 "zlarrv.f"
			newftt = wbegin + newfst - 1;
#line 643 "zlarrv.f"
		    } else {
#line 644 "zlarrv.f"
			if (wbegin + newfst - 1 < *dol) {
/*                       Store representation at the left end of Z array */
#line 646 "zlarrv.f"
			    newftt = *dol - 1;
#line 647 "zlarrv.f"
			} else if (wbegin + newfst - 1 > *dou) {
/*                       Store representation at the right end of Z array */
#line 649 "zlarrv.f"
			    newftt = *dou;
#line 650 "zlarrv.f"
			} else {
#line 651 "zlarrv.f"
			    newftt = wbegin + newfst - 1;
#line 652 "zlarrv.f"
			}
#line 653 "zlarrv.f"
		    }
#line 655 "zlarrv.f"
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
#line 670 "zlarrv.f"
			if (newfst == 1) {
/* Computing MAX */
#line 671 "zlarrv.f"
			    d__1 = 0., d__2 = w[wbegin] - werr[wbegin] - *vl;
#line 671 "zlarrv.f"
			    lgap = max(d__1,d__2);
#line 673 "zlarrv.f"
			} else {
#line 674 "zlarrv.f"
			    lgap = wgap[wbegin + newfst - 2];
#line 675 "zlarrv.f"
			}
#line 676 "zlarrv.f"
			rgap = wgap[wbegin + newlst - 1];

/*                    Compute left- and rightmost eigenvalue of child */
/*                    to high precision in order to shift as close */
/*                    as possible and obtain as large relative gaps */
/*                    as possible */

#line 683 "zlarrv.f"
			for (k = 1; k <= 2; ++k) {
#line 684 "zlarrv.f"
			    if (k == 1) {
#line 685 "zlarrv.f"
				p = indexw[wbegin - 1 + newfst];
#line 686 "zlarrv.f"
			    } else {
#line 687 "zlarrv.f"
				p = indexw[wbegin - 1 + newlst];
#line 688 "zlarrv.f"
			    }
#line 689 "zlarrv.f"
			    offset = indexw[wbegin] - 1;
#line 690 "zlarrv.f"
			    dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &p, &p, &rqtol, &rqtol, &offset, &
				    work[wbegin], &wgap[wbegin], &werr[wbegin]
				    , &work[indwrk], &iwork[iindwk], pivmin, &
				    spdiam, &in, &iinfo);
#line 697 "zlarrv.f"
/* L55: */
#line 697 "zlarrv.f"
			}

#line 699 "zlarrv.f"
			if (wbegin + newlst - 1 < *dol || wbegin + newfst - 1 
				> *dou) {
/*                       if the cluster contains no desired eigenvalues */
/*                       skip the computation of that branch of the rep. tree */

/*                       We could skip before the refinement of the extremal */
/*                       eigenvalues of the child, but then the representation */
/*                       tree could be different from the one when nothing is */
/*                       skipped. For this reason we skip at this place. */
#line 708 "zlarrv.f"
			    idone = idone + newlst - newfst + 1;
#line 709 "zlarrv.f"
			    goto L139;
#line 710 "zlarrv.f"
			}

/*                    Compute RRR of child cluster. */
/*                    Note that the new RRR is stored in Z */

/*                    DLARRF needs LWORK = 2*N */
#line 716 "zlarrv.f"
			dlarrf_(&in, &d__[ibegin], &l[ibegin], &work[indld + 
				ibegin - 1], &newfst, &newlst, &work[wbegin], 
				&wgap[wbegin], &werr[wbegin], &spdiam, &lgap, 
				&rgap, pivmin, &tau, &work[indin1], &work[
				indin2], &work[indwrk], &iinfo);
/*                    In the complex case, DLARRF cannot write */
/*                    the new RRR directly into Z and needs an intermediate */
/*                    workspace */
#line 726 "zlarrv.f"
			i__4 = in - 1;
#line 726 "zlarrv.f"
			for (k = 1; k <= i__4; ++k) {
#line 727 "zlarrv.f"
			    i__5 = ibegin + k - 1 + newftt * z_dim1;
#line 727 "zlarrv.f"
			    i__6 = indin1 + k - 1;
#line 727 "zlarrv.f"
			    z__1.r = work[i__6], z__1.i = 0.;
#line 727 "zlarrv.f"
			    z__[i__5].r = z__1.r, z__[i__5].i = z__1.i;
#line 729 "zlarrv.f"
			    i__5 = ibegin + k - 1 + (newftt + 1) * z_dim1;
#line 729 "zlarrv.f"
			    i__6 = indin2 + k - 1;
#line 729 "zlarrv.f"
			    z__1.r = work[i__6], z__1.i = 0.;
#line 729 "zlarrv.f"
			    z__[i__5].r = z__1.r, z__[i__5].i = z__1.i;
#line 731 "zlarrv.f"
/* L56: */
#line 731 "zlarrv.f"
			}
#line 732 "zlarrv.f"
			i__4 = iend + newftt * z_dim1;
#line 732 "zlarrv.f"
			i__5 = indin1 + in - 1;
#line 732 "zlarrv.f"
			z__1.r = work[i__5], z__1.i = 0.;
#line 732 "zlarrv.f"
			z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 734 "zlarrv.f"
			if (iinfo == 0) {
/*                       a new RRR for the cluster was found by DLARRF */
/*                       update shift and store it */
#line 737 "zlarrv.f"
			    ssigma = sigma + tau;
#line 738 "zlarrv.f"
			    i__4 = iend + (newftt + 1) * z_dim1;
#line 738 "zlarrv.f"
			    z__1.r = ssigma, z__1.i = 0.;
#line 738 "zlarrv.f"
			    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
/*                       WORK() are the midpoints and WERR() the semi-width */
/*                       Note that the entries in W are unchanged. */
#line 741 "zlarrv.f"
			    i__4 = newlst;
#line 741 "zlarrv.f"
			    for (k = newfst; k <= i__4; ++k) {
#line 742 "zlarrv.f"
				fudge = eps * 3. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
#line 744 "zlarrv.f"
				work[wbegin + k - 1] -= tau;
#line 746 "zlarrv.f"
				fudge += eps * 4. * (d__1 = work[wbegin + k - 
					1], abs(d__1));
/*                          Fudge errors */
#line 749 "zlarrv.f"
				werr[wbegin + k - 1] += fudge;
/*                          Gaps are not fudged. Provided that WERR is small */
/*                          when eigenvalues are close, a zero gap indicates */
/*                          that a new representation is needed for resolving */
/*                          the cluster. A fudge could lead to a wrong decision */
/*                          of judging eigenvalues 'separated' which in */
/*                          reality are not. This could have a negative impact */
/*                          on the orthogonality of the computed eigenvectors. */
#line 758 "zlarrv.f"
/* L116: */
#line 758 "zlarrv.f"
			    }
#line 760 "zlarrv.f"
			    ++nclus;
#line 761 "zlarrv.f"
			    k = newcls + (nclus << 1);
#line 762 "zlarrv.f"
			    iwork[k - 1] = newfst;
#line 763 "zlarrv.f"
			    iwork[k] = newlst;
#line 764 "zlarrv.f"
			} else {
#line 765 "zlarrv.f"
			    *info = -2;
#line 766 "zlarrv.f"
			    return 0;
#line 767 "zlarrv.f"
			}
#line 768 "zlarrv.f"
		    } else {

/*                    Compute eigenvector of singleton */

#line 772 "zlarrv.f"
			iter = 0;

#line 774 "zlarrv.f"
			tol = log((doublereal) in) * 4. * eps;

#line 776 "zlarrv.f"
			k = newfst;
#line 777 "zlarrv.f"
			windex = wbegin + k - 1;
/* Computing MAX */
#line 778 "zlarrv.f"
			i__4 = windex - 1;
#line 778 "zlarrv.f"
			windmn = max(i__4,1);
/* Computing MIN */
#line 779 "zlarrv.f"
			i__4 = windex + 1;
#line 779 "zlarrv.f"
			windpl = min(i__4,*m);
#line 780 "zlarrv.f"
			lambda = work[windex];
#line 781 "zlarrv.f"
			++done;
/*                    Check if eigenvector computation is to be skipped */
#line 783 "zlarrv.f"
			if (windex < *dol || windex > *dou) {
#line 785 "zlarrv.f"
			    eskip = TRUE_;
#line 786 "zlarrv.f"
			    goto L125;
#line 787 "zlarrv.f"
			} else {
#line 788 "zlarrv.f"
			    eskip = FALSE_;
#line 789 "zlarrv.f"
			}
#line 790 "zlarrv.f"
			left = work[windex] - werr[windex];
#line 791 "zlarrv.f"
			right = work[windex] + werr[windex];
#line 792 "zlarrv.f"
			indeig = indexw[windex];
/*                    Note that since we compute the eigenpairs for a child, */
/*                    all eigenvalue approximations are w.r.t the same shift. */
/*                    In this case, the entries in WORK should be used for */
/*                    computing the gaps since they exhibit even very small */
/*                    differences in the eigenvalues, as opposed to the */
/*                    entries in W which might "look" the same. */
#line 800 "zlarrv.f"
			if (k == 1) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VL, the formula */
/*                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA ) */
/*                       can lead to an overestimation of the left gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small left gap. */
/* Computing MAX */
#line 807 "zlarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 807 "zlarrv.f"
			    lgap = eps * max(d__1,d__2);
#line 808 "zlarrv.f"
			} else {
#line 809 "zlarrv.f"
			    lgap = wgap[windmn];
#line 810 "zlarrv.f"
			}
#line 811 "zlarrv.f"
			if (k == im) {
/*                       In the case RANGE='I' and with not much initial */
/*                       accuracy in LAMBDA and VU, the formula */
/*                       can lead to an overestimation of the right gap and */
/*                       thus to inadequately early RQI 'convergence'. */
/*                       Prevent this by forcing a small right gap. */
/* Computing MAX */
#line 817 "zlarrv.f"
			    d__1 = abs(left), d__2 = abs(right);
#line 817 "zlarrv.f"
			    rgap = eps * max(d__1,d__2);
#line 818 "zlarrv.f"
			} else {
#line 819 "zlarrv.f"
			    rgap = wgap[windex];
#line 820 "zlarrv.f"
			}
#line 821 "zlarrv.f"
			gap = min(lgap,rgap);
#line 822 "zlarrv.f"
			if (k == 1 || k == im) {
/*                       The eigenvector support can become wrong */
/*                       because significant entries could be cut off due to a */
/*                       large GAPTOL parameter in LAR1V. Prevent this. */
#line 826 "zlarrv.f"
			    gaptol = 0.;
#line 827 "zlarrv.f"
			} else {
#line 828 "zlarrv.f"
			    gaptol = gap * eps;
#line 829 "zlarrv.f"
			}
#line 830 "zlarrv.f"
			isupmn = in;
#line 831 "zlarrv.f"
			isupmx = 1;
/*                    Update WGAP so that it holds the minimum gap */
/*                    to the left or the right. This is crucial in the */
/*                    case where bisection is used to ensure that the */
/*                    eigenvalue is refined up to the required precision. */
/*                    The correct value is restored afterwards. */
#line 837 "zlarrv.f"
			savgap = wgap[windex];
#line 838 "zlarrv.f"
			wgap[windex] = gap;
/*                    We want to use the Rayleigh Quotient Correction */
/*                    as often as possible since it converges quadratically */
/*                    when we are close enough to the desired eigenvalue. */
/*                    However, the Rayleigh Quotient can have the wrong sign */
/*                    and lead us away from the desired eigenvalue. In this */
/*                    case, the best we can do is to use bisection. */
#line 845 "zlarrv.f"
			usedbs = FALSE_;
#line 846 "zlarrv.f"
			usedrq = FALSE_;
/*                    Bisection is initially turned off unless it is forced */
#line 848 "zlarrv.f"
			needbs = ! tryrqc;
#line 849 "zlarrv.f"
L120:
/*                    Check if bisection should be used to refine eigenvalue */
#line 851 "zlarrv.f"
			if (needbs) {
/*                       Take the bisection as new iterate */
#line 853 "zlarrv.f"
			    usedbs = TRUE_;
#line 854 "zlarrv.f"
			    itmp1 = iwork[iindr + windex];
#line 855 "zlarrv.f"
			    offset = indexw[wbegin] - 1;
#line 856 "zlarrv.f"
			    d__1 = eps * 2.;
#line 856 "zlarrv.f"
			    dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin 
				    - 1], &indeig, &indeig, &c_b28, &d__1, &
				    offset, &work[wbegin], &wgap[wbegin], &
				    werr[wbegin], &work[indwrk], &iwork[
				    iindwk], pivmin, &spdiam, &itmp1, &iinfo);
#line 863 "zlarrv.f"
			    if (iinfo != 0) {
#line 864 "zlarrv.f"
				*info = -3;
#line 865 "zlarrv.f"
				return 0;
#line 866 "zlarrv.f"
			    }
#line 867 "zlarrv.f"
			    lambda = work[windex];
/*                       Reset twist index from inaccurate LAMBDA to */
/*                       force computation of true MINGMA */
#line 870 "zlarrv.f"
			    iwork[iindr + windex] = 0;
#line 871 "zlarrv.f"
			}
/*                    Given LAMBDA, compute the eigenvector. */
#line 873 "zlarrv.f"
			L__1 = ! usedbs;
#line 873 "zlarrv.f"
			zlar1v_(&in, &c__1, &in, &lambda, &d__[ibegin], &l[
				ibegin], &work[indld + ibegin - 1], &work[
				indlld + ibegin - 1], pivmin, &gaptol, &z__[
				ibegin + windex * z_dim1], &L__1, &negcnt, &
				ztz, &mingma, &iwork[iindr + windex], &isuppz[
				(windex << 1) - 1], &nrminv, &resid, &rqcorr, 
				&work[indwrk]);
#line 880 "zlarrv.f"
			if (iter == 0) {
#line 881 "zlarrv.f"
			    bstres = resid;
#line 882 "zlarrv.f"
			    bstw = lambda;
#line 883 "zlarrv.f"
			} else if (resid < bstres) {
#line 884 "zlarrv.f"
			    bstres = resid;
#line 885 "zlarrv.f"
			    bstw = lambda;
#line 886 "zlarrv.f"
			}
/* Computing MIN */
#line 887 "zlarrv.f"
			i__4 = isupmn, i__5 = isuppz[(windex << 1) - 1];
#line 887 "zlarrv.f"
			isupmn = min(i__4,i__5);
/* Computing MAX */
#line 888 "zlarrv.f"
			i__4 = isupmx, i__5 = isuppz[windex * 2];
#line 888 "zlarrv.f"
			isupmx = max(i__4,i__5);
#line 889 "zlarrv.f"
			++iter;
/*                    sin alpha <= |resid|/gap */
/*                    Note that both the residual and the gap are */
/*                    proportional to the matrix, so ||T|| doesn't play */
/*                    a role in the quotient */

/*                    Convergence test for Rayleigh-Quotient iteration */
/*                    (omitted when Bisection has been used) */

#line 900 "zlarrv.f"
			if (resid > tol * gap && abs(rqcorr) > rqtol * abs(
				lambda) && ! usedbs) {
/*                       We need to check that the RQCORR update doesn't */
/*                       move the eigenvalue away from the desired one and */
/*                       towards a neighbor. -> protection with bisection */
#line 906 "zlarrv.f"
			    if (indeig <= negcnt) {
/*                          The wanted eigenvalue lies to the left */
#line 908 "zlarrv.f"
				sgndef = -1.;
#line 909 "zlarrv.f"
			    } else {
/*                          The wanted eigenvalue lies to the right */
#line 911 "zlarrv.f"
				sgndef = 1.;
#line 912 "zlarrv.f"
			    }
/*                       We only use the RQCORR if it improves the */
/*                       the iterate reasonably. */
#line 915 "zlarrv.f"
			    if (rqcorr * sgndef >= 0. && lambda + rqcorr <= 
				    right && lambda + rqcorr >= left) {
#line 919 "zlarrv.f"
				usedrq = TRUE_;
/*                          Store new midpoint of bisection interval in WORK */
#line 921 "zlarrv.f"
				if (sgndef == 1.) {
/*                             The current LAMBDA is on the left of the true */
/*                             eigenvalue */
#line 924 "zlarrv.f"
				    left = lambda;
/*                             We prefer to assume that the error estimate */
/*                             is correct. We could make the interval not */
/*                             as a bracket but to be modified if the RQCORR */
/*                             chooses to. In this case, the RIGHT side should */
/*                             be modified as follows: */
/*                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR) */
#line 931 "zlarrv.f"
				} else {
/*                             The current LAMBDA is on the right of the true */
/*                             eigenvalue */
#line 934 "zlarrv.f"
				    right = lambda;
/*                             See comment about assuming the error estimate is */
/*                             correct above. */
/*                              LEFT = MIN(LEFT, LAMBDA + RQCORR) */
#line 938 "zlarrv.f"
				}
#line 939 "zlarrv.f"
				work[windex] = (right + left) * .5;
/*                          Take RQCORR since it has the correct sign and */
/*                          improves the iterate reasonably */
#line 943 "zlarrv.f"
				lambda += rqcorr;
/*                          Update width of error interval */
#line 945 "zlarrv.f"
				werr[windex] = (right - left) * .5;
#line 947 "zlarrv.f"
			    } else {
#line 948 "zlarrv.f"
				needbs = TRUE_;
#line 949 "zlarrv.f"
			    }
#line 950 "zlarrv.f"
			    if (right - left < rqtol * abs(lambda)) {
/*                             The eigenvalue is computed to bisection accuracy */
/*                             compute eigenvector and stop */
#line 953 "zlarrv.f"
				usedbs = TRUE_;
#line 954 "zlarrv.f"
				goto L120;
#line 955 "zlarrv.f"
			    } else if (iter < 10) {
#line 956 "zlarrv.f"
				goto L120;
#line 957 "zlarrv.f"
			    } else if (iter == 10) {
#line 958 "zlarrv.f"
				needbs = TRUE_;
#line 959 "zlarrv.f"
				goto L120;
#line 960 "zlarrv.f"
			    } else {
#line 961 "zlarrv.f"
				*info = 5;
#line 962 "zlarrv.f"
				return 0;
#line 963 "zlarrv.f"
			    }
#line 964 "zlarrv.f"
			} else {
#line 965 "zlarrv.f"
			    stp2ii = FALSE_;
#line 966 "zlarrv.f"
			    if (usedrq && usedbs && bstres <= resid) {
#line 968 "zlarrv.f"
				lambda = bstw;
#line 969 "zlarrv.f"
				stp2ii = TRUE_;
#line 970 "zlarrv.f"
			    }
#line 971 "zlarrv.f"
			    if (stp2ii) {
/*                          improve error angle by second step */
#line 973 "zlarrv.f"
				L__1 = ! usedbs;
#line 973 "zlarrv.f"
				zlar1v_(&in, &c__1, &in, &lambda, &d__[ibegin]
					, &l[ibegin], &work[indld + ibegin - 
					1], &work[indlld + ibegin - 1], 
					pivmin, &gaptol, &z__[ibegin + windex 
					* z_dim1], &L__1, &negcnt, &ztz, &
					mingma, &iwork[iindr + windex], &
					isuppz[(windex << 1) - 1], &nrminv, &
					resid, &rqcorr, &work[indwrk]);
#line 982 "zlarrv.f"
			    }
#line 983 "zlarrv.f"
			    work[windex] = lambda;
#line 984 "zlarrv.f"
			}

/*                    Compute FP-vector support w.r.t. whole matrix */

#line 988 "zlarrv.f"
			isuppz[(windex << 1) - 1] += oldien;
#line 989 "zlarrv.f"
			isuppz[windex * 2] += oldien;
#line 990 "zlarrv.f"
			zfrom = isuppz[(windex << 1) - 1];
#line 991 "zlarrv.f"
			zto = isuppz[windex * 2];
#line 992 "zlarrv.f"
			isupmn += oldien;
#line 993 "zlarrv.f"
			isupmx += oldien;
/*                    Ensure vector is ok if support in the RQI has changed */
#line 995 "zlarrv.f"
			if (isupmn < zfrom) {
#line 996 "zlarrv.f"
			    i__4 = zfrom - 1;
#line 996 "zlarrv.f"
			    for (ii = isupmn; ii <= i__4; ++ii) {
#line 997 "zlarrv.f"
				i__5 = ii + windex * z_dim1;
#line 997 "zlarrv.f"
				z__[i__5].r = 0., z__[i__5].i = 0.;
#line 998 "zlarrv.f"
/* L122: */
#line 998 "zlarrv.f"
			    }
#line 999 "zlarrv.f"
			}
#line 1000 "zlarrv.f"
			if (isupmx > zto) {
#line 1001 "zlarrv.f"
			    i__4 = isupmx;
#line 1001 "zlarrv.f"
			    for (ii = zto + 1; ii <= i__4; ++ii) {
#line 1002 "zlarrv.f"
				i__5 = ii + windex * z_dim1;
#line 1002 "zlarrv.f"
				z__[i__5].r = 0., z__[i__5].i = 0.;
#line 1003 "zlarrv.f"
/* L123: */
#line 1003 "zlarrv.f"
			    }
#line 1004 "zlarrv.f"
			}
#line 1005 "zlarrv.f"
			i__4 = zto - zfrom + 1;
#line 1005 "zlarrv.f"
			zdscal_(&i__4, &nrminv, &z__[zfrom + windex * z_dim1],
				 &c__1);
#line 1007 "zlarrv.f"
L125:
/*                    Update W */
#line 1009 "zlarrv.f"
			w[windex] = lambda + sigma;
/*                    Recompute the gaps on the left and right */
/*                    But only allow them to become larger and not */
/*                    smaller (which can only happen through "bad" */
/*                    cancellation and doesn't reflect the theory */
/*                    where the initial gaps are underestimated due */
/*                    to WERR being too crude.) */
#line 1016 "zlarrv.f"
			if (! eskip) {
#line 1017 "zlarrv.f"
			    if (k > 1) {
/* Computing MAX */
#line 1018 "zlarrv.f"
				d__1 = wgap[windmn], d__2 = w[windex] - werr[
					windex] - w[windmn] - werr[windmn];
#line 1018 "zlarrv.f"
				wgap[windmn] = max(d__1,d__2);
#line 1021 "zlarrv.f"
			    }
#line 1022 "zlarrv.f"
			    if (windex < wend) {
/* Computing MAX */
#line 1023 "zlarrv.f"
				d__1 = savgap, d__2 = w[windpl] - werr[windpl]
					 - w[windex] - werr[windex];
#line 1023 "zlarrv.f"
				wgap[windex] = max(d__1,d__2);
#line 1026 "zlarrv.f"
			    }
#line 1027 "zlarrv.f"
			}
#line 1028 "zlarrv.f"
			++idone;
#line 1029 "zlarrv.f"
		    }
/*                 here ends the code for the current child */

#line 1032 "zlarrv.f"
L139:
/*                 Proceed to any remaining child nodes */
#line 1034 "zlarrv.f"
		    newfst = j + 1;
#line 1035 "zlarrv.f"
L140:
#line 1035 "zlarrv.f"
		    ;
#line 1035 "zlarrv.f"
		}
#line 1036 "zlarrv.f"
/* L150: */
#line 1036 "zlarrv.f"
	    }
#line 1037 "zlarrv.f"
	    ++ndepth;
#line 1038 "zlarrv.f"
	    goto L40;
#line 1039 "zlarrv.f"
	}
#line 1040 "zlarrv.f"
	ibegin = iend + 1;
#line 1041 "zlarrv.f"
	wbegin = wend + 1;
#line 1042 "zlarrv.f"
L170:
#line 1042 "zlarrv.f"
	;
#line 1042 "zlarrv.f"
    }

#line 1045 "zlarrv.f"
    return 0;

/*     End of ZLARRV */

} /* zlarrv_ */

