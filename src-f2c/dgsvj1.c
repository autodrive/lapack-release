#line 1 "dgsvj1.f"
/* dgsvj1.f -- translated by f2c (version 20100827).
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

#line 1 "dgsvj1.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b35 = 1.;

/* > \brief \b DGSVJ1 pre-processor for the routine sgesvj, applies Jacobi rotations targeting only particular
 pivots. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGSVJ1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgsvj1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgsvj1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgsvj1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, */
/*                          EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   EPS, SFMIN, TOL */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP */
/*       CHARACTER*1        JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), D( N ), SVA( N ), V( LDV, * ), */
/*      $                   WORK( LWORK ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGSVJ1 is called from DGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as DGESVJ does, but */
/* > it targets only particular pivots and it does not check convergence */
/* > (stopping criterion). Few tunning parameters (marked by [TP]) are */
/* > available for the implementer. */
/* > */
/* > Further Details */
/* > ~~~~~~~~~~~~~~~ */
/* > DGSVJ1 applies few sweeps of Jacobi rotations in the column space of */
/* > the input M-by-N matrix A. The pivot pairs are taken from the (1,2) */
/* > off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The */
/* > block-entries (tiles) of the (1,2) off-diagonal block are marked by the */
/* > [x]'s in the following scheme: */
/* > */
/* >    | *  *  * [x] [x] [x]| */
/* >    | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks. */
/* >    | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block. */
/* >    |[x] [x] [x] *  *  * | */
/* >    |[x] [x] [x] *  *  * | */
/* >    |[x] [x] [x] *  *  * | */
/* > */
/* > In terms of the columns of A, the first N1 columns are rotated 'against' */
/* > the remaining N-N1 columns, trying to increase the angle between the */
/* > corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is */
/* > tiled using quadratic tiles of side KBL. Here, KBL is a tunning parmeter. */
/* > The number of sweeps is given in NSWEEP and the orthogonality threshold */
/* > is given in TOL. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >          Specifies whether the output from this procedure is used */
/* >          to compute the matrix V: */
/* >          = 'V': the product of the Jacobi rotations is accumulated */
/* >                 by postmulyiplying the N-by-N array V. */
/* >                (See the description of V.) */
/* >          = 'A': the product of the Jacobi rotations is accumulated */
/* >                 by postmulyiplying the MV-by-N array V. */
/* >                (See the descriptions of MV and V.) */
/* >          = 'N': the Jacobi rotations are not accumulated. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the input matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the input matrix A. */
/* >          M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* >          N1 specifies the 2 x 2 block partition, the first N1 columns are */
/* >          rotated 'against' the remaining N-N1 columns of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, M-by-N matrix A, such that A*diag(D) represents */
/* >          the input matrix. */
/* >          On exit, */
/* >          A_onexit * D_onexit represents the input matrix A*diag(D) */
/* >          post-multiplied by a sequence of Jacobi rotations, where the */
/* >          rotation threshold and the total number of sweeps are given in */
/* >          TOL and NSWEEP, respectively. */
/* >          (See the descriptions of N1, D, TOL and NSWEEP.) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The array D accumulates the scaling factors from the fast scaled */
/* >          Jacobi rotations. */
/* >          On entry, A*diag(D) represents the input matrix. */
/* >          On exit, A_onexit*diag(D_onexit) represents the input matrix */
/* >          post-multiplied by a sequence of Jacobi rotations, where the */
/* >          rotation threshold and the total number of sweeps are given in */
/* >          TOL and NSWEEP, respectively. */
/* >          (See the descriptions of N1, A, TOL and NSWEEP.) */
/* > \endverbatim */
/* > */
/* > \param[in,out] SVA */
/* > \verbatim */
/* >          SVA is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, SVA contains the Euclidean norms of the columns of */
/* >          the matrix A*diag(D). */
/* >          On exit, SVA contains the Euclidean norms of the columns of */
/* >          the matrix onexit*diag(D_onexit). */
/* > \endverbatim */
/* > */
/* > \param[in] MV */
/* > \verbatim */
/* >          MV is INTEGER */
/* >          If JOBV .EQ. 'A', then MV rows of V are post-multipled by a */
/* >                           sequence of Jacobi rotations. */
/* >          If JOBV = 'N',   then MV is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension (LDV,N) */
/* >          If JOBV .EQ. 'V' then N rows of V are post-multipled by a */
/* >                           sequence of Jacobi rotations. */
/* >          If JOBV .EQ. 'A' then MV rows of V are post-multipled by a */
/* >                           sequence of Jacobi rotations. */
/* >          If JOBV = 'N',   then V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V,  LDV >= 1. */
/* >          If JOBV = 'V', LDV .GE. N. */
/* >          If JOBV = 'A', LDV .GE. MV. */
/* > \endverbatim */
/* > */
/* > \param[in] EPS */
/* > \verbatim */
/* >          EPS is DOUBLE PRECISION */
/* >          EPS = DLAMCH('Epsilon') */
/* > \endverbatim */
/* > */
/* > \param[in] SFMIN */
/* > \verbatim */
/* >          SFMIN is DOUBLE PRECISION */
/* >          SFMIN = DLAMCH('Safe Minimum') */
/* > \endverbatim */
/* > */
/* > \param[in] TOL */
/* > \verbatim */
/* >          TOL is DOUBLE PRECISION */
/* >          TOL is the threshold for Jacobi rotations. For a pair */
/* >          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is */
/* >          applied only if DABS(COS(angle(A(:,p),A(:,q)))) .GT. TOL. */
/* > \endverbatim */
/* > */
/* > \param[in] NSWEEP */
/* > \verbatim */
/* >          NSWEEP is INTEGER */
/* >          NSWEEP is the number of sweeps of Jacobi rotations to be */
/* >          performed. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          LWORK is the dimension of WORK. LWORK .GE. M. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0 : successful exit. */
/* >          < 0 : if INFO = -i, then the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */

/*  ===================================================================== */
/* Subroutine */ int dgsvj1_(char *jobv, integer *m, integer *n, integer *n1, 
	doublereal *a, integer *lda, doublereal *d__, doublereal *sva, 
	integer *mv, doublereal *v, integer *ldv, doublereal *eps, doublereal 
	*sfmin, doublereal *tol, integer *nsweep, doublereal *work, integer *
	lwork, integer *info, ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal bigtheta;
    static integer pskipped, i__, p, q;
    static doublereal t, rootsfmin, cs, sn;
    static integer jbc;
    static doublereal big;
    static integer kbl, igl, ibr, jgl, mvl, nblc;
    static doublereal aapp, aapq, aaqq;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nblr, ierr;
    static doublereal aapp0;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp1, large, apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta, small;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal fastr[5];
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical applv, rsvec;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), drotm_(integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *);
    static logical rotok;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ijblsk, swband, blskip;
    static doublereal mxaapq;
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal thsign, mxsinj;
    static integer emptsw, notrot, iswrot;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 289 "dgsvj1.f"
    /* Parameter adjustments */
#line 289 "dgsvj1.f"
    --sva;
#line 289 "dgsvj1.f"
    --d__;
#line 289 "dgsvj1.f"
    a_dim1 = *lda;
#line 289 "dgsvj1.f"
    a_offset = 1 + a_dim1;
#line 289 "dgsvj1.f"
    a -= a_offset;
#line 289 "dgsvj1.f"
    v_dim1 = *ldv;
#line 289 "dgsvj1.f"
    v_offset = 1 + v_dim1;
#line 289 "dgsvj1.f"
    v -= v_offset;
#line 289 "dgsvj1.f"
    --work;
#line 289 "dgsvj1.f"

#line 289 "dgsvj1.f"
    /* Function Body */
#line 289 "dgsvj1.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 290 "dgsvj1.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 291 "dgsvj1.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 292 "dgsvj1.f"
	*info = -1;
#line 293 "dgsvj1.f"
    } else if (*m < 0) {
#line 294 "dgsvj1.f"
	*info = -2;
#line 295 "dgsvj1.f"
    } else if (*n < 0 || *n > *m) {
#line 296 "dgsvj1.f"
	*info = -3;
#line 297 "dgsvj1.f"
    } else if (*n1 < 0) {
#line 298 "dgsvj1.f"
	*info = -4;
#line 299 "dgsvj1.f"
    } else if (*lda < *m) {
#line 300 "dgsvj1.f"
	*info = -6;
#line 301 "dgsvj1.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 302 "dgsvj1.f"
	*info = -9;
#line 303 "dgsvj1.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 305 "dgsvj1.f"
	*info = -11;
#line 306 "dgsvj1.f"
    } else if (*tol <= *eps) {
#line 307 "dgsvj1.f"
	*info = -14;
#line 308 "dgsvj1.f"
    } else if (*nsweep < 0) {
#line 309 "dgsvj1.f"
	*info = -15;
#line 310 "dgsvj1.f"
    } else if (*lwork < *m) {
#line 311 "dgsvj1.f"
	*info = -17;
#line 312 "dgsvj1.f"
    } else {
#line 313 "dgsvj1.f"
	*info = 0;
#line 314 "dgsvj1.f"
    }

/*     #:( */
#line 317 "dgsvj1.f"
    if (*info != 0) {
#line 318 "dgsvj1.f"
	i__1 = -(*info);
#line 318 "dgsvj1.f"
	xerbla_("DGSVJ1", &i__1, (ftnlen)6);
#line 319 "dgsvj1.f"
	return 0;
#line 320 "dgsvj1.f"
    }

#line 322 "dgsvj1.f"
    if (rsvec) {
#line 323 "dgsvj1.f"
	mvl = *n;
#line 324 "dgsvj1.f"
    } else if (applv) {
#line 325 "dgsvj1.f"
	mvl = *mv;
#line 326 "dgsvj1.f"
    }
#line 327 "dgsvj1.f"
    rsvec = rsvec || applv;
#line 329 "dgsvj1.f"
    rooteps = sqrt(*eps);
#line 330 "dgsvj1.f"
    rootsfmin = sqrt(*sfmin);
#line 331 "dgsvj1.f"
    small = *sfmin / *eps;
#line 332 "dgsvj1.f"
    big = 1. / *sfmin;
#line 333 "dgsvj1.f"
    rootbig = 1. / rootsfmin;
#line 334 "dgsvj1.f"
    large = big / sqrt((doublereal) (*m * *n));
#line 335 "dgsvj1.f"
    bigtheta = 1. / rooteps;
#line 336 "dgsvj1.f"
    roottol = sqrt(*tol);

/*     .. Initialize the right singular vector matrix .. */

/*     RSVEC = LSAME( JOBV, 'Y' ) */

#line 342 "dgsvj1.f"
    emptsw = *n1 * (*n - *n1);
#line 343 "dgsvj1.f"
    notrot = 0;
#line 344 "dgsvj1.f"
    fastr[0] = 0.;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 348 "dgsvj1.f"
    kbl = min(8,*n);
#line 349 "dgsvj1.f"
    nblr = *n1 / kbl;
#line 350 "dgsvj1.f"
    if (nblr * kbl != *n1) {
#line 350 "dgsvj1.f"
	++nblr;
#line 350 "dgsvj1.f"
    }
/*     .. the tiling is nblr-by-nblc [tiles] */
#line 354 "dgsvj1.f"
    nblc = (*n - *n1) / kbl;
#line 355 "dgsvj1.f"
    if (nblc * kbl != *n - *n1) {
#line 355 "dgsvj1.f"
	++nblc;
#line 355 "dgsvj1.f"
    }
/* Computing 2nd power */
#line 356 "dgsvj1.f"
    i__1 = kbl;
#line 356 "dgsvj1.f"
    blskip = i__1 * i__1 + 1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
#line 359 "dgsvj1.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */
#line 361 "dgsvj1.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter. It is meaningful and effective */
/*     if SGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm SGESVJ. */


/*     | *   *   * [x] [x] [x]| */
/*     | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks. */
/*     | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block. */
/*     |[x] [x] [x] *   *   * | */
/*     |[x] [x] [x] *   *   * | */
/*     |[x] [x] [x] *   *   * | */


#line 375 "dgsvj1.f"
    i__1 = *nsweep;
#line 375 "dgsvj1.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     .. go go go ... */

#line 378 "dgsvj1.f"
	mxaapq = 0.;
#line 379 "dgsvj1.f"
	mxsinj = 0.;
#line 380 "dgsvj1.f"
	iswrot = 0;

#line 382 "dgsvj1.f"
	notrot = 0;
#line 383 "dgsvj1.f"
	pskipped = 0;

#line 385 "dgsvj1.f"
	i__2 = nblr;
#line 385 "dgsvj1.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {
#line 387 "dgsvj1.f"
	    igl = (ibr - 1) * kbl + 1;


/* ........................................................ */
/* ... go to the off diagonal blocks */
#line 393 "dgsvj1.f"
	    igl = (ibr - 1) * kbl + 1;
#line 395 "dgsvj1.f"
	    i__3 = nblc;
#line 395 "dgsvj1.f"
	    for (jbc = 1; jbc <= i__3; ++jbc) {
#line 397 "dgsvj1.f"
		jgl = *n1 + (jbc - 1) * kbl + 1;
/*        doing the block at ( ibr, jbc ) */
#line 401 "dgsvj1.f"
		ijblsk = 0;
/* Computing MIN */
#line 402 "dgsvj1.f"
		i__5 = igl + kbl - 1;
#line 402 "dgsvj1.f"
		i__4 = min(i__5,*n1);
#line 402 "dgsvj1.f"
		for (p = igl; p <= i__4; ++p) {
#line 404 "dgsvj1.f"
		    aapp = sva[p];
#line 406 "dgsvj1.f"
		    if (aapp > 0.) {
#line 408 "dgsvj1.f"
			pskipped = 0;
/* Computing MIN */
#line 410 "dgsvj1.f"
			i__6 = jgl + kbl - 1;
#line 410 "dgsvj1.f"
			i__5 = min(i__6,*n);
#line 410 "dgsvj1.f"
			for (q = jgl; q <= i__5; ++q) {

#line 412 "dgsvj1.f"
			    aaqq = sva[q];
#line 414 "dgsvj1.f"
			    if (aaqq > 0.) {
#line 415 "dgsvj1.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        .. Safe Gram matrix computation .. */

#line 421 "dgsvj1.f"
				if (aaqq >= 1.) {
#line 422 "dgsvj1.f"
				    if (aapp >= aaqq) {
#line 423 "dgsvj1.f"
					rotok = small * aapp <= aaqq;
#line 424 "dgsvj1.f"
				    } else {
#line 425 "dgsvj1.f"
					rotok = small * aaqq <= aapp;
#line 426 "dgsvj1.f"
				    }
#line 427 "dgsvj1.f"
				    if (aapp < big / aaqq) {
#line 428 "dgsvj1.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 431 "dgsvj1.f"
				    } else {
#line 432 "dgsvj1.f"
					dcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 433 "dgsvj1.f"
					dlascl_("G", &c__0, &c__0, &aapp, &
						d__[p], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 435 "dgsvj1.f"
					aapq = ddot_(m, &work[1], &c__1, &a[q 
						* a_dim1 + 1], &c__1) * d__[q]
						 / aaqq;
#line 437 "dgsvj1.f"
				    }
#line 438 "dgsvj1.f"
				} else {
#line 439 "dgsvj1.f"
				    if (aapp >= aaqq) {
#line 440 "dgsvj1.f"
					rotok = aapp <= aaqq / small;
#line 441 "dgsvj1.f"
				    } else {
#line 442 "dgsvj1.f"
					rotok = aaqq <= aapp / small;
#line 443 "dgsvj1.f"
				    }
#line 444 "dgsvj1.f"
				    if (aapp > small / aaqq) {
#line 445 "dgsvj1.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 448 "dgsvj1.f"
				    } else {
#line 449 "dgsvj1.f"
					dcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 450 "dgsvj1.f"
					dlascl_("G", &c__0, &c__0, &aaqq, &
						d__[q], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 452 "dgsvj1.f"
					aapq = ddot_(m, &work[1], &c__1, &a[p 
						* a_dim1 + 1], &c__1) * d__[p]
						 / aapp;
#line 454 "dgsvj1.f"
				    }
#line 455 "dgsvj1.f"
				}
/* Computing MAX */
#line 457 "dgsvj1.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 457 "dgsvj1.f"
				mxaapq = max(d__1,d__2);
/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 461 "dgsvj1.f"
				if (abs(aapq) > *tol) {
#line 462 "dgsvj1.f"
				    notrot = 0;
/*           ROTATED  = ROTATED + 1 */
#line 464 "dgsvj1.f"
				    pskipped = 0;
#line 465 "dgsvj1.f"
				    ++iswrot;

#line 467 "dgsvj1.f"
				    if (rotok) {

#line 469 "dgsvj1.f"
					aqoap = aaqq / aapp;
#line 470 "dgsvj1.f"
					apoaq = aapp / aaqq;
#line 471 "dgsvj1.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;
#line 472 "dgsvj1.f"
					if (aaqq > aapp0) {
#line 472 "dgsvj1.f"
					    theta = -theta;
#line 472 "dgsvj1.f"
					}
#line 474 "dgsvj1.f"
					if (abs(theta) > bigtheta) {
#line 475 "dgsvj1.f"
					    t = .5 / theta;
#line 476 "dgsvj1.f"
					    fastr[2] = t * d__[p] / d__[q];
#line 477 "dgsvj1.f"
					    fastr[3] = -t * d__[q] / d__[p];
#line 478 "dgsvj1.f"
					    drotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 480 "dgsvj1.f"
					    if (rsvec) {
#line 480 "dgsvj1.f"
			  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 480 "dgsvj1.f"
					    }
/* Computing MAX */
#line 484 "dgsvj1.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 484 "dgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 486 "dgsvj1.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 486 "dgsvj1.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 488 "dgsvj1.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 488 "dgsvj1.f"
					    mxsinj = max(d__1,d__2);
#line 489 "dgsvj1.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 493 "dgsvj1.f"
					    thsign = -d_sign(&c_b35, &aapq);
#line 494 "dgsvj1.f"
					    if (aaqq > aapp0) {
#line 494 "dgsvj1.f"
			  thsign = -thsign;
#line 494 "dgsvj1.f"
					    }
#line 495 "dgsvj1.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 497 "dgsvj1.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 498 "dgsvj1.f"
					    sn = t * cs;
/* Computing MAX */
#line 499 "dgsvj1.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 499 "dgsvj1.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 500 "dgsvj1.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 500 "dgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 502 "dgsvj1.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 502 "dgsvj1.f"
					    aapp *= sqrt((max(d__1,d__2)));
#line 505 "dgsvj1.f"
					    apoaq = d__[p] / d__[q];
#line 506 "dgsvj1.f"
					    aqoap = d__[q] / d__[p];
#line 507 "dgsvj1.f"
					    if (d__[p] >= 1.) {

#line 509 "dgsvj1.f"
			  if (d__[q] >= 1.) {
#line 510 "dgsvj1.f"
			      fastr[2] = t * apoaq;
#line 511 "dgsvj1.f"
			      fastr[3] = -t * aqoap;
#line 512 "dgsvj1.f"
			      d__[p] *= cs;
#line 513 "dgsvj1.f"
			      d__[q] *= cs;
#line 514 "dgsvj1.f"
			      drotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 517 "dgsvj1.f"
			      if (rsvec) {
#line 517 "dgsvj1.f"
				  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 517 "dgsvj1.f"
			      }
#line 520 "dgsvj1.f"
			  } else {
#line 521 "dgsvj1.f"
			      d__1 = -t * aqoap;
#line 521 "dgsvj1.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 524 "dgsvj1.f"
			      d__1 = cs * sn * apoaq;
#line 524 "dgsvj1.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 527 "dgsvj1.f"
			      if (rsvec) {
#line 528 "dgsvj1.f"
				  d__1 = -t * aqoap;
#line 528 "dgsvj1.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 531 "dgsvj1.f"
				  d__1 = cs * sn * apoaq;
#line 531 "dgsvj1.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 535 "dgsvj1.f"
			      }
#line 536 "dgsvj1.f"
			      d__[p] *= cs;
#line 537 "dgsvj1.f"
			      d__[q] /= cs;
#line 538 "dgsvj1.f"
			  }
#line 539 "dgsvj1.f"
					    } else {
#line 540 "dgsvj1.f"
			  if (d__[q] >= 1.) {
#line 541 "dgsvj1.f"
			      d__1 = t * apoaq;
#line 541 "dgsvj1.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 544 "dgsvj1.f"
			      d__1 = -cs * sn * aqoap;
#line 544 "dgsvj1.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 547 "dgsvj1.f"
			      if (rsvec) {
#line 548 "dgsvj1.f"
				  d__1 = t * apoaq;
#line 548 "dgsvj1.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 551 "dgsvj1.f"
				  d__1 = -cs * sn * aqoap;
#line 551 "dgsvj1.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 555 "dgsvj1.f"
			      }
#line 556 "dgsvj1.f"
			      d__[p] /= cs;
#line 557 "dgsvj1.f"
			      d__[q] *= cs;
#line 558 "dgsvj1.f"
			  } else {
#line 559 "dgsvj1.f"
			      if (d__[p] >= d__[q]) {
#line 560 "dgsvj1.f"
				  d__1 = -t * aqoap;
#line 560 "dgsvj1.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 563 "dgsvj1.f"
				  d__1 = cs * sn * apoaq;
#line 563 "dgsvj1.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 566 "dgsvj1.f"
				  d__[p] *= cs;
#line 567 "dgsvj1.f"
				  d__[q] /= cs;
#line 568 "dgsvj1.f"
				  if (rsvec) {
#line 569 "dgsvj1.f"
				      d__1 = -t * aqoap;
#line 569 "dgsvj1.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 573 "dgsvj1.f"
				      d__1 = cs * sn * apoaq;
#line 573 "dgsvj1.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 577 "dgsvj1.f"
				  }
#line 578 "dgsvj1.f"
			      } else {
#line 579 "dgsvj1.f"
				  d__1 = t * apoaq;
#line 579 "dgsvj1.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 582 "dgsvj1.f"
				  d__1 = -cs * sn * aqoap;
#line 582 "dgsvj1.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 586 "dgsvj1.f"
				  d__[p] /= cs;
#line 587 "dgsvj1.f"
				  d__[q] *= cs;
#line 588 "dgsvj1.f"
				  if (rsvec) {
#line 589 "dgsvj1.f"
				      d__1 = t * apoaq;
#line 589 "dgsvj1.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 592 "dgsvj1.f"
				      d__1 = -cs * sn * aqoap;
#line 592 "dgsvj1.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 596 "dgsvj1.f"
				  }
#line 597 "dgsvj1.f"
			      }
#line 598 "dgsvj1.f"
			  }
#line 599 "dgsvj1.f"
					    }
#line 600 "dgsvj1.f"
					}
#line 602 "dgsvj1.f"
				    } else {
#line 603 "dgsvj1.f"
					if (aapp > aaqq) {
#line 604 "dgsvj1.f"
					    dcopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 606 "dgsvj1.f"
					    dlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b35, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 608 "dgsvj1.f"
					    dlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b35, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 611 "dgsvj1.f"
					    temp1 = -aapq * d__[p] / d__[q];
#line 612 "dgsvj1.f"
					    daxpy_(m, &temp1, &work[1], &c__1,
						     &a[q * a_dim1 + 1], &
						    c__1);
#line 614 "dgsvj1.f"
					    dlascl_("G", &c__0, &c__0, &c_b35,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 617 "dgsvj1.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 617 "dgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 619 "dgsvj1.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 620 "dgsvj1.f"
					} else {
#line 621 "dgsvj1.f"
					    dcopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 623 "dgsvj1.f"
					    dlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b35, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 625 "dgsvj1.f"
					    dlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b35, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 628 "dgsvj1.f"
					    temp1 = -aapq * d__[q] / d__[p];
#line 629 "dgsvj1.f"
					    daxpy_(m, &temp1, &work[1], &c__1,
						     &a[p * a_dim1 + 1], &
						    c__1);
#line 631 "dgsvj1.f"
					    dlascl_("G", &c__0, &c__0, &c_b35,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 634 "dgsvj1.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 634 "dgsvj1.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 636 "dgsvj1.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 637 "dgsvj1.f"
					}
#line 638 "dgsvj1.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q) */
/*           .. recompute SVA(q) */
/* Computing 2nd power */
#line 643 "dgsvj1.f"
				    d__1 = sva[q] / aaqq;
#line 643 "dgsvj1.f"
				    if (d__1 * d__1 <= rooteps) {
#line 645 "dgsvj1.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 647 "dgsvj1.f"
					    sva[q] = dnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * d__[q];
#line 649 "dgsvj1.f"
					} else {
#line 650 "dgsvj1.f"
					    t = 0.;
#line 651 "dgsvj1.f"
					    aaqq = 1.;
#line 652 "dgsvj1.f"
					    dlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 654 "dgsvj1.f"
					    sva[q] = t * sqrt(aaqq) * d__[q];
#line 655 "dgsvj1.f"
					}
#line 656 "dgsvj1.f"
				    }
/* Computing 2nd power */
#line 657 "dgsvj1.f"
				    d__1 = aapp / aapp0;
#line 657 "dgsvj1.f"
				    if (d__1 * d__1 <= rooteps) {
#line 658 "dgsvj1.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 660 "dgsvj1.f"
					    aapp = dnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * d__[p];
#line 662 "dgsvj1.f"
					} else {
#line 663 "dgsvj1.f"
					    t = 0.;
#line 664 "dgsvj1.f"
					    aapp = 1.;
#line 665 "dgsvj1.f"
					    dlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 667 "dgsvj1.f"
					    aapp = t * sqrt(aapp) * d__[p];
#line 668 "dgsvj1.f"
					}
#line 669 "dgsvj1.f"
					sva[p] = aapp;
#line 670 "dgsvj1.f"
				    }
/*              end of OK rotation */
#line 672 "dgsvj1.f"
				} else {
#line 673 "dgsvj1.f"
				    ++notrot;
/*           SKIPPED  = SKIPPED  + 1 */
#line 675 "dgsvj1.f"
				    ++pskipped;
#line 676 "dgsvj1.f"
				    ++ijblsk;
#line 677 "dgsvj1.f"
				}
#line 678 "dgsvj1.f"
			    } else {
#line 679 "dgsvj1.f"
				++notrot;
#line 680 "dgsvj1.f"
				++pskipped;
#line 681 "dgsvj1.f"
				++ijblsk;
#line 682 "dgsvj1.f"
			    }
/*      IF ( NOTROT .GE. EMPTSW )  GO TO 2011 */
#line 685 "dgsvj1.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 687 "dgsvj1.f"
				sva[p] = aapp;
#line 688 "dgsvj1.f"
				notrot = 0;
#line 689 "dgsvj1.f"
				goto L2011;
#line 690 "dgsvj1.f"
			    }
#line 691 "dgsvj1.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 693 "dgsvj1.f"
				aapp = -aapp;
#line 694 "dgsvj1.f"
				notrot = 0;
#line 695 "dgsvj1.f"
				goto L2203;
#line 696 "dgsvj1.f"
			    }

#line 699 "dgsvj1.f"
/* L2200: */
#line 699 "dgsvj1.f"
			}
/*        end of the q-loop */
#line 701 "dgsvj1.f"
L2203:
#line 703 "dgsvj1.f"
			sva[p] = aapp;

#line 705 "dgsvj1.f"
		    } else {
#line 706 "dgsvj1.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 706 "dgsvj1.f"
			    i__5 = jgl + kbl - 1;
#line 706 "dgsvj1.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 706 "dgsvj1.f"
			}
#line 708 "dgsvj1.f"
			if (aapp < 0.) {
#line 708 "dgsvj1.f"
			    notrot = 0;
#line 708 "dgsvj1.f"
			}
/* **      IF ( NOTROT .GE. EMPTSW )  GO TO 2011 */
#line 710 "dgsvj1.f"
		    }
#line 712 "dgsvj1.f"
/* L2100: */
#line 712 "dgsvj1.f"
		}
/*     end of the p-loop */
#line 714 "dgsvj1.f"
/* L2010: */
#line 714 "dgsvj1.f"
	    }
/*     end of the jbc-loop */
#line 716 "dgsvj1.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 718 "dgsvj1.f"
	    i__4 = igl + kbl - 1;
#line 718 "dgsvj1.f"
	    i__3 = min(i__4,*n);
#line 718 "dgsvj1.f"
	    for (p = igl; p <= i__3; ++p) {
#line 719 "dgsvj1.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 720 "dgsvj1.f"
/* L2012: */
#line 720 "dgsvj1.f"
	    }
/* **   IF ( NOTROT .GE. EMPTSW ) GO TO 1994 */
#line 722 "dgsvj1.f"
/* L2000: */
#line 722 "dgsvj1.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 726 "dgsvj1.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 728 "dgsvj1.f"
	    sva[*n] = dnrm2_(m, &a[*n * a_dim1 + 1], &c__1) * d__[*n];
#line 729 "dgsvj1.f"
	} else {
#line 730 "dgsvj1.f"
	    t = 0.;
#line 731 "dgsvj1.f"
	    aapp = 1.;
#line 732 "dgsvj1.f"
	    dlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 733 "dgsvj1.f"
	    sva[*n] = t * sqrt(aapp) * d__[*n];
#line 734 "dgsvj1.f"
	}

/*     Additional steering devices */

#line 738 "dgsvj1.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 738 "dgsvj1.f"
	    swband = i__;
#line 738 "dgsvj1.f"
	}
#line 741 "dgsvj1.f"
	if (i__ > swband + 1 && mxaapq < (doublereal) (*n) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 743 "dgsvj1.f"
	    goto L1994;
#line 744 "dgsvj1.f"
	}

#line 747 "dgsvj1.f"
	if (notrot >= emptsw) {
#line 747 "dgsvj1.f"
	    goto L1994;
#line 747 "dgsvj1.f"
	}
#line 749 "dgsvj1.f"
/* L1993: */
#line 749 "dgsvj1.f"
    }
/*     end i=1:NSWEEP loop */
/* #:) Reaching this point means that the procedure has completed the given */
/*     number of sweeps. */
#line 753 "dgsvj1.f"
    *info = *nsweep - 1;
#line 754 "dgsvj1.f"
    goto L1995;
#line 755 "dgsvj1.f"
L1994:
/* #:) Reaching this point means that during the i-th sweep all pivots were */
/*     below the given threshold, causing early exit. */
#line 759 "dgsvj1.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 761 "dgsvj1.f"
L1995:

/*     Sort the vector D */

#line 765 "dgsvj1.f"
    i__1 = *n - 1;
#line 765 "dgsvj1.f"
    for (p = 1; p <= i__1; ++p) {
#line 766 "dgsvj1.f"
	i__2 = *n - p + 1;
#line 766 "dgsvj1.f"
	q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 767 "dgsvj1.f"
	if (p != q) {
#line 768 "dgsvj1.f"
	    temp1 = sva[p];
#line 769 "dgsvj1.f"
	    sva[p] = sva[q];
#line 770 "dgsvj1.f"
	    sva[q] = temp1;
#line 771 "dgsvj1.f"
	    temp1 = d__[p];
#line 772 "dgsvj1.f"
	    d__[p] = d__[q];
#line 773 "dgsvj1.f"
	    d__[q] = temp1;
#line 774 "dgsvj1.f"
	    dswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 775 "dgsvj1.f"
	    if (rsvec) {
#line 775 "dgsvj1.f"
		dswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 775 "dgsvj1.f"
	    }
#line 776 "dgsvj1.f"
	}
#line 777 "dgsvj1.f"
/* L5991: */
#line 777 "dgsvj1.f"
    }

#line 779 "dgsvj1.f"
    return 0;
/*     .. */
/*     .. END OF DGSVJ1 */
/*     .. */
} /* dgsvj1_ */

