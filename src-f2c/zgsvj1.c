#line 1 "zgsvj1.f"
/* zgsvj1.f -- translated by f2c (version 20100827).
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

#line 1 "zgsvj1.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b18 = 1.;

/* > \brief \b ZGSVJ1 pre-processor for the routine zgesvj, applies Jacobi rotations targeting only particular
 pivots. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGSVJ1 + dependencies */
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

/*       SUBROUTINE ZGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, */
/*                          EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   EPS, SFMIN, TOL */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP */
/*       CHARACTER*1        JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK ) */
/*       DOUBLE PRECISION   SVA( N ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGSVJ1 is called from ZGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as ZGESVJ does, but */
/* > it targets only particular pivots and it does not check convergence */
/* > (stopping criterion). Few tunning parameters (marked by [TP]) are */
/* > available for the implementer. */
/* > */
/* > Further Details */
/* > ~~~~~~~~~~~~~~~ */
/* > ZGSVJ1 applies few sweeps of Jacobi rotations in the column space of */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          D is COMPLEX*16 array, dimension (N) */
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
/* >          V is COMPLEX*16 array, dimension (LDV,N) */
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
/* >          applied only if ABS(COS(angle(A(:,p),A(:,q)))) .GT. TOL. */
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
/* >          WORK is COMPLEX*16 array, dimension (LWORK) */
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

/* > \date June 2016 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Contributor: */
/*  ================== */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) */

/*  ===================================================================== */
/* Subroutine */ int zgsvj1_(char *jobv, integer *m, integer *n, integer *n1, 
	doublecomplex *a, integer *lda, doublecomplex *d__, doublereal *sva, 
	integer *mv, doublecomplex *v, integer *ldv, doublereal *eps, 
	doublereal *sfmin, doublereal *tol, integer *nsweep, doublecomplex *
	work, integer *lwork, integer *info, ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal bigtheta;
    static integer pskipped, i__, p, q;
    static doublereal t, rootsfmin, cs, sn;
    static integer jbc;
    static doublereal big;
    static integer kbl, igl, ibr, jgl, mvl, nblc;
    static doublereal aapp;
    static doublecomplex aapq;
    static doublereal aaqq;
    static integer nblr, ierr;
    static doublecomplex ompq;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublereal aapp0, aapq1, temp1, apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta, small;
    static logical applv, rsvec;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical rotok;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zaxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ijblsk, swband, blskip;
    static doublereal mxaapq;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static doublereal thsign, mxsinj;
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    static integer emptsw, notrot, iswrot;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. from BLAS */
/*     .. from LAPACK */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 293 "zgsvj1.f"
    /* Parameter adjustments */
#line 293 "zgsvj1.f"
    --sva;
#line 293 "zgsvj1.f"
    --d__;
#line 293 "zgsvj1.f"
    a_dim1 = *lda;
#line 293 "zgsvj1.f"
    a_offset = 1 + a_dim1;
#line 293 "zgsvj1.f"
    a -= a_offset;
#line 293 "zgsvj1.f"
    v_dim1 = *ldv;
#line 293 "zgsvj1.f"
    v_offset = 1 + v_dim1;
#line 293 "zgsvj1.f"
    v -= v_offset;
#line 293 "zgsvj1.f"
    --work;
#line 293 "zgsvj1.f"

#line 293 "zgsvj1.f"
    /* Function Body */
#line 293 "zgsvj1.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 294 "zgsvj1.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 295 "zgsvj1.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 296 "zgsvj1.f"
	*info = -1;
#line 297 "zgsvj1.f"
    } else if (*m < 0) {
#line 298 "zgsvj1.f"
	*info = -2;
#line 299 "zgsvj1.f"
    } else if (*n < 0 || *n > *m) {
#line 300 "zgsvj1.f"
	*info = -3;
#line 301 "zgsvj1.f"
    } else if (*n1 < 0) {
#line 302 "zgsvj1.f"
	*info = -4;
#line 303 "zgsvj1.f"
    } else if (*lda < *m) {
#line 304 "zgsvj1.f"
	*info = -6;
#line 305 "zgsvj1.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 306 "zgsvj1.f"
	*info = -9;
#line 307 "zgsvj1.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 309 "zgsvj1.f"
	*info = -11;
#line 310 "zgsvj1.f"
    } else if (*tol <= *eps) {
#line 311 "zgsvj1.f"
	*info = -14;
#line 312 "zgsvj1.f"
    } else if (*nsweep < 0) {
#line 313 "zgsvj1.f"
	*info = -15;
#line 314 "zgsvj1.f"
    } else if (*lwork < *m) {
#line 315 "zgsvj1.f"
	*info = -17;
#line 316 "zgsvj1.f"
    } else {
#line 317 "zgsvj1.f"
	*info = 0;
#line 318 "zgsvj1.f"
    }

/*     #:( */
#line 321 "zgsvj1.f"
    if (*info != 0) {
#line 322 "zgsvj1.f"
	i__1 = -(*info);
#line 322 "zgsvj1.f"
	xerbla_("ZGSVJ1", &i__1, (ftnlen)6);
#line 323 "zgsvj1.f"
	return 0;
#line 324 "zgsvj1.f"
    }

#line 326 "zgsvj1.f"
    if (rsvec) {
#line 327 "zgsvj1.f"
	mvl = *n;
#line 328 "zgsvj1.f"
    } else if (applv) {
#line 329 "zgsvj1.f"
	mvl = *mv;
#line 330 "zgsvj1.f"
    }
#line 331 "zgsvj1.f"
    rsvec = rsvec || applv;
#line 333 "zgsvj1.f"
    rooteps = sqrt(*eps);
#line 334 "zgsvj1.f"
    rootsfmin = sqrt(*sfmin);
#line 335 "zgsvj1.f"
    small = *sfmin / *eps;
#line 336 "zgsvj1.f"
    big = 1. / *sfmin;
#line 337 "zgsvj1.f"
    rootbig = 1. / rootsfmin;
/*     LARGE = BIG / SQRT( DBLE( M*N ) ) */
#line 339 "zgsvj1.f"
    bigtheta = 1. / rooteps;
#line 340 "zgsvj1.f"
    roottol = sqrt(*tol);

/*     .. Initialize the right singular vector matrix .. */

/*     RSVEC = LSAME( JOBV, 'Y' ) */

#line 346 "zgsvj1.f"
    emptsw = *n1 * (*n - *n1);
#line 347 "zgsvj1.f"
    notrot = 0;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 351 "zgsvj1.f"
    kbl = min(8,*n);
#line 352 "zgsvj1.f"
    nblr = *n1 / kbl;
#line 353 "zgsvj1.f"
    if (nblr * kbl != *n1) {
#line 353 "zgsvj1.f"
	++nblr;
#line 353 "zgsvj1.f"
    }
/*     .. the tiling is nblr-by-nblc [tiles] */
#line 357 "zgsvj1.f"
    nblc = (*n - *n1) / kbl;
#line 358 "zgsvj1.f"
    if (nblc * kbl != *n - *n1) {
#line 358 "zgsvj1.f"
	++nblc;
#line 358 "zgsvj1.f"
    }
/* Computing 2nd power */
#line 359 "zgsvj1.f"
    i__1 = kbl;
#line 359 "zgsvj1.f"
    blskip = i__1 * i__1 + 1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
#line 362 "zgsvj1.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */
#line 364 "zgsvj1.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter. It is meaningful and effective */
/*     if ZGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm ZGEJSV. */


/*     | *   *   * [x] [x] [x]| */
/*     | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks. */
/*     | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block. */
/*     |[x] [x] [x] *   *   * | */
/*     |[x] [x] [x] *   *   * | */
/*     |[x] [x] [x] *   *   * | */


#line 378 "zgsvj1.f"
    i__1 = *nsweep;
#line 378 "zgsvj1.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     .. go go go ... */

#line 382 "zgsvj1.f"
	mxaapq = 0.;
#line 383 "zgsvj1.f"
	mxsinj = 0.;
#line 384 "zgsvj1.f"
	iswrot = 0;

#line 386 "zgsvj1.f"
	notrot = 0;
#line 387 "zgsvj1.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 394 "zgsvj1.f"
	i__2 = nblr;
#line 394 "zgsvj1.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {

#line 396 "zgsvj1.f"
	    igl = (ibr - 1) * kbl + 1;


/* ... go to the off diagonal blocks */

#line 402 "zgsvj1.f"
	    igl = (ibr - 1) * kbl + 1;

/*            DO 2010 jbc = ibr + 1, NBL */
#line 405 "zgsvj1.f"
	    i__3 = nblc;
#line 405 "zgsvj1.f"
	    for (jbc = 1; jbc <= i__3; ++jbc) {

#line 407 "zgsvj1.f"
		jgl = (jbc - 1) * kbl + *n1 + 1;

/*        doing the block at ( ibr, jbc ) */

#line 411 "zgsvj1.f"
		ijblsk = 0;
/* Computing MIN */
#line 412 "zgsvj1.f"
		i__5 = igl + kbl - 1;
#line 412 "zgsvj1.f"
		i__4 = min(i__5,*n1);
#line 412 "zgsvj1.f"
		for (p = igl; p <= i__4; ++p) {

#line 414 "zgsvj1.f"
		    aapp = sva[p];
#line 415 "zgsvj1.f"
		    if (aapp > 0.) {

#line 417 "zgsvj1.f"
			pskipped = 0;

/* Computing MIN */
#line 419 "zgsvj1.f"
			i__6 = jgl + kbl - 1;
#line 419 "zgsvj1.f"
			i__5 = min(i__6,*n);
#line 419 "zgsvj1.f"
			for (q = jgl; q <= i__5; ++q) {

#line 421 "zgsvj1.f"
			    aaqq = sva[q];
#line 422 "zgsvj1.f"
			    if (aaqq > 0.) {
#line 423 "zgsvj1.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 429 "zgsvj1.f"
				if (aaqq >= 1.) {
#line 430 "zgsvj1.f"
				    if (aapp >= aaqq) {
#line 431 "zgsvj1.f"
					rotok = small * aapp <= aaqq;
#line 432 "zgsvj1.f"
				    } else {
#line 433 "zgsvj1.f"
					rotok = small * aaqq <= aapp;
#line 434 "zgsvj1.f"
				    }
#line 435 "zgsvj1.f"
				    if (aapp < big / aaqq) {
#line 436 "zgsvj1.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 436 "zgsvj1.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 436 "zgsvj1.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 436 "zgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 438 "zgsvj1.f"
				    } else {
#line 439 "zgsvj1.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 441 "zgsvj1.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b18, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 444 "zgsvj1.f"
					zdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 444 "zgsvj1.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 444 "zgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 446 "zgsvj1.f"
				    }
#line 447 "zgsvj1.f"
				} else {
#line 448 "zgsvj1.f"
				    if (aapp >= aaqq) {
#line 449 "zgsvj1.f"
					rotok = aapp <= aaqq / small;
#line 450 "zgsvj1.f"
				    } else {
#line 451 "zgsvj1.f"
					rotok = aaqq <= aapp / small;
#line 452 "zgsvj1.f"
				    }
#line 453 "zgsvj1.f"
				    if (aapp > small / aaqq) {
#line 454 "zgsvj1.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 454 "zgsvj1.f"
					d__1 = max(aaqq,aapp);
#line 454 "zgsvj1.f"
					z__2.r = z__3.r / d__1, z__2.i = 
						z__3.i / d__1;
#line 454 "zgsvj1.f"
					d__2 = min(aaqq,aapp);
#line 454 "zgsvj1.f"
					z__1.r = z__2.r / d__2, z__1.i = 
						z__2.i / d__2;
#line 454 "zgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 457 "zgsvj1.f"
				    } else {
#line 458 "zgsvj1.f"
					zcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 460 "zgsvj1.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b18, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 463 "zgsvj1.f"
					zdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 463 "zgsvj1.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 463 "zgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 465 "zgsvj1.f"
				    }
#line 466 "zgsvj1.f"
				}

/*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
#line 469 "zgsvj1.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 470 "zgsvj1.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 470 "zgsvj1.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 474 "zgsvj1.f"
				if (abs(aapq1) > *tol) {
#line 475 "zgsvj1.f"
				    d__1 = z_abs(&aapq);
#line 475 "zgsvj1.f"
				    z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					    d__1;
#line 475 "zgsvj1.f"
				    ompq.r = z__1.r, ompq.i = z__1.i;
#line 476 "zgsvj1.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 478 "zgsvj1.f"
				    pskipped = 0;
#line 479 "zgsvj1.f"
				    ++iswrot;

#line 481 "zgsvj1.f"
				    if (rotok) {

#line 483 "zgsvj1.f"
					aqoap = aaqq / aapp;
#line 484 "zgsvj1.f"
					apoaq = aapp / aaqq;
#line 485 "zgsvj1.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 486 "zgsvj1.f"
					if (aaqq > aapp0) {
#line 486 "zgsvj1.f"
					    theta = -theta;
#line 486 "zgsvj1.f"
					}

#line 488 "zgsvj1.f"
					if (abs(theta) > bigtheta) {
#line 489 "zgsvj1.f"
					    t = .5 / theta;
#line 490 "zgsvj1.f"
					    cs = 1.;
#line 491 "zgsvj1.f"
					    d_cnjg(&z__2, &ompq);
#line 491 "zgsvj1.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 491 "zgsvj1.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 493 "zgsvj1.f"
					    if (rsvec) {
#line 494 "zgsvj1.f"
			  d_cnjg(&z__2, &ompq);
#line 494 "zgsvj1.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 494 "zgsvj1.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 496 "zgsvj1.f"
					    }
/* Computing MAX */
#line 497 "zgsvj1.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 497 "zgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 499 "zgsvj1.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 499 "zgsvj1.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 501 "zgsvj1.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 501 "zgsvj1.f"
					    mxsinj = max(d__1,d__2);
#line 502 "zgsvj1.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 506 "zgsvj1.f"
					    thsign = -d_sign(&c_b18, &aapq1);
#line 507 "zgsvj1.f"
					    if (aaqq > aapp0) {
#line 507 "zgsvj1.f"
			  thsign = -thsign;
#line 507 "zgsvj1.f"
					    }
#line 508 "zgsvj1.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 510 "zgsvj1.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 511 "zgsvj1.f"
					    sn = t * cs;
/* Computing MAX */
#line 512 "zgsvj1.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 512 "zgsvj1.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 513 "zgsvj1.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 513 "zgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 515 "zgsvj1.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 515 "zgsvj1.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 518 "zgsvj1.f"
					    d_cnjg(&z__2, &ompq);
#line 518 "zgsvj1.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 518 "zgsvj1.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 520 "zgsvj1.f"
					    if (rsvec) {
#line 521 "zgsvj1.f"
			  d_cnjg(&z__2, &ompq);
#line 521 "zgsvj1.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 521 "zgsvj1.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 523 "zgsvj1.f"
					    }
#line 524 "zgsvj1.f"
					}
#line 525 "zgsvj1.f"
					i__6 = p;
#line 525 "zgsvj1.f"
					i__7 = q;
#line 525 "zgsvj1.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 525 "zgsvj1.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 525 "zgsvj1.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 527 "zgsvj1.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 529 "zgsvj1.f"
					if (aapp > aaqq) {
#line 530 "zgsvj1.f"
					    zcopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 532 "zgsvj1.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 535 "zgsvj1.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 538 "zgsvj1.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 538 "zgsvj1.f"
					    zaxpy_(m, &z__1, &work[1], &c__1, 
						    &a[q * a_dim1 + 1], &c__1)
						    ;
#line 540 "zgsvj1.f"
					    zlascl_("G", &c__0, &c__0, &c_b18,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 543 "zgsvj1.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 543 "zgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 545 "zgsvj1.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 546 "zgsvj1.f"
					} else {
#line 547 "zgsvj1.f"
					    zcopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 549 "zgsvj1.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 552 "zgsvj1.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 555 "zgsvj1.f"
					    d_cnjg(&z__2, &aapq);
#line 555 "zgsvj1.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 555 "zgsvj1.f"
					    zaxpy_(m, &z__1, &work[1], &c__1, 
						    &a[p * a_dim1 + 1], &c__1)
						    ;
#line 557 "zgsvj1.f"
					    zlascl_("G", &c__0, &c__0, &c_b18,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 560 "zgsvj1.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 560 "zgsvj1.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 562 "zgsvj1.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 563 "zgsvj1.f"
					}
#line 564 "zgsvj1.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 569 "zgsvj1.f"
				    d__1 = sva[q] / aaqq;
#line 569 "zgsvj1.f"
				    if (d__1 * d__1 <= rooteps) {
#line 571 "zgsvj1.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 573 "zgsvj1.f"
					    sva[q] = dznrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 574 "zgsvj1.f"
					} else {
#line 575 "zgsvj1.f"
					    t = 0.;
#line 576 "zgsvj1.f"
					    aaqq = 1.;
#line 577 "zgsvj1.f"
					    zlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 579 "zgsvj1.f"
					    sva[q] = t * sqrt(aaqq);
#line 580 "zgsvj1.f"
					}
#line 581 "zgsvj1.f"
				    }
/* Computing 2nd power */
#line 582 "zgsvj1.f"
				    d__1 = aapp / aapp0;
#line 582 "zgsvj1.f"
				    if (d__1 * d__1 <= rooteps) {
#line 583 "zgsvj1.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 585 "zgsvj1.f"
					    aapp = dznrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 586 "zgsvj1.f"
					} else {
#line 587 "zgsvj1.f"
					    t = 0.;
#line 588 "zgsvj1.f"
					    aapp = 1.;
#line 589 "zgsvj1.f"
					    zlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 591 "zgsvj1.f"
					    aapp = t * sqrt(aapp);
#line 592 "zgsvj1.f"
					}
#line 593 "zgsvj1.f"
					sva[p] = aapp;
#line 594 "zgsvj1.f"
				    }
/*              end of OK rotation */
#line 596 "zgsvj1.f"
				} else {
#line 597 "zgsvj1.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 599 "zgsvj1.f"
				    ++pskipped;
#line 600 "zgsvj1.f"
				    ++ijblsk;
#line 601 "zgsvj1.f"
				}
#line 602 "zgsvj1.f"
			    } else {
#line 603 "zgsvj1.f"
				++notrot;
#line 604 "zgsvj1.f"
				++pskipped;
#line 605 "zgsvj1.f"
				++ijblsk;
#line 606 "zgsvj1.f"
			    }

#line 608 "zgsvj1.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 610 "zgsvj1.f"
				sva[p] = aapp;
#line 611 "zgsvj1.f"
				notrot = 0;
#line 612 "zgsvj1.f"
				goto L2011;
#line 613 "zgsvj1.f"
			    }
#line 614 "zgsvj1.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 616 "zgsvj1.f"
				aapp = -aapp;
#line 617 "zgsvj1.f"
				notrot = 0;
#line 618 "zgsvj1.f"
				goto L2203;
#line 619 "zgsvj1.f"
			    }

#line 621 "zgsvj1.f"
/* L2200: */
#line 621 "zgsvj1.f"
			}
/*        end of the q-loop */
#line 623 "zgsvj1.f"
L2203:

#line 625 "zgsvj1.f"
			sva[p] = aapp;

#line 627 "zgsvj1.f"
		    } else {

#line 629 "zgsvj1.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 629 "zgsvj1.f"
			    i__5 = jgl + kbl - 1;
#line 629 "zgsvj1.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 629 "zgsvj1.f"
			}
#line 631 "zgsvj1.f"
			if (aapp < 0.) {
#line 631 "zgsvj1.f"
			    notrot = 0;
#line 631 "zgsvj1.f"
			}

#line 633 "zgsvj1.f"
		    }

#line 635 "zgsvj1.f"
/* L2100: */
#line 635 "zgsvj1.f"
		}
/*     end of the p-loop */
#line 637 "zgsvj1.f"
/* L2010: */
#line 637 "zgsvj1.f"
	    }
/*     end of the jbc-loop */
#line 639 "zgsvj1.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 641 "zgsvj1.f"
	    i__4 = igl + kbl - 1;
#line 641 "zgsvj1.f"
	    i__3 = min(i__4,*n);
#line 641 "zgsvj1.f"
	    for (p = igl; p <= i__3; ++p) {
#line 642 "zgsvj1.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 643 "zgsvj1.f"
/* L2012: */
#line 643 "zgsvj1.f"
	    }
/* ** */
#line 645 "zgsvj1.f"
/* L2000: */
#line 645 "zgsvj1.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 649 "zgsvj1.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 651 "zgsvj1.f"
	    sva[*n] = dznrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 652 "zgsvj1.f"
	} else {
#line 653 "zgsvj1.f"
	    t = 0.;
#line 654 "zgsvj1.f"
	    aapp = 1.;
#line 655 "zgsvj1.f"
	    zlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 656 "zgsvj1.f"
	    sva[*n] = t * sqrt(aapp);
#line 657 "zgsvj1.f"
	}

/*     Additional steering devices */

#line 661 "zgsvj1.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 661 "zgsvj1.f"
	    swband = i__;
#line 661 "zgsvj1.f"
	}

#line 664 "zgsvj1.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 666 "zgsvj1.f"
	    goto L1994;
#line 667 "zgsvj1.f"
	}

#line 669 "zgsvj1.f"
	if (notrot >= emptsw) {
#line 669 "zgsvj1.f"
	    goto L1994;
#line 669 "zgsvj1.f"
	}

#line 671 "zgsvj1.f"
/* L1993: */
#line 671 "zgsvj1.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 675 "zgsvj1.f"
    *info = *nsweep - 1;
#line 676 "zgsvj1.f"
    goto L1995;

#line 678 "zgsvj1.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 682 "zgsvj1.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 684 "zgsvj1.f"
L1995:

/*     Sort the vector SVA() of column norms. */
#line 687 "zgsvj1.f"
    i__1 = *n - 1;
#line 687 "zgsvj1.f"
    for (p = 1; p <= i__1; ++p) {
#line 688 "zgsvj1.f"
	i__2 = *n - p + 1;
#line 688 "zgsvj1.f"
	q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 689 "zgsvj1.f"
	if (p != q) {
#line 690 "zgsvj1.f"
	    temp1 = sva[p];
#line 691 "zgsvj1.f"
	    sva[p] = sva[q];
#line 692 "zgsvj1.f"
	    sva[q] = temp1;
#line 693 "zgsvj1.f"
	    i__2 = p;
#line 693 "zgsvj1.f"
	    aapq.r = d__[i__2].r, aapq.i = d__[i__2].i;
#line 694 "zgsvj1.f"
	    i__2 = p;
#line 694 "zgsvj1.f"
	    i__3 = q;
#line 694 "zgsvj1.f"
	    d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
#line 695 "zgsvj1.f"
	    i__2 = q;
#line 695 "zgsvj1.f"
	    d__[i__2].r = aapq.r, d__[i__2].i = aapq.i;
#line 696 "zgsvj1.f"
	    zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 697 "zgsvj1.f"
	    if (rsvec) {
#line 697 "zgsvj1.f"
		zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 697 "zgsvj1.f"
	    }
#line 698 "zgsvj1.f"
	}
#line 699 "zgsvj1.f"
/* L5991: */
#line 699 "zgsvj1.f"
    }


#line 702 "zgsvj1.f"
    return 0;
/*     .. */
/*     .. END OF ZGSVJ1 */
/*     .. */
} /* zgsvj1_ */

