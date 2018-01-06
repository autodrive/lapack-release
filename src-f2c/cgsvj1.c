#line 1 "cgsvj1.f"
/* cgsvj1.f -- translated by f2c (version 20100827).
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

#line 1 "cgsvj1.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b18 = 1.;

/* > \brief \b CGSVJ1 pre-processor for the routine sgesvj, applies Jacobi rotations targeting only particular
 pivots. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGSVJ1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgsvj1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgsvj1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgsvj1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, */
/*                          EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       REAL               EPS, SFMIN, TOL */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP */
/*       CHARACTER*1        JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX        A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK ) */
/*       REAL           SVA( N ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGSVJ1 is called from CGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as CGESVJ does, but */
/* > it targets only particular pivots and it does not check convergence */
/* > (stopping criterion). Few tunning parameters (marked by [TP]) are */
/* > available for the implementer. */
/* > */
/* > Further Details */
/* > ~~~~~~~~~~~~~~~ */
/* > CGSVJ1 applies few sweeps of Jacobi rotations in the column space of */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          D is REAL array, dimension (N) */
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
/* >          SVA is REAL array, dimension (N) */
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
/* >          V is REAL array, dimension (LDV,N) */
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
/* >          EPS is REAL */
/* >          EPS = SLAMCH('Epsilon') */
/* > \endverbatim */
/* > */
/* > \param[in] SFMIN */
/* > \verbatim */
/* >          SFMIN is REAL */
/* >          SFMIN = SLAMCH('Safe Minimum') */
/* > \endverbatim */
/* > */
/* > \param[in] TOL */
/* > \verbatim */
/* >          TOL is REAL */
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
/* >         WORK is COMPLEX array, dimension LWORK. */
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

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */

/*  ===================================================================== */
/* Subroutine */ int cgsvj1_(char *jobv, integer *m, integer *n, integer *n1, 
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
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublecomplex ompq;
    static doublereal aapp0, aapq1, temp1, large;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta, small;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static logical applv, rsvec;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical rotok;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen)
	    ;
    static integer ijblsk, swband;
    extern integer isamax_(integer *, doublereal *, integer *);
    static integer blskip;
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    static doublereal mxaapq, thsign, mxsinj;
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

#line 292 "cgsvj1.f"
    /* Parameter adjustments */
#line 292 "cgsvj1.f"
    --sva;
#line 292 "cgsvj1.f"
    --d__;
#line 292 "cgsvj1.f"
    a_dim1 = *lda;
#line 292 "cgsvj1.f"
    a_offset = 1 + a_dim1;
#line 292 "cgsvj1.f"
    a -= a_offset;
#line 292 "cgsvj1.f"
    v_dim1 = *ldv;
#line 292 "cgsvj1.f"
    v_offset = 1 + v_dim1;
#line 292 "cgsvj1.f"
    v -= v_offset;
#line 292 "cgsvj1.f"
    --work;
#line 292 "cgsvj1.f"

#line 292 "cgsvj1.f"
    /* Function Body */
#line 292 "cgsvj1.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 293 "cgsvj1.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 294 "cgsvj1.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 295 "cgsvj1.f"
	*info = -1;
#line 296 "cgsvj1.f"
    } else if (*m < 0) {
#line 297 "cgsvj1.f"
	*info = -2;
#line 298 "cgsvj1.f"
    } else if (*n < 0 || *n > *m) {
#line 299 "cgsvj1.f"
	*info = -3;
#line 300 "cgsvj1.f"
    } else if (*n1 < 0) {
#line 301 "cgsvj1.f"
	*info = -4;
#line 302 "cgsvj1.f"
    } else if (*lda < *m) {
#line 303 "cgsvj1.f"
	*info = -6;
#line 304 "cgsvj1.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 305 "cgsvj1.f"
	*info = -9;
#line 306 "cgsvj1.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 308 "cgsvj1.f"
	*info = -11;
#line 309 "cgsvj1.f"
    } else if (*tol <= *eps) {
#line 310 "cgsvj1.f"
	*info = -14;
#line 311 "cgsvj1.f"
    } else if (*nsweep < 0) {
#line 312 "cgsvj1.f"
	*info = -15;
#line 313 "cgsvj1.f"
    } else if (*lwork < *m) {
#line 314 "cgsvj1.f"
	*info = -17;
#line 315 "cgsvj1.f"
    } else {
#line 316 "cgsvj1.f"
	*info = 0;
#line 317 "cgsvj1.f"
    }

/*     #:( */
#line 320 "cgsvj1.f"
    if (*info != 0) {
#line 321 "cgsvj1.f"
	i__1 = -(*info);
#line 321 "cgsvj1.f"
	xerbla_("CGSVJ1", &i__1, (ftnlen)6);
#line 322 "cgsvj1.f"
	return 0;
#line 323 "cgsvj1.f"
    }

#line 325 "cgsvj1.f"
    if (rsvec) {
#line 326 "cgsvj1.f"
	mvl = *n;
#line 327 "cgsvj1.f"
    } else if (applv) {
#line 328 "cgsvj1.f"
	mvl = *mv;
#line 329 "cgsvj1.f"
    }
#line 330 "cgsvj1.f"
    rsvec = rsvec || applv;
#line 332 "cgsvj1.f"
    rooteps = sqrt(*eps);
#line 333 "cgsvj1.f"
    rootsfmin = sqrt(*sfmin);
#line 334 "cgsvj1.f"
    small = *sfmin / *eps;
#line 335 "cgsvj1.f"
    big = 1. / *sfmin;
#line 336 "cgsvj1.f"
    rootbig = 1. / rootsfmin;
#line 337 "cgsvj1.f"
    large = big / sqrt((doublereal) (*m * *n));
#line 338 "cgsvj1.f"
    bigtheta = 1. / rooteps;
#line 339 "cgsvj1.f"
    roottol = sqrt(*tol);

/*     .. Initialize the right singular vector matrix .. */

/*     RSVEC = LSAME( JOBV, 'Y' ) */

#line 345 "cgsvj1.f"
    emptsw = *n1 * (*n - *n1);
#line 346 "cgsvj1.f"
    notrot = 0;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 350 "cgsvj1.f"
    kbl = min(8,*n);
#line 351 "cgsvj1.f"
    nblr = *n1 / kbl;
#line 352 "cgsvj1.f"
    if (nblr * kbl != *n1) {
#line 352 "cgsvj1.f"
	++nblr;
#line 352 "cgsvj1.f"
    }
/*     .. the tiling is nblr-by-nblc [tiles] */
#line 356 "cgsvj1.f"
    nblc = (*n - *n1) / kbl;
#line 357 "cgsvj1.f"
    if (nblc * kbl != *n - *n1) {
#line 357 "cgsvj1.f"
	++nblc;
#line 357 "cgsvj1.f"
    }
/* Computing 2nd power */
#line 358 "cgsvj1.f"
    i__1 = kbl;
#line 358 "cgsvj1.f"
    blskip = i__1 * i__1 + 1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
#line 361 "cgsvj1.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */
#line 363 "cgsvj1.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter. It is meaningful and effective */
/*     if CGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm CGEJSV. */


/*     | *   *   * [x] [x] [x]| */
/*     | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks. */
/*     | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block. */
/*     |[x] [x] [x] *   *   * | */
/*     |[x] [x] [x] *   *   * | */
/*     |[x] [x] [x] *   *   * | */


#line 377 "cgsvj1.f"
    i__1 = *nsweep;
#line 377 "cgsvj1.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     .. go go go ... */

#line 381 "cgsvj1.f"
	mxaapq = 0.;
#line 382 "cgsvj1.f"
	mxsinj = 0.;
#line 383 "cgsvj1.f"
	iswrot = 0;

#line 385 "cgsvj1.f"
	notrot = 0;
#line 386 "cgsvj1.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 393 "cgsvj1.f"
	i__2 = nblr;
#line 393 "cgsvj1.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {

#line 395 "cgsvj1.f"
	    igl = (ibr - 1) * kbl + 1;


/* ... go to the off diagonal blocks */

#line 401 "cgsvj1.f"
	    igl = (ibr - 1) * kbl + 1;

/*            DO 2010 jbc = ibr + 1, NBL */
#line 404 "cgsvj1.f"
	    i__3 = nblc;
#line 404 "cgsvj1.f"
	    for (jbc = 1; jbc <= i__3; ++jbc) {

#line 406 "cgsvj1.f"
		jgl = (jbc - 1) * kbl + *n1 + 1;

/*        doing the block at ( ibr, jbc ) */

#line 410 "cgsvj1.f"
		ijblsk = 0;
/* Computing MIN */
#line 411 "cgsvj1.f"
		i__5 = igl + kbl - 1;
#line 411 "cgsvj1.f"
		i__4 = min(i__5,*n1);
#line 411 "cgsvj1.f"
		for (p = igl; p <= i__4; ++p) {

#line 413 "cgsvj1.f"
		    aapp = sva[p];
#line 414 "cgsvj1.f"
		    if (aapp > 0.) {

#line 416 "cgsvj1.f"
			pskipped = 0;

/* Computing MIN */
#line 418 "cgsvj1.f"
			i__6 = jgl + kbl - 1;
#line 418 "cgsvj1.f"
			i__5 = min(i__6,*n);
#line 418 "cgsvj1.f"
			for (q = jgl; q <= i__5; ++q) {

#line 420 "cgsvj1.f"
			    aaqq = sva[q];
#line 421 "cgsvj1.f"
			    if (aaqq > 0.) {
#line 422 "cgsvj1.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 428 "cgsvj1.f"
				if (aaqq >= 1.) {
#line 429 "cgsvj1.f"
				    if (aapp >= aaqq) {
#line 430 "cgsvj1.f"
					rotok = small * aapp <= aaqq;
#line 431 "cgsvj1.f"
				    } else {
#line 432 "cgsvj1.f"
					rotok = small * aaqq <= aapp;
#line 433 "cgsvj1.f"
				    }
#line 434 "cgsvj1.f"
				    if (aapp < big / aaqq) {
#line 435 "cgsvj1.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 435 "cgsvj1.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 435 "cgsvj1.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 435 "cgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 437 "cgsvj1.f"
				    } else {
#line 438 "cgsvj1.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 440 "cgsvj1.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b18, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 443 "cgsvj1.f"
					cdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 443 "cgsvj1.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 443 "cgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 445 "cgsvj1.f"
				    }
#line 446 "cgsvj1.f"
				} else {
#line 447 "cgsvj1.f"
				    if (aapp >= aaqq) {
#line 448 "cgsvj1.f"
					rotok = aapp <= aaqq / small;
#line 449 "cgsvj1.f"
				    } else {
#line 450 "cgsvj1.f"
					rotok = aaqq <= aapp / small;
#line 451 "cgsvj1.f"
				    }
#line 452 "cgsvj1.f"
				    if (aapp > small / aaqq) {
#line 453 "cgsvj1.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 453 "cgsvj1.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 453 "cgsvj1.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 453 "cgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 455 "cgsvj1.f"
				    } else {
#line 456 "cgsvj1.f"
					ccopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 458 "cgsvj1.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b18, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 461 "cgsvj1.f"
					cdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 461 "cgsvj1.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 461 "cgsvj1.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 463 "cgsvj1.f"
				    }
#line 464 "cgsvj1.f"
				}

#line 466 "cgsvj1.f"
				d__1 = z_abs(&aapq);
#line 466 "cgsvj1.f"
				z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					d__1;
#line 466 "cgsvj1.f"
				ompq.r = z__1.r, ompq.i = z__1.i;
/*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
#line 468 "cgsvj1.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 469 "cgsvj1.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 469 "cgsvj1.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 473 "cgsvj1.f"
				if (abs(aapq1) > *tol) {
#line 474 "cgsvj1.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 476 "cgsvj1.f"
				    pskipped = 0;
#line 477 "cgsvj1.f"
				    ++iswrot;

#line 479 "cgsvj1.f"
				    if (rotok) {

#line 481 "cgsvj1.f"
					aqoap = aaqq / aapp;
#line 482 "cgsvj1.f"
					apoaq = aapp / aaqq;
#line 483 "cgsvj1.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 484 "cgsvj1.f"
					if (aaqq > aapp0) {
#line 484 "cgsvj1.f"
					    theta = -theta;
#line 484 "cgsvj1.f"
					}

#line 486 "cgsvj1.f"
					if (abs(theta) > bigtheta) {
#line 487 "cgsvj1.f"
					    t = .5 / theta;
#line 488 "cgsvj1.f"
					    cs = 1.;
#line 489 "cgsvj1.f"
					    d_cnjg(&z__2, &ompq);
#line 489 "cgsvj1.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 489 "cgsvj1.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 491 "cgsvj1.f"
					    if (rsvec) {
#line 492 "cgsvj1.f"
			  d_cnjg(&z__2, &ompq);
#line 492 "cgsvj1.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 492 "cgsvj1.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 494 "cgsvj1.f"
					    }
/* Computing MAX */
#line 495 "cgsvj1.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 495 "cgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 497 "cgsvj1.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 497 "cgsvj1.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 499 "cgsvj1.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 499 "cgsvj1.f"
					    mxsinj = max(d__1,d__2);
#line 500 "cgsvj1.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 504 "cgsvj1.f"
					    thsign = -d_sign(&c_b18, &aapq1);
#line 505 "cgsvj1.f"
					    if (aaqq > aapp0) {
#line 505 "cgsvj1.f"
			  thsign = -thsign;
#line 505 "cgsvj1.f"
					    }
#line 506 "cgsvj1.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 508 "cgsvj1.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 509 "cgsvj1.f"
					    sn = t * cs;
/* Computing MAX */
#line 510 "cgsvj1.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 510 "cgsvj1.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 511 "cgsvj1.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 511 "cgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 513 "cgsvj1.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 513 "cgsvj1.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 516 "cgsvj1.f"
					    d_cnjg(&z__2, &ompq);
#line 516 "cgsvj1.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 516 "cgsvj1.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 518 "cgsvj1.f"
					    if (rsvec) {
#line 519 "cgsvj1.f"
			  d_cnjg(&z__2, &ompq);
#line 519 "cgsvj1.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 519 "cgsvj1.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 521 "cgsvj1.f"
					    }
#line 522 "cgsvj1.f"
					}
#line 523 "cgsvj1.f"
					i__6 = p;
#line 523 "cgsvj1.f"
					i__7 = q;
#line 523 "cgsvj1.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 523 "cgsvj1.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 523 "cgsvj1.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 525 "cgsvj1.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 527 "cgsvj1.f"
					if (aapp > aaqq) {
#line 528 "cgsvj1.f"
					    ccopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 530 "cgsvj1.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 533 "cgsvj1.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 536 "cgsvj1.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 536 "cgsvj1.f"
					    caxpy_(m, &z__1, &work[1], &c__1, 
						    &a[q * a_dim1 + 1], &c__1)
						    ;
#line 538 "cgsvj1.f"
					    clascl_("G", &c__0, &c__0, &c_b18,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 541 "cgsvj1.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 541 "cgsvj1.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 543 "cgsvj1.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 544 "cgsvj1.f"
					} else {
#line 545 "cgsvj1.f"
					    ccopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 547 "cgsvj1.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 550 "cgsvj1.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 553 "cgsvj1.f"
					    d_cnjg(&z__2, &aapq);
#line 553 "cgsvj1.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 553 "cgsvj1.f"
					    caxpy_(m, &z__1, &work[1], &c__1, 
						    &a[p * a_dim1 + 1], &c__1)
						    ;
#line 555 "cgsvj1.f"
					    clascl_("G", &c__0, &c__0, &c_b18,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 558 "cgsvj1.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 558 "cgsvj1.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 560 "cgsvj1.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 561 "cgsvj1.f"
					}
#line 562 "cgsvj1.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 567 "cgsvj1.f"
				    d__1 = sva[q] / aaqq;
#line 567 "cgsvj1.f"
				    if (d__1 * d__1 <= rooteps) {
#line 569 "cgsvj1.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 571 "cgsvj1.f"
					    sva[q] = scnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 572 "cgsvj1.f"
					} else {
#line 573 "cgsvj1.f"
					    t = 0.;
#line 574 "cgsvj1.f"
					    aaqq = 1.;
#line 575 "cgsvj1.f"
					    classq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 577 "cgsvj1.f"
					    sva[q] = t * sqrt(aaqq);
#line 578 "cgsvj1.f"
					}
#line 579 "cgsvj1.f"
				    }
/* Computing 2nd power */
#line 580 "cgsvj1.f"
				    d__1 = aapp / aapp0;
#line 580 "cgsvj1.f"
				    if (d__1 * d__1 <= rooteps) {
#line 581 "cgsvj1.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 583 "cgsvj1.f"
					    aapp = scnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 584 "cgsvj1.f"
					} else {
#line 585 "cgsvj1.f"
					    t = 0.;
#line 586 "cgsvj1.f"
					    aapp = 1.;
#line 587 "cgsvj1.f"
					    classq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 589 "cgsvj1.f"
					    aapp = t * sqrt(aapp);
#line 590 "cgsvj1.f"
					}
#line 591 "cgsvj1.f"
					sva[p] = aapp;
#line 592 "cgsvj1.f"
				    }
/*              end of OK rotation */
#line 594 "cgsvj1.f"
				} else {
#line 595 "cgsvj1.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 597 "cgsvj1.f"
				    ++pskipped;
#line 598 "cgsvj1.f"
				    ++ijblsk;
#line 599 "cgsvj1.f"
				}
#line 600 "cgsvj1.f"
			    } else {
#line 601 "cgsvj1.f"
				++notrot;
#line 602 "cgsvj1.f"
				++pskipped;
#line 603 "cgsvj1.f"
				++ijblsk;
#line 604 "cgsvj1.f"
			    }

#line 606 "cgsvj1.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 608 "cgsvj1.f"
				sva[p] = aapp;
#line 609 "cgsvj1.f"
				notrot = 0;
#line 610 "cgsvj1.f"
				goto L2011;
#line 611 "cgsvj1.f"
			    }
#line 612 "cgsvj1.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 614 "cgsvj1.f"
				aapp = -aapp;
#line 615 "cgsvj1.f"
				notrot = 0;
#line 616 "cgsvj1.f"
				goto L2203;
#line 617 "cgsvj1.f"
			    }

#line 619 "cgsvj1.f"
/* L2200: */
#line 619 "cgsvj1.f"
			}
/*        end of the q-loop */
#line 621 "cgsvj1.f"
L2203:

#line 623 "cgsvj1.f"
			sva[p] = aapp;

#line 625 "cgsvj1.f"
		    } else {

#line 627 "cgsvj1.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 627 "cgsvj1.f"
			    i__5 = jgl + kbl - 1;
#line 627 "cgsvj1.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 627 "cgsvj1.f"
			}
#line 629 "cgsvj1.f"
			if (aapp < 0.) {
#line 629 "cgsvj1.f"
			    notrot = 0;
#line 629 "cgsvj1.f"
			}

#line 631 "cgsvj1.f"
		    }

#line 633 "cgsvj1.f"
/* L2100: */
#line 633 "cgsvj1.f"
		}
/*     end of the p-loop */
#line 635 "cgsvj1.f"
/* L2010: */
#line 635 "cgsvj1.f"
	    }
/*     end of the jbc-loop */
#line 637 "cgsvj1.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 639 "cgsvj1.f"
	    i__4 = igl + kbl - 1;
#line 639 "cgsvj1.f"
	    i__3 = min(i__4,*n);
#line 639 "cgsvj1.f"
	    for (p = igl; p <= i__3; ++p) {
#line 640 "cgsvj1.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 641 "cgsvj1.f"
/* L2012: */
#line 641 "cgsvj1.f"
	    }
/* ** */
#line 643 "cgsvj1.f"
/* L2000: */
#line 643 "cgsvj1.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 647 "cgsvj1.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 649 "cgsvj1.f"
	    sva[*n] = scnrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 650 "cgsvj1.f"
	} else {
#line 651 "cgsvj1.f"
	    t = 0.;
#line 652 "cgsvj1.f"
	    aapp = 1.;
#line 653 "cgsvj1.f"
	    classq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 654 "cgsvj1.f"
	    sva[*n] = t * sqrt(aapp);
#line 655 "cgsvj1.f"
	}

/*     Additional steering devices */

#line 659 "cgsvj1.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 659 "cgsvj1.f"
	    swband = i__;
#line 659 "cgsvj1.f"
	}

#line 662 "cgsvj1.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 664 "cgsvj1.f"
	    goto L1994;
#line 665 "cgsvj1.f"
	}

#line 667 "cgsvj1.f"
	if (notrot >= emptsw) {
#line 667 "cgsvj1.f"
	    goto L1994;
#line 667 "cgsvj1.f"
	}

#line 669 "cgsvj1.f"
/* L1993: */
#line 669 "cgsvj1.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 673 "cgsvj1.f"
    *info = *nsweep - 1;
#line 674 "cgsvj1.f"
    goto L1995;

#line 676 "cgsvj1.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 680 "cgsvj1.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 682 "cgsvj1.f"
L1995:

/*     Sort the vector SVA() of column norms. */
#line 685 "cgsvj1.f"
    i__1 = *n - 1;
#line 685 "cgsvj1.f"
    for (p = 1; p <= i__1; ++p) {
#line 686 "cgsvj1.f"
	i__2 = *n - p + 1;
#line 686 "cgsvj1.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 687 "cgsvj1.f"
	if (p != q) {
#line 688 "cgsvj1.f"
	    temp1 = sva[p];
#line 689 "cgsvj1.f"
	    sva[p] = sva[q];
#line 690 "cgsvj1.f"
	    sva[q] = temp1;
#line 691 "cgsvj1.f"
	    i__2 = p;
#line 691 "cgsvj1.f"
	    aapq.r = d__[i__2].r, aapq.i = d__[i__2].i;
#line 692 "cgsvj1.f"
	    i__2 = p;
#line 692 "cgsvj1.f"
	    i__3 = q;
#line 692 "cgsvj1.f"
	    d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
#line 693 "cgsvj1.f"
	    i__2 = q;
#line 693 "cgsvj1.f"
	    d__[i__2].r = aapq.r, d__[i__2].i = aapq.i;
#line 694 "cgsvj1.f"
	    cswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 695 "cgsvj1.f"
	    if (rsvec) {
#line 695 "cgsvj1.f"
		cswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 695 "cgsvj1.f"
	    }
#line 696 "cgsvj1.f"
	}
#line 697 "cgsvj1.f"
/* L5991: */
#line 697 "cgsvj1.f"
    }


#line 700 "cgsvj1.f"
    return 0;
/*     .. */
/*     .. END OF CGSVJ1 */
/*     .. */
} /* cgsvj1_ */

