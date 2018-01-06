#line 1 "zgsvj0.f"
/* zgsvj0.f -- translated by f2c (version 20100827).
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

#line 1 "zgsvj0.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b27 = 1.;

/* > \brief \b ZGSVJ0 pre-processor for the routine dgesvj. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGSVJ0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgsvj0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgsvj0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgsvj0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, */
/*                          SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP */
/*       DOUBLE PRECISION   EPS, SFMIN, TOL */
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
/* > ZGSVJ0 is called from ZGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as ZGESVJ does, but */
/* > it does not check convergence (stopping criterion). Few tuning */
/* > parameters (marked by [TP]) are available for the implementer. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, M-by-N matrix A, such that A*diag(D) represents */
/* >          the input matrix. */
/* >          On exit, */
/* >          A_onexit * diag(D_onexit) represents the input matrix A*diag(D) */
/* >          post-multiplied by a sequence of Jacobi rotations, where the */
/* >          rotation threshold and the total number of sweeps are given in */
/* >          TOL and NSWEEP, respectively. */
/* >          (See the descriptions of D, TOL and NSWEEP.) */
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
/* >          The array D accumulates the scaling factors from the complex scaled */
/* >          Jacobi rotations. */
/* >          On entry, A*diag(D) represents the input matrix. */
/* >          On exit, A_onexit*diag(D_onexit) represents the input matrix */
/* >          post-multiplied by a sequence of Jacobi rotations, where the */
/* >          rotation threshold and the total number of sweeps are given in */
/* >          TOL and NSWEEP, respectively. */
/* >          (See the descriptions of A, TOL and NSWEEP.) */
/* > \endverbatim */
/* > */
/* > \param[in,out] SVA */
/* > \verbatim */
/* >          SVA is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, SVA contains the Euclidean norms of the columns of */
/* >          the matrix A*diag(D). */
/* >          On exit, SVA contains the Euclidean norms of the columns of */
/* >          the matrix A_onexit*diag(D_onexit). */
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
/* >          WORK is COMPLEX*16 array, dimension LWORK. */
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

/* > \ingroup complex16OTHERcomputational */
/* > */
/* > \par Further Details: */
/*  ===================== */
/* > */
/* > ZGSVJ0 is used just to enable ZGESVJ to call a simplified version of */
/* > itself to work on a submatrix of the original matrix. */
/* > */
/* > Contributors: */
/* ============= */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */
/* > */
/* > Bugs, Examples and Comments: */
/* ============================ */
/* > */
/* > Please report all bugs and send interesting test examples and comments to */
/* > drmac@math.hr. Thank you. */

/*  ===================================================================== */
/* Subroutine */ int zgsvj0_(char *jobv, integer *m, integer *n, 
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
    static integer ir1, jbc;
    static doublereal big;
    static integer kbl, igl, ibr, jgl, nbl, mvl;
    static doublereal aapp;
    static doublecomplex aapq;
    static doublereal aaqq;
    static integer ierr;
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
    static integer emptsw, notrot, iswrot, lkahead;
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
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     from BLAS */
/*     from LAPACK */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 278 "zgsvj0.f"
    /* Parameter adjustments */
#line 278 "zgsvj0.f"
    --sva;
#line 278 "zgsvj0.f"
    --d__;
#line 278 "zgsvj0.f"
    a_dim1 = *lda;
#line 278 "zgsvj0.f"
    a_offset = 1 + a_dim1;
#line 278 "zgsvj0.f"
    a -= a_offset;
#line 278 "zgsvj0.f"
    v_dim1 = *ldv;
#line 278 "zgsvj0.f"
    v_offset = 1 + v_dim1;
#line 278 "zgsvj0.f"
    v -= v_offset;
#line 278 "zgsvj0.f"
    --work;
#line 278 "zgsvj0.f"

#line 278 "zgsvj0.f"
    /* Function Body */
#line 278 "zgsvj0.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 279 "zgsvj0.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 280 "zgsvj0.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 281 "zgsvj0.f"
	*info = -1;
#line 282 "zgsvj0.f"
    } else if (*m < 0) {
#line 283 "zgsvj0.f"
	*info = -2;
#line 284 "zgsvj0.f"
    } else if (*n < 0 || *n > *m) {
#line 285 "zgsvj0.f"
	*info = -3;
#line 286 "zgsvj0.f"
    } else if (*lda < *m) {
#line 287 "zgsvj0.f"
	*info = -5;
#line 288 "zgsvj0.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 289 "zgsvj0.f"
	*info = -8;
#line 290 "zgsvj0.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 292 "zgsvj0.f"
	*info = -10;
#line 293 "zgsvj0.f"
    } else if (*tol <= *eps) {
#line 294 "zgsvj0.f"
	*info = -13;
#line 295 "zgsvj0.f"
    } else if (*nsweep < 0) {
#line 296 "zgsvj0.f"
	*info = -14;
#line 297 "zgsvj0.f"
    } else if (*lwork < *m) {
#line 298 "zgsvj0.f"
	*info = -16;
#line 299 "zgsvj0.f"
    } else {
#line 300 "zgsvj0.f"
	*info = 0;
#line 301 "zgsvj0.f"
    }

/*     #:( */
#line 304 "zgsvj0.f"
    if (*info != 0) {
#line 305 "zgsvj0.f"
	i__1 = -(*info);
#line 305 "zgsvj0.f"
	xerbla_("ZGSVJ0", &i__1, (ftnlen)6);
#line 306 "zgsvj0.f"
	return 0;
#line 307 "zgsvj0.f"
    }

#line 309 "zgsvj0.f"
    if (rsvec) {
#line 310 "zgsvj0.f"
	mvl = *n;
#line 311 "zgsvj0.f"
    } else if (applv) {
#line 312 "zgsvj0.f"
	mvl = *mv;
#line 313 "zgsvj0.f"
    }
#line 314 "zgsvj0.f"
    rsvec = rsvec || applv;
#line 316 "zgsvj0.f"
    rooteps = sqrt(*eps);
#line 317 "zgsvj0.f"
    rootsfmin = sqrt(*sfmin);
#line 318 "zgsvj0.f"
    small = *sfmin / *eps;
#line 319 "zgsvj0.f"
    big = 1. / *sfmin;
#line 320 "zgsvj0.f"
    rootbig = 1. / rootsfmin;
#line 321 "zgsvj0.f"
    bigtheta = 1. / rooteps;
#line 322 "zgsvj0.f"
    roottol = sqrt(*tol);

/*     .. Row-cyclic Jacobi SVD algorithm with column pivoting .. */

#line 326 "zgsvj0.f"
    emptsw = *n * (*n - 1) / 2;
#line 327 "zgsvj0.f"
    notrot = 0;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 332 "zgsvj0.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if ZGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 340 "zgsvj0.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 346 "zgsvj0.f"
    nbl = *n / kbl;
#line 347 "zgsvj0.f"
    if (nbl * kbl != *n) {
#line 347 "zgsvj0.f"
	++nbl;
#line 347 "zgsvj0.f"
    }

/* Computing 2nd power */
#line 349 "zgsvj0.f"
    i__1 = kbl;
#line 349 "zgsvj0.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 352 "zgsvj0.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 355 "zgsvj0.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */


/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 366 "zgsvj0.f"
    i__1 = *nsweep;
#line 366 "zgsvj0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     .. go go go ... */

#line 370 "zgsvj0.f"
	mxaapq = 0.;
#line 371 "zgsvj0.f"
	mxsinj = 0.;
#line 372 "zgsvj0.f"
	iswrot = 0;

#line 374 "zgsvj0.f"
	notrot = 0;
#line 375 "zgsvj0.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 382 "zgsvj0.f"
	i__2 = nbl;
#line 382 "zgsvj0.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {

#line 384 "zgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 386 "zgsvj0.f"
	    i__4 = lkahead, i__5 = nbl - ibr;
#line 386 "zgsvj0.f"
	    i__3 = min(i__4,i__5);
#line 386 "zgsvj0.f"
	    for (ir1 = 0; ir1 <= i__3; ++ir1) {

#line 388 "zgsvj0.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 390 "zgsvj0.f"
		i__5 = igl + kbl - 1, i__6 = *n - 1;
#line 390 "zgsvj0.f"
		i__4 = min(i__5,i__6);
#line 390 "zgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

/*     .. de Rijk's pivoting */

#line 394 "zgsvj0.f"
		    i__5 = *n - p + 1;
#line 394 "zgsvj0.f"
		    q = idamax_(&i__5, &sva[p], &c__1) + p - 1;
#line 395 "zgsvj0.f"
		    if (p != q) {
#line 396 "zgsvj0.f"
			zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 397 "zgsvj0.f"
			if (rsvec) {
#line 397 "zgsvj0.f"
			    zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 397 "zgsvj0.f"
			}
#line 399 "zgsvj0.f"
			temp1 = sva[p];
#line 400 "zgsvj0.f"
			sva[p] = sva[q];
#line 401 "zgsvj0.f"
			sva[q] = temp1;
#line 402 "zgsvj0.f"
			i__5 = p;
#line 402 "zgsvj0.f"
			aapq.r = d__[i__5].r, aapq.i = d__[i__5].i;
#line 403 "zgsvj0.f"
			i__5 = p;
#line 403 "zgsvj0.f"
			i__6 = q;
#line 403 "zgsvj0.f"
			d__[i__5].r = d__[i__6].r, d__[i__5].i = d__[i__6].i;
#line 404 "zgsvj0.f"
			i__5 = q;
#line 404 "zgsvj0.f"
			d__[i__5].r = aapq.r, d__[i__5].i = aapq.i;
#line 405 "zgsvj0.f"
		    }

#line 407 "zgsvj0.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/*        Caveat: */
/*        Unfortunately, some BLAS implementations compute SNCRM2(M,A(1,p),1) */
/*        as SQRT(S=ZDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to */
/*        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to */
/*        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
/*        Hence, DZNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented DZNRM2 is available, the IF-THEN-ELSE-END IF */
/*        below should be replaced with "AAPP = DZNRM2( M, A(1,p), 1 )". */

#line 421 "zgsvj0.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 423 "zgsvj0.f"
			    sva[p] = dznrm2_(m, &a[p * a_dim1 + 1], &c__1);
#line 424 "zgsvj0.f"
			} else {
#line 425 "zgsvj0.f"
			    temp1 = 0.;
#line 426 "zgsvj0.f"
			    aapp = 1.;
#line 427 "zgsvj0.f"
			    zlassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 428 "zgsvj0.f"
			    sva[p] = temp1 * sqrt(aapp);
#line 429 "zgsvj0.f"
			}
#line 430 "zgsvj0.f"
			aapp = sva[p];
#line 431 "zgsvj0.f"
		    } else {
#line 432 "zgsvj0.f"
			aapp = sva[p];
#line 433 "zgsvj0.f"
		    }

#line 435 "zgsvj0.f"
		    if (aapp > 0.) {

#line 437 "zgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 439 "zgsvj0.f"
			i__6 = igl + kbl - 1;
#line 439 "zgsvj0.f"
			i__5 = min(i__6,*n);
#line 439 "zgsvj0.f"
			for (q = p + 1; q <= i__5; ++q) {

#line 441 "zgsvj0.f"
			    aaqq = sva[q];

#line 443 "zgsvj0.f"
			    if (aaqq > 0.) {

#line 445 "zgsvj0.f"
				aapp0 = aapp;
#line 446 "zgsvj0.f"
				if (aaqq >= 1.) {
#line 447 "zgsvj0.f"
				    rotok = small * aapp <= aaqq;
#line 448 "zgsvj0.f"
				    if (aapp < big / aaqq) {
#line 449 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 449 "zgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 449 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 449 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 451 "zgsvj0.f"
				    } else {
#line 452 "zgsvj0.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 454 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 456 "zgsvj0.f"
					zdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 456 "zgsvj0.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 456 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 458 "zgsvj0.f"
				    }
#line 459 "zgsvj0.f"
				} else {
#line 460 "zgsvj0.f"
				    rotok = aapp <= aaqq / small;
#line 461 "zgsvj0.f"
				    if (aapp > small / aaqq) {
#line 462 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 462 "zgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 462 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 462 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 464 "zgsvj0.f"
				    } else {
#line 465 "zgsvj0.f"
					zcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 467 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 470 "zgsvj0.f"
					zdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 470 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 470 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 472 "zgsvj0.f"
				    }
#line 473 "zgsvj0.f"
				}

#line 475 "zgsvj0.f"
				d__1 = z_abs(&aapq);
#line 475 "zgsvj0.f"
				z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					d__1;
#line 475 "zgsvj0.f"
				ompq.r = z__1.r, ompq.i = z__1.i;
/*                           AAPQ = AAPQ * DCONJG( CWORK(p) ) * CWORK(q) */
#line 477 "zgsvj0.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 478 "zgsvj0.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 478 "zgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 482 "zgsvj0.f"
				if (abs(aapq1) > *tol) {

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 487 "zgsvj0.f"
				    if (ir1 == 0) {
#line 488 "zgsvj0.f"
					notrot = 0;
#line 489 "zgsvj0.f"
					pskipped = 0;
#line 490 "zgsvj0.f"
					++iswrot;
#line 491 "zgsvj0.f"
				    }

#line 493 "zgsvj0.f"
				    if (rotok) {

#line 495 "zgsvj0.f"
					aqoap = aaqq / aapp;
#line 496 "zgsvj0.f"
					apoaq = aapp / aaqq;
#line 497 "zgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;

#line 499 "zgsvj0.f"
					if (abs(theta) > bigtheta) {

#line 501 "zgsvj0.f"
					    t = .5 / theta;
#line 502 "zgsvj0.f"
					    cs = 1.;
#line 504 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 504 "zgsvj0.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 504 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 506 "zgsvj0.f"
					    if (rsvec) {
#line 507 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 507 "zgsvj0.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 507 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 509 "zgsvj0.f"
					    }
/* Computing MAX */
#line 511 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 511 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 513 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 513 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 515 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 515 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);

#line 517 "zgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 521 "zgsvj0.f"
					    thsign = -d_sign(&c_b27, &aapq1);
#line 522 "zgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 524 "zgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 525 "zgsvj0.f"
					    sn = t * cs;

/* Computing MAX */
#line 527 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 527 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 528 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 528 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 530 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 530 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 533 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 533 "zgsvj0.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 533 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 535 "zgsvj0.f"
					    if (rsvec) {
#line 536 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 536 "zgsvj0.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 536 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 538 "zgsvj0.f"
					    }
#line 539 "zgsvj0.f"
					}
#line 540 "zgsvj0.f"
					i__6 = p;
#line 540 "zgsvj0.f"
					i__7 = q;
#line 540 "zgsvj0.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 540 "zgsvj0.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 540 "zgsvj0.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 542 "zgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 544 "zgsvj0.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 546 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 549 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 551 "zgsvj0.f"
					z__1.r = -aapq.r, z__1.i = -aapq.i;
#line 551 "zgsvj0.f"
					zaxpy_(m, &z__1, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 553 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &c_b27, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 555 "zgsvj0.f"
					d__1 = 0., d__2 = 1. - aapq1 * aapq1;
#line 555 "zgsvj0.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 557 "zgsvj0.f"
					mxsinj = max(mxsinj,*sfmin);
#line 558 "zgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 564 "zgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 564 "zgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 566 "zgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 568 "zgsvj0.f"
					    sva[q] = dznrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 569 "zgsvj0.f"
					} else {
#line 570 "zgsvj0.f"
					    t = 0.;
#line 571 "zgsvj0.f"
					    aaqq = 1.;
#line 572 "zgsvj0.f"
					    zlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 574 "zgsvj0.f"
					    sva[q] = t * sqrt(aaqq);
#line 575 "zgsvj0.f"
					}
#line 576 "zgsvj0.f"
				    }
#line 577 "zgsvj0.f"
				    if (aapp / aapp0 <= rooteps) {
#line 578 "zgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 580 "zgsvj0.f"
					    aapp = dznrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 581 "zgsvj0.f"
					} else {
#line 582 "zgsvj0.f"
					    t = 0.;
#line 583 "zgsvj0.f"
					    aapp = 1.;
#line 584 "zgsvj0.f"
					    zlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 586 "zgsvj0.f"
					    aapp = t * sqrt(aapp);
#line 587 "zgsvj0.f"
					}
#line 588 "zgsvj0.f"
					sva[p] = aapp;
#line 589 "zgsvj0.f"
				    }

#line 591 "zgsvj0.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 593 "zgsvj0.f"
				    if (ir1 == 0) {
#line 593 "zgsvj0.f"
					++notrot;
#line 593 "zgsvj0.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 595 "zgsvj0.f"
				    ++pskipped;
#line 596 "zgsvj0.f"
				}
#line 597 "zgsvj0.f"
			    } else {
/*        A(:,q) is zero column */
#line 599 "zgsvj0.f"
				if (ir1 == 0) {
#line 599 "zgsvj0.f"
				    ++notrot;
#line 599 "zgsvj0.f"
				}
#line 600 "zgsvj0.f"
				++pskipped;
#line 601 "zgsvj0.f"
			    }

#line 603 "zgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 605 "zgsvj0.f"
				if (ir1 == 0) {
#line 605 "zgsvj0.f"
				    aapp = -aapp;
#line 605 "zgsvj0.f"
				}
#line 606 "zgsvj0.f"
				notrot = 0;
#line 607 "zgsvj0.f"
				goto L2103;
#line 608 "zgsvj0.f"
			    }

#line 610 "zgsvj0.f"
/* L2002: */
#line 610 "zgsvj0.f"
			}
/*     END q-LOOP */

#line 613 "zgsvj0.f"
L2103:
/*     bailed out of q-loop */

#line 616 "zgsvj0.f"
			sva[p] = aapp;

#line 618 "zgsvj0.f"
		    } else {
#line 619 "zgsvj0.f"
			sva[p] = aapp;
#line 620 "zgsvj0.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 620 "zgsvj0.f"
			    i__5 = igl + kbl - 1;
#line 620 "zgsvj0.f"
			    notrot = notrot + min(i__5,*n) - p;
#line 620 "zgsvj0.f"
			}
#line 622 "zgsvj0.f"
		    }

#line 624 "zgsvj0.f"
/* L2001: */
#line 624 "zgsvj0.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 627 "zgsvj0.f"
/* L1002: */
#line 627 "zgsvj0.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 632 "zgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

#line 634 "zgsvj0.f"
	    i__3 = nbl;
#line 634 "zgsvj0.f"
	    for (jbc = ibr + 1; jbc <= i__3; ++jbc) {

#line 636 "zgsvj0.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 640 "zgsvj0.f"
		ijblsk = 0;
/* Computing MIN */
#line 641 "zgsvj0.f"
		i__5 = igl + kbl - 1;
#line 641 "zgsvj0.f"
		i__4 = min(i__5,*n);
#line 641 "zgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

#line 643 "zgsvj0.f"
		    aapp = sva[p];
#line 644 "zgsvj0.f"
		    if (aapp > 0.) {

#line 646 "zgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 648 "zgsvj0.f"
			i__6 = jgl + kbl - 1;
#line 648 "zgsvj0.f"
			i__5 = min(i__6,*n);
#line 648 "zgsvj0.f"
			for (q = jgl; q <= i__5; ++q) {

#line 650 "zgsvj0.f"
			    aaqq = sva[q];
#line 651 "zgsvj0.f"
			    if (aaqq > 0.) {
#line 652 "zgsvj0.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 658 "zgsvj0.f"
				if (aaqq >= 1.) {
#line 659 "zgsvj0.f"
				    if (aapp >= aaqq) {
#line 660 "zgsvj0.f"
					rotok = small * aapp <= aaqq;
#line 661 "zgsvj0.f"
				    } else {
#line 662 "zgsvj0.f"
					rotok = small * aaqq <= aapp;
#line 663 "zgsvj0.f"
				    }
#line 664 "zgsvj0.f"
				    if (aapp < big / aaqq) {
#line 665 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 665 "zgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 665 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 665 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 667 "zgsvj0.f"
				    } else {
#line 668 "zgsvj0.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 670 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 673 "zgsvj0.f"
					zdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 673 "zgsvj0.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 673 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 675 "zgsvj0.f"
				    }
#line 676 "zgsvj0.f"
				} else {
#line 677 "zgsvj0.f"
				    if (aapp >= aaqq) {
#line 678 "zgsvj0.f"
					rotok = aapp <= aaqq / small;
#line 679 "zgsvj0.f"
				    } else {
#line 680 "zgsvj0.f"
					rotok = aaqq <= aapp / small;
#line 681 "zgsvj0.f"
				    }
#line 682 "zgsvj0.f"
				    if (aapp > small / aaqq) {
#line 683 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 683 "zgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 683 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 683 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 685 "zgsvj0.f"
				    } else {
#line 686 "zgsvj0.f"
					zcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 688 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 691 "zgsvj0.f"
					zdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 691 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 691 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 693 "zgsvj0.f"
				    }
#line 694 "zgsvj0.f"
				}

#line 696 "zgsvj0.f"
				d__1 = z_abs(&aapq);
#line 696 "zgsvj0.f"
				z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					d__1;
#line 696 "zgsvj0.f"
				ompq.r = z__1.r, ompq.i = z__1.i;
/*                           AAPQ = AAPQ * DCONJG(CWORK(p))*CWORK(q) */
#line 698 "zgsvj0.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 699 "zgsvj0.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 699 "zgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 703 "zgsvj0.f"
				if (abs(aapq1) > *tol) {
#line 704 "zgsvj0.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 706 "zgsvj0.f"
				    pskipped = 0;
#line 707 "zgsvj0.f"
				    ++iswrot;

#line 709 "zgsvj0.f"
				    if (rotok) {

#line 711 "zgsvj0.f"
					aqoap = aaqq / aapp;
#line 712 "zgsvj0.f"
					apoaq = aapp / aaqq;
#line 713 "zgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 714 "zgsvj0.f"
					if (aaqq > aapp0) {
#line 714 "zgsvj0.f"
					    theta = -theta;
#line 714 "zgsvj0.f"
					}

#line 716 "zgsvj0.f"
					if (abs(theta) > bigtheta) {
#line 717 "zgsvj0.f"
					    t = .5 / theta;
#line 718 "zgsvj0.f"
					    cs = 1.;
#line 719 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 719 "zgsvj0.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 719 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 721 "zgsvj0.f"
					    if (rsvec) {
#line 722 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 722 "zgsvj0.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 722 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 724 "zgsvj0.f"
					    }
/* Computing MAX */
#line 725 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 725 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 727 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 727 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 729 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 729 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);
#line 730 "zgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 734 "zgsvj0.f"
					    thsign = -d_sign(&c_b27, &aapq1);
#line 735 "zgsvj0.f"
					    if (aaqq > aapp0) {
#line 735 "zgsvj0.f"
			  thsign = -thsign;
#line 735 "zgsvj0.f"
					    }
#line 736 "zgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 738 "zgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 739 "zgsvj0.f"
					    sn = t * cs;
/* Computing MAX */
#line 740 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 740 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 741 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 741 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 743 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 743 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 746 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 746 "zgsvj0.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 746 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 748 "zgsvj0.f"
					    if (rsvec) {
#line 749 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 749 "zgsvj0.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 749 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 751 "zgsvj0.f"
					    }
#line 752 "zgsvj0.f"
					}
#line 753 "zgsvj0.f"
					i__6 = p;
#line 753 "zgsvj0.f"
					i__7 = q;
#line 753 "zgsvj0.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 753 "zgsvj0.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 753 "zgsvj0.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 755 "zgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 757 "zgsvj0.f"
					if (aapp > aaqq) {
#line 758 "zgsvj0.f"
					    zcopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 760 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b27, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 763 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b27, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 766 "zgsvj0.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 766 "zgsvj0.f"
					    zaxpy_(m, &z__1, &work[1], &c__1, 
						    &a[q * a_dim1 + 1], &c__1)
						    ;
#line 768 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &c_b27,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 771 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 771 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 773 "zgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 774 "zgsvj0.f"
					} else {
#line 775 "zgsvj0.f"
					    zcopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 777 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b27, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 780 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b27, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 783 "zgsvj0.f"
					    d_cnjg(&z__2, &aapq);
#line 783 "zgsvj0.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 783 "zgsvj0.f"
					    zaxpy_(m, &z__1, &work[1], &c__1, 
						    &a[p * a_dim1 + 1], &c__1)
						    ;
#line 785 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &c_b27,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 788 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 788 "zgsvj0.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 790 "zgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 791 "zgsvj0.f"
					}
#line 792 "zgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 797 "zgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 797 "zgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 799 "zgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 801 "zgsvj0.f"
					    sva[q] = dznrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 802 "zgsvj0.f"
					} else {
#line 803 "zgsvj0.f"
					    t = 0.;
#line 804 "zgsvj0.f"
					    aaqq = 1.;
#line 805 "zgsvj0.f"
					    zlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 807 "zgsvj0.f"
					    sva[q] = t * sqrt(aaqq);
#line 808 "zgsvj0.f"
					}
#line 809 "zgsvj0.f"
				    }
/* Computing 2nd power */
#line 810 "zgsvj0.f"
				    d__1 = aapp / aapp0;
#line 810 "zgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 811 "zgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 813 "zgsvj0.f"
					    aapp = dznrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 814 "zgsvj0.f"
					} else {
#line 815 "zgsvj0.f"
					    t = 0.;
#line 816 "zgsvj0.f"
					    aapp = 1.;
#line 817 "zgsvj0.f"
					    zlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 819 "zgsvj0.f"
					    aapp = t * sqrt(aapp);
#line 820 "zgsvj0.f"
					}
#line 821 "zgsvj0.f"
					sva[p] = aapp;
#line 822 "zgsvj0.f"
				    }
/*              end of OK rotation */
#line 824 "zgsvj0.f"
				} else {
#line 825 "zgsvj0.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 827 "zgsvj0.f"
				    ++pskipped;
#line 828 "zgsvj0.f"
				    ++ijblsk;
#line 829 "zgsvj0.f"
				}
#line 830 "zgsvj0.f"
			    } else {
#line 831 "zgsvj0.f"
				++notrot;
#line 832 "zgsvj0.f"
				++pskipped;
#line 833 "zgsvj0.f"
				++ijblsk;
#line 834 "zgsvj0.f"
			    }

#line 836 "zgsvj0.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 838 "zgsvj0.f"
				sva[p] = aapp;
#line 839 "zgsvj0.f"
				notrot = 0;
#line 840 "zgsvj0.f"
				goto L2011;
#line 841 "zgsvj0.f"
			    }
#line 842 "zgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 844 "zgsvj0.f"
				aapp = -aapp;
#line 845 "zgsvj0.f"
				notrot = 0;
#line 846 "zgsvj0.f"
				goto L2203;
#line 847 "zgsvj0.f"
			    }

#line 849 "zgsvj0.f"
/* L2200: */
#line 849 "zgsvj0.f"
			}
/*        end of the q-loop */
#line 851 "zgsvj0.f"
L2203:

#line 853 "zgsvj0.f"
			sva[p] = aapp;

#line 855 "zgsvj0.f"
		    } else {

#line 857 "zgsvj0.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 857 "zgsvj0.f"
			    i__5 = jgl + kbl - 1;
#line 857 "zgsvj0.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 857 "zgsvj0.f"
			}
#line 859 "zgsvj0.f"
			if (aapp < 0.) {
#line 859 "zgsvj0.f"
			    notrot = 0;
#line 859 "zgsvj0.f"
			}

#line 861 "zgsvj0.f"
		    }

#line 863 "zgsvj0.f"
/* L2100: */
#line 863 "zgsvj0.f"
		}
/*     end of the p-loop */
#line 865 "zgsvj0.f"
/* L2010: */
#line 865 "zgsvj0.f"
	    }
/*     end of the jbc-loop */
#line 867 "zgsvj0.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 869 "zgsvj0.f"
	    i__4 = igl + kbl - 1;
#line 869 "zgsvj0.f"
	    i__3 = min(i__4,*n);
#line 869 "zgsvj0.f"
	    for (p = igl; p <= i__3; ++p) {
#line 870 "zgsvj0.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 871 "zgsvj0.f"
/* L2012: */
#line 871 "zgsvj0.f"
	    }
/* ** */
#line 873 "zgsvj0.f"
/* L2000: */
#line 873 "zgsvj0.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 877 "zgsvj0.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 879 "zgsvj0.f"
	    sva[*n] = dznrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 880 "zgsvj0.f"
	} else {
#line 881 "zgsvj0.f"
	    t = 0.;
#line 882 "zgsvj0.f"
	    aapp = 1.;
#line 883 "zgsvj0.f"
	    zlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 884 "zgsvj0.f"
	    sva[*n] = t * sqrt(aapp);
#line 885 "zgsvj0.f"
	}

/*     Additional steering devices */

#line 889 "zgsvj0.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 889 "zgsvj0.f"
	    swband = i__;
#line 889 "zgsvj0.f"
	}

#line 892 "zgsvj0.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 894 "zgsvj0.f"
	    goto L1994;
#line 895 "zgsvj0.f"
	}

#line 897 "zgsvj0.f"
	if (notrot >= emptsw) {
#line 897 "zgsvj0.f"
	    goto L1994;
#line 897 "zgsvj0.f"
	}

#line 899 "zgsvj0.f"
/* L1993: */
#line 899 "zgsvj0.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 903 "zgsvj0.f"
    *info = *nsweep - 1;
#line 904 "zgsvj0.f"
    goto L1995;

#line 906 "zgsvj0.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 910 "zgsvj0.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 912 "zgsvj0.f"
L1995:

/*     Sort the vector SVA() of column norms. */
#line 915 "zgsvj0.f"
    i__1 = *n - 1;
#line 915 "zgsvj0.f"
    for (p = 1; p <= i__1; ++p) {
#line 916 "zgsvj0.f"
	i__2 = *n - p + 1;
#line 916 "zgsvj0.f"
	q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 917 "zgsvj0.f"
	if (p != q) {
#line 918 "zgsvj0.f"
	    temp1 = sva[p];
#line 919 "zgsvj0.f"
	    sva[p] = sva[q];
#line 920 "zgsvj0.f"
	    sva[q] = temp1;
#line 921 "zgsvj0.f"
	    i__2 = p;
#line 921 "zgsvj0.f"
	    aapq.r = d__[i__2].r, aapq.i = d__[i__2].i;
#line 922 "zgsvj0.f"
	    i__2 = p;
#line 922 "zgsvj0.f"
	    i__3 = q;
#line 922 "zgsvj0.f"
	    d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
#line 923 "zgsvj0.f"
	    i__2 = q;
#line 923 "zgsvj0.f"
	    d__[i__2].r = aapq.r, d__[i__2].i = aapq.i;
#line 924 "zgsvj0.f"
	    zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 925 "zgsvj0.f"
	    if (rsvec) {
#line 925 "zgsvj0.f"
		zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 925 "zgsvj0.f"
	    }
#line 926 "zgsvj0.f"
	}
#line 927 "zgsvj0.f"
/* L5991: */
#line 927 "zgsvj0.f"
    }

#line 929 "zgsvj0.f"
    return 0;
/*     .. */
/*     .. END OF ZGSVJ0 */
/*     .. */
} /* zgsvj0_ */

