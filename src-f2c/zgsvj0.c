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

/* > \brief <b> ZGSVJ0 pre-processor for the routine zgesvj. </b> */

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

/* > \date June 2016 */

/* > \ingroup complex16OTHERcomputational */
/* > */
/* > \par Further Details: */
/*  ===================== */
/* > */
/* > ZGSVJ0 is used just to enable ZGESVJ to call a simplified version of */
/* > itself to work on a submatrix of the original matrix. */
/* > */
/* > Contributor: */
/* ============= */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) */
/* > */
/* > \par Bugs, Examples and Comments: */
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     from BLAS */
/*     from LAPACK */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 279 "zgsvj0.f"
    /* Parameter adjustments */
#line 279 "zgsvj0.f"
    --sva;
#line 279 "zgsvj0.f"
    --d__;
#line 279 "zgsvj0.f"
    a_dim1 = *lda;
#line 279 "zgsvj0.f"
    a_offset = 1 + a_dim1;
#line 279 "zgsvj0.f"
    a -= a_offset;
#line 279 "zgsvj0.f"
    v_dim1 = *ldv;
#line 279 "zgsvj0.f"
    v_offset = 1 + v_dim1;
#line 279 "zgsvj0.f"
    v -= v_offset;
#line 279 "zgsvj0.f"
    --work;
#line 279 "zgsvj0.f"

#line 279 "zgsvj0.f"
    /* Function Body */
#line 279 "zgsvj0.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 280 "zgsvj0.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 281 "zgsvj0.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 282 "zgsvj0.f"
	*info = -1;
#line 283 "zgsvj0.f"
    } else if (*m < 0) {
#line 284 "zgsvj0.f"
	*info = -2;
#line 285 "zgsvj0.f"
    } else if (*n < 0 || *n > *m) {
#line 286 "zgsvj0.f"
	*info = -3;
#line 287 "zgsvj0.f"
    } else if (*lda < *m) {
#line 288 "zgsvj0.f"
	*info = -5;
#line 289 "zgsvj0.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 290 "zgsvj0.f"
	*info = -8;
#line 291 "zgsvj0.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 293 "zgsvj0.f"
	*info = -10;
#line 294 "zgsvj0.f"
    } else if (*tol <= *eps) {
#line 295 "zgsvj0.f"
	*info = -13;
#line 296 "zgsvj0.f"
    } else if (*nsweep < 0) {
#line 297 "zgsvj0.f"
	*info = -14;
#line 298 "zgsvj0.f"
    } else if (*lwork < *m) {
#line 299 "zgsvj0.f"
	*info = -16;
#line 300 "zgsvj0.f"
    } else {
#line 301 "zgsvj0.f"
	*info = 0;
#line 302 "zgsvj0.f"
    }

/*     #:( */
#line 305 "zgsvj0.f"
    if (*info != 0) {
#line 306 "zgsvj0.f"
	i__1 = -(*info);
#line 306 "zgsvj0.f"
	xerbla_("ZGSVJ0", &i__1, (ftnlen)6);
#line 307 "zgsvj0.f"
	return 0;
#line 308 "zgsvj0.f"
    }

#line 310 "zgsvj0.f"
    if (rsvec) {
#line 311 "zgsvj0.f"
	mvl = *n;
#line 312 "zgsvj0.f"
    } else if (applv) {
#line 313 "zgsvj0.f"
	mvl = *mv;
#line 314 "zgsvj0.f"
    }
#line 315 "zgsvj0.f"
    rsvec = rsvec || applv;
#line 317 "zgsvj0.f"
    rooteps = sqrt(*eps);
#line 318 "zgsvj0.f"
    rootsfmin = sqrt(*sfmin);
#line 319 "zgsvj0.f"
    small = *sfmin / *eps;
#line 320 "zgsvj0.f"
    big = 1. / *sfmin;
#line 321 "zgsvj0.f"
    rootbig = 1. / rootsfmin;
#line 322 "zgsvj0.f"
    bigtheta = 1. / rooteps;
#line 323 "zgsvj0.f"
    roottol = sqrt(*tol);

/*     .. Row-cyclic Jacobi SVD algorithm with column pivoting .. */

#line 327 "zgsvj0.f"
    emptsw = *n * (*n - 1) / 2;
#line 328 "zgsvj0.f"
    notrot = 0;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 333 "zgsvj0.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if ZGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 341 "zgsvj0.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 347 "zgsvj0.f"
    nbl = *n / kbl;
#line 348 "zgsvj0.f"
    if (nbl * kbl != *n) {
#line 348 "zgsvj0.f"
	++nbl;
#line 348 "zgsvj0.f"
    }

/* Computing 2nd power */
#line 350 "zgsvj0.f"
    i__1 = kbl;
#line 350 "zgsvj0.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 353 "zgsvj0.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 356 "zgsvj0.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */


/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 367 "zgsvj0.f"
    i__1 = *nsweep;
#line 367 "zgsvj0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     .. go go go ... */

#line 371 "zgsvj0.f"
	mxaapq = 0.;
#line 372 "zgsvj0.f"
	mxsinj = 0.;
#line 373 "zgsvj0.f"
	iswrot = 0;

#line 375 "zgsvj0.f"
	notrot = 0;
#line 376 "zgsvj0.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 383 "zgsvj0.f"
	i__2 = nbl;
#line 383 "zgsvj0.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {

#line 385 "zgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 387 "zgsvj0.f"
	    i__4 = lkahead, i__5 = nbl - ibr;
#line 387 "zgsvj0.f"
	    i__3 = min(i__4,i__5);
#line 387 "zgsvj0.f"
	    for (ir1 = 0; ir1 <= i__3; ++ir1) {

#line 389 "zgsvj0.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 391 "zgsvj0.f"
		i__5 = igl + kbl - 1, i__6 = *n - 1;
#line 391 "zgsvj0.f"
		i__4 = min(i__5,i__6);
#line 391 "zgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

/*     .. de Rijk's pivoting */

#line 395 "zgsvj0.f"
		    i__5 = *n - p + 1;
#line 395 "zgsvj0.f"
		    q = idamax_(&i__5, &sva[p], &c__1) + p - 1;
#line 396 "zgsvj0.f"
		    if (p != q) {
#line 397 "zgsvj0.f"
			zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 398 "zgsvj0.f"
			if (rsvec) {
#line 398 "zgsvj0.f"
			    zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 398 "zgsvj0.f"
			}
#line 400 "zgsvj0.f"
			temp1 = sva[p];
#line 401 "zgsvj0.f"
			sva[p] = sva[q];
#line 402 "zgsvj0.f"
			sva[q] = temp1;
#line 403 "zgsvj0.f"
			i__5 = p;
#line 403 "zgsvj0.f"
			aapq.r = d__[i__5].r, aapq.i = d__[i__5].i;
#line 404 "zgsvj0.f"
			i__5 = p;
#line 404 "zgsvj0.f"
			i__6 = q;
#line 404 "zgsvj0.f"
			d__[i__5].r = d__[i__6].r, d__[i__5].i = d__[i__6].i;
#line 405 "zgsvj0.f"
			i__5 = q;
#line 405 "zgsvj0.f"
			d__[i__5].r = aapq.r, d__[i__5].i = aapq.i;
#line 406 "zgsvj0.f"
		    }

#line 408 "zgsvj0.f"
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

#line 422 "zgsvj0.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 424 "zgsvj0.f"
			    sva[p] = dznrm2_(m, &a[p * a_dim1 + 1], &c__1);
#line 425 "zgsvj0.f"
			} else {
#line 426 "zgsvj0.f"
			    temp1 = 0.;
#line 427 "zgsvj0.f"
			    aapp = 1.;
#line 428 "zgsvj0.f"
			    zlassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 429 "zgsvj0.f"
			    sva[p] = temp1 * sqrt(aapp);
#line 430 "zgsvj0.f"
			}
#line 431 "zgsvj0.f"
			aapp = sva[p];
#line 432 "zgsvj0.f"
		    } else {
#line 433 "zgsvj0.f"
			aapp = sva[p];
#line 434 "zgsvj0.f"
		    }

#line 436 "zgsvj0.f"
		    if (aapp > 0.) {

#line 438 "zgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 440 "zgsvj0.f"
			i__6 = igl + kbl - 1;
#line 440 "zgsvj0.f"
			i__5 = min(i__6,*n);
#line 440 "zgsvj0.f"
			for (q = p + 1; q <= i__5; ++q) {

#line 442 "zgsvj0.f"
			    aaqq = sva[q];

#line 444 "zgsvj0.f"
			    if (aaqq > 0.) {

#line 446 "zgsvj0.f"
				aapp0 = aapp;
#line 447 "zgsvj0.f"
				if (aaqq >= 1.) {
#line 448 "zgsvj0.f"
				    rotok = small * aapp <= aaqq;
#line 449 "zgsvj0.f"
				    if (aapp < big / aaqq) {
#line 450 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 450 "zgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 450 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 450 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 452 "zgsvj0.f"
				    } else {
#line 453 "zgsvj0.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 455 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 457 "zgsvj0.f"
					zdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 457 "zgsvj0.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 457 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 459 "zgsvj0.f"
				    }
#line 460 "zgsvj0.f"
				} else {
#line 461 "zgsvj0.f"
				    rotok = aapp <= aaqq / small;
#line 462 "zgsvj0.f"
				    if (aapp > small / aaqq) {
#line 463 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 463 "zgsvj0.f"
					z__2.r = z__3.r / aapp, z__2.i = 
						z__3.i / aapp;
#line 463 "zgsvj0.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 463 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 465 "zgsvj0.f"
				    } else {
#line 466 "zgsvj0.f"
					zcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 468 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 471 "zgsvj0.f"
					zdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 471 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 471 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 473 "zgsvj0.f"
				    }
#line 474 "zgsvj0.f"
				}

/*                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q) */
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
#line 483 "zgsvj0.f"
				    d__1 = z_abs(&aapq);
#line 483 "zgsvj0.f"
				    z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					    d__1;
#line 483 "zgsvj0.f"
				    ompq.r = z__1.r, ompq.i = z__1.i;

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 488 "zgsvj0.f"
				    if (ir1 == 0) {
#line 489 "zgsvj0.f"
					notrot = 0;
#line 490 "zgsvj0.f"
					pskipped = 0;
#line 491 "zgsvj0.f"
					++iswrot;
#line 492 "zgsvj0.f"
				    }

#line 494 "zgsvj0.f"
				    if (rotok) {

#line 496 "zgsvj0.f"
					aqoap = aaqq / aapp;
#line 497 "zgsvj0.f"
					apoaq = aapp / aaqq;
#line 498 "zgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;

#line 500 "zgsvj0.f"
					if (abs(theta) > bigtheta) {

#line 502 "zgsvj0.f"
					    t = .5 / theta;
#line 503 "zgsvj0.f"
					    cs = 1.;
#line 505 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 505 "zgsvj0.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 505 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 507 "zgsvj0.f"
					    if (rsvec) {
#line 508 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 508 "zgsvj0.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 508 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 510 "zgsvj0.f"
					    }
/* Computing MAX */
#line 512 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 512 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 514 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 514 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 516 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 516 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);

#line 518 "zgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 522 "zgsvj0.f"
					    thsign = -d_sign(&c_b27, &aapq1);
#line 523 "zgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 525 "zgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 526 "zgsvj0.f"
					    sn = t * cs;

/* Computing MAX */
#line 528 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 528 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 529 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 529 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 531 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 531 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 534 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 534 "zgsvj0.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 534 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 536 "zgsvj0.f"
					    if (rsvec) {
#line 537 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 537 "zgsvj0.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 537 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 539 "zgsvj0.f"
					    }
#line 540 "zgsvj0.f"
					}
#line 541 "zgsvj0.f"
					i__6 = p;
#line 541 "zgsvj0.f"
					i__7 = q;
#line 541 "zgsvj0.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 541 "zgsvj0.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 541 "zgsvj0.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 543 "zgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 545 "zgsvj0.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 547 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 550 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 552 "zgsvj0.f"
					z__1.r = -aapq.r, z__1.i = -aapq.i;
#line 552 "zgsvj0.f"
					zaxpy_(m, &z__1, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 554 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &c_b27, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 556 "zgsvj0.f"
					d__1 = 0., d__2 = 1. - aapq1 * aapq1;
#line 556 "zgsvj0.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 558 "zgsvj0.f"
					mxsinj = max(mxsinj,*sfmin);
#line 559 "zgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 565 "zgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 565 "zgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 567 "zgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 569 "zgsvj0.f"
					    sva[q] = dznrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 570 "zgsvj0.f"
					} else {
#line 571 "zgsvj0.f"
					    t = 0.;
#line 572 "zgsvj0.f"
					    aaqq = 1.;
#line 573 "zgsvj0.f"
					    zlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 575 "zgsvj0.f"
					    sva[q] = t * sqrt(aaqq);
#line 576 "zgsvj0.f"
					}
#line 577 "zgsvj0.f"
				    }
#line 578 "zgsvj0.f"
				    if (aapp / aapp0 <= rooteps) {
#line 579 "zgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 581 "zgsvj0.f"
					    aapp = dznrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 582 "zgsvj0.f"
					} else {
#line 583 "zgsvj0.f"
					    t = 0.;
#line 584 "zgsvj0.f"
					    aapp = 1.;
#line 585 "zgsvj0.f"
					    zlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 587 "zgsvj0.f"
					    aapp = t * sqrt(aapp);
#line 588 "zgsvj0.f"
					}
#line 589 "zgsvj0.f"
					sva[p] = aapp;
#line 590 "zgsvj0.f"
				    }

#line 592 "zgsvj0.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 594 "zgsvj0.f"
				    if (ir1 == 0) {
#line 594 "zgsvj0.f"
					++notrot;
#line 594 "zgsvj0.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 596 "zgsvj0.f"
				    ++pskipped;
#line 597 "zgsvj0.f"
				}
#line 598 "zgsvj0.f"
			    } else {
/*        A(:,q) is zero column */
#line 600 "zgsvj0.f"
				if (ir1 == 0) {
#line 600 "zgsvj0.f"
				    ++notrot;
#line 600 "zgsvj0.f"
				}
#line 601 "zgsvj0.f"
				++pskipped;
#line 602 "zgsvj0.f"
			    }

#line 604 "zgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 606 "zgsvj0.f"
				if (ir1 == 0) {
#line 606 "zgsvj0.f"
				    aapp = -aapp;
#line 606 "zgsvj0.f"
				}
#line 607 "zgsvj0.f"
				notrot = 0;
#line 608 "zgsvj0.f"
				goto L2103;
#line 609 "zgsvj0.f"
			    }

#line 611 "zgsvj0.f"
/* L2002: */
#line 611 "zgsvj0.f"
			}
/*     END q-LOOP */

#line 614 "zgsvj0.f"
L2103:
/*     bailed out of q-loop */

#line 617 "zgsvj0.f"
			sva[p] = aapp;

#line 619 "zgsvj0.f"
		    } else {
#line 620 "zgsvj0.f"
			sva[p] = aapp;
#line 621 "zgsvj0.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 621 "zgsvj0.f"
			    i__5 = igl + kbl - 1;
#line 621 "zgsvj0.f"
			    notrot = notrot + min(i__5,*n) - p;
#line 621 "zgsvj0.f"
			}
#line 623 "zgsvj0.f"
		    }

#line 625 "zgsvj0.f"
/* L2001: */
#line 625 "zgsvj0.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 628 "zgsvj0.f"
/* L1002: */
#line 628 "zgsvj0.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 633 "zgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

#line 635 "zgsvj0.f"
	    i__3 = nbl;
#line 635 "zgsvj0.f"
	    for (jbc = ibr + 1; jbc <= i__3; ++jbc) {

#line 637 "zgsvj0.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 641 "zgsvj0.f"
		ijblsk = 0;
/* Computing MIN */
#line 642 "zgsvj0.f"
		i__5 = igl + kbl - 1;
#line 642 "zgsvj0.f"
		i__4 = min(i__5,*n);
#line 642 "zgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

#line 644 "zgsvj0.f"
		    aapp = sva[p];
#line 645 "zgsvj0.f"
		    if (aapp > 0.) {

#line 647 "zgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 649 "zgsvj0.f"
			i__6 = jgl + kbl - 1;
#line 649 "zgsvj0.f"
			i__5 = min(i__6,*n);
#line 649 "zgsvj0.f"
			for (q = jgl; q <= i__5; ++q) {

#line 651 "zgsvj0.f"
			    aaqq = sva[q];
#line 652 "zgsvj0.f"
			    if (aaqq > 0.) {
#line 653 "zgsvj0.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 659 "zgsvj0.f"
				if (aaqq >= 1.) {
#line 660 "zgsvj0.f"
				    if (aapp >= aaqq) {
#line 661 "zgsvj0.f"
					rotok = small * aapp <= aaqq;
#line 662 "zgsvj0.f"
				    } else {
#line 663 "zgsvj0.f"
					rotok = small * aaqq <= aapp;
#line 664 "zgsvj0.f"
				    }
#line 665 "zgsvj0.f"
				    if (aapp < big / aaqq) {
#line 666 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 666 "zgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 666 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 666 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 668 "zgsvj0.f"
				    } else {
#line 669 "zgsvj0.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 671 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 674 "zgsvj0.f"
					zdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 674 "zgsvj0.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 674 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 676 "zgsvj0.f"
				    }
#line 677 "zgsvj0.f"
				} else {
#line 678 "zgsvj0.f"
				    if (aapp >= aaqq) {
#line 679 "zgsvj0.f"
					rotok = aapp <= aaqq / small;
#line 680 "zgsvj0.f"
				    } else {
#line 681 "zgsvj0.f"
					rotok = aaqq <= aapp / small;
#line 682 "zgsvj0.f"
				    }
#line 683 "zgsvj0.f"
				    if (aapp > small / aaqq) {
#line 684 "zgsvj0.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 684 "zgsvj0.f"
					d__1 = max(aaqq,aapp);
#line 684 "zgsvj0.f"
					z__2.r = z__3.r / d__1, z__2.i = 
						z__3.i / d__1;
#line 684 "zgsvj0.f"
					d__2 = min(aaqq,aapp);
#line 684 "zgsvj0.f"
					z__1.r = z__2.r / d__2, z__1.i = 
						z__2.i / d__2;
#line 684 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 687 "zgsvj0.f"
				    } else {
#line 688 "zgsvj0.f"
					zcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 690 "zgsvj0.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 693 "zgsvj0.f"
					zdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 693 "zgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 693 "zgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 695 "zgsvj0.f"
				    }
#line 696 "zgsvj0.f"
				}

/*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
#line 699 "zgsvj0.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 700 "zgsvj0.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 700 "zgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 704 "zgsvj0.f"
				if (abs(aapq1) > *tol) {
#line 705 "zgsvj0.f"
				    d__1 = z_abs(&aapq);
#line 705 "zgsvj0.f"
				    z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					    d__1;
#line 705 "zgsvj0.f"
				    ompq.r = z__1.r, ompq.i = z__1.i;
#line 706 "zgsvj0.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 708 "zgsvj0.f"
				    pskipped = 0;
#line 709 "zgsvj0.f"
				    ++iswrot;

#line 711 "zgsvj0.f"
				    if (rotok) {

#line 713 "zgsvj0.f"
					aqoap = aaqq / aapp;
#line 714 "zgsvj0.f"
					apoaq = aapp / aaqq;
#line 715 "zgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 716 "zgsvj0.f"
					if (aaqq > aapp0) {
#line 716 "zgsvj0.f"
					    theta = -theta;
#line 716 "zgsvj0.f"
					}

#line 718 "zgsvj0.f"
					if (abs(theta) > bigtheta) {
#line 719 "zgsvj0.f"
					    t = .5 / theta;
#line 720 "zgsvj0.f"
					    cs = 1.;
#line 721 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 721 "zgsvj0.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 721 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 723 "zgsvj0.f"
					    if (rsvec) {
#line 724 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 724 "zgsvj0.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 724 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 726 "zgsvj0.f"
					    }
/* Computing MAX */
#line 727 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 727 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 729 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 729 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 731 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 731 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);
#line 732 "zgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 736 "zgsvj0.f"
					    thsign = -d_sign(&c_b27, &aapq1);
#line 737 "zgsvj0.f"
					    if (aaqq > aapp0) {
#line 737 "zgsvj0.f"
			  thsign = -thsign;
#line 737 "zgsvj0.f"
					    }
#line 738 "zgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 740 "zgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 741 "zgsvj0.f"
					    sn = t * cs;
/* Computing MAX */
#line 742 "zgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 742 "zgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 743 "zgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 743 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 745 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 745 "zgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 748 "zgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 748 "zgsvj0.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 748 "zgsvj0.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 750 "zgsvj0.f"
					    if (rsvec) {
#line 751 "zgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 751 "zgsvj0.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 751 "zgsvj0.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 753 "zgsvj0.f"
					    }
#line 754 "zgsvj0.f"
					}
#line 755 "zgsvj0.f"
					i__6 = p;
#line 755 "zgsvj0.f"
					i__7 = q;
#line 755 "zgsvj0.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 755 "zgsvj0.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 755 "zgsvj0.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 757 "zgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 759 "zgsvj0.f"
					if (aapp > aaqq) {
#line 760 "zgsvj0.f"
					    zcopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 762 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b27, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 765 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b27, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 768 "zgsvj0.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 768 "zgsvj0.f"
					    zaxpy_(m, &z__1, &work[1], &c__1, 
						    &a[q * a_dim1 + 1], &c__1)
						    ;
#line 770 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &c_b27,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 773 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 773 "zgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 775 "zgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 776 "zgsvj0.f"
					} else {
#line 777 "zgsvj0.f"
					    zcopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 779 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b27, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 782 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b27, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 785 "zgsvj0.f"
					    d_cnjg(&z__2, &aapq);
#line 785 "zgsvj0.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 785 "zgsvj0.f"
					    zaxpy_(m, &z__1, &work[1], &c__1, 
						    &a[p * a_dim1 + 1], &c__1)
						    ;
#line 787 "zgsvj0.f"
					    zlascl_("G", &c__0, &c__0, &c_b27,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 790 "zgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 790 "zgsvj0.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 792 "zgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 793 "zgsvj0.f"
					}
#line 794 "zgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 799 "zgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 799 "zgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 801 "zgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 803 "zgsvj0.f"
					    sva[q] = dznrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 804 "zgsvj0.f"
					} else {
#line 805 "zgsvj0.f"
					    t = 0.;
#line 806 "zgsvj0.f"
					    aaqq = 1.;
#line 807 "zgsvj0.f"
					    zlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 809 "zgsvj0.f"
					    sva[q] = t * sqrt(aaqq);
#line 810 "zgsvj0.f"
					}
#line 811 "zgsvj0.f"
				    }
/* Computing 2nd power */
#line 812 "zgsvj0.f"
				    d__1 = aapp / aapp0;
#line 812 "zgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 813 "zgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 815 "zgsvj0.f"
					    aapp = dznrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 816 "zgsvj0.f"
					} else {
#line 817 "zgsvj0.f"
					    t = 0.;
#line 818 "zgsvj0.f"
					    aapp = 1.;
#line 819 "zgsvj0.f"
					    zlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 821 "zgsvj0.f"
					    aapp = t * sqrt(aapp);
#line 822 "zgsvj0.f"
					}
#line 823 "zgsvj0.f"
					sva[p] = aapp;
#line 824 "zgsvj0.f"
				    }
/*              end of OK rotation */
#line 826 "zgsvj0.f"
				} else {
#line 827 "zgsvj0.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 829 "zgsvj0.f"
				    ++pskipped;
#line 830 "zgsvj0.f"
				    ++ijblsk;
#line 831 "zgsvj0.f"
				}
#line 832 "zgsvj0.f"
			    } else {
#line 833 "zgsvj0.f"
				++notrot;
#line 834 "zgsvj0.f"
				++pskipped;
#line 835 "zgsvj0.f"
				++ijblsk;
#line 836 "zgsvj0.f"
			    }

#line 838 "zgsvj0.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 840 "zgsvj0.f"
				sva[p] = aapp;
#line 841 "zgsvj0.f"
				notrot = 0;
#line 842 "zgsvj0.f"
				goto L2011;
#line 843 "zgsvj0.f"
			    }
#line 844 "zgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 846 "zgsvj0.f"
				aapp = -aapp;
#line 847 "zgsvj0.f"
				notrot = 0;
#line 848 "zgsvj0.f"
				goto L2203;
#line 849 "zgsvj0.f"
			    }

#line 851 "zgsvj0.f"
/* L2200: */
#line 851 "zgsvj0.f"
			}
/*        end of the q-loop */
#line 853 "zgsvj0.f"
L2203:

#line 855 "zgsvj0.f"
			sva[p] = aapp;

#line 857 "zgsvj0.f"
		    } else {

#line 859 "zgsvj0.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 859 "zgsvj0.f"
			    i__5 = jgl + kbl - 1;
#line 859 "zgsvj0.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 859 "zgsvj0.f"
			}
#line 861 "zgsvj0.f"
			if (aapp < 0.) {
#line 861 "zgsvj0.f"
			    notrot = 0;
#line 861 "zgsvj0.f"
			}

#line 863 "zgsvj0.f"
		    }

#line 865 "zgsvj0.f"
/* L2100: */
#line 865 "zgsvj0.f"
		}
/*     end of the p-loop */
#line 867 "zgsvj0.f"
/* L2010: */
#line 867 "zgsvj0.f"
	    }
/*     end of the jbc-loop */
#line 869 "zgsvj0.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 871 "zgsvj0.f"
	    i__4 = igl + kbl - 1;
#line 871 "zgsvj0.f"
	    i__3 = min(i__4,*n);
#line 871 "zgsvj0.f"
	    for (p = igl; p <= i__3; ++p) {
#line 872 "zgsvj0.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 873 "zgsvj0.f"
/* L2012: */
#line 873 "zgsvj0.f"
	    }
/* ** */
#line 875 "zgsvj0.f"
/* L2000: */
#line 875 "zgsvj0.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 879 "zgsvj0.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 881 "zgsvj0.f"
	    sva[*n] = dznrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 882 "zgsvj0.f"
	} else {
#line 883 "zgsvj0.f"
	    t = 0.;
#line 884 "zgsvj0.f"
	    aapp = 1.;
#line 885 "zgsvj0.f"
	    zlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 886 "zgsvj0.f"
	    sva[*n] = t * sqrt(aapp);
#line 887 "zgsvj0.f"
	}

/*     Additional steering devices */

#line 891 "zgsvj0.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 891 "zgsvj0.f"
	    swband = i__;
#line 891 "zgsvj0.f"
	}

#line 894 "zgsvj0.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 896 "zgsvj0.f"
	    goto L1994;
#line 897 "zgsvj0.f"
	}

#line 899 "zgsvj0.f"
	if (notrot >= emptsw) {
#line 899 "zgsvj0.f"
	    goto L1994;
#line 899 "zgsvj0.f"
	}

#line 901 "zgsvj0.f"
/* L1993: */
#line 901 "zgsvj0.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 905 "zgsvj0.f"
    *info = *nsweep - 1;
#line 906 "zgsvj0.f"
    goto L1995;

#line 908 "zgsvj0.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 912 "zgsvj0.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 914 "zgsvj0.f"
L1995:

/*     Sort the vector SVA() of column norms. */
#line 917 "zgsvj0.f"
    i__1 = *n - 1;
#line 917 "zgsvj0.f"
    for (p = 1; p <= i__1; ++p) {
#line 918 "zgsvj0.f"
	i__2 = *n - p + 1;
#line 918 "zgsvj0.f"
	q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 919 "zgsvj0.f"
	if (p != q) {
#line 920 "zgsvj0.f"
	    temp1 = sva[p];
#line 921 "zgsvj0.f"
	    sva[p] = sva[q];
#line 922 "zgsvj0.f"
	    sva[q] = temp1;
#line 923 "zgsvj0.f"
	    i__2 = p;
#line 923 "zgsvj0.f"
	    aapq.r = d__[i__2].r, aapq.i = d__[i__2].i;
#line 924 "zgsvj0.f"
	    i__2 = p;
#line 924 "zgsvj0.f"
	    i__3 = q;
#line 924 "zgsvj0.f"
	    d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
#line 925 "zgsvj0.f"
	    i__2 = q;
#line 925 "zgsvj0.f"
	    d__[i__2].r = aapq.r, d__[i__2].i = aapq.i;
#line 926 "zgsvj0.f"
	    zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 927 "zgsvj0.f"
	    if (rsvec) {
#line 927 "zgsvj0.f"
		zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 927 "zgsvj0.f"
	    }
#line 928 "zgsvj0.f"
	}
#line 929 "zgsvj0.f"
/* L5991: */
#line 929 "zgsvj0.f"
    }

#line 931 "zgsvj0.f"
    return 0;
/*     .. */
/*     .. END OF ZGSVJ0 */
/*     .. */
} /* zgsvj0_ */

