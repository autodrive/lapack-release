#line 1 "dgsvj0.f"
/* dgsvj0.f -- translated by f2c (version 20100827).
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

#line 1 "dgsvj0.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b42 = 1.;

/* > \brief \b DGSVJ0 pre-processor for the routine dgesvj. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGSVJ0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgsvj0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgsvj0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgsvj0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, */
/*                          SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP */
/*       DOUBLE PRECISION   EPS, SFMIN, TOL */
/*       CHARACTER*1        JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), SVA( N ), D( N ), V( LDV, * ), */
/*      $                   WORK( LWORK ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGSVJ0 is called from DGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as DGESVJ does, but */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, M-by-N matrix A, such that A*diag(D) represents */
/* >          the input matrix. */
/* >          On exit, */
/* >          A_onexit * D_onexit represents the input matrix A*diag(D) */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The array D accumulates the scaling factors from the fast scaled */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > DGSVJ0 is used just to enable DGESVJ to call a simplified version of */
/* > itself to work on a submatrix of the original matrix. */
/* > */
/* > \par Contributors: */
/*  ================== */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */
/* > */
/* > \par Bugs, Examples and Comments: */
/*  ================================= */
/* > */
/* > Please report all bugs and send interesting test examples and comments to */
/* > drmac@math.hr. Thank you. */

/*  ===================================================================== */
/* Subroutine */ int dgsvj0_(char *jobv, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *d__, doublereal *sva, integer *mv, 
	doublereal *v, integer *ldv, doublereal *eps, doublereal *sfmin, 
	doublereal *tol, integer *nsweep, doublereal *work, integer *lwork, 
	integer *info, ftnlen jobv_len)
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
    static integer ir1, jbc;
    static doublereal big;
    static integer kbl, igl, ibr, jgl, nbl, mvl;
    static doublereal aapp, aapq, aaqq;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal aapp0;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp1, apoaq, aqoap;
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
    static integer emptsw, notrot, iswrot, lkahead;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 271 "dgsvj0.f"
    /* Parameter adjustments */
#line 271 "dgsvj0.f"
    --sva;
#line 271 "dgsvj0.f"
    --d__;
#line 271 "dgsvj0.f"
    a_dim1 = *lda;
#line 271 "dgsvj0.f"
    a_offset = 1 + a_dim1;
#line 271 "dgsvj0.f"
    a -= a_offset;
#line 271 "dgsvj0.f"
    v_dim1 = *ldv;
#line 271 "dgsvj0.f"
    v_offset = 1 + v_dim1;
#line 271 "dgsvj0.f"
    v -= v_offset;
#line 271 "dgsvj0.f"
    --work;
#line 271 "dgsvj0.f"

#line 271 "dgsvj0.f"
    /* Function Body */
#line 271 "dgsvj0.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 272 "dgsvj0.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 273 "dgsvj0.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 274 "dgsvj0.f"
	*info = -1;
#line 275 "dgsvj0.f"
    } else if (*m < 0) {
#line 276 "dgsvj0.f"
	*info = -2;
#line 277 "dgsvj0.f"
    } else if (*n < 0 || *n > *m) {
#line 278 "dgsvj0.f"
	*info = -3;
#line 279 "dgsvj0.f"
    } else if (*lda < *m) {
#line 280 "dgsvj0.f"
	*info = -5;
#line 281 "dgsvj0.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 282 "dgsvj0.f"
	*info = -8;
#line 283 "dgsvj0.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 285 "dgsvj0.f"
	*info = -10;
#line 286 "dgsvj0.f"
    } else if (*tol <= *eps) {
#line 287 "dgsvj0.f"
	*info = -13;
#line 288 "dgsvj0.f"
    } else if (*nsweep < 0) {
#line 289 "dgsvj0.f"
	*info = -14;
#line 290 "dgsvj0.f"
    } else if (*lwork < *m) {
#line 291 "dgsvj0.f"
	*info = -16;
#line 292 "dgsvj0.f"
    } else {
#line 293 "dgsvj0.f"
	*info = 0;
#line 294 "dgsvj0.f"
    }

/*     #:( */
#line 297 "dgsvj0.f"
    if (*info != 0) {
#line 298 "dgsvj0.f"
	i__1 = -(*info);
#line 298 "dgsvj0.f"
	xerbla_("DGSVJ0", &i__1, (ftnlen)6);
#line 299 "dgsvj0.f"
	return 0;
#line 300 "dgsvj0.f"
    }

#line 302 "dgsvj0.f"
    if (rsvec) {
#line 303 "dgsvj0.f"
	mvl = *n;
#line 304 "dgsvj0.f"
    } else if (applv) {
#line 305 "dgsvj0.f"
	mvl = *mv;
#line 306 "dgsvj0.f"
    }
#line 307 "dgsvj0.f"
    rsvec = rsvec || applv;
#line 309 "dgsvj0.f"
    rooteps = sqrt(*eps);
#line 310 "dgsvj0.f"
    rootsfmin = sqrt(*sfmin);
#line 311 "dgsvj0.f"
    small = *sfmin / *eps;
#line 312 "dgsvj0.f"
    big = 1. / *sfmin;
#line 313 "dgsvj0.f"
    rootbig = 1. / rootsfmin;
#line 314 "dgsvj0.f"
    bigtheta = 1. / rooteps;
#line 315 "dgsvj0.f"
    roottol = sqrt(*tol);

/*     -#- Row-cyclic Jacobi SVD algorithm with column pivoting -#- */

#line 319 "dgsvj0.f"
    emptsw = *n * (*n - 1) / 2;
#line 320 "dgsvj0.f"
    notrot = 0;
#line 321 "dgsvj0.f"
    fastr[0] = 0.;

/*     -#- Row-cyclic pivot strategy with de Rijk's pivoting -#- */

#line 326 "dgsvj0.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter. It is meaningful and effective */
/*     if SGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure */
/*     ...... */
#line 332 "dgsvj0.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 338 "dgsvj0.f"
    nbl = *n / kbl;
#line 339 "dgsvj0.f"
    if (nbl * kbl != *n) {
#line 339 "dgsvj0.f"
	++nbl;
#line 339 "dgsvj0.f"
    }
/* Computing 2nd power */
#line 341 "dgsvj0.f"
    i__1 = kbl;
#line 341 "dgsvj0.f"
    blskip = i__1 * i__1 + 1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
#line 344 "dgsvj0.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */
#line 347 "dgsvj0.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */
#line 349 "dgsvj0.f"
    swband = 0;
#line 350 "dgsvj0.f"
    pskipped = 0;

#line 352 "dgsvj0.f"
    i__1 = *nsweep;
#line 352 "dgsvj0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     .. go go go ... */

#line 355 "dgsvj0.f"
	mxaapq = 0.;
#line 356 "dgsvj0.f"
	mxsinj = 0.;
#line 357 "dgsvj0.f"
	iswrot = 0;

#line 359 "dgsvj0.f"
	notrot = 0;
#line 360 "dgsvj0.f"
	pskipped = 0;

#line 362 "dgsvj0.f"
	i__2 = nbl;
#line 362 "dgsvj0.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {
#line 364 "dgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 366 "dgsvj0.f"
	    i__4 = lkahead, i__5 = nbl - ibr;
#line 366 "dgsvj0.f"
	    i__3 = min(i__4,i__5);
#line 366 "dgsvj0.f"
	    for (ir1 = 0; ir1 <= i__3; ++ir1) {

#line 368 "dgsvj0.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 370 "dgsvj0.f"
		i__5 = igl + kbl - 1, i__6 = *n - 1;
#line 370 "dgsvj0.f"
		i__4 = min(i__5,i__6);
#line 370 "dgsvj0.f"
		for (p = igl; p <= i__4; ++p) {
/*     .. de Rijk's pivoting */
#line 373 "dgsvj0.f"
		    i__5 = *n - p + 1;
#line 373 "dgsvj0.f"
		    q = idamax_(&i__5, &sva[p], &c__1) + p - 1;
#line 374 "dgsvj0.f"
		    if (p != q) {
#line 375 "dgsvj0.f"
			dswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 376 "dgsvj0.f"
			if (rsvec) {
#line 376 "dgsvj0.f"
			    dswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 376 "dgsvj0.f"
			}
#line 378 "dgsvj0.f"
			temp1 = sva[p];
#line 379 "dgsvj0.f"
			sva[p] = sva[q];
#line 380 "dgsvj0.f"
			sva[q] = temp1;
#line 381 "dgsvj0.f"
			temp1 = d__[p];
#line 382 "dgsvj0.f"
			d__[p] = d__[q];
#line 383 "dgsvj0.f"
			d__[q] = temp1;
#line 384 "dgsvj0.f"
		    }

#line 386 "dgsvj0.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/*        Caveat: */
/*        Some BLAS implementations compute DNRM2(M,A(1,p),1) */
/*        as DSQRT(DDOT(M,A(1,p),1,A(1,p),1)), which may result in */
/*        overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and */
/*        undeflow for ||A(:,p)||_2 < DSQRT(underflow_threshold). */
/*        Hence, DNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented DNRM2 is available, the IF-THEN-ELSE */
/*        below should read "AAPP = DNRM2( M, A(1,p), 1 ) * D(p)". */

#line 400 "dgsvj0.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 402 "dgsvj0.f"
			    sva[p] = dnrm2_(m, &a[p * a_dim1 + 1], &c__1) * 
				    d__[p];
#line 403 "dgsvj0.f"
			} else {
#line 404 "dgsvj0.f"
			    temp1 = 0.;
#line 405 "dgsvj0.f"
			    aapp = 1.;
#line 406 "dgsvj0.f"
			    dlassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 407 "dgsvj0.f"
			    sva[p] = temp1 * sqrt(aapp) * d__[p];
#line 408 "dgsvj0.f"
			}
#line 409 "dgsvj0.f"
			aapp = sva[p];
#line 410 "dgsvj0.f"
		    } else {
#line 411 "dgsvj0.f"
			aapp = sva[p];
#line 412 "dgsvj0.f"
		    }

#line 415 "dgsvj0.f"
		    if (aapp > 0.) {

#line 417 "dgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 419 "dgsvj0.f"
			i__6 = igl + kbl - 1;
#line 419 "dgsvj0.f"
			i__5 = min(i__6,*n);
#line 419 "dgsvj0.f"
			for (q = p + 1; q <= i__5; ++q) {

#line 421 "dgsvj0.f"
			    aaqq = sva[q];
#line 423 "dgsvj0.f"
			    if (aaqq > 0.) {

#line 425 "dgsvj0.f"
				aapp0 = aapp;
#line 426 "dgsvj0.f"
				if (aaqq >= 1.) {
#line 427 "dgsvj0.f"
				    rotok = small * aapp <= aaqq;
#line 428 "dgsvj0.f"
				    if (aapp < big / aaqq) {
#line 429 "dgsvj0.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 432 "dgsvj0.f"
				    } else {
#line 433 "dgsvj0.f"
					dcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 434 "dgsvj0.f"
					dlascl_("G", &c__0, &c__0, &aapp, &
						d__[p], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 436 "dgsvj0.f"
					aapq = ddot_(m, &work[1], &c__1, &a[q 
						* a_dim1 + 1], &c__1) * d__[q]
						 / aaqq;
#line 438 "dgsvj0.f"
				    }
#line 439 "dgsvj0.f"
				} else {
#line 440 "dgsvj0.f"
				    rotok = aapp <= aaqq / small;
#line 441 "dgsvj0.f"
				    if (aapp > small / aaqq) {
#line 442 "dgsvj0.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 445 "dgsvj0.f"
				    } else {
#line 446 "dgsvj0.f"
					dcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 447 "dgsvj0.f"
					dlascl_("G", &c__0, &c__0, &aaqq, &
						d__[q], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 449 "dgsvj0.f"
					aapq = ddot_(m, &work[1], &c__1, &a[p 
						* a_dim1 + 1], &c__1) * d__[p]
						 / aapp;
#line 451 "dgsvj0.f"
				    }
#line 452 "dgsvj0.f"
				}

/* Computing MAX */
#line 454 "dgsvj0.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 454 "dgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 458 "dgsvj0.f"
				if (abs(aapq) > *tol) {

/*           .. rotate */
/*           ROTATED = ROTATED + ONE */

#line 463 "dgsvj0.f"
				    if (ir1 == 0) {
#line 464 "dgsvj0.f"
					notrot = 0;
#line 465 "dgsvj0.f"
					pskipped = 0;
#line 466 "dgsvj0.f"
					++iswrot;
#line 467 "dgsvj0.f"
				    }

#line 469 "dgsvj0.f"
				    if (rotok) {

#line 471 "dgsvj0.f"
					aqoap = aaqq / aapp;
#line 472 "dgsvj0.f"
					apoaq = aapp / aaqq;
#line 473 "dgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;

#line 475 "dgsvj0.f"
					if (abs(theta) > bigtheta) {

#line 477 "dgsvj0.f"
					    t = .5 / theta;
#line 478 "dgsvj0.f"
					    fastr[2] = t * d__[p] / d__[q];
#line 479 "dgsvj0.f"
					    fastr[3] = -t * d__[q] / d__[p];
#line 480 "dgsvj0.f"
					    drotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 482 "dgsvj0.f"
					    if (rsvec) {
#line 482 "dgsvj0.f"
			  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 482 "dgsvj0.f"
					    }
/* Computing MAX */
#line 486 "dgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 486 "dgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 488 "dgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 488 "dgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 490 "dgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 490 "dgsvj0.f"
					    mxsinj = max(d__1,d__2);

#line 492 "dgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 496 "dgsvj0.f"
					    thsign = -d_sign(&c_b42, &aapq);
#line 497 "dgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 499 "dgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 500 "dgsvj0.f"
					    sn = t * cs;

/* Computing MAX */
#line 502 "dgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 502 "dgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 503 "dgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 503 "dgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 505 "dgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 505 "dgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 508 "dgsvj0.f"
					    apoaq = d__[p] / d__[q];
#line 509 "dgsvj0.f"
					    aqoap = d__[q] / d__[p];
#line 510 "dgsvj0.f"
					    if (d__[p] >= 1.) {
#line 511 "dgsvj0.f"
			  if (d__[q] >= 1.) {
#line 512 "dgsvj0.f"
			      fastr[2] = t * apoaq;
#line 513 "dgsvj0.f"
			      fastr[3] = -t * aqoap;
#line 514 "dgsvj0.f"
			      d__[p] *= cs;
#line 515 "dgsvj0.f"
			      d__[q] *= cs;
#line 516 "dgsvj0.f"
			      drotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 519 "dgsvj0.f"
			      if (rsvec) {
#line 519 "dgsvj0.f"
				  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 519 "dgsvj0.f"
			      }
#line 522 "dgsvj0.f"
			  } else {
#line 523 "dgsvj0.f"
			      d__1 = -t * aqoap;
#line 523 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 526 "dgsvj0.f"
			      d__1 = cs * sn * apoaq;
#line 526 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 529 "dgsvj0.f"
			      d__[p] *= cs;
#line 530 "dgsvj0.f"
			      d__[q] /= cs;
#line 531 "dgsvj0.f"
			      if (rsvec) {
#line 532 "dgsvj0.f"
				  d__1 = -t * aqoap;
#line 532 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 535 "dgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 535 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 539 "dgsvj0.f"
			      }
#line 540 "dgsvj0.f"
			  }
#line 541 "dgsvj0.f"
					    } else {
#line 542 "dgsvj0.f"
			  if (d__[q] >= 1.) {
#line 543 "dgsvj0.f"
			      d__1 = t * apoaq;
#line 543 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 546 "dgsvj0.f"
			      d__1 = -cs * sn * aqoap;
#line 546 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 549 "dgsvj0.f"
			      d__[p] /= cs;
#line 550 "dgsvj0.f"
			      d__[q] *= cs;
#line 551 "dgsvj0.f"
			      if (rsvec) {
#line 552 "dgsvj0.f"
				  d__1 = t * apoaq;
#line 552 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 555 "dgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 555 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 559 "dgsvj0.f"
			      }
#line 560 "dgsvj0.f"
			  } else {
#line 561 "dgsvj0.f"
			      if (d__[p] >= d__[q]) {
#line 562 "dgsvj0.f"
				  d__1 = -t * aqoap;
#line 562 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 565 "dgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 565 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 568 "dgsvj0.f"
				  d__[p] *= cs;
#line 569 "dgsvj0.f"
				  d__[q] /= cs;
#line 570 "dgsvj0.f"
				  if (rsvec) {
#line 571 "dgsvj0.f"
				      d__1 = -t * aqoap;
#line 571 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 575 "dgsvj0.f"
				      d__1 = cs * sn * apoaq;
#line 575 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 579 "dgsvj0.f"
				  }
#line 580 "dgsvj0.f"
			      } else {
#line 581 "dgsvj0.f"
				  d__1 = t * apoaq;
#line 581 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 584 "dgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 584 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 588 "dgsvj0.f"
				  d__[p] /= cs;
#line 589 "dgsvj0.f"
				  d__[q] *= cs;
#line 590 "dgsvj0.f"
				  if (rsvec) {
#line 591 "dgsvj0.f"
				      d__1 = t * apoaq;
#line 591 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 594 "dgsvj0.f"
				      d__1 = -cs * sn * aqoap;
#line 594 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 598 "dgsvj0.f"
				  }
#line 599 "dgsvj0.f"
			      }
#line 600 "dgsvj0.f"
			  }
#line 601 "dgsvj0.f"
					    }
#line 602 "dgsvj0.f"
					}

#line 604 "dgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 606 "dgsvj0.f"
					dcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 607 "dgsvj0.f"
					dlascl_("G", &c__0, &c__0, &aapp, &
						c_b42, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 609 "dgsvj0.f"
					dlascl_("G", &c__0, &c__0, &aaqq, &
						c_b42, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 611 "dgsvj0.f"
					temp1 = -aapq * d__[p] / d__[q];
#line 612 "dgsvj0.f"
					daxpy_(m, &temp1, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 614 "dgsvj0.f"
					dlascl_("G", &c__0, &c__0, &c_b42, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 616 "dgsvj0.f"
					d__1 = 0., d__2 = 1. - aapq * aapq;
#line 616 "dgsvj0.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 618 "dgsvj0.f"
					mxsinj = max(mxsinj,*sfmin);
#line 619 "dgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */
/* Computing 2nd power */
#line 624 "dgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 624 "dgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 626 "dgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 628 "dgsvj0.f"
					    sva[q] = dnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * d__[q];
#line 630 "dgsvj0.f"
					} else {
#line 631 "dgsvj0.f"
					    t = 0.;
#line 632 "dgsvj0.f"
					    aaqq = 1.;
#line 633 "dgsvj0.f"
					    dlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 635 "dgsvj0.f"
					    sva[q] = t * sqrt(aaqq) * d__[q];
#line 636 "dgsvj0.f"
					}
#line 637 "dgsvj0.f"
				    }
#line 638 "dgsvj0.f"
				    if (aapp / aapp0 <= rooteps) {
#line 639 "dgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 641 "dgsvj0.f"
					    aapp = dnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * d__[p];
#line 643 "dgsvj0.f"
					} else {
#line 644 "dgsvj0.f"
					    t = 0.;
#line 645 "dgsvj0.f"
					    aapp = 1.;
#line 646 "dgsvj0.f"
					    dlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 648 "dgsvj0.f"
					    aapp = t * sqrt(aapp) * d__[p];
#line 649 "dgsvj0.f"
					}
#line 650 "dgsvj0.f"
					sva[p] = aapp;
#line 651 "dgsvj0.f"
				    }

#line 653 "dgsvj0.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 655 "dgsvj0.f"
				    if (ir1 == 0) {
#line 655 "dgsvj0.f"
					++notrot;
#line 655 "dgsvj0.f"
				    }
#line 656 "dgsvj0.f"
				    ++pskipped;
#line 657 "dgsvj0.f"
				}
#line 658 "dgsvj0.f"
			    } else {
/*        A(:,q) is zero column */
#line 660 "dgsvj0.f"
				if (ir1 == 0) {
#line 660 "dgsvj0.f"
				    ++notrot;
#line 660 "dgsvj0.f"
				}
#line 661 "dgsvj0.f"
				++pskipped;
#line 662 "dgsvj0.f"
			    }

#line 664 "dgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 666 "dgsvj0.f"
				if (ir1 == 0) {
#line 666 "dgsvj0.f"
				    aapp = -aapp;
#line 666 "dgsvj0.f"
				}
#line 667 "dgsvj0.f"
				notrot = 0;
#line 668 "dgsvj0.f"
				goto L2103;
#line 669 "dgsvj0.f"
			    }

#line 671 "dgsvj0.f"
/* L2002: */
#line 671 "dgsvj0.f"
			}
/*     END q-LOOP */

#line 674 "dgsvj0.f"
L2103:
/*     bailed out of q-loop */
#line 677 "dgsvj0.f"
			sva[p] = aapp;
#line 679 "dgsvj0.f"
		    } else {
#line 680 "dgsvj0.f"
			sva[p] = aapp;
#line 681 "dgsvj0.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 681 "dgsvj0.f"
			    i__5 = igl + kbl - 1;
#line 681 "dgsvj0.f"
			    notrot = notrot + min(i__5,*n) - p;
#line 681 "dgsvj0.f"
			}
#line 683 "dgsvj0.f"
		    }

#line 685 "dgsvj0.f"
/* L2001: */
#line 685 "dgsvj0.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 688 "dgsvj0.f"
/* L1002: */
#line 688 "dgsvj0.f"
	    }
/*     end of ir1-loop */

/* ........................................................ */
/* ... go to the off diagonal blocks */

#line 694 "dgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

#line 696 "dgsvj0.f"
	    i__3 = nbl;
#line 696 "dgsvj0.f"
	    for (jbc = ibr + 1; jbc <= i__3; ++jbc) {

#line 698 "dgsvj0.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 702 "dgsvj0.f"
		ijblsk = 0;
/* Computing MIN */
#line 703 "dgsvj0.f"
		i__5 = igl + kbl - 1;
#line 703 "dgsvj0.f"
		i__4 = min(i__5,*n);
#line 703 "dgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

#line 705 "dgsvj0.f"
		    aapp = sva[p];

#line 707 "dgsvj0.f"
		    if (aapp > 0.) {

#line 709 "dgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 711 "dgsvj0.f"
			i__6 = jgl + kbl - 1;
#line 711 "dgsvj0.f"
			i__5 = min(i__6,*n);
#line 711 "dgsvj0.f"
			for (q = jgl; q <= i__5; ++q) {

#line 713 "dgsvj0.f"
			    aaqq = sva[q];

#line 715 "dgsvj0.f"
			    if (aaqq > 0.) {
#line 716 "dgsvj0.f"
				aapp0 = aapp;

/*     -#- M x 2 Jacobi SVD -#- */

/*        -#- Safe Gram matrix computation -#- */

#line 722 "dgsvj0.f"
				if (aaqq >= 1.) {
#line 723 "dgsvj0.f"
				    if (aapp >= aaqq) {
#line 724 "dgsvj0.f"
					rotok = small * aapp <= aaqq;
#line 725 "dgsvj0.f"
				    } else {
#line 726 "dgsvj0.f"
					rotok = small * aaqq <= aapp;
#line 727 "dgsvj0.f"
				    }
#line 728 "dgsvj0.f"
				    if (aapp < big / aaqq) {
#line 729 "dgsvj0.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 732 "dgsvj0.f"
				    } else {
#line 733 "dgsvj0.f"
					dcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 734 "dgsvj0.f"
					dlascl_("G", &c__0, &c__0, &aapp, &
						d__[p], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 736 "dgsvj0.f"
					aapq = ddot_(m, &work[1], &c__1, &a[q 
						* a_dim1 + 1], &c__1) * d__[q]
						 / aaqq;
#line 738 "dgsvj0.f"
				    }
#line 739 "dgsvj0.f"
				} else {
#line 740 "dgsvj0.f"
				    if (aapp >= aaqq) {
#line 741 "dgsvj0.f"
					rotok = aapp <= aaqq / small;
#line 742 "dgsvj0.f"
				    } else {
#line 743 "dgsvj0.f"
					rotok = aaqq <= aapp / small;
#line 744 "dgsvj0.f"
				    }
#line 745 "dgsvj0.f"
				    if (aapp > small / aaqq) {
#line 746 "dgsvj0.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 749 "dgsvj0.f"
				    } else {
#line 750 "dgsvj0.f"
					dcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 751 "dgsvj0.f"
					dlascl_("G", &c__0, &c__0, &aaqq, &
						d__[q], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 753 "dgsvj0.f"
					aapq = ddot_(m, &work[1], &c__1, &a[p 
						* a_dim1 + 1], &c__1) * d__[p]
						 / aapp;
#line 755 "dgsvj0.f"
				    }
#line 756 "dgsvj0.f"
				}

/* Computing MAX */
#line 758 "dgsvj0.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 758 "dgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 762 "dgsvj0.f"
				if (abs(aapq) > *tol) {
#line 763 "dgsvj0.f"
				    notrot = 0;
/*           ROTATED  = ROTATED + 1 */
#line 765 "dgsvj0.f"
				    pskipped = 0;
#line 766 "dgsvj0.f"
				    ++iswrot;

#line 768 "dgsvj0.f"
				    if (rotok) {

#line 770 "dgsvj0.f"
					aqoap = aaqq / aapp;
#line 771 "dgsvj0.f"
					apoaq = aapp / aaqq;
#line 772 "dgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;
#line 773 "dgsvj0.f"
					if (aaqq > aapp0) {
#line 773 "dgsvj0.f"
					    theta = -theta;
#line 773 "dgsvj0.f"
					}

#line 775 "dgsvj0.f"
					if (abs(theta) > bigtheta) {
#line 776 "dgsvj0.f"
					    t = .5 / theta;
#line 777 "dgsvj0.f"
					    fastr[2] = t * d__[p] / d__[q];
#line 778 "dgsvj0.f"
					    fastr[3] = -t * d__[q] / d__[p];
#line 779 "dgsvj0.f"
					    drotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 781 "dgsvj0.f"
					    if (rsvec) {
#line 781 "dgsvj0.f"
			  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 781 "dgsvj0.f"
					    }
/* Computing MAX */
#line 785 "dgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 785 "dgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 787 "dgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 787 "dgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 789 "dgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 789 "dgsvj0.f"
					    mxsinj = max(d__1,d__2);
#line 790 "dgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 794 "dgsvj0.f"
					    thsign = -d_sign(&c_b42, &aapq);
#line 795 "dgsvj0.f"
					    if (aaqq > aapp0) {
#line 795 "dgsvj0.f"
			  thsign = -thsign;
#line 795 "dgsvj0.f"
					    }
#line 796 "dgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 798 "dgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 799 "dgsvj0.f"
					    sn = t * cs;
/* Computing MAX */
#line 800 "dgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 800 "dgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 801 "dgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 801 "dgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 803 "dgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 803 "dgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 806 "dgsvj0.f"
					    apoaq = d__[p] / d__[q];
#line 807 "dgsvj0.f"
					    aqoap = d__[q] / d__[p];
#line 808 "dgsvj0.f"
					    if (d__[p] >= 1.) {

#line 810 "dgsvj0.f"
			  if (d__[q] >= 1.) {
#line 811 "dgsvj0.f"
			      fastr[2] = t * apoaq;
#line 812 "dgsvj0.f"
			      fastr[3] = -t * aqoap;
#line 813 "dgsvj0.f"
			      d__[p] *= cs;
#line 814 "dgsvj0.f"
			      d__[q] *= cs;
#line 815 "dgsvj0.f"
			      drotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 818 "dgsvj0.f"
			      if (rsvec) {
#line 818 "dgsvj0.f"
				  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 818 "dgsvj0.f"
			      }
#line 821 "dgsvj0.f"
			  } else {
#line 822 "dgsvj0.f"
			      d__1 = -t * aqoap;
#line 822 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 825 "dgsvj0.f"
			      d__1 = cs * sn * apoaq;
#line 825 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 828 "dgsvj0.f"
			      if (rsvec) {
#line 829 "dgsvj0.f"
				  d__1 = -t * aqoap;
#line 829 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 832 "dgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 832 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 836 "dgsvj0.f"
			      }
#line 837 "dgsvj0.f"
			      d__[p] *= cs;
#line 838 "dgsvj0.f"
			      d__[q] /= cs;
#line 839 "dgsvj0.f"
			  }
#line 840 "dgsvj0.f"
					    } else {
#line 841 "dgsvj0.f"
			  if (d__[q] >= 1.) {
#line 842 "dgsvj0.f"
			      d__1 = t * apoaq;
#line 842 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 845 "dgsvj0.f"
			      d__1 = -cs * sn * aqoap;
#line 845 "dgsvj0.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 848 "dgsvj0.f"
			      if (rsvec) {
#line 849 "dgsvj0.f"
				  d__1 = t * apoaq;
#line 849 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 852 "dgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 852 "dgsvj0.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 856 "dgsvj0.f"
			      }
#line 857 "dgsvj0.f"
			      d__[p] /= cs;
#line 858 "dgsvj0.f"
			      d__[q] *= cs;
#line 859 "dgsvj0.f"
			  } else {
#line 860 "dgsvj0.f"
			      if (d__[p] >= d__[q]) {
#line 861 "dgsvj0.f"
				  d__1 = -t * aqoap;
#line 861 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 864 "dgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 864 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 867 "dgsvj0.f"
				  d__[p] *= cs;
#line 868 "dgsvj0.f"
				  d__[q] /= cs;
#line 869 "dgsvj0.f"
				  if (rsvec) {
#line 870 "dgsvj0.f"
				      d__1 = -t * aqoap;
#line 870 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 874 "dgsvj0.f"
				      d__1 = cs * sn * apoaq;
#line 874 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 878 "dgsvj0.f"
				  }
#line 879 "dgsvj0.f"
			      } else {
#line 880 "dgsvj0.f"
				  d__1 = t * apoaq;
#line 880 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 883 "dgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 883 "dgsvj0.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 887 "dgsvj0.f"
				  d__[p] /= cs;
#line 888 "dgsvj0.f"
				  d__[q] *= cs;
#line 889 "dgsvj0.f"
				  if (rsvec) {
#line 890 "dgsvj0.f"
				      d__1 = t * apoaq;
#line 890 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 893 "dgsvj0.f"
				      d__1 = -cs * sn * aqoap;
#line 893 "dgsvj0.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 897 "dgsvj0.f"
				  }
#line 898 "dgsvj0.f"
			      }
#line 899 "dgsvj0.f"
			  }
#line 900 "dgsvj0.f"
					    }
#line 901 "dgsvj0.f"
					}

#line 903 "dgsvj0.f"
				    } else {
#line 904 "dgsvj0.f"
					if (aapp > aaqq) {
#line 905 "dgsvj0.f"
					    dcopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 907 "dgsvj0.f"
					    dlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 909 "dgsvj0.f"
					    dlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 912 "dgsvj0.f"
					    temp1 = -aapq * d__[p] / d__[q];
#line 913 "dgsvj0.f"
					    daxpy_(m, &temp1, &work[1], &c__1,
						     &a[q * a_dim1 + 1], &
						    c__1);
#line 915 "dgsvj0.f"
					    dlascl_("G", &c__0, &c__0, &c_b42,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 918 "dgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 918 "dgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 920 "dgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 921 "dgsvj0.f"
					} else {
#line 922 "dgsvj0.f"
					    dcopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 924 "dgsvj0.f"
					    dlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 926 "dgsvj0.f"
					    dlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 929 "dgsvj0.f"
					    temp1 = -aapq * d__[q] / d__[p];
#line 930 "dgsvj0.f"
					    daxpy_(m, &temp1, &work[1], &c__1,
						     &a[p * a_dim1 + 1], &
						    c__1);
#line 932 "dgsvj0.f"
					    dlascl_("G", &c__0, &c__0, &c_b42,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 935 "dgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 935 "dgsvj0.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 937 "dgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 938 "dgsvj0.f"
					}
#line 939 "dgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q) */
/*           .. recompute SVA(q) */
/* Computing 2nd power */
#line 944 "dgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 944 "dgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 946 "dgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 948 "dgsvj0.f"
					    sva[q] = dnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * d__[q];
#line 950 "dgsvj0.f"
					} else {
#line 951 "dgsvj0.f"
					    t = 0.;
#line 952 "dgsvj0.f"
					    aaqq = 1.;
#line 953 "dgsvj0.f"
					    dlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 955 "dgsvj0.f"
					    sva[q] = t * sqrt(aaqq) * d__[q];
#line 956 "dgsvj0.f"
					}
#line 957 "dgsvj0.f"
				    }
/* Computing 2nd power */
#line 958 "dgsvj0.f"
				    d__1 = aapp / aapp0;
#line 958 "dgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 959 "dgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 961 "dgsvj0.f"
					    aapp = dnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * d__[p];
#line 963 "dgsvj0.f"
					} else {
#line 964 "dgsvj0.f"
					    t = 0.;
#line 965 "dgsvj0.f"
					    aapp = 1.;
#line 966 "dgsvj0.f"
					    dlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 968 "dgsvj0.f"
					    aapp = t * sqrt(aapp) * d__[p];
#line 969 "dgsvj0.f"
					}
#line 970 "dgsvj0.f"
					sva[p] = aapp;
#line 971 "dgsvj0.f"
				    }
/*              end of OK rotation */
#line 973 "dgsvj0.f"
				} else {
#line 974 "dgsvj0.f"
				    ++notrot;
#line 975 "dgsvj0.f"
				    ++pskipped;
#line 976 "dgsvj0.f"
				    ++ijblsk;
#line 977 "dgsvj0.f"
				}
#line 978 "dgsvj0.f"
			    } else {
#line 979 "dgsvj0.f"
				++notrot;
#line 980 "dgsvj0.f"
				++pskipped;
#line 981 "dgsvj0.f"
				++ijblsk;
#line 982 "dgsvj0.f"
			    }

#line 984 "dgsvj0.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 986 "dgsvj0.f"
				sva[p] = aapp;
#line 987 "dgsvj0.f"
				notrot = 0;
#line 988 "dgsvj0.f"
				goto L2011;
#line 989 "dgsvj0.f"
			    }
#line 990 "dgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 992 "dgsvj0.f"
				aapp = -aapp;
#line 993 "dgsvj0.f"
				notrot = 0;
#line 994 "dgsvj0.f"
				goto L2203;
#line 995 "dgsvj0.f"
			    }

#line 997 "dgsvj0.f"
/* L2200: */
#line 997 "dgsvj0.f"
			}
/*        end of the q-loop */
#line 999 "dgsvj0.f"
L2203:

#line 1001 "dgsvj0.f"
			sva[p] = aapp;

#line 1003 "dgsvj0.f"
		    } else {
#line 1004 "dgsvj0.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1004 "dgsvj0.f"
			    i__5 = jgl + kbl - 1;
#line 1004 "dgsvj0.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 1004 "dgsvj0.f"
			}
#line 1006 "dgsvj0.f"
			if (aapp < 0.) {
#line 1006 "dgsvj0.f"
			    notrot = 0;
#line 1006 "dgsvj0.f"
			}
#line 1007 "dgsvj0.f"
		    }
#line 1009 "dgsvj0.f"
/* L2100: */
#line 1009 "dgsvj0.f"
		}
/*     end of the p-loop */
#line 1011 "dgsvj0.f"
/* L2010: */
#line 1011 "dgsvj0.f"
	    }
/*     end of the jbc-loop */
#line 1013 "dgsvj0.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1015 "dgsvj0.f"
	    i__4 = igl + kbl - 1;
#line 1015 "dgsvj0.f"
	    i__3 = min(i__4,*n);
#line 1015 "dgsvj0.f"
	    for (p = igl; p <= i__3; ++p) {
#line 1016 "dgsvj0.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1017 "dgsvj0.f"
/* L2012: */
#line 1017 "dgsvj0.f"
	    }

#line 1019 "dgsvj0.f"
/* L2000: */
#line 1019 "dgsvj0.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1023 "dgsvj0.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1025 "dgsvj0.f"
	    sva[*n] = dnrm2_(m, &a[*n * a_dim1 + 1], &c__1) * d__[*n];
#line 1026 "dgsvj0.f"
	} else {
#line 1027 "dgsvj0.f"
	    t = 0.;
#line 1028 "dgsvj0.f"
	    aapp = 1.;
#line 1029 "dgsvj0.f"
	    dlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1030 "dgsvj0.f"
	    sva[*n] = t * sqrt(aapp) * d__[*n];
#line 1031 "dgsvj0.f"
	}

/*     Additional steering devices */

#line 1035 "dgsvj0.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1035 "dgsvj0.f"
	    swband = i__;
#line 1035 "dgsvj0.f"
	}

#line 1038 "dgsvj0.f"
	if (i__ > swband + 1 && mxaapq < (doublereal) (*n) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 1040 "dgsvj0.f"
	    goto L1994;
#line 1041 "dgsvj0.f"
	}

#line 1043 "dgsvj0.f"
	if (notrot >= emptsw) {
#line 1043 "dgsvj0.f"
	    goto L1994;
#line 1043 "dgsvj0.f"
	}
#line 1045 "dgsvj0.f"
/* L1993: */
#line 1045 "dgsvj0.f"
    }
/*     end i=1:NSWEEP loop */
/* #:) Reaching this point means that the procedure has comleted the given */
/*     number of iterations. */
#line 1049 "dgsvj0.f"
    *info = *nsweep - 1;
#line 1050 "dgsvj0.f"
    goto L1995;
#line 1051 "dgsvj0.f"
L1994:
/* #:) Reaching this point means that during the i-th sweep all pivots were */
/*     below the given tolerance, causing early exit. */

#line 1055 "dgsvj0.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1057 "dgsvj0.f"
L1995:

/*     Sort the vector D. */
#line 1060 "dgsvj0.f"
    i__1 = *n - 1;
#line 1060 "dgsvj0.f"
    for (p = 1; p <= i__1; ++p) {
#line 1061 "dgsvj0.f"
	i__2 = *n - p + 1;
#line 1061 "dgsvj0.f"
	q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1062 "dgsvj0.f"
	if (p != q) {
#line 1063 "dgsvj0.f"
	    temp1 = sva[p];
#line 1064 "dgsvj0.f"
	    sva[p] = sva[q];
#line 1065 "dgsvj0.f"
	    sva[q] = temp1;
#line 1066 "dgsvj0.f"
	    temp1 = d__[p];
#line 1067 "dgsvj0.f"
	    d__[p] = d__[q];
#line 1068 "dgsvj0.f"
	    d__[q] = temp1;
#line 1069 "dgsvj0.f"
	    dswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1070 "dgsvj0.f"
	    if (rsvec) {
#line 1070 "dgsvj0.f"
		dswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1070 "dgsvj0.f"
	    }
#line 1071 "dgsvj0.f"
	}
#line 1072 "dgsvj0.f"
/* L5991: */
#line 1072 "dgsvj0.f"
    }

#line 1074 "dgsvj0.f"
    return 0;
/*     .. */
/*     .. END OF DGSVJ0 */
/*     .. */
} /* dgsvj0_ */

