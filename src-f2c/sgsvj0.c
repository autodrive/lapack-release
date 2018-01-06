#line 1 "sgsvj0.f"
/* sgsvj0.f -- translated by f2c (version 20100827).
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

#line 1 "sgsvj0.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b42 = 1.;

/* > \brief \b SGSVJ0 pre-processor for the routine sgesvj. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGSVJ0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgsvj0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgsvj0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgsvj0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, */
/*                          SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP */
/*       REAL               EPS, SFMIN, TOL */
/*       CHARACTER*1        JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), SVA( N ), D( N ), V( LDV, * ), */
/*      $                   WORK( LWORK ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGSVJ0 is called from SGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as SGESVJ does, but */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          D is REAL array, dimension (N) */
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
/* >          WORK is REAL array, dimension (LWORK) */
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

/* > \date November 2017 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > SGSVJ0 is used just to enable SGESVJ to call a simplified version of */
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
/* Subroutine */ int sgsvj0_(char *jobv, integer *m, integer *n, doublereal *
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
    static integer ierr;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal aapp0, temp1;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta, small, fastr[5];
    static logical applv, rsvec;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical rotok;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), srotm_(integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *)
	    , xerbla_(char *, integer *, ftnlen);
    static integer ijblsk, swband;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static integer blskip;
    static doublereal mxaapq, thsign;
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal mxsinj;
    static integer emptsw, notrot, iswrot, lkahead;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 272 "sgsvj0.f"
    /* Parameter adjustments */
#line 272 "sgsvj0.f"
    --sva;
#line 272 "sgsvj0.f"
    --d__;
#line 272 "sgsvj0.f"
    a_dim1 = *lda;
#line 272 "sgsvj0.f"
    a_offset = 1 + a_dim1;
#line 272 "sgsvj0.f"
    a -= a_offset;
#line 272 "sgsvj0.f"
    v_dim1 = *ldv;
#line 272 "sgsvj0.f"
    v_offset = 1 + v_dim1;
#line 272 "sgsvj0.f"
    v -= v_offset;
#line 272 "sgsvj0.f"
    --work;
#line 272 "sgsvj0.f"

#line 272 "sgsvj0.f"
    /* Function Body */
#line 272 "sgsvj0.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 273 "sgsvj0.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 274 "sgsvj0.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 275 "sgsvj0.f"
	*info = -1;
#line 276 "sgsvj0.f"
    } else if (*m < 0) {
#line 277 "sgsvj0.f"
	*info = -2;
#line 278 "sgsvj0.f"
    } else if (*n < 0 || *n > *m) {
#line 279 "sgsvj0.f"
	*info = -3;
#line 280 "sgsvj0.f"
    } else if (*lda < *m) {
#line 281 "sgsvj0.f"
	*info = -5;
#line 282 "sgsvj0.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 283 "sgsvj0.f"
	*info = -8;
#line 284 "sgsvj0.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 286 "sgsvj0.f"
	*info = -10;
#line 287 "sgsvj0.f"
    } else if (*tol <= *eps) {
#line 288 "sgsvj0.f"
	*info = -13;
#line 289 "sgsvj0.f"
    } else if (*nsweep < 0) {
#line 290 "sgsvj0.f"
	*info = -14;
#line 291 "sgsvj0.f"
    } else if (*lwork < *m) {
#line 292 "sgsvj0.f"
	*info = -16;
#line 293 "sgsvj0.f"
    } else {
#line 294 "sgsvj0.f"
	*info = 0;
#line 295 "sgsvj0.f"
    }

/*     #:( */
#line 298 "sgsvj0.f"
    if (*info != 0) {
#line 299 "sgsvj0.f"
	i__1 = -(*info);
#line 299 "sgsvj0.f"
	xerbla_("SGSVJ0", &i__1, (ftnlen)6);
#line 300 "sgsvj0.f"
	return 0;
#line 301 "sgsvj0.f"
    }

#line 303 "sgsvj0.f"
    if (rsvec) {
#line 304 "sgsvj0.f"
	mvl = *n;
#line 305 "sgsvj0.f"
    } else if (applv) {
#line 306 "sgsvj0.f"
	mvl = *mv;
#line 307 "sgsvj0.f"
    }
#line 308 "sgsvj0.f"
    rsvec = rsvec || applv;
#line 310 "sgsvj0.f"
    rooteps = sqrt(*eps);
#line 311 "sgsvj0.f"
    rootsfmin = sqrt(*sfmin);
#line 312 "sgsvj0.f"
    small = *sfmin / *eps;
#line 313 "sgsvj0.f"
    big = 1. / *sfmin;
#line 314 "sgsvj0.f"
    rootbig = 1. / rootsfmin;
#line 315 "sgsvj0.f"
    bigtheta = 1. / rooteps;
#line 316 "sgsvj0.f"
    roottol = sqrt(*tol);

/*     .. Row-cyclic Jacobi SVD algorithm with column pivoting .. */

#line 320 "sgsvj0.f"
    emptsw = *n * (*n - 1) / 2;
#line 321 "sgsvj0.f"
    notrot = 0;
#line 322 "sgsvj0.f"
    fastr[0] = 0.;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 327 "sgsvj0.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter. It is meaningful and effective */
/*     if SGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure */
/*     ...... */
#line 333 "sgsvj0.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 339 "sgsvj0.f"
    nbl = *n / kbl;
#line 340 "sgsvj0.f"
    if (nbl * kbl != *n) {
#line 340 "sgsvj0.f"
	++nbl;
#line 340 "sgsvj0.f"
    }
/* Computing 2nd power */
#line 342 "sgsvj0.f"
    i__1 = kbl;
#line 342 "sgsvj0.f"
    blskip = i__1 * i__1 + 1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
#line 345 "sgsvj0.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */
#line 348 "sgsvj0.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */
#line 350 "sgsvj0.f"
    swband = 0;
#line 351 "sgsvj0.f"
    pskipped = 0;

#line 353 "sgsvj0.f"
    i__1 = *nsweep;
#line 353 "sgsvj0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     .. go go go ... */

#line 356 "sgsvj0.f"
	mxaapq = 0.;
#line 357 "sgsvj0.f"
	mxsinj = 0.;
#line 358 "sgsvj0.f"
	iswrot = 0;

#line 360 "sgsvj0.f"
	notrot = 0;
#line 361 "sgsvj0.f"
	pskipped = 0;

#line 363 "sgsvj0.f"
	i__2 = nbl;
#line 363 "sgsvj0.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {
#line 365 "sgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 367 "sgsvj0.f"
	    i__4 = lkahead, i__5 = nbl - ibr;
#line 367 "sgsvj0.f"
	    i__3 = min(i__4,i__5);
#line 367 "sgsvj0.f"
	    for (ir1 = 0; ir1 <= i__3; ++ir1) {

#line 369 "sgsvj0.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 371 "sgsvj0.f"
		i__5 = igl + kbl - 1, i__6 = *n - 1;
#line 371 "sgsvj0.f"
		i__4 = min(i__5,i__6);
#line 371 "sgsvj0.f"
		for (p = igl; p <= i__4; ++p) {
/*     .. de Rijk's pivoting */
#line 374 "sgsvj0.f"
		    i__5 = *n - p + 1;
#line 374 "sgsvj0.f"
		    q = isamax_(&i__5, &sva[p], &c__1) + p - 1;
#line 375 "sgsvj0.f"
		    if (p != q) {
#line 376 "sgsvj0.f"
			sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 377 "sgsvj0.f"
			if (rsvec) {
#line 377 "sgsvj0.f"
			    sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 377 "sgsvj0.f"
			}
#line 379 "sgsvj0.f"
			temp1 = sva[p];
#line 380 "sgsvj0.f"
			sva[p] = sva[q];
#line 381 "sgsvj0.f"
			sva[q] = temp1;
#line 382 "sgsvj0.f"
			temp1 = d__[p];
#line 383 "sgsvj0.f"
			d__[p] = d__[q];
#line 384 "sgsvj0.f"
			d__[q] = temp1;
#line 385 "sgsvj0.f"
		    }

#line 387 "sgsvj0.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/*        Caveat: */
/*        Some BLAS implementations compute SNRM2(M,A(1,p),1) */
/*        as SQRT(SDOT(M,A(1,p),1,A(1,p),1)), which may result in */
/*        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and */
/*        undeflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
/*        Hence, SNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented SNRM2 is available, the IF-THEN-ELSE */
/*        below should read "AAPP = SNRM2( M, A(1,p), 1 ) * D(p)". */

#line 401 "sgsvj0.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 403 "sgsvj0.f"
			    sva[p] = snrm2_(m, &a[p * a_dim1 + 1], &c__1) * 
				    d__[p];
#line 404 "sgsvj0.f"
			} else {
#line 405 "sgsvj0.f"
			    temp1 = 0.;
#line 406 "sgsvj0.f"
			    aapp = 1.;
#line 407 "sgsvj0.f"
			    slassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 408 "sgsvj0.f"
			    sva[p] = temp1 * sqrt(aapp) * d__[p];
#line 409 "sgsvj0.f"
			}
#line 410 "sgsvj0.f"
			aapp = sva[p];
#line 411 "sgsvj0.f"
		    } else {
#line 412 "sgsvj0.f"
			aapp = sva[p];
#line 413 "sgsvj0.f"
		    }

#line 416 "sgsvj0.f"
		    if (aapp > 0.) {

#line 418 "sgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 420 "sgsvj0.f"
			i__6 = igl + kbl - 1;
#line 420 "sgsvj0.f"
			i__5 = min(i__6,*n);
#line 420 "sgsvj0.f"
			for (q = p + 1; q <= i__5; ++q) {

#line 422 "sgsvj0.f"
			    aaqq = sva[q];
#line 424 "sgsvj0.f"
			    if (aaqq > 0.) {

#line 426 "sgsvj0.f"
				aapp0 = aapp;
#line 427 "sgsvj0.f"
				if (aaqq >= 1.) {
#line 428 "sgsvj0.f"
				    rotok = small * aapp <= aaqq;
#line 429 "sgsvj0.f"
				    if (aapp < big / aaqq) {
#line 430 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 433 "sgsvj0.f"
				    } else {
#line 434 "sgsvj0.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 435 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						d__[p], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 437 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[q 
						* a_dim1 + 1], &c__1) * d__[q]
						 / aaqq;
#line 439 "sgsvj0.f"
				    }
#line 440 "sgsvj0.f"
				} else {
#line 441 "sgsvj0.f"
				    rotok = aapp <= aaqq / small;
#line 442 "sgsvj0.f"
				    if (aapp > small / aaqq) {
#line 443 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 446 "sgsvj0.f"
				    } else {
#line 447 "sgsvj0.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 448 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						d__[q], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 450 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[p 
						* a_dim1 + 1], &c__1) * d__[p]
						 / aapp;
#line 452 "sgsvj0.f"
				    }
#line 453 "sgsvj0.f"
				}

/* Computing MAX */
#line 455 "sgsvj0.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 455 "sgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 459 "sgsvj0.f"
				if (abs(aapq) > *tol) {

/*           .. rotate */
/*           ROTATED = ROTATED + ONE */

#line 464 "sgsvj0.f"
				    if (ir1 == 0) {
#line 465 "sgsvj0.f"
					notrot = 0;
#line 466 "sgsvj0.f"
					pskipped = 0;
#line 467 "sgsvj0.f"
					++iswrot;
#line 468 "sgsvj0.f"
				    }

#line 470 "sgsvj0.f"
				    if (rotok) {

#line 472 "sgsvj0.f"
					aqoap = aaqq / aapp;
#line 473 "sgsvj0.f"
					apoaq = aapp / aaqq;
#line 474 "sgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;

#line 476 "sgsvj0.f"
					if (abs(theta) > bigtheta) {

#line 478 "sgsvj0.f"
					    t = .5 / theta;
#line 479 "sgsvj0.f"
					    fastr[2] = t * d__[p] / d__[q];
#line 480 "sgsvj0.f"
					    fastr[3] = -t * d__[q] / d__[p];
#line 481 "sgsvj0.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 483 "sgsvj0.f"
					    if (rsvec) {
#line 483 "sgsvj0.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 483 "sgsvj0.f"
					    }
/* Computing MAX */
#line 487 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 487 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 489 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 489 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 491 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 491 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);

#line 493 "sgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 497 "sgsvj0.f"
					    thsign = -d_sign(&c_b42, &aapq);
#line 498 "sgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 500 "sgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 501 "sgsvj0.f"
					    sn = t * cs;

/* Computing MAX */
#line 503 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 503 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 504 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 504 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 506 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 506 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 509 "sgsvj0.f"
					    apoaq = d__[p] / d__[q];
#line 510 "sgsvj0.f"
					    aqoap = d__[q] / d__[p];
#line 511 "sgsvj0.f"
					    if (d__[p] >= 1.) {
#line 512 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 513 "sgsvj0.f"
			      fastr[2] = t * apoaq;
#line 514 "sgsvj0.f"
			      fastr[3] = -t * aqoap;
#line 515 "sgsvj0.f"
			      d__[p] *= cs;
#line 516 "sgsvj0.f"
			      d__[q] *= cs;
#line 517 "sgsvj0.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 520 "sgsvj0.f"
			      if (rsvec) {
#line 520 "sgsvj0.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 520 "sgsvj0.f"
			      }
#line 523 "sgsvj0.f"
			  } else {
#line 524 "sgsvj0.f"
			      d__1 = -t * aqoap;
#line 524 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 527 "sgsvj0.f"
			      d__1 = cs * sn * apoaq;
#line 527 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 530 "sgsvj0.f"
			      d__[p] *= cs;
#line 531 "sgsvj0.f"
			      d__[q] /= cs;
#line 532 "sgsvj0.f"
			      if (rsvec) {
#line 533 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 533 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 536 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 536 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 540 "sgsvj0.f"
			      }
#line 541 "sgsvj0.f"
			  }
#line 542 "sgsvj0.f"
					    } else {
#line 543 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 544 "sgsvj0.f"
			      d__1 = t * apoaq;
#line 544 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 547 "sgsvj0.f"
			      d__1 = -cs * sn * aqoap;
#line 547 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 550 "sgsvj0.f"
			      d__[p] /= cs;
#line 551 "sgsvj0.f"
			      d__[q] *= cs;
#line 552 "sgsvj0.f"
			      if (rsvec) {
#line 553 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 553 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 556 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 556 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 560 "sgsvj0.f"
			      }
#line 561 "sgsvj0.f"
			  } else {
#line 562 "sgsvj0.f"
			      if (d__[p] >= d__[q]) {
#line 563 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 563 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 566 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 566 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 569 "sgsvj0.f"
				  d__[p] *= cs;
#line 570 "sgsvj0.f"
				  d__[q] /= cs;
#line 571 "sgsvj0.f"
				  if (rsvec) {
#line 572 "sgsvj0.f"
				      d__1 = -t * aqoap;
#line 572 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 576 "sgsvj0.f"
				      d__1 = cs * sn * apoaq;
#line 576 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 580 "sgsvj0.f"
				  }
#line 581 "sgsvj0.f"
			      } else {
#line 582 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 582 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 585 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 585 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 589 "sgsvj0.f"
				  d__[p] /= cs;
#line 590 "sgsvj0.f"
				  d__[q] *= cs;
#line 591 "sgsvj0.f"
				  if (rsvec) {
#line 592 "sgsvj0.f"
				      d__1 = t * apoaq;
#line 592 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 595 "sgsvj0.f"
				      d__1 = -cs * sn * aqoap;
#line 595 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 599 "sgsvj0.f"
				  }
#line 600 "sgsvj0.f"
			      }
#line 601 "sgsvj0.f"
			  }
#line 602 "sgsvj0.f"
					    }
#line 603 "sgsvj0.f"
					}

#line 605 "sgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 607 "sgsvj0.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 608 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						c_b42, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 610 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						c_b42, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 612 "sgsvj0.f"
					temp1 = -aapq * d__[p] / d__[q];
#line 613 "sgsvj0.f"
					saxpy_(m, &temp1, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 615 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &c_b42, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 617 "sgsvj0.f"
					d__1 = 0., d__2 = 1. - aapq * aapq;
#line 617 "sgsvj0.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 619 "sgsvj0.f"
					mxsinj = max(mxsinj,*sfmin);
#line 620 "sgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */
/* Computing 2nd power */
#line 625 "sgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 625 "sgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 627 "sgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 629 "sgsvj0.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * d__[q];
#line 631 "sgsvj0.f"
					} else {
#line 632 "sgsvj0.f"
					    t = 0.;
#line 633 "sgsvj0.f"
					    aaqq = 1.;
#line 634 "sgsvj0.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 636 "sgsvj0.f"
					    sva[q] = t * sqrt(aaqq) * d__[q];
#line 637 "sgsvj0.f"
					}
#line 638 "sgsvj0.f"
				    }
#line 639 "sgsvj0.f"
				    if (aapp / aapp0 <= rooteps) {
#line 640 "sgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 642 "sgsvj0.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * d__[p];
#line 644 "sgsvj0.f"
					} else {
#line 645 "sgsvj0.f"
					    t = 0.;
#line 646 "sgsvj0.f"
					    aapp = 1.;
#line 647 "sgsvj0.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 649 "sgsvj0.f"
					    aapp = t * sqrt(aapp) * d__[p];
#line 650 "sgsvj0.f"
					}
#line 651 "sgsvj0.f"
					sva[p] = aapp;
#line 652 "sgsvj0.f"
				    }

#line 654 "sgsvj0.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 656 "sgsvj0.f"
				    if (ir1 == 0) {
#line 656 "sgsvj0.f"
					++notrot;
#line 656 "sgsvj0.f"
				    }
#line 657 "sgsvj0.f"
				    ++pskipped;
#line 658 "sgsvj0.f"
				}
#line 659 "sgsvj0.f"
			    } else {
/*        A(:,q) is zero column */
#line 661 "sgsvj0.f"
				if (ir1 == 0) {
#line 661 "sgsvj0.f"
				    ++notrot;
#line 661 "sgsvj0.f"
				}
#line 662 "sgsvj0.f"
				++pskipped;
#line 663 "sgsvj0.f"
			    }

#line 665 "sgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 667 "sgsvj0.f"
				if (ir1 == 0) {
#line 667 "sgsvj0.f"
				    aapp = -aapp;
#line 667 "sgsvj0.f"
				}
#line 668 "sgsvj0.f"
				notrot = 0;
#line 669 "sgsvj0.f"
				goto L2103;
#line 670 "sgsvj0.f"
			    }

#line 672 "sgsvj0.f"
/* L2002: */
#line 672 "sgsvj0.f"
			}
/*     END q-LOOP */

#line 675 "sgsvj0.f"
L2103:
/*     bailed out of q-loop */
#line 678 "sgsvj0.f"
			sva[p] = aapp;
#line 680 "sgsvj0.f"
		    } else {
#line 681 "sgsvj0.f"
			sva[p] = aapp;
#line 682 "sgsvj0.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 682 "sgsvj0.f"
			    i__5 = igl + kbl - 1;
#line 682 "sgsvj0.f"
			    notrot = notrot + min(i__5,*n) - p;
#line 682 "sgsvj0.f"
			}
#line 684 "sgsvj0.f"
		    }

#line 686 "sgsvj0.f"
/* L2001: */
#line 686 "sgsvj0.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 689 "sgsvj0.f"
/* L1002: */
#line 689 "sgsvj0.f"
	    }
/*     end of ir1-loop */

/* ........................................................ */
/* ... go to the off diagonal blocks */

#line 695 "sgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

#line 697 "sgsvj0.f"
	    i__3 = nbl;
#line 697 "sgsvj0.f"
	    for (jbc = ibr + 1; jbc <= i__3; ++jbc) {

#line 699 "sgsvj0.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 703 "sgsvj0.f"
		ijblsk = 0;
/* Computing MIN */
#line 704 "sgsvj0.f"
		i__5 = igl + kbl - 1;
#line 704 "sgsvj0.f"
		i__4 = min(i__5,*n);
#line 704 "sgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

#line 706 "sgsvj0.f"
		    aapp = sva[p];

#line 708 "sgsvj0.f"
		    if (aapp > 0.) {

#line 710 "sgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 712 "sgsvj0.f"
			i__6 = jgl + kbl - 1;
#line 712 "sgsvj0.f"
			i__5 = min(i__6,*n);
#line 712 "sgsvj0.f"
			for (q = jgl; q <= i__5; ++q) {

#line 714 "sgsvj0.f"
			    aaqq = sva[q];

#line 716 "sgsvj0.f"
			    if (aaqq > 0.) {
#line 717 "sgsvj0.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        .. Safe Gram matrix computation .. */

#line 723 "sgsvj0.f"
				if (aaqq >= 1.) {
#line 724 "sgsvj0.f"
				    if (aapp >= aaqq) {
#line 725 "sgsvj0.f"
					rotok = small * aapp <= aaqq;
#line 726 "sgsvj0.f"
				    } else {
#line 727 "sgsvj0.f"
					rotok = small * aaqq <= aapp;
#line 728 "sgsvj0.f"
				    }
#line 729 "sgsvj0.f"
				    if (aapp < big / aaqq) {
#line 730 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 733 "sgsvj0.f"
				    } else {
#line 734 "sgsvj0.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 735 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						d__[p], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 737 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[q 
						* a_dim1 + 1], &c__1) * d__[q]
						 / aaqq;
#line 739 "sgsvj0.f"
				    }
#line 740 "sgsvj0.f"
				} else {
#line 741 "sgsvj0.f"
				    if (aapp >= aaqq) {
#line 742 "sgsvj0.f"
					rotok = aapp <= aaqq / small;
#line 743 "sgsvj0.f"
				    } else {
#line 744 "sgsvj0.f"
					rotok = aaqq <= aapp / small;
#line 745 "sgsvj0.f"
				    }
#line 746 "sgsvj0.f"
				    if (aapp > small / aaqq) {
#line 747 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 750 "sgsvj0.f"
				    } else {
#line 751 "sgsvj0.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 752 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						d__[q], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 754 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[p 
						* a_dim1 + 1], &c__1) * d__[p]
						 / aapp;
#line 756 "sgsvj0.f"
				    }
#line 757 "sgsvj0.f"
				}

/* Computing MAX */
#line 759 "sgsvj0.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 759 "sgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 763 "sgsvj0.f"
				if (abs(aapq) > *tol) {
#line 764 "sgsvj0.f"
				    notrot = 0;
/*           ROTATED  = ROTATED + 1 */
#line 766 "sgsvj0.f"
				    pskipped = 0;
#line 767 "sgsvj0.f"
				    ++iswrot;

#line 769 "sgsvj0.f"
				    if (rotok) {

#line 771 "sgsvj0.f"
					aqoap = aaqq / aapp;
#line 772 "sgsvj0.f"
					apoaq = aapp / aaqq;
#line 773 "sgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;
#line 774 "sgsvj0.f"
					if (aaqq > aapp0) {
#line 774 "sgsvj0.f"
					    theta = -theta;
#line 774 "sgsvj0.f"
					}

#line 776 "sgsvj0.f"
					if (abs(theta) > bigtheta) {
#line 777 "sgsvj0.f"
					    t = .5 / theta;
#line 778 "sgsvj0.f"
					    fastr[2] = t * d__[p] / d__[q];
#line 779 "sgsvj0.f"
					    fastr[3] = -t * d__[q] / d__[p];
#line 780 "sgsvj0.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 782 "sgsvj0.f"
					    if (rsvec) {
#line 782 "sgsvj0.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 782 "sgsvj0.f"
					    }
/* Computing MAX */
#line 786 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 786 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 788 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 788 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 790 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 790 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);
#line 791 "sgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 795 "sgsvj0.f"
					    thsign = -d_sign(&c_b42, &aapq);
#line 796 "sgsvj0.f"
					    if (aaqq > aapp0) {
#line 796 "sgsvj0.f"
			  thsign = -thsign;
#line 796 "sgsvj0.f"
					    }
#line 797 "sgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 799 "sgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 800 "sgsvj0.f"
					    sn = t * cs;
/* Computing MAX */
#line 801 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 801 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 802 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 802 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 804 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 804 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 807 "sgsvj0.f"
					    apoaq = d__[p] / d__[q];
#line 808 "sgsvj0.f"
					    aqoap = d__[q] / d__[p];
#line 809 "sgsvj0.f"
					    if (d__[p] >= 1.) {

#line 811 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 812 "sgsvj0.f"
			      fastr[2] = t * apoaq;
#line 813 "sgsvj0.f"
			      fastr[3] = -t * aqoap;
#line 814 "sgsvj0.f"
			      d__[p] *= cs;
#line 815 "sgsvj0.f"
			      d__[q] *= cs;
#line 816 "sgsvj0.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 819 "sgsvj0.f"
			      if (rsvec) {
#line 819 "sgsvj0.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 819 "sgsvj0.f"
			      }
#line 822 "sgsvj0.f"
			  } else {
#line 823 "sgsvj0.f"
			      d__1 = -t * aqoap;
#line 823 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 826 "sgsvj0.f"
			      d__1 = cs * sn * apoaq;
#line 826 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 829 "sgsvj0.f"
			      if (rsvec) {
#line 830 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 830 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 833 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 833 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 837 "sgsvj0.f"
			      }
#line 838 "sgsvj0.f"
			      d__[p] *= cs;
#line 839 "sgsvj0.f"
			      d__[q] /= cs;
#line 840 "sgsvj0.f"
			  }
#line 841 "sgsvj0.f"
					    } else {
#line 842 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 843 "sgsvj0.f"
			      d__1 = t * apoaq;
#line 843 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 846 "sgsvj0.f"
			      d__1 = -cs * sn * aqoap;
#line 846 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 849 "sgsvj0.f"
			      if (rsvec) {
#line 850 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 850 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 853 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 853 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 857 "sgsvj0.f"
			      }
#line 858 "sgsvj0.f"
			      d__[p] /= cs;
#line 859 "sgsvj0.f"
			      d__[q] *= cs;
#line 860 "sgsvj0.f"
			  } else {
#line 861 "sgsvj0.f"
			      if (d__[p] >= d__[q]) {
#line 862 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 862 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 865 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 865 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 868 "sgsvj0.f"
				  d__[p] *= cs;
#line 869 "sgsvj0.f"
				  d__[q] /= cs;
#line 870 "sgsvj0.f"
				  if (rsvec) {
#line 871 "sgsvj0.f"
				      d__1 = -t * aqoap;
#line 871 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 875 "sgsvj0.f"
				      d__1 = cs * sn * apoaq;
#line 875 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 879 "sgsvj0.f"
				  }
#line 880 "sgsvj0.f"
			      } else {
#line 881 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 881 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 884 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 884 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 888 "sgsvj0.f"
				  d__[p] /= cs;
#line 889 "sgsvj0.f"
				  d__[q] *= cs;
#line 890 "sgsvj0.f"
				  if (rsvec) {
#line 891 "sgsvj0.f"
				      d__1 = t * apoaq;
#line 891 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 894 "sgsvj0.f"
				      d__1 = -cs * sn * aqoap;
#line 894 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 898 "sgsvj0.f"
				  }
#line 899 "sgsvj0.f"
			      }
#line 900 "sgsvj0.f"
			  }
#line 901 "sgsvj0.f"
					    }
#line 902 "sgsvj0.f"
					}

#line 904 "sgsvj0.f"
				    } else {
#line 905 "sgsvj0.f"
					if (aapp > aaqq) {
#line 906 "sgsvj0.f"
					    scopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 908 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 910 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 913 "sgsvj0.f"
					    temp1 = -aapq * d__[p] / d__[q];
#line 914 "sgsvj0.f"
					    saxpy_(m, &temp1, &work[1], &c__1,
						     &a[q * a_dim1 + 1], &
						    c__1);
#line 916 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &c_b42,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 919 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 919 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 921 "sgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 922 "sgsvj0.f"
					} else {
#line 923 "sgsvj0.f"
					    scopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 925 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 927 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 930 "sgsvj0.f"
					    temp1 = -aapq * d__[q] / d__[p];
#line 931 "sgsvj0.f"
					    saxpy_(m, &temp1, &work[1], &c__1,
						     &a[p * a_dim1 + 1], &
						    c__1);
#line 933 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &c_b42,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 936 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 936 "sgsvj0.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 938 "sgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 939 "sgsvj0.f"
					}
#line 940 "sgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q) */
/*           .. recompute SVA(q) */
/* Computing 2nd power */
#line 945 "sgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 945 "sgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 947 "sgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 949 "sgsvj0.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * d__[q];
#line 951 "sgsvj0.f"
					} else {
#line 952 "sgsvj0.f"
					    t = 0.;
#line 953 "sgsvj0.f"
					    aaqq = 1.;
#line 954 "sgsvj0.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 956 "sgsvj0.f"
					    sva[q] = t * sqrt(aaqq) * d__[q];
#line 957 "sgsvj0.f"
					}
#line 958 "sgsvj0.f"
				    }
/* Computing 2nd power */
#line 959 "sgsvj0.f"
				    d__1 = aapp / aapp0;
#line 959 "sgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 960 "sgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 962 "sgsvj0.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * d__[p];
#line 964 "sgsvj0.f"
					} else {
#line 965 "sgsvj0.f"
					    t = 0.;
#line 966 "sgsvj0.f"
					    aapp = 1.;
#line 967 "sgsvj0.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 969 "sgsvj0.f"
					    aapp = t * sqrt(aapp) * d__[p];
#line 970 "sgsvj0.f"
					}
#line 971 "sgsvj0.f"
					sva[p] = aapp;
#line 972 "sgsvj0.f"
				    }
/*              end of OK rotation */
#line 974 "sgsvj0.f"
				} else {
#line 975 "sgsvj0.f"
				    ++notrot;
#line 976 "sgsvj0.f"
				    ++pskipped;
#line 977 "sgsvj0.f"
				    ++ijblsk;
#line 978 "sgsvj0.f"
				}
#line 979 "sgsvj0.f"
			    } else {
#line 980 "sgsvj0.f"
				++notrot;
#line 981 "sgsvj0.f"
				++pskipped;
#line 982 "sgsvj0.f"
				++ijblsk;
#line 983 "sgsvj0.f"
			    }

#line 985 "sgsvj0.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 987 "sgsvj0.f"
				sva[p] = aapp;
#line 988 "sgsvj0.f"
				notrot = 0;
#line 989 "sgsvj0.f"
				goto L2011;
#line 990 "sgsvj0.f"
			    }
#line 991 "sgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 993 "sgsvj0.f"
				aapp = -aapp;
#line 994 "sgsvj0.f"
				notrot = 0;
#line 995 "sgsvj0.f"
				goto L2203;
#line 996 "sgsvj0.f"
			    }

#line 998 "sgsvj0.f"
/* L2200: */
#line 998 "sgsvj0.f"
			}
/*        end of the q-loop */
#line 1000 "sgsvj0.f"
L2203:

#line 1002 "sgsvj0.f"
			sva[p] = aapp;

#line 1004 "sgsvj0.f"
		    } else {
#line 1005 "sgsvj0.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1005 "sgsvj0.f"
			    i__5 = jgl + kbl - 1;
#line 1005 "sgsvj0.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 1005 "sgsvj0.f"
			}
#line 1007 "sgsvj0.f"
			if (aapp < 0.) {
#line 1007 "sgsvj0.f"
			    notrot = 0;
#line 1007 "sgsvj0.f"
			}
#line 1008 "sgsvj0.f"
		    }
#line 1010 "sgsvj0.f"
/* L2100: */
#line 1010 "sgsvj0.f"
		}
/*     end of the p-loop */
#line 1012 "sgsvj0.f"
/* L2010: */
#line 1012 "sgsvj0.f"
	    }
/*     end of the jbc-loop */
#line 1014 "sgsvj0.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1016 "sgsvj0.f"
	    i__4 = igl + kbl - 1;
#line 1016 "sgsvj0.f"
	    i__3 = min(i__4,*n);
#line 1016 "sgsvj0.f"
	    for (p = igl; p <= i__3; ++p) {
#line 1017 "sgsvj0.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1018 "sgsvj0.f"
/* L2012: */
#line 1018 "sgsvj0.f"
	    }

#line 1020 "sgsvj0.f"
/* L2000: */
#line 1020 "sgsvj0.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1024 "sgsvj0.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1026 "sgsvj0.f"
	    sva[*n] = snrm2_(m, &a[*n * a_dim1 + 1], &c__1) * d__[*n];
#line 1027 "sgsvj0.f"
	} else {
#line 1028 "sgsvj0.f"
	    t = 0.;
#line 1029 "sgsvj0.f"
	    aapp = 1.;
#line 1030 "sgsvj0.f"
	    slassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1031 "sgsvj0.f"
	    sva[*n] = t * sqrt(aapp) * d__[*n];
#line 1032 "sgsvj0.f"
	}

/*     Additional steering devices */

#line 1036 "sgsvj0.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1036 "sgsvj0.f"
	    swband = i__;
#line 1036 "sgsvj0.f"
	}

#line 1039 "sgsvj0.f"
	if (i__ > swband + 1 && mxaapq < (doublereal) (*n) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 1041 "sgsvj0.f"
	    goto L1994;
#line 1042 "sgsvj0.f"
	}

#line 1044 "sgsvj0.f"
	if (notrot >= emptsw) {
#line 1044 "sgsvj0.f"
	    goto L1994;
#line 1044 "sgsvj0.f"
	}
#line 1046 "sgsvj0.f"
/* L1993: */
#line 1046 "sgsvj0.f"
    }
/*     end i=1:NSWEEP loop */
/* #:) Reaching this point means that the procedure has comleted the given */
/*     number of iterations. */
#line 1050 "sgsvj0.f"
    *info = *nsweep - 1;
#line 1051 "sgsvj0.f"
    goto L1995;
#line 1052 "sgsvj0.f"
L1994:
/* #:) Reaching this point means that during the i-th sweep all pivots were */
/*     below the given tolerance, causing early exit. */

#line 1056 "sgsvj0.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1058 "sgsvj0.f"
L1995:

/*     Sort the vector D. */
#line 1061 "sgsvj0.f"
    i__1 = *n - 1;
#line 1061 "sgsvj0.f"
    for (p = 1; p <= i__1; ++p) {
#line 1062 "sgsvj0.f"
	i__2 = *n - p + 1;
#line 1062 "sgsvj0.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1063 "sgsvj0.f"
	if (p != q) {
#line 1064 "sgsvj0.f"
	    temp1 = sva[p];
#line 1065 "sgsvj0.f"
	    sva[p] = sva[q];
#line 1066 "sgsvj0.f"
	    sva[q] = temp1;
#line 1067 "sgsvj0.f"
	    temp1 = d__[p];
#line 1068 "sgsvj0.f"
	    d__[p] = d__[q];
#line 1069 "sgsvj0.f"
	    d__[q] = temp1;
#line 1070 "sgsvj0.f"
	    sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1071 "sgsvj0.f"
	    if (rsvec) {
#line 1071 "sgsvj0.f"
		sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1071 "sgsvj0.f"
	    }
#line 1072 "sgsvj0.f"
	}
#line 1073 "sgsvj0.f"
/* L5991: */
#line 1073 "sgsvj0.f"
    }

#line 1075 "sgsvj0.f"
    return 0;
/*     .. */
/*     .. END OF SGSVJ0 */
/*     .. */
} /* sgsvj0_ */

