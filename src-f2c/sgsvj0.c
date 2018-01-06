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
/* >          WORK is REAL array, dimension LWORK. */
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

#line 271 "sgsvj0.f"
    /* Parameter adjustments */
#line 271 "sgsvj0.f"
    --sva;
#line 271 "sgsvj0.f"
    --d__;
#line 271 "sgsvj0.f"
    a_dim1 = *lda;
#line 271 "sgsvj0.f"
    a_offset = 1 + a_dim1;
#line 271 "sgsvj0.f"
    a -= a_offset;
#line 271 "sgsvj0.f"
    v_dim1 = *ldv;
#line 271 "sgsvj0.f"
    v_offset = 1 + v_dim1;
#line 271 "sgsvj0.f"
    v -= v_offset;
#line 271 "sgsvj0.f"
    --work;
#line 271 "sgsvj0.f"

#line 271 "sgsvj0.f"
    /* Function Body */
#line 271 "sgsvj0.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 272 "sgsvj0.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 273 "sgsvj0.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 274 "sgsvj0.f"
	*info = -1;
#line 275 "sgsvj0.f"
    } else if (*m < 0) {
#line 276 "sgsvj0.f"
	*info = -2;
#line 277 "sgsvj0.f"
    } else if (*n < 0 || *n > *m) {
#line 278 "sgsvj0.f"
	*info = -3;
#line 279 "sgsvj0.f"
    } else if (*lda < *m) {
#line 280 "sgsvj0.f"
	*info = -5;
#line 281 "sgsvj0.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 282 "sgsvj0.f"
	*info = -8;
#line 283 "sgsvj0.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 285 "sgsvj0.f"
	*info = -10;
#line 286 "sgsvj0.f"
    } else if (*tol <= *eps) {
#line 287 "sgsvj0.f"
	*info = -13;
#line 288 "sgsvj0.f"
    } else if (*nsweep < 0) {
#line 289 "sgsvj0.f"
	*info = -14;
#line 290 "sgsvj0.f"
    } else if (*lwork < *m) {
#line 291 "sgsvj0.f"
	*info = -16;
#line 292 "sgsvj0.f"
    } else {
#line 293 "sgsvj0.f"
	*info = 0;
#line 294 "sgsvj0.f"
    }

/*     #:( */
#line 297 "sgsvj0.f"
    if (*info != 0) {
#line 298 "sgsvj0.f"
	i__1 = -(*info);
#line 298 "sgsvj0.f"
	xerbla_("SGSVJ0", &i__1, (ftnlen)6);
#line 299 "sgsvj0.f"
	return 0;
#line 300 "sgsvj0.f"
    }

#line 302 "sgsvj0.f"
    if (rsvec) {
#line 303 "sgsvj0.f"
	mvl = *n;
#line 304 "sgsvj0.f"
    } else if (applv) {
#line 305 "sgsvj0.f"
	mvl = *mv;
#line 306 "sgsvj0.f"
    }
#line 307 "sgsvj0.f"
    rsvec = rsvec || applv;
#line 309 "sgsvj0.f"
    rooteps = sqrt(*eps);
#line 310 "sgsvj0.f"
    rootsfmin = sqrt(*sfmin);
#line 311 "sgsvj0.f"
    small = *sfmin / *eps;
#line 312 "sgsvj0.f"
    big = 1. / *sfmin;
#line 313 "sgsvj0.f"
    rootbig = 1. / rootsfmin;
#line 314 "sgsvj0.f"
    bigtheta = 1. / rooteps;
#line 315 "sgsvj0.f"
    roottol = sqrt(*tol);

/*     .. Row-cyclic Jacobi SVD algorithm with column pivoting .. */

#line 319 "sgsvj0.f"
    emptsw = *n * (*n - 1) / 2;
#line 320 "sgsvj0.f"
    notrot = 0;
#line 321 "sgsvj0.f"
    fastr[0] = 0.;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 326 "sgsvj0.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter. It is meaningful and effective */
/*     if SGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure */
/*     ...... */
#line 332 "sgsvj0.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 338 "sgsvj0.f"
    nbl = *n / kbl;
#line 339 "sgsvj0.f"
    if (nbl * kbl != *n) {
#line 339 "sgsvj0.f"
	++nbl;
#line 339 "sgsvj0.f"
    }
/* Computing 2nd power */
#line 341 "sgsvj0.f"
    i__1 = kbl;
#line 341 "sgsvj0.f"
    blskip = i__1 * i__1 + 1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
#line 344 "sgsvj0.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */
#line 347 "sgsvj0.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */
#line 349 "sgsvj0.f"
    swband = 0;
#line 350 "sgsvj0.f"
    pskipped = 0;

#line 352 "sgsvj0.f"
    i__1 = *nsweep;
#line 352 "sgsvj0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     .. go go go ... */

#line 355 "sgsvj0.f"
	mxaapq = 0.;
#line 356 "sgsvj0.f"
	mxsinj = 0.;
#line 357 "sgsvj0.f"
	iswrot = 0;

#line 359 "sgsvj0.f"
	notrot = 0;
#line 360 "sgsvj0.f"
	pskipped = 0;

#line 362 "sgsvj0.f"
	i__2 = nbl;
#line 362 "sgsvj0.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {
#line 364 "sgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 366 "sgsvj0.f"
	    i__4 = lkahead, i__5 = nbl - ibr;
#line 366 "sgsvj0.f"
	    i__3 = min(i__4,i__5);
#line 366 "sgsvj0.f"
	    for (ir1 = 0; ir1 <= i__3; ++ir1) {

#line 368 "sgsvj0.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 370 "sgsvj0.f"
		i__5 = igl + kbl - 1, i__6 = *n - 1;
#line 370 "sgsvj0.f"
		i__4 = min(i__5,i__6);
#line 370 "sgsvj0.f"
		for (p = igl; p <= i__4; ++p) {
/*     .. de Rijk's pivoting */
#line 373 "sgsvj0.f"
		    i__5 = *n - p + 1;
#line 373 "sgsvj0.f"
		    q = isamax_(&i__5, &sva[p], &c__1) + p - 1;
#line 374 "sgsvj0.f"
		    if (p != q) {
#line 375 "sgsvj0.f"
			sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 376 "sgsvj0.f"
			if (rsvec) {
#line 376 "sgsvj0.f"
			    sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 376 "sgsvj0.f"
			}
#line 378 "sgsvj0.f"
			temp1 = sva[p];
#line 379 "sgsvj0.f"
			sva[p] = sva[q];
#line 380 "sgsvj0.f"
			sva[q] = temp1;
#line 381 "sgsvj0.f"
			temp1 = d__[p];
#line 382 "sgsvj0.f"
			d__[p] = d__[q];
#line 383 "sgsvj0.f"
			d__[q] = temp1;
#line 384 "sgsvj0.f"
		    }

#line 386 "sgsvj0.f"
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

#line 400 "sgsvj0.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 402 "sgsvj0.f"
			    sva[p] = snrm2_(m, &a[p * a_dim1 + 1], &c__1) * 
				    d__[p];
#line 403 "sgsvj0.f"
			} else {
#line 404 "sgsvj0.f"
			    temp1 = 0.;
#line 405 "sgsvj0.f"
			    aapp = 1.;
#line 406 "sgsvj0.f"
			    slassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 407 "sgsvj0.f"
			    sva[p] = temp1 * sqrt(aapp) * d__[p];
#line 408 "sgsvj0.f"
			}
#line 409 "sgsvj0.f"
			aapp = sva[p];
#line 410 "sgsvj0.f"
		    } else {
#line 411 "sgsvj0.f"
			aapp = sva[p];
#line 412 "sgsvj0.f"
		    }

#line 415 "sgsvj0.f"
		    if (aapp > 0.) {

#line 417 "sgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 419 "sgsvj0.f"
			i__6 = igl + kbl - 1;
#line 419 "sgsvj0.f"
			i__5 = min(i__6,*n);
#line 419 "sgsvj0.f"
			for (q = p + 1; q <= i__5; ++q) {

#line 421 "sgsvj0.f"
			    aaqq = sva[q];
#line 423 "sgsvj0.f"
			    if (aaqq > 0.) {

#line 425 "sgsvj0.f"
				aapp0 = aapp;
#line 426 "sgsvj0.f"
				if (aaqq >= 1.) {
#line 427 "sgsvj0.f"
				    rotok = small * aapp <= aaqq;
#line 428 "sgsvj0.f"
				    if (aapp < big / aaqq) {
#line 429 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 432 "sgsvj0.f"
				    } else {
#line 433 "sgsvj0.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 434 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						d__[p], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 436 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[q 
						* a_dim1 + 1], &c__1) * d__[q]
						 / aaqq;
#line 438 "sgsvj0.f"
				    }
#line 439 "sgsvj0.f"
				} else {
#line 440 "sgsvj0.f"
				    rotok = aapp <= aaqq / small;
#line 441 "sgsvj0.f"
				    if (aapp > small / aaqq) {
#line 442 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 445 "sgsvj0.f"
				    } else {
#line 446 "sgsvj0.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 447 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						d__[q], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 449 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[p 
						* a_dim1 + 1], &c__1) * d__[p]
						 / aapp;
#line 451 "sgsvj0.f"
				    }
#line 452 "sgsvj0.f"
				}

/* Computing MAX */
#line 454 "sgsvj0.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 454 "sgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 458 "sgsvj0.f"
				if (abs(aapq) > *tol) {

/*           .. rotate */
/*           ROTATED = ROTATED + ONE */

#line 463 "sgsvj0.f"
				    if (ir1 == 0) {
#line 464 "sgsvj0.f"
					notrot = 0;
#line 465 "sgsvj0.f"
					pskipped = 0;
#line 466 "sgsvj0.f"
					++iswrot;
#line 467 "sgsvj0.f"
				    }

#line 469 "sgsvj0.f"
				    if (rotok) {

#line 471 "sgsvj0.f"
					aqoap = aaqq / aapp;
#line 472 "sgsvj0.f"
					apoaq = aapp / aaqq;
#line 473 "sgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;

#line 475 "sgsvj0.f"
					if (abs(theta) > bigtheta) {

#line 477 "sgsvj0.f"
					    t = .5 / theta;
#line 478 "sgsvj0.f"
					    fastr[2] = t * d__[p] / d__[q];
#line 479 "sgsvj0.f"
					    fastr[3] = -t * d__[q] / d__[p];
#line 480 "sgsvj0.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 482 "sgsvj0.f"
					    if (rsvec) {
#line 482 "sgsvj0.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 482 "sgsvj0.f"
					    }
/* Computing MAX */
#line 486 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 486 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 488 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 488 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 490 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 490 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);

#line 492 "sgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 496 "sgsvj0.f"
					    thsign = -d_sign(&c_b42, &aapq);
#line 497 "sgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 499 "sgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 500 "sgsvj0.f"
					    sn = t * cs;

/* Computing MAX */
#line 502 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 502 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 503 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 503 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 505 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 505 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 508 "sgsvj0.f"
					    apoaq = d__[p] / d__[q];
#line 509 "sgsvj0.f"
					    aqoap = d__[q] / d__[p];
#line 510 "sgsvj0.f"
					    if (d__[p] >= 1.) {
#line 511 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 512 "sgsvj0.f"
			      fastr[2] = t * apoaq;
#line 513 "sgsvj0.f"
			      fastr[3] = -t * aqoap;
#line 514 "sgsvj0.f"
			      d__[p] *= cs;
#line 515 "sgsvj0.f"
			      d__[q] *= cs;
#line 516 "sgsvj0.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 519 "sgsvj0.f"
			      if (rsvec) {
#line 519 "sgsvj0.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 519 "sgsvj0.f"
			      }
#line 522 "sgsvj0.f"
			  } else {
#line 523 "sgsvj0.f"
			      d__1 = -t * aqoap;
#line 523 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 526 "sgsvj0.f"
			      d__1 = cs * sn * apoaq;
#line 526 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 529 "sgsvj0.f"
			      d__[p] *= cs;
#line 530 "sgsvj0.f"
			      d__[q] /= cs;
#line 531 "sgsvj0.f"
			      if (rsvec) {
#line 532 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 532 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 535 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 535 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 539 "sgsvj0.f"
			      }
#line 540 "sgsvj0.f"
			  }
#line 541 "sgsvj0.f"
					    } else {
#line 542 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 543 "sgsvj0.f"
			      d__1 = t * apoaq;
#line 543 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 546 "sgsvj0.f"
			      d__1 = -cs * sn * aqoap;
#line 546 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 549 "sgsvj0.f"
			      d__[p] /= cs;
#line 550 "sgsvj0.f"
			      d__[q] *= cs;
#line 551 "sgsvj0.f"
			      if (rsvec) {
#line 552 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 552 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 555 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 555 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 559 "sgsvj0.f"
			      }
#line 560 "sgsvj0.f"
			  } else {
#line 561 "sgsvj0.f"
			      if (d__[p] >= d__[q]) {
#line 562 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 562 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 565 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 565 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 568 "sgsvj0.f"
				  d__[p] *= cs;
#line 569 "sgsvj0.f"
				  d__[q] /= cs;
#line 570 "sgsvj0.f"
				  if (rsvec) {
#line 571 "sgsvj0.f"
				      d__1 = -t * aqoap;
#line 571 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 575 "sgsvj0.f"
				      d__1 = cs * sn * apoaq;
#line 575 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 579 "sgsvj0.f"
				  }
#line 580 "sgsvj0.f"
			      } else {
#line 581 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 581 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 584 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 584 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 588 "sgsvj0.f"
				  d__[p] /= cs;
#line 589 "sgsvj0.f"
				  d__[q] *= cs;
#line 590 "sgsvj0.f"
				  if (rsvec) {
#line 591 "sgsvj0.f"
				      d__1 = t * apoaq;
#line 591 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 594 "sgsvj0.f"
				      d__1 = -cs * sn * aqoap;
#line 594 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 598 "sgsvj0.f"
				  }
#line 599 "sgsvj0.f"
			      }
#line 600 "sgsvj0.f"
			  }
#line 601 "sgsvj0.f"
					    }
#line 602 "sgsvj0.f"
					}

#line 604 "sgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 606 "sgsvj0.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 607 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						c_b42, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 609 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						c_b42, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 611 "sgsvj0.f"
					temp1 = -aapq * d__[p] / d__[q];
#line 612 "sgsvj0.f"
					saxpy_(m, &temp1, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 614 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &c_b42, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 616 "sgsvj0.f"
					d__1 = 0., d__2 = 1. - aapq * aapq;
#line 616 "sgsvj0.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 618 "sgsvj0.f"
					mxsinj = max(mxsinj,*sfmin);
#line 619 "sgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */
/* Computing 2nd power */
#line 624 "sgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 624 "sgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 626 "sgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 628 "sgsvj0.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * d__[q];
#line 630 "sgsvj0.f"
					} else {
#line 631 "sgsvj0.f"
					    t = 0.;
#line 632 "sgsvj0.f"
					    aaqq = 1.;
#line 633 "sgsvj0.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 635 "sgsvj0.f"
					    sva[q] = t * sqrt(aaqq) * d__[q];
#line 636 "sgsvj0.f"
					}
#line 637 "sgsvj0.f"
				    }
#line 638 "sgsvj0.f"
				    if (aapp / aapp0 <= rooteps) {
#line 639 "sgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 641 "sgsvj0.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * d__[p];
#line 643 "sgsvj0.f"
					} else {
#line 644 "sgsvj0.f"
					    t = 0.;
#line 645 "sgsvj0.f"
					    aapp = 1.;
#line 646 "sgsvj0.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 648 "sgsvj0.f"
					    aapp = t * sqrt(aapp) * d__[p];
#line 649 "sgsvj0.f"
					}
#line 650 "sgsvj0.f"
					sva[p] = aapp;
#line 651 "sgsvj0.f"
				    }

#line 653 "sgsvj0.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 655 "sgsvj0.f"
				    if (ir1 == 0) {
#line 655 "sgsvj0.f"
					++notrot;
#line 655 "sgsvj0.f"
				    }
#line 656 "sgsvj0.f"
				    ++pskipped;
#line 657 "sgsvj0.f"
				}
#line 658 "sgsvj0.f"
			    } else {
/*        A(:,q) is zero column */
#line 660 "sgsvj0.f"
				if (ir1 == 0) {
#line 660 "sgsvj0.f"
				    ++notrot;
#line 660 "sgsvj0.f"
				}
#line 661 "sgsvj0.f"
				++pskipped;
#line 662 "sgsvj0.f"
			    }

#line 664 "sgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 666 "sgsvj0.f"
				if (ir1 == 0) {
#line 666 "sgsvj0.f"
				    aapp = -aapp;
#line 666 "sgsvj0.f"
				}
#line 667 "sgsvj0.f"
				notrot = 0;
#line 668 "sgsvj0.f"
				goto L2103;
#line 669 "sgsvj0.f"
			    }

#line 671 "sgsvj0.f"
/* L2002: */
#line 671 "sgsvj0.f"
			}
/*     END q-LOOP */

#line 674 "sgsvj0.f"
L2103:
/*     bailed out of q-loop */
#line 677 "sgsvj0.f"
			sva[p] = aapp;
#line 679 "sgsvj0.f"
		    } else {
#line 680 "sgsvj0.f"
			sva[p] = aapp;
#line 681 "sgsvj0.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 681 "sgsvj0.f"
			    i__5 = igl + kbl - 1;
#line 681 "sgsvj0.f"
			    notrot = notrot + min(i__5,*n) - p;
#line 681 "sgsvj0.f"
			}
#line 683 "sgsvj0.f"
		    }

#line 685 "sgsvj0.f"
/* L2001: */
#line 685 "sgsvj0.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 688 "sgsvj0.f"
/* L1002: */
#line 688 "sgsvj0.f"
	    }
/*     end of ir1-loop */

/* ........................................................ */
/* ... go to the off diagonal blocks */

#line 694 "sgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

#line 696 "sgsvj0.f"
	    i__3 = nbl;
#line 696 "sgsvj0.f"
	    for (jbc = ibr + 1; jbc <= i__3; ++jbc) {

#line 698 "sgsvj0.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 702 "sgsvj0.f"
		ijblsk = 0;
/* Computing MIN */
#line 703 "sgsvj0.f"
		i__5 = igl + kbl - 1;
#line 703 "sgsvj0.f"
		i__4 = min(i__5,*n);
#line 703 "sgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

#line 705 "sgsvj0.f"
		    aapp = sva[p];

#line 707 "sgsvj0.f"
		    if (aapp > 0.) {

#line 709 "sgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 711 "sgsvj0.f"
			i__6 = jgl + kbl - 1;
#line 711 "sgsvj0.f"
			i__5 = min(i__6,*n);
#line 711 "sgsvj0.f"
			for (q = jgl; q <= i__5; ++q) {

#line 713 "sgsvj0.f"
			    aaqq = sva[q];

#line 715 "sgsvj0.f"
			    if (aaqq > 0.) {
#line 716 "sgsvj0.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        .. Safe Gram matrix computation .. */

#line 722 "sgsvj0.f"
				if (aaqq >= 1.) {
#line 723 "sgsvj0.f"
				    if (aapp >= aaqq) {
#line 724 "sgsvj0.f"
					rotok = small * aapp <= aaqq;
#line 725 "sgsvj0.f"
				    } else {
#line 726 "sgsvj0.f"
					rotok = small * aaqq <= aapp;
#line 727 "sgsvj0.f"
				    }
#line 728 "sgsvj0.f"
				    if (aapp < big / aaqq) {
#line 729 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 732 "sgsvj0.f"
				    } else {
#line 733 "sgsvj0.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 734 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						d__[p], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 736 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[q 
						* a_dim1 + 1], &c__1) * d__[q]
						 / aaqq;
#line 738 "sgsvj0.f"
				    }
#line 739 "sgsvj0.f"
				} else {
#line 740 "sgsvj0.f"
				    if (aapp >= aaqq) {
#line 741 "sgsvj0.f"
					rotok = aapp <= aaqq / small;
#line 742 "sgsvj0.f"
				    } else {
#line 743 "sgsvj0.f"
					rotok = aaqq <= aapp / small;
#line 744 "sgsvj0.f"
				    }
#line 745 "sgsvj0.f"
				    if (aapp > small / aaqq) {
#line 746 "sgsvj0.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * d__[p] * d__[q] / 
						aaqq / aapp;
#line 749 "sgsvj0.f"
				    } else {
#line 750 "sgsvj0.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 751 "sgsvj0.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						d__[q], m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 753 "sgsvj0.f"
					aapq = sdot_(m, &work[1], &c__1, &a[p 
						* a_dim1 + 1], &c__1) * d__[p]
						 / aapp;
#line 755 "sgsvj0.f"
				    }
#line 756 "sgsvj0.f"
				}

/* Computing MAX */
#line 758 "sgsvj0.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 758 "sgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 762 "sgsvj0.f"
				if (abs(aapq) > *tol) {
#line 763 "sgsvj0.f"
				    notrot = 0;
/*           ROTATED  = ROTATED + 1 */
#line 765 "sgsvj0.f"
				    pskipped = 0;
#line 766 "sgsvj0.f"
				    ++iswrot;

#line 768 "sgsvj0.f"
				    if (rotok) {

#line 770 "sgsvj0.f"
					aqoap = aaqq / aapp;
#line 771 "sgsvj0.f"
					apoaq = aapp / aaqq;
#line 772 "sgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;
#line 773 "sgsvj0.f"
					if (aaqq > aapp0) {
#line 773 "sgsvj0.f"
					    theta = -theta;
#line 773 "sgsvj0.f"
					}

#line 775 "sgsvj0.f"
					if (abs(theta) > bigtheta) {
#line 776 "sgsvj0.f"
					    t = .5 / theta;
#line 777 "sgsvj0.f"
					    fastr[2] = t * d__[p] / d__[q];
#line 778 "sgsvj0.f"
					    fastr[3] = -t * d__[q] / d__[p];
#line 779 "sgsvj0.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 781 "sgsvj0.f"
					    if (rsvec) {
#line 781 "sgsvj0.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 781 "sgsvj0.f"
					    }
/* Computing MAX */
#line 785 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 785 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 787 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 787 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 789 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 789 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);
#line 790 "sgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 794 "sgsvj0.f"
					    thsign = -d_sign(&c_b42, &aapq);
#line 795 "sgsvj0.f"
					    if (aaqq > aapp0) {
#line 795 "sgsvj0.f"
			  thsign = -thsign;
#line 795 "sgsvj0.f"
					    }
#line 796 "sgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 798 "sgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 799 "sgsvj0.f"
					    sn = t * cs;
/* Computing MAX */
#line 800 "sgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 800 "sgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 801 "sgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 801 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 803 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 803 "sgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 806 "sgsvj0.f"
					    apoaq = d__[p] / d__[q];
#line 807 "sgsvj0.f"
					    aqoap = d__[q] / d__[p];
#line 808 "sgsvj0.f"
					    if (d__[p] >= 1.) {

#line 810 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 811 "sgsvj0.f"
			      fastr[2] = t * apoaq;
#line 812 "sgsvj0.f"
			      fastr[3] = -t * aqoap;
#line 813 "sgsvj0.f"
			      d__[p] *= cs;
#line 814 "sgsvj0.f"
			      d__[q] *= cs;
#line 815 "sgsvj0.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 818 "sgsvj0.f"
			      if (rsvec) {
#line 818 "sgsvj0.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 818 "sgsvj0.f"
			      }
#line 821 "sgsvj0.f"
			  } else {
#line 822 "sgsvj0.f"
			      d__1 = -t * aqoap;
#line 822 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 825 "sgsvj0.f"
			      d__1 = cs * sn * apoaq;
#line 825 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 828 "sgsvj0.f"
			      if (rsvec) {
#line 829 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 829 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 832 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 832 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 836 "sgsvj0.f"
			      }
#line 837 "sgsvj0.f"
			      d__[p] *= cs;
#line 838 "sgsvj0.f"
			      d__[q] /= cs;
#line 839 "sgsvj0.f"
			  }
#line 840 "sgsvj0.f"
					    } else {
#line 841 "sgsvj0.f"
			  if (d__[q] >= 1.) {
#line 842 "sgsvj0.f"
			      d__1 = t * apoaq;
#line 842 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 845 "sgsvj0.f"
			      d__1 = -cs * sn * aqoap;
#line 845 "sgsvj0.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 848 "sgsvj0.f"
			      if (rsvec) {
#line 849 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 849 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 852 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 852 "sgsvj0.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 856 "sgsvj0.f"
			      }
#line 857 "sgsvj0.f"
			      d__[p] /= cs;
#line 858 "sgsvj0.f"
			      d__[q] *= cs;
#line 859 "sgsvj0.f"
			  } else {
#line 860 "sgsvj0.f"
			      if (d__[p] >= d__[q]) {
#line 861 "sgsvj0.f"
				  d__1 = -t * aqoap;
#line 861 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 864 "sgsvj0.f"
				  d__1 = cs * sn * apoaq;
#line 864 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 867 "sgsvj0.f"
				  d__[p] *= cs;
#line 868 "sgsvj0.f"
				  d__[q] /= cs;
#line 869 "sgsvj0.f"
				  if (rsvec) {
#line 870 "sgsvj0.f"
				      d__1 = -t * aqoap;
#line 870 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 874 "sgsvj0.f"
				      d__1 = cs * sn * apoaq;
#line 874 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 878 "sgsvj0.f"
				  }
#line 879 "sgsvj0.f"
			      } else {
#line 880 "sgsvj0.f"
				  d__1 = t * apoaq;
#line 880 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 883 "sgsvj0.f"
				  d__1 = -cs * sn * aqoap;
#line 883 "sgsvj0.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 887 "sgsvj0.f"
				  d__[p] /= cs;
#line 888 "sgsvj0.f"
				  d__[q] *= cs;
#line 889 "sgsvj0.f"
				  if (rsvec) {
#line 890 "sgsvj0.f"
				      d__1 = t * apoaq;
#line 890 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 893 "sgsvj0.f"
				      d__1 = -cs * sn * aqoap;
#line 893 "sgsvj0.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 897 "sgsvj0.f"
				  }
#line 898 "sgsvj0.f"
			      }
#line 899 "sgsvj0.f"
			  }
#line 900 "sgsvj0.f"
					    }
#line 901 "sgsvj0.f"
					}

#line 903 "sgsvj0.f"
				    } else {
#line 904 "sgsvj0.f"
					if (aapp > aaqq) {
#line 905 "sgsvj0.f"
					    scopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 907 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 909 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 912 "sgsvj0.f"
					    temp1 = -aapq * d__[p] / d__[q];
#line 913 "sgsvj0.f"
					    saxpy_(m, &temp1, &work[1], &c__1,
						     &a[q * a_dim1 + 1], &
						    c__1);
#line 915 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &c_b42,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 918 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 918 "sgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 920 "sgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 921 "sgsvj0.f"
					} else {
#line 922 "sgsvj0.f"
					    scopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 924 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 926 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 929 "sgsvj0.f"
					    temp1 = -aapq * d__[q] / d__[p];
#line 930 "sgsvj0.f"
					    saxpy_(m, &temp1, &work[1], &c__1,
						     &a[p * a_dim1 + 1], &
						    c__1);
#line 932 "sgsvj0.f"
					    slascl_("G", &c__0, &c__0, &c_b42,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 935 "sgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 935 "sgsvj0.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 937 "sgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 938 "sgsvj0.f"
					}
#line 939 "sgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q) */
/*           .. recompute SVA(q) */
/* Computing 2nd power */
#line 944 "sgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 944 "sgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 946 "sgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 948 "sgsvj0.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * d__[q];
#line 950 "sgsvj0.f"
					} else {
#line 951 "sgsvj0.f"
					    t = 0.;
#line 952 "sgsvj0.f"
					    aaqq = 1.;
#line 953 "sgsvj0.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 955 "sgsvj0.f"
					    sva[q] = t * sqrt(aaqq) * d__[q];
#line 956 "sgsvj0.f"
					}
#line 957 "sgsvj0.f"
				    }
/* Computing 2nd power */
#line 958 "sgsvj0.f"
				    d__1 = aapp / aapp0;
#line 958 "sgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 959 "sgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 961 "sgsvj0.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * d__[p];
#line 963 "sgsvj0.f"
					} else {
#line 964 "sgsvj0.f"
					    t = 0.;
#line 965 "sgsvj0.f"
					    aapp = 1.;
#line 966 "sgsvj0.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 968 "sgsvj0.f"
					    aapp = t * sqrt(aapp) * d__[p];
#line 969 "sgsvj0.f"
					}
#line 970 "sgsvj0.f"
					sva[p] = aapp;
#line 971 "sgsvj0.f"
				    }
/*              end of OK rotation */
#line 973 "sgsvj0.f"
				} else {
#line 974 "sgsvj0.f"
				    ++notrot;
#line 975 "sgsvj0.f"
				    ++pskipped;
#line 976 "sgsvj0.f"
				    ++ijblsk;
#line 977 "sgsvj0.f"
				}
#line 978 "sgsvj0.f"
			    } else {
#line 979 "sgsvj0.f"
				++notrot;
#line 980 "sgsvj0.f"
				++pskipped;
#line 981 "sgsvj0.f"
				++ijblsk;
#line 982 "sgsvj0.f"
			    }

#line 984 "sgsvj0.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 986 "sgsvj0.f"
				sva[p] = aapp;
#line 987 "sgsvj0.f"
				notrot = 0;
#line 988 "sgsvj0.f"
				goto L2011;
#line 989 "sgsvj0.f"
			    }
#line 990 "sgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 992 "sgsvj0.f"
				aapp = -aapp;
#line 993 "sgsvj0.f"
				notrot = 0;
#line 994 "sgsvj0.f"
				goto L2203;
#line 995 "sgsvj0.f"
			    }

#line 997 "sgsvj0.f"
/* L2200: */
#line 997 "sgsvj0.f"
			}
/*        end of the q-loop */
#line 999 "sgsvj0.f"
L2203:

#line 1001 "sgsvj0.f"
			sva[p] = aapp;

#line 1003 "sgsvj0.f"
		    } else {
#line 1004 "sgsvj0.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1004 "sgsvj0.f"
			    i__5 = jgl + kbl - 1;
#line 1004 "sgsvj0.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 1004 "sgsvj0.f"
			}
#line 1006 "sgsvj0.f"
			if (aapp < 0.) {
#line 1006 "sgsvj0.f"
			    notrot = 0;
#line 1006 "sgsvj0.f"
			}
#line 1007 "sgsvj0.f"
		    }
#line 1009 "sgsvj0.f"
/* L2100: */
#line 1009 "sgsvj0.f"
		}
/*     end of the p-loop */
#line 1011 "sgsvj0.f"
/* L2010: */
#line 1011 "sgsvj0.f"
	    }
/*     end of the jbc-loop */
#line 1013 "sgsvj0.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1015 "sgsvj0.f"
	    i__4 = igl + kbl - 1;
#line 1015 "sgsvj0.f"
	    i__3 = min(i__4,*n);
#line 1015 "sgsvj0.f"
	    for (p = igl; p <= i__3; ++p) {
#line 1016 "sgsvj0.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1017 "sgsvj0.f"
/* L2012: */
#line 1017 "sgsvj0.f"
	    }

#line 1019 "sgsvj0.f"
/* L2000: */
#line 1019 "sgsvj0.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1023 "sgsvj0.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1025 "sgsvj0.f"
	    sva[*n] = snrm2_(m, &a[*n * a_dim1 + 1], &c__1) * d__[*n];
#line 1026 "sgsvj0.f"
	} else {
#line 1027 "sgsvj0.f"
	    t = 0.;
#line 1028 "sgsvj0.f"
	    aapp = 1.;
#line 1029 "sgsvj0.f"
	    slassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1030 "sgsvj0.f"
	    sva[*n] = t * sqrt(aapp) * d__[*n];
#line 1031 "sgsvj0.f"
	}

/*     Additional steering devices */

#line 1035 "sgsvj0.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1035 "sgsvj0.f"
	    swband = i__;
#line 1035 "sgsvj0.f"
	}

#line 1038 "sgsvj0.f"
	if (i__ > swband + 1 && mxaapq < (doublereal) (*n) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 1040 "sgsvj0.f"
	    goto L1994;
#line 1041 "sgsvj0.f"
	}

#line 1043 "sgsvj0.f"
	if (notrot >= emptsw) {
#line 1043 "sgsvj0.f"
	    goto L1994;
#line 1043 "sgsvj0.f"
	}
#line 1045 "sgsvj0.f"
/* L1993: */
#line 1045 "sgsvj0.f"
    }
/*     end i=1:NSWEEP loop */
/* #:) Reaching this point means that the procedure has comleted the given */
/*     number of iterations. */
#line 1049 "sgsvj0.f"
    *info = *nsweep - 1;
#line 1050 "sgsvj0.f"
    goto L1995;
#line 1051 "sgsvj0.f"
L1994:
/* #:) Reaching this point means that during the i-th sweep all pivots were */
/*     below the given tolerance, causing early exit. */

#line 1055 "sgsvj0.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1057 "sgsvj0.f"
L1995:

/*     Sort the vector D. */
#line 1060 "sgsvj0.f"
    i__1 = *n - 1;
#line 1060 "sgsvj0.f"
    for (p = 1; p <= i__1; ++p) {
#line 1061 "sgsvj0.f"
	i__2 = *n - p + 1;
#line 1061 "sgsvj0.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1062 "sgsvj0.f"
	if (p != q) {
#line 1063 "sgsvj0.f"
	    temp1 = sva[p];
#line 1064 "sgsvj0.f"
	    sva[p] = sva[q];
#line 1065 "sgsvj0.f"
	    sva[q] = temp1;
#line 1066 "sgsvj0.f"
	    temp1 = d__[p];
#line 1067 "sgsvj0.f"
	    d__[p] = d__[q];
#line 1068 "sgsvj0.f"
	    d__[q] = temp1;
#line 1069 "sgsvj0.f"
	    sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1070 "sgsvj0.f"
	    if (rsvec) {
#line 1070 "sgsvj0.f"
		sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1070 "sgsvj0.f"
	    }
#line 1071 "sgsvj0.f"
	}
#line 1072 "sgsvj0.f"
/* L5991: */
#line 1072 "sgsvj0.f"
    }

#line 1074 "sgsvj0.f"
    return 0;
/*     .. */
/*     .. END OF SGSVJ0 */
/*     .. */
} /* sgsvj0_ */

