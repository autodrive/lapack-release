#line 1 "cgsvj0.f"
/* cgsvj0.f -- translated by f2c (version 20100827).
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

#line 1 "cgsvj0.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b27 = 1.;

/* > \brief \b CGSVJ0 pre-processor for the routine sgesvj. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGSVJ0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgsvj0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgsvj0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgsvj0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, */
/*                          SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP */
/*       REAL               EPS, SFMIN, TOL */
/*       CHARACTER*1        JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK ) */
/*       REAL               SVA( N ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGSVJ0 is called from CGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as CGESVJ does, but */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          D is COMPLEX array, dimension (N) */
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
/* >          SVA is REAL array, dimension (N) */
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
/* >          V is COMPLEX array, dimension (LDV,N) */
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
/* >          WORK is COMPLEX array, dimension LWORK. */
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

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > CGSVJ0 is used just to enable CGESVJ to call a simplified version of */
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
/* Subroutine */ int cgsvj0_(char *jobv, integer *m, integer *n, 
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
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublecomplex ompq;
    static doublereal aapp0, aapq1, temp1;
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

#line 279 "cgsvj0.f"
    /* Parameter adjustments */
#line 279 "cgsvj0.f"
    --sva;
#line 279 "cgsvj0.f"
    --d__;
#line 279 "cgsvj0.f"
    a_dim1 = *lda;
#line 279 "cgsvj0.f"
    a_offset = 1 + a_dim1;
#line 279 "cgsvj0.f"
    a -= a_offset;
#line 279 "cgsvj0.f"
    v_dim1 = *ldv;
#line 279 "cgsvj0.f"
    v_offset = 1 + v_dim1;
#line 279 "cgsvj0.f"
    v -= v_offset;
#line 279 "cgsvj0.f"
    --work;
#line 279 "cgsvj0.f"

#line 279 "cgsvj0.f"
    /* Function Body */
#line 279 "cgsvj0.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 280 "cgsvj0.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 281 "cgsvj0.f"
    if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 282 "cgsvj0.f"
	*info = -1;
#line 283 "cgsvj0.f"
    } else if (*m < 0) {
#line 284 "cgsvj0.f"
	*info = -2;
#line 285 "cgsvj0.f"
    } else if (*n < 0 || *n > *m) {
#line 286 "cgsvj0.f"
	*info = -3;
#line 287 "cgsvj0.f"
    } else if (*lda < *m) {
#line 288 "cgsvj0.f"
	*info = -5;
#line 289 "cgsvj0.f"
    } else if ((rsvec || applv) && *mv < 0) {
#line 290 "cgsvj0.f"
	*info = -8;
#line 291 "cgsvj0.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 293 "cgsvj0.f"
	*info = -10;
#line 294 "cgsvj0.f"
    } else if (*tol <= *eps) {
#line 295 "cgsvj0.f"
	*info = -13;
#line 296 "cgsvj0.f"
    } else if (*nsweep < 0) {
#line 297 "cgsvj0.f"
	*info = -14;
#line 298 "cgsvj0.f"
    } else if (*lwork < *m) {
#line 299 "cgsvj0.f"
	*info = -16;
#line 300 "cgsvj0.f"
    } else {
#line 301 "cgsvj0.f"
	*info = 0;
#line 302 "cgsvj0.f"
    }

/*     #:( */
#line 305 "cgsvj0.f"
    if (*info != 0) {
#line 306 "cgsvj0.f"
	i__1 = -(*info);
#line 306 "cgsvj0.f"
	xerbla_("CGSVJ0", &i__1, (ftnlen)6);
#line 307 "cgsvj0.f"
	return 0;
#line 308 "cgsvj0.f"
    }

#line 310 "cgsvj0.f"
    if (rsvec) {
#line 311 "cgsvj0.f"
	mvl = *n;
#line 312 "cgsvj0.f"
    } else if (applv) {
#line 313 "cgsvj0.f"
	mvl = *mv;
#line 314 "cgsvj0.f"
    }
#line 315 "cgsvj0.f"
    rsvec = rsvec || applv;
#line 317 "cgsvj0.f"
    rooteps = sqrt(*eps);
#line 318 "cgsvj0.f"
    rootsfmin = sqrt(*sfmin);
#line 319 "cgsvj0.f"
    small = *sfmin / *eps;
#line 320 "cgsvj0.f"
    big = 1. / *sfmin;
#line 321 "cgsvj0.f"
    rootbig = 1. / rootsfmin;
#line 322 "cgsvj0.f"
    bigtheta = 1. / rooteps;
#line 323 "cgsvj0.f"
    roottol = sqrt(*tol);

/*     .. Row-cyclic Jacobi SVD algorithm with column pivoting .. */

#line 327 "cgsvj0.f"
    emptsw = *n * (*n - 1) / 2;
#line 328 "cgsvj0.f"
    notrot = 0;

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 333 "cgsvj0.f"
    swband = 0;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if CGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm CGEJSV. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 341 "cgsvj0.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 347 "cgsvj0.f"
    nbl = *n / kbl;
#line 348 "cgsvj0.f"
    if (nbl * kbl != *n) {
#line 348 "cgsvj0.f"
	++nbl;
#line 348 "cgsvj0.f"
    }

/* Computing 2nd power */
#line 350 "cgsvj0.f"
    i__1 = kbl;
#line 350 "cgsvj0.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 353 "cgsvj0.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 356 "cgsvj0.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */


/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 367 "cgsvj0.f"
    i__1 = *nsweep;
#line 367 "cgsvj0.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     .. go go go ... */

#line 371 "cgsvj0.f"
	mxaapq = 0.;
#line 372 "cgsvj0.f"
	mxsinj = 0.;
#line 373 "cgsvj0.f"
	iswrot = 0;

#line 375 "cgsvj0.f"
	notrot = 0;
#line 376 "cgsvj0.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 383 "cgsvj0.f"
	i__2 = nbl;
#line 383 "cgsvj0.f"
	for (ibr = 1; ibr <= i__2; ++ibr) {

#line 385 "cgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 387 "cgsvj0.f"
	    i__4 = lkahead, i__5 = nbl - ibr;
#line 387 "cgsvj0.f"
	    i__3 = min(i__4,i__5);
#line 387 "cgsvj0.f"
	    for (ir1 = 0; ir1 <= i__3; ++ir1) {

#line 389 "cgsvj0.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 391 "cgsvj0.f"
		i__5 = igl + kbl - 1, i__6 = *n - 1;
#line 391 "cgsvj0.f"
		i__4 = min(i__5,i__6);
#line 391 "cgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

/*     .. de Rijk's pivoting */

#line 395 "cgsvj0.f"
		    i__5 = *n - p + 1;
#line 395 "cgsvj0.f"
		    q = isamax_(&i__5, &sva[p], &c__1) + p - 1;
#line 396 "cgsvj0.f"
		    if (p != q) {
#line 397 "cgsvj0.f"
			cswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 398 "cgsvj0.f"
			if (rsvec) {
#line 398 "cgsvj0.f"
			    cswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 398 "cgsvj0.f"
			}
#line 400 "cgsvj0.f"
			temp1 = sva[p];
#line 401 "cgsvj0.f"
			sva[p] = sva[q];
#line 402 "cgsvj0.f"
			sva[q] = temp1;
#line 403 "cgsvj0.f"
			i__5 = p;
#line 403 "cgsvj0.f"
			aapq.r = d__[i__5].r, aapq.i = d__[i__5].i;
#line 404 "cgsvj0.f"
			i__5 = p;
#line 404 "cgsvj0.f"
			i__6 = q;
#line 404 "cgsvj0.f"
			d__[i__5].r = d__[i__6].r, d__[i__5].i = d__[i__6].i;
#line 405 "cgsvj0.f"
			i__5 = q;
#line 405 "cgsvj0.f"
			d__[i__5].r = aapq.r, d__[i__5].i = aapq.i;
#line 406 "cgsvj0.f"
		    }

#line 408 "cgsvj0.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/*        Caveat: */
/*        Unfortunately, some BLAS implementations compute SNCRM2(M,A(1,p),1) */
/*        as SQRT(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to */
/*        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to */
/*        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
/*        Hence, SCNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF */
/*        below should be replaced with "AAPP = SCNRM2( M, A(1,p), 1 )". */

#line 422 "cgsvj0.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 424 "cgsvj0.f"
			    sva[p] = scnrm2_(m, &a[p * a_dim1 + 1], &c__1);
#line 425 "cgsvj0.f"
			} else {
#line 426 "cgsvj0.f"
			    temp1 = 0.;
#line 427 "cgsvj0.f"
			    aapp = 1.;
#line 428 "cgsvj0.f"
			    classq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 429 "cgsvj0.f"
			    sva[p] = temp1 * sqrt(aapp);
#line 430 "cgsvj0.f"
			}
#line 431 "cgsvj0.f"
			aapp = sva[p];
#line 432 "cgsvj0.f"
		    } else {
#line 433 "cgsvj0.f"
			aapp = sva[p];
#line 434 "cgsvj0.f"
		    }

#line 436 "cgsvj0.f"
		    if (aapp > 0.) {

#line 438 "cgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 440 "cgsvj0.f"
			i__6 = igl + kbl - 1;
#line 440 "cgsvj0.f"
			i__5 = min(i__6,*n);
#line 440 "cgsvj0.f"
			for (q = p + 1; q <= i__5; ++q) {

#line 442 "cgsvj0.f"
			    aaqq = sva[q];

#line 444 "cgsvj0.f"
			    if (aaqq > 0.) {

#line 446 "cgsvj0.f"
				aapp0 = aapp;
#line 447 "cgsvj0.f"
				if (aaqq >= 1.) {
#line 448 "cgsvj0.f"
				    rotok = small * aapp <= aaqq;
#line 449 "cgsvj0.f"
				    if (aapp < big / aaqq) {
#line 450 "cgsvj0.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 450 "cgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 450 "cgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 450 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 452 "cgsvj0.f"
				    } else {
#line 453 "cgsvj0.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 455 "cgsvj0.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 457 "cgsvj0.f"
					cdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 457 "cgsvj0.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 457 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 459 "cgsvj0.f"
				    }
#line 460 "cgsvj0.f"
				} else {
#line 461 "cgsvj0.f"
				    rotok = aapp <= aaqq / small;
#line 462 "cgsvj0.f"
				    if (aapp > small / aaqq) {
#line 463 "cgsvj0.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 463 "cgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 463 "cgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 463 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 465 "cgsvj0.f"
				    } else {
#line 466 "cgsvj0.f"
					ccopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 468 "cgsvj0.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 471 "cgsvj0.f"
					cdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 471 "cgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 471 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 473 "cgsvj0.f"
				    }
#line 474 "cgsvj0.f"
				}

#line 476 "cgsvj0.f"
				d__1 = z_abs(&aapq);
#line 476 "cgsvj0.f"
				z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					d__1;
#line 476 "cgsvj0.f"
				ompq.r = z__1.r, ompq.i = z__1.i;
/*                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q) */
#line 478 "cgsvj0.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 479 "cgsvj0.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 479 "cgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 483 "cgsvj0.f"
				if (abs(aapq1) > *tol) {

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 488 "cgsvj0.f"
				    if (ir1 == 0) {
#line 489 "cgsvj0.f"
					notrot = 0;
#line 490 "cgsvj0.f"
					pskipped = 0;
#line 491 "cgsvj0.f"
					++iswrot;
#line 492 "cgsvj0.f"
				    }

#line 494 "cgsvj0.f"
				    if (rotok) {

#line 496 "cgsvj0.f"
					aqoap = aaqq / aapp;
#line 497 "cgsvj0.f"
					apoaq = aapp / aaqq;
#line 498 "cgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;

#line 500 "cgsvj0.f"
					if (abs(theta) > bigtheta) {

#line 502 "cgsvj0.f"
					    t = .5 / theta;
#line 503 "cgsvj0.f"
					    cs = 1.;
#line 505 "cgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 505 "cgsvj0.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 505 "cgsvj0.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 507 "cgsvj0.f"
					    if (rsvec) {
#line 508 "cgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 508 "cgsvj0.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 508 "cgsvj0.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 510 "cgsvj0.f"
					    }
/* Computing MAX */
#line 512 "cgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 512 "cgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 514 "cgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 514 "cgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 516 "cgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 516 "cgsvj0.f"
					    mxsinj = max(d__1,d__2);

#line 518 "cgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 522 "cgsvj0.f"
					    thsign = -d_sign(&c_b27, &aapq1);
#line 523 "cgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 525 "cgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 526 "cgsvj0.f"
					    sn = t * cs;

/* Computing MAX */
#line 528 "cgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 528 "cgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 529 "cgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 529 "cgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 531 "cgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 531 "cgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 534 "cgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 534 "cgsvj0.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 534 "cgsvj0.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 536 "cgsvj0.f"
					    if (rsvec) {
#line 537 "cgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 537 "cgsvj0.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 537 "cgsvj0.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 539 "cgsvj0.f"
					    }
#line 540 "cgsvj0.f"
					}
#line 541 "cgsvj0.f"
					i__6 = p;
#line 541 "cgsvj0.f"
					i__7 = q;
#line 541 "cgsvj0.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 541 "cgsvj0.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 541 "cgsvj0.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 543 "cgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 545 "cgsvj0.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 547 "cgsvj0.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 550 "cgsvj0.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 552 "cgsvj0.f"
					z__1.r = -aapq.r, z__1.i = -aapq.i;
#line 552 "cgsvj0.f"
					caxpy_(m, &z__1, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 554 "cgsvj0.f"
					clascl_("G", &c__0, &c__0, &c_b27, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 556 "cgsvj0.f"
					d__1 = 0., d__2 = 1. - aapq1 * aapq1;
#line 556 "cgsvj0.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 558 "cgsvj0.f"
					mxsinj = max(mxsinj,*sfmin);
#line 559 "cgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 565 "cgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 565 "cgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 567 "cgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 569 "cgsvj0.f"
					    sva[q] = scnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 570 "cgsvj0.f"
					} else {
#line 571 "cgsvj0.f"
					    t = 0.;
#line 572 "cgsvj0.f"
					    aaqq = 1.;
#line 573 "cgsvj0.f"
					    classq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 575 "cgsvj0.f"
					    sva[q] = t * sqrt(aaqq);
#line 576 "cgsvj0.f"
					}
#line 577 "cgsvj0.f"
				    }
#line 578 "cgsvj0.f"
				    if (aapp / aapp0 <= rooteps) {
#line 579 "cgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 581 "cgsvj0.f"
					    aapp = scnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 582 "cgsvj0.f"
					} else {
#line 583 "cgsvj0.f"
					    t = 0.;
#line 584 "cgsvj0.f"
					    aapp = 1.;
#line 585 "cgsvj0.f"
					    classq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 587 "cgsvj0.f"
					    aapp = t * sqrt(aapp);
#line 588 "cgsvj0.f"
					}
#line 589 "cgsvj0.f"
					sva[p] = aapp;
#line 590 "cgsvj0.f"
				    }

#line 592 "cgsvj0.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 594 "cgsvj0.f"
				    if (ir1 == 0) {
#line 594 "cgsvj0.f"
					++notrot;
#line 594 "cgsvj0.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 596 "cgsvj0.f"
				    ++pskipped;
#line 597 "cgsvj0.f"
				}
#line 598 "cgsvj0.f"
			    } else {
/*        A(:,q) is zero column */
#line 600 "cgsvj0.f"
				if (ir1 == 0) {
#line 600 "cgsvj0.f"
				    ++notrot;
#line 600 "cgsvj0.f"
				}
#line 601 "cgsvj0.f"
				++pskipped;
#line 602 "cgsvj0.f"
			    }

#line 604 "cgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 606 "cgsvj0.f"
				if (ir1 == 0) {
#line 606 "cgsvj0.f"
				    aapp = -aapp;
#line 606 "cgsvj0.f"
				}
#line 607 "cgsvj0.f"
				notrot = 0;
#line 608 "cgsvj0.f"
				goto L2103;
#line 609 "cgsvj0.f"
			    }

#line 611 "cgsvj0.f"
/* L2002: */
#line 611 "cgsvj0.f"
			}
/*     END q-LOOP */

#line 614 "cgsvj0.f"
L2103:
/*     bailed out of q-loop */

#line 617 "cgsvj0.f"
			sva[p] = aapp;

#line 619 "cgsvj0.f"
		    } else {
#line 620 "cgsvj0.f"
			sva[p] = aapp;
#line 621 "cgsvj0.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 621 "cgsvj0.f"
			    i__5 = igl + kbl - 1;
#line 621 "cgsvj0.f"
			    notrot = notrot + min(i__5,*n) - p;
#line 621 "cgsvj0.f"
			}
#line 623 "cgsvj0.f"
		    }

#line 625 "cgsvj0.f"
/* L2001: */
#line 625 "cgsvj0.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 628 "cgsvj0.f"
/* L1002: */
#line 628 "cgsvj0.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 633 "cgsvj0.f"
	    igl = (ibr - 1) * kbl + 1;

#line 635 "cgsvj0.f"
	    i__3 = nbl;
#line 635 "cgsvj0.f"
	    for (jbc = ibr + 1; jbc <= i__3; ++jbc) {

#line 637 "cgsvj0.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 641 "cgsvj0.f"
		ijblsk = 0;
/* Computing MIN */
#line 642 "cgsvj0.f"
		i__5 = igl + kbl - 1;
#line 642 "cgsvj0.f"
		i__4 = min(i__5,*n);
#line 642 "cgsvj0.f"
		for (p = igl; p <= i__4; ++p) {

#line 644 "cgsvj0.f"
		    aapp = sva[p];
#line 645 "cgsvj0.f"
		    if (aapp > 0.) {

#line 647 "cgsvj0.f"
			pskipped = 0;

/* Computing MIN */
#line 649 "cgsvj0.f"
			i__6 = jgl + kbl - 1;
#line 649 "cgsvj0.f"
			i__5 = min(i__6,*n);
#line 649 "cgsvj0.f"
			for (q = jgl; q <= i__5; ++q) {

#line 651 "cgsvj0.f"
			    aaqq = sva[q];
#line 652 "cgsvj0.f"
			    if (aaqq > 0.) {
#line 653 "cgsvj0.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 659 "cgsvj0.f"
				if (aaqq >= 1.) {
#line 660 "cgsvj0.f"
				    if (aapp >= aaqq) {
#line 661 "cgsvj0.f"
					rotok = small * aapp <= aaqq;
#line 662 "cgsvj0.f"
				    } else {
#line 663 "cgsvj0.f"
					rotok = small * aaqq <= aapp;
#line 664 "cgsvj0.f"
				    }
#line 665 "cgsvj0.f"
				    if (aapp < big / aaqq) {
#line 666 "cgsvj0.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 666 "cgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 666 "cgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 666 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 668 "cgsvj0.f"
				    } else {
#line 669 "cgsvj0.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 671 "cgsvj0.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 674 "cgsvj0.f"
					cdotc_(&z__2, m, &work[1], &c__1, &a[
						q * a_dim1 + 1], &c__1);
#line 674 "cgsvj0.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 674 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 676 "cgsvj0.f"
				    }
#line 677 "cgsvj0.f"
				} else {
#line 678 "cgsvj0.f"
				    if (aapp >= aaqq) {
#line 679 "cgsvj0.f"
					rotok = aapp <= aaqq / small;
#line 680 "cgsvj0.f"
				    } else {
#line 681 "cgsvj0.f"
					rotok = aaqq <= aapp / small;
#line 682 "cgsvj0.f"
				    }
#line 683 "cgsvj0.f"
				    if (aapp > small / aaqq) {
#line 684 "cgsvj0.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 684 "cgsvj0.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 684 "cgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 684 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 686 "cgsvj0.f"
				    } else {
#line 687 "cgsvj0.f"
					ccopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[1], &c__1);
#line 689 "cgsvj0.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b27, m, &c__1, &work[1], 
						lda, &ierr, (ftnlen)1);
#line 692 "cgsvj0.f"
					cdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &work[1], &c__1);
#line 692 "cgsvj0.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 692 "cgsvj0.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 694 "cgsvj0.f"
				    }
#line 695 "cgsvj0.f"
				}

#line 697 "cgsvj0.f"
				d__1 = z_abs(&aapq);
#line 697 "cgsvj0.f"
				z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					d__1;
#line 697 "cgsvj0.f"
				ompq.r = z__1.r, ompq.i = z__1.i;
/*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
#line 699 "cgsvj0.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 700 "cgsvj0.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 700 "cgsvj0.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 704 "cgsvj0.f"
				if (abs(aapq1) > *tol) {
#line 705 "cgsvj0.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 707 "cgsvj0.f"
				    pskipped = 0;
#line 708 "cgsvj0.f"
				    ++iswrot;

#line 710 "cgsvj0.f"
				    if (rotok) {

#line 712 "cgsvj0.f"
					aqoap = aaqq / aapp;
#line 713 "cgsvj0.f"
					apoaq = aapp / aaqq;
#line 714 "cgsvj0.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 715 "cgsvj0.f"
					if (aaqq > aapp0) {
#line 715 "cgsvj0.f"
					    theta = -theta;
#line 715 "cgsvj0.f"
					}

#line 717 "cgsvj0.f"
					if (abs(theta) > bigtheta) {
#line 718 "cgsvj0.f"
					    t = .5 / theta;
#line 719 "cgsvj0.f"
					    cs = 1.;
#line 720 "cgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 720 "cgsvj0.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 720 "cgsvj0.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 722 "cgsvj0.f"
					    if (rsvec) {
#line 723 "cgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 723 "cgsvj0.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 723 "cgsvj0.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 725 "cgsvj0.f"
					    }
/* Computing MAX */
#line 726 "cgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 726 "cgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 728 "cgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 728 "cgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 730 "cgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 730 "cgsvj0.f"
					    mxsinj = max(d__1,d__2);
#line 731 "cgsvj0.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 735 "cgsvj0.f"
					    thsign = -d_sign(&c_b27, &aapq1);
#line 736 "cgsvj0.f"
					    if (aaqq > aapp0) {
#line 736 "cgsvj0.f"
			  thsign = -thsign;
#line 736 "cgsvj0.f"
					    }
#line 737 "cgsvj0.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 739 "cgsvj0.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 740 "cgsvj0.f"
					    sn = t * cs;
/* Computing MAX */
#line 741 "cgsvj0.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 741 "cgsvj0.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 742 "cgsvj0.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 742 "cgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 744 "cgsvj0.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 744 "cgsvj0.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 747 "cgsvj0.f"
					    d_cnjg(&z__2, &ompq);
#line 747 "cgsvj0.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 747 "cgsvj0.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 749 "cgsvj0.f"
					    if (rsvec) {
#line 750 "cgsvj0.f"
			  d_cnjg(&z__2, &ompq);
#line 750 "cgsvj0.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 750 "cgsvj0.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 752 "cgsvj0.f"
					    }
#line 753 "cgsvj0.f"
					}
#line 754 "cgsvj0.f"
					i__6 = p;
#line 754 "cgsvj0.f"
					i__7 = q;
#line 754 "cgsvj0.f"
					z__2.r = -d__[i__7].r, z__2.i = -d__[
						i__7].i;
#line 754 "cgsvj0.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 754 "cgsvj0.f"
					d__[i__6].r = z__1.r, d__[i__6].i = 
						z__1.i;

#line 756 "cgsvj0.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 758 "cgsvj0.f"
					if (aapp > aaqq) {
#line 759 "cgsvj0.f"
					    ccopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 761 "cgsvj0.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b27, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 764 "cgsvj0.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b27, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 767 "cgsvj0.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 767 "cgsvj0.f"
					    caxpy_(m, &z__1, &work[1], &c__1, 
						    &a[q * a_dim1 + 1], &c__1)
						    ;
#line 769 "cgsvj0.f"
					    clascl_("G", &c__0, &c__0, &c_b27,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 772 "cgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 772 "cgsvj0.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 774 "cgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 775 "cgsvj0.f"
					} else {
#line 776 "cgsvj0.f"
					    ccopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[1], &c__1);
#line 778 "cgsvj0.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b27, m, &c__1, &work[1]
						    , lda, &ierr, (ftnlen)1);
#line 781 "cgsvj0.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b27, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 784 "cgsvj0.f"
					    d_cnjg(&z__2, &aapq);
#line 784 "cgsvj0.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 784 "cgsvj0.f"
					    caxpy_(m, &z__1, &work[1], &c__1, 
						    &a[p * a_dim1 + 1], &c__1)
						    ;
#line 786 "cgsvj0.f"
					    clascl_("G", &c__0, &c__0, &c_b27,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 789 "cgsvj0.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 789 "cgsvj0.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 791 "cgsvj0.f"
					    mxsinj = max(mxsinj,*sfmin);
#line 792 "cgsvj0.f"
					}
#line 793 "cgsvj0.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 798 "cgsvj0.f"
				    d__1 = sva[q] / aaqq;
#line 798 "cgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 800 "cgsvj0.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 802 "cgsvj0.f"
					    sva[q] = scnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 803 "cgsvj0.f"
					} else {
#line 804 "cgsvj0.f"
					    t = 0.;
#line 805 "cgsvj0.f"
					    aaqq = 1.;
#line 806 "cgsvj0.f"
					    classq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 808 "cgsvj0.f"
					    sva[q] = t * sqrt(aaqq);
#line 809 "cgsvj0.f"
					}
#line 810 "cgsvj0.f"
				    }
/* Computing 2nd power */
#line 811 "cgsvj0.f"
				    d__1 = aapp / aapp0;
#line 811 "cgsvj0.f"
				    if (d__1 * d__1 <= rooteps) {
#line 812 "cgsvj0.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 814 "cgsvj0.f"
					    aapp = scnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 815 "cgsvj0.f"
					} else {
#line 816 "cgsvj0.f"
					    t = 0.;
#line 817 "cgsvj0.f"
					    aapp = 1.;
#line 818 "cgsvj0.f"
					    classq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 820 "cgsvj0.f"
					    aapp = t * sqrt(aapp);
#line 821 "cgsvj0.f"
					}
#line 822 "cgsvj0.f"
					sva[p] = aapp;
#line 823 "cgsvj0.f"
				    }
/*              end of OK rotation */
#line 825 "cgsvj0.f"
				} else {
#line 826 "cgsvj0.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 828 "cgsvj0.f"
				    ++pskipped;
#line 829 "cgsvj0.f"
				    ++ijblsk;
#line 830 "cgsvj0.f"
				}
#line 831 "cgsvj0.f"
			    } else {
#line 832 "cgsvj0.f"
				++notrot;
#line 833 "cgsvj0.f"
				++pskipped;
#line 834 "cgsvj0.f"
				++ijblsk;
#line 835 "cgsvj0.f"
			    }

#line 837 "cgsvj0.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 839 "cgsvj0.f"
				sva[p] = aapp;
#line 840 "cgsvj0.f"
				notrot = 0;
#line 841 "cgsvj0.f"
				goto L2011;
#line 842 "cgsvj0.f"
			    }
#line 843 "cgsvj0.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 845 "cgsvj0.f"
				aapp = -aapp;
#line 846 "cgsvj0.f"
				notrot = 0;
#line 847 "cgsvj0.f"
				goto L2203;
#line 848 "cgsvj0.f"
			    }

#line 850 "cgsvj0.f"
/* L2200: */
#line 850 "cgsvj0.f"
			}
/*        end of the q-loop */
#line 852 "cgsvj0.f"
L2203:

#line 854 "cgsvj0.f"
			sva[p] = aapp;

#line 856 "cgsvj0.f"
		    } else {

#line 858 "cgsvj0.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 858 "cgsvj0.f"
			    i__5 = jgl + kbl - 1;
#line 858 "cgsvj0.f"
			    notrot = notrot + min(i__5,*n) - jgl + 1;
#line 858 "cgsvj0.f"
			}
#line 860 "cgsvj0.f"
			if (aapp < 0.) {
#line 860 "cgsvj0.f"
			    notrot = 0;
#line 860 "cgsvj0.f"
			}

#line 862 "cgsvj0.f"
		    }

#line 864 "cgsvj0.f"
/* L2100: */
#line 864 "cgsvj0.f"
		}
/*     end of the p-loop */
#line 866 "cgsvj0.f"
/* L2010: */
#line 866 "cgsvj0.f"
	    }
/*     end of the jbc-loop */
#line 868 "cgsvj0.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 870 "cgsvj0.f"
	    i__4 = igl + kbl - 1;
#line 870 "cgsvj0.f"
	    i__3 = min(i__4,*n);
#line 870 "cgsvj0.f"
	    for (p = igl; p <= i__3; ++p) {
#line 871 "cgsvj0.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 872 "cgsvj0.f"
/* L2012: */
#line 872 "cgsvj0.f"
	    }
/* ** */
#line 874 "cgsvj0.f"
/* L2000: */
#line 874 "cgsvj0.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 878 "cgsvj0.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 880 "cgsvj0.f"
	    sva[*n] = scnrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 881 "cgsvj0.f"
	} else {
#line 882 "cgsvj0.f"
	    t = 0.;
#line 883 "cgsvj0.f"
	    aapp = 1.;
#line 884 "cgsvj0.f"
	    classq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 885 "cgsvj0.f"
	    sva[*n] = t * sqrt(aapp);
#line 886 "cgsvj0.f"
	}

/*     Additional steering devices */

#line 890 "cgsvj0.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 890 "cgsvj0.f"
	    swband = i__;
#line 890 "cgsvj0.f"
	}

#line 893 "cgsvj0.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * *tol && (
		doublereal) (*n) * mxaapq * mxsinj < *tol) {
#line 895 "cgsvj0.f"
	    goto L1994;
#line 896 "cgsvj0.f"
	}

#line 898 "cgsvj0.f"
	if (notrot >= emptsw) {
#line 898 "cgsvj0.f"
	    goto L1994;
#line 898 "cgsvj0.f"
	}

#line 900 "cgsvj0.f"
/* L1993: */
#line 900 "cgsvj0.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 904 "cgsvj0.f"
    *info = *nsweep - 1;
#line 905 "cgsvj0.f"
    goto L1995;

#line 907 "cgsvj0.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 911 "cgsvj0.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 913 "cgsvj0.f"
L1995:

/*     Sort the vector SVA() of column norms. */
#line 916 "cgsvj0.f"
    i__1 = *n - 1;
#line 916 "cgsvj0.f"
    for (p = 1; p <= i__1; ++p) {
#line 917 "cgsvj0.f"
	i__2 = *n - p + 1;
#line 917 "cgsvj0.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 918 "cgsvj0.f"
	if (p != q) {
#line 919 "cgsvj0.f"
	    temp1 = sva[p];
#line 920 "cgsvj0.f"
	    sva[p] = sva[q];
#line 921 "cgsvj0.f"
	    sva[q] = temp1;
#line 922 "cgsvj0.f"
	    i__2 = p;
#line 922 "cgsvj0.f"
	    aapq.r = d__[i__2].r, aapq.i = d__[i__2].i;
#line 923 "cgsvj0.f"
	    i__2 = p;
#line 923 "cgsvj0.f"
	    i__3 = q;
#line 923 "cgsvj0.f"
	    d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
#line 924 "cgsvj0.f"
	    i__2 = q;
#line 924 "cgsvj0.f"
	    d__[i__2].r = aapq.r, d__[i__2].i = aapq.i;
#line 925 "cgsvj0.f"
	    cswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 926 "cgsvj0.f"
	    if (rsvec) {
#line 926 "cgsvj0.f"
		cswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 926 "cgsvj0.f"
	    }
#line 927 "cgsvj0.f"
	}
#line 928 "cgsvj0.f"
/* L5991: */
#line 928 "cgsvj0.f"
    }

#line 930 "cgsvj0.f"
    return 0;
/*     .. */
/*     .. END OF CGSVJ0 */
/*     .. */
} /* cgsvj0_ */

