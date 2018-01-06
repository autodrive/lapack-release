#line 1 "dgesdd.f"
/* dgesdd.f -- translated by f2c (version 20100827).
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

#line 1 "dgesdd.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b63 = 0.;
static integer c__1 = 1;
static doublereal c_b84 = 1.;

/* > \brief \b DGESDD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGESDD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesdd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesdd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesdd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGESDD computes the singular value decomposition (SVD) of a real */
/* > M-by-N matrix A, optionally computing the left and right singular */
/* > vectors.  If singular vectors are desired, it uses a */
/* > divide-and-conquer algorithm. */
/* > */
/* > The SVD is written */
/* > */
/* >      A = U * SIGMA * transpose(V) */
/* > */
/* > where SIGMA is an M-by-N matrix which is zero except for its */
/* > min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and */
/* > V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA */
/* > are the singular values of A; they are real and non-negative, and */
/* > are returned in descending order.  The first min(m,n) columns of */
/* > U and V are the left and right singular vectors of A. */
/* > */
/* > Note that the routine returns VT = V**T, not V. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          Specifies options for computing all or part of the matrix U: */
/* >          = 'A':  all M columns of U and all N rows of V**T are */
/* >                  returned in the arrays U and VT; */
/* >          = 'S':  the first min(M,N) columns of U and the first */
/* >                  min(M,N) rows of V**T are returned in the arrays U */
/* >                  and VT; */
/* >          = 'O':  If M >= N, the first N columns of U are overwritten */
/* >                  on the array A and all rows of V**T are returned in */
/* >                  the array VT; */
/* >                  otherwise, all columns of U are returned in the */
/* >                  array U and the first M rows of V**T are overwritten */
/* >                  in the array A; */
/* >          = 'N':  no columns of U or rows of V**T are computed. */
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
/* >          The number of columns of the input matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          if JOBZ = 'O',  A is overwritten with the first N columns */
/* >                          of U (the left singular vectors, stored */
/* >                          columnwise) if M >= N; */
/* >                          A is overwritten with the first M rows */
/* >                          of V**T (the right singular vectors, stored */
/* >                          rowwise) otherwise. */
/* >          if JOBZ .ne. 'O', the contents of A are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension (LDU,UCOL) */
/* >          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N; */
/* >          UCOL = min(M,N) if JOBZ = 'S'. */
/* >          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M */
/* >          orthogonal matrix U; */
/* >          if JOBZ = 'S', U contains the first min(M,N) columns of U */
/* >          (the left singular vectors, stored columnwise); */
/* >          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U.  LDU >= 1; if */
/* >          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension (LDVT,N) */
/* >          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the */
/* >          N-by-N orthogonal matrix V**T; */
/* >          if JOBZ = 'S', VT contains the first min(M,N) rows of */
/* >          V**T (the right singular vectors, stored rowwise); */
/* >          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >          The leading dimension of the array VT.  LDVT >= 1; */
/* >          if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N; */
/* >          if JOBZ = 'S', LDVT >= min(M,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK; */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= 1. */
/* >          If LWORK = -1, a workspace query is assumed.  The optimal */
/* >          size for the WORK array is calculated and stored in WORK(1), */
/* >          and no other work except argument checking is performed. */
/* > */
/* >          Let mx = max(M,N) and mn = min(M,N). */
/* >          If JOBZ = 'N', LWORK >= 3*mn + max( mx, 7*mn ). */
/* >          If JOBZ = 'O', LWORK >= 3*mn + max( mx, 5*mn*mn + 4*mn ). */
/* >          If JOBZ = 'S', LWORK >= 4*mn*mn + 7*mn. */
/* >          If JOBZ = 'A', LWORK >= 4*mn*mn + 6*mn + mx. */
/* >          These are not tight minimums in all cases; see comments inside code. */
/* >          For good performance, LWORK should generally be larger; */
/* >          a query is recommended. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (8*min(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  DBDSDC did not converge, updating process failed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup doubleGEsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *s, doublereal *u, integer *ldu, 
	doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *iwork, integer *info, ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer lwork_dorglq_mn__, lwork_dorglq_nn__, lwork_dorgqr_mm__, 
	    lwork_dorgqr_mn__, i__, ie, lwork_dorgbr_p_mm__, il, 
	    lwork_dorgbr_q_nn__, ir, iu, blk;
    static doublereal dum[1], eps;
    static integer ivt, iscl;
    static doublereal anrm;
    static integer idum[1], ierr, itau, lwork_dormbr_qln_mm__, 
	    lwork_dormbr_qln_mn__, lwork_dormbr_qln_nn__, 
	    lwork_dormbr_prt_mm__, lwork_dormbr_prt_mn__, 
	    lwork_dormbr_prt_nn__;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn, wrkbl, itaup, itauq, mnthr;
    static logical wntqa;
    static integer nwork;
    static logical wntqn, wntqo, wntqs;
    extern /* Subroutine */ int dbdsdc_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static integer bdspac;
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dlacpy_(char *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dorgbr_(char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), dorglq_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dorgqr_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *);
    static integer ldwrkl, ldwrkr, minwrk, ldwrku, maxwrk, ldwkvt;
    static doublereal smlnum;
    static logical wntqas, lquery;
    static integer lwork_dgebrd_mm__, lwork_dgebrd_mn__, lwork_dgebrd_nn__, 
	    lwork_dgelqf_mn__, lwork_dgeqrf_mn__;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 281 "dgesdd.f"
    /* Parameter adjustments */
#line 281 "dgesdd.f"
    a_dim1 = *lda;
#line 281 "dgesdd.f"
    a_offset = 1 + a_dim1;
#line 281 "dgesdd.f"
    a -= a_offset;
#line 281 "dgesdd.f"
    --s;
#line 281 "dgesdd.f"
    u_dim1 = *ldu;
#line 281 "dgesdd.f"
    u_offset = 1 + u_dim1;
#line 281 "dgesdd.f"
    u -= u_offset;
#line 281 "dgesdd.f"
    vt_dim1 = *ldvt;
#line 281 "dgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 281 "dgesdd.f"
    vt -= vt_offset;
#line 281 "dgesdd.f"
    --work;
#line 281 "dgesdd.f"
    --iwork;
#line 281 "dgesdd.f"

#line 281 "dgesdd.f"
    /* Function Body */
#line 281 "dgesdd.f"
    *info = 0;
#line 282 "dgesdd.f"
    minmn = min(*m,*n);
#line 283 "dgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 284 "dgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 285 "dgesdd.f"
    wntqas = wntqa || wntqs;
#line 286 "dgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 287 "dgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 288 "dgesdd.f"
    lquery = *lwork == -1;

#line 290 "dgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 291 "dgesdd.f"
	*info = -1;
#line 292 "dgesdd.f"
    } else if (*m < 0) {
#line 293 "dgesdd.f"
	*info = -2;
#line 294 "dgesdd.f"
    } else if (*n < 0) {
#line 295 "dgesdd.f"
	*info = -3;
#line 296 "dgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 297 "dgesdd.f"
	*info = -5;
#line 298 "dgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 300 "dgesdd.f"
	*info = -8;
#line 301 "dgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 304 "dgesdd.f"
	*info = -10;
#line 305 "dgesdd.f"
    }

/*     Compute workspace */
/*       Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace allocated at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */

#line 314 "dgesdd.f"
    if (*info == 0) {
#line 315 "dgesdd.f"
	minwrk = 1;
#line 316 "dgesdd.f"
	maxwrk = 1;
#line 317 "dgesdd.f"
	bdspac = 0;
#line 318 "dgesdd.f"
	mnthr = (integer) (minmn * 11. / 6.);
#line 319 "dgesdd.f"
	if (*m >= *n && minmn > 0) {

/*           Compute space needed for DBDSDC */

#line 323 "dgesdd.f"
	    if (wntqn) {
/*              dbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6) */
/*              keep 7*N for backwards compatability. */
#line 326 "dgesdd.f"
		bdspac = *n * 7;
#line 327 "dgesdd.f"
	    } else {
#line 328 "dgesdd.f"
		bdspac = *n * 3 * *n + (*n << 2);
#line 329 "dgesdd.f"
	    }

/*           Compute space preferred for each routine */
#line 332 "dgesdd.f"
	    dgebrd_(m, n, dum, m, dum, dum, dum, dum, dum, &c_n1, &ierr);
#line 334 "dgesdd.f"
	    lwork_dgebrd_mn__ = (integer) dum[0];

#line 336 "dgesdd.f"
	    dgebrd_(n, n, dum, n, dum, dum, dum, dum, dum, &c_n1, &ierr);
#line 338 "dgesdd.f"
	    lwork_dgebrd_nn__ = (integer) dum[0];

#line 340 "dgesdd.f"
	    dgeqrf_(m, n, dum, m, dum, dum, &c_n1, &ierr);
#line 341 "dgesdd.f"
	    lwork_dgeqrf_mn__ = (integer) dum[0];

#line 343 "dgesdd.f"
	    dorgbr_("Q", n, n, n, dum, n, dum, dum, &c_n1, &ierr, (ftnlen)1);
#line 345 "dgesdd.f"
	    lwork_dorgbr_q_nn__ = (integer) dum[0];

#line 347 "dgesdd.f"
	    dorgqr_(m, m, n, dum, m, dum, dum, &c_n1, &ierr);
#line 348 "dgesdd.f"
	    lwork_dorgqr_mm__ = (integer) dum[0];

#line 350 "dgesdd.f"
	    dorgqr_(m, n, n, dum, m, dum, dum, &c_n1, &ierr);
#line 351 "dgesdd.f"
	    lwork_dorgqr_mn__ = (integer) dum[0];

#line 353 "dgesdd.f"
	    dormbr_("P", "R", "T", n, n, n, dum, n, dum, dum, n, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 355 "dgesdd.f"
	    lwork_dormbr_prt_nn__ = (integer) dum[0];

#line 357 "dgesdd.f"
	    dormbr_("Q", "L", "N", n, n, n, dum, n, dum, dum, n, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 359 "dgesdd.f"
	    lwork_dormbr_qln_nn__ = (integer) dum[0];

#line 361 "dgesdd.f"
	    dormbr_("Q", "L", "N", m, n, n, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 363 "dgesdd.f"
	    lwork_dormbr_qln_mn__ = (integer) dum[0];

#line 365 "dgesdd.f"
	    dormbr_("Q", "L", "N", m, m, n, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 367 "dgesdd.f"
	    lwork_dormbr_qln_mm__ = (integer) dum[0];

#line 369 "dgesdd.f"
	    if (*m >= mnthr) {
#line 370 "dgesdd.f"
		if (wntqn) {

/*                 Path 1 (M >> N, JOBZ='N') */

#line 374 "dgesdd.f"
		    wrkbl = *n + lwork_dgeqrf_mn__;
/* Computing MAX */
#line 375 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
#line 375 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 376 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n;
#line 376 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 377 "dgesdd.f"
		    minwrk = bdspac + *n;
#line 378 "dgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M >> N, JOBZ='O') */

#line 382 "dgesdd.f"
		    wrkbl = *n + lwork_dgeqrf_mn__;
/* Computing MAX */
#line 383 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_dorgqr_mn__;
#line 383 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 384 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
#line 384 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 385 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_nn__;
#line 385 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 386 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
#line 386 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 387 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 387 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 388 "dgesdd.f"
		    maxwrk = wrkbl + (*n << 1) * *n;
#line 389 "dgesdd.f"
		    minwrk = bdspac + (*n << 1) * *n + *n * 3;
#line 390 "dgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M >> N, JOBZ='S') */

#line 394 "dgesdd.f"
		    wrkbl = *n + lwork_dgeqrf_mn__;
/* Computing MAX */
#line 395 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_dorgqr_mn__;
#line 395 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 396 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
#line 396 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 397 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_nn__;
#line 397 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 398 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
#line 398 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 399 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 399 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 400 "dgesdd.f"
		    maxwrk = wrkbl + *n * *n;
#line 401 "dgesdd.f"
		    minwrk = bdspac + *n * *n + *n * 3;
#line 402 "dgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M >> N, JOBZ='A') */

#line 406 "dgesdd.f"
		    wrkbl = *n + lwork_dgeqrf_mn__;
/* Computing MAX */
#line 407 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_dorgqr_mm__;
#line 407 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 408 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
#line 408 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 409 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_nn__;
#line 409 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 410 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
#line 410 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 411 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 411 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 412 "dgesdd.f"
		    maxwrk = wrkbl + *n * *n;
/* Computing MAX */
#line 413 "dgesdd.f"
		    i__1 = *n * 3 + bdspac, i__2 = *n + *m;
#line 413 "dgesdd.f"
		    minwrk = *n * *n + max(i__1,i__2);
#line 414 "dgesdd.f"
		}
#line 415 "dgesdd.f"
	    } else {

/*              Path 5 (M >= N, but not much larger) */

#line 419 "dgesdd.f"
		wrkbl = *n * 3 + lwork_dgebrd_mn__;
#line 420 "dgesdd.f"
		if (wntqn) {
/*                 Path 5n (M >= N, jobz='N') */
/* Computing MAX */
#line 422 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 422 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 423 "dgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 424 "dgesdd.f"
		} else if (wntqo) {
/*                 Path 5o (M >= N, jobz='O') */
/* Computing MAX */
#line 426 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
#line 426 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 427 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_mn__;
#line 427 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 428 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 428 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 429 "dgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 430 "dgesdd.f"
		    i__1 = *m, i__2 = *n * *n + bdspac;
#line 430 "dgesdd.f"
		    minwrk = *n * 3 + max(i__1,i__2);
#line 431 "dgesdd.f"
		} else if (wntqs) {
/*                 Path 5s (M >= N, jobz='S') */
/* Computing MAX */
#line 433 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_mn__;
#line 433 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 434 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
#line 434 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 435 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 435 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 436 "dgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 437 "dgesdd.f"
		} else if (wntqa) {
/*                 Path 5a (M >= N, jobz='A') */
/* Computing MAX */
#line 439 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_mm__;
#line 439 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 440 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
#line 440 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 441 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 441 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 442 "dgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 443 "dgesdd.f"
		}
#line 444 "dgesdd.f"
	    }
#line 445 "dgesdd.f"
	} else if (minmn > 0) {

/*           Compute space needed for DBDSDC */

#line 449 "dgesdd.f"
	    if (wntqn) {
/*              dbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6) */
/*              keep 7*N for backwards compatability. */
#line 452 "dgesdd.f"
		bdspac = *m * 7;
#line 453 "dgesdd.f"
	    } else {
#line 454 "dgesdd.f"
		bdspac = *m * 3 * *m + (*m << 2);
#line 455 "dgesdd.f"
	    }

/*           Compute space preferred for each routine */
#line 458 "dgesdd.f"
	    dgebrd_(m, n, dum, m, dum, dum, dum, dum, dum, &c_n1, &ierr);
#line 460 "dgesdd.f"
	    lwork_dgebrd_mn__ = (integer) dum[0];

#line 462 "dgesdd.f"
	    dgebrd_(m, m, &a[a_offset], m, &s[1], dum, dum, dum, dum, &c_n1, &
		    ierr);
#line 464 "dgesdd.f"
	    lwork_dgebrd_mm__ = (integer) dum[0];

#line 466 "dgesdd.f"
	    dgelqf_(m, n, &a[a_offset], m, dum, dum, &c_n1, &ierr);
#line 467 "dgesdd.f"
	    lwork_dgelqf_mn__ = (integer) dum[0];

#line 469 "dgesdd.f"
	    dorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
#line 470 "dgesdd.f"
	    lwork_dorglq_nn__ = (integer) dum[0];

#line 472 "dgesdd.f"
	    dorglq_(m, n, m, &a[a_offset], m, dum, dum, &c_n1, &ierr);
#line 473 "dgesdd.f"
	    lwork_dorglq_mn__ = (integer) dum[0];

#line 475 "dgesdd.f"
	    dorgbr_("P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 476 "dgesdd.f"
	    lwork_dorgbr_p_mm__ = (integer) dum[0];

#line 478 "dgesdd.f"
	    dormbr_("P", "R", "T", m, m, m, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 480 "dgesdd.f"
	    lwork_dormbr_prt_mm__ = (integer) dum[0];

#line 482 "dgesdd.f"
	    dormbr_("P", "R", "T", m, n, m, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 484 "dgesdd.f"
	    lwork_dormbr_prt_mn__ = (integer) dum[0];

#line 486 "dgesdd.f"
	    dormbr_("P", "R", "T", n, n, m, dum, n, dum, dum, n, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 488 "dgesdd.f"
	    lwork_dormbr_prt_nn__ = (integer) dum[0];

#line 490 "dgesdd.f"
	    dormbr_("Q", "L", "N", m, m, m, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 492 "dgesdd.f"
	    lwork_dormbr_qln_mm__ = (integer) dum[0];

#line 494 "dgesdd.f"
	    if (*n >= mnthr) {
#line 495 "dgesdd.f"
		if (wntqn) {

/*                 Path 1t (N >> M, JOBZ='N') */

#line 499 "dgesdd.f"
		    wrkbl = *m + lwork_dgelqf_mn__;
/* Computing MAX */
#line 500 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
#line 500 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 501 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m;
#line 501 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 502 "dgesdd.f"
		    minwrk = bdspac + *m;
#line 503 "dgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N >> M, JOBZ='O') */

#line 507 "dgesdd.f"
		    wrkbl = *m + lwork_dgelqf_mn__;
/* Computing MAX */
#line 508 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_dorglq_mn__;
#line 508 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 509 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
#line 509 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 510 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
#line 510 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 511 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mm__;
#line 511 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 512 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 512 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 513 "dgesdd.f"
		    maxwrk = wrkbl + (*m << 1) * *m;
#line 514 "dgesdd.f"
		    minwrk = bdspac + (*m << 1) * *m + *m * 3;
#line 515 "dgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N >> M, JOBZ='S') */

#line 519 "dgesdd.f"
		    wrkbl = *m + lwork_dgelqf_mn__;
/* Computing MAX */
#line 520 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_dorglq_mn__;
#line 520 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 521 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
#line 521 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 522 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
#line 522 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 523 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mm__;
#line 523 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 524 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 524 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 525 "dgesdd.f"
		    maxwrk = wrkbl + *m * *m;
#line 526 "dgesdd.f"
		    minwrk = bdspac + *m * *m + *m * 3;
#line 527 "dgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N >> M, JOBZ='A') */

#line 531 "dgesdd.f"
		    wrkbl = *m + lwork_dgelqf_mn__;
/* Computing MAX */
#line 532 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_dorglq_nn__;
#line 532 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 533 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
#line 533 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 534 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
#line 534 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 535 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mm__;
#line 535 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 536 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 536 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 537 "dgesdd.f"
		    maxwrk = wrkbl + *m * *m;
/* Computing MAX */
#line 538 "dgesdd.f"
		    i__1 = *m * 3 + bdspac, i__2 = *m + *n;
#line 538 "dgesdd.f"
		    minwrk = *m * *m + max(i__1,i__2);
#line 539 "dgesdd.f"
		}
#line 540 "dgesdd.f"
	    } else {

/*              Path 5t (N > M, but not much larger) */

#line 544 "dgesdd.f"
		wrkbl = *m * 3 + lwork_dgebrd_mn__;
#line 545 "dgesdd.f"
		if (wntqn) {
/*                 Path 5tn (N > M, jobz='N') */
/* Computing MAX */
#line 547 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 547 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 548 "dgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 549 "dgesdd.f"
		} else if (wntqo) {
/*                 Path 5to (N > M, jobz='O') */
/* Computing MAX */
#line 551 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
#line 551 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 552 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mn__;
#line 552 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 553 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 553 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 554 "dgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 555 "dgesdd.f"
		    i__1 = *n, i__2 = *m * *m + bdspac;
#line 555 "dgesdd.f"
		    minwrk = *m * 3 + max(i__1,i__2);
#line 556 "dgesdd.f"
		} else if (wntqs) {
/*                 Path 5ts (N > M, jobz='S') */
/* Computing MAX */
#line 558 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
#line 558 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 559 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mn__;
#line 559 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 560 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 560 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 561 "dgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 562 "dgesdd.f"
		} else if (wntqa) {
/*                 Path 5ta (N > M, jobz='A') */
/* Computing MAX */
#line 564 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
#line 564 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 565 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_nn__;
#line 565 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 566 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 566 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 567 "dgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 568 "dgesdd.f"
		}
#line 569 "dgesdd.f"
	    }
#line 570 "dgesdd.f"
	}
#line 572 "dgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 573 "dgesdd.f"
	work[1] = (doublereal) maxwrk;

#line 575 "dgesdd.f"
	if (*lwork < minwrk && ! lquery) {
#line 576 "dgesdd.f"
	    *info = -12;
#line 577 "dgesdd.f"
	}
#line 578 "dgesdd.f"
    }

#line 580 "dgesdd.f"
    if (*info != 0) {
#line 581 "dgesdd.f"
	i__1 = -(*info);
#line 581 "dgesdd.f"
	xerbla_("DGESDD", &i__1, (ftnlen)6);
#line 582 "dgesdd.f"
	return 0;
#line 583 "dgesdd.f"
    } else if (lquery) {
#line 584 "dgesdd.f"
	return 0;
#line 585 "dgesdd.f"
    }

/*     Quick return if possible */

#line 589 "dgesdd.f"
    if (*m == 0 || *n == 0) {
#line 590 "dgesdd.f"
	return 0;
#line 591 "dgesdd.f"
    }

/*     Get machine constants */

#line 595 "dgesdd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 596 "dgesdd.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 597 "dgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 601 "dgesdd.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 602 "dgesdd.f"
    iscl = 0;
#line 603 "dgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 604 "dgesdd.f"
	iscl = 1;
#line 605 "dgesdd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 606 "dgesdd.f"
    } else if (anrm > bignum) {
#line 607 "dgesdd.f"
	iscl = 1;
#line 608 "dgesdd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 609 "dgesdd.f"
    }

#line 611 "dgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 617 "dgesdd.f"
	if (*m >= mnthr) {

#line 619 "dgesdd.f"
	    if (wntqn) {

/*              Path 1 (M >> N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 624 "dgesdd.f"
		itau = 1;
#line 625 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              Workspace: need   N [tau] + N    [work] */
/*              Workspace: prefer N [tau] + N*NB [work] */

#line 631 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 631 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 636 "dgesdd.f"
		i__1 = *n - 1;
#line 636 "dgesdd.f"
		i__2 = *n - 1;
#line 636 "dgesdd.f"
		dlaset_("L", &i__1, &i__2, &c_b63, &c_b63, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 637 "dgesdd.f"
		ie = 1;
#line 638 "dgesdd.f"
		itauq = ie + *n;
#line 639 "dgesdd.f"
		itaup = itauq + *n;
#line 640 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              Workspace: need   3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 646 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 646 "dgesdd.f"
		dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 649 "dgesdd.f"
		nwork = ie + *n;

/*              Perform bidiagonal SVD, computing singular values only */
/*              Workspace: need   N [e] + BDSPAC */

#line 654 "dgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 657 "dgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M >> N, JOBZ = 'O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 663 "dgesdd.f"
		ir = 1;

/*              WORK(IR) is LDWRKR by N */

#line 667 "dgesdd.f"
		if (*lwork >= *lda * *n + *n * *n + *n * 3 + bdspac) {
#line 668 "dgesdd.f"
		    ldwrkr = *lda;
#line 669 "dgesdd.f"
		} else {
#line 670 "dgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3 - bdspac) / *n;
#line 671 "dgesdd.f"
		}
#line 672 "dgesdd.f"
		itau = ir + ldwrkr * *n;
#line 673 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 679 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 679 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 684 "dgesdd.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 685 "dgesdd.f"
		i__1 = *n - 1;
#line 685 "dgesdd.f"
		i__2 = *n - 1;
#line 685 "dgesdd.f"
		dlaset_("L", &i__1, &i__2, &c_b63, &c_b63, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 692 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 692 "dgesdd.f"
		dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 694 "dgesdd.f"
		ie = itau;
#line 695 "dgesdd.f"
		itauq = ie + *n;
#line 696 "dgesdd.f"
		itaup = itauq + *n;
#line 697 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 703 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 703 "dgesdd.f"
		dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              WORK(IU) is N by N */

#line 709 "dgesdd.f"
		iu = nwork;
#line 710 "dgesdd.f"
		nwork = iu + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + BDSPAC */

#line 717 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R */
/*              and VT by right singular vectors of R */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N    [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N*NB [work] */

#line 726 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 726 "dgesdd.f"
		dormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], n, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 729 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 729 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] */
/*              Workspace: prefer M*N [R] + 3*N [e, tauq, taup] + N*N [U] */

#line 738 "dgesdd.f"
		i__1 = *m;
#line 738 "dgesdd.f"
		i__2 = ldwrkr;
#line 738 "dgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 739 "dgesdd.f"
		    i__3 = *m - i__ + 1;
#line 739 "dgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 740 "dgesdd.f"
		    dgemm_("N", "N", &chunk, n, n, &c_b84, &a[i__ + a_dim1], 
			    lda, &work[iu], n, &c_b63, &work[ir], &ldwrkr, (
			    ftnlen)1, (ftnlen)1);
#line 743 "dgesdd.f"
		    dlacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 745 "dgesdd.f"
/* L10: */
#line 745 "dgesdd.f"
		}

#line 747 "dgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M >> N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 753 "dgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 757 "dgesdd.f"
		ldwrkr = *n;
#line 758 "dgesdd.f"
		itau = ir + ldwrkr * *n;
#line 759 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 765 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 765 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 770 "dgesdd.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 771 "dgesdd.f"
		i__2 = *n - 1;
#line 771 "dgesdd.f"
		i__1 = *n - 1;
#line 771 "dgesdd.f"
		dlaset_("L", &i__2, &i__1, &c_b63, &c_b63, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 778 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 778 "dgesdd.f"
		dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 780 "dgesdd.f"
		ie = itau;
#line 781 "dgesdd.f"
		itauq = ie + *n;
#line 782 "dgesdd.f"
		itaup = itauq + *n;
#line 783 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 789 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 789 "dgesdd.f"
		dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagoal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + BDSPAC */

#line 798 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N    [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*NB [work] */

#line 807 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 807 "dgesdd.f"
		dormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 811 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 811 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              Workspace: need   N*N [R] */

#line 819 "dgesdd.f"
		dlacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 820 "dgesdd.f"
		dgemm_("N", "N", m, n, n, &c_b84, &a[a_offset], lda, &work[ir]
			, &ldwrkr, &c_b63, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 823 "dgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M >> N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 829 "dgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 833 "dgesdd.f"
		ldwrku = *n;
#line 834 "dgesdd.f"
		itau = iu + ldwrku * *n;
#line 835 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              Workspace: need   N*N [U] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [U] + N [tau] + N*NB [work] */

#line 841 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 841 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 843 "dgesdd.f"
		dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              Workspace: need   N*N [U] + N [tau] + M    [work] */
/*              Workspace: prefer N*N [U] + N [tau] + M*NB [work] */
#line 848 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 848 "dgesdd.f"
		dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out other entries */

#line 853 "dgesdd.f"
		i__2 = *n - 1;
#line 853 "dgesdd.f"
		i__1 = *n - 1;
#line 853 "dgesdd.f"
		dlaset_("L", &i__2, &i__1, &c_b63, &c_b63, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 854 "dgesdd.f"
		ie = itau;
#line 855 "dgesdd.f"
		itauq = ie + *n;
#line 856 "dgesdd.f"
		itaup = itauq + *n;
#line 857 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 863 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 863 "dgesdd.f"
		dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + BDSPAC */

#line 872 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N    [work] */
/*              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + N*NB [work] */

#line 881 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 881 "dgesdd.f"
		dormbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 884 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 884 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              Workspace: need   N*N [U] */

#line 892 "dgesdd.f"
		dgemm_("N", "N", m, n, n, &c_b84, &u[u_offset], ldu, &work[iu]
			, &ldwrku, &c_b63, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 897 "dgesdd.f"
		dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 899 "dgesdd.f"
	    }

#line 901 "dgesdd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 5 (M >= N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 908 "dgesdd.f"
	    ie = 1;
#line 909 "dgesdd.f"
	    itauq = ie + *n;
#line 910 "dgesdd.f"
	    itaup = itauq + *n;
#line 911 "dgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           Workspace: need   3*N [e, tauq, taup] + M        [work] */
/*           Workspace: prefer 3*N [e, tauq, taup] + (M+N)*NB [work] */

#line 917 "dgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 917 "dgesdd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 920 "dgesdd.f"
	    if (wntqn) {

/*              Path 5n (M >= N, JOBZ='N') */
/*              Perform bidiagonal SVD, only computing singular values */
/*              Workspace: need   3*N [e, tauq, taup] + BDSPAC */

#line 926 "dgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 928 "dgesdd.f"
	    } else if (wntqo) {
/*              Path 5o (M >= N, JOBZ='O') */
#line 930 "dgesdd.f"
		iu = nwork;
#line 931 "dgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 WORK( IU ) is M by N */

#line 935 "dgesdd.f"
		    ldwrku = *m;
#line 936 "dgesdd.f"
		    nwork = iu + ldwrku * *n;
#line 937 "dgesdd.f"
		    dlaset_("F", m, n, &c_b63, &c_b63, &work[iu], &ldwrku, (
			    ftnlen)1);
/*                 IR is unused; silence compile warnings */
#line 940 "dgesdd.f"
		    ir = -1;
#line 941 "dgesdd.f"
		} else {

/*                 WORK( IU ) is N by N */

#line 945 "dgesdd.f"
		    ldwrku = *n;
#line 946 "dgesdd.f"
		    nwork = iu + ldwrku * *n;

/*                 WORK(IR) is LDWRKR by N */

#line 950 "dgesdd.f"
		    ir = nwork;
#line 951 "dgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 952 "dgesdd.f"
		}
#line 953 "dgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*N [e, tauq, taup] + N*N [U] + BDSPAC */

#line 960 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], &ldwrku, &
			vt[vt_offset], ldvt, dum, idum, &work[nwork], &iwork[
			1], info, (ftnlen)1, (ftnlen)1);

/*              Overwrite VT by right singular vectors of A */
/*              Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work] */

#line 968 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 968 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 972 "dgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 Path 5o-fast */
/*                 Overwrite WORK(IU) by left singular vectors of A */
/*                 Workspace: need   3*N [e, tauq, taup] + M*N [U] + N    [work] */
/*                 Workspace: prefer 3*N [e, tauq, taup] + M*N [U] + N*NB [work] */

#line 979 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 979 "dgesdd.f"
		    dormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy left singular vectors of A from WORK(IU) to A */

#line 985 "dgesdd.f"
		    dlacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 986 "dgesdd.f"
		} else {

/*                 Path 5o-slow */
/*                 Generate Q in A */
/*                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work] */
/*                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work] */

#line 993 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 993 "dgesdd.f"
		    dorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by left singular vectors of */
/*                 bidiagonal matrix in WORK(IU), storing result in */
/*                 WORK(IR) and copying to A */
/*                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + NB*N [R] */
/*                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + M*N  [R] */

#line 1002 "dgesdd.f"
		    i__2 = *m;
#line 1002 "dgesdd.f"
		    i__1 = ldwrkr;
#line 1002 "dgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1003 "dgesdd.f"
			i__3 = *m - i__ + 1;
#line 1003 "dgesdd.f"
			chunk = min(i__3,ldwrkr);
#line 1004 "dgesdd.f"
			dgemm_("N", "N", &chunk, n, n, &c_b84, &a[i__ + 
				a_dim1], lda, &work[iu], &ldwrku, &c_b63, &
				work[ir], &ldwrkr, (ftnlen)1, (ftnlen)1);
#line 1007 "dgesdd.f"
			dlacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 1009 "dgesdd.f"
/* L20: */
#line 1009 "dgesdd.f"
		    }
#line 1010 "dgesdd.f"
		}

#line 1012 "dgesdd.f"
	    } else if (wntqs) {

/*              Path 5s (M >= N, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*N [e, tauq, taup] + BDSPAC */

#line 1020 "dgesdd.f"
		dlaset_("F", m, n, &c_b63, &c_b63, &u[u_offset], ldu, (ftnlen)
			1);
#line 1021 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*N [e, tauq, taup] + N    [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + N*NB [work] */

#line 1030 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1030 "dgesdd.f"
		dormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1033 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1033 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1036 "dgesdd.f"
	    } else if (wntqa) {

/*              Path 5a (M >= N, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*N [e, tauq, taup] + BDSPAC */

#line 1044 "dgesdd.f"
		dlaset_("F", m, m, &c_b63, &c_b63, &u[u_offset], ldu, (ftnlen)
			1);
#line 1045 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 1051 "dgesdd.f"
		if (*m > *n) {
#line 1052 "dgesdd.f"
		    i__1 = *m - *n;
#line 1052 "dgesdd.f"
		    i__2 = *m - *n;
#line 1052 "dgesdd.f"
		    dlaset_("F", &i__1, &i__2, &c_b63, &c_b84, &u[*n + 1 + (*
			    n + 1) * u_dim1], ldu, (ftnlen)1);
#line 1054 "dgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*N [e, tauq, taup] + M    [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + M*NB [work] */

#line 1061 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1061 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1064 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1064 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1067 "dgesdd.f"
	    }

#line 1069 "dgesdd.f"
	}

#line 1071 "dgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 1077 "dgesdd.f"
	if (*n >= mnthr) {

#line 1079 "dgesdd.f"
	    if (wntqn) {

/*              Path 1t (N >> M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 1084 "dgesdd.f"
		itau = 1;
#line 1085 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              Workspace: need   M [tau] + M [work] */
/*              Workspace: prefer M [tau] + M*NB [work] */

#line 1091 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1091 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out above L */

#line 1096 "dgesdd.f"
		i__1 = *m - 1;
#line 1096 "dgesdd.f"
		i__2 = *m - 1;
#line 1096 "dgesdd.f"
		dlaset_("U", &i__1, &i__2, &c_b63, &c_b63, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 1097 "dgesdd.f"
		ie = 1;
#line 1098 "dgesdd.f"
		itauq = ie + *m;
#line 1099 "dgesdd.f"
		itaup = itauq + *m;
#line 1100 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              Workspace: need   3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1106 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1106 "dgesdd.f"
		dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 1109 "dgesdd.f"
		nwork = ie + *m;

/*              Perform bidiagonal SVD, computing singular values only */
/*              Workspace: need   M [e] + BDSPAC */

#line 1114 "dgesdd.f"
		dbdsdc_("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 1117 "dgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N >> M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1123 "dgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */
/*              WORK(IL)  is M by M; it is later resized to M by chunk for gemm */

#line 1128 "dgesdd.f"
		il = ivt + *m * *m;
#line 1129 "dgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3 + bdspac) {
#line 1130 "dgesdd.f"
		    ldwrkl = *m;
#line 1131 "dgesdd.f"
		    chunk = *n;
#line 1132 "dgesdd.f"
		} else {
#line 1133 "dgesdd.f"
		    ldwrkl = *m;
#line 1134 "dgesdd.f"
		    chunk = (*lwork - *m * *m) / *m;
#line 1135 "dgesdd.f"
		}
#line 1136 "dgesdd.f"
		itau = il + ldwrkl * *m;
#line 1137 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */

#line 1143 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1143 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1148 "dgesdd.f"
		dlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1149 "dgesdd.f"
		i__1 = *m - 1;
#line 1149 "dgesdd.f"
		i__2 = *m - 1;
#line 1149 "dgesdd.f"
		dlaset_("U", &i__1, &i__2, &c_b63, &c_b63, &work[il + ldwrkl],
			 &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */

#line 1156 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1156 "dgesdd.f"
		dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1158 "dgesdd.f"
		ie = itau;
#line 1159 "dgesdd.f"
		itauq = ie + *m;
#line 1160 "dgesdd.f"
		itaup = itauq + *m;
#line 1161 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1167 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1167 "dgesdd.f"
		dgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U, and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + BDSPAC */

#line 1176 "dgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], m, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M    [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M*NB [work] */

#line 1185 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1185 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1188 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1188 "dgesdd.f"
		dormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], m, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              Workspace: need   M*M [VT] + M*M [L] */
/*              Workspace: prefer M*M [VT] + M*N [L] */
/*              At this point, L is resized as M by chunk. */

#line 1198 "dgesdd.f"
		i__1 = *n;
#line 1198 "dgesdd.f"
		i__2 = chunk;
#line 1198 "dgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1199 "dgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1199 "dgesdd.f"
		    blk = min(i__3,chunk);
#line 1200 "dgesdd.f"
		    dgemm_("N", "N", m, &blk, m, &c_b84, &work[ivt], m, &a[
			    i__ * a_dim1 + 1], lda, &c_b63, &work[il], &
			    ldwrkl, (ftnlen)1, (ftnlen)1);
#line 1202 "dgesdd.f"
		    dlacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1204 "dgesdd.f"
/* L30: */
#line 1204 "dgesdd.f"
		}

#line 1206 "dgesdd.f"
	    } else if (wntqs) {

/*              Path 3t (N >> M, JOBZ='S') */
/*              M right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1212 "dgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1216 "dgesdd.f"
		ldwrkl = *m;
#line 1217 "dgesdd.f"
		itau = il + ldwrkl * *m;
#line 1218 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              Workspace: need   M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [L] + M [tau] + M*NB [work] */

#line 1224 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1224 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1229 "dgesdd.f"
		dlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1230 "dgesdd.f"
		i__2 = *m - 1;
#line 1230 "dgesdd.f"
		i__1 = *m - 1;
#line 1230 "dgesdd.f"
		dlaset_("U", &i__2, &i__1, &c_b63, &c_b63, &work[il + ldwrkl],
			 &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [L] + M [tau] + M*NB [work] */

#line 1237 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1237 "dgesdd.f"
		dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1239 "dgesdd.f"
		ie = itau;
#line 1240 "dgesdd.f"
		itauq = ie + *m;
#line 1241 "dgesdd.f"
		itaup = itauq + *m;
#line 1242 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IU). */
/*              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1248 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1248 "dgesdd.f"
		dgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + BDSPAC */

#line 1257 "dgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and VT */
/*              by right singular vectors of L */
/*              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M    [work] */
/*              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + M*NB [work] */

#line 1266 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1266 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1269 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1269 "dgesdd.f"
		dormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by */
/*              Q in A, storing result in VT */
/*              Workspace: need   M*M [L] */

#line 1277 "dgesdd.f"
		dlacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1278 "dgesdd.f"
		dgemm_("N", "N", m, n, m, &c_b84, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b63, &vt[vt_offset], ldvt, (ftnlen)
			1, (ftnlen)1);

#line 1281 "dgesdd.f"
	    } else if (wntqa) {

/*              Path 4t (N >> M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1287 "dgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1291 "dgesdd.f"
		ldwkvt = *m;
#line 1292 "dgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1293 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              Workspace: need   M*M [VT] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [VT] + M [tau] + M*NB [work] */

#line 1299 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1299 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 1301 "dgesdd.f"
		dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              Workspace: need   M*M [VT] + M [tau] + N    [work] */
/*              Workspace: prefer M*M [VT] + M [tau] + N*NB [work] */

#line 1307 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1307 "dgesdd.f"
		dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__2, &ierr);

/*              Produce L in A, zeroing out other entries */

#line 1312 "dgesdd.f"
		i__2 = *m - 1;
#line 1312 "dgesdd.f"
		i__1 = *m - 1;
#line 1312 "dgesdd.f"
		dlaset_("U", &i__2, &i__1, &c_b63, &c_b63, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 1313 "dgesdd.f"
		ie = itau;
#line 1314 "dgesdd.f"
		itauq = ie + *m;
#line 1315 "dgesdd.f"
		itaup = itauq + *m;
#line 1316 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1322 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1322 "dgesdd.f"
		dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + BDSPAC */

#line 1331 "dgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              Workspace: need   M*M [VT] + 3*M [e, tauq, taup]+ M    [work] */
/*              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup]+ M*NB [work] */

#line 1340 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1340 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1343 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1343 "dgesdd.f"
		dormbr_("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              Workspace: need   M*M [VT] */

#line 1351 "dgesdd.f"
		dgemm_("N", "N", m, n, m, &c_b84, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b63, &a[a_offset], lda, (ftnlen)
			1, (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1356 "dgesdd.f"
		dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1358 "dgesdd.f"
	    }

#line 1360 "dgesdd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 5t (N > M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 1367 "dgesdd.f"
	    ie = 1;
#line 1368 "dgesdd.f"
	    itauq = ie + *m;
#line 1369 "dgesdd.f"
	    itaup = itauq + *m;
#line 1370 "dgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           Workspace: need   3*M [e, tauq, taup] + N        [work] */
/*           Workspace: prefer 3*M [e, tauq, taup] + (M+N)*NB [work] */

#line 1376 "dgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1376 "dgesdd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 1379 "dgesdd.f"
	    if (wntqn) {

/*              Path 5tn (N > M, JOBZ='N') */
/*              Perform bidiagonal SVD, only computing singular values */
/*              Workspace: need   3*M [e, tauq, taup] + BDSPAC */

#line 1385 "dgesdd.f"
		dbdsdc_("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 1387 "dgesdd.f"
	    } else if (wntqo) {
/*              Path 5to (N > M, JOBZ='O') */
#line 1389 "dgesdd.f"
		ldwkvt = *m;
#line 1390 "dgesdd.f"
		ivt = nwork;
#line 1391 "dgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 WORK( IVT ) is M by N */

#line 1395 "dgesdd.f"
		    dlaset_("F", m, n, &c_b63, &c_b63, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 1397 "dgesdd.f"
		    nwork = ivt + ldwkvt * *n;
/*                 IL is unused; silence compile warnings */
#line 1399 "dgesdd.f"
		    il = -1;
#line 1400 "dgesdd.f"
		} else {

/*                 WORK( IVT ) is M by M */

#line 1404 "dgesdd.f"
		    nwork = ivt + ldwkvt * *m;
#line 1405 "dgesdd.f"
		    il = nwork;

/*                 WORK(IL) is M by CHUNK */

#line 1409 "dgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1410 "dgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + BDSPAC */

#line 1417 "dgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A */
/*              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work] */

#line 1425 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1425 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1429 "dgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 Path 5to-fast */
/*                 Overwrite WORK(IVT) by left singular vectors of A */
/*                 Workspace: need   3*M [e, tauq, taup] + M*N [VT] + M    [work] */
/*                 Workspace: prefer 3*M [e, tauq, taup] + M*N [VT] + M*NB [work] */

#line 1436 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1436 "dgesdd.f"
		    dormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy right singular vectors of A from WORK(IVT) to A */

#line 1442 "dgesdd.f"
		    dlacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 1443 "dgesdd.f"
		} else {

/*                 Path 5to-slow */
/*                 Generate P**T in A */
/*                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work] */
/*                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work] */

#line 1450 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1450 "dgesdd.f"
		    dorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by right singular vectors of */
/*                 bidiagonal matrix in WORK(IVT), storing result in */
/*                 WORK(IL) and copying to A */
/*                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M*NB [L] */
/*                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*N  [L] */

#line 1459 "dgesdd.f"
		    i__2 = *n;
#line 1459 "dgesdd.f"
		    i__1 = chunk;
#line 1459 "dgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1460 "dgesdd.f"
			i__3 = *n - i__ + 1;
#line 1460 "dgesdd.f"
			blk = min(i__3,chunk);
#line 1461 "dgesdd.f"
			dgemm_("N", "N", m, &blk, m, &c_b84, &work[ivt], &
				ldwkvt, &a[i__ * a_dim1 + 1], lda, &c_b63, &
				work[il], m, (ftnlen)1, (ftnlen)1);
#line 1464 "dgesdd.f"
			dlacpy_("F", m, &blk, &work[il], m, &a[i__ * a_dim1 + 
				1], lda, (ftnlen)1);
#line 1466 "dgesdd.f"
/* L40: */
#line 1466 "dgesdd.f"
		    }
#line 1467 "dgesdd.f"
		}
#line 1468 "dgesdd.f"
	    } else if (wntqs) {

/*              Path 5ts (N > M, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*M [e, tauq, taup] + BDSPAC */

#line 1476 "dgesdd.f"
		dlaset_("F", m, n, &c_b63, &c_b63, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1477 "dgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*M [e, tauq, taup] + M    [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + M*NB [work] */

#line 1486 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1486 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1489 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1489 "dgesdd.f"
		dormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1492 "dgesdd.f"
	    } else if (wntqa) {

/*              Path 5ta (N > M, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*M [e, tauq, taup] + BDSPAC */

#line 1500 "dgesdd.f"
		dlaset_("F", n, n, &c_b63, &c_b63, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1501 "dgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of VT to identity matrix */

#line 1507 "dgesdd.f"
		if (*n > *m) {
#line 1508 "dgesdd.f"
		    i__1 = *n - *m;
#line 1508 "dgesdd.f"
		    i__2 = *n - *m;
#line 1508 "dgesdd.f"
		    dlaset_("F", &i__1, &i__2, &c_b63, &c_b84, &vt[*m + 1 + (*
			    m + 1) * vt_dim1], ldvt, (ftnlen)1);
#line 1510 "dgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*M [e, tauq, taup] + N    [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + N*NB [work] */

#line 1517 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1517 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1520 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1520 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1523 "dgesdd.f"
	    }

#line 1525 "dgesdd.f"
	}

#line 1527 "dgesdd.f"
    }

/*     Undo scaling if necessary */

#line 1531 "dgesdd.f"
    if (iscl == 1) {
#line 1532 "dgesdd.f"
	if (anrm > bignum) {
#line 1532 "dgesdd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1532 "dgesdd.f"
	}
#line 1535 "dgesdd.f"
	if (anrm < smlnum) {
#line 1535 "dgesdd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1535 "dgesdd.f"
	}
#line 1538 "dgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 1542 "dgesdd.f"
    work[1] = (doublereal) maxwrk;

#line 1544 "dgesdd.f"
    return 0;

/*     End of DGESDD */

} /* dgesdd_ */

