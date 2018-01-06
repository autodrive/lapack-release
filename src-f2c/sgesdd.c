#line 1 "sgesdd.f"
/* sgesdd.f -- translated by f2c (version 20100827).
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

#line 1 "sgesdd.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b63 = 0.;
static integer c__1 = 1;
static doublereal c_b84 = 1.;

/* > \brief \b SGESDD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGESDD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesdd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesdd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesdd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL   A( LDA, * ), S( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGESDD computes the singular value decomposition (SVD) of a real */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          S is REAL array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array, dimension (LDU,UCOL) */
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
/* >          VT is REAL array, dimension (LDVT,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* >          > 0:  SBDSDC did not converge, updating process failed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup realGEsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgesdd_(char *jobz, integer *m, integer *n, doublereal *
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
    static integer lwork_sgelqf_mn__, lwork_sgeqrf_mn__, lwork_sorglq_mn__, 
	    lwork_sorglq_nn__, lwork_sorgqr_mm__, lwork_sorgqr_mn__, i__, ie, 
	    il, ir, iu, lwork_sorgbr_p_mm__, lwork_sorgbr_q_nn__, blk;
    static doublereal dum[1], eps;
    static integer ivt, iscl;
    static doublereal anrm;
    static integer idum[1], ierr, itau, lwork_sormbr_qln_mm__, 
	    lwork_sormbr_qln_mn__, lwork_sormbr_qln_nn__, 
	    lwork_sormbr_prt_mm__, lwork_sormbr_prt_mn__, 
	    lwork_sormbr_prt_nn__;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer minmn, wrkbl, itaup, itauq, mnthr;
    static logical wntqa;
    static integer nwork;
    static logical wntqn, wntqo, wntqs;
    static integer bdspac;
    extern /* Subroutine */ int sbdsdc_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), sgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int sgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    slascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     sgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), slacpy_(char *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    sorgbr_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static integer ldwrkl;
    extern /* Subroutine */ int sormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer ldwrkr, minwrk, ldwrku, maxwrk;
    extern /* Subroutine */ int sorglq_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static integer ldwkvt;
    static doublereal smlnum;
    static logical wntqas;
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical lquery;
    static integer lwork_sgebrd_mm__, lwork_sgebrd_mn__, lwork_sgebrd_nn__;


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

#line 281 "sgesdd.f"
    /* Parameter adjustments */
#line 281 "sgesdd.f"
    a_dim1 = *lda;
#line 281 "sgesdd.f"
    a_offset = 1 + a_dim1;
#line 281 "sgesdd.f"
    a -= a_offset;
#line 281 "sgesdd.f"
    --s;
#line 281 "sgesdd.f"
    u_dim1 = *ldu;
#line 281 "sgesdd.f"
    u_offset = 1 + u_dim1;
#line 281 "sgesdd.f"
    u -= u_offset;
#line 281 "sgesdd.f"
    vt_dim1 = *ldvt;
#line 281 "sgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 281 "sgesdd.f"
    vt -= vt_offset;
#line 281 "sgesdd.f"
    --work;
#line 281 "sgesdd.f"
    --iwork;
#line 281 "sgesdd.f"

#line 281 "sgesdd.f"
    /* Function Body */
#line 281 "sgesdd.f"
    *info = 0;
#line 282 "sgesdd.f"
    minmn = min(*m,*n);
#line 283 "sgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 284 "sgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 285 "sgesdd.f"
    wntqas = wntqa || wntqs;
#line 286 "sgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 287 "sgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 288 "sgesdd.f"
    lquery = *lwork == -1;

#line 290 "sgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 291 "sgesdd.f"
	*info = -1;
#line 292 "sgesdd.f"
    } else if (*m < 0) {
#line 293 "sgesdd.f"
	*info = -2;
#line 294 "sgesdd.f"
    } else if (*n < 0) {
#line 295 "sgesdd.f"
	*info = -3;
#line 296 "sgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 297 "sgesdd.f"
	*info = -5;
#line 298 "sgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 300 "sgesdd.f"
	*info = -8;
#line 301 "sgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 304 "sgesdd.f"
	*info = -10;
#line 305 "sgesdd.f"
    }

/*     Compute workspace */
/*       Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace allocated at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */

#line 314 "sgesdd.f"
    if (*info == 0) {
#line 315 "sgesdd.f"
	minwrk = 1;
#line 316 "sgesdd.f"
	maxwrk = 1;
#line 317 "sgesdd.f"
	bdspac = 0;
#line 318 "sgesdd.f"
	mnthr = (integer) (minmn * 11. / 6.);
#line 319 "sgesdd.f"
	if (*m >= *n && minmn > 0) {

/*           Compute space needed for SBDSDC */

#line 323 "sgesdd.f"
	    if (wntqn) {
/*              sbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6) */
/*              keep 7*N for backwards compatability. */
#line 326 "sgesdd.f"
		bdspac = *n * 7;
#line 327 "sgesdd.f"
	    } else {
#line 328 "sgesdd.f"
		bdspac = *n * 3 * *n + (*n << 2);
#line 329 "sgesdd.f"
	    }

/*           Compute space preferred for each routine */
#line 332 "sgesdd.f"
	    sgebrd_(m, n, dum, m, dum, dum, dum, dum, dum, &c_n1, &ierr);
#line 334 "sgesdd.f"
	    lwork_sgebrd_mn__ = (integer) dum[0];

#line 336 "sgesdd.f"
	    sgebrd_(n, n, dum, n, dum, dum, dum, dum, dum, &c_n1, &ierr);
#line 338 "sgesdd.f"
	    lwork_sgebrd_nn__ = (integer) dum[0];

#line 340 "sgesdd.f"
	    sgeqrf_(m, n, dum, m, dum, dum, &c_n1, &ierr);
#line 341 "sgesdd.f"
	    lwork_sgeqrf_mn__ = (integer) dum[0];

#line 343 "sgesdd.f"
	    sorgbr_("Q", n, n, n, dum, n, dum, dum, &c_n1, &ierr, (ftnlen)1);
#line 345 "sgesdd.f"
	    lwork_sorgbr_q_nn__ = (integer) dum[0];

#line 347 "sgesdd.f"
	    sorgqr_(m, m, n, dum, m, dum, dum, &c_n1, &ierr);
#line 348 "sgesdd.f"
	    lwork_sorgqr_mm__ = (integer) dum[0];

#line 350 "sgesdd.f"
	    sorgqr_(m, n, n, dum, m, dum, dum, &c_n1, &ierr);
#line 351 "sgesdd.f"
	    lwork_sorgqr_mn__ = (integer) dum[0];

#line 353 "sgesdd.f"
	    sormbr_("P", "R", "T", n, n, n, dum, n, dum, dum, n, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 355 "sgesdd.f"
	    lwork_sormbr_prt_nn__ = (integer) dum[0];

#line 357 "sgesdd.f"
	    sormbr_("Q", "L", "N", n, n, n, dum, n, dum, dum, n, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 359 "sgesdd.f"
	    lwork_sormbr_qln_nn__ = (integer) dum[0];

#line 361 "sgesdd.f"
	    sormbr_("Q", "L", "N", m, n, n, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 363 "sgesdd.f"
	    lwork_sormbr_qln_mn__ = (integer) dum[0];

#line 365 "sgesdd.f"
	    sormbr_("Q", "L", "N", m, m, n, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 367 "sgesdd.f"
	    lwork_sormbr_qln_mm__ = (integer) dum[0];

#line 369 "sgesdd.f"
	    if (*m >= mnthr) {
#line 370 "sgesdd.f"
		if (wntqn) {

/*                 Path 1 (M >> N, JOBZ='N') */

#line 374 "sgesdd.f"
		    wrkbl = *n + lwork_sgeqrf_mn__;
/* Computing MAX */
#line 375 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sgebrd_nn__;
#line 375 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 376 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n;
#line 376 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 377 "sgesdd.f"
		    minwrk = bdspac + *n;
#line 378 "sgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M >> N, JOBZ='O') */

#line 382 "sgesdd.f"
		    wrkbl = *n + lwork_sgeqrf_mn__;
/* Computing MAX */
#line 383 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_sorgqr_mn__;
#line 383 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 384 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sgebrd_nn__;
#line 384 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 385 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_qln_nn__;
#line 385 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 386 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_prt_nn__;
#line 386 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 387 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 387 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 388 "sgesdd.f"
		    maxwrk = wrkbl + (*n << 1) * *n;
#line 389 "sgesdd.f"
		    minwrk = bdspac + (*n << 1) * *n + *n * 3;
#line 390 "sgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M >> N, JOBZ='S') */

#line 394 "sgesdd.f"
		    wrkbl = *n + lwork_sgeqrf_mn__;
/* Computing MAX */
#line 395 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_sorgqr_mn__;
#line 395 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 396 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sgebrd_nn__;
#line 396 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 397 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_qln_nn__;
#line 397 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 398 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_prt_nn__;
#line 398 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 399 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 399 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 400 "sgesdd.f"
		    maxwrk = wrkbl + *n * *n;
#line 401 "sgesdd.f"
		    minwrk = bdspac + *n * *n + *n * 3;
#line 402 "sgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M >> N, JOBZ='A') */

#line 406 "sgesdd.f"
		    wrkbl = *n + lwork_sgeqrf_mn__;
/* Computing MAX */
#line 407 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_sorgqr_mm__;
#line 407 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 408 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sgebrd_nn__;
#line 408 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 409 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_qln_nn__;
#line 409 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 410 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_prt_nn__;
#line 410 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 411 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 411 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 412 "sgesdd.f"
		    maxwrk = wrkbl + *n * *n;
/* Computing MAX */
#line 413 "sgesdd.f"
		    i__1 = *n * 3 + bdspac, i__2 = *n + *m;
#line 413 "sgesdd.f"
		    minwrk = *n * *n + max(i__1,i__2);
#line 414 "sgesdd.f"
		}
#line 415 "sgesdd.f"
	    } else {

/*              Path 5 (M >= N, but not much larger) */

#line 419 "sgesdd.f"
		wrkbl = *n * 3 + lwork_sgebrd_mn__;
#line 420 "sgesdd.f"
		if (wntqn) {
/*                 Path 5n (M >= N, jobz='N') */
/* Computing MAX */
#line 422 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 422 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 423 "sgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 424 "sgesdd.f"
		} else if (wntqo) {
/*                 Path 5o (M >= N, jobz='O') */
/* Computing MAX */
#line 426 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_prt_nn__;
#line 426 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 427 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_qln_mn__;
#line 427 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 428 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 428 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 429 "sgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 430 "sgesdd.f"
		    i__1 = *m, i__2 = *n * *n + bdspac;
#line 430 "sgesdd.f"
		    minwrk = *n * 3 + max(i__1,i__2);
#line 431 "sgesdd.f"
		} else if (wntqs) {
/*                 Path 5s (M >= N, jobz='S') */
/* Computing MAX */
#line 433 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_qln_mn__;
#line 433 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 434 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_prt_nn__;
#line 434 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 435 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 435 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 436 "sgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 437 "sgesdd.f"
		} else if (wntqa) {
/*                 Path 5a (M >= N, jobz='A') */
/* Computing MAX */
#line 439 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_qln_mm__;
#line 439 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 440 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + lwork_sormbr_prt_nn__;
#line 440 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 441 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
#line 441 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 442 "sgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 443 "sgesdd.f"
		}
#line 444 "sgesdd.f"
	    }
#line 445 "sgesdd.f"
	} else if (minmn > 0) {

/*           Compute space needed for SBDSDC */

#line 449 "sgesdd.f"
	    if (wntqn) {
/*              sbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6) */
/*              keep 7*N for backwards compatability. */
#line 452 "sgesdd.f"
		bdspac = *m * 7;
#line 453 "sgesdd.f"
	    } else {
#line 454 "sgesdd.f"
		bdspac = *m * 3 * *m + (*m << 2);
#line 455 "sgesdd.f"
	    }

/*           Compute space preferred for each routine */
#line 458 "sgesdd.f"
	    sgebrd_(m, n, dum, m, dum, dum, dum, dum, dum, &c_n1, &ierr);
#line 460 "sgesdd.f"
	    lwork_sgebrd_mn__ = (integer) dum[0];

#line 462 "sgesdd.f"
	    sgebrd_(m, m, &a[a_offset], m, &s[1], dum, dum, dum, dum, &c_n1, &
		    ierr);
#line 464 "sgesdd.f"
	    lwork_sgebrd_mm__ = (integer) dum[0];

#line 466 "sgesdd.f"
	    sgelqf_(m, n, &a[a_offset], m, dum, dum, &c_n1, &ierr);
#line 467 "sgesdd.f"
	    lwork_sgelqf_mn__ = (integer) dum[0];

#line 469 "sgesdd.f"
	    sorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
#line 470 "sgesdd.f"
	    lwork_sorglq_nn__ = (integer) dum[0];

#line 472 "sgesdd.f"
	    sorglq_(m, n, m, &a[a_offset], m, dum, dum, &c_n1, &ierr);
#line 473 "sgesdd.f"
	    lwork_sorglq_mn__ = (integer) dum[0];

#line 475 "sgesdd.f"
	    sorgbr_("P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 476 "sgesdd.f"
	    lwork_sorgbr_p_mm__ = (integer) dum[0];

#line 478 "sgesdd.f"
	    sormbr_("P", "R", "T", m, m, m, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 480 "sgesdd.f"
	    lwork_sormbr_prt_mm__ = (integer) dum[0];

#line 482 "sgesdd.f"
	    sormbr_("P", "R", "T", m, n, m, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 484 "sgesdd.f"
	    lwork_sormbr_prt_mn__ = (integer) dum[0];

#line 486 "sgesdd.f"
	    sormbr_("P", "R", "T", n, n, m, dum, n, dum, dum, n, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 488 "sgesdd.f"
	    lwork_sormbr_prt_nn__ = (integer) dum[0];

#line 490 "sgesdd.f"
	    sormbr_("Q", "L", "N", m, m, m, dum, m, dum, dum, m, dum, &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 492 "sgesdd.f"
	    lwork_sormbr_qln_mm__ = (integer) dum[0];

#line 494 "sgesdd.f"
	    if (*n >= mnthr) {
#line 495 "sgesdd.f"
		if (wntqn) {

/*                 Path 1t (N >> M, JOBZ='N') */

#line 499 "sgesdd.f"
		    wrkbl = *m + lwork_sgelqf_mn__;
/* Computing MAX */
#line 500 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sgebrd_mm__;
#line 500 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 501 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m;
#line 501 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 502 "sgesdd.f"
		    minwrk = bdspac + *m;
#line 503 "sgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N >> M, JOBZ='O') */

#line 507 "sgesdd.f"
		    wrkbl = *m + lwork_sgelqf_mn__;
/* Computing MAX */
#line 508 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_sorglq_mn__;
#line 508 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 509 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sgebrd_mm__;
#line 509 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 510 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_qln_mm__;
#line 510 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 511 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_prt_mm__;
#line 511 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 512 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 512 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 513 "sgesdd.f"
		    maxwrk = wrkbl + (*m << 1) * *m;
#line 514 "sgesdd.f"
		    minwrk = bdspac + (*m << 1) * *m + *m * 3;
#line 515 "sgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N >> M, JOBZ='S') */

#line 519 "sgesdd.f"
		    wrkbl = *m + lwork_sgelqf_mn__;
/* Computing MAX */
#line 520 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_sorglq_mn__;
#line 520 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 521 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sgebrd_mm__;
#line 521 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 522 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_qln_mm__;
#line 522 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 523 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_prt_mm__;
#line 523 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 524 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 524 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 525 "sgesdd.f"
		    maxwrk = wrkbl + *m * *m;
#line 526 "sgesdd.f"
		    minwrk = bdspac + *m * *m + *m * 3;
#line 527 "sgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N >> M, JOBZ='A') */

#line 531 "sgesdd.f"
		    wrkbl = *m + lwork_sgelqf_mn__;
/* Computing MAX */
#line 532 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_sorglq_nn__;
#line 532 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 533 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sgebrd_mm__;
#line 533 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 534 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_qln_mm__;
#line 534 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 535 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_prt_mm__;
#line 535 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 536 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 536 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 537 "sgesdd.f"
		    maxwrk = wrkbl + *m * *m;
/* Computing MAX */
#line 538 "sgesdd.f"
		    i__1 = *m * 3 + bdspac, i__2 = *m + *n;
#line 538 "sgesdd.f"
		    minwrk = *m * *m + max(i__1,i__2);
#line 539 "sgesdd.f"
		}
#line 540 "sgesdd.f"
	    } else {

/*              Path 5t (N > M, but not much larger) */

#line 544 "sgesdd.f"
		wrkbl = *m * 3 + lwork_sgebrd_mn__;
#line 545 "sgesdd.f"
		if (wntqn) {
/*                 Path 5tn (N > M, jobz='N') */
/* Computing MAX */
#line 547 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 547 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 548 "sgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 549 "sgesdd.f"
		} else if (wntqo) {
/*                 Path 5to (N > M, jobz='O') */
/* Computing MAX */
#line 551 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_qln_mm__;
#line 551 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 552 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_prt_mn__;
#line 552 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 553 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 553 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 554 "sgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 555 "sgesdd.f"
		    i__1 = *n, i__2 = *m * *m + bdspac;
#line 555 "sgesdd.f"
		    minwrk = *m * 3 + max(i__1,i__2);
#line 556 "sgesdd.f"
		} else if (wntqs) {
/*                 Path 5ts (N > M, jobz='S') */
/* Computing MAX */
#line 558 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_qln_mm__;
#line 558 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 559 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_prt_mn__;
#line 559 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 560 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 560 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 561 "sgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 562 "sgesdd.f"
		} else if (wntqa) {
/*                 Path 5ta (N > M, jobz='A') */
/* Computing MAX */
#line 564 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_qln_mm__;
#line 564 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 565 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + lwork_sormbr_prt_nn__;
#line 565 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 566 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
#line 566 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 567 "sgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 568 "sgesdd.f"
		}
#line 569 "sgesdd.f"
	    }
#line 570 "sgesdd.f"
	}
#line 572 "sgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 573 "sgesdd.f"
	work[1] = (doublereal) maxwrk;

#line 575 "sgesdd.f"
	if (*lwork < minwrk && ! lquery) {
#line 576 "sgesdd.f"
	    *info = -12;
#line 577 "sgesdd.f"
	}
#line 578 "sgesdd.f"
    }

#line 580 "sgesdd.f"
    if (*info != 0) {
#line 581 "sgesdd.f"
	i__1 = -(*info);
#line 581 "sgesdd.f"
	xerbla_("SGESDD", &i__1, (ftnlen)6);
#line 582 "sgesdd.f"
	return 0;
#line 583 "sgesdd.f"
    } else if (lquery) {
#line 584 "sgesdd.f"
	return 0;
#line 585 "sgesdd.f"
    }

/*     Quick return if possible */

#line 589 "sgesdd.f"
    if (*m == 0 || *n == 0) {
#line 590 "sgesdd.f"
	return 0;
#line 591 "sgesdd.f"
    }

/*     Get machine constants */

#line 595 "sgesdd.f"
    eps = slamch_("P", (ftnlen)1);
#line 596 "sgesdd.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 597 "sgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 601 "sgesdd.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 602 "sgesdd.f"
    iscl = 0;
#line 603 "sgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 604 "sgesdd.f"
	iscl = 1;
#line 605 "sgesdd.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 606 "sgesdd.f"
    } else if (anrm > bignum) {
#line 607 "sgesdd.f"
	iscl = 1;
#line 608 "sgesdd.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 609 "sgesdd.f"
    }

#line 611 "sgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 617 "sgesdd.f"
	if (*m >= mnthr) {

#line 619 "sgesdd.f"
	    if (wntqn) {

/*              Path 1 (M >> N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 624 "sgesdd.f"
		itau = 1;
#line 625 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              Workspace: need   N [tau] + N    [work] */
/*              Workspace: prefer N [tau] + N*NB [work] */

#line 631 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 631 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 636 "sgesdd.f"
		i__1 = *n - 1;
#line 636 "sgesdd.f"
		i__2 = *n - 1;
#line 636 "sgesdd.f"
		slaset_("L", &i__1, &i__2, &c_b63, &c_b63, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 637 "sgesdd.f"
		ie = 1;
#line 638 "sgesdd.f"
		itauq = ie + *n;
#line 639 "sgesdd.f"
		itaup = itauq + *n;
#line 640 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              Workspace: need   3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 646 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 646 "sgesdd.f"
		sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 649 "sgesdd.f"
		nwork = ie + *n;

/*              Perform bidiagonal SVD, computing singular values only */
/*              Workspace: need   N [e] + BDSPAC */

#line 654 "sgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 657 "sgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M >> N, JOBZ = 'O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 663 "sgesdd.f"
		ir = 1;

/*              WORK(IR) is LDWRKR by N */

#line 667 "sgesdd.f"
		if (*lwork >= *lda * *n + *n * *n + *n * 3 + bdspac) {
#line 668 "sgesdd.f"
		    ldwrkr = *lda;
#line 669 "sgesdd.f"
		} else {
#line 670 "sgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3 - bdspac) / *n;
#line 671 "sgesdd.f"
		}
#line 672 "sgesdd.f"
		itau = ir + ldwrkr * *n;
#line 673 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 679 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 679 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 684 "sgesdd.f"
		slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 685 "sgesdd.f"
		i__1 = *n - 1;
#line 685 "sgesdd.f"
		i__2 = *n - 1;
#line 685 "sgesdd.f"
		slaset_("L", &i__1, &i__2, &c_b63, &c_b63, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 692 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 692 "sgesdd.f"
		sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 694 "sgesdd.f"
		ie = itau;
#line 695 "sgesdd.f"
		itauq = ie + *n;
#line 696 "sgesdd.f"
		itaup = itauq + *n;
#line 697 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 703 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 703 "sgesdd.f"
		sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              WORK(IU) is N by N */

#line 709 "sgesdd.f"
		iu = nwork;
#line 710 "sgesdd.f"
		nwork = iu + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + BDSPAC */

#line 717 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R */
/*              and VT by right singular vectors of R */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N    [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N*NB [work] */

#line 726 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 726 "sgesdd.f"
		sormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], n, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 729 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 729 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] */
/*              Workspace: prefer M*N [R] + 3*N [e, tauq, taup] + N*N [U] */

#line 738 "sgesdd.f"
		i__1 = *m;
#line 738 "sgesdd.f"
		i__2 = ldwrkr;
#line 738 "sgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 739 "sgesdd.f"
		    i__3 = *m - i__ + 1;
#line 739 "sgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 740 "sgesdd.f"
		    sgemm_("N", "N", &chunk, n, n, &c_b84, &a[i__ + a_dim1], 
			    lda, &work[iu], n, &c_b63, &work[ir], &ldwrkr, (
			    ftnlen)1, (ftnlen)1);
#line 743 "sgesdd.f"
		    slacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 745 "sgesdd.f"
/* L10: */
#line 745 "sgesdd.f"
		}

#line 747 "sgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M >> N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 753 "sgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 757 "sgesdd.f"
		ldwrkr = *n;
#line 758 "sgesdd.f"
		itau = ir + ldwrkr * *n;
#line 759 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 765 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 765 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 770 "sgesdd.f"
		slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 771 "sgesdd.f"
		i__2 = *n - 1;
#line 771 "sgesdd.f"
		i__1 = *n - 1;
#line 771 "sgesdd.f"
		slaset_("L", &i__2, &i__1, &c_b63, &c_b63, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   N*N [R] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [R] + N [tau] + N*NB [work] */

#line 778 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 778 "sgesdd.f"
		sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 780 "sgesdd.f"
		ie = itau;
#line 781 "sgesdd.f"
		itauq = ie + *n;
#line 782 "sgesdd.f"
		itaup = itauq + *n;
#line 783 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 789 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 789 "sgesdd.f"
		sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagoal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + BDSPAC */

#line 798 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N    [work] */
/*              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*NB [work] */

#line 807 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 807 "sgesdd.f"
		sormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 811 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 811 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              Workspace: need   N*N [R] */

#line 819 "sgesdd.f"
		slacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 820 "sgesdd.f"
		sgemm_("N", "N", m, n, n, &c_b84, &a[a_offset], lda, &work[ir]
			, &ldwrkr, &c_b63, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 823 "sgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M >> N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 829 "sgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 833 "sgesdd.f"
		ldwrku = *n;
#line 834 "sgesdd.f"
		itau = iu + ldwrku * *n;
#line 835 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              Workspace: need   N*N [U] + N [tau] + N    [work] */
/*              Workspace: prefer N*N [U] + N [tau] + N*NB [work] */

#line 841 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 841 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 843 "sgesdd.f"
		slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              Workspace: need   N*N [U] + N [tau] + M    [work] */
/*              Workspace: prefer N*N [U] + N [tau] + M*NB [work] */
#line 848 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 848 "sgesdd.f"
		sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out other entries */

#line 853 "sgesdd.f"
		i__2 = *n - 1;
#line 853 "sgesdd.f"
		i__1 = *n - 1;
#line 853 "sgesdd.f"
		slaset_("L", &i__2, &i__1, &c_b63, &c_b63, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 854 "sgesdd.f"
		ie = itau;
#line 855 "sgesdd.f"
		itauq = ie + *n;
#line 856 "sgesdd.f"
		itaup = itauq + *n;
#line 857 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N      [work] */
/*              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + 2*N*NB [work] */

#line 863 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 863 "sgesdd.f"
		sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + BDSPAC */

#line 872 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N    [work] */
/*              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + N*NB [work] */

#line 881 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 881 "sgesdd.f"
		sormbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 884 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 884 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              Workspace: need   N*N [U] */

#line 892 "sgesdd.f"
		sgemm_("N", "N", m, n, n, &c_b84, &u[u_offset], ldu, &work[iu]
			, &ldwrku, &c_b63, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 897 "sgesdd.f"
		slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 899 "sgesdd.f"
	    }

#line 901 "sgesdd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 5 (M >= N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 908 "sgesdd.f"
	    ie = 1;
#line 909 "sgesdd.f"
	    itauq = ie + *n;
#line 910 "sgesdd.f"
	    itaup = itauq + *n;
#line 911 "sgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           Workspace: need   3*N [e, tauq, taup] + M        [work] */
/*           Workspace: prefer 3*N [e, tauq, taup] + (M+N)*NB [work] */

#line 917 "sgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 917 "sgesdd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 920 "sgesdd.f"
	    if (wntqn) {

/*              Path 5n (M >= N, JOBZ='N') */
/*              Perform bidiagonal SVD, only computing singular values */
/*              Workspace: need   3*N [e, tauq, taup] + BDSPAC */

#line 926 "sgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 928 "sgesdd.f"
	    } else if (wntqo) {
/*              Path 5o (M >= N, JOBZ='O') */
#line 930 "sgesdd.f"
		iu = nwork;
#line 931 "sgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 WORK( IU ) is M by N */

#line 935 "sgesdd.f"
		    ldwrku = *m;
#line 936 "sgesdd.f"
		    nwork = iu + ldwrku * *n;
#line 937 "sgesdd.f"
		    slaset_("F", m, n, &c_b63, &c_b63, &work[iu], &ldwrku, (
			    ftnlen)1);
/*                 IR is unused; silence compile warnings */
#line 940 "sgesdd.f"
		    ir = -1;
#line 941 "sgesdd.f"
		} else {

/*                 WORK( IU ) is N by N */

#line 945 "sgesdd.f"
		    ldwrku = *n;
#line 946 "sgesdd.f"
		    nwork = iu + ldwrku * *n;

/*                 WORK(IR) is LDWRKR by N */

#line 950 "sgesdd.f"
		    ir = nwork;
#line 951 "sgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 952 "sgesdd.f"
		}
#line 953 "sgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*N [e, tauq, taup] + N*N [U] + BDSPAC */

#line 960 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], &ldwrku, &
			vt[vt_offset], ldvt, dum, idum, &work[nwork], &iwork[
			1], info, (ftnlen)1, (ftnlen)1);

/*              Overwrite VT by right singular vectors of A */
/*              Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work] */

#line 968 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 968 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 972 "sgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 Path 5o-fast */
/*                 Overwrite WORK(IU) by left singular vectors of A */
/*                 Workspace: need   3*N [e, tauq, taup] + M*N [U] + N    [work] */
/*                 Workspace: prefer 3*N [e, tauq, taup] + M*N [U] + N*NB [work] */

#line 979 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 979 "sgesdd.f"
		    sormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy left singular vectors of A from WORK(IU) to A */

#line 985 "sgesdd.f"
		    slacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 986 "sgesdd.f"
		} else {

/*                 Path 5o-slow */
/*                 Generate Q in A */
/*                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work] */
/*                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work] */

#line 993 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 993 "sgesdd.f"
		    sorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by left singular vectors of */
/*                 bidiagonal matrix in WORK(IU), storing result in */
/*                 WORK(IR) and copying to A */
/*                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + NB*N [R] */
/*                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + M*N  [R] */

#line 1002 "sgesdd.f"
		    i__2 = *m;
#line 1002 "sgesdd.f"
		    i__1 = ldwrkr;
#line 1002 "sgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1003 "sgesdd.f"
			i__3 = *m - i__ + 1;
#line 1003 "sgesdd.f"
			chunk = min(i__3,ldwrkr);
#line 1004 "sgesdd.f"
			sgemm_("N", "N", &chunk, n, n, &c_b84, &a[i__ + 
				a_dim1], lda, &work[iu], &ldwrku, &c_b63, &
				work[ir], &ldwrkr, (ftnlen)1, (ftnlen)1);
#line 1007 "sgesdd.f"
			slacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 1009 "sgesdd.f"
/* L20: */
#line 1009 "sgesdd.f"
		    }
#line 1010 "sgesdd.f"
		}

#line 1012 "sgesdd.f"
	    } else if (wntqs) {

/*              Path 5s (M >= N, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*N [e, tauq, taup] + BDSPAC */

#line 1020 "sgesdd.f"
		slaset_("F", m, n, &c_b63, &c_b63, &u[u_offset], ldu, (ftnlen)
			1);
#line 1021 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*N [e, tauq, taup] + N    [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + N*NB [work] */

#line 1030 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1030 "sgesdd.f"
		sormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1033 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1033 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1036 "sgesdd.f"
	    } else if (wntqa) {

/*              Path 5a (M >= N, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*N [e, tauq, taup] + BDSPAC */

#line 1044 "sgesdd.f"
		slaset_("F", m, m, &c_b63, &c_b63, &u[u_offset], ldu, (ftnlen)
			1);
#line 1045 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 1051 "sgesdd.f"
		if (*m > *n) {
#line 1052 "sgesdd.f"
		    i__1 = *m - *n;
#line 1052 "sgesdd.f"
		    i__2 = *m - *n;
#line 1052 "sgesdd.f"
		    slaset_("F", &i__1, &i__2, &c_b63, &c_b84, &u[*n + 1 + (*
			    n + 1) * u_dim1], ldu, (ftnlen)1);
#line 1054 "sgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*N [e, tauq, taup] + M    [work] */
/*              Workspace: prefer 3*N [e, tauq, taup] + M*NB [work] */

#line 1061 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1061 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1064 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1064 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1067 "sgesdd.f"
	    }

#line 1069 "sgesdd.f"
	}

#line 1071 "sgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 1077 "sgesdd.f"
	if (*n >= mnthr) {

#line 1079 "sgesdd.f"
	    if (wntqn) {

/*              Path 1t (N >> M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 1084 "sgesdd.f"
		itau = 1;
#line 1085 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              Workspace: need   M [tau] + M [work] */
/*              Workspace: prefer M [tau] + M*NB [work] */

#line 1091 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1091 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out above L */

#line 1096 "sgesdd.f"
		i__1 = *m - 1;
#line 1096 "sgesdd.f"
		i__2 = *m - 1;
#line 1096 "sgesdd.f"
		slaset_("U", &i__1, &i__2, &c_b63, &c_b63, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 1097 "sgesdd.f"
		ie = 1;
#line 1098 "sgesdd.f"
		itauq = ie + *m;
#line 1099 "sgesdd.f"
		itaup = itauq + *m;
#line 1100 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              Workspace: need   3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1106 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1106 "sgesdd.f"
		sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 1109 "sgesdd.f"
		nwork = ie + *m;

/*              Perform bidiagonal SVD, computing singular values only */
/*              Workspace: need   M [e] + BDSPAC */

#line 1114 "sgesdd.f"
		sbdsdc_("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 1117 "sgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N >> M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1123 "sgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */
/*              WORK(IL)  is M by M; it is later resized to M by chunk for gemm */

#line 1128 "sgesdd.f"
		il = ivt + *m * *m;
#line 1129 "sgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3 + bdspac) {
#line 1130 "sgesdd.f"
		    ldwrkl = *m;
#line 1131 "sgesdd.f"
		    chunk = *n;
#line 1132 "sgesdd.f"
		} else {
#line 1133 "sgesdd.f"
		    ldwrkl = *m;
#line 1134 "sgesdd.f"
		    chunk = (*lwork - *m * *m) / *m;
#line 1135 "sgesdd.f"
		}
#line 1136 "sgesdd.f"
		itau = il + ldwrkl * *m;
#line 1137 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */

#line 1143 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1143 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1148 "sgesdd.f"
		slacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1149 "sgesdd.f"
		i__1 = *m - 1;
#line 1149 "sgesdd.f"
		i__2 = *m - 1;
#line 1149 "sgesdd.f"
		slaset_("U", &i__1, &i__2, &c_b63, &c_b63, &work[il + ldwrkl],
			 &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */

#line 1156 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1156 "sgesdd.f"
		sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1158 "sgesdd.f"
		ie = itau;
#line 1159 "sgesdd.f"
		itauq = ie + *m;
#line 1160 "sgesdd.f"
		itaup = itauq + *m;
#line 1161 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1167 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1167 "sgesdd.f"
		sgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U, and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + BDSPAC */

#line 1176 "sgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], m, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M    [work] */
/*              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M*NB [work] */

#line 1185 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1185 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1188 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1188 "sgesdd.f"
		sormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], m, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              Workspace: need   M*M [VT] + M*M [L] */
/*              Workspace: prefer M*M [VT] + M*N [L] */
/*              At this point, L is resized as M by chunk. */

#line 1198 "sgesdd.f"
		i__1 = *n;
#line 1198 "sgesdd.f"
		i__2 = chunk;
#line 1198 "sgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1199 "sgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1199 "sgesdd.f"
		    blk = min(i__3,chunk);
#line 1200 "sgesdd.f"
		    sgemm_("N", "N", m, &blk, m, &c_b84, &work[ivt], m, &a[
			    i__ * a_dim1 + 1], lda, &c_b63, &work[il], &
			    ldwrkl, (ftnlen)1, (ftnlen)1);
#line 1202 "sgesdd.f"
		    slacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1204 "sgesdd.f"
/* L30: */
#line 1204 "sgesdd.f"
		}

#line 1206 "sgesdd.f"
	    } else if (wntqs) {

/*              Path 3t (N >> M, JOBZ='S') */
/*              M right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1212 "sgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1216 "sgesdd.f"
		ldwrkl = *m;
#line 1217 "sgesdd.f"
		itau = il + ldwrkl * *m;
#line 1218 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              Workspace: need   M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [L] + M [tau] + M*NB [work] */

#line 1224 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1224 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1229 "sgesdd.f"
		slacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1230 "sgesdd.f"
		i__2 = *m - 1;
#line 1230 "sgesdd.f"
		i__1 = *m - 1;
#line 1230 "sgesdd.f"
		slaset_("U", &i__2, &i__1, &c_b63, &c_b63, &work[il + ldwrkl],
			 &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              Workspace: need   M*M [L] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [L] + M [tau] + M*NB [work] */

#line 1237 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1237 "sgesdd.f"
		sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1239 "sgesdd.f"
		ie = itau;
#line 1240 "sgesdd.f"
		itauq = ie + *m;
#line 1241 "sgesdd.f"
		itaup = itauq + *m;
#line 1242 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IU). */
/*              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1248 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1248 "sgesdd.f"
		sgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + BDSPAC */

#line 1257 "sgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and VT */
/*              by right singular vectors of L */
/*              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M    [work] */
/*              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + M*NB [work] */

#line 1266 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1266 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1269 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1269 "sgesdd.f"
		sormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by */
/*              Q in A, storing result in VT */
/*              Workspace: need   M*M [L] */

#line 1277 "sgesdd.f"
		slacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1278 "sgesdd.f"
		sgemm_("N", "N", m, n, m, &c_b84, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b63, &vt[vt_offset], ldvt, (ftnlen)
			1, (ftnlen)1);

#line 1281 "sgesdd.f"
	    } else if (wntqa) {

/*              Path 4t (N >> M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1287 "sgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1291 "sgesdd.f"
		ldwkvt = *m;
#line 1292 "sgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1293 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              Workspace: need   M*M [VT] + M [tau] + M    [work] */
/*              Workspace: prefer M*M [VT] + M [tau] + M*NB [work] */

#line 1299 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1299 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 1301 "sgesdd.f"
		slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              Workspace: need   M*M [VT] + M [tau] + N    [work] */
/*              Workspace: prefer M*M [VT] + M [tau] + N*NB [work] */

#line 1307 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1307 "sgesdd.f"
		sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__2, &ierr);

/*              Produce L in A, zeroing out other entries */

#line 1312 "sgesdd.f"
		i__2 = *m - 1;
#line 1312 "sgesdd.f"
		i__1 = *m - 1;
#line 1312 "sgesdd.f"
		slaset_("U", &i__2, &i__1, &c_b63, &c_b63, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 1313 "sgesdd.f"
		ie = itau;
#line 1314 "sgesdd.f"
		itauq = ie + *m;
#line 1315 "sgesdd.f"
		itaup = itauq + *m;
#line 1316 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + M      [work] */
/*              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup] + 2*M*NB [work] */

#line 1322 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1322 "sgesdd.f"
		sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + BDSPAC */

#line 1331 "sgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              Workspace: need   M*M [VT] + 3*M [e, tauq, taup]+ M    [work] */
/*              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup]+ M*NB [work] */

#line 1340 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1340 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1343 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1343 "sgesdd.f"
		sormbr_("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              Workspace: need   M*M [VT] */

#line 1351 "sgesdd.f"
		sgemm_("N", "N", m, n, m, &c_b84, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b63, &a[a_offset], lda, (ftnlen)
			1, (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1356 "sgesdd.f"
		slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1358 "sgesdd.f"
	    }

#line 1360 "sgesdd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 5t (N > M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 1367 "sgesdd.f"
	    ie = 1;
#line 1368 "sgesdd.f"
	    itauq = ie + *m;
#line 1369 "sgesdd.f"
	    itaup = itauq + *m;
#line 1370 "sgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           Workspace: need   3*M [e, tauq, taup] + N        [work] */
/*           Workspace: prefer 3*M [e, tauq, taup] + (M+N)*NB [work] */

#line 1376 "sgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1376 "sgesdd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 1379 "sgesdd.f"
	    if (wntqn) {

/*              Path 5tn (N > M, JOBZ='N') */
/*              Perform bidiagonal SVD, only computing singular values */
/*              Workspace: need   3*M [e, tauq, taup] + BDSPAC */

#line 1385 "sgesdd.f"
		sbdsdc_("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 1387 "sgesdd.f"
	    } else if (wntqo) {
/*              Path 5to (N > M, JOBZ='O') */
#line 1389 "sgesdd.f"
		ldwkvt = *m;
#line 1390 "sgesdd.f"
		ivt = nwork;
#line 1391 "sgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 WORK( IVT ) is M by N */

#line 1395 "sgesdd.f"
		    slaset_("F", m, n, &c_b63, &c_b63, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 1397 "sgesdd.f"
		    nwork = ivt + ldwkvt * *n;
/*                 IL is unused; silence compile warnings */
#line 1399 "sgesdd.f"
		    il = -1;
#line 1400 "sgesdd.f"
		} else {

/*                 WORK( IVT ) is M by M */

#line 1404 "sgesdd.f"
		    nwork = ivt + ldwkvt * *m;
#line 1405 "sgesdd.f"
		    il = nwork;

/*                 WORK(IL) is M by CHUNK */

#line 1409 "sgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1410 "sgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + BDSPAC */

#line 1417 "sgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A */
/*              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work] */

#line 1425 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1425 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1429 "sgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 Path 5to-fast */
/*                 Overwrite WORK(IVT) by left singular vectors of A */
/*                 Workspace: need   3*M [e, tauq, taup] + M*N [VT] + M    [work] */
/*                 Workspace: prefer 3*M [e, tauq, taup] + M*N [VT] + M*NB [work] */

#line 1436 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1436 "sgesdd.f"
		    sormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy right singular vectors of A from WORK(IVT) to A */

#line 1442 "sgesdd.f"
		    slacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 1443 "sgesdd.f"
		} else {

/*                 Path 5to-slow */
/*                 Generate P**T in A */
/*                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work] */
/*                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work] */

#line 1450 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1450 "sgesdd.f"
		    sorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by right singular vectors of */
/*                 bidiagonal matrix in WORK(IVT), storing result in */
/*                 WORK(IL) and copying to A */
/*                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M*NB [L] */
/*                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*N  [L] */

#line 1459 "sgesdd.f"
		    i__2 = *n;
#line 1459 "sgesdd.f"
		    i__1 = chunk;
#line 1459 "sgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1460 "sgesdd.f"
			i__3 = *n - i__ + 1;
#line 1460 "sgesdd.f"
			blk = min(i__3,chunk);
#line 1461 "sgesdd.f"
			sgemm_("N", "N", m, &blk, m, &c_b84, &work[ivt], &
				ldwkvt, &a[i__ * a_dim1 + 1], lda, &c_b63, &
				work[il], m, (ftnlen)1, (ftnlen)1);
#line 1464 "sgesdd.f"
			slacpy_("F", m, &blk, &work[il], m, &a[i__ * a_dim1 + 
				1], lda, (ftnlen)1);
#line 1466 "sgesdd.f"
/* L40: */
#line 1466 "sgesdd.f"
		    }
#line 1467 "sgesdd.f"
		}
#line 1468 "sgesdd.f"
	    } else if (wntqs) {

/*              Path 5ts (N > M, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*M [e, tauq, taup] + BDSPAC */

#line 1476 "sgesdd.f"
		slaset_("F", m, n, &c_b63, &c_b63, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1477 "sgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*M [e, tauq, taup] + M    [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + M*NB [work] */

#line 1486 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1486 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1489 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1489 "sgesdd.f"
		sormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1492 "sgesdd.f"
	    } else if (wntqa) {

/*              Path 5ta (N > M, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              Workspace: need   3*M [e, tauq, taup] + BDSPAC */

#line 1500 "sgesdd.f"
		slaset_("F", n, n, &c_b63, &c_b63, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1501 "sgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of VT to identity matrix */

#line 1507 "sgesdd.f"
		if (*n > *m) {
#line 1508 "sgesdd.f"
		    i__1 = *n - *m;
#line 1508 "sgesdd.f"
		    i__2 = *n - *m;
#line 1508 "sgesdd.f"
		    slaset_("F", &i__1, &i__2, &c_b63, &c_b84, &vt[*m + 1 + (*
			    m + 1) * vt_dim1], ldvt, (ftnlen)1);
#line 1510 "sgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              Workspace: need   3*M [e, tauq, taup] + N    [work] */
/*              Workspace: prefer 3*M [e, tauq, taup] + N*NB [work] */

#line 1517 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1517 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1520 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1520 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1523 "sgesdd.f"
	    }

#line 1525 "sgesdd.f"
	}

#line 1527 "sgesdd.f"
    }

/*     Undo scaling if necessary */

#line 1531 "sgesdd.f"
    if (iscl == 1) {
#line 1532 "sgesdd.f"
	if (anrm > bignum) {
#line 1532 "sgesdd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1532 "sgesdd.f"
	}
#line 1535 "sgesdd.f"
	if (anrm < smlnum) {
#line 1535 "sgesdd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1535 "sgesdd.f"
	}
#line 1538 "sgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 1542 "sgesdd.f"
    work[1] = (doublereal) maxwrk;

#line 1544 "sgesdd.f"
    return 0;

/*     End of SGESDD */

} /* sgesdd_ */

