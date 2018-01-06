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

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b227 = 0.;
static doublereal c_b248 = 1.;

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

/*       SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, */
/*                          LWORK, IWORK, INFO ) */

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
/* >          The leading dimension of the array VT.  LDVT >= 1; if */
/* >          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N; */
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
/* >          If JOBZ = 'N', */
/* >            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)). */
/* >          If JOBZ = 'O', */
/* >            LWORK >= 3*min(M,N) + */
/* >                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)). */
/* >          If JOBZ = 'S' or 'A' */
/* >            LWORK >= min(M,N)*(7+4*min(M,N)) */
/* >          For good performance, LWORK should generally be larger. */
/* >          If LWORK = -1 but other input arguments are legal, WORK(1) */
/* >          returns the optimal LWORK. */
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

/* > \date November 2015 */

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
    static integer i__, ie, il, ir, iu, blk;
    static doublereal dum[1], eps;
    static integer ivt, iscl;
    static doublereal anrm;
    static integer idum[1], ierr, itau;
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
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 270 "dgesdd.f"
    /* Parameter adjustments */
#line 270 "dgesdd.f"
    a_dim1 = *lda;
#line 270 "dgesdd.f"
    a_offset = 1 + a_dim1;
#line 270 "dgesdd.f"
    a -= a_offset;
#line 270 "dgesdd.f"
    --s;
#line 270 "dgesdd.f"
    u_dim1 = *ldu;
#line 270 "dgesdd.f"
    u_offset = 1 + u_dim1;
#line 270 "dgesdd.f"
    u -= u_offset;
#line 270 "dgesdd.f"
    vt_dim1 = *ldvt;
#line 270 "dgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 270 "dgesdd.f"
    vt -= vt_offset;
#line 270 "dgesdd.f"
    --work;
#line 270 "dgesdd.f"
    --iwork;
#line 270 "dgesdd.f"

#line 270 "dgesdd.f"
    /* Function Body */
#line 270 "dgesdd.f"
    *info = 0;
#line 271 "dgesdd.f"
    minmn = min(*m,*n);
#line 272 "dgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 273 "dgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 274 "dgesdd.f"
    wntqas = wntqa || wntqs;
#line 275 "dgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 276 "dgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 277 "dgesdd.f"
    lquery = *lwork == -1;

#line 279 "dgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 280 "dgesdd.f"
	*info = -1;
#line 281 "dgesdd.f"
    } else if (*m < 0) {
#line 282 "dgesdd.f"
	*info = -2;
#line 283 "dgesdd.f"
    } else if (*n < 0) {
#line 284 "dgesdd.f"
	*info = -3;
#line 285 "dgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 286 "dgesdd.f"
	*info = -5;
#line 287 "dgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 289 "dgesdd.f"
	*info = -8;
#line 290 "dgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 293 "dgesdd.f"
	*info = -10;
#line 294 "dgesdd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 303 "dgesdd.f"
    if (*info == 0) {
#line 304 "dgesdd.f"
	minwrk = 1;
#line 305 "dgesdd.f"
	maxwrk = 1;
#line 306 "dgesdd.f"
	if (*m >= *n && minmn > 0) {

/*           Compute space needed for DBDSDC */

#line 310 "dgesdd.f"
	    mnthr = (integer) (minmn * 11. / 6.);
#line 311 "dgesdd.f"
	    if (wntqn) {
#line 312 "dgesdd.f"
		bdspac = *n * 7;
#line 313 "dgesdd.f"
	    } else {
#line 314 "dgesdd.f"
		bdspac = *n * 3 * *n + (*n << 2);
#line 315 "dgesdd.f"
	    }
#line 316 "dgesdd.f"
	    if (*m >= mnthr) {
#line 317 "dgesdd.f"
		if (wntqn) {

/*                 Path 1 (M much larger than N, JOBZ='N') */

#line 321 "dgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 323 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 323 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 325 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n;
#line 325 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 326 "dgesdd.f"
		    minwrk = bdspac + *n;
#line 327 "dgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M much larger than N, JOBZ='O') */

#line 331 "dgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 332 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 332 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 334 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 334 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 336 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "QLN", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 336 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 338 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 338 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 340 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 340 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 341 "dgesdd.f"
		    maxwrk = wrkbl + (*n << 1) * *n;
#line 342 "dgesdd.f"
		    minwrk = bdspac + (*n << 1) * *n + *n * 3;
#line 343 "dgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M much larger than N, JOBZ='S') */

#line 347 "dgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 348 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 348 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 350 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 350 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 352 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "QLN", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 352 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 354 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 354 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 356 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 356 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 357 "dgesdd.f"
		    maxwrk = wrkbl + *n * *n;
#line 358 "dgesdd.f"
		    minwrk = bdspac + *n * *n + *n * 3;
#line 359 "dgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M much larger than N, JOBZ='A') */

#line 363 "dgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 364 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *m * ilaenv_(&c__1, "DORGQR", 
			    " ", m, m, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 364 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 366 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 366 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 368 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "QLN", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 368 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 370 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 370 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 372 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 372 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 373 "dgesdd.f"
		    maxwrk = wrkbl + *n * *n;
#line 374 "dgesdd.f"
		    minwrk = bdspac + *n * *n + (*n << 1) + *m;
#line 375 "dgesdd.f"
		}
#line 376 "dgesdd.f"
	    } else {

/*              Path 5 (M at least N, but not much larger) */

#line 380 "dgesdd.f"
		wrkbl = *n * 3 + (*m + *n) * ilaenv_(&c__1, "DGEBRD", " ", m, 
			n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 382 "dgesdd.f"
		if (wntqn) {
/* Computing MAX */
#line 383 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 383 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 384 "dgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 385 "dgesdd.f"
		} else if (wntqo) {
/* Computing MAX */
#line 386 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 386 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 388 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 388 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 390 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 390 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 391 "dgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 392 "dgesdd.f"
		    i__1 = *m, i__2 = *n * *n + bdspac;
#line 392 "dgesdd.f"
		    minwrk = *n * 3 + max(i__1,i__2);
#line 393 "dgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 394 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 394 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 396 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 396 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 398 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 398 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 399 "dgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 400 "dgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 401 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 401 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 403 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "DORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 403 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 405 "dgesdd.f"
		    i__1 = maxwrk, i__2 = bdspac + *n * 3;
#line 405 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 406 "dgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 407 "dgesdd.f"
		}
#line 408 "dgesdd.f"
	    }
#line 409 "dgesdd.f"
	} else if (minmn > 0) {

/*           Compute space needed for DBDSDC */

#line 413 "dgesdd.f"
	    mnthr = (integer) (minmn * 11. / 6.);
#line 414 "dgesdd.f"
	    if (wntqn) {
#line 415 "dgesdd.f"
		bdspac = *m * 7;
#line 416 "dgesdd.f"
	    } else {
#line 417 "dgesdd.f"
		bdspac = *m * 3 * *m + (*m << 2);
#line 418 "dgesdd.f"
	    }
#line 419 "dgesdd.f"
	    if (*n >= mnthr) {
#line 420 "dgesdd.f"
		if (wntqn) {

/*                 Path 1t (N much larger than M, JOBZ='N') */

#line 424 "dgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 426 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 426 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 428 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m;
#line 428 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 429 "dgesdd.f"
		    minwrk = bdspac + *m;
#line 430 "dgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N much larger than M, JOBZ='O') */

#line 434 "dgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 435 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "DORGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 435 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 437 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 437 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 439 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 439 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 441 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "PRT", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 441 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 443 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 443 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 444 "dgesdd.f"
		    maxwrk = wrkbl + (*m << 1) * *m;
#line 445 "dgesdd.f"
		    minwrk = bdspac + (*m << 1) * *m + *m * 3;
#line 446 "dgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N much larger than M, JOBZ='S') */

#line 450 "dgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 451 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "DORGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 451 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 453 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 453 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 455 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 455 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 457 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "PRT", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 457 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 459 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 459 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 460 "dgesdd.f"
		    maxwrk = wrkbl + *m * *m;
#line 461 "dgesdd.f"
		    minwrk = bdspac + *m * *m + *m * 3;
#line 462 "dgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N much larger than M, JOBZ='A') */

#line 466 "dgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 467 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *n * ilaenv_(&c__1, "DORGLQ", 
			    " ", n, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 467 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 469 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "DGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 469 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 471 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 471 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 473 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "PRT", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 473 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 475 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 475 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 476 "dgesdd.f"
		    maxwrk = wrkbl + *m * *m;
#line 477 "dgesdd.f"
		    minwrk = bdspac + *m * *m + *m * 3;
#line 478 "dgesdd.f"
		}
#line 479 "dgesdd.f"
	    } else {

/*              Path 5t (N greater than M, but not much larger) */

#line 483 "dgesdd.f"
		wrkbl = *m * 3 + (*m + *n) * ilaenv_(&c__1, "DGEBRD", " ", m, 
			n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 485 "dgesdd.f"
		if (wntqn) {
/* Computing MAX */
#line 486 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 486 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 487 "dgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 488 "dgesdd.f"
		} else if (wntqo) {
/* Computing MAX */
#line 489 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 489 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 491 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "PRT", m, n, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 491 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 493 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 493 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 494 "dgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 495 "dgesdd.f"
		    i__1 = *n, i__2 = *m * *m + bdspac;
#line 495 "dgesdd.f"
		    minwrk = *m * 3 + max(i__1,i__2);
#line 496 "dgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 497 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 497 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 499 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "PRT", m, n, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 499 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 501 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 501 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 502 "dgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 503 "dgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 504 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 504 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 506 "dgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR"
			    , "PRT", n, n, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 506 "dgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 508 "dgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 508 "dgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 509 "dgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 510 "dgesdd.f"
		}
#line 511 "dgesdd.f"
	    }
#line 512 "dgesdd.f"
	}
#line 513 "dgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 514 "dgesdd.f"
	work[1] = (doublereal) maxwrk;

#line 516 "dgesdd.f"
	if (*lwork < minwrk && ! lquery) {
#line 517 "dgesdd.f"
	    *info = -12;
#line 518 "dgesdd.f"
	}
#line 519 "dgesdd.f"
    }

#line 521 "dgesdd.f"
    if (*info != 0) {
#line 522 "dgesdd.f"
	i__1 = -(*info);
#line 522 "dgesdd.f"
	xerbla_("DGESDD", &i__1, (ftnlen)6);
#line 523 "dgesdd.f"
	return 0;
#line 524 "dgesdd.f"
    } else if (lquery) {
#line 525 "dgesdd.f"
	return 0;
#line 526 "dgesdd.f"
    }

/*     Quick return if possible */

#line 530 "dgesdd.f"
    if (*m == 0 || *n == 0) {
#line 531 "dgesdd.f"
	return 0;
#line 532 "dgesdd.f"
    }

/*     Get machine constants */

#line 536 "dgesdd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 537 "dgesdd.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 538 "dgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 542 "dgesdd.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 543 "dgesdd.f"
    iscl = 0;
#line 544 "dgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 545 "dgesdd.f"
	iscl = 1;
#line 546 "dgesdd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 547 "dgesdd.f"
    } else if (anrm > bignum) {
#line 548 "dgesdd.f"
	iscl = 1;
#line 549 "dgesdd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 550 "dgesdd.f"
    }

#line 552 "dgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 558 "dgesdd.f"
	if (*m >= mnthr) {

#line 560 "dgesdd.f"
	    if (wntqn) {

/*              Path 1 (M much larger than N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 565 "dgesdd.f"
		itau = 1;
#line 566 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need 2*N, prefer N+N*NB) */

#line 571 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 571 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 576 "dgesdd.f"
		i__1 = *n - 1;
#line 576 "dgesdd.f"
		i__2 = *n - 1;
#line 576 "dgesdd.f"
		dlaset_("L", &i__1, &i__2, &c_b227, &c_b227, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 577 "dgesdd.f"
		ie = 1;
#line 578 "dgesdd.f"
		itauq = ie + *n;
#line 579 "dgesdd.f"
		itaup = itauq + *n;
#line 580 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 585 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 585 "dgesdd.f"
		dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 588 "dgesdd.f"
		nwork = ie + *n;

/*              Perform bidiagonal SVD, computing singular values only */
/*              (Workspace: need N+BDSPAC) */

#line 593 "dgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 596 "dgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M much larger than N, JOBZ = 'O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 602 "dgesdd.f"
		ir = 1;

/*              WORK(IR) is LDWRKR by N */

#line 606 "dgesdd.f"
		if (*lwork >= *lda * *n + *n * *n + *n * 3 + bdspac) {
#line 607 "dgesdd.f"
		    ldwrkr = *lda;
#line 608 "dgesdd.f"
		} else {
#line 609 "dgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3 - bdspac) / *n;
#line 610 "dgesdd.f"
		}
#line 611 "dgesdd.f"
		itau = ir + ldwrkr * *n;
#line 612 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 617 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 617 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 622 "dgesdd.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 623 "dgesdd.f"
		i__1 = *n - 1;
#line 623 "dgesdd.f"
		i__2 = *n - 1;
#line 623 "dgesdd.f"
		dlaset_("L", &i__1, &i__2, &c_b227, &c_b227, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 629 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 629 "dgesdd.f"
		dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 631 "dgesdd.f"
		ie = itau;
#line 632 "dgesdd.f"
		itauq = ie + *n;
#line 633 "dgesdd.f"
		itaup = itauq + *n;
#line 634 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in VT, copying result to WORK(IR) */
/*              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 639 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 639 "dgesdd.f"
		dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              WORK(IU) is N by N */

#line 645 "dgesdd.f"
		iu = nwork;
#line 646 "dgesdd.f"
		nwork = iu + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+N*N+BDSPAC) */

#line 653 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R */
/*              and VT by right singular vectors of R */
/*              (Workspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */

#line 661 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 661 "dgesdd.f"
		dormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], n, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 664 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 664 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              (Workspace: need 2*N*N, prefer N*N+M*N) */

#line 672 "dgesdd.f"
		i__1 = *m;
#line 672 "dgesdd.f"
		i__2 = ldwrkr;
#line 672 "dgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 673 "dgesdd.f"
		    i__3 = *m - i__ + 1;
#line 673 "dgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 674 "dgesdd.f"
		    dgemm_("N", "N", &chunk, n, n, &c_b248, &a[i__ + a_dim1], 
			    lda, &work[iu], n, &c_b227, &work[ir], &ldwrkr, (
			    ftnlen)1, (ftnlen)1);
#line 677 "dgesdd.f"
		    dlacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 679 "dgesdd.f"
/* L10: */
#line 679 "dgesdd.f"
		}

#line 681 "dgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M much larger than N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 687 "dgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 691 "dgesdd.f"
		ldwrkr = *n;
#line 692 "dgesdd.f"
		itau = ir + ldwrkr * *n;
#line 693 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 698 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 698 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 703 "dgesdd.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 704 "dgesdd.f"
		i__2 = *n - 1;
#line 704 "dgesdd.f"
		i__1 = *n - 1;
#line 704 "dgesdd.f"
		dlaset_("L", &i__2, &i__1, &c_b227, &c_b227, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 710 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 710 "dgesdd.f"
		dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 712 "dgesdd.f"
		ie = itau;
#line 713 "dgesdd.f"
		itauq = ie + *n;
#line 714 "dgesdd.f"
		itaup = itauq + *n;
#line 715 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 720 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 720 "dgesdd.f"
		dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagoal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+BDSPAC) */

#line 729 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB) */

#line 737 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 737 "dgesdd.f"
		dormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 741 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 741 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              (Workspace: need N*N) */

#line 749 "dgesdd.f"
		dlacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 750 "dgesdd.f"
		dgemm_("N", "N", m, n, n, &c_b248, &a[a_offset], lda, &work[
			ir], &ldwrkr, &c_b227, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 753 "dgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M much larger than N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 759 "dgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 763 "dgesdd.f"
		ldwrku = *n;
#line 764 "dgesdd.f"
		itau = iu + ldwrku * *n;
#line 765 "dgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 770 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 770 "dgesdd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 772 "dgesdd.f"
		dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */
#line 776 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 776 "dgesdd.f"
		dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out other entries */

#line 781 "dgesdd.f"
		i__2 = *n - 1;
#line 781 "dgesdd.f"
		i__1 = *n - 1;
#line 781 "dgesdd.f"
		dlaset_("L", &i__2, &i__1, &c_b227, &c_b227, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 782 "dgesdd.f"
		ie = itau;
#line 783 "dgesdd.f"
		itauq = ie + *n;
#line 784 "dgesdd.f"
		itaup = itauq + *n;
#line 785 "dgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 790 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 790 "dgesdd.f"
		dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+N*N+BDSPAC) */

#line 799 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB) */

#line 807 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 807 "dgesdd.f"
		dormbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 810 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 810 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              (Workspace: need N*N) */

#line 818 "dgesdd.f"
		dgemm_("N", "N", m, n, n, &c_b248, &u[u_offset], ldu, &work[
			iu], &ldwrku, &c_b227, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 823 "dgesdd.f"
		dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 825 "dgesdd.f"
	    }

#line 827 "dgesdd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 5 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 834 "dgesdd.f"
	    ie = 1;
#line 835 "dgesdd.f"
	    itauq = ie + *n;
#line 836 "dgesdd.f"
	    itaup = itauq + *n;
#line 837 "dgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 842 "dgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 842 "dgesdd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 845 "dgesdd.f"
	    if (wntqn) {

/*              Perform bidiagonal SVD, only computing singular values */
/*              (Workspace: need N+BDSPAC) */

#line 850 "dgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 852 "dgesdd.f"
	    } else if (wntqo) {
#line 853 "dgesdd.f"
		iu = nwork;
#line 854 "dgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 WORK( IU ) is M by N */

#line 858 "dgesdd.f"
		    ldwrku = *m;
#line 859 "dgesdd.f"
		    nwork = iu + ldwrku * *n;
#line 860 "dgesdd.f"
		    dlaset_("F", m, n, &c_b227, &c_b227, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 862 "dgesdd.f"
		} else {

/*                 WORK( IU ) is N by N */

#line 866 "dgesdd.f"
		    ldwrku = *n;
#line 867 "dgesdd.f"
		    nwork = iu + ldwrku * *n;

/*                 WORK(IR) is LDWRKR by N */

#line 871 "dgesdd.f"
		    ir = nwork;
#line 872 "dgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 873 "dgesdd.f"
		}
#line 874 "dgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+N*N+BDSPAC) */

#line 881 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], &ldwrku, &
			vt[vt_offset], ldvt, dum, idum, &work[nwork], &iwork[
			1], info, (ftnlen)1, (ftnlen)1);

/*              Overwrite VT by right singular vectors of A */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 888 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 888 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 892 "dgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 Overwrite WORK(IU) by left singular vectors of A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 897 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 897 "dgesdd.f"
		    dormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy left singular vectors of A from WORK(IU) to A */

#line 903 "dgesdd.f"
		    dlacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 904 "dgesdd.f"
		} else {

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 909 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 909 "dgesdd.f"
		    dorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by left singular vectors of */
/*                 bidiagonal matrix in WORK(IU), storing result in */
/*                 WORK(IR) and copying to A */
/*                 (Workspace: need 2*N*N, prefer N*N+M*N) */

#line 917 "dgesdd.f"
		    i__2 = *m;
#line 917 "dgesdd.f"
		    i__1 = ldwrkr;
#line 917 "dgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 918 "dgesdd.f"
			i__3 = *m - i__ + 1;
#line 918 "dgesdd.f"
			chunk = min(i__3,ldwrkr);
#line 919 "dgesdd.f"
			dgemm_("N", "N", &chunk, n, n, &c_b248, &a[i__ + 
				a_dim1], lda, &work[iu], &ldwrku, &c_b227, &
				work[ir], &ldwrkr, (ftnlen)1, (ftnlen)1);
#line 922 "dgesdd.f"
			dlacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 924 "dgesdd.f"
/* L20: */
#line 924 "dgesdd.f"
		    }
#line 925 "dgesdd.f"
		}

#line 927 "dgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+BDSPAC) */

#line 934 "dgesdd.f"
		dlaset_("F", m, n, &c_b227, &c_b227, &u[u_offset], ldu, (
			ftnlen)1);
#line 935 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need 3*N, prefer 2*N+N*NB) */

#line 943 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 943 "dgesdd.f"
		dormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 946 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 946 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 949 "dgesdd.f"
	    } else if (wntqa) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+BDSPAC) */

#line 956 "dgesdd.f"
		dlaset_("F", m, m, &c_b227, &c_b227, &u[u_offset], ldu, (
			ftnlen)1);
#line 957 "dgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 963 "dgesdd.f"
		if (*m > *n) {
#line 964 "dgesdd.f"
		    i__1 = *m - *n;
#line 964 "dgesdd.f"
		    i__2 = *m - *n;
#line 964 "dgesdd.f"
		    dlaset_("F", &i__1, &i__2, &c_b227, &c_b248, &u[*n + 1 + (
			    *n + 1) * u_dim1], ldu, (ftnlen)1);
#line 966 "dgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need N*N+2*N+M, prefer N*N+2*N+M*NB) */

#line 972 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 972 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 975 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 975 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 978 "dgesdd.f"
	    }

#line 980 "dgesdd.f"
	}

#line 982 "dgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 988 "dgesdd.f"
	if (*n >= mnthr) {

#line 990 "dgesdd.f"
	    if (wntqn) {

/*              Path 1t (N much larger than M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 995 "dgesdd.f"
		itau = 1;
#line 996 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need 2*M, prefer M+M*NB) */

#line 1001 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1001 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out above L */

#line 1006 "dgesdd.f"
		i__1 = *m - 1;
#line 1006 "dgesdd.f"
		i__2 = *m - 1;
#line 1006 "dgesdd.f"
		dlaset_("U", &i__1, &i__2, &c_b227, &c_b227, &a[(a_dim1 << 1) 
			+ 1], lda, (ftnlen)1);
#line 1007 "dgesdd.f"
		ie = 1;
#line 1008 "dgesdd.f"
		itauq = ie + *m;
#line 1009 "dgesdd.f"
		itaup = itauq + *m;
#line 1010 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 1015 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1015 "dgesdd.f"
		dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 1018 "dgesdd.f"
		nwork = ie + *m;

/*              Perform bidiagonal SVD, computing singular values only */
/*              (Workspace: need M+BDSPAC) */

#line 1023 "dgesdd.f"
		dbdsdc_("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 1026 "dgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N much larger than M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1032 "dgesdd.f"
		ivt = 1;

/*              IVT is M by M */

#line 1036 "dgesdd.f"
		il = ivt + *m * *m;
#line 1037 "dgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3 + bdspac) {

/*                 WORK(IL) is M by N */

#line 1041 "dgesdd.f"
		    ldwrkl = *m;
#line 1042 "dgesdd.f"
		    chunk = *n;
#line 1043 "dgesdd.f"
		} else {
#line 1044 "dgesdd.f"
		    ldwrkl = *m;
#line 1045 "dgesdd.f"
		    chunk = (*lwork - *m * *m) / *m;
#line 1046 "dgesdd.f"
		}
#line 1047 "dgesdd.f"
		itau = il + ldwrkl * *m;
#line 1048 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1053 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1053 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1058 "dgesdd.f"
		dlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1059 "dgesdd.f"
		i__1 = *m - 1;
#line 1059 "dgesdd.f"
		i__2 = *m - 1;
#line 1059 "dgesdd.f"
		dlaset_("U", &i__1, &i__2, &c_b227, &c_b227, &work[il + 
			ldwrkl], &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1065 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1065 "dgesdd.f"
		dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1067 "dgesdd.f"
		ie = itau;
#line 1068 "dgesdd.f"
		itauq = ie + *m;
#line 1069 "dgesdd.f"
		itaup = itauq + *m;
#line 1070 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 1075 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1075 "dgesdd.f"
		dgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U, and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              (Workspace: need M+M*M+BDSPAC) */

#line 1084 "dgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], m, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              (Workspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */

#line 1092 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1092 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1095 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1095 "dgesdd.f"
		dormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], m, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              (Workspace: need 2*M*M, prefer M*M+M*N) */

#line 1103 "dgesdd.f"
		i__1 = *n;
#line 1103 "dgesdd.f"
		i__2 = chunk;
#line 1103 "dgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1104 "dgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1104 "dgesdd.f"
		    blk = min(i__3,chunk);
#line 1105 "dgesdd.f"
		    dgemm_("N", "N", m, &blk, m, &c_b248, &work[ivt], m, &a[
			    i__ * a_dim1 + 1], lda, &c_b227, &work[il], &
			    ldwrkl, (ftnlen)1, (ftnlen)1);
#line 1107 "dgesdd.f"
		    dlacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1109 "dgesdd.f"
/* L30: */
#line 1109 "dgesdd.f"
		}

#line 1111 "dgesdd.f"
	    } else if (wntqs) {

/*              Path 3t (N much larger than M, JOBZ='S') */
/*              M right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1117 "dgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1121 "dgesdd.f"
		ldwrkl = *m;
#line 1122 "dgesdd.f"
		itau = il + ldwrkl * *m;
#line 1123 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1128 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1128 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1133 "dgesdd.f"
		dlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1134 "dgesdd.f"
		i__2 = *m - 1;
#line 1134 "dgesdd.f"
		i__1 = *m - 1;
#line 1134 "dgesdd.f"
		dlaset_("U", &i__2, &i__1, &c_b227, &c_b227, &work[il + 
			ldwrkl], &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1140 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1140 "dgesdd.f"
		dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1142 "dgesdd.f"
		ie = itau;
#line 1143 "dgesdd.f"
		itauq = ie + *m;
#line 1144 "dgesdd.f"
		itaup = itauq + *m;
#line 1145 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IU), copying result to U */
/*              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 1150 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1150 "dgesdd.f"
		dgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need M+BDSPAC) */

#line 1159 "dgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and VT */
/*              by right singular vectors of L */
/*              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB) */

#line 1167 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1167 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1170 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1170 "dgesdd.f"
		dormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by */
/*              Q in A, storing result in VT */
/*              (Workspace: need M*M) */

#line 1178 "dgesdd.f"
		dlacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1179 "dgesdd.f"
		dgemm_("N", "N", m, n, m, &c_b248, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b227, &vt[vt_offset], ldvt, (
			ftnlen)1, (ftnlen)1);

#line 1182 "dgesdd.f"
	    } else if (wntqa) {

/*              Path 4t (N much larger than M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1188 "dgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1192 "dgesdd.f"
		ldwkvt = *m;
#line 1193 "dgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1194 "dgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1199 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1199 "dgesdd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 1201 "dgesdd.f"
		dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1206 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1206 "dgesdd.f"
		dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__2, &ierr);

/*              Produce L in A, zeroing out other entries */

#line 1211 "dgesdd.f"
		i__2 = *m - 1;
#line 1211 "dgesdd.f"
		i__1 = *m - 1;
#line 1211 "dgesdd.f"
		dlaset_("U", &i__2, &i__1, &c_b227, &c_b227, &a[(a_dim1 << 1) 
			+ 1], lda, (ftnlen)1);
#line 1212 "dgesdd.f"
		ie = itau;
#line 1213 "dgesdd.f"
		itauq = ie + *m;
#line 1214 "dgesdd.f"
		itaup = itauq + *m;
#line 1215 "dgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 1220 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1220 "dgesdd.f"
		dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              (Workspace: need M+M*M+BDSPAC) */

#line 1229 "dgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB) */

#line 1237 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1237 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1240 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1240 "dgesdd.f"
		dormbr_("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              (Workspace: need M*M) */

#line 1248 "dgesdd.f"
		dgemm_("N", "N", m, n, m, &c_b248, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b227, &a[a_offset], lda, (ftnlen)
			1, (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1253 "dgesdd.f"
		dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1255 "dgesdd.f"
	    }

#line 1257 "dgesdd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 5t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 1264 "dgesdd.f"
	    ie = 1;
#line 1265 "dgesdd.f"
	    itauq = ie + *m;
#line 1266 "dgesdd.f"
	    itaup = itauq + *m;
#line 1267 "dgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 1272 "dgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1272 "dgesdd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 1275 "dgesdd.f"
	    if (wntqn) {

/*              Perform bidiagonal SVD, only computing singular values */
/*              (Workspace: need M+BDSPAC) */

#line 1280 "dgesdd.f"
		dbdsdc_("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 1282 "dgesdd.f"
	    } else if (wntqo) {
#line 1283 "dgesdd.f"
		ldwkvt = *m;
#line 1284 "dgesdd.f"
		ivt = nwork;
#line 1285 "dgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 WORK( IVT ) is M by N */

#line 1289 "dgesdd.f"
		    dlaset_("F", m, n, &c_b227, &c_b227, &work[ivt], &ldwkvt, 
			    (ftnlen)1);
#line 1291 "dgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1292 "dgesdd.f"
		} else {

/*                 WORK( IVT ) is M by M */

#line 1296 "dgesdd.f"
		    nwork = ivt + ldwkvt * *m;
#line 1297 "dgesdd.f"
		    il = nwork;

/*                 WORK(IL) is M by CHUNK */

#line 1301 "dgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1302 "dgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              (Workspace: need M*M+BDSPAC) */

#line 1309 "dgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1316 "dgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1316 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1320 "dgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 Overwrite WORK(IVT) by left singular vectors of A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1325 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1325 "dgesdd.f"
		    dormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy right singular vectors of A from WORK(IVT) to A */

#line 1331 "dgesdd.f"
		    dlacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 1332 "dgesdd.f"
		} else {

/*                 Generate P**T in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1337 "dgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1337 "dgesdd.f"
		    dorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by right singular vectors of */
/*                 bidiagonal matrix in WORK(IVT), storing result in */
/*                 WORK(IL) and copying to A */
/*                 (Workspace: need 2*M*M, prefer M*M+M*N) */

#line 1345 "dgesdd.f"
		    i__2 = *n;
#line 1345 "dgesdd.f"
		    i__1 = chunk;
#line 1345 "dgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1346 "dgesdd.f"
			i__3 = *n - i__ + 1;
#line 1346 "dgesdd.f"
			blk = min(i__3,chunk);
#line 1347 "dgesdd.f"
			dgemm_("N", "N", m, &blk, m, &c_b248, &work[ivt], &
				ldwkvt, &a[i__ * a_dim1 + 1], lda, &c_b227, &
				work[il], m, (ftnlen)1, (ftnlen)1);
#line 1350 "dgesdd.f"
			dlacpy_("F", m, &blk, &work[il], m, &a[i__ * a_dim1 + 
				1], lda, (ftnlen)1);
#line 1352 "dgesdd.f"
/* L40: */
#line 1352 "dgesdd.f"
		    }
#line 1353 "dgesdd.f"
		}
#line 1354 "dgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need M+BDSPAC) */

#line 1361 "dgesdd.f"
		dlaset_("F", m, n, &c_b227, &c_b227, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1362 "dgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need 3*M, prefer 2*M+M*NB) */

#line 1370 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1370 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1373 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1373 "dgesdd.f"
		dormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1376 "dgesdd.f"
	    } else if (wntqa) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need M+BDSPAC) */

#line 1383 "dgesdd.f"
		dlaset_("F", n, n, &c_b227, &c_b227, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1384 "dgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of VT to identity matrix */

#line 1390 "dgesdd.f"
		if (*n > *m) {
#line 1391 "dgesdd.f"
		    i__1 = *n - *m;
#line 1391 "dgesdd.f"
		    i__2 = *n - *m;
#line 1391 "dgesdd.f"
		    dlaset_("F", &i__1, &i__2, &c_b227, &c_b248, &vt[*m + 1 + 
			    (*m + 1) * vt_dim1], ldvt, (ftnlen)1);
#line 1393 "dgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need 2*M+N, prefer 2*M+N*NB) */

#line 1399 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1399 "dgesdd.f"
		dormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1402 "dgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1402 "dgesdd.f"
		dormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1405 "dgesdd.f"
	    }

#line 1407 "dgesdd.f"
	}

#line 1409 "dgesdd.f"
    }

/*     Undo scaling if necessary */

#line 1413 "dgesdd.f"
    if (iscl == 1) {
#line 1414 "dgesdd.f"
	if (anrm > bignum) {
#line 1414 "dgesdd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1414 "dgesdd.f"
	}
#line 1417 "dgesdd.f"
	if (anrm < smlnum) {
#line 1417 "dgesdd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1417 "dgesdd.f"
	}
#line 1420 "dgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 1424 "dgesdd.f"
    work[1] = (doublereal) maxwrk;

#line 1426 "dgesdd.f"
    return 0;

/*     End of DGESDD */

} /* dgesdd_ */

