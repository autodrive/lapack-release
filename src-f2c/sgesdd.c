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

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b227 = 0.;
static doublereal c_b248 = 1.;

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

/*       SUBROUTINE SGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, */
/*                          LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), S( * ), U( LDU, * ), */
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
/* >          The leading dimension of the array VT.  LDVT >= 1; if */
/* >          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N; */
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
/* >          If JOBZ = 'N', */
/* >            LWORK >= 3*min(M,N) + max(max(M,N),6*min(M,N)). */
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
/* >          > 0:  SBDSDC did not converge, updating process failed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

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
    static integer i__, ie, il, ir, iu, blk;
    static doublereal dum[1], eps;
    static integer ivt, iscl;
    static doublereal anrm;
    static integer idum[1], ierr, itau;
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
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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

#line 270 "sgesdd.f"
    /* Parameter adjustments */
#line 270 "sgesdd.f"
    a_dim1 = *lda;
#line 270 "sgesdd.f"
    a_offset = 1 + a_dim1;
#line 270 "sgesdd.f"
    a -= a_offset;
#line 270 "sgesdd.f"
    --s;
#line 270 "sgesdd.f"
    u_dim1 = *ldu;
#line 270 "sgesdd.f"
    u_offset = 1 + u_dim1;
#line 270 "sgesdd.f"
    u -= u_offset;
#line 270 "sgesdd.f"
    vt_dim1 = *ldvt;
#line 270 "sgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 270 "sgesdd.f"
    vt -= vt_offset;
#line 270 "sgesdd.f"
    --work;
#line 270 "sgesdd.f"
    --iwork;
#line 270 "sgesdd.f"

#line 270 "sgesdd.f"
    /* Function Body */
#line 270 "sgesdd.f"
    *info = 0;
#line 271 "sgesdd.f"
    minmn = min(*m,*n);
#line 272 "sgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 273 "sgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 274 "sgesdd.f"
    wntqas = wntqa || wntqs;
#line 275 "sgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 276 "sgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 277 "sgesdd.f"
    lquery = *lwork == -1;

#line 279 "sgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 280 "sgesdd.f"
	*info = -1;
#line 281 "sgesdd.f"
    } else if (*m < 0) {
#line 282 "sgesdd.f"
	*info = -2;
#line 283 "sgesdd.f"
    } else if (*n < 0) {
#line 284 "sgesdd.f"
	*info = -3;
#line 285 "sgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 286 "sgesdd.f"
	*info = -5;
#line 287 "sgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 289 "sgesdd.f"
	*info = -8;
#line 290 "sgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 293 "sgesdd.f"
	*info = -10;
#line 294 "sgesdd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 303 "sgesdd.f"
    if (*info == 0) {
#line 304 "sgesdd.f"
	minwrk = 1;
#line 305 "sgesdd.f"
	maxwrk = 1;
#line 306 "sgesdd.f"
	if (*m >= *n && minmn > 0) {

/*           Compute space needed for SBDSDC */

#line 310 "sgesdd.f"
	    mnthr = (integer) (minmn * 11. / 6.);
#line 311 "sgesdd.f"
	    if (wntqn) {
#line 312 "sgesdd.f"
		bdspac = *n * 7;
#line 313 "sgesdd.f"
	    } else {
#line 314 "sgesdd.f"
		bdspac = *n * 3 * *n + (*n << 2);
#line 315 "sgesdd.f"
	    }
#line 316 "sgesdd.f"
	    if (*m >= mnthr) {
#line 317 "sgesdd.f"
		if (wntqn) {

/*                 Path 1 (M much larger than N, JOBZ='N') */

#line 321 "sgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 323 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 323 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 325 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n;
#line 325 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 326 "sgesdd.f"
		    minwrk = bdspac + *n;
#line 327 "sgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M much larger than N, JOBZ='O') */

#line 331 "sgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 332 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "SORGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 332 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 334 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 334 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 336 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "QLN", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 336 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 338 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 338 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 340 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 340 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 341 "sgesdd.f"
		    maxwrk = wrkbl + (*n << 1) * *n;
#line 342 "sgesdd.f"
		    minwrk = bdspac + (*n << 1) * *n + *n * 3;
#line 343 "sgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M much larger than N, JOBZ='S') */

#line 347 "sgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 348 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "SORGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 348 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 350 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 350 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 352 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "QLN", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 352 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 354 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 354 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 356 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 356 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 357 "sgesdd.f"
		    maxwrk = wrkbl + *n * *n;
#line 358 "sgesdd.f"
		    minwrk = bdspac + *n * *n + *n * 3;
#line 359 "sgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M much larger than N, JOBZ='A') */

#line 363 "sgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 364 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *m * ilaenv_(&c__1, "SORGQR", 
			    " ", m, m, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 364 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 366 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + (*n << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 366 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 368 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "QLN", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 368 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 370 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 370 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 372 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 372 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 373 "sgesdd.f"
		    maxwrk = wrkbl + *n * *n;
#line 374 "sgesdd.f"
		    minwrk = bdspac + *n * *n + (*n << 1) + *m;
#line 375 "sgesdd.f"
		}
#line 376 "sgesdd.f"
	    } else {

/*              Path 5 (M at least N, but not much larger) */

#line 380 "sgesdd.f"
		wrkbl = *n * 3 + (*m + *n) * ilaenv_(&c__1, "SGEBRD", " ", m, 
			n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 382 "sgesdd.f"
		if (wntqn) {
/* Computing MAX */
#line 383 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 383 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 384 "sgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 385 "sgesdd.f"
		} else if (wntqo) {
/* Computing MAX */
#line 386 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 386 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 388 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 388 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 390 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 390 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 391 "sgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 392 "sgesdd.f"
		    i__1 = *m, i__2 = *n * *n + bdspac;
#line 392 "sgesdd.f"
		    minwrk = *n * 3 + max(i__1,i__2);
#line 393 "sgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 394 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 394 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 396 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 396 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 398 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
#line 398 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 399 "sgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 400 "sgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 401 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 401 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 403 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *n * 3 + *n * ilaenv_(&c__1, "SORMBR"
			    , "PRT", n, n, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 403 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 405 "sgesdd.f"
		    i__1 = maxwrk, i__2 = bdspac + *n * 3;
#line 405 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 406 "sgesdd.f"
		    minwrk = *n * 3 + max(*m,bdspac);
#line 407 "sgesdd.f"
		}
#line 408 "sgesdd.f"
	    }
#line 409 "sgesdd.f"
	} else if (minmn > 0) {

/*           Compute space needed for SBDSDC */

#line 413 "sgesdd.f"
	    mnthr = (integer) (minmn * 11. / 6.);
#line 414 "sgesdd.f"
	    if (wntqn) {
#line 415 "sgesdd.f"
		bdspac = *m * 7;
#line 416 "sgesdd.f"
	    } else {
#line 417 "sgesdd.f"
		bdspac = *m * 3 * *m + (*m << 2);
#line 418 "sgesdd.f"
	    }
#line 419 "sgesdd.f"
	    if (*n >= mnthr) {
#line 420 "sgesdd.f"
		if (wntqn) {

/*                 Path 1t (N much larger than M, JOBZ='N') */

#line 424 "sgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 426 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 426 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 428 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m;
#line 428 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 429 "sgesdd.f"
		    minwrk = bdspac + *m;
#line 430 "sgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N much larger than M, JOBZ='O') */

#line 434 "sgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 435 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "SORGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 435 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 437 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 437 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 439 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 439 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 441 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "PRT", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 441 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 443 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 443 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 444 "sgesdd.f"
		    maxwrk = wrkbl + (*m << 1) * *m;
#line 445 "sgesdd.f"
		    minwrk = bdspac + (*m << 1) * *m + *m * 3;
#line 446 "sgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N much larger than M, JOBZ='S') */

#line 450 "sgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 451 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "SORGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 451 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 453 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 453 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 455 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 455 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 457 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "PRT", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 457 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 459 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 459 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 460 "sgesdd.f"
		    maxwrk = wrkbl + *m * *m;
#line 461 "sgesdd.f"
		    minwrk = bdspac + *m * *m + *m * 3;
#line 462 "sgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N much larger than M, JOBZ='A') */

#line 466 "sgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 467 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *n * ilaenv_(&c__1, "SORGLQ", 
			    " ", n, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 467 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 469 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + (*m << 1) * ilaenv_(&c__1, 
			    "SGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
#line 469 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 471 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 471 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 473 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "PRT", m, m, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 473 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 475 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 475 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 476 "sgesdd.f"
		    maxwrk = wrkbl + *m * *m;
#line 477 "sgesdd.f"
		    minwrk = bdspac + *m * *m + *m * 3;
#line 478 "sgesdd.f"
		}
#line 479 "sgesdd.f"
	    } else {

/*              Path 5t (N greater than M, but not much larger) */

#line 483 "sgesdd.f"
		wrkbl = *m * 3 + (*m + *n) * ilaenv_(&c__1, "SGEBRD", " ", m, 
			n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 485 "sgesdd.f"
		if (wntqn) {
/* Computing MAX */
#line 486 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 486 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 487 "sgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 488 "sgesdd.f"
		} else if (wntqo) {
/* Computing MAX */
#line 489 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 489 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 491 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "PRT", m, n, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 491 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 493 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 493 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 494 "sgesdd.f"
		    maxwrk = wrkbl + *m * *n;
/* Computing MAX */
#line 495 "sgesdd.f"
		    i__1 = *n, i__2 = *m * *m + bdspac;
#line 495 "sgesdd.f"
		    minwrk = *m * 3 + max(i__1,i__2);
#line 496 "sgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 497 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 497 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 499 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "PRT", m, n, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 499 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 501 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 501 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 502 "sgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 503 "sgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 504 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "QLN", m, m, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 504 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 506 "sgesdd.f"
		    i__1 = wrkbl, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORMBR"
			    , "PRT", n, n, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 506 "sgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 508 "sgesdd.f"
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
#line 508 "sgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 509 "sgesdd.f"
		    minwrk = *m * 3 + max(*n,bdspac);
#line 510 "sgesdd.f"
		}
#line 511 "sgesdd.f"
	    }
#line 512 "sgesdd.f"
	}
#line 513 "sgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 514 "sgesdd.f"
	work[1] = (doublereal) maxwrk;

#line 516 "sgesdd.f"
	if (*lwork < minwrk && ! lquery) {
#line 517 "sgesdd.f"
	    *info = -12;
#line 518 "sgesdd.f"
	}
#line 519 "sgesdd.f"
    }

#line 521 "sgesdd.f"
    if (*info != 0) {
#line 522 "sgesdd.f"
	i__1 = -(*info);
#line 522 "sgesdd.f"
	xerbla_("SGESDD", &i__1, (ftnlen)6);
#line 523 "sgesdd.f"
	return 0;
#line 524 "sgesdd.f"
    } else if (lquery) {
#line 525 "sgesdd.f"
	return 0;
#line 526 "sgesdd.f"
    }

/*     Quick return if possible */

#line 530 "sgesdd.f"
    if (*m == 0 || *n == 0) {
#line 531 "sgesdd.f"
	return 0;
#line 532 "sgesdd.f"
    }

/*     Get machine constants */

#line 536 "sgesdd.f"
    eps = slamch_("P", (ftnlen)1);
#line 537 "sgesdd.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 538 "sgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 542 "sgesdd.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 543 "sgesdd.f"
    iscl = 0;
#line 544 "sgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 545 "sgesdd.f"
	iscl = 1;
#line 546 "sgesdd.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 547 "sgesdd.f"
    } else if (anrm > bignum) {
#line 548 "sgesdd.f"
	iscl = 1;
#line 549 "sgesdd.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 550 "sgesdd.f"
    }

#line 552 "sgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 558 "sgesdd.f"
	if (*m >= mnthr) {

#line 560 "sgesdd.f"
	    if (wntqn) {

/*              Path 1 (M much larger than N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 565 "sgesdd.f"
		itau = 1;
#line 566 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need 2*N, prefer N+N*NB) */

#line 571 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 571 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 576 "sgesdd.f"
		i__1 = *n - 1;
#line 576 "sgesdd.f"
		i__2 = *n - 1;
#line 576 "sgesdd.f"
		slaset_("L", &i__1, &i__2, &c_b227, &c_b227, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 577 "sgesdd.f"
		ie = 1;
#line 578 "sgesdd.f"
		itauq = ie + *n;
#line 579 "sgesdd.f"
		itaup = itauq + *n;
#line 580 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 585 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 585 "sgesdd.f"
		sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 588 "sgesdd.f"
		nwork = ie + *n;

/*              Perform bidiagonal SVD, computing singular values only */
/*              (Workspace: need N+BDSPAC) */

#line 593 "sgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 596 "sgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M much larger than N, JOBZ = 'O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 602 "sgesdd.f"
		ir = 1;

/*              WORK(IR) is LDWRKR by N */

#line 606 "sgesdd.f"
		if (*lwork >= *lda * *n + *n * *n + *n * 3 + bdspac) {
#line 607 "sgesdd.f"
		    ldwrkr = *lda;
#line 608 "sgesdd.f"
		} else {
#line 609 "sgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3 - bdspac) / *n;
#line 610 "sgesdd.f"
		}
#line 611 "sgesdd.f"
		itau = ir + ldwrkr * *n;
#line 612 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 617 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 617 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 622 "sgesdd.f"
		slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 623 "sgesdd.f"
		i__1 = *n - 1;
#line 623 "sgesdd.f"
		i__2 = *n - 1;
#line 623 "sgesdd.f"
		slaset_("L", &i__1, &i__2, &c_b227, &c_b227, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 629 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 629 "sgesdd.f"
		sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 631 "sgesdd.f"
		ie = itau;
#line 632 "sgesdd.f"
		itauq = ie + *n;
#line 633 "sgesdd.f"
		itaup = itauq + *n;
#line 634 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in VT, copying result to WORK(IR) */
/*              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 639 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 639 "sgesdd.f"
		sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              WORK(IU) is N by N */

#line 645 "sgesdd.f"
		iu = nwork;
#line 646 "sgesdd.f"
		nwork = iu + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+N*N+BDSPAC) */

#line 653 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R */
/*              and VT by right singular vectors of R */
/*              (Workspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */

#line 661 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 661 "sgesdd.f"
		sormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], n, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 664 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 664 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              (Workspace: need 2*N*N, prefer N*N+M*N) */

#line 672 "sgesdd.f"
		i__1 = *m;
#line 672 "sgesdd.f"
		i__2 = ldwrkr;
#line 672 "sgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 673 "sgesdd.f"
		    i__3 = *m - i__ + 1;
#line 673 "sgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 674 "sgesdd.f"
		    sgemm_("N", "N", &chunk, n, n, &c_b248, &a[i__ + a_dim1], 
			    lda, &work[iu], n, &c_b227, &work[ir], &ldwrkr, (
			    ftnlen)1, (ftnlen)1);
#line 677 "sgesdd.f"
		    slacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 679 "sgesdd.f"
/* L10: */
#line 679 "sgesdd.f"
		}

#line 681 "sgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M much larger than N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 687 "sgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 691 "sgesdd.f"
		ldwrkr = *n;
#line 692 "sgesdd.f"
		itau = ir + ldwrkr * *n;
#line 693 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 698 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 698 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 703 "sgesdd.f"
		slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 704 "sgesdd.f"
		i__2 = *n - 1;
#line 704 "sgesdd.f"
		i__1 = *n - 1;
#line 704 "sgesdd.f"
		slaset_("L", &i__2, &i__1, &c_b227, &c_b227, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 710 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 710 "sgesdd.f"
		sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 712 "sgesdd.f"
		ie = itau;
#line 713 "sgesdd.f"
		itauq = ie + *n;
#line 714 "sgesdd.f"
		itaup = itauq + *n;
#line 715 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 720 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 720 "sgesdd.f"
		sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagoal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+BDSPAC) */

#line 729 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB) */

#line 737 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 737 "sgesdd.f"
		sormbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 741 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 741 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              (Workspace: need N*N) */

#line 749 "sgesdd.f"
		slacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 750 "sgesdd.f"
		sgemm_("N", "N", m, n, n, &c_b248, &a[a_offset], lda, &work[
			ir], &ldwrkr, &c_b227, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 753 "sgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M much larger than N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 759 "sgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 763 "sgesdd.f"
		ldwrku = *n;
#line 764 "sgesdd.f"
		itau = iu + ldwrku * *n;
#line 765 "sgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 770 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 770 "sgesdd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 772 "sgesdd.f"
		slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */
#line 776 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 776 "sgesdd.f"
		sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out other entries */

#line 781 "sgesdd.f"
		i__2 = *n - 1;
#line 781 "sgesdd.f"
		i__1 = *n - 1;
#line 781 "sgesdd.f"
		slaset_("L", &i__2, &i__1, &c_b227, &c_b227, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 782 "sgesdd.f"
		ie = itau;
#line 783 "sgesdd.f"
		itauq = ie + *n;
#line 784 "sgesdd.f"
		itaup = itauq + *n;
#line 785 "sgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 790 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 790 "sgesdd.f"
		sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+N*N+BDSPAC) */

#line 799 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite WORK(IU) by left singular vectors of R and VT */
/*              by right singular vectors of R */
/*              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB) */

#line 807 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 807 "sgesdd.f"
		sormbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 810 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 810 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              (Workspace: need N*N) */

#line 818 "sgesdd.f"
		sgemm_("N", "N", m, n, n, &c_b248, &u[u_offset], ldu, &work[
			iu], &ldwrku, &c_b227, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 823 "sgesdd.f"
		slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 825 "sgesdd.f"
	    }

#line 827 "sgesdd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 5 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 834 "sgesdd.f"
	    ie = 1;
#line 835 "sgesdd.f"
	    itauq = ie + *n;
#line 836 "sgesdd.f"
	    itaup = itauq + *n;
#line 837 "sgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 842 "sgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 842 "sgesdd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 845 "sgesdd.f"
	    if (wntqn) {

/*              Perform bidiagonal SVD, only computing singular values */
/*              (Workspace: need N+BDSPAC) */

#line 850 "sgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 852 "sgesdd.f"
	    } else if (wntqo) {
#line 853 "sgesdd.f"
		iu = nwork;
#line 854 "sgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 WORK( IU ) is M by N */

#line 858 "sgesdd.f"
		    ldwrku = *m;
#line 859 "sgesdd.f"
		    nwork = iu + ldwrku * *n;
#line 860 "sgesdd.f"
		    slaset_("F", m, n, &c_b227, &c_b227, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 862 "sgesdd.f"
		} else {

/*                 WORK( IU ) is N by N */

#line 866 "sgesdd.f"
		    ldwrku = *n;
#line 867 "sgesdd.f"
		    nwork = iu + ldwrku * *n;

/*                 WORK(IR) is LDWRKR by N */

#line 871 "sgesdd.f"
		    ir = nwork;
#line 872 "sgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 873 "sgesdd.f"
		}
#line 874 "sgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in WORK(IU) and computing right */
/*              singular vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+N*N+BDSPAC) */

#line 881 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &work[iu], &ldwrku, &
			vt[vt_offset], ldvt, dum, idum, &work[nwork], &iwork[
			1], info, (ftnlen)1, (ftnlen)1);

/*              Overwrite VT by right singular vectors of A */
/*              (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 888 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 888 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 892 "sgesdd.f"
		if (*lwork >= *m * *n + *n * 3 + bdspac) {

/*                 Overwrite WORK(IU) by left singular vectors of A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 897 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 897 "sgesdd.f"
		    sormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy left singular vectors of A from WORK(IU) to A */

#line 903 "sgesdd.f"
		    slacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 904 "sgesdd.f"
		} else {

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 909 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 909 "sgesdd.f"
		    sorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by left singular vectors of */
/*                 bidiagonal matrix in WORK(IU), storing result in */
/*                 WORK(IR) and copying to A */
/*                 (Workspace: need 2*N*N, prefer N*N+M*N) */

#line 917 "sgesdd.f"
		    i__2 = *m;
#line 917 "sgesdd.f"
		    i__1 = ldwrkr;
#line 917 "sgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 918 "sgesdd.f"
			i__3 = *m - i__ + 1;
#line 918 "sgesdd.f"
			chunk = min(i__3,ldwrkr);
#line 919 "sgesdd.f"
			sgemm_("N", "N", &chunk, n, n, &c_b248, &a[i__ + 
				a_dim1], lda, &work[iu], &ldwrku, &c_b227, &
				work[ir], &ldwrkr, (ftnlen)1, (ftnlen)1);
#line 922 "sgesdd.f"
			slacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 924 "sgesdd.f"
/* L20: */
#line 924 "sgesdd.f"
		    }
#line 925 "sgesdd.f"
		}

#line 927 "sgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+BDSPAC) */

#line 934 "sgesdd.f"
		slaset_("F", m, n, &c_b227, &c_b227, &u[u_offset], ldu, (
			ftnlen)1);
#line 935 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need 3*N, prefer 2*N+N*NB) */

#line 943 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 943 "sgesdd.f"
		sormbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 946 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 946 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 949 "sgesdd.f"
	    } else if (wntqa) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need N+BDSPAC) */

#line 956 "sgesdd.f"
		slaset_("F", m, m, &c_b227, &c_b227, &u[u_offset], ldu, (
			ftnlen)1);
#line 957 "sgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 963 "sgesdd.f"
		if (*m > *n) {
#line 964 "sgesdd.f"
		    i__1 = *m - *n;
#line 964 "sgesdd.f"
		    i__2 = *m - *n;
#line 964 "sgesdd.f"
		    slaset_("F", &i__1, &i__2, &c_b227, &c_b248, &u[*n + 1 + (
			    *n + 1) * u_dim1], ldu, (ftnlen)1);
#line 966 "sgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need N*N+2*N+M, prefer N*N+2*N+M*NB) */

#line 972 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 972 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 975 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 975 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 978 "sgesdd.f"
	    }

#line 980 "sgesdd.f"
	}

#line 982 "sgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 988 "sgesdd.f"
	if (*n >= mnthr) {

#line 990 "sgesdd.f"
	    if (wntqn) {

/*              Path 1t (N much larger than M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 995 "sgesdd.f"
		itau = 1;
#line 996 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need 2*M, prefer M+M*NB) */

#line 1001 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1001 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out above L */

#line 1006 "sgesdd.f"
		i__1 = *m - 1;
#line 1006 "sgesdd.f"
		i__2 = *m - 1;
#line 1006 "sgesdd.f"
		slaset_("U", &i__1, &i__2, &c_b227, &c_b227, &a[(a_dim1 << 1) 
			+ 1], lda, (ftnlen)1);
#line 1007 "sgesdd.f"
		ie = 1;
#line 1008 "sgesdd.f"
		itauq = ie + *m;
#line 1009 "sgesdd.f"
		itaup = itauq + *m;
#line 1010 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 1015 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1015 "sgesdd.f"
		sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 1018 "sgesdd.f"
		nwork = ie + *m;

/*              Perform bidiagonal SVD, computing singular values only */
/*              (Workspace: need M+BDSPAC) */

#line 1023 "sgesdd.f"
		sbdsdc_("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);

#line 1026 "sgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N much larger than M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1032 "sgesdd.f"
		ivt = 1;

/*              IVT is M by M */

#line 1036 "sgesdd.f"
		il = ivt + *m * *m;
#line 1037 "sgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3 + bdspac) {

/*                 WORK(IL) is M by N */

#line 1041 "sgesdd.f"
		    ldwrkl = *m;
#line 1042 "sgesdd.f"
		    chunk = *n;
#line 1043 "sgesdd.f"
		} else {
#line 1044 "sgesdd.f"
		    ldwrkl = *m;
#line 1045 "sgesdd.f"
		    chunk = (*lwork - *m * *m) / *m;
#line 1046 "sgesdd.f"
		}
#line 1047 "sgesdd.f"
		itau = il + ldwrkl * *m;
#line 1048 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1053 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1053 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1058 "sgesdd.f"
		slacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1059 "sgesdd.f"
		i__1 = *m - 1;
#line 1059 "sgesdd.f"
		i__2 = *m - 1;
#line 1059 "sgesdd.f"
		slaset_("U", &i__1, &i__2, &c_b227, &c_b227, &work[il + 
			ldwrkl], &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1065 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1065 "sgesdd.f"
		sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1067 "sgesdd.f"
		ie = itau;
#line 1068 "sgesdd.f"
		itauq = ie + *m;
#line 1069 "sgesdd.f"
		itaup = itauq + *m;
#line 1070 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 1075 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1075 "sgesdd.f"
		sgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U, and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              (Workspace: need M+M*M+BDSPAC) */

#line 1084 "sgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], m, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              (Workspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */

#line 1092 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1092 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1095 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1095 "sgesdd.f"
		sormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], m, &work[nwork], &i__1, &ierr, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              (Workspace: need 2*M*M, prefer M*M+M*N) */

#line 1103 "sgesdd.f"
		i__1 = *n;
#line 1103 "sgesdd.f"
		i__2 = chunk;
#line 1103 "sgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1104 "sgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1104 "sgesdd.f"
		    blk = min(i__3,chunk);
#line 1105 "sgesdd.f"
		    sgemm_("N", "N", m, &blk, m, &c_b248, &work[ivt], m, &a[
			    i__ * a_dim1 + 1], lda, &c_b227, &work[il], &
			    ldwrkl, (ftnlen)1, (ftnlen)1);
#line 1107 "sgesdd.f"
		    slacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1109 "sgesdd.f"
/* L30: */
#line 1109 "sgesdd.f"
		}

#line 1111 "sgesdd.f"
	    } else if (wntqs) {

/*              Path 3t (N much larger than M, JOBZ='S') */
/*              M right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1117 "sgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1121 "sgesdd.f"
		ldwrkl = *m;
#line 1122 "sgesdd.f"
		itau = il + ldwrkl * *m;
#line 1123 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1128 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1128 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1133 "sgesdd.f"
		slacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1134 "sgesdd.f"
		i__2 = *m - 1;
#line 1134 "sgesdd.f"
		i__1 = *m - 1;
#line 1134 "sgesdd.f"
		slaset_("U", &i__2, &i__1, &c_b227, &c_b227, &work[il + 
			ldwrkl], &ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1140 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1140 "sgesdd.f"
		sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1142 "sgesdd.f"
		ie = itau;
#line 1143 "sgesdd.f"
		itauq = ie + *m;
#line 1144 "sgesdd.f"
		itaup = itauq + *m;
#line 1145 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IU), copying result to U */
/*              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 1150 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1150 "sgesdd.f"
		sgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need M+BDSPAC) */

#line 1159 "sgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and VT */
/*              by right singular vectors of L */
/*              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB) */

#line 1167 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1167 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1170 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1170 "sgesdd.f"
		sormbr_("P", "R", "T", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by */
/*              Q in A, storing result in VT */
/*              (Workspace: need M*M) */

#line 1178 "sgesdd.f"
		slacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1179 "sgesdd.f"
		sgemm_("N", "N", m, n, m, &c_b248, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b227, &vt[vt_offset], ldvt, (
			ftnlen)1, (ftnlen)1);

#line 1182 "sgesdd.f"
	    } else if (wntqa) {

/*              Path 4t (N much larger than M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1188 "sgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1192 "sgesdd.f"
		ldwkvt = *m;
#line 1193 "sgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1194 "sgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1199 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1199 "sgesdd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 1201 "sgesdd.f"
		slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1206 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1206 "sgesdd.f"
		sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__2, &ierr);

/*              Produce L in A, zeroing out other entries */

#line 1211 "sgesdd.f"
		i__2 = *m - 1;
#line 1211 "sgesdd.f"
		i__1 = *m - 1;
#line 1211 "sgesdd.f"
		slaset_("U", &i__2, &i__1, &c_b227, &c_b227, &a[(a_dim1 << 1) 
			+ 1], lda, (ftnlen)1);
#line 1212 "sgesdd.f"
		ie = itau;
#line 1213 "sgesdd.f"
		itauq = ie + *m;
#line 1214 "sgesdd.f"
		itaup = itauq + *m;
#line 1215 "sgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 1220 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1220 "sgesdd.f"
		sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              (Workspace: need M+M*M+BDSPAC) */

#line 1229 "sgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of L and WORK(IVT) */
/*              by right singular vectors of L */
/*              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB) */

#line 1237 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1237 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1240 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1240 "sgesdd.f"
		sormbr_("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              (Workspace: need M*M) */

#line 1248 "sgesdd.f"
		sgemm_("N", "N", m, n, m, &c_b248, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b227, &a[a_offset], lda, (ftnlen)
			1, (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1253 "sgesdd.f"
		slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1255 "sgesdd.f"
	    }

#line 1257 "sgesdd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 5t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 1264 "sgesdd.f"
	    ie = 1;
#line 1265 "sgesdd.f"
	    itauq = ie + *m;
#line 1266 "sgesdd.f"
	    itaup = itauq + *m;
#line 1267 "sgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 1272 "sgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1272 "sgesdd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__2, &ierr);
#line 1275 "sgesdd.f"
	    if (wntqn) {

/*              Perform bidiagonal SVD, only computing singular values */
/*              (Workspace: need M+BDSPAC) */

#line 1280 "sgesdd.f"
		sbdsdc_("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, 
			(ftnlen)1);
#line 1282 "sgesdd.f"
	    } else if (wntqo) {
#line 1283 "sgesdd.f"
		ldwkvt = *m;
#line 1284 "sgesdd.f"
		ivt = nwork;
#line 1285 "sgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 WORK( IVT ) is M by N */

#line 1289 "sgesdd.f"
		    slaset_("F", m, n, &c_b227, &c_b227, &work[ivt], &ldwkvt, 
			    (ftnlen)1);
#line 1291 "sgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1292 "sgesdd.f"
		} else {

/*                 WORK( IVT ) is M by M */

#line 1296 "sgesdd.f"
		    nwork = ivt + ldwkvt * *m;
#line 1297 "sgesdd.f"
		    il = nwork;

/*                 WORK(IL) is M by CHUNK */

#line 1301 "sgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1302 "sgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in WORK(IVT) */
/*              (Workspace: need M*M+BDSPAC) */

#line 1309 "sgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A */
/*              (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1316 "sgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1316 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1320 "sgesdd.f"
		if (*lwork >= *m * *n + *m * 3 + bdspac) {

/*                 Overwrite WORK(IVT) by left singular vectors of A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1325 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1325 "sgesdd.f"
		    sormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Copy right singular vectors of A from WORK(IVT) to A */

#line 1331 "sgesdd.f"
		    slacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 1332 "sgesdd.f"
		} else {

/*                 Generate P**T in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 1337 "sgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1337 "sgesdd.f"
		    sorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by right singular vectors of */
/*                 bidiagonal matrix in WORK(IVT), storing result in */
/*                 WORK(IL) and copying to A */
/*                 (Workspace: need 2*M*M, prefer M*M+M*N) */

#line 1345 "sgesdd.f"
		    i__2 = *n;
#line 1345 "sgesdd.f"
		    i__1 = chunk;
#line 1345 "sgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1346 "sgesdd.f"
			i__3 = *n - i__ + 1;
#line 1346 "sgesdd.f"
			blk = min(i__3,chunk);
#line 1347 "sgesdd.f"
			sgemm_("N", "N", m, &blk, m, &c_b248, &work[ivt], &
				ldwkvt, &a[i__ * a_dim1 + 1], lda, &c_b227, &
				work[il], m, (ftnlen)1, (ftnlen)1);
#line 1350 "sgesdd.f"
			slacpy_("F", m, &blk, &work[il], m, &a[i__ * a_dim1 + 
				1], lda, (ftnlen)1);
#line 1352 "sgesdd.f"
/* L40: */
#line 1352 "sgesdd.f"
		    }
#line 1353 "sgesdd.f"
		}
#line 1354 "sgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need M+BDSPAC) */

#line 1361 "sgesdd.f"
		slaset_("F", m, n, &c_b227, &c_b227, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1362 "sgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need 3*M, prefer 2*M+M*NB) */

#line 1370 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1370 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1373 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1373 "sgesdd.f"
		sormbr_("P", "R", "T", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1376 "sgesdd.f"
	    } else if (wntqa) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in U and computing right singular */
/*              vectors of bidiagonal matrix in VT */
/*              (Workspace: need M+BDSPAC) */

#line 1383 "sgesdd.f"
		slaset_("F", n, n, &c_b227, &c_b227, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1384 "sgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of VT to identity matrix */

#line 1390 "sgesdd.f"
		if (*n > *m) {
#line 1391 "sgesdd.f"
		    i__1 = *n - *m;
#line 1391 "sgesdd.f"
		    i__2 = *n - *m;
#line 1391 "sgesdd.f"
		    slaset_("F", &i__1, &i__2, &c_b227, &c_b248, &vt[*m + 1 + 
			    (*m + 1) * vt_dim1], ldvt, (ftnlen)1);
#line 1393 "sgesdd.f"
		}

/*              Overwrite U by left singular vectors of A and VT */
/*              by right singular vectors of A */
/*              (Workspace: need 2*M+N, prefer 2*M+N*NB) */

#line 1399 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1399 "sgesdd.f"
		sormbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1402 "sgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1402 "sgesdd.f"
		sormbr_("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1405 "sgesdd.f"
	    }

#line 1407 "sgesdd.f"
	}

#line 1409 "sgesdd.f"
    }

/*     Undo scaling if necessary */

#line 1413 "sgesdd.f"
    if (iscl == 1) {
#line 1414 "sgesdd.f"
	if (anrm > bignum) {
#line 1414 "sgesdd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1414 "sgesdd.f"
	}
#line 1417 "sgesdd.f"
	if (anrm < smlnum) {
#line 1417 "sgesdd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 1417 "sgesdd.f"
	}
#line 1420 "sgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 1424 "sgesdd.f"
    work[1] = (doublereal) maxwrk;

#line 1426 "sgesdd.f"
    return 0;

/*     End of SGESDD */

} /* sgesdd_ */

