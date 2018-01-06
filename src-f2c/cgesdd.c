#line 1 "cgesdd.f"
/* cgesdd.f -- translated by f2c (version 20100827).
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

#line 1 "cgesdd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;

/* > \brief \b CGESDD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESDD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesdd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesdd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesdd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), S( * ) */
/*       COMPLEX            A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESDD computes the singular value decomposition (SVD) of a complex */
/* > M-by-N matrix A, optionally computing the left and/or right singular */
/* > vectors, by using divide-and-conquer method. The SVD is written */
/* > */
/* >      A = U * SIGMA * conjugate-transpose(V) */
/* > */
/* > where SIGMA is an M-by-N matrix which is zero except for its */
/* > min(m,n) diagonal elements, U is an M-by-M unitary matrix, and */
/* > V is an N-by-N unitary matrix.  The diagonal elements of SIGMA */
/* > are the singular values of A; they are real and non-negative, and */
/* > are returned in descending order.  The first min(m,n) columns of */
/* > U and V are the left and right singular vectors of A. */
/* > */
/* > Note that the routine returns VT = V**H, not V. */
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
/* >          = 'A':  all M columns of U and all N rows of V**H are */
/* >                  returned in the arrays U and VT; */
/* >          = 'S':  the first min(M,N) columns of U and the first */
/* >                  min(M,N) rows of V**H are returned in the arrays U */
/* >                  and VT; */
/* >          = 'O':  If M >= N, the first N columns of U are overwritten */
/* >                  in the array A and all rows of V**H are returned in */
/* >                  the array VT; */
/* >                  otherwise, all columns of U are returned in the */
/* >                  array U and the first M rows of V**H are overwritten */
/* >                  in the array A; */
/* >          = 'N':  no columns of U or rows of V**H are computed. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          if JOBZ = 'O',  A is overwritten with the first N columns */
/* >                          of U (the left singular vectors, stored */
/* >                          columnwise) if M >= N; */
/* >                          A is overwritten with the first M rows */
/* >                          of V**H (the right singular vectors, stored */
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
/* >          U is COMPLEX array, dimension (LDU,UCOL) */
/* >          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N; */
/* >          UCOL = min(M,N) if JOBZ = 'S'. */
/* >          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M */
/* >          unitary matrix U; */
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
/* >          VT is COMPLEX array, dimension (LDVT,N) */
/* >          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the */
/* >          N-by-N unitary matrix V**H; */
/* >          if JOBZ = 'S', VT contains the first min(M,N) rows of */
/* >          V**H (the right singular vectors, stored rowwise); */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= 1. */
/* >          if JOBZ = 'N', LWORK >= 2*min(M,N)+max(M,N). */
/* >          if JOBZ = 'O', */
/* >                LWORK >= 2*min(M,N)*min(M,N)+2*min(M,N)+max(M,N). */
/* >          if JOBZ = 'S' or 'A', */
/* >                LWORK >= min(M,N)*min(M,N)+2*min(M,N)+max(M,N). */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, a workspace query is assumed.  The optimal */
/* >          size for the WORK array is calculated and stored in WORK(1), */
/* >          and no other work except argument checking is performed. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* >          If JOBZ = 'N', LRWORK >= 5*min(M,N). */
/* >          Otherwise, */
/* >          LRWORK >= min(M,N)*max(5*min(M,N)+7,2*max(M,N)+2*min(M,N)+1) */
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
/* >          > 0:  The updating process of SBDSDC did not converge. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexGEsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgesdd_(char *jobz, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, 
	integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info, 
	ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ie, il, ir, iu, blk;
    static doublereal dum[1], eps;
    static integer iru, ivt, iscl;
    static doublereal anrm;
    static integer idum[1], ierr, itau, irvt;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn, wrkbl, itaup, itauq;
    static logical wntqa;
    static integer nwork;
    extern /* Subroutine */ int clacp2_(char *, integer *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, ftnlen);
    static logical wntqn, wntqo, wntqs;
    static integer mnthr1, mnthr2;
    extern /* Subroutine */ int cgebrd_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clacrm_(integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, doublereal *)
	    , clarcm_(integer *, integer *, doublereal *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *), clascl_(char *, integer *, integer *, doublereal *,
	     doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen), sbdsdc_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), cgeqrf_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int cungbr_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), cunmbr_(char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), cunglq_(integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, integer *);
    static integer ldwrkl;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer ldwrkr, minwrk, ldwrku, maxwrk, ldwkvt;
    static doublereal smlnum;
    static logical wntqas;
    static integer nrwork;


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 282 "cgesdd.f"
    /* Parameter adjustments */
#line 282 "cgesdd.f"
    a_dim1 = *lda;
#line 282 "cgesdd.f"
    a_offset = 1 + a_dim1;
#line 282 "cgesdd.f"
    a -= a_offset;
#line 282 "cgesdd.f"
    --s;
#line 282 "cgesdd.f"
    u_dim1 = *ldu;
#line 282 "cgesdd.f"
    u_offset = 1 + u_dim1;
#line 282 "cgesdd.f"
    u -= u_offset;
#line 282 "cgesdd.f"
    vt_dim1 = *ldvt;
#line 282 "cgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 282 "cgesdd.f"
    vt -= vt_offset;
#line 282 "cgesdd.f"
    --work;
#line 282 "cgesdd.f"
    --rwork;
#line 282 "cgesdd.f"
    --iwork;
#line 282 "cgesdd.f"

#line 282 "cgesdd.f"
    /* Function Body */
#line 282 "cgesdd.f"
    *info = 0;
#line 283 "cgesdd.f"
    minmn = min(*m,*n);
#line 284 "cgesdd.f"
    mnthr1 = (integer) (minmn * 17. / 9.);
#line 285 "cgesdd.f"
    mnthr2 = (integer) (minmn * 5. / 3.);
#line 286 "cgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 287 "cgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 288 "cgesdd.f"
    wntqas = wntqa || wntqs;
#line 289 "cgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 290 "cgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 291 "cgesdd.f"
    minwrk = 1;
#line 292 "cgesdd.f"
    maxwrk = 1;

#line 294 "cgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 295 "cgesdd.f"
	*info = -1;
#line 296 "cgesdd.f"
    } else if (*m < 0) {
#line 297 "cgesdd.f"
	*info = -2;
#line 298 "cgesdd.f"
    } else if (*n < 0) {
#line 299 "cgesdd.f"
	*info = -3;
#line 300 "cgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 301 "cgesdd.f"
	*info = -5;
#line 302 "cgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 304 "cgesdd.f"
	*info = -8;
#line 305 "cgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 308 "cgesdd.f"
	*info = -10;
#line 309 "cgesdd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to */
/*       real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 319 "cgesdd.f"
    if (*info == 0 && *m > 0 && *n > 0) {
#line 320 "cgesdd.f"
	if (*m >= *n) {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD is BDSPAC */
/*           for computing singular values and singular vectors; BDSPAN */
/*           for computing singular values only. */
/*           BDSPAC = 5*N*N + 7*N */
/*           BDSPAN = MAX(7*N+4, 3*N+2+SMLSIZ*(SMLSIZ+8)) */

#line 329 "cgesdd.f"
	    if (*m >= mnthr1) {
#line 330 "cgesdd.f"
		if (wntqn) {

/*                 Path 1 (M much larger than N, JOBZ='N') */

#line 334 "cgesdd.f"
		    maxwrk = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 336 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 336 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 338 "cgesdd.f"
		    minwrk = *n * 3;
#line 339 "cgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M much larger than N, JOBZ='O') */

#line 343 "cgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 344 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "CUNGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 344 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 346 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 346 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 348 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 348 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 350 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 350 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 352 "cgesdd.f"
		    maxwrk = *m * *n + *n * *n + wrkbl;
#line 353 "cgesdd.f"
		    minwrk = (*n << 1) * *n + *n * 3;
#line 354 "cgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M much larger than N, JOBZ='S') */

#line 358 "cgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 359 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "CUNGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 359 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 361 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 361 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 363 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 363 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 365 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 365 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 367 "cgesdd.f"
		    maxwrk = *n * *n + wrkbl;
#line 368 "cgesdd.f"
		    minwrk = *n * *n + *n * 3;
#line 369 "cgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M much larger than N, JOBZ='A') */

#line 373 "cgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 374 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *m * ilaenv_(&c__1, "CUNGQR", 
			    " ", m, m, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 374 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 376 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 376 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 378 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 378 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 380 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 380 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 382 "cgesdd.f"
		    maxwrk = *n * *n + wrkbl;
#line 383 "cgesdd.f"
		    minwrk = *n * *n + (*n << 1) + *m;
#line 384 "cgesdd.f"
		}
#line 385 "cgesdd.f"
	    } else if (*m >= mnthr2) {

/*              Path 5 (M much larger than N, but not as much as MNTHR1) */

#line 389 "cgesdd.f"
		maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 391 "cgesdd.f"
		minwrk = (*n << 1) + *m;
#line 392 "cgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 393 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "P", n, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 393 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 395 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "Q", m, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 395 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 397 "cgesdd.f"
		    maxwrk += *m * *n;
#line 398 "cgesdd.f"
		    minwrk += *n * *n;
#line 399 "cgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 400 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "P", n, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 400 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 402 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "Q", m, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 402 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 404 "cgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 405 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "P", n, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 405 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 407 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 407 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 409 "cgesdd.f"
		}
#line 410 "cgesdd.f"
	    } else {

/*              Path 6 (M at least N, but not much larger) */

#line 414 "cgesdd.f"
		maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 416 "cgesdd.f"
		minwrk = (*n << 1) + *m;
#line 417 "cgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 418 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 418 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 420 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", m, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 420 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 422 "cgesdd.f"
		    maxwrk += *m * *n;
#line 423 "cgesdd.f"
		    minwrk += *n * *n;
#line 424 "cgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 425 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 425 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 427 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", m, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 427 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 429 "cgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 430 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 430 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 432 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 432 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 434 "cgesdd.f"
		}
#line 435 "cgesdd.f"
	    }
#line 436 "cgesdd.f"
	} else {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD is BDSPAC */
/*           for computing singular values and singular vectors; BDSPAN */
/*           for computing singular values only. */
/*           BDSPAC = 5*M*M + 7*M */
/*           BDSPAN = MAX(7*M+4, 3*M+2+SMLSIZ*(SMLSIZ+8)) */

#line 445 "cgesdd.f"
	    if (*n >= mnthr1) {
#line 446 "cgesdd.f"
		if (wntqn) {

/*                 Path 1t (N much larger than M, JOBZ='N') */

#line 450 "cgesdd.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "CGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 452 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 452 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 454 "cgesdd.f"
		    minwrk = *m * 3;
#line 455 "cgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N much larger than M, JOBZ='O') */

#line 459 "cgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "CGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 460 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "CUNGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 460 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 462 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 462 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 464 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 464 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 466 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 466 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 468 "cgesdd.f"
		    maxwrk = *m * *n + *m * *m + wrkbl;
#line 469 "cgesdd.f"
		    minwrk = (*m << 1) * *m + *m * 3;
#line 470 "cgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N much larger than M, JOBZ='S') */

#line 474 "cgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "CGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 475 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "CUNGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 475 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 477 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 477 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 479 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 479 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 481 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 481 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 483 "cgesdd.f"
		    maxwrk = *m * *m + wrkbl;
#line 484 "cgesdd.f"
		    minwrk = *m * *m + *m * 3;
#line 485 "cgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N much larger than M, JOBZ='A') */

#line 489 "cgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "CGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 490 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *n * ilaenv_(&c__1, "CUNGLQ", 
			    " ", n, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 490 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 492 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "CGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 492 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 494 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 494 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 496 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 496 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 498 "cgesdd.f"
		    maxwrk = *m * *m + wrkbl;
#line 499 "cgesdd.f"
		    minwrk = *m * *m + (*m << 1) + *n;
#line 500 "cgesdd.f"
		}
#line 501 "cgesdd.f"
	    } else if (*n >= mnthr2) {

/*              Path 5t (N much larger than M, but not as much as MNTHR1) */

#line 505 "cgesdd.f"
		maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 507 "cgesdd.f"
		minwrk = (*m << 1) + *n;
#line 508 "cgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 509 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "P", m, n, m, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 509 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 511 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 511 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 513 "cgesdd.f"
		    maxwrk += *m * *n;
#line 514 "cgesdd.f"
		    minwrk += *m * *m;
#line 515 "cgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 516 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "P", m, n, m, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 516 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 518 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 518 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 520 "cgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 521 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "P", n, n, m, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 521 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 523 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 523 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 525 "cgesdd.f"
		}
#line 526 "cgesdd.f"
	    } else {

/*              Path 6t (N greater than M, but not much larger) */

#line 530 "cgesdd.f"
		maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 532 "cgesdd.f"
		minwrk = (*m << 1) + *n;
#line 533 "cgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 534 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "PRC", m, n, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 534 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 536 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 536 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 538 "cgesdd.f"
		    maxwrk += *m * *n;
#line 539 "cgesdd.f"
		    minwrk += *m * *m;
#line 540 "cgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 541 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "PRC", m, n, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 541 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 543 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 543 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 545 "cgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 546 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *n * ilaenv_(&c__1, 
			    "CUNGBR", "PRC", n, n, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 546 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 548 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNGBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 548 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 550 "cgesdd.f"
		}
#line 551 "cgesdd.f"
	    }
#line 552 "cgesdd.f"
	}
#line 553 "cgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 554 "cgesdd.f"
    }
#line 555 "cgesdd.f"
    if (*info == 0) {
#line 556 "cgesdd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 557 "cgesdd.f"
	if (*lwork < minwrk && *lwork != -1) {
#line 557 "cgesdd.f"
	    *info = -13;
#line 557 "cgesdd.f"
	}
#line 559 "cgesdd.f"
    }

/*     Quick returns */

#line 563 "cgesdd.f"
    if (*info != 0) {
#line 564 "cgesdd.f"
	i__1 = -(*info);
#line 564 "cgesdd.f"
	xerbla_("CGESDD", &i__1, (ftnlen)6);
#line 565 "cgesdd.f"
	return 0;
#line 566 "cgesdd.f"
    }
#line 567 "cgesdd.f"
    if (*lwork == -1) {
#line 567 "cgesdd.f"
	return 0;
#line 567 "cgesdd.f"
    }
#line 569 "cgesdd.f"
    if (*m == 0 || *n == 0) {
#line 570 "cgesdd.f"
	return 0;
#line 571 "cgesdd.f"
    }

/*     Get machine constants */

#line 575 "cgesdd.f"
    eps = slamch_("P", (ftnlen)1);
#line 576 "cgesdd.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 577 "cgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 581 "cgesdd.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 582 "cgesdd.f"
    iscl = 0;
#line 583 "cgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 584 "cgesdd.f"
	iscl = 1;
#line 585 "cgesdd.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 586 "cgesdd.f"
    } else if (anrm > bignum) {
#line 587 "cgesdd.f"
	iscl = 1;
#line 588 "cgesdd.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 589 "cgesdd.f"
    }

#line 591 "cgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 597 "cgesdd.f"
	if (*m >= mnthr1) {

#line 599 "cgesdd.f"
	    if (wntqn) {

/*              Path 1 (M much larger than N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 604 "cgesdd.f"
		itau = 1;
#line 605 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: need 0) */

#line 611 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 611 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 616 "cgesdd.f"
		i__1 = *n - 1;
#line 616 "cgesdd.f"
		i__2 = *n - 1;
#line 616 "cgesdd.f"
		claset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 618 "cgesdd.f"
		ie = 1;
#line 619 "cgesdd.f"
		itauq = 1;
#line 620 "cgesdd.f"
		itaup = itauq + *n;
#line 621 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 627 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 627 "cgesdd.f"
		cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 630 "cgesdd.f"
		nrwork = ie + *n;

/*              Perform bidiagonal SVD, compute singular values only */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAN) */

#line 636 "cgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 639 "cgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M much larger than N, JOBZ='O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 645 "cgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 649 "cgesdd.f"
		ldwrku = *n;
#line 650 "cgesdd.f"
		ir = iu + ldwrku * *n;
#line 651 "cgesdd.f"
		if (*lwork >= *m * *n + *n * *n + *n * 3) {

/*                 WORK(IR) is M by N */

#line 655 "cgesdd.f"
		    ldwrkr = *m;
#line 656 "cgesdd.f"
		} else {
#line 657 "cgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 658 "cgesdd.f"
		}
#line 659 "cgesdd.f"
		itau = ir + ldwrkr * *n;
#line 660 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need N*N+2*N, prefer M*N+N+N*NB) */
/*              (RWorkspace: 0) */

#line 666 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 666 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK( IR ), zeroing out below it */

#line 671 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 672 "cgesdd.f"
		i__1 = *n - 1;
#line 672 "cgesdd.f"
		i__2 = *n - 1;
#line 672 "cgesdd.f"
		claset_("L", &i__1, &i__2, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 679 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 679 "cgesdd.f"
		cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 681 "cgesdd.f"
		ie = 1;
#line 682 "cgesdd.f"
		itauq = itau;
#line 683 "cgesdd.f"
		itaup = itauq + *n;
#line 684 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 690 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 690 "cgesdd.f"
		cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of R in WORK(IRU) and computing right singular vectors */
/*              of R in WORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 700 "cgesdd.f"
		iru = ie + *n;
#line 701 "cgesdd.f"
		irvt = iru + *n * *n;
#line 702 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 703 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of R */
/*              (CWorkspace: need 2*N*N+3*N, prefer M*N+N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 712 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 714 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 714 "cgesdd.f"
		cunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by the right singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 723 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 724 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 724 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              (CWorkspace: need 2*N*N, prefer N*N+M*N) */
/*              (RWorkspace: 0) */

#line 733 "cgesdd.f"
		i__1 = *m;
#line 733 "cgesdd.f"
		i__2 = ldwrkr;
#line 733 "cgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 734 "cgesdd.f"
		    i__3 = *m - i__ + 1;
#line 734 "cgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 735 "cgesdd.f"
		    cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1], 
			    lda, &work[iu], &ldwrku, &c_b1, &work[ir], &
			    ldwrkr, (ftnlen)1, (ftnlen)1);
#line 738 "cgesdd.f"
		    clacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 740 "cgesdd.f"
/* L10: */
#line 740 "cgesdd.f"
		}

#line 742 "cgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M much larger than N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 748 "cgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 752 "cgesdd.f"
		ldwrkr = *n;
#line 753 "cgesdd.f"
		itau = ir + ldwrkr * *n;
#line 754 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*              (RWorkspace: 0) */

#line 760 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 760 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 765 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 766 "cgesdd.f"
		i__2 = *n - 1;
#line 766 "cgesdd.f"
		i__1 = *n - 1;
#line 766 "cgesdd.f"
		claset_("L", &i__2, &i__1, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 773 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 773 "cgesdd.f"
		cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 775 "cgesdd.f"
		ie = 1;
#line 776 "cgesdd.f"
		itauq = itau;
#line 777 "cgesdd.f"
		itaup = itauq + *n;
#line 778 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 784 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 784 "cgesdd.f"
		cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 794 "cgesdd.f"
		iru = ie + *n;
#line 795 "cgesdd.f"
		irvt = iru + *n * *n;
#line 796 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 797 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 806 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 807 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 807 "cgesdd.f"
		cunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 816 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 817 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 817 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              (CWorkspace: need N*N) */
/*              (RWorkspace: 0) */

#line 826 "cgesdd.f"
		clacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 827 "cgesdd.f"
		cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &work[ir],
			 &ldwrkr, &c_b1, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 830 "cgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M much larger than N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 836 "cgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 840 "cgesdd.f"
		ldwrku = *n;
#line 841 "cgesdd.f"
		itau = iu + ldwrku * *n;
#line 842 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 848 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 848 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 850 "cgesdd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              (CWorkspace: need N+M, prefer N+M*NB) */
/*              (RWorkspace: 0) */

#line 856 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 856 "cgesdd.f"
		cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out below it */

#line 861 "cgesdd.f"
		i__2 = *n - 1;
#line 861 "cgesdd.f"
		i__1 = *n - 1;
#line 861 "cgesdd.f"
		claset_("L", &i__2, &i__1, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 863 "cgesdd.f"
		ie = 1;
#line 864 "cgesdd.f"
		itauq = itau;
#line 865 "cgesdd.f"
		itaup = itauq + *n;
#line 866 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 872 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 872 "cgesdd.f"
		cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 875 "cgesdd.f"
		iru = ie + *n;
#line 876 "cgesdd.f"
		irvt = iru + *n * *n;
#line 877 "cgesdd.f"
		nrwork = irvt + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 885 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by left singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 894 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 896 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 896 "cgesdd.f"
		cunmbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 905 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 906 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 906 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              (CWorkspace: need N*N) */
/*              (RWorkspace: 0) */

#line 915 "cgesdd.f"
		cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &work[iu],
			 &ldwrku, &c_b1, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 920 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 922 "cgesdd.f"
	    }

#line 924 "cgesdd.f"
	} else if (*m >= mnthr2) {

/*           MNTHR2 <= M < MNTHR1 */

/*           Path 5 (M much larger than N, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           CUNGBR and matrix multiplication to compute singular vectors */

#line 932 "cgesdd.f"
	    ie = 1;
#line 933 "cgesdd.f"
	    nrwork = ie + *n;
#line 934 "cgesdd.f"
	    itauq = 1;
#line 935 "cgesdd.f"
	    itaup = itauq + *n;
#line 936 "cgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 942 "cgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 942 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 945 "cgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 951 "cgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 953 "cgesdd.f"
	    } else if (wntqo) {
#line 954 "cgesdd.f"
		iu = nwork;
#line 955 "cgesdd.f"
		iru = nrwork;
#line 956 "cgesdd.f"
		irvt = iru + *n * *n;
#line 957 "cgesdd.f"
		nrwork = irvt + *n * *n;

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 963 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 964 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 964 "cgesdd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 971 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 971 "cgesdd.f"
		cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

#line 974 "cgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 978 "cgesdd.f"
		    ldwrku = *m;
#line 979 "cgesdd.f"
		} else {

/*                 WORK(IU) is LDWRKU by N */

#line 983 "cgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 984 "cgesdd.f"
		}
#line 985 "cgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 993 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in WORK(IU), copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1002 "cgesdd.f"
		clarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &work[iu]
			, &ldwrku, &rwork[nrwork]);
#line 1004 "cgesdd.f"
		clacpy_("F", n, n, &work[iu], &ldwrku, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in A by real matrix RWORK(IRU), storing the */
/*              result in WORK(IU), copying to A */
/*              (CWorkspace: need N*N, prefer M*N) */
/*              (Rworkspace: need 3*N*N, prefer N*N+2*M*N) */

#line 1011 "cgesdd.f"
		nrwork = irvt;
#line 1012 "cgesdd.f"
		i__2 = *m;
#line 1012 "cgesdd.f"
		i__1 = ldwrku;
#line 1012 "cgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1013 "cgesdd.f"
		    i__3 = *m - i__ + 1;
#line 1013 "cgesdd.f"
		    chunk = min(i__3,ldwrku);
#line 1014 "cgesdd.f"
		    clacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru], n, 
			    &work[iu], &ldwrku, &rwork[nrwork]);
#line 1016 "cgesdd.f"
		    clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 1018 "cgesdd.f"
/* L20: */
#line 1018 "cgesdd.f"
		}

#line 1020 "cgesdd.f"
	    } else if (wntqs) {

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1026 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1027 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1027 "cgesdd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1034 "cgesdd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1035 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1035 "cgesdd.f"
		cungbr_("Q", m, n, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1044 "cgesdd.f"
		iru = nrwork;
#line 1045 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1046 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1047 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1056 "cgesdd.f"
		clarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1058 "cgesdd.f"
		clacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: need 0) */
/*              (Rworkspace: need N*N+2*M*N) */

#line 1065 "cgesdd.f"
		nrwork = irvt;
#line 1066 "cgesdd.f"
		clacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1068 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1069 "cgesdd.f"
	    } else {

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1075 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1076 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1076 "cgesdd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1083 "cgesdd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1084 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1084 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1093 "cgesdd.f"
		iru = nrwork;
#line 1094 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1095 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1096 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1105 "cgesdd.f"
		clarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1107 "cgesdd.f"
		clacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1114 "cgesdd.f"
		nrwork = irvt;
#line 1115 "cgesdd.f"
		clacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1117 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1118 "cgesdd.f"
	    }

#line 1120 "cgesdd.f"
	} else {

/*           M .LT. MNTHR2 */

/*           Path 6 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */
/*           Use CUNMBR to compute singular vectors */

#line 1128 "cgesdd.f"
	    ie = 1;
#line 1129 "cgesdd.f"
	    nrwork = ie + *n;
#line 1130 "cgesdd.f"
	    itauq = 1;
#line 1131 "cgesdd.f"
	    itaup = itauq + *n;
#line 1132 "cgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 1138 "cgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1138 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);
#line 1141 "cgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 1147 "cgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1149 "cgesdd.f"
	    } else if (wntqo) {
#line 1150 "cgesdd.f"
		iu = nwork;
#line 1151 "cgesdd.f"
		iru = nrwork;
#line 1152 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1153 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1154 "cgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 1158 "cgesdd.f"
		    ldwrku = *m;
#line 1159 "cgesdd.f"
		} else {

/*                 WORK( IU ) is LDWRKU by N */

#line 1163 "cgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 1164 "cgesdd.f"
		}
#line 1165 "cgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1173 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: need 0) */

#line 1182 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1183 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1183 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1187 "cgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by left singular vectors of A, copying */
/*              to A */
/*              (Cworkspace: need M*N+2*N, prefer M*N+N+N*NB) */
/*              (Rworkspace: need 0) */

#line 1195 "cgesdd.f"
		    claset_("F", m, n, &c_b1, &c_b1, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1197 "cgesdd.f"
		    clacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1199 "cgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1199 "cgesdd.f"
		    cunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1202 "cgesdd.f"
		    clacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 1203 "cgesdd.f"
		} else {

/*                 Generate Q in A */
/*                 (Cworkspace: need 2*N, prefer N+N*NB) */
/*                 (Rworkspace: need 0) */

#line 1209 "cgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1209 "cgesdd.f"
		    cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__1, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 (CWorkspace: need N*N, prefer M*N) */
/*                 (Rworkspace: need 3*N*N, prefer N*N+2*M*N) */

#line 1217 "cgesdd.f"
		    nrwork = irvt;
#line 1218 "cgesdd.f"
		    i__1 = *m;
#line 1218 "cgesdd.f"
		    i__2 = ldwrku;
#line 1218 "cgesdd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ +=
			     i__2) {
/* Computing MIN */
#line 1219 "cgesdd.f"
			i__3 = *m - i__ + 1;
#line 1219 "cgesdd.f"
			chunk = min(i__3,ldwrku);
#line 1220 "cgesdd.f"
			clacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru],
				 n, &work[iu], &ldwrku, &rwork[nrwork]);
#line 1223 "cgesdd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 1225 "cgesdd.f"
/* L30: */
#line 1225 "cgesdd.f"
		    }
#line 1226 "cgesdd.f"
		}

#line 1228 "cgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1236 "cgesdd.f"
		iru = nrwork;
#line 1237 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1238 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1239 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1248 "cgesdd.f"
		claset_("F", m, n, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1249 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1250 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1250 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1259 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1260 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1260 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1263 "cgesdd.f"
	    } else {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1271 "cgesdd.f"
		iru = nrwork;
#line 1272 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1273 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1274 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 1280 "cgesdd.f"
		claset_("F", m, m, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1281 "cgesdd.f"
		if (*m > *n) {
#line 1282 "cgesdd.f"
		    i__2 = *m - *n;
#line 1282 "cgesdd.f"
		    i__1 = *m - *n;
#line 1282 "cgesdd.f"
		    claset_("F", &i__2, &i__1, &c_b1, &c_b2, &u[*n + 1 + (*n 
			    + 1) * u_dim1], ldu, (ftnlen)1);
#line 1284 "cgesdd.f"
		}

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*              (RWorkspace: 0) */

#line 1291 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1292 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1292 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1301 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1302 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1302 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1305 "cgesdd.f"
	    }

#line 1307 "cgesdd.f"
	}

#line 1309 "cgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 1315 "cgesdd.f"
	if (*n >= mnthr1) {

#line 1317 "cgesdd.f"
	    if (wntqn) {

/*              Path 1t (N much larger than M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 1322 "cgesdd.f"
		itau = 1;
#line 1323 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1329 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1329 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 1334 "cgesdd.f"
		i__2 = *m - 1;
#line 1334 "cgesdd.f"
		i__1 = *m - 1;
#line 1334 "cgesdd.f"
		claset_("U", &i__2, &i__1, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1336 "cgesdd.f"
		ie = 1;
#line 1337 "cgesdd.f"
		itauq = 1;
#line 1338 "cgesdd.f"
		itaup = itauq + *m;
#line 1339 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1345 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1345 "cgesdd.f"
		cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 1348 "cgesdd.f"
		nrwork = ie + *m;

/*              Perform bidiagonal SVD, compute singular values only */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAN) */

#line 1354 "cgesdd.f"
		sbdsdc_("U", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 1357 "cgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N much larger than M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1363 "cgesdd.f"
		ivt = 1;
#line 1364 "cgesdd.f"
		ldwkvt = *m;

/*              WORK(IVT) is M by M */

#line 1368 "cgesdd.f"
		il = ivt + ldwkvt * *m;
#line 1369 "cgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3) {

/*                 WORK(IL) M by N */

#line 1373 "cgesdd.f"
		    ldwrkl = *m;
#line 1374 "cgesdd.f"
		    chunk = *n;
#line 1375 "cgesdd.f"
		} else {

/*                 WORK(IL) is M by CHUNK */

#line 1379 "cgesdd.f"
		    ldwrkl = *m;
#line 1380 "cgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1381 "cgesdd.f"
		}
#line 1382 "cgesdd.f"
		itau = il + ldwrkl * chunk;
#line 1383 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1389 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1389 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1394 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1395 "cgesdd.f"
		i__2 = *m - 1;
#line 1395 "cgesdd.f"
		i__1 = *m - 1;
#line 1395 "cgesdd.f"
		claset_("U", &i__2, &i__1, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*              (RWorkspace: 0) */

#line 1402 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1402 "cgesdd.f"
		cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1404 "cgesdd.f"
		ie = 1;
#line 1405 "cgesdd.f"
		itauq = itau;
#line 1406 "cgesdd.f"
		itaup = itauq + *m;
#line 1407 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1413 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1413 "cgesdd.f"
		cgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1423 "cgesdd.f"
		iru = ie + *m;
#line 1424 "cgesdd.f"
		irvt = iru + *m * *m;
#line 1425 "cgesdd.f"
		nrwork = irvt + *m * *m;
#line 1426 "cgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of L */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1435 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1436 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1436 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by the right singular vectors of L */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1445 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1447 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1447 "cgesdd.f"
		cunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              (CWorkspace: need 2*M*M, prefer M*M+M*N)) */
/*              (RWorkspace: 0) */

#line 1456 "cgesdd.f"
		i__2 = *n;
#line 1456 "cgesdd.f"
		i__1 = chunk;
#line 1456 "cgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1457 "cgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1457 "cgesdd.f"
		    blk = min(i__3,chunk);
#line 1458 "cgesdd.f"
		    cgemm_("N", "N", m, &blk, m, &c_b2, &work[ivt], m, &a[i__ 
			    * a_dim1 + 1], lda, &c_b1, &work[il], &ldwrkl, (
			    ftnlen)1, (ftnlen)1);
#line 1461 "cgesdd.f"
		    clacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1463 "cgesdd.f"
/* L40: */
#line 1463 "cgesdd.f"
		}

#line 1465 "cgesdd.f"
	    } else if (wntqs) {

/*             Path 3t (N much larger than M, JOBZ='S') */
/*             M right singular vectors to be computed in VT and */
/*             M left singular vectors to be computed in U */

#line 1471 "cgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1475 "cgesdd.f"
		ldwrkl = *m;
#line 1476 "cgesdd.f"
		itau = il + ldwrkl * *m;
#line 1477 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1483 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1483 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1488 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1489 "cgesdd.f"
		i__1 = *m - 1;
#line 1489 "cgesdd.f"
		i__2 = *m - 1;
#line 1489 "cgesdd.f"
		claset_("U", &i__1, &i__2, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*              (RWorkspace: 0) */

#line 1496 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1496 "cgesdd.f"
		cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1498 "cgesdd.f"
		ie = 1;
#line 1499 "cgesdd.f"
		itauq = itau;
#line 1500 "cgesdd.f"
		itaup = itauq + *m;
#line 1501 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1507 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1507 "cgesdd.f"
		cgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1517 "cgesdd.f"
		iru = ie + *m;
#line 1518 "cgesdd.f"
		irvt = iru + *m * *m;
#line 1519 "cgesdd.f"
		nrwork = irvt + *m * *m;
#line 1520 "cgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1529 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1530 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1530 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by left singular vectors of L */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1539 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1540 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1540 "cgesdd.f"
		cunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy VT to WORK(IL), multiply right singular vectors of L */
/*              in WORK(IL) by Q in A, storing result in VT */
/*              (CWorkspace: need M*M) */
/*              (RWorkspace: 0) */

#line 1549 "cgesdd.f"
		clacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1550 "cgesdd.f"
		cgemm_("N", "N", m, n, m, &c_b2, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b1, &vt[vt_offset], ldvt, (ftnlen)
			1, (ftnlen)1);

#line 1553 "cgesdd.f"
	    } else if (wntqa) {

/*              Path 9t (N much larger than M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1559 "cgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1563 "cgesdd.f"
		ldwkvt = *m;
#line 1564 "cgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1565 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1571 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1571 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
#line 1573 "cgesdd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              (CWorkspace: need M+N, prefer M+N*NB) */
/*              (RWorkspace: 0) */

#line 1579 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1579 "cgesdd.f"
		cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

/*              Produce L in A, zeroing out above it */

#line 1584 "cgesdd.f"
		i__1 = *m - 1;
#line 1584 "cgesdd.f"
		i__2 = *m - 1;
#line 1584 "cgesdd.f"
		claset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1586 "cgesdd.f"
		ie = 1;
#line 1587 "cgesdd.f"
		itauq = itau;
#line 1588 "cgesdd.f"
		itaup = itauq + *m;
#line 1589 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1595 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1595 "cgesdd.f"
		cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1605 "cgesdd.f"
		iru = ie + *m;
#line 1606 "cgesdd.f"
		irvt = iru + *m * *m;
#line 1607 "cgesdd.f"
		nrwork = irvt + *m * *m;
#line 1608 "cgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1617 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1618 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1618 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by right singular vectors of L */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1627 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1629 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1629 "cgesdd.f"
		cunmbr_("P", "R", "C", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              (CWorkspace: need M*M) */
/*              (RWorkspace: 0) */

#line 1638 "cgesdd.f"
		cgemm_("N", "N", m, n, m, &c_b2, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b1, &a[a_offset], lda, (ftnlen)1,
			 (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1643 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1645 "cgesdd.f"
	    }

#line 1647 "cgesdd.f"
	} else if (*n >= mnthr2) {

/*           MNTHR2 <= N < MNTHR1 */

/*           Path 5t (N much larger than M, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           CUNGBR and matrix multiplication to compute singular vectors */


#line 1656 "cgesdd.f"
	    ie = 1;
#line 1657 "cgesdd.f"
	    nrwork = ie + *m;
#line 1658 "cgesdd.f"
	    itauq = 1;
#line 1659 "cgesdd.f"
	    itaup = itauq + *m;
#line 1660 "cgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 1666 "cgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1666 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);

#line 1670 "cgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 1676 "cgesdd.f"
		sbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1678 "cgesdd.f"
	    } else if (wntqo) {
#line 1679 "cgesdd.f"
		irvt = nrwork;
#line 1680 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1681 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1682 "cgesdd.f"
		ivt = nwork;

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1688 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1689 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1689 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Generate P**H in A */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1696 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1696 "cgesdd.f"
		cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

#line 1699 "cgesdd.f"
		ldwkvt = *m;
#line 1700 "cgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 1704 "cgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1705 "cgesdd.f"
		    chunk = *n;
#line 1706 "cgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 1710 "cgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 1711 "cgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 1712 "cgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1720 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRVT) */
/*              storing the result in WORK(IVT), copying to U */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 2*M*M) */

#line 1729 "cgesdd.f"
		clacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &work[ivt], &
			ldwkvt, &rwork[nrwork]);
#line 1731 "cgesdd.f"
		clacpy_("F", m, m, &work[ivt], &ldwkvt, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply RWORK(IRVT) by P**H in A, storing the */
/*              result in WORK(IVT), copying to A */
/*              (CWorkspace: need M*M, prefer M*N) */
/*              (Rworkspace: need 2*M*M, prefer 2*M*N) */

#line 1738 "cgesdd.f"
		nrwork = iru;
#line 1739 "cgesdd.f"
		i__1 = *n;
#line 1739 "cgesdd.f"
		i__2 = chunk;
#line 1739 "cgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1740 "cgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1740 "cgesdd.f"
		    blk = min(i__3,chunk);
#line 1741 "cgesdd.f"
		    clarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1], 
			    lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 1743 "cgesdd.f"
		    clacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
			    a_dim1 + 1], lda, (ftnlen)1);
#line 1745 "cgesdd.f"
/* L50: */
#line 1745 "cgesdd.f"
		}
#line 1746 "cgesdd.f"
	    } else if (wntqs) {

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1752 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1753 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1753 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1760 "cgesdd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1761 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1761 "cgesdd.f"
		cungbr_("P", m, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1770 "cgesdd.f"
		irvt = nrwork;
#line 1771 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1772 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1773 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: need 0) */
/*              (Rworkspace: need 3*M*M) */

#line 1782 "cgesdd.f"
		clacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1784 "cgesdd.f"
		clacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need M*M+2*M*N) */

#line 1791 "cgesdd.f"
		nrwork = iru;
#line 1792 "cgesdd.f"
		clarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1794 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1795 "cgesdd.f"
	    } else {

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1801 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1802 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1802 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1809 "cgesdd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1810 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1810 "cgesdd.f"
		cungbr_("P", n, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1819 "cgesdd.f"
		irvt = nrwork;
#line 1820 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1821 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1822 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: need 0) */
/*              (Rworkspace: need 3*M*M) */

#line 1831 "cgesdd.f"
		clacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1833 "cgesdd.f"
		clacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need M*M+2*M*N) */

#line 1840 "cgesdd.f"
		clarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1842 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1843 "cgesdd.f"
	    }

#line 1845 "cgesdd.f"
	} else {

/*           N .LT. MNTHR2 */

/*           Path 6t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           Use CUNMBR to compute singular vectors */

#line 1853 "cgesdd.f"
	    ie = 1;
#line 1854 "cgesdd.f"
	    nrwork = ie + *m;
#line 1855 "cgesdd.f"
	    itauq = 1;
#line 1856 "cgesdd.f"
	    itaup = itauq + *m;
#line 1857 "cgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 1863 "cgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1863 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 1866 "cgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 1872 "cgesdd.f"
		sbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1874 "cgesdd.f"
	    } else if (wntqo) {
#line 1875 "cgesdd.f"
		ldwkvt = *m;
#line 1876 "cgesdd.f"
		ivt = nwork;
#line 1877 "cgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 1881 "cgesdd.f"
		    claset_("F", m, n, &c_b1, &c_b1, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 1883 "cgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1884 "cgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 1888 "cgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 1889 "cgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 1890 "cgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1898 "cgesdd.f"
		irvt = nrwork;
#line 1899 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1900 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1901 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: need 0) */

#line 1910 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1911 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1911 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1915 "cgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by right singular vectors of A, */
/*              copying to A */
/*              (Cworkspace: need M*N+2*M, prefer M*N+M+M*NB) */
/*              (Rworkspace: need 0) */

#line 1923 "cgesdd.f"
		    clacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 1925 "cgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1925 "cgesdd.f"
		    cunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1928 "cgesdd.f"
		    clacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 1929 "cgesdd.f"
		} else {

/*                 Generate P**H in A */
/*                 (Cworkspace: need 2*M, prefer M+M*NB) */
/*                 (Rworkspace: need 0) */

#line 1935 "cgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1935 "cgesdd.f"
		    cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 (CWorkspace: need M*M, prefer M*N) */
/*                 (Rworkspace: need 3*M*M, prefer M*M+2*M*N) */

#line 1943 "cgesdd.f"
		    nrwork = iru;
#line 1944 "cgesdd.f"
		    i__2 = *n;
#line 1944 "cgesdd.f"
		    i__1 = chunk;
#line 1944 "cgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1945 "cgesdd.f"
			i__3 = *n - i__ + 1;
#line 1945 "cgesdd.f"
			blk = min(i__3,chunk);
#line 1946 "cgesdd.f"
			clarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1]
				, lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 1949 "cgesdd.f"
			clacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 1951 "cgesdd.f"
/* L60: */
#line 1951 "cgesdd.f"
		    }
#line 1952 "cgesdd.f"
		}
#line 1953 "cgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1961 "cgesdd.f"
		irvt = nrwork;
#line 1962 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1963 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1964 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: M*M) */

#line 1973 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1974 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1974 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: M*M) */

#line 1983 "cgesdd.f"
		claset_("F", m, n, &c_b1, &c_b1, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1984 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1985 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1985 "cgesdd.f"
		cunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1988 "cgesdd.f"
	    } else {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1996 "cgesdd.f"
		irvt = nrwork;
#line 1997 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1998 "cgesdd.f"
		nrwork = iru + *m * *m;

#line 2000 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: M*M) */

#line 2009 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2010 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2010 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Set all of VT to identity matrix */

#line 2016 "cgesdd.f"
		claset_("F", n, n, &c_b1, &c_b2, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*              (RWorkspace: M*M) */

#line 2023 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2024 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2024 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2027 "cgesdd.f"
	    }

#line 2029 "cgesdd.f"
	}

#line 2031 "cgesdd.f"
    }

/*     Undo scaling if necessary */

#line 2035 "cgesdd.f"
    if (iscl == 1) {
#line 2036 "cgesdd.f"
	if (anrm > bignum) {
#line 2036 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2036 "cgesdd.f"
	}
#line 2039 "cgesdd.f"
	if (*info != 0 && anrm > bignum) {
#line 2039 "cgesdd.f"
	    i__1 = minmn - 1;
#line 2039 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2039 "cgesdd.f"
	}
#line 2042 "cgesdd.f"
	if (anrm < smlnum) {
#line 2042 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2042 "cgesdd.f"
	}
#line 2045 "cgesdd.f"
	if (*info != 0 && anrm < smlnum) {
#line 2045 "cgesdd.f"
	    i__1 = minmn - 1;
#line 2045 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2045 "cgesdd.f"
	}
#line 2048 "cgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 2052 "cgesdd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 2054 "cgesdd.f"
    return 0;

/*     End of CGESDD */

} /* cgesdd_ */

