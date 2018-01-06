#line 1 "zgesdd.f"
/* zgesdd.f -- translated by f2c (version 20100827).
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

#line 1 "zgesdd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;

/* > \brief \b ZGESDD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGESDD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesdd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesdd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesdd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, */
/*                          LWORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGESDD computes the singular value decomposition (SVD) of a complex */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX*16 array, dimension (LDU,UCOL) */
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
/* >          VT is COMPLEX*16 array, dimension (LDVT,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
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
/* >          > 0:  The updating process of DBDSDC did not converge. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16GEsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgesdd_(char *jobz, integer *m, integer *n, 
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer wrkbl, itaup, itauq;
    static logical wntqa;
    static integer nwork;
    static logical wntqn, wntqo, wntqs;
    extern /* Subroutine */ int zlacp2_(char *, integer *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, ftnlen);
    static integer mnthr1, mnthr2;
    extern /* Subroutine */ int dbdsdc_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen),
	     zgebrd_(integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zlacrm_(integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, doublereal *)
	    , zlarcm_(integer *, integer *, doublereal *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *), zlascl_(char *, integer *, integer *, doublereal *,
	     doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen), zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static integer ldwrkl;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer ldwrkr, minwrk, ldwrku, maxwrk;
    extern /* Subroutine */ int zungbr_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *, ftnlen);
    static integer ldwkvt;
    static doublereal smlnum;
    static logical wntqas;
    extern /* Subroutine */ int zunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen), zunglq_(integer *, integer *, integer *
	    , doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer nrwork;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);


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

#line 282 "zgesdd.f"
    /* Parameter adjustments */
#line 282 "zgesdd.f"
    a_dim1 = *lda;
#line 282 "zgesdd.f"
    a_offset = 1 + a_dim1;
#line 282 "zgesdd.f"
    a -= a_offset;
#line 282 "zgesdd.f"
    --s;
#line 282 "zgesdd.f"
    u_dim1 = *ldu;
#line 282 "zgesdd.f"
    u_offset = 1 + u_dim1;
#line 282 "zgesdd.f"
    u -= u_offset;
#line 282 "zgesdd.f"
    vt_dim1 = *ldvt;
#line 282 "zgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 282 "zgesdd.f"
    vt -= vt_offset;
#line 282 "zgesdd.f"
    --work;
#line 282 "zgesdd.f"
    --rwork;
#line 282 "zgesdd.f"
    --iwork;
#line 282 "zgesdd.f"

#line 282 "zgesdd.f"
    /* Function Body */
#line 282 "zgesdd.f"
    *info = 0;
#line 283 "zgesdd.f"
    minmn = min(*m,*n);
#line 284 "zgesdd.f"
    mnthr1 = (integer) (minmn * 17. / 9.);
#line 285 "zgesdd.f"
    mnthr2 = (integer) (minmn * 5. / 3.);
#line 286 "zgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 287 "zgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 288 "zgesdd.f"
    wntqas = wntqa || wntqs;
#line 289 "zgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 290 "zgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 291 "zgesdd.f"
    minwrk = 1;
#line 292 "zgesdd.f"
    maxwrk = 1;

#line 294 "zgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 295 "zgesdd.f"
	*info = -1;
#line 296 "zgesdd.f"
    } else if (*m < 0) {
#line 297 "zgesdd.f"
	*info = -2;
#line 298 "zgesdd.f"
    } else if (*n < 0) {
#line 299 "zgesdd.f"
	*info = -3;
#line 300 "zgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 301 "zgesdd.f"
	*info = -5;
#line 302 "zgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 304 "zgesdd.f"
	*info = -8;
#line 305 "zgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 308 "zgesdd.f"
	*info = -10;
#line 309 "zgesdd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to */
/*       real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 319 "zgesdd.f"
    if (*info == 0 && *m > 0 && *n > 0) {
#line 320 "zgesdd.f"
	if (*m >= *n) {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD is BDSPAC */
/*           for computing singular values and singular vectors; BDSPAN */
/*           for computing singular values only. */
/*           BDSPAC = 5*N*N + 7*N */
/*           BDSPAN = MAX(7*N+4, 3*N+2+SMLSIZ*(SMLSIZ+8)) */

#line 329 "zgesdd.f"
	    if (*m >= mnthr1) {
#line 330 "zgesdd.f"
		if (wntqn) {

/*                 Path 1 (M much larger than N, JOBZ='N') */

#line 334 "zgesdd.f"
		    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 336 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 336 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 338 "zgesdd.f"
		    minwrk = *n * 3;
#line 339 "zgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M much larger than N, JOBZ='O') */

#line 343 "zgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 344 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 344 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 346 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 346 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 348 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 348 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 350 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 350 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 352 "zgesdd.f"
		    maxwrk = *m * *n + *n * *n + wrkbl;
#line 353 "zgesdd.f"
		    minwrk = (*n << 1) * *n + *n * 3;
#line 354 "zgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M much larger than N, JOBZ='S') */

#line 358 "zgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 359 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", 
			    " ", m, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 359 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 361 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 361 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 363 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 363 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 365 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 365 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 367 "zgesdd.f"
		    maxwrk = *n * *n + wrkbl;
#line 368 "zgesdd.f"
		    minwrk = *n * *n + *n * 3;
#line 369 "zgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M much larger than N, JOBZ='A') */

#line 373 "zgesdd.f"
		    wrkbl = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 374 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *n + *m * ilaenv_(&c__1, "ZUNGQR", 
			    " ", m, m, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 374 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 376 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + (*n << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 376 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 378 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 378 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 380 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 380 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 382 "zgesdd.f"
		    maxwrk = *n * *n + wrkbl;
#line 383 "zgesdd.f"
		    minwrk = *n * *n + (*n << 1) + *m;
#line 384 "zgesdd.f"
		}
#line 385 "zgesdd.f"
	    } else if (*m >= mnthr2) {

/*              Path 5 (M much larger than N, but not as much as MNTHR1) */

#line 389 "zgesdd.f"
		maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 391 "zgesdd.f"
		minwrk = (*n << 1) + *m;
#line 392 "zgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 393 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "P", n, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 393 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 395 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "Q", m, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 395 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 397 "zgesdd.f"
		    maxwrk += *m * *n;
#line 398 "zgesdd.f"
		    minwrk += *n * *n;
#line 399 "zgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 400 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "P", n, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 400 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 402 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "Q", m, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 402 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 404 "zgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 405 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "P", n, n, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 405 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 407 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 407 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 409 "zgesdd.f"
		}
#line 410 "zgesdd.f"
	    } else {

/*              Path 6 (M at least N, but not much larger) */

#line 414 "zgesdd.f"
		maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 416 "zgesdd.f"
		minwrk = (*n << 1) + *m;
#line 417 "zgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 418 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 418 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 420 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", m, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 420 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 422 "zgesdd.f"
		    maxwrk += *m * *n;
#line 423 "zgesdd.f"
		    minwrk += *n * *n;
#line 424 "zgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 425 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 425 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 427 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", m, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 427 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 429 "zgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 430 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "PRC", n, n, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 430 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 432 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 432 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 434 "zgesdd.f"
		}
#line 435 "zgesdd.f"
	    }
#line 436 "zgesdd.f"
	} else {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD is BDSPAC */
/*           for computing singular values and singular vectors; BDSPAN */
/*           for computing singular values only. */
/*           BDSPAC = 5*M*M + 7*M */
/*           BDSPAN = MAX(7*M+4, 3*M+2+SMLSIZ*(SMLSIZ+8)) */

#line 445 "zgesdd.f"
	    if (*n >= mnthr1) {
#line 446 "zgesdd.f"
		if (wntqn) {

/*                 Path 1t (N much larger than M, JOBZ='N') */

#line 450 "zgesdd.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 452 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 452 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 454 "zgesdd.f"
		    minwrk = *m * 3;
#line 455 "zgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N much larger than M, JOBZ='O') */

#line 459 "zgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 460 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "ZUNGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 460 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 462 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 462 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 464 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 464 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 466 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 466 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 468 "zgesdd.f"
		    maxwrk = *m * *n + *m * *m + wrkbl;
#line 469 "zgesdd.f"
		    minwrk = (*m << 1) * *m + *m * 3;
#line 470 "zgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N much larger than M, JOBZ='S') */

#line 474 "zgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 475 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *m * ilaenv_(&c__1, "ZUNGLQ", 
			    " ", m, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 475 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 477 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 477 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 479 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 479 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 481 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 481 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 483 "zgesdd.f"
		    maxwrk = *m * *m + wrkbl;
#line 484 "zgesdd.f"
		    minwrk = *m * *m + *m * 3;
#line 485 "zgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N much larger than M, JOBZ='A') */

#line 489 "zgesdd.f"
		    wrkbl = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 490 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *m + *n * ilaenv_(&c__1, "ZUNGLQ", 
			    " ", n, n, m, &c_n1, (ftnlen)6, (ftnlen)1);
#line 490 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 492 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + (*m << 1) * ilaenv_(&
			    c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 492 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 494 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 494 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 496 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", m, m, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 496 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 498 "zgesdd.f"
		    maxwrk = *m * *m + wrkbl;
#line 499 "zgesdd.f"
		    minwrk = *m * *m + (*m << 1) + *n;
#line 500 "zgesdd.f"
		}
#line 501 "zgesdd.f"
	    } else if (*n >= mnthr2) {

/*              Path 5t (N much larger than M, but not as much as MNTHR1) */

#line 505 "zgesdd.f"
		maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 507 "zgesdd.f"
		minwrk = (*m << 1) + *n;
#line 508 "zgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 509 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "P", m, n, m, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 509 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 511 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 511 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 513 "zgesdd.f"
		    maxwrk += *m * *n;
#line 514 "zgesdd.f"
		    minwrk += *m * *m;
#line 515 "zgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 516 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "P", m, n, m, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 516 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 518 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 518 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 520 "zgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 521 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "P", n, n, m, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 521 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 523 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "Q", m, m, n, &c_n1, (ftnlen)6, (ftnlen)
			    1);
#line 523 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 525 "zgesdd.f"
		}
#line 526 "zgesdd.f"
	    } else {

/*              Path 6t (N greater than M, but not much larger) */

#line 530 "zgesdd.f"
		maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 532 "zgesdd.f"
		minwrk = (*m << 1) + *n;
#line 533 "zgesdd.f"
		if (wntqo) {
/* Computing MAX */
#line 534 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "PRC", m, n, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 534 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 536 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 536 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 538 "zgesdd.f"
		    maxwrk += *m * *n;
#line 539 "zgesdd.f"
		    minwrk += *m * *m;
#line 540 "zgesdd.f"
		} else if (wntqs) {
/* Computing MAX */
#line 541 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "PRC", m, n, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 541 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 543 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 543 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 545 "zgesdd.f"
		} else if (wntqa) {
/* Computing MAX */
#line 546 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *n * ilaenv_(&c__1, 
			    "ZUNGBR", "PRC", n, n, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 546 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 548 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNGBR", "QLN", m, m, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 548 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 550 "zgesdd.f"
		}
#line 551 "zgesdd.f"
	    }
#line 552 "zgesdd.f"
	}
#line 553 "zgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 554 "zgesdd.f"
    }
#line 555 "zgesdd.f"
    if (*info == 0) {
#line 556 "zgesdd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 557 "zgesdd.f"
	if (*lwork < minwrk && *lwork != -1) {
#line 557 "zgesdd.f"
	    *info = -13;
#line 557 "zgesdd.f"
	}
#line 559 "zgesdd.f"
    }

/*     Quick returns */

#line 563 "zgesdd.f"
    if (*info != 0) {
#line 564 "zgesdd.f"
	i__1 = -(*info);
#line 564 "zgesdd.f"
	xerbla_("ZGESDD", &i__1, (ftnlen)6);
#line 565 "zgesdd.f"
	return 0;
#line 566 "zgesdd.f"
    }
#line 567 "zgesdd.f"
    if (*lwork == -1) {
#line 567 "zgesdd.f"
	return 0;
#line 567 "zgesdd.f"
    }
#line 569 "zgesdd.f"
    if (*m == 0 || *n == 0) {
#line 570 "zgesdd.f"
	return 0;
#line 571 "zgesdd.f"
    }

/*     Get machine constants */

#line 575 "zgesdd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 576 "zgesdd.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 577 "zgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 581 "zgesdd.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 582 "zgesdd.f"
    iscl = 0;
#line 583 "zgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 584 "zgesdd.f"
	iscl = 1;
#line 585 "zgesdd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 586 "zgesdd.f"
    } else if (anrm > bignum) {
#line 587 "zgesdd.f"
	iscl = 1;
#line 588 "zgesdd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 589 "zgesdd.f"
    }

#line 591 "zgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 597 "zgesdd.f"
	if (*m >= mnthr1) {

#line 599 "zgesdd.f"
	    if (wntqn) {

/*              Path 1 (M much larger than N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 604 "zgesdd.f"
		itau = 1;
#line 605 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: need 0) */

#line 611 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 611 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 616 "zgesdd.f"
		i__1 = *n - 1;
#line 616 "zgesdd.f"
		i__2 = *n - 1;
#line 616 "zgesdd.f"
		zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 618 "zgesdd.f"
		ie = 1;
#line 619 "zgesdd.f"
		itauq = 1;
#line 620 "zgesdd.f"
		itaup = itauq + *n;
#line 621 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 627 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 627 "zgesdd.f"
		zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 630 "zgesdd.f"
		nrwork = ie + *n;

/*              Perform bidiagonal SVD, compute singular values only */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAN) */

#line 636 "zgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 639 "zgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M much larger than N, JOBZ='O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 645 "zgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 649 "zgesdd.f"
		ldwrku = *n;
#line 650 "zgesdd.f"
		ir = iu + ldwrku * *n;
#line 651 "zgesdd.f"
		if (*lwork >= *m * *n + *n * *n + *n * 3) {

/*                 WORK(IR) is M by N */

#line 655 "zgesdd.f"
		    ldwrkr = *m;
#line 656 "zgesdd.f"
		} else {
#line 657 "zgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 658 "zgesdd.f"
		}
#line 659 "zgesdd.f"
		itau = ir + ldwrkr * *n;
#line 660 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need N*N+2*N, prefer M*N+N+N*NB) */
/*              (RWorkspace: 0) */

#line 666 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 666 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK( IR ), zeroing out below it */

#line 671 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 672 "zgesdd.f"
		i__1 = *n - 1;
#line 672 "zgesdd.f"
		i__2 = *n - 1;
#line 672 "zgesdd.f"
		zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 679 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 679 "zgesdd.f"
		zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 681 "zgesdd.f"
		ie = 1;
#line 682 "zgesdd.f"
		itauq = itau;
#line 683 "zgesdd.f"
		itaup = itauq + *n;
#line 684 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 690 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 690 "zgesdd.f"
		zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of R in WORK(IRU) and computing right singular vectors */
/*              of R in WORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 700 "zgesdd.f"
		iru = ie + *n;
#line 701 "zgesdd.f"
		irvt = iru + *n * *n;
#line 702 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 703 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of R */
/*              (CWorkspace: need 2*N*N+3*N, prefer M*N+N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 712 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 714 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 714 "zgesdd.f"
		zunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by the right singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 723 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 724 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 724 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              (CWorkspace: need 2*N*N, prefer N*N+M*N) */
/*              (RWorkspace: 0) */

#line 733 "zgesdd.f"
		i__1 = *m;
#line 733 "zgesdd.f"
		i__2 = ldwrkr;
#line 733 "zgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 734 "zgesdd.f"
		    i__3 = *m - i__ + 1;
#line 734 "zgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 735 "zgesdd.f"
		    zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1], 
			    lda, &work[iu], &ldwrku, &c_b1, &work[ir], &
			    ldwrkr, (ftnlen)1, (ftnlen)1);
#line 738 "zgesdd.f"
		    zlacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 740 "zgesdd.f"
/* L10: */
#line 740 "zgesdd.f"
		}

#line 742 "zgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M much larger than N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 748 "zgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 752 "zgesdd.f"
		ldwrkr = *n;
#line 753 "zgesdd.f"
		itau = ir + ldwrkr * *n;
#line 754 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*              (RWorkspace: 0) */

#line 760 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 760 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 765 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 766 "zgesdd.f"
		i__2 = *n - 1;
#line 766 "zgesdd.f"
		i__1 = *n - 1;
#line 766 "zgesdd.f"
		zlaset_("L", &i__2, &i__1, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 773 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 773 "zgesdd.f"
		zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 775 "zgesdd.f"
		ie = 1;
#line 776 "zgesdd.f"
		itauq = itau;
#line 777 "zgesdd.f"
		itaup = itauq + *n;
#line 778 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 784 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 784 "zgesdd.f"
		zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 794 "zgesdd.f"
		iru = ie + *n;
#line 795 "zgesdd.f"
		irvt = iru + *n * *n;
#line 796 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 797 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 806 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 807 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 807 "zgesdd.f"
		zunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 816 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 817 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 817 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              (CWorkspace: need N*N) */
/*              (RWorkspace: 0) */

#line 826 "zgesdd.f"
		zlacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 827 "zgesdd.f"
		zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &work[ir],
			 &ldwrkr, &c_b1, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 830 "zgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M much larger than N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 836 "zgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 840 "zgesdd.f"
		ldwrku = *n;
#line 841 "zgesdd.f"
		itau = iu + ldwrku * *n;
#line 842 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 848 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 848 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 850 "zgesdd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              (CWorkspace: need N+M, prefer N+M*NB) */
/*              (RWorkspace: 0) */

#line 856 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 856 "zgesdd.f"
		zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out below it */

#line 861 "zgesdd.f"
		i__2 = *n - 1;
#line 861 "zgesdd.f"
		i__1 = *n - 1;
#line 861 "zgesdd.f"
		zlaset_("L", &i__2, &i__1, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 863 "zgesdd.f"
		ie = 1;
#line 864 "zgesdd.f"
		itauq = itau;
#line 865 "zgesdd.f"
		itaup = itauq + *n;
#line 866 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 872 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 872 "zgesdd.f"
		zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 875 "zgesdd.f"
		iru = ie + *n;
#line 876 "zgesdd.f"
		irvt = iru + *n * *n;
#line 877 "zgesdd.f"
		nrwork = irvt + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 885 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by left singular vectors of R */
/*              (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 894 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 896 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 896 "zgesdd.f"
		zunmbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 905 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 906 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 906 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              (CWorkspace: need N*N) */
/*              (RWorkspace: 0) */

#line 915 "zgesdd.f"
		zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &work[iu],
			 &ldwrku, &c_b1, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 920 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 922 "zgesdd.f"
	    }

#line 924 "zgesdd.f"
	} else if (*m >= mnthr2) {

/*           MNTHR2 <= M < MNTHR1 */

/*           Path 5 (M much larger than N, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           ZUNGBR and matrix multiplication to compute singular vectors */

#line 932 "zgesdd.f"
	    ie = 1;
#line 933 "zgesdd.f"
	    nrwork = ie + *n;
#line 934 "zgesdd.f"
	    itauq = 1;
#line 935 "zgesdd.f"
	    itaup = itauq + *n;
#line 936 "zgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 942 "zgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 942 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 945 "zgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 951 "zgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 953 "zgesdd.f"
	    } else if (wntqo) {
#line 954 "zgesdd.f"
		iu = nwork;
#line 955 "zgesdd.f"
		iru = nrwork;
#line 956 "zgesdd.f"
		irvt = iru + *n * *n;
#line 957 "zgesdd.f"
		nrwork = irvt + *n * *n;

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 963 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 964 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 964 "zgesdd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: 0) */

#line 971 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 971 "zgesdd.f"
		zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

#line 974 "zgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 978 "zgesdd.f"
		    ldwrku = *m;
#line 979 "zgesdd.f"
		} else {

/*                 WORK(IU) is LDWRKU by N */

#line 983 "zgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 984 "zgesdd.f"
		}
#line 985 "zgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 993 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in WORK(IU), copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1002 "zgesdd.f"
		zlarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &work[iu]
			, &ldwrku, &rwork[nrwork]);
#line 1004 "zgesdd.f"
		zlacpy_("F", n, n, &work[iu], &ldwrku, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in A by real matrix RWORK(IRU), storing the */
/*              result in WORK(IU), copying to A */
/*              (CWorkspace: need N*N, prefer M*N) */
/*              (Rworkspace: need 3*N*N, prefer N*N+2*M*N) */

#line 1011 "zgesdd.f"
		nrwork = irvt;
#line 1012 "zgesdd.f"
		i__2 = *m;
#line 1012 "zgesdd.f"
		i__1 = ldwrku;
#line 1012 "zgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1013 "zgesdd.f"
		    i__3 = *m - i__ + 1;
#line 1013 "zgesdd.f"
		    chunk = min(i__3,ldwrku);
#line 1014 "zgesdd.f"
		    zlacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru], n, 
			    &work[iu], &ldwrku, &rwork[nrwork]);
#line 1016 "zgesdd.f"
		    zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 1018 "zgesdd.f"
/* L20: */
#line 1018 "zgesdd.f"
		}

#line 1020 "zgesdd.f"
	    } else if (wntqs) {

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1026 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1027 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1027 "zgesdd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1034 "zgesdd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1035 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1035 "zgesdd.f"
		zungbr_("Q", m, n, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1044 "zgesdd.f"
		iru = nrwork;
#line 1045 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1046 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1047 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1056 "zgesdd.f"
		zlarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1058 "zgesdd.f"
		zlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: need 0) */
/*              (Rworkspace: need N*N+2*M*N) */

#line 1065 "zgesdd.f"
		nrwork = irvt;
#line 1066 "zgesdd.f"
		zlacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1068 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1069 "zgesdd.f"
	    } else {

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1075 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1076 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1076 "zgesdd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: 0) */

#line 1083 "zgesdd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1084 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1084 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1093 "zgesdd.f"
		iru = nrwork;
#line 1094 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1095 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1096 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1105 "zgesdd.f"
		zlarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1107 "zgesdd.f"
		zlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: 0) */
/*              (Rworkspace: need 3*N*N) */

#line 1114 "zgesdd.f"
		nrwork = irvt;
#line 1115 "zgesdd.f"
		zlacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1117 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1118 "zgesdd.f"
	    }

#line 1120 "zgesdd.f"
	} else {

/*           M .LT. MNTHR2 */

/*           Path 6 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */
/*           Use ZUNMBR to compute singular vectors */

#line 1128 "zgesdd.f"
	    ie = 1;
#line 1129 "zgesdd.f"
	    nrwork = ie + *n;
#line 1130 "zgesdd.f"
	    itauq = 1;
#line 1131 "zgesdd.f"
	    itaup = itauq + *n;
#line 1132 "zgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 1138 "zgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1138 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);
#line 1141 "zgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 1147 "zgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1149 "zgesdd.f"
	    } else if (wntqo) {
#line 1150 "zgesdd.f"
		iu = nwork;
#line 1151 "zgesdd.f"
		iru = nrwork;
#line 1152 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1153 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1154 "zgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 1158 "zgesdd.f"
		    ldwrku = *m;
#line 1159 "zgesdd.f"
		} else {

/*                 WORK( IU ) is LDWRKU by N */

#line 1163 "zgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 1164 "zgesdd.f"
		}
#line 1165 "zgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1173 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (Cworkspace: need 2*N, prefer N+N*NB) */
/*              (Rworkspace: need 0) */

#line 1182 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1183 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1183 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1187 "zgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by left singular vectors of A, copying */
/*              to A */
/*              (Cworkspace: need M*N+2*N, prefer M*N+N+N*NB) */
/*              (Rworkspace: need 0) */

#line 1195 "zgesdd.f"
		    zlaset_("F", m, n, &c_b1, &c_b1, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1197 "zgesdd.f"
		    zlacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1199 "zgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1199 "zgesdd.f"
		    zunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1202 "zgesdd.f"
		    zlacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 1203 "zgesdd.f"
		} else {

/*                 Generate Q in A */
/*                 (Cworkspace: need 2*N, prefer N+N*NB) */
/*                 (Rworkspace: need 0) */

#line 1209 "zgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1209 "zgesdd.f"
		    zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__1, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 (CWorkspace: need N*N, prefer M*N) */
/*                 (Rworkspace: need 3*N*N, prefer N*N+2*M*N) */

#line 1217 "zgesdd.f"
		    nrwork = irvt;
#line 1218 "zgesdd.f"
		    i__1 = *m;
#line 1218 "zgesdd.f"
		    i__2 = ldwrku;
#line 1218 "zgesdd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ +=
			     i__2) {
/* Computing MIN */
#line 1219 "zgesdd.f"
			i__3 = *m - i__ + 1;
#line 1219 "zgesdd.f"
			chunk = min(i__3,ldwrku);
#line 1220 "zgesdd.f"
			zlacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru],
				 n, &work[iu], &ldwrku, &rwork[nrwork]);
#line 1223 "zgesdd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 1225 "zgesdd.f"
/* L30: */
#line 1225 "zgesdd.f"
		    }
#line 1226 "zgesdd.f"
		}

#line 1228 "zgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1236 "zgesdd.f"
		iru = nrwork;
#line 1237 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1238 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1239 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1248 "zgesdd.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1249 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1250 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1250 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1259 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1260 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1260 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1263 "zgesdd.f"
	    } else {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1271 "zgesdd.f"
		iru = nrwork;
#line 1272 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1273 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1274 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 1280 "zgesdd.f"
		zlaset_("F", m, m, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1281 "zgesdd.f"
		if (*m > *n) {
#line 1282 "zgesdd.f"
		    i__2 = *m - *n;
#line 1282 "zgesdd.f"
		    i__1 = *m - *n;
#line 1282 "zgesdd.f"
		    zlaset_("F", &i__2, &i__1, &c_b1, &c_b2, &u[*n + 1 + (*n 
			    + 1) * u_dim1], ldu, (ftnlen)1);
#line 1284 "zgesdd.f"
		}

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*              (RWorkspace: 0) */

#line 1291 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1292 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1292 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1301 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1302 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1302 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1305 "zgesdd.f"
	    }

#line 1307 "zgesdd.f"
	}

#line 1309 "zgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 1315 "zgesdd.f"
	if (*n >= mnthr1) {

#line 1317 "zgesdd.f"
	    if (wntqn) {

/*              Path 1t (N much larger than M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 1322 "zgesdd.f"
		itau = 1;
#line 1323 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1329 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1329 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 1334 "zgesdd.f"
		i__2 = *m - 1;
#line 1334 "zgesdd.f"
		i__1 = *m - 1;
#line 1334 "zgesdd.f"
		zlaset_("U", &i__2, &i__1, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1336 "zgesdd.f"
		ie = 1;
#line 1337 "zgesdd.f"
		itauq = 1;
#line 1338 "zgesdd.f"
		itaup = itauq + *m;
#line 1339 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1345 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1345 "zgesdd.f"
		zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 1348 "zgesdd.f"
		nrwork = ie + *m;

/*              Perform bidiagonal SVD, compute singular values only */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAN) */

#line 1354 "zgesdd.f"
		dbdsdc_("U", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 1357 "zgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N much larger than M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1363 "zgesdd.f"
		ivt = 1;
#line 1364 "zgesdd.f"
		ldwkvt = *m;

/*              WORK(IVT) is M by M */

#line 1368 "zgesdd.f"
		il = ivt + ldwkvt * *m;
#line 1369 "zgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3) {

/*                 WORK(IL) M by N */

#line 1373 "zgesdd.f"
		    ldwrkl = *m;
#line 1374 "zgesdd.f"
		    chunk = *n;
#line 1375 "zgesdd.f"
		} else {

/*                 WORK(IL) is M by CHUNK */

#line 1379 "zgesdd.f"
		    ldwrkl = *m;
#line 1380 "zgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1381 "zgesdd.f"
		}
#line 1382 "zgesdd.f"
		itau = il + ldwrkl * chunk;
#line 1383 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1389 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1389 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1394 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1395 "zgesdd.f"
		i__2 = *m - 1;
#line 1395 "zgesdd.f"
		i__1 = *m - 1;
#line 1395 "zgesdd.f"
		zlaset_("U", &i__2, &i__1, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*              (RWorkspace: 0) */

#line 1402 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1402 "zgesdd.f"
		zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1404 "zgesdd.f"
		ie = 1;
#line 1405 "zgesdd.f"
		itauq = itau;
#line 1406 "zgesdd.f"
		itaup = itauq + *m;
#line 1407 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1413 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1413 "zgesdd.f"
		zgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1423 "zgesdd.f"
		iru = ie + *m;
#line 1424 "zgesdd.f"
		irvt = iru + *m * *m;
#line 1425 "zgesdd.f"
		nrwork = irvt + *m * *m;
#line 1426 "zgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of L */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1435 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1436 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1436 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by the right singular vectors of L */
/*              (CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 1445 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1447 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1447 "zgesdd.f"
		zunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              (CWorkspace: need 2*M*M, prefer M*M+M*N)) */
/*              (RWorkspace: 0) */

#line 1456 "zgesdd.f"
		i__2 = *n;
#line 1456 "zgesdd.f"
		i__1 = chunk;
#line 1456 "zgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1457 "zgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1457 "zgesdd.f"
		    blk = min(i__3,chunk);
#line 1458 "zgesdd.f"
		    zgemm_("N", "N", m, &blk, m, &c_b2, &work[ivt], m, &a[i__ 
			    * a_dim1 + 1], lda, &c_b1, &work[il], &ldwrkl, (
			    ftnlen)1, (ftnlen)1);
#line 1461 "zgesdd.f"
		    zlacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1463 "zgesdd.f"
/* L40: */
#line 1463 "zgesdd.f"
		}

#line 1465 "zgesdd.f"
	    } else if (wntqs) {

/*             Path 3t (N much larger than M, JOBZ='S') */
/*             M right singular vectors to be computed in VT and */
/*             M left singular vectors to be computed in U */

#line 1471 "zgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1475 "zgesdd.f"
		ldwrkl = *m;
#line 1476 "zgesdd.f"
		itau = il + ldwrkl * *m;
#line 1477 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1483 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1483 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1488 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1489 "zgesdd.f"
		i__1 = *m - 1;
#line 1489 "zgesdd.f"
		i__2 = *m - 1;
#line 1489 "zgesdd.f"
		zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*              (RWorkspace: 0) */

#line 1496 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1496 "zgesdd.f"
		zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1498 "zgesdd.f"
		ie = 1;
#line 1499 "zgesdd.f"
		itauq = itau;
#line 1500 "zgesdd.f"
		itaup = itauq + *m;
#line 1501 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1507 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1507 "zgesdd.f"
		zgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1517 "zgesdd.f"
		iru = ie + *m;
#line 1518 "zgesdd.f"
		irvt = iru + *m * *m;
#line 1519 "zgesdd.f"
		nrwork = irvt + *m * *m;
#line 1520 "zgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1529 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1530 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1530 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by left singular vectors of L */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1539 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1540 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1540 "zgesdd.f"
		zunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy VT to WORK(IL), multiply right singular vectors of L */
/*              in WORK(IL) by Q in A, storing result in VT */
/*              (CWorkspace: need M*M) */
/*              (RWorkspace: 0) */

#line 1549 "zgesdd.f"
		zlacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1550 "zgesdd.f"
		zgemm_("N", "N", m, n, m, &c_b2, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b1, &vt[vt_offset], ldvt, (ftnlen)
			1, (ftnlen)1);

#line 1553 "zgesdd.f"
	    } else if (wntqa) {

/*              Path 9t (N much larger than M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1559 "zgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1563 "zgesdd.f"
		ldwkvt = *m;
#line 1564 "zgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1565 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 1571 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1571 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
#line 1573 "zgesdd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              (CWorkspace: need M+N, prefer M+N*NB) */
/*              (RWorkspace: 0) */

#line 1579 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1579 "zgesdd.f"
		zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

/*              Produce L in A, zeroing out above it */

#line 1584 "zgesdd.f"
		i__1 = *m - 1;
#line 1584 "zgesdd.f"
		i__2 = *m - 1;
#line 1584 "zgesdd.f"
		zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1586 "zgesdd.f"
		ie = 1;
#line 1587 "zgesdd.f"
		itauq = itau;
#line 1588 "zgesdd.f"
		itaup = itauq + *m;
#line 1589 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 1595 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1595 "zgesdd.f"
		zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1605 "zgesdd.f"
		iru = ie + *m;
#line 1606 "zgesdd.f"
		irvt = iru + *m * *m;
#line 1607 "zgesdd.f"
		nrwork = irvt + *m * *m;
#line 1608 "zgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1617 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1618 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1618 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by right singular vectors of L */
/*              (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 1627 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1629 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1629 "zgesdd.f"
		zunmbr_("P", "R", "C", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              (CWorkspace: need M*M) */
/*              (RWorkspace: 0) */

#line 1638 "zgesdd.f"
		zgemm_("N", "N", m, n, m, &c_b2, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b1, &a[a_offset], lda, (ftnlen)1,
			 (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1643 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1645 "zgesdd.f"
	    }

#line 1647 "zgesdd.f"
	} else if (*n >= mnthr2) {

/*           MNTHR2 <= N < MNTHR1 */

/*           Path 5t (N much larger than M, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           ZUNGBR and matrix multiplication to compute singular vectors */


#line 1656 "zgesdd.f"
	    ie = 1;
#line 1657 "zgesdd.f"
	    nrwork = ie + *m;
#line 1658 "zgesdd.f"
	    itauq = 1;
#line 1659 "zgesdd.f"
	    itaup = itauq + *m;
#line 1660 "zgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 1666 "zgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1666 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);

#line 1670 "zgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 1676 "zgesdd.f"
		dbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1678 "zgesdd.f"
	    } else if (wntqo) {
#line 1679 "zgesdd.f"
		irvt = nrwork;
#line 1680 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1681 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1682 "zgesdd.f"
		ivt = nwork;

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1688 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1689 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1689 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Generate P**H in A */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1696 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1696 "zgesdd.f"
		zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

#line 1699 "zgesdd.f"
		ldwkvt = *m;
#line 1700 "zgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 1704 "zgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1705 "zgesdd.f"
		    chunk = *n;
#line 1706 "zgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 1710 "zgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 1711 "zgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 1712 "zgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1720 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRVT) */
/*              storing the result in WORK(IVT), copying to U */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need 2*M*M) */

#line 1729 "zgesdd.f"
		zlacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &work[ivt], &
			ldwkvt, &rwork[nrwork]);
#line 1731 "zgesdd.f"
		zlacpy_("F", m, m, &work[ivt], &ldwkvt, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply RWORK(IRVT) by P**H in A, storing the */
/*              result in WORK(IVT), copying to A */
/*              (CWorkspace: need M*M, prefer M*N) */
/*              (Rworkspace: need 2*M*M, prefer 2*M*N) */

#line 1738 "zgesdd.f"
		nrwork = iru;
#line 1739 "zgesdd.f"
		i__1 = *n;
#line 1739 "zgesdd.f"
		i__2 = chunk;
#line 1739 "zgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1740 "zgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1740 "zgesdd.f"
		    blk = min(i__3,chunk);
#line 1741 "zgesdd.f"
		    zlarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1], 
			    lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 1743 "zgesdd.f"
		    zlacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
			    a_dim1 + 1], lda, (ftnlen)1);
#line 1745 "zgesdd.f"
/* L50: */
#line 1745 "zgesdd.f"
		}
#line 1746 "zgesdd.f"
	    } else if (wntqs) {

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1752 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1753 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1753 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1760 "zgesdd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1761 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1761 "zgesdd.f"
		zungbr_("P", m, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1770 "zgesdd.f"
		irvt = nrwork;
#line 1771 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1772 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1773 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: need 0) */
/*              (Rworkspace: need 3*M*M) */

#line 1782 "zgesdd.f"
		zlacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1784 "zgesdd.f"
		zlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need M*M+2*M*N) */

#line 1791 "zgesdd.f"
		nrwork = iru;
#line 1792 "zgesdd.f"
		zlarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1794 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1795 "zgesdd.f"
	    } else {

/*              Copy A to U, generate Q */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1801 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1802 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1802 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: 0) */

#line 1809 "zgesdd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1810 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1810 "zgesdd.f"
		zungbr_("P", n, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1819 "zgesdd.f"
		irvt = nrwork;
#line 1820 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1821 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1822 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              (CWorkspace: need 0) */
/*              (Rworkspace: need 3*M*M) */

#line 1831 "zgesdd.f"
		zlacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1833 "zgesdd.f"
		zlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              (Cworkspace: need 0) */
/*              (Rworkspace: need M*M+2*M*N) */

#line 1840 "zgesdd.f"
		zlarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1842 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1843 "zgesdd.f"
	    }

#line 1845 "zgesdd.f"
	} else {

/*           N .LT. MNTHR2 */

/*           Path 6t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           Use ZUNMBR to compute singular vectors */

#line 1853 "zgesdd.f"
	    ie = 1;
#line 1854 "zgesdd.f"
	    nrwork = ie + *m;
#line 1855 "zgesdd.f"
	    itauq = 1;
#line 1856 "zgesdd.f"
	    itaup = itauq + *m;
#line 1857 "zgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 1863 "zgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1863 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 1866 "zgesdd.f"
	    if (wntqn) {

/*              Compute singular values only */
/*              (Cworkspace: 0) */
/*              (Rworkspace: need BDSPAN) */

#line 1872 "zgesdd.f"
		dbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1874 "zgesdd.f"
	    } else if (wntqo) {
#line 1875 "zgesdd.f"
		ldwkvt = *m;
#line 1876 "zgesdd.f"
		ivt = nwork;
#line 1877 "zgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 1881 "zgesdd.f"
		    zlaset_("F", m, n, &c_b1, &c_b1, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 1883 "zgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1884 "zgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 1888 "zgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 1889 "zgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 1890 "zgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1898 "zgesdd.f"
		irvt = nrwork;
#line 1899 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1900 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1901 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (Cworkspace: need 2*M, prefer M+M*NB) */
/*              (Rworkspace: need 0) */

#line 1910 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1911 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1911 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1915 "zgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by right singular vectors of A, */
/*              copying to A */
/*              (Cworkspace: need M*N+2*M, prefer M*N+M+M*NB) */
/*              (Rworkspace: need 0) */

#line 1923 "zgesdd.f"
		    zlacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 1925 "zgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1925 "zgesdd.f"
		    zunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1928 "zgesdd.f"
		    zlacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 1929 "zgesdd.f"
		} else {

/*                 Generate P**H in A */
/*                 (Cworkspace: need 2*M, prefer M+M*NB) */
/*                 (Rworkspace: need 0) */

#line 1935 "zgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 1935 "zgesdd.f"
		    zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 (CWorkspace: need M*M, prefer M*N) */
/*                 (Rworkspace: need 3*M*M, prefer M*M+2*M*N) */

#line 1943 "zgesdd.f"
		    nrwork = iru;
#line 1944 "zgesdd.f"
		    i__2 = *n;
#line 1944 "zgesdd.f"
		    i__1 = chunk;
#line 1944 "zgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 1945 "zgesdd.f"
			i__3 = *n - i__ + 1;
#line 1945 "zgesdd.f"
			blk = min(i__3,chunk);
#line 1946 "zgesdd.f"
			zlarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1]
				, lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 1949 "zgesdd.f"
			zlacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 1951 "zgesdd.f"
/* L60: */
#line 1951 "zgesdd.f"
		    }
#line 1952 "zgesdd.f"
		}
#line 1953 "zgesdd.f"
	    } else if (wntqs) {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1961 "zgesdd.f"
		irvt = nrwork;
#line 1962 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1963 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1964 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: M*M) */

#line 1973 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1974 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1974 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: M*M) */

#line 1983 "zgesdd.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1984 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1985 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1985 "zgesdd.f"
		zunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1988 "zgesdd.f"
	    } else {

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              (CWorkspace: need 0) */
/*              (RWorkspace: need BDSPAC) */

#line 1996 "zgesdd.f"
		irvt = nrwork;
#line 1997 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1998 "zgesdd.f"
		nrwork = iru + *m * *m;

#line 2000 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: M*M) */

#line 2009 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2010 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2010 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Set all of VT to identity matrix */

#line 2016 "zgesdd.f"
		zlaset_("F", n, n, &c_b1, &c_b2, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*              (RWorkspace: M*M) */

#line 2023 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2024 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2024 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2027 "zgesdd.f"
	    }

#line 2029 "zgesdd.f"
	}

#line 2031 "zgesdd.f"
    }

/*     Undo scaling if necessary */

#line 2035 "zgesdd.f"
    if (iscl == 1) {
#line 2036 "zgesdd.f"
	if (anrm > bignum) {
#line 2036 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2036 "zgesdd.f"
	}
#line 2039 "zgesdd.f"
	if (*info != 0 && anrm > bignum) {
#line 2039 "zgesdd.f"
	    i__1 = minmn - 1;
#line 2039 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2039 "zgesdd.f"
	}
#line 2042 "zgesdd.f"
	if (anrm < smlnum) {
#line 2042 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2042 "zgesdd.f"
	}
#line 2045 "zgesdd.f"
	if (*info != 0 && anrm < smlnum) {
#line 2045 "zgesdd.f"
	    i__1 = minmn - 1;
#line 2045 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2045 "zgesdd.f"
	}
#line 2048 "zgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 2052 "zgesdd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 2054 "zgesdd.f"
    return 0;

/*     End of ZGESDD */

} /* zgesdd_ */

