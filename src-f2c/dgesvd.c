#line 1 "dgesvd.f"
/* dgesvd.f -- translated by f2c (version 20100827).
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

#line 1 "dgesvd.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c_n1 = -1;
static doublereal c_b57 = 0.;
static integer c__1 = 1;
static doublereal c_b79 = 1.;

/* > \brief <b> DGESVD computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGESVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU, JOBVT */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGESVD computes the singular value decomposition (SVD) of a real */
/* > M-by-N matrix A, optionally computing the left and/or right singular */
/* > vectors. The SVD is written */
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
/* > Note that the routine returns V**T, not V. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          Specifies options for computing all or part of the matrix U: */
/* >          = 'A':  all M columns of U are returned in array U: */
/* >          = 'S':  the first min(m,n) columns of U (the left singular */
/* >                  vectors) are returned in the array U; */
/* >          = 'O':  the first min(m,n) columns of U (the left singular */
/* >                  vectors) are overwritten on the array A; */
/* >          = 'N':  no columns of U (no left singular vectors) are */
/* >                  computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVT */
/* > \verbatim */
/* >          JOBVT is CHARACTER*1 */
/* >          Specifies options for computing all or part of the matrix */
/* >          V**T: */
/* >          = 'A':  all N rows of V**T are returned in the array VT; */
/* >          = 'S':  the first min(m,n) rows of V**T (the right singular */
/* >                  vectors) are returned in the array VT; */
/* >          = 'O':  the first min(m,n) rows of V**T (the right singular */
/* >                  vectors) are overwritten on the array A; */
/* >          = 'N':  no rows of V**T (no right singular vectors) are */
/* >                  computed. */
/* > */
/* >          JOBVT and JOBU cannot both be 'O'. */
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
/* >          if JOBU = 'O',  A is overwritten with the first min(m,n) */
/* >                          columns of U (the left singular vectors, */
/* >                          stored columnwise); */
/* >          if JOBVT = 'O', A is overwritten with the first min(m,n) */
/* >                          rows of V**T (the right singular vectors, */
/* >                          stored rowwise); */
/* >          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A */
/* >                          are destroyed. */
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
/* >          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'. */
/* >          If JOBU = 'A', U contains the M-by-M orthogonal matrix U; */
/* >          if JOBU = 'S', U contains the first min(m,n) columns of U */
/* >          (the left singular vectors, stored columnwise); */
/* >          if JOBU = 'N' or 'O', U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U.  LDU >= 1; if */
/* >          JOBU = 'S' or 'A', LDU >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension (LDVT,N) */
/* >          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix */
/* >          V**T; */
/* >          if JOBVT = 'S', VT contains the first min(m,n) rows of */
/* >          V**T (the right singular vectors, stored rowwise); */
/* >          if JOBVT = 'N' or 'O', VT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >          The leading dimension of the array VT.  LDVT >= 1; if */
/* >          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK; */
/* >          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged */
/* >          superdiagonal elements of an upper bidiagonal matrix B */
/* >          whose diagonal is in S (not necessarily sorted). B */
/* >          satisfies A = U * B * VT, so it has the same singular values */
/* >          as A, and singular vectors related by U and VT. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code): */
/* >             - PATH 1  (M much larger than N, JOBU='N') */
/* >             - PATH 1t (N much larger than M, JOBVT='N') */
/* >          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) for the other paths */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if DBDSQR did not converge, INFO specifies how many */
/* >                superdiagonals of an intermediate bidiagonal form B */
/* >                did not converge to zero. See the description of WORK */
/* >                above for details. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup doubleGEsing */

/*  ===================================================================== */
/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *info, ftnlen jobu_len, ftnlen jobvt_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1[2], 
	    i__2, i__3, i__4;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ie, ir, iu, blk, ncu;
    static doublereal dum[1], eps;
    static integer nru, iscl;
    static doublereal anrm;
    static integer ierr, itau, ncvt, nrvt, lwork_dgebrd__, lwork_dgelqf__, 
	    lwork_dgeqrf__;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn, wrkbl, itaup, itauq, mnthr, iwork;
    static logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    extern /* Subroutine */ int dgebrd_(integer *, integer *, doublereal *, 
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
	    dbdsqr_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dorgbr_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), dorglq_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dorgqr_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *);
    static integer ldwrkr, minwrk, ldwrku, maxwrk;
    static doublereal smlnum;
    static logical lquery, wntuas, wntvas;
    static integer lwork_dorgbr_p__, lwork_dorgbr_q__, lwork_dorglq_m__, 
	    lwork_dorglq_n__, lwork_dorgqr_m__, lwork_dorgqr_n__;


/*  -- LAPACK driver routine (version 3.4.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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

#line 267 "dgesvd.f"
    /* Parameter adjustments */
#line 267 "dgesvd.f"
    a_dim1 = *lda;
#line 267 "dgesvd.f"
    a_offset = 1 + a_dim1;
#line 267 "dgesvd.f"
    a -= a_offset;
#line 267 "dgesvd.f"
    --s;
#line 267 "dgesvd.f"
    u_dim1 = *ldu;
#line 267 "dgesvd.f"
    u_offset = 1 + u_dim1;
#line 267 "dgesvd.f"
    u -= u_offset;
#line 267 "dgesvd.f"
    vt_dim1 = *ldvt;
#line 267 "dgesvd.f"
    vt_offset = 1 + vt_dim1;
#line 267 "dgesvd.f"
    vt -= vt_offset;
#line 267 "dgesvd.f"
    --work;
#line 267 "dgesvd.f"

#line 267 "dgesvd.f"
    /* Function Body */
#line 267 "dgesvd.f"
    *info = 0;
#line 268 "dgesvd.f"
    minmn = min(*m,*n);
#line 269 "dgesvd.f"
    wntua = lsame_(jobu, "A", (ftnlen)1, (ftnlen)1);
#line 270 "dgesvd.f"
    wntus = lsame_(jobu, "S", (ftnlen)1, (ftnlen)1);
#line 271 "dgesvd.f"
    wntuas = wntua || wntus;
#line 272 "dgesvd.f"
    wntuo = lsame_(jobu, "O", (ftnlen)1, (ftnlen)1);
#line 273 "dgesvd.f"
    wntun = lsame_(jobu, "N", (ftnlen)1, (ftnlen)1);
#line 274 "dgesvd.f"
    wntva = lsame_(jobvt, "A", (ftnlen)1, (ftnlen)1);
#line 275 "dgesvd.f"
    wntvs = lsame_(jobvt, "S", (ftnlen)1, (ftnlen)1);
#line 276 "dgesvd.f"
    wntvas = wntva || wntvs;
#line 277 "dgesvd.f"
    wntvo = lsame_(jobvt, "O", (ftnlen)1, (ftnlen)1);
#line 278 "dgesvd.f"
    wntvn = lsame_(jobvt, "N", (ftnlen)1, (ftnlen)1);
#line 279 "dgesvd.f"
    lquery = *lwork == -1;

#line 281 "dgesvd.f"
    if (! (wntua || wntus || wntuo || wntun)) {
#line 282 "dgesvd.f"
	*info = -1;
#line 283 "dgesvd.f"
    } else if (! (wntva || wntvs || wntvo || wntvn) || wntvo && wntuo) {
#line 285 "dgesvd.f"
	*info = -2;
#line 286 "dgesvd.f"
    } else if (*m < 0) {
#line 287 "dgesvd.f"
	*info = -3;
#line 288 "dgesvd.f"
    } else if (*n < 0) {
#line 289 "dgesvd.f"
	*info = -4;
#line 290 "dgesvd.f"
    } else if (*lda < max(1,*m)) {
#line 291 "dgesvd.f"
	*info = -6;
#line 292 "dgesvd.f"
    } else if (*ldu < 1 || wntuas && *ldu < *m) {
#line 293 "dgesvd.f"
	*info = -9;
#line 294 "dgesvd.f"
    } else if (*ldvt < 1 || wntva && *ldvt < *n || wntvs && *ldvt < minmn) {
#line 296 "dgesvd.f"
	*info = -11;
#line 297 "dgesvd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 306 "dgesvd.f"
    if (*info == 0) {
#line 307 "dgesvd.f"
	minwrk = 1;
#line 308 "dgesvd.f"
	maxwrk = 1;
#line 309 "dgesvd.f"
	if (*m >= *n && minmn > 0) {

/*           Compute space needed for DBDSQR */

/* Writing concatenation */
#line 313 "dgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 313 "dgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 313 "dgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 313 "dgesvd.f"
	    mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
#line 314 "dgesvd.f"
	    bdspac = *n * 5;
/*           Compute space needed for DGEQRF */
#line 316 "dgesvd.f"
	    dgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 317 "dgesvd.f"
	    lwork_dgeqrf__ = (integer) dum[0];
/*           Compute space needed for DORGQR */
#line 319 "dgesvd.f"
	    dorgqr_(m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 320 "dgesvd.f"
	    lwork_dorgqr_n__ = (integer) dum[0];
#line 321 "dgesvd.f"
	    dorgqr_(m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 322 "dgesvd.f"
	    lwork_dorgqr_m__ = (integer) dum[0];
/*           Compute space needed for DGEBRD */
#line 324 "dgesvd.f"
	    dgebrd_(n, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 326 "dgesvd.f"
	    lwork_dgebrd__ = (integer) dum[0];
/*           Compute space needed for DORGBR P */
#line 328 "dgesvd.f"
	    dorgbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 330 "dgesvd.f"
	    lwork_dorgbr_p__ = (integer) dum[0];
/*           Compute space needed for DORGBR Q */
#line 332 "dgesvd.f"
	    dorgbr_("Q", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 334 "dgesvd.f"
	    lwork_dorgbr_q__ = (integer) dum[0];

#line 336 "dgesvd.f"
	    if (*m >= mnthr) {
#line 337 "dgesvd.f"
		if (wntun) {

/*                 Path 1 (M much larger than N, JOBU='N') */

#line 341 "dgesvd.f"
		    maxwrk = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 342 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_dgebrd__;
#line 342 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 343 "dgesvd.f"
		    if (wntvo || wntvas) {
/* Computing MAX */
#line 343 "dgesvd.f"
			i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_p__;
#line 343 "dgesvd.f"
			maxwrk = max(i__2,i__3);
#line 343 "dgesvd.f"
		    }
#line 345 "dgesvd.f"
		    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 346 "dgesvd.f"
		    i__2 = *n << 2;
#line 346 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 347 "dgesvd.f"
		} else if (wntuo && wntvn) {

/*                 Path 2 (M much larger than N, JOBU='O', JOBVT='N') */

#line 351 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 352 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
#line 352 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 353 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 353 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 354 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 354 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 355 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 356 "dgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
#line 356 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 357 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 357 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 358 "dgesvd.f"
		} else if (wntuo && wntvas) {

/*                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or */
/*                 'A') */

#line 363 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 364 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
#line 364 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 365 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 365 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 366 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 366 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 367 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
#line 367 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 368 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 369 "dgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
#line 369 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 370 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 370 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 371 "dgesvd.f"
		} else if (wntus && wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */

#line 375 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 376 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
#line 376 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 377 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 377 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 378 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 378 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 379 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 380 "dgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 381 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 381 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 382 "dgesvd.f"
		} else if (wntus && wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */

#line 386 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 387 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
#line 387 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 388 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 388 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 389 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 389 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 390 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
#line 390 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 391 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 392 "dgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
#line 393 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 393 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 394 "dgesvd.f"
		} else if (wntus && wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or */
/*                 'A') */

#line 399 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 400 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
#line 400 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 401 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 401 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 402 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 402 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 403 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
#line 403 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 404 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 405 "dgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 406 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 406 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 407 "dgesvd.f"
		} else if (wntua && wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */

#line 411 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 412 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
#line 412 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 413 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 413 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 414 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 414 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 415 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 416 "dgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 417 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 417 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 418 "dgesvd.f"
		} else if (wntua && wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */

#line 422 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 423 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
#line 423 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 424 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 424 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 425 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 425 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 426 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
#line 426 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 427 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 428 "dgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
#line 429 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 429 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 430 "dgesvd.f"
		} else if (wntua && wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or */
/*                 'A') */

#line 435 "dgesvd.f"
		    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
#line 436 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
#line 436 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 437 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
#line 437 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 438 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 438 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 439 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
#line 439 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 440 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 441 "dgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 442 "dgesvd.f"
		    i__2 = *n * 3 + *m;
#line 442 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 443 "dgesvd.f"
		}
#line 444 "dgesvd.f"
	    } else {

/*              Path 10 (M at least N, but not much larger) */

#line 448 "dgesvd.f"
		dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 450 "dgesvd.f"
		lwork_dgebrd__ = (integer) dum[0];
#line 451 "dgesvd.f"
		maxwrk = *n * 3 + lwork_dgebrd__;
#line 452 "dgesvd.f"
		if (wntus || wntuo) {
#line 453 "dgesvd.f"
		    dorgbr_("Q", m, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 455 "dgesvd.f"
		    lwork_dorgbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 456 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 456 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 457 "dgesvd.f"
		}
#line 458 "dgesvd.f"
		if (wntua) {
#line 459 "dgesvd.f"
		    dorgbr_("Q", m, m, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 461 "dgesvd.f"
		    lwork_dorgbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 462 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_q__;
#line 462 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 463 "dgesvd.f"
		}
#line 464 "dgesvd.f"
		if (! wntvn) {
/* Computing MAX */
#line 465 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_p__;
#line 465 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 466 "dgesvd.f"
		}
#line 467 "dgesvd.f"
		maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 468 "dgesvd.f"
		i__2 = *n * 3 + *m;
#line 468 "dgesvd.f"
		minwrk = max(i__2,bdspac);
#line 469 "dgesvd.f"
	    }
#line 470 "dgesvd.f"
	} else if (minmn > 0) {

/*           Compute space needed for DBDSQR */

/* Writing concatenation */
#line 474 "dgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 474 "dgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 474 "dgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 474 "dgesvd.f"
	    mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
#line 475 "dgesvd.f"
	    bdspac = *m * 5;
/*           Compute space needed for DGELQF */
#line 477 "dgesvd.f"
	    dgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 478 "dgesvd.f"
	    lwork_dgelqf__ = (integer) dum[0];
/*           Compute space needed for DORGLQ */
#line 480 "dgesvd.f"
	    dorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
#line 481 "dgesvd.f"
	    lwork_dorglq_n__ = (integer) dum[0];
#line 482 "dgesvd.f"
	    dorglq_(m, n, m, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 483 "dgesvd.f"
	    lwork_dorglq_m__ = (integer) dum[0];
/*           Compute space needed for DGEBRD */
#line 485 "dgesvd.f"
	    dgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 487 "dgesvd.f"
	    lwork_dgebrd__ = (integer) dum[0];
/*            Compute space needed for DORGBR P */
#line 489 "dgesvd.f"
	    dorgbr_("P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 491 "dgesvd.f"
	    lwork_dorgbr_p__ = (integer) dum[0];
/*           Compute space needed for DORGBR Q */
#line 493 "dgesvd.f"
	    dorgbr_("Q", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 495 "dgesvd.f"
	    lwork_dorgbr_q__ = (integer) dum[0];
#line 496 "dgesvd.f"
	    if (*n >= mnthr) {
#line 497 "dgesvd.f"
		if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N') */

#line 501 "dgesvd.f"
		    maxwrk = *m + lwork_dgelqf__;
/* Computing MAX */
#line 502 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_dgebrd__;
#line 502 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 503 "dgesvd.f"
		    if (wntuo || wntuas) {
/* Computing MAX */
#line 503 "dgesvd.f"
			i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_q__;
#line 503 "dgesvd.f"
			maxwrk = max(i__2,i__3);
#line 503 "dgesvd.f"
		    }
#line 505 "dgesvd.f"
		    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 506 "dgesvd.f"
		    i__2 = *m << 2;
#line 506 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 507 "dgesvd.f"
		} else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O') */

#line 511 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 512 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
#line 512 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 513 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 513 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 514 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 514 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 515 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 516 "dgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
#line 516 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 517 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 517 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 518 "dgesvd.f"
		} else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='O') */

#line 523 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 524 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
#line 524 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 525 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 525 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 526 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 526 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 527 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
#line 527 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 528 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 529 "dgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
#line 529 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 530 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 530 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 531 "dgesvd.f"
		} else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */

#line 535 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 536 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
#line 536 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 537 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 537 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 538 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 538 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 539 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 540 "dgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 541 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 541 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 542 "dgesvd.f"
		} else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */

#line 546 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 547 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
#line 547 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 548 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 548 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 549 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 549 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 550 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
#line 550 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 551 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 552 "dgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
#line 553 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 553 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 554 "dgesvd.f"
		} else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='S') */

#line 559 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 560 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
#line 560 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 561 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 561 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 562 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 562 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 563 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
#line 563 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 564 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 565 "dgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 566 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 566 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 567 "dgesvd.f"
		} else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */

#line 571 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 572 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
#line 572 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 573 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 573 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 574 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 574 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 575 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 576 "dgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 577 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 577 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 578 "dgesvd.f"
		} else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */

#line 582 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 583 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
#line 583 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 584 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 584 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 585 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 585 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 586 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
#line 586 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 587 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 588 "dgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
#line 589 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 589 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 590 "dgesvd.f"
		} else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='A') */

#line 595 "dgesvd.f"
		    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
#line 596 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
#line 596 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 597 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
#line 597 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 598 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 598 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 599 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
#line 599 "dgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 600 "dgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 601 "dgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 602 "dgesvd.f"
		    i__2 = *m * 3 + *n;
#line 602 "dgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 603 "dgesvd.f"
		}
#line 604 "dgesvd.f"
	    } else {

/*              Path 10t(N greater than M, but not much larger) */

#line 608 "dgesvd.f"
		dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 610 "dgesvd.f"
		lwork_dgebrd__ = (integer) dum[0];
#line 611 "dgesvd.f"
		maxwrk = *m * 3 + lwork_dgebrd__;
#line 612 "dgesvd.f"
		if (wntvs || wntvo) {
/*                Compute space needed for DORGBR P */
#line 614 "dgesvd.f"
		    dorgbr_("P", m, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 616 "dgesvd.f"
		    lwork_dorgbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 617 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 617 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 618 "dgesvd.f"
		}
#line 619 "dgesvd.f"
		if (wntva) {
#line 620 "dgesvd.f"
		    dorgbr_("P", n, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 622 "dgesvd.f"
		    lwork_dorgbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 623 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_p__;
#line 623 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 624 "dgesvd.f"
		}
#line 625 "dgesvd.f"
		if (! wntun) {
/* Computing MAX */
#line 626 "dgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_q__;
#line 626 "dgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 627 "dgesvd.f"
		}
#line 628 "dgesvd.f"
		maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 629 "dgesvd.f"
		i__2 = *m * 3 + *n;
#line 629 "dgesvd.f"
		minwrk = max(i__2,bdspac);
#line 630 "dgesvd.f"
	    }
#line 631 "dgesvd.f"
	}
#line 632 "dgesvd.f"
	maxwrk = max(maxwrk,minwrk);
#line 633 "dgesvd.f"
	work[1] = (doublereal) maxwrk;

#line 635 "dgesvd.f"
	if (*lwork < minwrk && ! lquery) {
#line 636 "dgesvd.f"
	    *info = -13;
#line 637 "dgesvd.f"
	}
#line 638 "dgesvd.f"
    }

#line 640 "dgesvd.f"
    if (*info != 0) {
#line 641 "dgesvd.f"
	i__2 = -(*info);
#line 641 "dgesvd.f"
	xerbla_("DGESVD", &i__2, (ftnlen)6);
#line 642 "dgesvd.f"
	return 0;
#line 643 "dgesvd.f"
    } else if (lquery) {
#line 644 "dgesvd.f"
	return 0;
#line 645 "dgesvd.f"
    }

/*     Quick return if possible */

#line 649 "dgesvd.f"
    if (*m == 0 || *n == 0) {
#line 650 "dgesvd.f"
	return 0;
#line 651 "dgesvd.f"
    }

/*     Get machine constants */

#line 655 "dgesvd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 656 "dgesvd.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 657 "dgesvd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 661 "dgesvd.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 662 "dgesvd.f"
    iscl = 0;
#line 663 "dgesvd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 664 "dgesvd.f"
	iscl = 1;
#line 665 "dgesvd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 666 "dgesvd.f"
    } else if (anrm > bignum) {
#line 667 "dgesvd.f"
	iscl = 1;
#line 668 "dgesvd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 669 "dgesvd.f"
    }

#line 671 "dgesvd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 677 "dgesvd.f"
	if (*m >= mnthr) {

#line 679 "dgesvd.f"
	    if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N') */
/*              No left singular vectors to be computed */

#line 684 "dgesvd.f"
		itau = 1;
#line 685 "dgesvd.f"
		iwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need 2*N, prefer N+N*NB) */

#line 690 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 690 "dgesvd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out below R */

#line 695 "dgesvd.f"
		i__2 = *n - 1;
#line 695 "dgesvd.f"
		i__3 = *n - 1;
#line 695 "dgesvd.f"
		dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 696 "dgesvd.f"
		ie = 1;
#line 697 "dgesvd.f"
		itauq = ie + *n;
#line 698 "dgesvd.f"
		itaup = itauq + *n;
#line 699 "dgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 704 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 704 "dgesvd.f"
		dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 707 "dgesvd.f"
		ncvt = 0;
#line 708 "dgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 713 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 713 "dgesvd.f"
		    dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 715 "dgesvd.f"
		    ncvt = *n;
#line 716 "dgesvd.f"
		}
#line 717 "dgesvd.f"
		iwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 723 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, dum, &c__1, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 728 "dgesvd.f"
		if (wntvas) {
#line 728 "dgesvd.f"
		    dlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 728 "dgesvd.f"
		}

#line 731 "dgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

/* Computing MAX */
#line 737 "dgesvd.f"
		i__2 = *n << 2;
#line 737 "dgesvd.f"
		if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 741 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 742 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 742 "dgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 746 "dgesvd.f"
			ldwrku = *lda;
#line 747 "dgesvd.f"
			ldwrkr = *lda;
#line 748 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 748 "dgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 748 "dgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 752 "dgesvd.f"
			    ldwrku = *lda;
#line 753 "dgesvd.f"
			    ldwrkr = *n;
#line 754 "dgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 758 "dgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 759 "dgesvd.f"
			    ldwrkr = *n;
#line 760 "dgesvd.f"
			}
#line 760 "dgesvd.f"
		    }
#line 761 "dgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 762 "dgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 767 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 767 "dgesvd.f"
		    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 772 "dgesvd.f"
		    dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 773 "dgesvd.f"
		    i__2 = *n - 1;
#line 773 "dgesvd.f"
		    i__3 = *n - 1;
#line 773 "dgesvd.f"
		    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1], 
			    &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 779 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 779 "dgesvd.f"
		    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 781 "dgesvd.f"
		    ie = itau;
#line 782 "dgesvd.f"
		    itauq = ie + *n;
#line 783 "dgesvd.f"
		    itaup = itauq + *n;
#line 784 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 789 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 789 "dgesvd.f"
		    dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 796 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 796 "dgesvd.f"
		    dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 799 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (Workspace: need N*N+BDSPAC) */

#line 805 "dgesvd.f"
		    dbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &
			    c__1, &work[ir], &ldwrkr, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 808 "dgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N+2*N, prefer N*N+M*N+N) */

#line 814 "dgesvd.f"
		    i__2 = *m;
#line 814 "dgesvd.f"
		    i__3 = ldwrku;
#line 814 "dgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 815 "dgesvd.f"
			i__4 = *m - i__ + 1;
#line 815 "dgesvd.f"
			chunk = min(i__4,ldwrku);
#line 816 "dgesvd.f"
			dgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 819 "dgesvd.f"
			dlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 821 "dgesvd.f"
/* L10: */
#line 821 "dgesvd.f"
		    }

#line 823 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 827 "dgesvd.f"
		    ie = 1;
#line 828 "dgesvd.f"
		    itauq = ie + *n;
#line 829 "dgesvd.f"
		    itaup = itauq + *n;
#line 830 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 835 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 835 "dgesvd.f"
		    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (Workspace: need 4*N, prefer 3*N+N*NB) */

#line 842 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 842 "dgesvd.f"
		    dorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 844 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 850 "dgesvd.f"
		    dbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &
			    c__1, &a[a_offset], lda, dum, &c__1, &work[iwork],
			     info, (ftnlen)1);

#line 853 "dgesvd.f"
		}

#line 855 "dgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

/* Computing MAX */
#line 861 "dgesvd.f"
		i__3 = *n << 2;
#line 861 "dgesvd.f"
		if (*lwork >= *n * *n + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 865 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 866 "dgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 866 "dgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 870 "dgesvd.f"
			ldwrku = *lda;
#line 871 "dgesvd.f"
			ldwrkr = *lda;
#line 872 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 872 "dgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 872 "dgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 876 "dgesvd.f"
			    ldwrku = *lda;
#line 877 "dgesvd.f"
			    ldwrkr = *n;
#line 878 "dgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 882 "dgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 883 "dgesvd.f"
			    ldwrkr = *n;
#line 884 "dgesvd.f"
			}
#line 884 "dgesvd.f"
		    }
#line 885 "dgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 886 "dgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 891 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 891 "dgesvd.f"
		    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 896 "dgesvd.f"
		    dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 897 "dgesvd.f"
		    if (*n > 1) {
#line 897 "dgesvd.f"
			i__3 = *n - 1;
#line 897 "dgesvd.f"
			i__2 = *n - 1;
#line 897 "dgesvd.f"
			dlaset_("L", &i__3, &i__2, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 897 "dgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 904 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 904 "dgesvd.f"
		    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 906 "dgesvd.f"
		    ie = itau;
#line 907 "dgesvd.f"
		    itauq = ie + *n;
#line 908 "dgesvd.f"
		    itaup = itauq + *n;
#line 909 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 914 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 914 "dgesvd.f"
		    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 917 "dgesvd.f"
		    dlacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 922 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 922 "dgesvd.f"
		    dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need N*N+4*N-1, prefer N*N+3*N+(N-1)*NB) */

#line 929 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 929 "dgesvd.f"
		    dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 931 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (Workspace: need N*N+BDSPAC) */

#line 938 "dgesvd.f"
		    dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, dum, &c__1, 
			    &work[iwork], info, (ftnlen)1);
#line 941 "dgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N+2*N, prefer N*N+M*N+N) */

#line 947 "dgesvd.f"
		    i__3 = *m;
#line 947 "dgesvd.f"
		    i__2 = ldwrku;
#line 947 "dgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 948 "dgesvd.f"
			i__4 = *m - i__ + 1;
#line 948 "dgesvd.f"
			chunk = min(i__4,ldwrku);
#line 949 "dgesvd.f"
			dgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 952 "dgesvd.f"
			dlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 954 "dgesvd.f"
/* L20: */
#line 954 "dgesvd.f"
		    }

#line 956 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 960 "dgesvd.f"
		    itau = 1;
#line 961 "dgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need 2*N, prefer N+N*NB) */

#line 966 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 966 "dgesvd.f"
		    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 971 "dgesvd.f"
		    dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 972 "dgesvd.f"
		    if (*n > 1) {
#line 972 "dgesvd.f"
			i__2 = *n - 1;
#line 972 "dgesvd.f"
			i__3 = *n - 1;
#line 972 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 972 "dgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need 2*N, prefer N+N*NB) */

#line 979 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 979 "dgesvd.f"
		    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 981 "dgesvd.f"
		    ie = itau;
#line 982 "dgesvd.f"
		    itauq = ie + *n;
#line 983 "dgesvd.f"
		    itaup = itauq + *n;
#line 984 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 989 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 989 "dgesvd.f"
		    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 996 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 996 "dgesvd.f"
		    dormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1003 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1003 "dgesvd.f"
		    dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1005 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (Workspace: need BDSPAC) */

#line 1012 "dgesvd.f"
		    dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 1015 "dgesvd.f"
		}

#line 1017 "dgesvd.f"
	    } else if (wntus) {

#line 1019 "dgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1025 "dgesvd.f"
		    i__2 = *n << 2;
#line 1025 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1029 "dgesvd.f"
			ir = 1;
#line 1030 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1034 "dgesvd.f"
			    ldwrkr = *lda;
#line 1035 "dgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1039 "dgesvd.f"
			    ldwrkr = *n;
#line 1040 "dgesvd.f"
			}
#line 1041 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1042 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1047 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1047 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1052 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1054 "dgesvd.f"
			i__2 = *n - 1;
#line 1054 "dgesvd.f"
			i__3 = *n - 1;
#line 1054 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1060 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1060 "dgesvd.f"
			dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1062 "dgesvd.f"
			ie = itau;
#line 1063 "dgesvd.f"
			itauq = ie + *n;
#line 1064 "dgesvd.f"
			itaup = itauq + *n;
#line 1065 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1070 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1070 "dgesvd.f"
			dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1078 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1078 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1081 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1087 "dgesvd.f"
			dbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (Workspace: need N*N) */

#line 1095 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1098 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1102 "dgesvd.f"
			itau = 1;
#line 1103 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1108 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1108 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1110 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1115 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1115 "dgesvd.f"
			dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1117 "dgesvd.f"
			ie = itau;
#line 1118 "dgesvd.f"
			itauq = ie + *n;
#line 1119 "dgesvd.f"
			itaup = itauq + *n;
#line 1120 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1124 "dgesvd.f"
			i__2 = *n - 1;
#line 1124 "dgesvd.f"
			i__3 = *n - 1;
#line 1124 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1130 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1130 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1137 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1137 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1140 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1146 "dgesvd.f"
			dbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1150 "dgesvd.f"
		    }

#line 1152 "dgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1158 "dgesvd.f"
		    i__2 = *n << 2;
#line 1158 "dgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1162 "dgesvd.f"
			iu = 1;
#line 1163 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1167 "dgesvd.f"
			    ldwrku = *lda;
#line 1168 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1169 "dgesvd.f"
			    ldwrkr = *lda;
#line 1170 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1174 "dgesvd.f"
			    ldwrku = *lda;
#line 1175 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1176 "dgesvd.f"
			    ldwrkr = *n;
#line 1177 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1181 "dgesvd.f"
			    ldwrku = *n;
#line 1182 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1183 "dgesvd.f"
			    ldwrkr = *n;
#line 1184 "dgesvd.f"
			}
#line 1185 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1186 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1191 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1191 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1196 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1198 "dgesvd.f"
			i__2 = *n - 1;
#line 1198 "dgesvd.f"
			i__3 = *n - 1;
#line 1198 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1204 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1204 "dgesvd.f"
			dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1206 "dgesvd.f"
			ie = itau;
#line 1207 "dgesvd.f"
			itauq = ie + *n;
#line 1208 "dgesvd.f"
			itaup = itauq + *n;
#line 1209 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1216 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1216 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1220 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB) */

#line 1226 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1226 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1234 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1234 "dgesvd.f"
			dorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1237 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N+BDSPAC) */

#line 1244 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1252 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (Workspace: need N*N) */

#line 1258 "dgesvd.f"
			dlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1261 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1265 "dgesvd.f"
			itau = 1;
#line 1266 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1271 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1271 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1273 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1278 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1278 "dgesvd.f"
			dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1280 "dgesvd.f"
			ie = itau;
#line 1281 "dgesvd.f"
			itauq = ie + *n;
#line 1282 "dgesvd.f"
			itaup = itauq + *n;
#line 1283 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1287 "dgesvd.f"
			i__2 = *n - 1;
#line 1287 "dgesvd.f"
			i__3 = *n - 1;
#line 1287 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1293 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1293 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1300 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1300 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1307 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1307 "dgesvd.f"
			dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1309 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1316 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1320 "dgesvd.f"
		    }

#line 1322 "dgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1329 "dgesvd.f"
		    i__2 = *n << 2;
#line 1329 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1333 "dgesvd.f"
			iu = 1;
#line 1334 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1338 "dgesvd.f"
			    ldwrku = *lda;
#line 1339 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1343 "dgesvd.f"
			    ldwrku = *n;
#line 1344 "dgesvd.f"
			}
#line 1345 "dgesvd.f"
			itau = iu + ldwrku * *n;
#line 1346 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1351 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1351 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1356 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1358 "dgesvd.f"
			i__2 = *n - 1;
#line 1358 "dgesvd.f"
			i__3 = *n - 1;
#line 1358 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1364 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1364 "dgesvd.f"
			dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1366 "dgesvd.f"
			ie = itau;
#line 1367 "dgesvd.f"
			itauq = ie + *n;
#line 1368 "dgesvd.f"
			itaup = itauq + *n;
#line 1369 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1374 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1374 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1378 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1384 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1384 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N+4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1392 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1392 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1394 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1401 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1409 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1412 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1416 "dgesvd.f"
			itau = 1;
#line 1417 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1422 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1422 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1424 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1429 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1429 "dgesvd.f"
			dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1434 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1435 "dgesvd.f"
			if (*n > 1) {
#line 1435 "dgesvd.f"
			    i__2 = *n - 1;
#line 1435 "dgesvd.f"
			    i__3 = *n - 1;
#line 1435 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1435 "dgesvd.f"
			}
#line 1438 "dgesvd.f"
			ie = itau;
#line 1439 "dgesvd.f"
			itauq = ie + *n;
#line 1440 "dgesvd.f"
			itaup = itauq + *n;
#line 1441 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1446 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1446 "dgesvd.f"
			dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1454 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1454 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1461 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1461 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1463 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1470 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1474 "dgesvd.f"
		    }

#line 1476 "dgesvd.f"
		}

#line 1478 "dgesvd.f"
	    } else if (wntua) {

#line 1480 "dgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1486 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1486 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1490 "dgesvd.f"
			ir = 1;
#line 1491 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1495 "dgesvd.f"
			    ldwrkr = *lda;
#line 1496 "dgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1500 "dgesvd.f"
			    ldwrkr = *n;
#line 1501 "dgesvd.f"
			}
#line 1502 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1503 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1508 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1508 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1510 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1514 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1516 "dgesvd.f"
			i__2 = *n - 1;
#line 1516 "dgesvd.f"
			i__3 = *n - 1;
#line 1516 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 1522 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1522 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1524 "dgesvd.f"
			ie = itau;
#line 1525 "dgesvd.f"
			itauq = ie + *n;
#line 1526 "dgesvd.f"
			itaup = itauq + *n;
#line 1527 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1532 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1532 "dgesvd.f"
			dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1540 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1540 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1543 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1549 "dgesvd.f"
			dbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (Workspace: need N*N) */

#line 1557 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1562 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1564 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1568 "dgesvd.f"
			itau = 1;
#line 1569 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1574 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1574 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1576 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1581 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1581 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1583 "dgesvd.f"
			ie = itau;
#line 1584 "dgesvd.f"
			itauq = ie + *n;
#line 1585 "dgesvd.f"
			itaup = itauq + *n;
#line 1586 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1590 "dgesvd.f"
			i__2 = *n - 1;
#line 1590 "dgesvd.f"
			i__3 = *n - 1;
#line 1590 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1596 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1596 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1604 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1604 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1607 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1613 "dgesvd.f"
			dbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1617 "dgesvd.f"
		    }

#line 1619 "dgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1625 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1625 "dgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1629 "dgesvd.f"
			iu = 1;
#line 1630 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1634 "dgesvd.f"
			    ldwrku = *lda;
#line 1635 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1636 "dgesvd.f"
			    ldwrkr = *lda;
#line 1637 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1641 "dgesvd.f"
			    ldwrku = *lda;
#line 1642 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1643 "dgesvd.f"
			    ldwrkr = *n;
#line 1644 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1648 "dgesvd.f"
			    ldwrku = *n;
#line 1649 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1650 "dgesvd.f"
			    ldwrkr = *n;
#line 1651 "dgesvd.f"
			}
#line 1652 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1653 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1658 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1658 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1660 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */

#line 1665 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1665 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1670 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1672 "dgesvd.f"
			i__2 = *n - 1;
#line 1672 "dgesvd.f"
			i__3 = *n - 1;
#line 1672 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1674 "dgesvd.f"
			ie = itau;
#line 1675 "dgesvd.f"
			itauq = ie + *n;
#line 1676 "dgesvd.f"
			itaup = itauq + *n;
#line 1677 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1684 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1684 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1688 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB) */

#line 1694 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1694 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1702 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1702 "dgesvd.f"
			dorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1705 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N+BDSPAC) */

#line 1712 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1720 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1725 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1729 "dgesvd.f"
			dlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1732 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1736 "dgesvd.f"
			itau = 1;
#line 1737 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1742 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1742 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1744 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1749 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1749 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1751 "dgesvd.f"
			ie = itau;
#line 1752 "dgesvd.f"
			itauq = ie + *n;
#line 1753 "dgesvd.f"
			itaup = itauq + *n;
#line 1754 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1758 "dgesvd.f"
			i__2 = *n - 1;
#line 1758 "dgesvd.f"
			i__3 = *n - 1;
#line 1758 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1764 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1764 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1772 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1772 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1779 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1779 "dgesvd.f"
			dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1781 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1788 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1792 "dgesvd.f"
		    }

#line 1794 "dgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1801 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1801 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1805 "dgesvd.f"
			iu = 1;
#line 1806 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1810 "dgesvd.f"
			    ldwrku = *lda;
#line 1811 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1815 "dgesvd.f"
			    ldwrku = *n;
#line 1816 "dgesvd.f"
			}
#line 1817 "dgesvd.f"
			itau = iu + ldwrku * *n;
#line 1818 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1823 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1823 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1825 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 1830 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1830 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1835 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1837 "dgesvd.f"
			i__2 = *n - 1;
#line 1837 "dgesvd.f"
			i__3 = *n - 1;
#line 1837 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1839 "dgesvd.f"
			ie = itau;
#line 1840 "dgesvd.f"
			itauq = ie + *n;
#line 1841 "dgesvd.f"
			itaup = itauq + *n;
#line 1842 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1847 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1847 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1851 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1857 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1857 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N+4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1865 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1865 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1867 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1874 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1882 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1887 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1889 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1893 "dgesvd.f"
			itau = 1;
#line 1894 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1899 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1899 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1901 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1906 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1906 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 1911 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1912 "dgesvd.f"
			if (*n > 1) {
#line 1912 "dgesvd.f"
			    i__2 = *n - 1;
#line 1912 "dgesvd.f"
			    i__3 = *n - 1;
#line 1912 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1912 "dgesvd.f"
			}
#line 1915 "dgesvd.f"
			ie = itau;
#line 1916 "dgesvd.f"
			itauq = ie + *n;
#line 1917 "dgesvd.f"
			itaup = itauq + *n;
#line 1918 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1923 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1923 "dgesvd.f"
			dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1931 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1931 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1938 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1938 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1940 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1947 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1951 "dgesvd.f"
		    }

#line 1953 "dgesvd.f"
		}

#line 1955 "dgesvd.f"
	    }

#line 1957 "dgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 1964 "dgesvd.f"
	    ie = 1;
#line 1965 "dgesvd.f"
	    itauq = ie + *n;
#line 1966 "dgesvd.f"
	    itaup = itauq + *n;
#line 1967 "dgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 1972 "dgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 1972 "dgesvd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 1975 "dgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 3*N+NCU, prefer 3*N+NCU*NB) */

#line 1981 "dgesvd.f"
		dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1982 "dgesvd.f"
		if (wntus) {
#line 1982 "dgesvd.f"
		    ncu = *n;
#line 1982 "dgesvd.f"
		}
#line 1984 "dgesvd.f"
		if (wntua) {
#line 1984 "dgesvd.f"
		    ncu = *m;
#line 1984 "dgesvd.f"
		}
#line 1986 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 1986 "dgesvd.f"
		dorgbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1988 "dgesvd.f"
	    }
#line 1989 "dgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1995 "dgesvd.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1996 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 1996 "dgesvd.f"
		dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1998 "dgesvd.f"
	    }
#line 1999 "dgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N, prefer 3*N+N*NB) */

#line 2005 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2005 "dgesvd.f"
		dorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2007 "dgesvd.f"
	    }
#line 2008 "dgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 2014 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2014 "dgesvd.f"
		dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2016 "dgesvd.f"
	    }
#line 2017 "dgesvd.f"
	    iwork = ie + *n;
#line 2018 "dgesvd.f"
	    if (wntuas || wntuo) {
#line 2018 "dgesvd.f"
		nru = *m;
#line 2018 "dgesvd.f"
	    }
#line 2020 "dgesvd.f"
	    if (wntun) {
#line 2020 "dgesvd.f"
		nru = 0;
#line 2020 "dgesvd.f"
	    }
#line 2022 "dgesvd.f"
	    if (wntvas || wntvo) {
#line 2022 "dgesvd.f"
		ncvt = *n;
#line 2022 "dgesvd.f"
	    }
#line 2024 "dgesvd.f"
	    if (wntvn) {
#line 2024 "dgesvd.f"
		ncvt = 0;
#line 2024 "dgesvd.f"
	    }
#line 2026 "dgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2033 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2035 "dgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 2042 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 2044 "dgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2051 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2053 "dgesvd.f"
	    }

#line 2055 "dgesvd.f"
	}

#line 2057 "dgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2063 "dgesvd.f"
	if (*n >= mnthr) {

#line 2065 "dgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2070 "dgesvd.f"
		itau = 1;
#line 2071 "dgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need 2*M, prefer M+M*NB) */

#line 2076 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2076 "dgesvd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2081 "dgesvd.f"
		i__2 = *m - 1;
#line 2081 "dgesvd.f"
		i__3 = *m - 1;
#line 2081 "dgesvd.f"
		dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 2082 "dgesvd.f"
		ie = 1;
#line 2083 "dgesvd.f"
		itauq = ie + *m;
#line 2084 "dgesvd.f"
		itaup = itauq + *m;
#line 2085 "dgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2090 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2090 "dgesvd.f"
		dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2093 "dgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2098 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2098 "dgesvd.f"
		    dorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2100 "dgesvd.f"
		}
#line 2101 "dgesvd.f"
		iwork = ie + *m;
#line 2102 "dgesvd.f"
		nru = 0;
#line 2103 "dgesvd.f"
		if (wntuo || wntuas) {
#line 2103 "dgesvd.f"
		    nru = *m;
#line 2103 "dgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 2110 "dgesvd.f"
		dbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &work[ie], dum, &
			c__1, &a[a_offset], lda, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2115 "dgesvd.f"
		if (wntuas) {
#line 2115 "dgesvd.f"
		    dlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2115 "dgesvd.f"
		}

#line 2118 "dgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

/* Computing MAX */
#line 2124 "dgesvd.f"
		i__2 = *m << 2;
#line 2124 "dgesvd.f"
		if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2128 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2129 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2129 "dgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2133 "dgesvd.f"
			ldwrku = *lda;
#line 2134 "dgesvd.f"
			chunk = *n;
#line 2135 "dgesvd.f"
			ldwrkr = *lda;
#line 2136 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2136 "dgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2136 "dgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2140 "dgesvd.f"
			    ldwrku = *lda;
#line 2141 "dgesvd.f"
			    chunk = *n;
#line 2142 "dgesvd.f"
			    ldwrkr = *m;
#line 2143 "dgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2147 "dgesvd.f"
			    ldwrku = *m;
#line 2148 "dgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2149 "dgesvd.f"
			    ldwrkr = *m;
#line 2150 "dgesvd.f"
			}
#line 2150 "dgesvd.f"
		    }
#line 2151 "dgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2152 "dgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2157 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2157 "dgesvd.f"
		    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2162 "dgesvd.f"
		    dlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2163 "dgesvd.f"
		    i__2 = *m - 1;
#line 2163 "dgesvd.f"
		    i__3 = *m - 1;
#line 2163 "dgesvd.f"
		    dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2169 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2169 "dgesvd.f"
		    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2171 "dgesvd.f"
		    ie = itau;
#line 2172 "dgesvd.f"
		    itauq = ie + *m;
#line 2173 "dgesvd.f"
		    itaup = itauq + *m;
#line 2174 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2179 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2179 "dgesvd.f"
		    dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */

#line 2186 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2186 "dgesvd.f"
		    dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2189 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M+BDSPAC) */

#line 2195 "dgesvd.f"
		    dbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[
			    ir], &ldwrkr, dum, &c__1, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 2198 "dgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M*N+M) */

#line 2204 "dgesvd.f"
		    i__2 = *n;
#line 2204 "dgesvd.f"
		    i__3 = chunk;
#line 2204 "dgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2205 "dgesvd.f"
			i__4 = *n - i__ + 1;
#line 2205 "dgesvd.f"
			blk = min(i__4,chunk);
#line 2206 "dgesvd.f"
			dgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2209 "dgesvd.f"
			dlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2211 "dgesvd.f"
/* L30: */
#line 2211 "dgesvd.f"
		    }

#line 2213 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2217 "dgesvd.f"
		    ie = 1;
#line 2218 "dgesvd.f"
		    itauq = ie + *m;
#line 2219 "dgesvd.f"
		    itaup = itauq + *m;
#line 2220 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 2225 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2225 "dgesvd.f"
		    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2232 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2232 "dgesvd.f"
		    dorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2234 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2240 "dgesvd.f"
		    dbdsqr_("L", m, n, &c__0, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, dum, &c__1, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);

#line 2243 "dgesvd.f"
		}

#line 2245 "dgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

/* Computing MAX */
#line 2251 "dgesvd.f"
		i__3 = *m << 2;
#line 2251 "dgesvd.f"
		if (*lwork >= *m * *m + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2255 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2256 "dgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2256 "dgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2260 "dgesvd.f"
			ldwrku = *lda;
#line 2261 "dgesvd.f"
			chunk = *n;
#line 2262 "dgesvd.f"
			ldwrkr = *lda;
#line 2263 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2263 "dgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2263 "dgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2267 "dgesvd.f"
			    ldwrku = *lda;
#line 2268 "dgesvd.f"
			    chunk = *n;
#line 2269 "dgesvd.f"
			    ldwrkr = *m;
#line 2270 "dgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2274 "dgesvd.f"
			    ldwrku = *m;
#line 2275 "dgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2276 "dgesvd.f"
			    ldwrkr = *m;
#line 2277 "dgesvd.f"
			}
#line 2277 "dgesvd.f"
		    }
#line 2278 "dgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2279 "dgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2284 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2284 "dgesvd.f"
		    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2289 "dgesvd.f"
		    dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2290 "dgesvd.f"
		    i__3 = *m - 1;
#line 2290 "dgesvd.f"
		    i__2 = *m - 1;
#line 2290 "dgesvd.f"
		    dlaset_("U", &i__3, &i__2, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2296 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2296 "dgesvd.f"
		    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2298 "dgesvd.f"
		    ie = itau;
#line 2299 "dgesvd.f"
		    itauq = ie + *m;
#line 2300 "dgesvd.f"
		    itaup = itauq + *m;
#line 2301 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2306 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2306 "dgesvd.f"
		    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2309 "dgesvd.f"
		    dlacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */

#line 2314 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2314 "dgesvd.f"
		    dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 2321 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2321 "dgesvd.f"
		    dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2323 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M+BDSPAC) */

#line 2330 "dgesvd.f"
		    dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[ir], 
			    &ldwrkr, &u[u_offset], ldu, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);
#line 2333 "dgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M*N+M)) */

#line 2339 "dgesvd.f"
		    i__3 = *n;
#line 2339 "dgesvd.f"
		    i__2 = chunk;
#line 2339 "dgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2340 "dgesvd.f"
			i__4 = *n - i__ + 1;
#line 2340 "dgesvd.f"
			blk = min(i__4,chunk);
#line 2341 "dgesvd.f"
			dgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2344 "dgesvd.f"
			dlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2346 "dgesvd.f"
/* L40: */
#line 2346 "dgesvd.f"
		    }

#line 2348 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2352 "dgesvd.f"
		    itau = 1;
#line 2353 "dgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need 2*M, prefer M+M*NB) */

#line 2358 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2358 "dgesvd.f"
		    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2363 "dgesvd.f"
		    dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2364 "dgesvd.f"
		    i__2 = *m - 1;
#line 2364 "dgesvd.f"
		    i__3 = *m - 1;
#line 2364 "dgesvd.f"
		    dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need 2*M, prefer M+M*NB) */

#line 2370 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2370 "dgesvd.f"
		    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2372 "dgesvd.f"
		    ie = itau;
#line 2373 "dgesvd.f"
		    itauq = ie + *m;
#line 2374 "dgesvd.f"
		    itaup = itauq + *m;
#line 2375 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2380 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2380 "dgesvd.f"
		    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2387 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2387 "dgesvd.f"
		    dormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2394 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2394 "dgesvd.f"
		    dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2396 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2403 "dgesvd.f"
		    dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 2406 "dgesvd.f"
		}

#line 2408 "dgesvd.f"
	    } else if (wntvs) {

#line 2410 "dgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2416 "dgesvd.f"
		    i__2 = *m << 2;
#line 2416 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2420 "dgesvd.f"
			ir = 1;
#line 2421 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2425 "dgesvd.f"
			    ldwrkr = *lda;
#line 2426 "dgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2430 "dgesvd.f"
			    ldwrkr = *m;
#line 2431 "dgesvd.f"
			}
#line 2432 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2433 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2438 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2438 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2443 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2445 "dgesvd.f"
			i__2 = *m - 1;
#line 2445 "dgesvd.f"
			i__3 = *m - 1;
#line 2445 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2451 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2451 "dgesvd.f"
			dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2453 "dgesvd.f"
			ie = itau;
#line 2454 "dgesvd.f"
			itauq = ie + *m;
#line 2455 "dgesvd.f"
			itaup = itauq + *m;
#line 2456 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2461 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2461 "dgesvd.f"
			dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB) */

#line 2470 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2470 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2473 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2479 "dgesvd.f"
			dbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2487 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2490 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2494 "dgesvd.f"
			itau = 1;
#line 2495 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2500 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2500 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2505 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2510 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2510 "dgesvd.f"
			dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2512 "dgesvd.f"
			ie = itau;
#line 2513 "dgesvd.f"
			itauq = ie + *m;
#line 2514 "dgesvd.f"
			itaup = itauq + *m;
#line 2515 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2519 "dgesvd.f"
			i__2 = *m - 1;
#line 2519 "dgesvd.f"
			i__3 = *m - 1;
#line 2519 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2525 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2525 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2532 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2532 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2535 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2541 "dgesvd.f"
			dbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 2545 "dgesvd.f"
		    }

#line 2547 "dgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 2553 "dgesvd.f"
		    i__2 = *m << 2;
#line 2553 "dgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2557 "dgesvd.f"
			iu = 1;
#line 2558 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2562 "dgesvd.f"
			    ldwrku = *lda;
#line 2563 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2564 "dgesvd.f"
			    ldwrkr = *lda;
#line 2565 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2569 "dgesvd.f"
			    ldwrku = *lda;
#line 2570 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2571 "dgesvd.f"
			    ldwrkr = *m;
#line 2572 "dgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2576 "dgesvd.f"
			    ldwrku = *m;
#line 2577 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2578 "dgesvd.f"
			    ldwrkr = *m;
#line 2579 "dgesvd.f"
			}
#line 2580 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2581 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 2586 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2586 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2591 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2593 "dgesvd.f"
			i__2 = *m - 1;
#line 2593 "dgesvd.f"
			i__3 = *m - 1;
#line 2593 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 2599 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2599 "dgesvd.f"
			dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2601 "dgesvd.f"
			ie = itau;
#line 2602 "dgesvd.f"
			itauq = ie + *m;
#line 2603 "dgesvd.f"
			itaup = itauq + *m;
#line 2604 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 2611 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2611 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2615 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M+4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 2622 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2622 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB) */

#line 2629 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2629 "dgesvd.f"
			dorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2632 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M+BDSPAC) */

#line 2639 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2647 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (Workspace: need M*M) */

#line 2653 "dgesvd.f"
			dlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2656 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2660 "dgesvd.f"
			itau = 1;
#line 2661 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2666 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2666 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2668 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2673 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2673 "dgesvd.f"
			dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2675 "dgesvd.f"
			ie = itau;
#line 2676 "dgesvd.f"
			itauq = ie + *m;
#line 2677 "dgesvd.f"
			itaup = itauq + *m;
#line 2678 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2682 "dgesvd.f"
			i__2 = *m - 1;
#line 2682 "dgesvd.f"
			i__3 = *m - 1;
#line 2682 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2688 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2688 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2695 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2695 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2702 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2702 "dgesvd.f"
			dorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2704 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, compute left */
/*                    singular vectors of A in A and compute right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2711 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2715 "dgesvd.f"
		    }

#line 2717 "dgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 2724 "dgesvd.f"
		    i__2 = *m << 2;
#line 2724 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2728 "dgesvd.f"
			iu = 1;
#line 2729 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2733 "dgesvd.f"
			    ldwrku = *lda;
#line 2734 "dgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2738 "dgesvd.f"
			    ldwrku = *m;
#line 2739 "dgesvd.f"
			}
#line 2740 "dgesvd.f"
			itau = iu + ldwrku * *m;
#line 2741 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2746 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2746 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2751 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2753 "dgesvd.f"
			i__2 = *m - 1;
#line 2753 "dgesvd.f"
			i__3 = *m - 1;
#line 2753 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2759 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2759 "dgesvd.f"
			dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2761 "dgesvd.f"
			ie = itau;
#line 2762 "dgesvd.f"
			itauq = ie + *m;
#line 2763 "dgesvd.f"
			itaup = itauq + *m;
#line 2764 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2769 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2769 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2773 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M+4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2780 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2780 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 2787 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2787 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2789 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2796 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2804 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2807 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2811 "dgesvd.f"
			itau = 1;
#line 2812 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2817 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2817 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2819 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2824 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2824 "dgesvd.f"
			dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2829 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2830 "dgesvd.f"
			i__2 = *m - 1;
#line 2830 "dgesvd.f"
			i__3 = *m - 1;
#line 2830 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 2832 "dgesvd.f"
			ie = itau;
#line 2833 "dgesvd.f"
			itauq = ie + *m;
#line 2834 "dgesvd.f"
			itaup = itauq + *m;
#line 2835 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2840 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2840 "dgesvd.f"
			dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2848 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2848 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2855 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2855 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2857 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2864 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2868 "dgesvd.f"
		    }

#line 2870 "dgesvd.f"
		}

#line 2872 "dgesvd.f"
	    } else if (wntva) {

#line 2874 "dgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2880 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 2880 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2884 "dgesvd.f"
			ir = 1;
#line 2885 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2889 "dgesvd.f"
			    ldwrkr = *lda;
#line 2890 "dgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2894 "dgesvd.f"
			    ldwrkr = *m;
#line 2895 "dgesvd.f"
			}
#line 2896 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2897 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2902 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2902 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2904 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2908 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2910 "dgesvd.f"
			i__2 = *m - 1;
#line 2910 "dgesvd.f"
			i__3 = *m - 1;
#line 2910 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB) */

#line 2916 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2916 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2918 "dgesvd.f"
			ie = itau;
#line 2919 "dgesvd.f"
			itauq = ie + *m;
#line 2920 "dgesvd.f"
			itaup = itauq + *m;
#line 2921 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2926 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2926 "dgesvd.f"
			dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need M*M+4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2935 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2935 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2938 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2944 "dgesvd.f"
			dbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 2952 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 2957 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 2959 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2963 "dgesvd.f"
			itau = 1;
#line 2964 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2969 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2969 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2971 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 2976 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2976 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2978 "dgesvd.f"
			ie = itau;
#line 2979 "dgesvd.f"
			itauq = ie + *m;
#line 2980 "dgesvd.f"
			itaup = itauq + *m;
#line 2981 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2985 "dgesvd.f"
			i__2 = *m - 1;
#line 2985 "dgesvd.f"
			i__3 = *m - 1;
#line 2985 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2991 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2991 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2999 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2999 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3002 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3008 "dgesvd.f"
			dbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 3012 "dgesvd.f"
		    }

#line 3014 "dgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3020 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3020 "dgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3024 "dgesvd.f"
			iu = 1;
#line 3025 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3029 "dgesvd.f"
			    ldwrku = *lda;
#line 3030 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3031 "dgesvd.f"
			    ldwrkr = *lda;
#line 3032 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3036 "dgesvd.f"
			    ldwrku = *lda;
#line 3037 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3038 "dgesvd.f"
			    ldwrkr = *m;
#line 3039 "dgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3043 "dgesvd.f"
			    ldwrku = *m;
#line 3044 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3045 "dgesvd.f"
			    ldwrkr = *m;
#line 3046 "dgesvd.f"
			}
#line 3047 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3048 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 3053 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3053 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3055 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */

#line 3060 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3060 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3065 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3067 "dgesvd.f"
			i__2 = *m - 1;
#line 3067 "dgesvd.f"
			i__3 = *m - 1;
#line 3067 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3069 "dgesvd.f"
			ie = itau;
#line 3070 "dgesvd.f"
			itauq = ie + *m;
#line 3071 "dgesvd.f"
			itaup = itauq + *m;
#line 3072 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 3079 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3079 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3083 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M+4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 3090 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3090 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB) */

#line 3097 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3097 "dgesvd.f"
			dorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3100 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M+BDSPAC) */

#line 3107 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3115 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3120 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3124 "dgesvd.f"
			dlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3127 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3131 "dgesvd.f"
			itau = 1;
#line 3132 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 3137 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3137 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3139 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 3144 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3144 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3146 "dgesvd.f"
			ie = itau;
#line 3147 "dgesvd.f"
			itauq = ie + *m;
#line 3148 "dgesvd.f"
			itaup = itauq + *m;
#line 3149 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3153 "dgesvd.f"
			i__2 = *m - 1;
#line 3153 "dgesvd.f"
			i__3 = *m - 1;
#line 3153 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 3159 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3159 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3167 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3167 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3174 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3174 "dgesvd.f"
			dorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3176 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3183 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3187 "dgesvd.f"
		    }

#line 3189 "dgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3196 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3196 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3200 "dgesvd.f"
			iu = 1;
#line 3201 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3205 "dgesvd.f"
			    ldwrku = *lda;
#line 3206 "dgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3210 "dgesvd.f"
			    ldwrku = *m;
#line 3211 "dgesvd.f"
			}
#line 3212 "dgesvd.f"
			itau = iu + ldwrku * *m;
#line 3213 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 3218 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3218 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3220 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB) */

#line 3225 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3225 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3230 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3232 "dgesvd.f"
			i__2 = *m - 1;
#line 3232 "dgesvd.f"
			i__3 = *m - 1;
#line 3232 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3234 "dgesvd.f"
			ie = itau;
#line 3235 "dgesvd.f"
			itauq = ie + *m;
#line 3236 "dgesvd.f"
			itaup = itauq + *m;
#line 3237 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 3242 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3242 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3246 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB) */

#line 3252 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3252 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 3259 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3259 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3261 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 3268 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3276 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3281 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3283 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3287 "dgesvd.f"
			itau = 1;
#line 3288 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 3293 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3293 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3295 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 3300 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3300 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3305 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3306 "dgesvd.f"
			i__2 = *m - 1;
#line 3306 "dgesvd.f"
			i__3 = *m - 1;
#line 3306 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 3308 "dgesvd.f"
			ie = itau;
#line 3309 "dgesvd.f"
			itauq = ie + *m;
#line 3310 "dgesvd.f"
			itaup = itauq + *m;
#line 3311 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 3316 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3316 "dgesvd.f"
			dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3324 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3324 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3331 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3331 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3333 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3340 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3344 "dgesvd.f"
		    }

#line 3346 "dgesvd.f"
		}

#line 3348 "dgesvd.f"
	    }

#line 3350 "dgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3357 "dgesvd.f"
	    ie = 1;
#line 3358 "dgesvd.f"
	    itauq = ie + *m;
#line 3359 "dgesvd.f"
	    itaup = itauq + *m;
#line 3360 "dgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 3365 "dgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3365 "dgesvd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 3368 "dgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) */

#line 3374 "dgesvd.f"
		dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3375 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3375 "dgesvd.f"
		dorgbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3377 "dgesvd.f"
	    }
#line 3378 "dgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB) */

#line 3384 "dgesvd.f"
		dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3385 "dgesvd.f"
		if (wntva) {
#line 3385 "dgesvd.f"
		    nrvt = *n;
#line 3385 "dgesvd.f"
		}
#line 3387 "dgesvd.f"
		if (wntvs) {
#line 3387 "dgesvd.f"
		    nrvt = *m;
#line 3387 "dgesvd.f"
		}
#line 3389 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3389 "dgesvd.f"
		dorgbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3391 "dgesvd.f"
	    }
#line 3392 "dgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) */

#line 3398 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3398 "dgesvd.f"
		dorgbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3400 "dgesvd.f"
	    }
#line 3401 "dgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3407 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3407 "dgesvd.f"
		dorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3409 "dgesvd.f"
	    }
#line 3410 "dgesvd.f"
	    iwork = ie + *m;
#line 3411 "dgesvd.f"
	    if (wntuas || wntuo) {
#line 3411 "dgesvd.f"
		nru = *m;
#line 3411 "dgesvd.f"
	    }
#line 3413 "dgesvd.f"
	    if (wntun) {
#line 3413 "dgesvd.f"
		nru = 0;
#line 3413 "dgesvd.f"
	    }
#line 3415 "dgesvd.f"
	    if (wntvas || wntvo) {
#line 3415 "dgesvd.f"
		ncvt = *n;
#line 3415 "dgesvd.f"
	    }
#line 3417 "dgesvd.f"
	    if (wntvn) {
#line 3417 "dgesvd.f"
		ncvt = 0;
#line 3417 "dgesvd.f"
	    }
#line 3419 "dgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3426 "dgesvd.f"
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3428 "dgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 3435 "dgesvd.f"
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 3437 "dgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3444 "dgesvd.f"
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3446 "dgesvd.f"
	    }

#line 3448 "dgesvd.f"
	}

#line 3450 "dgesvd.f"
    }

/*     If DBDSQR failed to converge, copy unconverged superdiagonals */
/*     to WORK( 2:MINMN ) */

#line 3455 "dgesvd.f"
    if (*info != 0) {
#line 3456 "dgesvd.f"
	if (ie > 2) {
#line 3457 "dgesvd.f"
	    i__2 = minmn - 1;
#line 3457 "dgesvd.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 3458 "dgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3459 "dgesvd.f"
/* L50: */
#line 3459 "dgesvd.f"
	    }
#line 3460 "dgesvd.f"
	}
#line 3461 "dgesvd.f"
	if (ie < 2) {
#line 3462 "dgesvd.f"
	    for (i__ = minmn - 1; i__ >= 1; --i__) {
#line 3463 "dgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3464 "dgesvd.f"
/* L60: */
#line 3464 "dgesvd.f"
	    }
#line 3465 "dgesvd.f"
	}
#line 3466 "dgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3470 "dgesvd.f"
    if (iscl == 1) {
#line 3471 "dgesvd.f"
	if (anrm > bignum) {
#line 3471 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3471 "dgesvd.f"
	}
#line 3474 "dgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3474 "dgesvd.f"
	    i__2 = minmn - 1;
#line 3474 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3474 "dgesvd.f"
	}
#line 3477 "dgesvd.f"
	if (anrm < smlnum) {
#line 3477 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3477 "dgesvd.f"
	}
#line 3480 "dgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3480 "dgesvd.f"
	    i__2 = minmn - 1;
#line 3480 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3480 "dgesvd.f"
	}
#line 3483 "dgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3487 "dgesvd.f"
    work[1] = (doublereal) maxwrk;

#line 3489 "dgesvd.f"
    return 0;

/*     End of DGESVD */

} /* dgesvd_ */

