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
/* >          LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths */
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


/*  -- LAPACK driver routine (version 3.7.0) -- */
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
/*              (Workspace: need 2*N, prefer N + N*NB) */

#line 690 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 690 "dgesvd.f"
		dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out below R */

#line 695 "dgesvd.f"
		if (*n > 1) {
#line 696 "dgesvd.f"
		    i__2 = *n - 1;
#line 696 "dgesvd.f"
		    i__3 = *n - 1;
#line 696 "dgesvd.f"
		    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2],
			     lda, (ftnlen)1);
#line 698 "dgesvd.f"
		}
#line 699 "dgesvd.f"
		ie = 1;
#line 700 "dgesvd.f"
		itauq = ie + *n;
#line 701 "dgesvd.f"
		itaup = itauq + *n;
#line 702 "dgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 707 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 707 "dgesvd.f"
		dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 710 "dgesvd.f"
		ncvt = 0;
#line 711 "dgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 716 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 716 "dgesvd.f"
		    dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 718 "dgesvd.f"
		    ncvt = *n;
#line 719 "dgesvd.f"
		}
#line 720 "dgesvd.f"
		iwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 726 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, dum, &c__1, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 731 "dgesvd.f"
		if (wntvas) {
#line 731 "dgesvd.f"
		    dlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 731 "dgesvd.f"
		}

#line 734 "dgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

/* Computing MAX */
#line 740 "dgesvd.f"
		i__2 = *n << 2;
#line 740 "dgesvd.f"
		if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 744 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 745 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 745 "dgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 749 "dgesvd.f"
			ldwrku = *lda;
#line 750 "dgesvd.f"
			ldwrkr = *lda;
#line 751 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 751 "dgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 751 "dgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 755 "dgesvd.f"
			    ldwrku = *lda;
#line 756 "dgesvd.f"
			    ldwrkr = *n;
#line 757 "dgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 761 "dgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 762 "dgesvd.f"
			    ldwrkr = *n;
#line 763 "dgesvd.f"
			}
#line 763 "dgesvd.f"
		    }
#line 764 "dgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 765 "dgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 770 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 770 "dgesvd.f"
		    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 775 "dgesvd.f"
		    dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 776 "dgesvd.f"
		    i__2 = *n - 1;
#line 776 "dgesvd.f"
		    i__3 = *n - 1;
#line 776 "dgesvd.f"
		    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1], 
			    &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 782 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 782 "dgesvd.f"
		    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 784 "dgesvd.f"
		    ie = itau;
#line 785 "dgesvd.f"
		    itauq = ie + *n;
#line 786 "dgesvd.f"
		    itaup = itauq + *n;
#line 787 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

#line 792 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 792 "dgesvd.f"
		    dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

#line 799 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 799 "dgesvd.f"
		    dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 802 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (Workspace: need N*N + BDSPAC) */

#line 808 "dgesvd.f"
		    dbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &
			    c__1, &work[ir], &ldwrkr, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 811 "dgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N) */

#line 817 "dgesvd.f"
		    i__2 = *m;
#line 817 "dgesvd.f"
		    i__3 = ldwrku;
#line 817 "dgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 818 "dgesvd.f"
			i__4 = *m - i__ + 1;
#line 818 "dgesvd.f"
			chunk = min(i__4,ldwrku);
#line 819 "dgesvd.f"
			dgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 822 "dgesvd.f"
			dlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 824 "dgesvd.f"
/* L10: */
#line 824 "dgesvd.f"
		    }

#line 826 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 830 "dgesvd.f"
		    ie = 1;
#line 831 "dgesvd.f"
		    itauq = ie + *n;
#line 832 "dgesvd.f"
		    itaup = itauq + *n;
#line 833 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB) */

#line 838 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 838 "dgesvd.f"
		    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (Workspace: need 4*N, prefer 3*N + N*NB) */

#line 845 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 845 "dgesvd.f"
		    dorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 847 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 853 "dgesvd.f"
		    dbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &
			    c__1, &a[a_offset], lda, dum, &c__1, &work[iwork],
			     info, (ftnlen)1);

#line 856 "dgesvd.f"
		}

#line 858 "dgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

/* Computing MAX */
#line 864 "dgesvd.f"
		i__3 = *n << 2;
#line 864 "dgesvd.f"
		if (*lwork >= *n * *n + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 868 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 869 "dgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 869 "dgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 873 "dgesvd.f"
			ldwrku = *lda;
#line 874 "dgesvd.f"
			ldwrkr = *lda;
#line 875 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 875 "dgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 875 "dgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 879 "dgesvd.f"
			    ldwrku = *lda;
#line 880 "dgesvd.f"
			    ldwrkr = *n;
#line 881 "dgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 885 "dgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 886 "dgesvd.f"
			    ldwrkr = *n;
#line 887 "dgesvd.f"
			}
#line 887 "dgesvd.f"
		    }
#line 888 "dgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 889 "dgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 894 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 894 "dgesvd.f"
		    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 899 "dgesvd.f"
		    dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 900 "dgesvd.f"
		    if (*n > 1) {
#line 900 "dgesvd.f"
			i__3 = *n - 1;
#line 900 "dgesvd.f"
			i__2 = *n - 1;
#line 900 "dgesvd.f"
			dlaset_("L", &i__3, &i__2, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 900 "dgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 907 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 907 "dgesvd.f"
		    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 909 "dgesvd.f"
		    ie = itau;
#line 910 "dgesvd.f"
		    itauq = ie + *n;
#line 911 "dgesvd.f"
		    itaup = itauq + *n;
#line 912 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

#line 917 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 917 "dgesvd.f"
		    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 920 "dgesvd.f"
		    dlacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

#line 925 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 925 "dgesvd.f"
		    dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need N*N + 4*N-1, prefer N*N + 3*N + (N-1)*NB) */

#line 932 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 932 "dgesvd.f"
		    dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 934 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (Workspace: need N*N + BDSPAC) */

#line 941 "dgesvd.f"
		    dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, dum, &c__1, 
			    &work[iwork], info, (ftnlen)1);
#line 944 "dgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N) */

#line 950 "dgesvd.f"
		    i__3 = *m;
#line 950 "dgesvd.f"
		    i__2 = ldwrku;
#line 950 "dgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 951 "dgesvd.f"
			i__4 = *m - i__ + 1;
#line 951 "dgesvd.f"
			chunk = min(i__4,ldwrku);
#line 952 "dgesvd.f"
			dgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 955 "dgesvd.f"
			dlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 957 "dgesvd.f"
/* L20: */
#line 957 "dgesvd.f"
		    }

#line 959 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 963 "dgesvd.f"
		    itau = 1;
#line 964 "dgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need 2*N, prefer N + N*NB) */

#line 969 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 969 "dgesvd.f"
		    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 974 "dgesvd.f"
		    dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 975 "dgesvd.f"
		    if (*n > 1) {
#line 975 "dgesvd.f"
			i__2 = *n - 1;
#line 975 "dgesvd.f"
			i__3 = *n - 1;
#line 975 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 975 "dgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need 2*N, prefer N + N*NB) */

#line 982 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 982 "dgesvd.f"
		    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 984 "dgesvd.f"
		    ie = itau;
#line 985 "dgesvd.f"
		    itauq = ie + *n;
#line 986 "dgesvd.f"
		    itaup = itauq + *n;
#line 987 "dgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 992 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 992 "dgesvd.f"
		    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (Workspace: need 3*N + M, prefer 3*N + M*NB) */

#line 999 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 999 "dgesvd.f"
		    dormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 1006 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1006 "dgesvd.f"
		    dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1008 "dgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (Workspace: need BDSPAC) */

#line 1015 "dgesvd.f"
		    dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 1018 "dgesvd.f"
		}

#line 1020 "dgesvd.f"
	    } else if (wntus) {

#line 1022 "dgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1028 "dgesvd.f"
		    i__2 = *n << 2;
#line 1028 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1032 "dgesvd.f"
			ir = 1;
#line 1033 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1037 "dgesvd.f"
			    ldwrkr = *lda;
#line 1038 "dgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1042 "dgesvd.f"
			    ldwrkr = *n;
#line 1043 "dgesvd.f"
			}
#line 1044 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1045 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 1050 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1050 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1055 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1057 "dgesvd.f"
			i__2 = *n - 1;
#line 1057 "dgesvd.f"
			i__3 = *n - 1;
#line 1057 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 1063 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1063 "dgesvd.f"
			dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1065 "dgesvd.f"
			ie = itau;
#line 1066 "dgesvd.f"
			itauq = ie + *n;
#line 1067 "dgesvd.f"
			itaup = itauq + *n;
#line 1068 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

#line 1073 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1073 "dgesvd.f"
			dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

#line 1081 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1081 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1084 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N + BDSPAC) */

#line 1090 "dgesvd.f"
			dbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (Workspace: need N*N) */

#line 1098 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1101 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1105 "dgesvd.f"
			itau = 1;
#line 1106 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1111 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1111 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1113 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1118 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1118 "dgesvd.f"
			dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1120 "dgesvd.f"
			ie = itau;
#line 1121 "dgesvd.f"
			itauq = ie + *n;
#line 1122 "dgesvd.f"
			itaup = itauq + *n;
#line 1123 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1127 "dgesvd.f"
			if (*n > 1) {
#line 1128 "dgesvd.f"
			    i__2 = *n - 1;
#line 1128 "dgesvd.f"
			    i__3 = *n - 1;
#line 1128 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1130 "dgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 1135 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1135 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

#line 1142 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1142 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1145 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1151 "dgesvd.f"
			dbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1155 "dgesvd.f"
		    }

#line 1157 "dgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1163 "dgesvd.f"
		    i__2 = *n << 2;
#line 1163 "dgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1167 "dgesvd.f"
			iu = 1;
#line 1168 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1172 "dgesvd.f"
			    ldwrku = *lda;
#line 1173 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1174 "dgesvd.f"
			    ldwrkr = *lda;
#line 1175 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1179 "dgesvd.f"
			    ldwrku = *lda;
#line 1180 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1181 "dgesvd.f"
			    ldwrkr = *n;
#line 1182 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1186 "dgesvd.f"
			    ldwrku = *n;
#line 1187 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1188 "dgesvd.f"
			    ldwrkr = *n;
#line 1189 "dgesvd.f"
			}
#line 1190 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1191 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB) */

#line 1196 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1196 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1201 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1203 "dgesvd.f"
			i__2 = *n - 1;
#line 1203 "dgesvd.f"
			i__3 = *n - 1;
#line 1203 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB) */

#line 1209 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1209 "dgesvd.f"
			dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1211 "dgesvd.f"
			ie = itau;
#line 1212 "dgesvd.f"
			itauq = ie + *n;
#line 1213 "dgesvd.f"
			itaup = itauq + *n;
#line 1214 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1221 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1221 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1225 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB) */

#line 1231 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1231 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1239 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1239 "dgesvd.f"
			dorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1242 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N + BDSPAC) */

#line 1249 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1257 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (Workspace: need N*N) */

#line 1263 "dgesvd.f"
			dlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1266 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1270 "dgesvd.f"
			itau = 1;
#line 1271 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1276 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1276 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1278 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1283 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1283 "dgesvd.f"
			dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1285 "dgesvd.f"
			ie = itau;
#line 1286 "dgesvd.f"
			itauq = ie + *n;
#line 1287 "dgesvd.f"
			itaup = itauq + *n;
#line 1288 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1292 "dgesvd.f"
			if (*n > 1) {
#line 1293 "dgesvd.f"
			    i__2 = *n - 1;
#line 1293 "dgesvd.f"
			    i__3 = *n - 1;
#line 1293 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1295 "dgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 1300 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1300 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

#line 1307 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1307 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 1314 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1314 "dgesvd.f"
			dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1316 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1323 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1327 "dgesvd.f"
		    }

#line 1329 "dgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1336 "dgesvd.f"
		    i__2 = *n << 2;
#line 1336 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1340 "dgesvd.f"
			iu = 1;
#line 1341 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1345 "dgesvd.f"
			    ldwrku = *lda;
#line 1346 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1350 "dgesvd.f"
			    ldwrku = *n;
#line 1351 "dgesvd.f"
			}
#line 1352 "dgesvd.f"
			itau = iu + ldwrku * *n;
#line 1353 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 1358 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1358 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1363 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1365 "dgesvd.f"
			i__2 = *n - 1;
#line 1365 "dgesvd.f"
			i__3 = *n - 1;
#line 1365 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 1371 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1371 "dgesvd.f"
			dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1373 "dgesvd.f"
			ie = itau;
#line 1374 "dgesvd.f"
			itauq = ie + *n;
#line 1375 "dgesvd.f"
			itaup = itauq + *n;
#line 1376 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

#line 1381 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1381 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1385 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

#line 1391 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1391 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N + 4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1399 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1399 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1401 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N + BDSPAC) */

#line 1408 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1416 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1419 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1423 "dgesvd.f"
			itau = 1;
#line 1424 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1429 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1429 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1431 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1436 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1436 "dgesvd.f"
			dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1441 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1442 "dgesvd.f"
			if (*n > 1) {
#line 1442 "dgesvd.f"
			    i__2 = *n - 1;
#line 1442 "dgesvd.f"
			    i__3 = *n - 1;
#line 1442 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1442 "dgesvd.f"
			}
#line 1445 "dgesvd.f"
			ie = itau;
#line 1446 "dgesvd.f"
			itauq = ie + *n;
#line 1447 "dgesvd.f"
			itaup = itauq + *n;
#line 1448 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 1453 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1453 "dgesvd.f"
			dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

#line 1461 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1461 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 1468 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1468 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1470 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1477 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1481 "dgesvd.f"
		    }

#line 1483 "dgesvd.f"
		}

#line 1485 "dgesvd.f"
	    } else if (wntua) {

#line 1487 "dgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1493 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1493 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1497 "dgesvd.f"
			ir = 1;
#line 1498 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1502 "dgesvd.f"
			    ldwrkr = *lda;
#line 1503 "dgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1507 "dgesvd.f"
			    ldwrkr = *n;
#line 1508 "dgesvd.f"
			}
#line 1509 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1510 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 1515 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1515 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1517 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1521 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1523 "dgesvd.f"
			i__2 = *n - 1;
#line 1523 "dgesvd.f"
			i__3 = *n - 1;
#line 1523 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB) */

#line 1529 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1529 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1531 "dgesvd.f"
			ie = itau;
#line 1532 "dgesvd.f"
			itauq = ie + *n;
#line 1533 "dgesvd.f"
			itaup = itauq + *n;
#line 1534 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

#line 1539 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1539 "dgesvd.f"
			dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

#line 1547 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1547 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1550 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N + BDSPAC) */

#line 1556 "dgesvd.f"
			dbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (Workspace: need N*N) */

#line 1564 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1569 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1571 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1575 "dgesvd.f"
			itau = 1;
#line 1576 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1581 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1581 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1583 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N + M, prefer N + M*NB) */

#line 1588 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1588 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1590 "dgesvd.f"
			ie = itau;
#line 1591 "dgesvd.f"
			itauq = ie + *n;
#line 1592 "dgesvd.f"
			itaup = itauq + *n;
#line 1593 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1597 "dgesvd.f"
			if (*n > 1) {
#line 1598 "dgesvd.f"
			    i__2 = *n - 1;
#line 1598 "dgesvd.f"
			    i__3 = *n - 1;
#line 1598 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1600 "dgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 1605 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1605 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

#line 1613 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1613 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1616 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1622 "dgesvd.f"
			dbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1626 "dgesvd.f"
		    }

#line 1628 "dgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1634 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1634 "dgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1638 "dgesvd.f"
			iu = 1;
#line 1639 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1643 "dgesvd.f"
			    ldwrku = *lda;
#line 1644 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1645 "dgesvd.f"
			    ldwrkr = *lda;
#line 1646 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1650 "dgesvd.f"
			    ldwrku = *lda;
#line 1651 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1652 "dgesvd.f"
			    ldwrkr = *n;
#line 1653 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1657 "dgesvd.f"
			    ldwrku = *n;
#line 1658 "dgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1659 "dgesvd.f"
			    ldwrkr = *n;
#line 1660 "dgesvd.f"
			}
#line 1661 "dgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1662 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB) */

#line 1667 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1667 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1669 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N*N + N + M, prefer 2*N*N + N + M*NB) */

#line 1674 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1674 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1679 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1681 "dgesvd.f"
			i__2 = *n - 1;
#line 1681 "dgesvd.f"
			i__3 = *n - 1;
#line 1681 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1683 "dgesvd.f"
			ie = itau;
#line 1684 "dgesvd.f"
			itauq = ie + *n;
#line 1685 "dgesvd.f"
			itaup = itauq + *n;
#line 1686 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1693 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1693 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1697 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB) */

#line 1703 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1703 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1711 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1711 "dgesvd.f"
			dorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1714 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N + BDSPAC) */

#line 1721 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1729 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1734 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1738 "dgesvd.f"
			dlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1741 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1745 "dgesvd.f"
			itau = 1;
#line 1746 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1751 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1751 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1753 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N + M, prefer N + M*NB) */

#line 1758 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1758 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1760 "dgesvd.f"
			ie = itau;
#line 1761 "dgesvd.f"
			itauq = ie + *n;
#line 1762 "dgesvd.f"
			itaup = itauq + *n;
#line 1763 "dgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1767 "dgesvd.f"
			if (*n > 1) {
#line 1768 "dgesvd.f"
			    i__2 = *n - 1;
#line 1768 "dgesvd.f"
			    i__3 = *n - 1;
#line 1768 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1770 "dgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 1775 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1775 "dgesvd.f"
			dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

#line 1783 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1783 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 1790 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1790 "dgesvd.f"
			dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1792 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1799 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1803 "dgesvd.f"
		    }

#line 1805 "dgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1812 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1812 "dgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1816 "dgesvd.f"
			iu = 1;
#line 1817 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1821 "dgesvd.f"
			    ldwrku = *lda;
#line 1822 "dgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1826 "dgesvd.f"
			    ldwrku = *n;
#line 1827 "dgesvd.f"
			}
#line 1828 "dgesvd.f"
			itau = iu + ldwrku * *n;
#line 1829 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

#line 1834 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1834 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1836 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB) */

#line 1841 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1841 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1846 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1848 "dgesvd.f"
			i__2 = *n - 1;
#line 1848 "dgesvd.f"
			i__3 = *n - 1;
#line 1848 "dgesvd.f"
			dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1850 "dgesvd.f"
			ie = itau;
#line 1851 "dgesvd.f"
			itauq = ie + *n;
#line 1852 "dgesvd.f"
			itaup = itauq + *n;
#line 1853 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

#line 1858 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1858 "dgesvd.f"
			dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1862 "dgesvd.f"
			dlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

#line 1868 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1868 "dgesvd.f"
			dorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N + 4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1876 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1876 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1878 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N + BDSPAC) */

#line 1885 "dgesvd.f"
			dbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1893 "dgesvd.f"
			dgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1898 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1900 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1904 "dgesvd.f"
			itau = 1;
#line 1905 "dgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

#line 1910 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1910 "dgesvd.f"
			dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1912 "dgesvd.f"
			dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N + M, prefer N + M*NB) */

#line 1917 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1917 "dgesvd.f"
			dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 1922 "dgesvd.f"
			dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1923 "dgesvd.f"
			if (*n > 1) {
#line 1923 "dgesvd.f"
			    i__2 = *n - 1;
#line 1923 "dgesvd.f"
			    i__3 = *n - 1;
#line 1923 "dgesvd.f"
			    dlaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1923 "dgesvd.f"
			}
#line 1926 "dgesvd.f"
			ie = itau;
#line 1927 "dgesvd.f"
			itauq = ie + *n;
#line 1928 "dgesvd.f"
			itaup = itauq + *n;
#line 1929 "dgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

#line 1934 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1934 "dgesvd.f"
			dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

#line 1942 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1942 "dgesvd.f"
			dormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 1949 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1949 "dgesvd.f"
			dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1951 "dgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1958 "dgesvd.f"
			dbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1962 "dgesvd.f"
		    }

#line 1964 "dgesvd.f"
		}

#line 1966 "dgesvd.f"
	    }

#line 1968 "dgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 1975 "dgesvd.f"
	    ie = 1;
#line 1976 "dgesvd.f"
	    itauq = ie + *n;
#line 1977 "dgesvd.f"
	    itaup = itauq + *n;
#line 1978 "dgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB) */

#line 1983 "dgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 1983 "dgesvd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 1986 "dgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 3*N + NCU, prefer 3*N + NCU*NB) */

#line 1992 "dgesvd.f"
		dlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1993 "dgesvd.f"
		if (wntus) {
#line 1993 "dgesvd.f"
		    ncu = *n;
#line 1993 "dgesvd.f"
		}
#line 1995 "dgesvd.f"
		if (wntua) {
#line 1995 "dgesvd.f"
		    ncu = *m;
#line 1995 "dgesvd.f"
		}
#line 1997 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 1997 "dgesvd.f"
		dorgbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1999 "dgesvd.f"
	    }
#line 2000 "dgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 2006 "dgesvd.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2007 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2007 "dgesvd.f"
		dorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2009 "dgesvd.f"
	    }
#line 2010 "dgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N, prefer 3*N + N*NB) */

#line 2016 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2016 "dgesvd.f"
		dorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2018 "dgesvd.f"
	    }
#line 2019 "dgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

#line 2025 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2025 "dgesvd.f"
		dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2027 "dgesvd.f"
	    }
#line 2028 "dgesvd.f"
	    iwork = ie + *n;
#line 2029 "dgesvd.f"
	    if (wntuas || wntuo) {
#line 2029 "dgesvd.f"
		nru = *m;
#line 2029 "dgesvd.f"
	    }
#line 2031 "dgesvd.f"
	    if (wntun) {
#line 2031 "dgesvd.f"
		nru = 0;
#line 2031 "dgesvd.f"
	    }
#line 2033 "dgesvd.f"
	    if (wntvas || wntvo) {
#line 2033 "dgesvd.f"
		ncvt = *n;
#line 2033 "dgesvd.f"
	    }
#line 2035 "dgesvd.f"
	    if (wntvn) {
#line 2035 "dgesvd.f"
		ncvt = 0;
#line 2035 "dgesvd.f"
	    }
#line 2037 "dgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2044 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2046 "dgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 2053 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 2055 "dgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2062 "dgesvd.f"
		dbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2064 "dgesvd.f"
	    }

#line 2066 "dgesvd.f"
	}

#line 2068 "dgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2074 "dgesvd.f"
	if (*n >= mnthr) {

#line 2076 "dgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2081 "dgesvd.f"
		itau = 1;
#line 2082 "dgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need 2*M, prefer M + M*NB) */

#line 2087 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2087 "dgesvd.f"
		dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2092 "dgesvd.f"
		i__2 = *m - 1;
#line 2092 "dgesvd.f"
		i__3 = *m - 1;
#line 2092 "dgesvd.f"
		dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 2093 "dgesvd.f"
		ie = 1;
#line 2094 "dgesvd.f"
		itauq = ie + *m;
#line 2095 "dgesvd.f"
		itaup = itauq + *m;
#line 2096 "dgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 2101 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2101 "dgesvd.f"
		dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2104 "dgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 2109 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2109 "dgesvd.f"
		    dorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2111 "dgesvd.f"
		}
#line 2112 "dgesvd.f"
		iwork = ie + *m;
#line 2113 "dgesvd.f"
		nru = 0;
#line 2114 "dgesvd.f"
		if (wntuo || wntuas) {
#line 2114 "dgesvd.f"
		    nru = *m;
#line 2114 "dgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 2121 "dgesvd.f"
		dbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &work[ie], dum, &
			c__1, &a[a_offset], lda, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2126 "dgesvd.f"
		if (wntuas) {
#line 2126 "dgesvd.f"
		    dlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2126 "dgesvd.f"
		}

#line 2129 "dgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

/* Computing MAX */
#line 2135 "dgesvd.f"
		i__2 = *m << 2;
#line 2135 "dgesvd.f"
		if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2139 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2140 "dgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2140 "dgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2144 "dgesvd.f"
			ldwrku = *lda;
#line 2145 "dgesvd.f"
			chunk = *n;
#line 2146 "dgesvd.f"
			ldwrkr = *lda;
#line 2147 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2147 "dgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2147 "dgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2151 "dgesvd.f"
			    ldwrku = *lda;
#line 2152 "dgesvd.f"
			    chunk = *n;
#line 2153 "dgesvd.f"
			    ldwrkr = *m;
#line 2154 "dgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2158 "dgesvd.f"
			    ldwrku = *m;
#line 2159 "dgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2160 "dgesvd.f"
			    ldwrkr = *m;
#line 2161 "dgesvd.f"
			}
#line 2161 "dgesvd.f"
		    }
#line 2162 "dgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2163 "dgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2168 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2168 "dgesvd.f"
		    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2173 "dgesvd.f"
		    dlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2174 "dgesvd.f"
		    i__2 = *m - 1;
#line 2174 "dgesvd.f"
		    i__3 = *m - 1;
#line 2174 "dgesvd.f"
		    dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2180 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2180 "dgesvd.f"
		    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2182 "dgesvd.f"
		    ie = itau;
#line 2183 "dgesvd.f"
		    itauq = ie + *m;
#line 2184 "dgesvd.f"
		    itaup = itauq + *m;
#line 2185 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

#line 2190 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2190 "dgesvd.f"
		    dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB) */

#line 2197 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2197 "dgesvd.f"
		    dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2200 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M + BDSPAC) */

#line 2206 "dgesvd.f"
		    dbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[
			    ir], &ldwrkr, dum, &c__1, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 2209 "dgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M) */

#line 2215 "dgesvd.f"
		    i__2 = *n;
#line 2215 "dgesvd.f"
		    i__3 = chunk;
#line 2215 "dgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2216 "dgesvd.f"
			i__4 = *n - i__ + 1;
#line 2216 "dgesvd.f"
			blk = min(i__4,chunk);
#line 2217 "dgesvd.f"
			dgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2220 "dgesvd.f"
			dlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2222 "dgesvd.f"
/* L30: */
#line 2222 "dgesvd.f"
		    }

#line 2224 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2228 "dgesvd.f"
		    ie = 1;
#line 2229 "dgesvd.f"
		    itauq = ie + *m;
#line 2230 "dgesvd.f"
		    itaup = itauq + *m;
#line 2231 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB) */

#line 2236 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2236 "dgesvd.f"
		    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 2243 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2243 "dgesvd.f"
		    dorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2245 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2251 "dgesvd.f"
		    dbdsqr_("L", m, n, &c__0, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, dum, &c__1, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);

#line 2254 "dgesvd.f"
		}

#line 2256 "dgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

/* Computing MAX */
#line 2262 "dgesvd.f"
		i__3 = *m << 2;
#line 2262 "dgesvd.f"
		if (*lwork >= *m * *m + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2266 "dgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2267 "dgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2267 "dgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2271 "dgesvd.f"
			ldwrku = *lda;
#line 2272 "dgesvd.f"
			chunk = *n;
#line 2273 "dgesvd.f"
			ldwrkr = *lda;
#line 2274 "dgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2274 "dgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2274 "dgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2278 "dgesvd.f"
			    ldwrku = *lda;
#line 2279 "dgesvd.f"
			    chunk = *n;
#line 2280 "dgesvd.f"
			    ldwrkr = *m;
#line 2281 "dgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2285 "dgesvd.f"
			    ldwrku = *m;
#line 2286 "dgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2287 "dgesvd.f"
			    ldwrkr = *m;
#line 2288 "dgesvd.f"
			}
#line 2288 "dgesvd.f"
		    }
#line 2289 "dgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2290 "dgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2295 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2295 "dgesvd.f"
		    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2300 "dgesvd.f"
		    dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2301 "dgesvd.f"
		    i__3 = *m - 1;
#line 2301 "dgesvd.f"
		    i__2 = *m - 1;
#line 2301 "dgesvd.f"
		    dlaset_("U", &i__3, &i__2, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2307 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2307 "dgesvd.f"
		    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2309 "dgesvd.f"
		    ie = itau;
#line 2310 "dgesvd.f"
		    itauq = ie + *m;
#line 2311 "dgesvd.f"
		    itaup = itauq + *m;
#line 2312 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

#line 2317 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2317 "dgesvd.f"
		    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2320 "dgesvd.f"
		    dlacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB) */

#line 2325 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2325 "dgesvd.f"
		    dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB) */

#line 2332 "dgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2332 "dgesvd.f"
		    dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2334 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M + BDSPAC) */

#line 2341 "dgesvd.f"
		    dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[ir], 
			    &ldwrkr, &u[u_offset], ldu, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);
#line 2344 "dgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M)) */

#line 2350 "dgesvd.f"
		    i__3 = *n;
#line 2350 "dgesvd.f"
		    i__2 = chunk;
#line 2350 "dgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2351 "dgesvd.f"
			i__4 = *n - i__ + 1;
#line 2351 "dgesvd.f"
			blk = min(i__4,chunk);
#line 2352 "dgesvd.f"
			dgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2355 "dgesvd.f"
			dlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2357 "dgesvd.f"
/* L40: */
#line 2357 "dgesvd.f"
		    }

#line 2359 "dgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2363 "dgesvd.f"
		    itau = 1;
#line 2364 "dgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need 2*M, prefer M + M*NB) */

#line 2369 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2369 "dgesvd.f"
		    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2374 "dgesvd.f"
		    dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2375 "dgesvd.f"
		    i__2 = *m - 1;
#line 2375 "dgesvd.f"
		    i__3 = *m - 1;
#line 2375 "dgesvd.f"
		    dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need 2*M, prefer M + M*NB) */

#line 2381 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2381 "dgesvd.f"
		    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2383 "dgesvd.f"
		    ie = itau;
#line 2384 "dgesvd.f"
		    itauq = ie + *m;
#line 2385 "dgesvd.f"
		    itaup = itauq + *m;
#line 2386 "dgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 2391 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2391 "dgesvd.f"
		    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (Workspace: need 3*M + N, prefer 3*M + N*NB) */

#line 2398 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2398 "dgesvd.f"
		    dormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 2405 "dgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2405 "dgesvd.f"
		    dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2407 "dgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2414 "dgesvd.f"
		    dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 2417 "dgesvd.f"
		}

#line 2419 "dgesvd.f"
	    } else if (wntvs) {

#line 2421 "dgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2427 "dgesvd.f"
		    i__2 = *m << 2;
#line 2427 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2431 "dgesvd.f"
			ir = 1;
#line 2432 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2436 "dgesvd.f"
			    ldwrkr = *lda;
#line 2437 "dgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2441 "dgesvd.f"
			    ldwrkr = *m;
#line 2442 "dgesvd.f"
			}
#line 2443 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2444 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2449 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2449 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2454 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2456 "dgesvd.f"
			i__2 = *m - 1;
#line 2456 "dgesvd.f"
			i__3 = *m - 1;
#line 2456 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2462 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2462 "dgesvd.f"
			dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2464 "dgesvd.f"
			ie = itau;
#line 2465 "dgesvd.f"
			itauq = ie + *m;
#line 2466 "dgesvd.f"
			itaup = itauq + *m;
#line 2467 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

#line 2472 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2472 "dgesvd.f"
			dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB) */

#line 2481 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2481 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2484 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M + BDSPAC) */

#line 2490 "dgesvd.f"
			dbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2498 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2501 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2505 "dgesvd.f"
			itau = 1;
#line 2506 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 2511 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2511 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2516 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 2521 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2521 "dgesvd.f"
			dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2523 "dgesvd.f"
			ie = itau;
#line 2524 "dgesvd.f"
			itauq = ie + *m;
#line 2525 "dgesvd.f"
			itaup = itauq + *m;
#line 2526 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2530 "dgesvd.f"
			i__2 = *m - 1;
#line 2530 "dgesvd.f"
			i__3 = *m - 1;
#line 2530 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 2536 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2536 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

#line 2543 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2543 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2546 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2552 "dgesvd.f"
			dbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 2556 "dgesvd.f"
		    }

#line 2558 "dgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 2564 "dgesvd.f"
		    i__2 = *m << 2;
#line 2564 "dgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2568 "dgesvd.f"
			iu = 1;
#line 2569 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2573 "dgesvd.f"
			    ldwrku = *lda;
#line 2574 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2575 "dgesvd.f"
			    ldwrkr = *lda;
#line 2576 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2580 "dgesvd.f"
			    ldwrku = *lda;
#line 2581 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2582 "dgesvd.f"
			    ldwrkr = *m;
#line 2583 "dgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2587 "dgesvd.f"
			    ldwrku = *m;
#line 2588 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2589 "dgesvd.f"
			    ldwrkr = *m;
#line 2590 "dgesvd.f"
			}
#line 2591 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2592 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB) */

#line 2597 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2597 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2602 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2604 "dgesvd.f"
			i__2 = *m - 1;
#line 2604 "dgesvd.f"
			i__3 = *m - 1;
#line 2604 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB) */

#line 2610 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2610 "dgesvd.f"
			dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2612 "dgesvd.f"
			ie = itau;
#line 2613 "dgesvd.f"
			itauq = ie + *m;
#line 2614 "dgesvd.f"
			itaup = itauq + *m;
#line 2615 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 2622 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2622 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2626 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M + 4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 2633 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2633 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB) */

#line 2640 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2640 "dgesvd.f"
			dorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2643 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M + BDSPAC) */

#line 2650 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2658 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (Workspace: need M*M) */

#line 2664 "dgesvd.f"
			dlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2667 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2671 "dgesvd.f"
			itau = 1;
#line 2672 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 2677 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2677 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2679 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 2684 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2684 "dgesvd.f"
			dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2686 "dgesvd.f"
			ie = itau;
#line 2687 "dgesvd.f"
			itauq = ie + *m;
#line 2688 "dgesvd.f"
			itaup = itauq + *m;
#line 2689 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2693 "dgesvd.f"
			i__2 = *m - 1;
#line 2693 "dgesvd.f"
			i__3 = *m - 1;
#line 2693 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 2699 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2699 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

#line 2706 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2706 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 2713 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2713 "dgesvd.f"
			dorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2715 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, compute left */
/*                    singular vectors of A in A and compute right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2722 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2726 "dgesvd.f"
		    }

#line 2728 "dgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 2735 "dgesvd.f"
		    i__2 = *m << 2;
#line 2735 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2739 "dgesvd.f"
			iu = 1;
#line 2740 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2744 "dgesvd.f"
			    ldwrku = *lda;
#line 2745 "dgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2749 "dgesvd.f"
			    ldwrku = *m;
#line 2750 "dgesvd.f"
			}
#line 2751 "dgesvd.f"
			itau = iu + ldwrku * *m;
#line 2752 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2757 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2757 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2762 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2764 "dgesvd.f"
			i__2 = *m - 1;
#line 2764 "dgesvd.f"
			i__3 = *m - 1;
#line 2764 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2770 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2770 "dgesvd.f"
			dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2772 "dgesvd.f"
			ie = itau;
#line 2773 "dgesvd.f"
			itauq = ie + *m;
#line 2774 "dgesvd.f"
			itaup = itauq + *m;
#line 2775 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

#line 2780 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2780 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2784 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M + 4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2791 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2791 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB) */

#line 2798 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2798 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2800 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M + BDSPAC) */

#line 2807 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2815 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2818 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2822 "dgesvd.f"
			itau = 1;
#line 2823 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 2828 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2828 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2830 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 2835 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2835 "dgesvd.f"
			dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2840 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2841 "dgesvd.f"
			i__2 = *m - 1;
#line 2841 "dgesvd.f"
			i__3 = *m - 1;
#line 2841 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 2843 "dgesvd.f"
			ie = itau;
#line 2844 "dgesvd.f"
			itauq = ie + *m;
#line 2845 "dgesvd.f"
			itaup = itauq + *m;
#line 2846 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 2851 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2851 "dgesvd.f"
			dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

#line 2859 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2859 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 2866 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2866 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2868 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2875 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2879 "dgesvd.f"
		    }

#line 2881 "dgesvd.f"
		}

#line 2883 "dgesvd.f"
	    } else if (wntva) {

#line 2885 "dgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2891 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 2891 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2895 "dgesvd.f"
			ir = 1;
#line 2896 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2900 "dgesvd.f"
			    ldwrkr = *lda;
#line 2901 "dgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2905 "dgesvd.f"
			    ldwrkr = *m;
#line 2906 "dgesvd.f"
			}
#line 2907 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2908 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 2913 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2913 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2915 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2919 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2921 "dgesvd.f"
			i__2 = *m - 1;
#line 2921 "dgesvd.f"
			i__3 = *m - 1;
#line 2921 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB) */

#line 2927 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2927 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2929 "dgesvd.f"
			ie = itau;
#line 2930 "dgesvd.f"
			itauq = ie + *m;
#line 2931 "dgesvd.f"
			itaup = itauq + *m;
#line 2932 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

#line 2937 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2937 "dgesvd.f"
			dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need M*M + 4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2946 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2946 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2949 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M + BDSPAC) */

#line 2955 "dgesvd.f"
			dbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 2963 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 2968 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 2970 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2974 "dgesvd.f"
			itau = 1;
#line 2975 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 2980 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2980 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2982 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M + N, prefer M + N*NB) */

#line 2987 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2987 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2989 "dgesvd.f"
			ie = itau;
#line 2990 "dgesvd.f"
			itauq = ie + *m;
#line 2991 "dgesvd.f"
			itaup = itauq + *m;
#line 2992 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2996 "dgesvd.f"
			i__2 = *m - 1;
#line 2996 "dgesvd.f"
			i__3 = *m - 1;
#line 2996 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 3002 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3002 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

#line 3010 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3010 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3013 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3019 "dgesvd.f"
			dbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 3023 "dgesvd.f"
		    }

#line 3025 "dgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3031 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3031 "dgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3035 "dgesvd.f"
			iu = 1;
#line 3036 "dgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3040 "dgesvd.f"
			    ldwrku = *lda;
#line 3041 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3042 "dgesvd.f"
			    ldwrkr = *lda;
#line 3043 "dgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3047 "dgesvd.f"
			    ldwrku = *lda;
#line 3048 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3049 "dgesvd.f"
			    ldwrkr = *m;
#line 3050 "dgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3054 "dgesvd.f"
			    ldwrku = *m;
#line 3055 "dgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3056 "dgesvd.f"
			    ldwrkr = *m;
#line 3057 "dgesvd.f"
			}
#line 3058 "dgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3059 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB) */

#line 3064 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3064 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3066 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M*M + M + N, prefer 2*M*M + M + N*NB) */

#line 3071 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3071 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3076 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3078 "dgesvd.f"
			i__2 = *m - 1;
#line 3078 "dgesvd.f"
			i__3 = *m - 1;
#line 3078 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3080 "dgesvd.f"
			ie = itau;
#line 3081 "dgesvd.f"
			itauq = ie + *m;
#line 3082 "dgesvd.f"
			itaup = itauq + *m;
#line 3083 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 3090 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3090 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3094 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M + 4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 3101 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3101 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB) */

#line 3108 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3108 "dgesvd.f"
			dorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3111 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M + BDSPAC) */

#line 3118 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3126 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3131 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3135 "dgesvd.f"
			dlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3138 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3142 "dgesvd.f"
			itau = 1;
#line 3143 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 3148 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3148 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3150 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M + N, prefer M + N*NB) */

#line 3155 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3155 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3157 "dgesvd.f"
			ie = itau;
#line 3158 "dgesvd.f"
			itauq = ie + *m;
#line 3159 "dgesvd.f"
			itaup = itauq + *m;
#line 3160 "dgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3164 "dgesvd.f"
			i__2 = *m - 1;
#line 3164 "dgesvd.f"
			i__3 = *m - 1;
#line 3164 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 3170 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3170 "dgesvd.f"
			dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

#line 3178 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3178 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 3185 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3185 "dgesvd.f"
			dorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3187 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3194 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3198 "dgesvd.f"
		    }

#line 3200 "dgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3207 "dgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3207 "dgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3211 "dgesvd.f"
			iu = 1;
#line 3212 "dgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3216 "dgesvd.f"
			    ldwrku = *lda;
#line 3217 "dgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3221 "dgesvd.f"
			    ldwrku = *m;
#line 3222 "dgesvd.f"
			}
#line 3223 "dgesvd.f"
			itau = iu + ldwrku * *m;
#line 3224 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

#line 3229 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3229 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3231 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB) */

#line 3236 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3236 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3241 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3243 "dgesvd.f"
			i__2 = *m - 1;
#line 3243 "dgesvd.f"
			i__3 = *m - 1;
#line 3243 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3245 "dgesvd.f"
			ie = itau;
#line 3246 "dgesvd.f"
			itauq = ie + *m;
#line 3247 "dgesvd.f"
			itaup = itauq + *m;
#line 3248 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

#line 3253 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3253 "dgesvd.f"
			dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3257 "dgesvd.f"
			dlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB) */

#line 3263 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3263 "dgesvd.f"
			dorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB) */

#line 3270 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3270 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3272 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M + BDSPAC) */

#line 3279 "dgesvd.f"
			dbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3287 "dgesvd.f"
			dgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3292 "dgesvd.f"
			dlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3294 "dgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3298 "dgesvd.f"
			itau = 1;
#line 3299 "dgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

#line 3304 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3304 "dgesvd.f"
			dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3306 "dgesvd.f"
			dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M + N, prefer M + N*NB) */

#line 3311 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3311 "dgesvd.f"
			dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3316 "dgesvd.f"
			dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3317 "dgesvd.f"
			i__2 = *m - 1;
#line 3317 "dgesvd.f"
			i__3 = *m - 1;
#line 3317 "dgesvd.f"
			dlaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 3319 "dgesvd.f"
			ie = itau;
#line 3320 "dgesvd.f"
			itauq = ie + *m;
#line 3321 "dgesvd.f"
			itaup = itauq + *m;
#line 3322 "dgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

#line 3327 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3327 "dgesvd.f"
			dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

#line 3335 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3335 "dgesvd.f"
			dormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 3342 "dgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3342 "dgesvd.f"
			dorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3344 "dgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3351 "dgesvd.f"
			dbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3355 "dgesvd.f"
		    }

#line 3357 "dgesvd.f"
		}

#line 3359 "dgesvd.f"
	    }

#line 3361 "dgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3368 "dgesvd.f"
	    ie = 1;
#line 3369 "dgesvd.f"
	    itauq = ie + *m;
#line 3370 "dgesvd.f"
	    itaup = itauq + *m;
#line 3371 "dgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB) */

#line 3376 "dgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3376 "dgesvd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 3379 "dgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB) */

#line 3385 "dgesvd.f"
		dlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3386 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3386 "dgesvd.f"
		dorgbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3388 "dgesvd.f"
	    }
#line 3389 "dgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 3*M + NRVT, prefer 3*M + NRVT*NB) */

#line 3395 "dgesvd.f"
		dlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3396 "dgesvd.f"
		if (wntva) {
#line 3396 "dgesvd.f"
		    nrvt = *n;
#line 3396 "dgesvd.f"
		}
#line 3398 "dgesvd.f"
		if (wntvs) {
#line 3398 "dgesvd.f"
		    nrvt = *m;
#line 3398 "dgesvd.f"
		}
#line 3400 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3400 "dgesvd.f"
		dorgbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3402 "dgesvd.f"
	    }
#line 3403 "dgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB) */

#line 3409 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3409 "dgesvd.f"
		dorgbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3411 "dgesvd.f"
	    }
#line 3412 "dgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M, prefer 3*M + M*NB) */

#line 3418 "dgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3418 "dgesvd.f"
		dorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3420 "dgesvd.f"
	    }
#line 3421 "dgesvd.f"
	    iwork = ie + *m;
#line 3422 "dgesvd.f"
	    if (wntuas || wntuo) {
#line 3422 "dgesvd.f"
		nru = *m;
#line 3422 "dgesvd.f"
	    }
#line 3424 "dgesvd.f"
	    if (wntun) {
#line 3424 "dgesvd.f"
		nru = 0;
#line 3424 "dgesvd.f"
	    }
#line 3426 "dgesvd.f"
	    if (wntvas || wntvo) {
#line 3426 "dgesvd.f"
		ncvt = *n;
#line 3426 "dgesvd.f"
	    }
#line 3428 "dgesvd.f"
	    if (wntvn) {
#line 3428 "dgesvd.f"
		ncvt = 0;
#line 3428 "dgesvd.f"
	    }
#line 3430 "dgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3437 "dgesvd.f"
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3439 "dgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 3446 "dgesvd.f"
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 3448 "dgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3455 "dgesvd.f"
		dbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3457 "dgesvd.f"
	    }

#line 3459 "dgesvd.f"
	}

#line 3461 "dgesvd.f"
    }

/*     If DBDSQR failed to converge, copy unconverged superdiagonals */
/*     to WORK( 2:MINMN ) */

#line 3466 "dgesvd.f"
    if (*info != 0) {
#line 3467 "dgesvd.f"
	if (ie > 2) {
#line 3468 "dgesvd.f"
	    i__2 = minmn - 1;
#line 3468 "dgesvd.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 3469 "dgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3470 "dgesvd.f"
/* L50: */
#line 3470 "dgesvd.f"
	    }
#line 3471 "dgesvd.f"
	}
#line 3472 "dgesvd.f"
	if (ie < 2) {
#line 3473 "dgesvd.f"
	    for (i__ = minmn - 1; i__ >= 1; --i__) {
#line 3474 "dgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3475 "dgesvd.f"
/* L60: */
#line 3475 "dgesvd.f"
	    }
#line 3476 "dgesvd.f"
	}
#line 3477 "dgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3481 "dgesvd.f"
    if (iscl == 1) {
#line 3482 "dgesvd.f"
	if (anrm > bignum) {
#line 3482 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3482 "dgesvd.f"
	}
#line 3485 "dgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3485 "dgesvd.f"
	    i__2 = minmn - 1;
#line 3485 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3485 "dgesvd.f"
	}
#line 3488 "dgesvd.f"
	if (anrm < smlnum) {
#line 3488 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3488 "dgesvd.f"
	}
#line 3491 "dgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3491 "dgesvd.f"
	    i__2 = minmn - 1;
#line 3491 "dgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3491 "dgesvd.f"
	}
#line 3494 "dgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3498 "dgesvd.f"
    work[1] = (doublereal) maxwrk;

#line 3500 "dgesvd.f"
    return 0;

/*     End of DGESVD */

} /* dgesvd_ */

