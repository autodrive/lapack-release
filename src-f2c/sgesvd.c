#line 1 "sgesvd.f"
/* sgesvd.f -- translated by f2c (version 20100827).
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

#line 1 "sgesvd.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c_n1 = -1;
static doublereal c_b57 = 0.;
static integer c__1 = 1;
static doublereal c_b79 = 1.;

/* > \brief <b> SGESVD computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGESVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU, JOBVT */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), S( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGESVD computes the singular value decomposition (SVD) of a real */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          S is REAL array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array, dimension (LDU,UCOL) */
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
/* >          VT is REAL array, dimension (LDVT,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* >          > 0:  if SBDSQR did not converge, INFO specifies how many */
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

/* > \ingroup realGEsing */

/*  ===================================================================== */
/* Subroutine */ int sgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
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
    static integer ierr, itau, ncvt, nrvt, lwork_sgebrd__, lwork_sgelqf__, 
	    lwork_sgeqrf__;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer minmn, wrkbl, itaup, itauq, mnthr, iwork;
    static logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    static integer bdspac;
    extern /* Subroutine */ int sgebrd_(integer *, integer *, doublereal *, 
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
	    sbdsqr_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), sorgbr_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), sormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer ldwrkr, minwrk, ldwrku, maxwrk;
    extern /* Subroutine */ int sorglq_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal smlnum;
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical lquery, wntuas, wntvas;
    static integer lwork_sorgbr_p__, lwork_sorgbr_q__, lwork_sorglq_m__, 
	    lwork_sorglq_n__, lwork_sorgqr_m__, lwork_sorgqr_n__;


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

#line 267 "sgesvd.f"
    /* Parameter adjustments */
#line 267 "sgesvd.f"
    a_dim1 = *lda;
#line 267 "sgesvd.f"
    a_offset = 1 + a_dim1;
#line 267 "sgesvd.f"
    a -= a_offset;
#line 267 "sgesvd.f"
    --s;
#line 267 "sgesvd.f"
    u_dim1 = *ldu;
#line 267 "sgesvd.f"
    u_offset = 1 + u_dim1;
#line 267 "sgesvd.f"
    u -= u_offset;
#line 267 "sgesvd.f"
    vt_dim1 = *ldvt;
#line 267 "sgesvd.f"
    vt_offset = 1 + vt_dim1;
#line 267 "sgesvd.f"
    vt -= vt_offset;
#line 267 "sgesvd.f"
    --work;
#line 267 "sgesvd.f"

#line 267 "sgesvd.f"
    /* Function Body */
#line 267 "sgesvd.f"
    *info = 0;
#line 268 "sgesvd.f"
    minmn = min(*m,*n);
#line 269 "sgesvd.f"
    wntua = lsame_(jobu, "A", (ftnlen)1, (ftnlen)1);
#line 270 "sgesvd.f"
    wntus = lsame_(jobu, "S", (ftnlen)1, (ftnlen)1);
#line 271 "sgesvd.f"
    wntuas = wntua || wntus;
#line 272 "sgesvd.f"
    wntuo = lsame_(jobu, "O", (ftnlen)1, (ftnlen)1);
#line 273 "sgesvd.f"
    wntun = lsame_(jobu, "N", (ftnlen)1, (ftnlen)1);
#line 274 "sgesvd.f"
    wntva = lsame_(jobvt, "A", (ftnlen)1, (ftnlen)1);
#line 275 "sgesvd.f"
    wntvs = lsame_(jobvt, "S", (ftnlen)1, (ftnlen)1);
#line 276 "sgesvd.f"
    wntvas = wntva || wntvs;
#line 277 "sgesvd.f"
    wntvo = lsame_(jobvt, "O", (ftnlen)1, (ftnlen)1);
#line 278 "sgesvd.f"
    wntvn = lsame_(jobvt, "N", (ftnlen)1, (ftnlen)1);
#line 279 "sgesvd.f"
    lquery = *lwork == -1;

#line 281 "sgesvd.f"
    if (! (wntua || wntus || wntuo || wntun)) {
#line 282 "sgesvd.f"
	*info = -1;
#line 283 "sgesvd.f"
    } else if (! (wntva || wntvs || wntvo || wntvn) || wntvo && wntuo) {
#line 285 "sgesvd.f"
	*info = -2;
#line 286 "sgesvd.f"
    } else if (*m < 0) {
#line 287 "sgesvd.f"
	*info = -3;
#line 288 "sgesvd.f"
    } else if (*n < 0) {
#line 289 "sgesvd.f"
	*info = -4;
#line 290 "sgesvd.f"
    } else if (*lda < max(1,*m)) {
#line 291 "sgesvd.f"
	*info = -6;
#line 292 "sgesvd.f"
    } else if (*ldu < 1 || wntuas && *ldu < *m) {
#line 293 "sgesvd.f"
	*info = -9;
#line 294 "sgesvd.f"
    } else if (*ldvt < 1 || wntva && *ldvt < *n || wntvs && *ldvt < minmn) {
#line 296 "sgesvd.f"
	*info = -11;
#line 297 "sgesvd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 306 "sgesvd.f"
    if (*info == 0) {
#line 307 "sgesvd.f"
	minwrk = 1;
#line 308 "sgesvd.f"
	maxwrk = 1;
#line 309 "sgesvd.f"
	if (*m >= *n && minmn > 0) {

/*           Compute space needed for SBDSQR */

/* Writing concatenation */
#line 313 "sgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 313 "sgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 313 "sgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 313 "sgesvd.f"
	    mnthr = ilaenv_(&c__6, "SGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
#line 314 "sgesvd.f"
	    bdspac = *n * 5;
/*           Compute space needed for SGEQRF */
#line 316 "sgesvd.f"
	    sgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 317 "sgesvd.f"
	    lwork_sgeqrf__ = (integer) dum[0];
/*           Compute space needed for SORGQR */
#line 319 "sgesvd.f"
	    sorgqr_(m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 320 "sgesvd.f"
	    lwork_sorgqr_n__ = (integer) dum[0];
#line 321 "sgesvd.f"
	    sorgqr_(m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 322 "sgesvd.f"
	    lwork_sorgqr_m__ = (integer) dum[0];
/*           Compute space needed for SGEBRD */
#line 324 "sgesvd.f"
	    sgebrd_(n, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 326 "sgesvd.f"
	    lwork_sgebrd__ = (integer) dum[0];
/*           Compute space needed for SORGBR P */
#line 328 "sgesvd.f"
	    sorgbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 330 "sgesvd.f"
	    lwork_sorgbr_p__ = (integer) dum[0];
/*           Compute space needed for SORGBR Q */
#line 332 "sgesvd.f"
	    sorgbr_("Q", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 334 "sgesvd.f"
	    lwork_sorgbr_q__ = (integer) dum[0];

#line 336 "sgesvd.f"
	    if (*m >= mnthr) {
#line 337 "sgesvd.f"
		if (wntun) {

/*                 Path 1 (M much larger than N, JOBU='N') */

#line 341 "sgesvd.f"
		    maxwrk = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 342 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_sgebrd__;
#line 342 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 343 "sgesvd.f"
		    if (wntvo || wntvas) {
/* Computing MAX */
#line 343 "sgesvd.f"
			i__2 = maxwrk, i__3 = *n * 3 + lwork_sorgbr_p__;
#line 343 "sgesvd.f"
			maxwrk = max(i__2,i__3);
#line 343 "sgesvd.f"
		    }
#line 345 "sgesvd.f"
		    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 346 "sgesvd.f"
		    i__2 = *n << 2;
#line 346 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 347 "sgesvd.f"
		} else if (wntuo && wntvn) {

/*                 Path 2 (M much larger than N, JOBU='O', JOBVT='N') */

#line 351 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 352 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_n__;
#line 352 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 353 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 353 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 354 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 354 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 355 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 356 "sgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
#line 356 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 357 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 357 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 358 "sgesvd.f"
		} else if (wntuo && wntvas) {

/*                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or */
/*                 'A') */

#line 363 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 364 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_n__;
#line 364 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 365 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 365 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 366 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 366 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 367 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_p__;
#line 367 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 368 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 369 "sgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
#line 369 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 370 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 370 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 371 "sgesvd.f"
		} else if (wntus && wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */

#line 375 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 376 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_n__;
#line 376 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 377 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 377 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 378 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 378 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 379 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 380 "sgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 381 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 381 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 382 "sgesvd.f"
		} else if (wntus && wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */

#line 386 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 387 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_n__;
#line 387 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 388 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 388 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 389 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 389 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 390 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_p__;
#line 390 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 391 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 392 "sgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
#line 393 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 393 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 394 "sgesvd.f"
		} else if (wntus && wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or */
/*                 'A') */

#line 399 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 400 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_n__;
#line 400 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 401 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 401 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 402 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 402 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 403 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_p__;
#line 403 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 404 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 405 "sgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 406 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 406 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 407 "sgesvd.f"
		} else if (wntua && wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */

#line 411 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 412 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_m__;
#line 412 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 413 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 413 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 414 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 414 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 415 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 416 "sgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 417 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 417 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 418 "sgesvd.f"
		} else if (wntua && wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */

#line 422 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 423 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_m__;
#line 423 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 424 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 424 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 425 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 425 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 426 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_p__;
#line 426 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 427 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 428 "sgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
#line 429 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 429 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 430 "sgesvd.f"
		} else if (wntua && wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or */
/*                 'A') */

#line 435 "sgesvd.f"
		    wrkbl = *n + lwork_sgeqrf__;
/* Computing MAX */
#line 436 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_sorgqr_m__;
#line 436 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 437 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sgebrd__;
#line 437 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 438 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 438 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 439 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *n * 3 + lwork_sorgbr_p__;
#line 439 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 440 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 441 "sgesvd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 442 "sgesvd.f"
		    i__2 = *n * 3 + *m;
#line 442 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 443 "sgesvd.f"
		}
#line 444 "sgesvd.f"
	    } else {

/*              Path 10 (M at least N, but not much larger) */

#line 448 "sgesvd.f"
		sgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 450 "sgesvd.f"
		lwork_sgebrd__ = (integer) dum[0];
#line 451 "sgesvd.f"
		maxwrk = *n * 3 + lwork_sgebrd__;
#line 452 "sgesvd.f"
		if (wntus || wntuo) {
#line 453 "sgesvd.f"
		    sorgbr_("Q", m, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 455 "sgesvd.f"
		    lwork_sorgbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 456 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 456 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 457 "sgesvd.f"
		}
#line 458 "sgesvd.f"
		if (wntua) {
#line 459 "sgesvd.f"
		    sorgbr_("Q", m, m, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 461 "sgesvd.f"
		    lwork_sorgbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 462 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_sorgbr_q__;
#line 462 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 463 "sgesvd.f"
		}
#line 464 "sgesvd.f"
		if (! wntvn) {
/* Computing MAX */
#line 465 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *n * 3 + lwork_sorgbr_p__;
#line 465 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 466 "sgesvd.f"
		}
#line 467 "sgesvd.f"
		maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 468 "sgesvd.f"
		i__2 = *n * 3 + *m;
#line 468 "sgesvd.f"
		minwrk = max(i__2,bdspac);
#line 469 "sgesvd.f"
	    }
#line 470 "sgesvd.f"
	} else if (minmn > 0) {

/*           Compute space needed for SBDSQR */

/* Writing concatenation */
#line 474 "sgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 474 "sgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 474 "sgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 474 "sgesvd.f"
	    mnthr = ilaenv_(&c__6, "SGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
#line 475 "sgesvd.f"
	    bdspac = *m * 5;
/*           Compute space needed for SGELQF */
#line 477 "sgesvd.f"
	    sgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 478 "sgesvd.f"
	    lwork_sgelqf__ = (integer) dum[0];
/*           Compute space needed for SORGLQ */
#line 480 "sgesvd.f"
	    sorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
#line 481 "sgesvd.f"
	    lwork_sorglq_n__ = (integer) dum[0];
#line 482 "sgesvd.f"
	    sorglq_(m, n, m, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 483 "sgesvd.f"
	    lwork_sorglq_m__ = (integer) dum[0];
/*           Compute space needed for SGEBRD */
#line 485 "sgesvd.f"
	    sgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 487 "sgesvd.f"
	    lwork_sgebrd__ = (integer) dum[0];
/*            Compute space needed for SORGBR P */
#line 489 "sgesvd.f"
	    sorgbr_("P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 491 "sgesvd.f"
	    lwork_sorgbr_p__ = (integer) dum[0];
/*           Compute space needed for SORGBR Q */
#line 493 "sgesvd.f"
	    sorgbr_("Q", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 495 "sgesvd.f"
	    lwork_sorgbr_q__ = (integer) dum[0];
#line 496 "sgesvd.f"
	    if (*n >= mnthr) {
#line 497 "sgesvd.f"
		if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N') */

#line 501 "sgesvd.f"
		    maxwrk = *m + lwork_sgelqf__;
/* Computing MAX */
#line 502 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_sgebrd__;
#line 502 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 503 "sgesvd.f"
		    if (wntuo || wntuas) {
/* Computing MAX */
#line 503 "sgesvd.f"
			i__2 = maxwrk, i__3 = *m * 3 + lwork_sorgbr_q__;
#line 503 "sgesvd.f"
			maxwrk = max(i__2,i__3);
#line 503 "sgesvd.f"
		    }
#line 505 "sgesvd.f"
		    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 506 "sgesvd.f"
		    i__2 = *m << 2;
#line 506 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 507 "sgesvd.f"
		} else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O') */

#line 511 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 512 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_m__;
#line 512 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 513 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 513 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 514 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 514 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 515 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 516 "sgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
#line 516 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 517 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 517 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 518 "sgesvd.f"
		} else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='O') */

#line 523 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 524 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_m__;
#line 524 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 525 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 525 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 526 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 526 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 527 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_q__;
#line 527 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 528 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
#line 529 "sgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
#line 529 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
/* Computing MAX */
#line 530 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 530 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 531 "sgesvd.f"
		} else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */

#line 535 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 536 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_m__;
#line 536 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 537 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 537 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 538 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 538 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 539 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 540 "sgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 541 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 541 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 542 "sgesvd.f"
		} else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */

#line 546 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 547 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_m__;
#line 547 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 548 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 548 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 549 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 549 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 550 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_q__;
#line 550 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 551 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 552 "sgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
#line 553 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 553 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 554 "sgesvd.f"
		    maxwrk = max(maxwrk,minwrk);
#line 555 "sgesvd.f"
		} else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='S') */

#line 560 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 561 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_m__;
#line 561 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 562 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 562 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 563 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 563 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 564 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_q__;
#line 564 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 565 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 566 "sgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 567 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 567 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 568 "sgesvd.f"
		} else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */

#line 572 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 573 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_n__;
#line 573 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 574 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 574 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 575 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 575 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 576 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 577 "sgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 578 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 578 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 579 "sgesvd.f"
		} else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */

#line 583 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 584 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_n__;
#line 584 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 585 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 585 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 586 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 586 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 587 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_q__;
#line 587 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 588 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 589 "sgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
#line 590 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 590 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 591 "sgesvd.f"
		} else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='A') */

#line 596 "sgesvd.f"
		    wrkbl = *m + lwork_sgelqf__;
/* Computing MAX */
#line 597 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_sorglq_n__;
#line 597 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 598 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sgebrd__;
#line 598 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 599 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 599 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 600 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *m * 3 + lwork_sorgbr_q__;
#line 600 "sgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 601 "sgesvd.f"
		    wrkbl = max(wrkbl,bdspac);
#line 602 "sgesvd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 603 "sgesvd.f"
		    i__2 = *m * 3 + *n;
#line 603 "sgesvd.f"
		    minwrk = max(i__2,bdspac);
#line 604 "sgesvd.f"
		}
#line 605 "sgesvd.f"
	    } else {

/*              Path 10t(N greater than M, but not much larger) */

#line 609 "sgesvd.f"
		sgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 611 "sgesvd.f"
		lwork_sgebrd__ = (integer) dum[0];
#line 612 "sgesvd.f"
		maxwrk = *m * 3 + lwork_sgebrd__;
#line 613 "sgesvd.f"
		if (wntvs || wntvo) {
/*                Compute space needed for SORGBR P */
#line 615 "sgesvd.f"
		    sorgbr_("P", m, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 617 "sgesvd.f"
		    lwork_sorgbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 618 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 618 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 619 "sgesvd.f"
		}
#line 620 "sgesvd.f"
		if (wntva) {
#line 621 "sgesvd.f"
		    sorgbr_("P", n, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 623 "sgesvd.f"
		    lwork_sorgbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 624 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_sorgbr_p__;
#line 624 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 625 "sgesvd.f"
		}
#line 626 "sgesvd.f"
		if (! wntun) {
/* Computing MAX */
#line 627 "sgesvd.f"
		    i__2 = maxwrk, i__3 = *m * 3 + lwork_sorgbr_q__;
#line 627 "sgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 628 "sgesvd.f"
		}
#line 629 "sgesvd.f"
		maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 630 "sgesvd.f"
		i__2 = *m * 3 + *n;
#line 630 "sgesvd.f"
		minwrk = max(i__2,bdspac);
#line 631 "sgesvd.f"
	    }
#line 632 "sgesvd.f"
	}
#line 633 "sgesvd.f"
	maxwrk = max(maxwrk,minwrk);
#line 634 "sgesvd.f"
	work[1] = (doublereal) maxwrk;

#line 636 "sgesvd.f"
	if (*lwork < minwrk && ! lquery) {
#line 637 "sgesvd.f"
	    *info = -13;
#line 638 "sgesvd.f"
	}
#line 639 "sgesvd.f"
    }

#line 641 "sgesvd.f"
    if (*info != 0) {
#line 642 "sgesvd.f"
	i__2 = -(*info);
#line 642 "sgesvd.f"
	xerbla_("SGESVD", &i__2, (ftnlen)6);
#line 643 "sgesvd.f"
	return 0;
#line 644 "sgesvd.f"
    } else if (lquery) {
#line 645 "sgesvd.f"
	return 0;
#line 646 "sgesvd.f"
    }

/*     Quick return if possible */

#line 650 "sgesvd.f"
    if (*m == 0 || *n == 0) {
#line 651 "sgesvd.f"
	return 0;
#line 652 "sgesvd.f"
    }

/*     Get machine constants */

#line 656 "sgesvd.f"
    eps = slamch_("P", (ftnlen)1);
#line 657 "sgesvd.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 658 "sgesvd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 662 "sgesvd.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 663 "sgesvd.f"
    iscl = 0;
#line 664 "sgesvd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 665 "sgesvd.f"
	iscl = 1;
#line 666 "sgesvd.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 667 "sgesvd.f"
    } else if (anrm > bignum) {
#line 668 "sgesvd.f"
	iscl = 1;
#line 669 "sgesvd.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 670 "sgesvd.f"
    }

#line 672 "sgesvd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 678 "sgesvd.f"
	if (*m >= mnthr) {

#line 680 "sgesvd.f"
	    if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N') */
/*              No left singular vectors to be computed */

#line 685 "sgesvd.f"
		itau = 1;
#line 686 "sgesvd.f"
		iwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need 2*N, prefer N+N*NB) */

#line 691 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 691 "sgesvd.f"
		sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out below R */

#line 696 "sgesvd.f"
		i__2 = *n - 1;
#line 696 "sgesvd.f"
		i__3 = *n - 1;
#line 696 "sgesvd.f"
		slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 697 "sgesvd.f"
		ie = 1;
#line 698 "sgesvd.f"
		itauq = ie + *n;
#line 699 "sgesvd.f"
		itaup = itauq + *n;
#line 700 "sgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 705 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 705 "sgesvd.f"
		sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 708 "sgesvd.f"
		ncvt = 0;
#line 709 "sgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 714 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 714 "sgesvd.f"
		    sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 716 "sgesvd.f"
		    ncvt = *n;
#line 717 "sgesvd.f"
		}
#line 718 "sgesvd.f"
		iwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 724 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, dum, &c__1, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 729 "sgesvd.f"
		if (wntvas) {
#line 729 "sgesvd.f"
		    slacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 729 "sgesvd.f"
		}

#line 732 "sgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

/* Computing MAX */
#line 738 "sgesvd.f"
		i__2 = *n << 2;
#line 738 "sgesvd.f"
		if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 742 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 743 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 743 "sgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 747 "sgesvd.f"
			ldwrku = *lda;
#line 748 "sgesvd.f"
			ldwrkr = *lda;
#line 749 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 749 "sgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 749 "sgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 753 "sgesvd.f"
			    ldwrku = *lda;
#line 754 "sgesvd.f"
			    ldwrkr = *n;
#line 755 "sgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 759 "sgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 760 "sgesvd.f"
			    ldwrkr = *n;
#line 761 "sgesvd.f"
			}
#line 761 "sgesvd.f"
		    }
#line 762 "sgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 763 "sgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 768 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 768 "sgesvd.f"
		    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 773 "sgesvd.f"
		    slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 774 "sgesvd.f"
		    i__2 = *n - 1;
#line 774 "sgesvd.f"
		    i__3 = *n - 1;
#line 774 "sgesvd.f"
		    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1], 
			    &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 780 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 780 "sgesvd.f"
		    sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 782 "sgesvd.f"
		    ie = itau;
#line 783 "sgesvd.f"
		    itauq = ie + *n;
#line 784 "sgesvd.f"
		    itaup = itauq + *n;
#line 785 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 790 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 790 "sgesvd.f"
		    sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 797 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 797 "sgesvd.f"
		    sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 800 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (Workspace: need N*N+BDSPAC) */

#line 806 "sgesvd.f"
		    sbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &
			    c__1, &work[ir], &ldwrkr, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 809 "sgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N+2*N, prefer N*N+M*N+N) */

#line 815 "sgesvd.f"
		    i__2 = *m;
#line 815 "sgesvd.f"
		    i__3 = ldwrku;
#line 815 "sgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 816 "sgesvd.f"
			i__4 = *m - i__ + 1;
#line 816 "sgesvd.f"
			chunk = min(i__4,ldwrku);
#line 817 "sgesvd.f"
			sgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 820 "sgesvd.f"
			slacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 822 "sgesvd.f"
/* L10: */
#line 822 "sgesvd.f"
		    }

#line 824 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 828 "sgesvd.f"
		    ie = 1;
#line 829 "sgesvd.f"
		    itauq = ie + *n;
#line 830 "sgesvd.f"
		    itaup = itauq + *n;
#line 831 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 836 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 836 "sgesvd.f"
		    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (Workspace: need 4*N, prefer 3*N+N*NB) */

#line 843 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 843 "sgesvd.f"
		    sorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 845 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 851 "sgesvd.f"
		    sbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &
			    c__1, &a[a_offset], lda, dum, &c__1, &work[iwork],
			     info, (ftnlen)1);

#line 854 "sgesvd.f"
		}

#line 856 "sgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

/* Computing MAX */
#line 862 "sgesvd.f"
		i__3 = *n << 2;
#line 862 "sgesvd.f"
		if (*lwork >= *n * *n + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 866 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 867 "sgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 867 "sgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 871 "sgesvd.f"
			ldwrku = *lda;
#line 872 "sgesvd.f"
			ldwrkr = *lda;
#line 873 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 873 "sgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 873 "sgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 877 "sgesvd.f"
			    ldwrku = *lda;
#line 878 "sgesvd.f"
			    ldwrkr = *n;
#line 879 "sgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 883 "sgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 884 "sgesvd.f"
			    ldwrkr = *n;
#line 885 "sgesvd.f"
			}
#line 885 "sgesvd.f"
		    }
#line 886 "sgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 887 "sgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 892 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 892 "sgesvd.f"
		    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 897 "sgesvd.f"
		    slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 898 "sgesvd.f"
		    if (*n > 1) {
#line 898 "sgesvd.f"
			i__3 = *n - 1;
#line 898 "sgesvd.f"
			i__2 = *n - 1;
#line 898 "sgesvd.f"
			slaset_("L", &i__3, &i__2, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 898 "sgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 905 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 905 "sgesvd.f"
		    sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 907 "sgesvd.f"
		    ie = itau;
#line 908 "sgesvd.f"
		    itauq = ie + *n;
#line 909 "sgesvd.f"
		    itaup = itauq + *n;
#line 910 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 915 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 915 "sgesvd.f"
		    sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 918 "sgesvd.f"
		    slacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 923 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 923 "sgesvd.f"
		    sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need N*N+4*N-1, prefer N*N+3*N+(N-1)*NB) */

#line 930 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 930 "sgesvd.f"
		    sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 932 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (Workspace: need N*N+BDSPAC) */

#line 939 "sgesvd.f"
		    sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, dum, &c__1, 
			    &work[iwork], info, (ftnlen)1);
#line 942 "sgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N+2*N, prefer N*N+M*N+N) */

#line 948 "sgesvd.f"
		    i__3 = *m;
#line 948 "sgesvd.f"
		    i__2 = ldwrku;
#line 948 "sgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 949 "sgesvd.f"
			i__4 = *m - i__ + 1;
#line 949 "sgesvd.f"
			chunk = min(i__4,ldwrku);
#line 950 "sgesvd.f"
			sgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 953 "sgesvd.f"
			slacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 955 "sgesvd.f"
/* L20: */
#line 955 "sgesvd.f"
		    }

#line 957 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 961 "sgesvd.f"
		    itau = 1;
#line 962 "sgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need 2*N, prefer N+N*NB) */

#line 967 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 967 "sgesvd.f"
		    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 972 "sgesvd.f"
		    slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 973 "sgesvd.f"
		    if (*n > 1) {
#line 973 "sgesvd.f"
			i__2 = *n - 1;
#line 973 "sgesvd.f"
			i__3 = *n - 1;
#line 973 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 973 "sgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need 2*N, prefer N+N*NB) */

#line 980 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 980 "sgesvd.f"
		    sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 982 "sgesvd.f"
		    ie = itau;
#line 983 "sgesvd.f"
		    itauq = ie + *n;
#line 984 "sgesvd.f"
		    itaup = itauq + *n;
#line 985 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 990 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 990 "sgesvd.f"
		    sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 997 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 997 "sgesvd.f"
		    sormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1004 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1004 "sgesvd.f"
		    sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1006 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (Workspace: need BDSPAC) */

#line 1013 "sgesvd.f"
		    sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 1016 "sgesvd.f"
		}

#line 1018 "sgesvd.f"
	    } else if (wntus) {

#line 1020 "sgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1026 "sgesvd.f"
		    i__2 = *n << 2;
#line 1026 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1030 "sgesvd.f"
			ir = 1;
#line 1031 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1035 "sgesvd.f"
			    ldwrkr = *lda;
#line 1036 "sgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1040 "sgesvd.f"
			    ldwrkr = *n;
#line 1041 "sgesvd.f"
			}
#line 1042 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1043 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1048 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1048 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1053 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1055 "sgesvd.f"
			i__2 = *n - 1;
#line 1055 "sgesvd.f"
			i__3 = *n - 1;
#line 1055 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1061 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1061 "sgesvd.f"
			sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1063 "sgesvd.f"
			ie = itau;
#line 1064 "sgesvd.f"
			itauq = ie + *n;
#line 1065 "sgesvd.f"
			itaup = itauq + *n;
#line 1066 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1071 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1071 "sgesvd.f"
			sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1079 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1079 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1082 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1088 "sgesvd.f"
			sbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (Workspace: need N*N) */

#line 1096 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1099 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1103 "sgesvd.f"
			itau = 1;
#line 1104 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1109 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1109 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1111 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1116 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1116 "sgesvd.f"
			sorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1118 "sgesvd.f"
			ie = itau;
#line 1119 "sgesvd.f"
			itauq = ie + *n;
#line 1120 "sgesvd.f"
			itaup = itauq + *n;
#line 1121 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1125 "sgesvd.f"
			i__2 = *n - 1;
#line 1125 "sgesvd.f"
			i__3 = *n - 1;
#line 1125 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1131 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1131 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1138 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1138 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1141 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1147 "sgesvd.f"
			sbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1151 "sgesvd.f"
		    }

#line 1153 "sgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1159 "sgesvd.f"
		    i__2 = *n << 2;
#line 1159 "sgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1163 "sgesvd.f"
			iu = 1;
#line 1164 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1168 "sgesvd.f"
			    ldwrku = *lda;
#line 1169 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1170 "sgesvd.f"
			    ldwrkr = *lda;
#line 1171 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1175 "sgesvd.f"
			    ldwrku = *lda;
#line 1176 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1177 "sgesvd.f"
			    ldwrkr = *n;
#line 1178 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1182 "sgesvd.f"
			    ldwrku = *n;
#line 1183 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1184 "sgesvd.f"
			    ldwrkr = *n;
#line 1185 "sgesvd.f"
			}
#line 1186 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1187 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1192 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1192 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1197 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1199 "sgesvd.f"
			i__2 = *n - 1;
#line 1199 "sgesvd.f"
			i__3 = *n - 1;
#line 1199 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1205 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1205 "sgesvd.f"
			sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1207 "sgesvd.f"
			ie = itau;
#line 1208 "sgesvd.f"
			itauq = ie + *n;
#line 1209 "sgesvd.f"
			itaup = itauq + *n;
#line 1210 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1217 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1217 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1221 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB) */

#line 1227 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1227 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1235 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1235 "sgesvd.f"
			sorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1238 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N+BDSPAC) */

#line 1245 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1253 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (Workspace: need N*N) */

#line 1259 "sgesvd.f"
			slacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1262 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1266 "sgesvd.f"
			itau = 1;
#line 1267 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1272 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1272 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1274 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1279 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1279 "sgesvd.f"
			sorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1281 "sgesvd.f"
			ie = itau;
#line 1282 "sgesvd.f"
			itauq = ie + *n;
#line 1283 "sgesvd.f"
			itaup = itauq + *n;
#line 1284 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1288 "sgesvd.f"
			i__2 = *n - 1;
#line 1288 "sgesvd.f"
			i__3 = *n - 1;
#line 1288 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1294 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1294 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1301 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1301 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1308 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1308 "sgesvd.f"
			sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1310 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1317 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1321 "sgesvd.f"
		    }

#line 1323 "sgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1330 "sgesvd.f"
		    i__2 = *n << 2;
#line 1330 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1334 "sgesvd.f"
			iu = 1;
#line 1335 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1339 "sgesvd.f"
			    ldwrku = *lda;
#line 1340 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1344 "sgesvd.f"
			    ldwrku = *n;
#line 1345 "sgesvd.f"
			}
#line 1346 "sgesvd.f"
			itau = iu + ldwrku * *n;
#line 1347 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1352 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1352 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1357 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1359 "sgesvd.f"
			i__2 = *n - 1;
#line 1359 "sgesvd.f"
			i__3 = *n - 1;
#line 1359 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1365 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1365 "sgesvd.f"
			sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1367 "sgesvd.f"
			ie = itau;
#line 1368 "sgesvd.f"
			itauq = ie + *n;
#line 1369 "sgesvd.f"
			itaup = itauq + *n;
#line 1370 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1375 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1375 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1379 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1385 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1385 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N+4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1393 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1393 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1395 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1402 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1410 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1413 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1417 "sgesvd.f"
			itau = 1;
#line 1418 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1423 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1423 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1425 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1430 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1430 "sgesvd.f"
			sorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1435 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1436 "sgesvd.f"
			if (*n > 1) {
#line 1436 "sgesvd.f"
			    i__2 = *n - 1;
#line 1436 "sgesvd.f"
			    i__3 = *n - 1;
#line 1436 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1436 "sgesvd.f"
			}
#line 1439 "sgesvd.f"
			ie = itau;
#line 1440 "sgesvd.f"
			itauq = ie + *n;
#line 1441 "sgesvd.f"
			itaup = itauq + *n;
#line 1442 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1447 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1447 "sgesvd.f"
			sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1455 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1455 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1462 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1462 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1464 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1471 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1475 "sgesvd.f"
		    }

#line 1477 "sgesvd.f"
		}

#line 1479 "sgesvd.f"
	    } else if (wntua) {

#line 1481 "sgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1487 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1487 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1491 "sgesvd.f"
			ir = 1;
#line 1492 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1496 "sgesvd.f"
			    ldwrkr = *lda;
#line 1497 "sgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1501 "sgesvd.f"
			    ldwrkr = *n;
#line 1502 "sgesvd.f"
			}
#line 1503 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1504 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1509 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1509 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1511 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1515 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1517 "sgesvd.f"
			i__2 = *n - 1;
#line 1517 "sgesvd.f"
			i__3 = *n - 1;
#line 1517 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 1523 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1523 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1525 "sgesvd.f"
			ie = itau;
#line 1526 "sgesvd.f"
			itauq = ie + *n;
#line 1527 "sgesvd.f"
			itaup = itauq + *n;
#line 1528 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1533 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1533 "sgesvd.f"
			sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1541 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1541 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1544 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1550 "sgesvd.f"
			sbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (Workspace: need N*N) */

#line 1558 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1563 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1565 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1569 "sgesvd.f"
			itau = 1;
#line 1570 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1575 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1575 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1577 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1582 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1582 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1584 "sgesvd.f"
			ie = itau;
#line 1585 "sgesvd.f"
			itauq = ie + *n;
#line 1586 "sgesvd.f"
			itaup = itauq + *n;
#line 1587 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1591 "sgesvd.f"
			i__2 = *n - 1;
#line 1591 "sgesvd.f"
			i__3 = *n - 1;
#line 1591 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1597 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1597 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1605 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1605 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1608 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1614 "sgesvd.f"
			sbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1618 "sgesvd.f"
		    }

#line 1620 "sgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1626 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1626 "sgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1630 "sgesvd.f"
			iu = 1;
#line 1631 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1635 "sgesvd.f"
			    ldwrku = *lda;
#line 1636 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1637 "sgesvd.f"
			    ldwrkr = *lda;
#line 1638 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1642 "sgesvd.f"
			    ldwrku = *lda;
#line 1643 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1644 "sgesvd.f"
			    ldwrkr = *n;
#line 1645 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1649 "sgesvd.f"
			    ldwrku = *n;
#line 1650 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1651 "sgesvd.f"
			    ldwrkr = *n;
#line 1652 "sgesvd.f"
			}
#line 1653 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1654 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1659 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1659 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1661 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */

#line 1666 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1666 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1671 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1673 "sgesvd.f"
			i__2 = *n - 1;
#line 1673 "sgesvd.f"
			i__3 = *n - 1;
#line 1673 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1675 "sgesvd.f"
			ie = itau;
#line 1676 "sgesvd.f"
			itauq = ie + *n;
#line 1677 "sgesvd.f"
			itaup = itauq + *n;
#line 1678 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1685 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1685 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1689 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB) */

#line 1695 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1695 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1703 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1703 "sgesvd.f"
			sorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1706 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N+BDSPAC) */

#line 1713 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1721 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1726 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1730 "sgesvd.f"
			slacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1733 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1737 "sgesvd.f"
			itau = 1;
#line 1738 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1743 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1743 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1745 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1750 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1750 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1752 "sgesvd.f"
			ie = itau;
#line 1753 "sgesvd.f"
			itauq = ie + *n;
#line 1754 "sgesvd.f"
			itaup = itauq + *n;
#line 1755 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1759 "sgesvd.f"
			i__2 = *n - 1;
#line 1759 "sgesvd.f"
			i__3 = *n - 1;
#line 1759 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 
				+ 2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1765 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1765 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1773 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1773 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1780 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1780 "sgesvd.f"
			sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1782 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1789 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1793 "sgesvd.f"
		    }

#line 1795 "sgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1802 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1802 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1806 "sgesvd.f"
			iu = 1;
#line 1807 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1811 "sgesvd.f"
			    ldwrku = *lda;
#line 1812 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1816 "sgesvd.f"
			    ldwrku = *n;
#line 1817 "sgesvd.f"
			}
#line 1818 "sgesvd.f"
			itau = iu + ldwrku * *n;
#line 1819 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1824 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1824 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1826 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 1831 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1831 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1836 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1838 "sgesvd.f"
			i__2 = *n - 1;
#line 1838 "sgesvd.f"
			i__3 = *n - 1;
#line 1838 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1840 "sgesvd.f"
			ie = itau;
#line 1841 "sgesvd.f"
			itauq = ie + *n;
#line 1842 "sgesvd.f"
			itaup = itauq + *n;
#line 1843 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1848 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1848 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1852 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1858 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1858 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N+4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1866 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1866 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1868 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1875 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1883 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1888 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1890 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1894 "sgesvd.f"
			itau = 1;
#line 1895 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1900 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1900 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1902 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1907 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1907 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 1912 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1913 "sgesvd.f"
			if (*n > 1) {
#line 1913 "sgesvd.f"
			    i__2 = *n - 1;
#line 1913 "sgesvd.f"
			    i__3 = *n - 1;
#line 1913 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1913 "sgesvd.f"
			}
#line 1916 "sgesvd.f"
			ie = itau;
#line 1917 "sgesvd.f"
			itauq = ie + *n;
#line 1918 "sgesvd.f"
			itaup = itauq + *n;
#line 1919 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1924 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1924 "sgesvd.f"
			sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1932 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1932 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1939 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1939 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1941 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1948 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1952 "sgesvd.f"
		    }

#line 1954 "sgesvd.f"
		}

#line 1956 "sgesvd.f"
	    }

#line 1958 "sgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 1965 "sgesvd.f"
	    ie = 1;
#line 1966 "sgesvd.f"
	    itauq = ie + *n;
#line 1967 "sgesvd.f"
	    itaup = itauq + *n;
#line 1968 "sgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 1973 "sgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 1973 "sgesvd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 1976 "sgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 3*N+NCU, prefer 3*N+NCU*NB) */

#line 1982 "sgesvd.f"
		slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1983 "sgesvd.f"
		if (wntus) {
#line 1983 "sgesvd.f"
		    ncu = *n;
#line 1983 "sgesvd.f"
		}
#line 1985 "sgesvd.f"
		if (wntua) {
#line 1985 "sgesvd.f"
		    ncu = *m;
#line 1985 "sgesvd.f"
		}
#line 1987 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 1987 "sgesvd.f"
		sorgbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1989 "sgesvd.f"
	    }
#line 1990 "sgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1996 "sgesvd.f"
		slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1997 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 1997 "sgesvd.f"
		sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1999 "sgesvd.f"
	    }
#line 2000 "sgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N, prefer 3*N+N*NB) */

#line 2006 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2006 "sgesvd.f"
		sorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2008 "sgesvd.f"
	    }
#line 2009 "sgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 2015 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2015 "sgesvd.f"
		sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2017 "sgesvd.f"
	    }
#line 2018 "sgesvd.f"
	    iwork = ie + *n;
#line 2019 "sgesvd.f"
	    if (wntuas || wntuo) {
#line 2019 "sgesvd.f"
		nru = *m;
#line 2019 "sgesvd.f"
	    }
#line 2021 "sgesvd.f"
	    if (wntun) {
#line 2021 "sgesvd.f"
		nru = 0;
#line 2021 "sgesvd.f"
	    }
#line 2023 "sgesvd.f"
	    if (wntvas || wntvo) {
#line 2023 "sgesvd.f"
		ncvt = *n;
#line 2023 "sgesvd.f"
	    }
#line 2025 "sgesvd.f"
	    if (wntvn) {
#line 2025 "sgesvd.f"
		ncvt = 0;
#line 2025 "sgesvd.f"
	    }
#line 2027 "sgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2034 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2036 "sgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 2043 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 2045 "sgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2052 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2054 "sgesvd.f"
	    }

#line 2056 "sgesvd.f"
	}

#line 2058 "sgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2064 "sgesvd.f"
	if (*n >= mnthr) {

#line 2066 "sgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2071 "sgesvd.f"
		itau = 1;
#line 2072 "sgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need 2*M, prefer M+M*NB) */

#line 2077 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2077 "sgesvd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2082 "sgesvd.f"
		i__2 = *m - 1;
#line 2082 "sgesvd.f"
		i__3 = *m - 1;
#line 2082 "sgesvd.f"
		slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 2083 "sgesvd.f"
		ie = 1;
#line 2084 "sgesvd.f"
		itauq = ie + *m;
#line 2085 "sgesvd.f"
		itaup = itauq + *m;
#line 2086 "sgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2091 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2091 "sgesvd.f"
		sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2094 "sgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2099 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2099 "sgesvd.f"
		    sorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2101 "sgesvd.f"
		}
#line 2102 "sgesvd.f"
		iwork = ie + *m;
#line 2103 "sgesvd.f"
		nru = 0;
#line 2104 "sgesvd.f"
		if (wntuo || wntuas) {
#line 2104 "sgesvd.f"
		    nru = *m;
#line 2104 "sgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 2111 "sgesvd.f"
		sbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &work[ie], dum, &
			c__1, &a[a_offset], lda, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2116 "sgesvd.f"
		if (wntuas) {
#line 2116 "sgesvd.f"
		    slacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2116 "sgesvd.f"
		}

#line 2119 "sgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

/* Computing MAX */
#line 2125 "sgesvd.f"
		i__2 = *m << 2;
#line 2125 "sgesvd.f"
		if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2129 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2130 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2130 "sgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2134 "sgesvd.f"
			ldwrku = *lda;
#line 2135 "sgesvd.f"
			chunk = *n;
#line 2136 "sgesvd.f"
			ldwrkr = *lda;
#line 2137 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2137 "sgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2137 "sgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2141 "sgesvd.f"
			    ldwrku = *lda;
#line 2142 "sgesvd.f"
			    chunk = *n;
#line 2143 "sgesvd.f"
			    ldwrkr = *m;
#line 2144 "sgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2148 "sgesvd.f"
			    ldwrku = *m;
#line 2149 "sgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2150 "sgesvd.f"
			    ldwrkr = *m;
#line 2151 "sgesvd.f"
			}
#line 2151 "sgesvd.f"
		    }
#line 2152 "sgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2153 "sgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2158 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2158 "sgesvd.f"
		    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2163 "sgesvd.f"
		    slacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2164 "sgesvd.f"
		    i__2 = *m - 1;
#line 2164 "sgesvd.f"
		    i__3 = *m - 1;
#line 2164 "sgesvd.f"
		    slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2170 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2170 "sgesvd.f"
		    sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2172 "sgesvd.f"
		    ie = itau;
#line 2173 "sgesvd.f"
		    itauq = ie + *m;
#line 2174 "sgesvd.f"
		    itaup = itauq + *m;
#line 2175 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2180 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2180 "sgesvd.f"
		    sgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */

#line 2187 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2187 "sgesvd.f"
		    sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2190 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M+BDSPAC) */

#line 2196 "sgesvd.f"
		    sbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[
			    ir], &ldwrkr, dum, &c__1, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 2199 "sgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M*N+M) */

#line 2205 "sgesvd.f"
		    i__2 = *n;
#line 2205 "sgesvd.f"
		    i__3 = chunk;
#line 2205 "sgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2206 "sgesvd.f"
			i__4 = *n - i__ + 1;
#line 2206 "sgesvd.f"
			blk = min(i__4,chunk);
#line 2207 "sgesvd.f"
			sgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2210 "sgesvd.f"
			slacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2212 "sgesvd.f"
/* L30: */
#line 2212 "sgesvd.f"
		    }

#line 2214 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2218 "sgesvd.f"
		    ie = 1;
#line 2219 "sgesvd.f"
		    itauq = ie + *m;
#line 2220 "sgesvd.f"
		    itaup = itauq + *m;
#line 2221 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 2226 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2226 "sgesvd.f"
		    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2233 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2233 "sgesvd.f"
		    sorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2235 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2241 "sgesvd.f"
		    sbdsqr_("L", m, n, &c__0, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, dum, &c__1, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);

#line 2244 "sgesvd.f"
		}

#line 2246 "sgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

/* Computing MAX */
#line 2252 "sgesvd.f"
		i__3 = *m << 2;
#line 2252 "sgesvd.f"
		if (*lwork >= *m * *m + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2256 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2257 "sgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2257 "sgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2261 "sgesvd.f"
			ldwrku = *lda;
#line 2262 "sgesvd.f"
			chunk = *n;
#line 2263 "sgesvd.f"
			ldwrkr = *lda;
#line 2264 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2264 "sgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2264 "sgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2268 "sgesvd.f"
			    ldwrku = *lda;
#line 2269 "sgesvd.f"
			    chunk = *n;
#line 2270 "sgesvd.f"
			    ldwrkr = *m;
#line 2271 "sgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2275 "sgesvd.f"
			    ldwrku = *m;
#line 2276 "sgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2277 "sgesvd.f"
			    ldwrkr = *m;
#line 2278 "sgesvd.f"
			}
#line 2278 "sgesvd.f"
		    }
#line 2279 "sgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2280 "sgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2285 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2285 "sgesvd.f"
		    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2290 "sgesvd.f"
		    slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2291 "sgesvd.f"
		    i__3 = *m - 1;
#line 2291 "sgesvd.f"
		    i__2 = *m - 1;
#line 2291 "sgesvd.f"
		    slaset_("U", &i__3, &i__2, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2297 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2297 "sgesvd.f"
		    sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2299 "sgesvd.f"
		    ie = itau;
#line 2300 "sgesvd.f"
		    itauq = ie + *m;
#line 2301 "sgesvd.f"
		    itaup = itauq + *m;
#line 2302 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2307 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2307 "sgesvd.f"
		    sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2310 "sgesvd.f"
		    slacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */

#line 2315 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2315 "sgesvd.f"
		    sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 2322 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2322 "sgesvd.f"
		    sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2324 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M+BDSPAC) */

#line 2331 "sgesvd.f"
		    sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[ir], 
			    &ldwrkr, &u[u_offset], ldu, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);
#line 2334 "sgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M*N+M)) */

#line 2340 "sgesvd.f"
		    i__3 = *n;
#line 2340 "sgesvd.f"
		    i__2 = chunk;
#line 2340 "sgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2341 "sgesvd.f"
			i__4 = *n - i__ + 1;
#line 2341 "sgesvd.f"
			blk = min(i__4,chunk);
#line 2342 "sgesvd.f"
			sgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2345 "sgesvd.f"
			slacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2347 "sgesvd.f"
/* L40: */
#line 2347 "sgesvd.f"
		    }

#line 2349 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2353 "sgesvd.f"
		    itau = 1;
#line 2354 "sgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need 2*M, prefer M+M*NB) */

#line 2359 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2359 "sgesvd.f"
		    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2364 "sgesvd.f"
		    slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2365 "sgesvd.f"
		    i__2 = *m - 1;
#line 2365 "sgesvd.f"
		    i__3 = *m - 1;
#line 2365 "sgesvd.f"
		    slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need 2*M, prefer M+M*NB) */

#line 2371 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2371 "sgesvd.f"
		    sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2373 "sgesvd.f"
		    ie = itau;
#line 2374 "sgesvd.f"
		    itauq = ie + *m;
#line 2375 "sgesvd.f"
		    itaup = itauq + *m;
#line 2376 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2381 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2381 "sgesvd.f"
		    sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2388 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2388 "sgesvd.f"
		    sormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2395 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2395 "sgesvd.f"
		    sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2397 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2404 "sgesvd.f"
		    sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 2407 "sgesvd.f"
		}

#line 2409 "sgesvd.f"
	    } else if (wntvs) {

#line 2411 "sgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2417 "sgesvd.f"
		    i__2 = *m << 2;
#line 2417 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2421 "sgesvd.f"
			ir = 1;
#line 2422 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2426 "sgesvd.f"
			    ldwrkr = *lda;
#line 2427 "sgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2431 "sgesvd.f"
			    ldwrkr = *m;
#line 2432 "sgesvd.f"
			}
#line 2433 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2434 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2439 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2439 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2444 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2446 "sgesvd.f"
			i__2 = *m - 1;
#line 2446 "sgesvd.f"
			i__3 = *m - 1;
#line 2446 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2452 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2452 "sgesvd.f"
			sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2454 "sgesvd.f"
			ie = itau;
#line 2455 "sgesvd.f"
			itauq = ie + *m;
#line 2456 "sgesvd.f"
			itaup = itauq + *m;
#line 2457 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2462 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2462 "sgesvd.f"
			sgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB) */

#line 2471 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2471 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2474 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2480 "sgesvd.f"
			sbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2488 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2491 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2495 "sgesvd.f"
			itau = 1;
#line 2496 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2501 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2501 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2506 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2511 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2511 "sgesvd.f"
			sorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2513 "sgesvd.f"
			ie = itau;
#line 2514 "sgesvd.f"
			itauq = ie + *m;
#line 2515 "sgesvd.f"
			itaup = itauq + *m;
#line 2516 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2520 "sgesvd.f"
			i__2 = *m - 1;
#line 2520 "sgesvd.f"
			i__3 = *m - 1;
#line 2520 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2526 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2526 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2533 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2533 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2536 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2542 "sgesvd.f"
			sbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 2546 "sgesvd.f"
		    }

#line 2548 "sgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 2554 "sgesvd.f"
		    i__2 = *m << 2;
#line 2554 "sgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2558 "sgesvd.f"
			iu = 1;
#line 2559 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2563 "sgesvd.f"
			    ldwrku = *lda;
#line 2564 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2565 "sgesvd.f"
			    ldwrkr = *lda;
#line 2566 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2570 "sgesvd.f"
			    ldwrku = *lda;
#line 2571 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2572 "sgesvd.f"
			    ldwrkr = *m;
#line 2573 "sgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2577 "sgesvd.f"
			    ldwrku = *m;
#line 2578 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2579 "sgesvd.f"
			    ldwrkr = *m;
#line 2580 "sgesvd.f"
			}
#line 2581 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2582 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 2587 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2587 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2592 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2594 "sgesvd.f"
			i__2 = *m - 1;
#line 2594 "sgesvd.f"
			i__3 = *m - 1;
#line 2594 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 2600 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2600 "sgesvd.f"
			sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2602 "sgesvd.f"
			ie = itau;
#line 2603 "sgesvd.f"
			itauq = ie + *m;
#line 2604 "sgesvd.f"
			itaup = itauq + *m;
#line 2605 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 2612 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2612 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2616 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M+4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 2623 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2623 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB) */

#line 2630 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2630 "sgesvd.f"
			sorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2633 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M+BDSPAC) */

#line 2640 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2648 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (Workspace: need M*M) */

#line 2654 "sgesvd.f"
			slacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2657 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2661 "sgesvd.f"
			itau = 1;
#line 2662 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2667 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2667 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2669 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2674 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2674 "sgesvd.f"
			sorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2676 "sgesvd.f"
			ie = itau;
#line 2677 "sgesvd.f"
			itauq = ie + *m;
#line 2678 "sgesvd.f"
			itaup = itauq + *m;
#line 2679 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2683 "sgesvd.f"
			i__2 = *m - 1;
#line 2683 "sgesvd.f"
			i__3 = *m - 1;
#line 2683 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2689 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2689 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2696 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2696 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2703 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2703 "sgesvd.f"
			sorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2705 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, compute left */
/*                    singular vectors of A in A and compute right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2712 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2716 "sgesvd.f"
		    }

#line 2718 "sgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 2725 "sgesvd.f"
		    i__2 = *m << 2;
#line 2725 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2729 "sgesvd.f"
			iu = 1;
#line 2730 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2734 "sgesvd.f"
			    ldwrku = *lda;
#line 2735 "sgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2739 "sgesvd.f"
			    ldwrku = *m;
#line 2740 "sgesvd.f"
			}
#line 2741 "sgesvd.f"
			itau = iu + ldwrku * *m;
#line 2742 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2747 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2747 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2752 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2754 "sgesvd.f"
			i__2 = *m - 1;
#line 2754 "sgesvd.f"
			i__3 = *m - 1;
#line 2754 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2760 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2760 "sgesvd.f"
			sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2762 "sgesvd.f"
			ie = itau;
#line 2763 "sgesvd.f"
			itauq = ie + *m;
#line 2764 "sgesvd.f"
			itaup = itauq + *m;
#line 2765 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2770 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2770 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2774 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M+4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2781 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2781 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 2788 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2788 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2790 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2797 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2805 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2808 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2812 "sgesvd.f"
			itau = 1;
#line 2813 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2818 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2818 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2820 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2825 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2825 "sgesvd.f"
			sorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2830 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2831 "sgesvd.f"
			i__2 = *m - 1;
#line 2831 "sgesvd.f"
			i__3 = *m - 1;
#line 2831 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 2833 "sgesvd.f"
			ie = itau;
#line 2834 "sgesvd.f"
			itauq = ie + *m;
#line 2835 "sgesvd.f"
			itaup = itauq + *m;
#line 2836 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2841 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2841 "sgesvd.f"
			sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2849 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2849 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2856 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2856 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2858 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2865 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2869 "sgesvd.f"
		    }

#line 2871 "sgesvd.f"
		}

#line 2873 "sgesvd.f"
	    } else if (wntva) {

#line 2875 "sgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2881 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 2881 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2885 "sgesvd.f"
			ir = 1;
#line 2886 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2890 "sgesvd.f"
			    ldwrkr = *lda;
#line 2891 "sgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2895 "sgesvd.f"
			    ldwrkr = *m;
#line 2896 "sgesvd.f"
			}
#line 2897 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2898 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2903 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2903 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2905 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2909 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2911 "sgesvd.f"
			i__2 = *m - 1;
#line 2911 "sgesvd.f"
			i__3 = *m - 1;
#line 2911 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB) */

#line 2917 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2917 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2919 "sgesvd.f"
			ie = itau;
#line 2920 "sgesvd.f"
			itauq = ie + *m;
#line 2921 "sgesvd.f"
			itaup = itauq + *m;
#line 2922 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2927 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2927 "sgesvd.f"
			sgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need M*M+4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2936 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2936 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2939 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2945 "sgesvd.f"
			sbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 2953 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 2958 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 2960 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2964 "sgesvd.f"
			itau = 1;
#line 2965 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2970 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2970 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2972 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 2977 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2977 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2979 "sgesvd.f"
			ie = itau;
#line 2980 "sgesvd.f"
			itauq = ie + *m;
#line 2981 "sgesvd.f"
			itaup = itauq + *m;
#line 2982 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2986 "sgesvd.f"
			i__2 = *m - 1;
#line 2986 "sgesvd.f"
			i__3 = *m - 1;
#line 2986 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2992 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2992 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3000 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3000 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3003 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3009 "sgesvd.f"
			sbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 3013 "sgesvd.f"
		    }

#line 3015 "sgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3021 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3021 "sgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3025 "sgesvd.f"
			iu = 1;
#line 3026 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3030 "sgesvd.f"
			    ldwrku = *lda;
#line 3031 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3032 "sgesvd.f"
			    ldwrkr = *lda;
#line 3033 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3037 "sgesvd.f"
			    ldwrku = *lda;
#line 3038 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3039 "sgesvd.f"
			    ldwrkr = *m;
#line 3040 "sgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3044 "sgesvd.f"
			    ldwrku = *m;
#line 3045 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3046 "sgesvd.f"
			    ldwrkr = *m;
#line 3047 "sgesvd.f"
			}
#line 3048 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3049 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 3054 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3054 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3056 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */

#line 3061 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3061 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3066 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3068 "sgesvd.f"
			i__2 = *m - 1;
#line 3068 "sgesvd.f"
			i__3 = *m - 1;
#line 3068 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3070 "sgesvd.f"
			ie = itau;
#line 3071 "sgesvd.f"
			itauq = ie + *m;
#line 3072 "sgesvd.f"
			itaup = itauq + *m;
#line 3073 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 3080 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3080 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3084 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M+4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 3091 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3091 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB) */

#line 3098 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3098 "sgesvd.f"
			sorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3101 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M+BDSPAC) */

#line 3108 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3116 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3121 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3125 "sgesvd.f"
			slacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3128 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3132 "sgesvd.f"
			itau = 1;
#line 3133 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 3138 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3138 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3140 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 3145 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3145 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3147 "sgesvd.f"
			ie = itau;
#line 3148 "sgesvd.f"
			itauq = ie + *m;
#line 3149 "sgesvd.f"
			itaup = itauq + *m;
#line 3150 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3154 "sgesvd.f"
			i__2 = *m - 1;
#line 3154 "sgesvd.f"
			i__3 = *m - 1;
#line 3154 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 3160 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3160 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3168 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3168 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3175 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3175 "sgesvd.f"
			sorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3177 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3184 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3188 "sgesvd.f"
		    }

#line 3190 "sgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3197 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3197 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3201 "sgesvd.f"
			iu = 1;
#line 3202 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3206 "sgesvd.f"
			    ldwrku = *lda;
#line 3207 "sgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3211 "sgesvd.f"
			    ldwrku = *m;
#line 3212 "sgesvd.f"
			}
#line 3213 "sgesvd.f"
			itau = iu + ldwrku * *m;
#line 3214 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 3219 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3219 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3221 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB) */

#line 3226 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3226 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3231 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3233 "sgesvd.f"
			i__2 = *m - 1;
#line 3233 "sgesvd.f"
			i__3 = *m - 1;
#line 3233 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3235 "sgesvd.f"
			ie = itau;
#line 3236 "sgesvd.f"
			itauq = ie + *m;
#line 3237 "sgesvd.f"
			itaup = itauq + *m;
#line 3238 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 3243 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3243 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3247 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB) */

#line 3253 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3253 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 3260 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3260 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3262 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 3269 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3277 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3282 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3284 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3288 "sgesvd.f"
			itau = 1;
#line 3289 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 3294 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3294 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3296 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 3301 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3301 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3306 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3307 "sgesvd.f"
			i__2 = *m - 1;
#line 3307 "sgesvd.f"
			i__3 = *m - 1;
#line 3307 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 3309 "sgesvd.f"
			ie = itau;
#line 3310 "sgesvd.f"
			itauq = ie + *m;
#line 3311 "sgesvd.f"
			itaup = itauq + *m;
#line 3312 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 3317 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3317 "sgesvd.f"
			sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3325 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3325 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3332 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3332 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3334 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3341 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3345 "sgesvd.f"
		    }

#line 3347 "sgesvd.f"
		}

#line 3349 "sgesvd.f"
	    }

#line 3351 "sgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3358 "sgesvd.f"
	    ie = 1;
#line 3359 "sgesvd.f"
	    itauq = ie + *m;
#line 3360 "sgesvd.f"
	    itaup = itauq + *m;
#line 3361 "sgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 3366 "sgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3366 "sgesvd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 3369 "sgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) */

#line 3375 "sgesvd.f"
		slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3376 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3376 "sgesvd.f"
		sorgbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3378 "sgesvd.f"
	    }
#line 3379 "sgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB) */

#line 3385 "sgesvd.f"
		slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3386 "sgesvd.f"
		if (wntva) {
#line 3386 "sgesvd.f"
		    nrvt = *n;
#line 3386 "sgesvd.f"
		}
#line 3388 "sgesvd.f"
		if (wntvs) {
#line 3388 "sgesvd.f"
		    nrvt = *m;
#line 3388 "sgesvd.f"
		}
#line 3390 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3390 "sgesvd.f"
		sorgbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3392 "sgesvd.f"
	    }
#line 3393 "sgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) */

#line 3399 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3399 "sgesvd.f"
		sorgbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3401 "sgesvd.f"
	    }
#line 3402 "sgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3408 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3408 "sgesvd.f"
		sorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3410 "sgesvd.f"
	    }
#line 3411 "sgesvd.f"
	    iwork = ie + *m;
#line 3412 "sgesvd.f"
	    if (wntuas || wntuo) {
#line 3412 "sgesvd.f"
		nru = *m;
#line 3412 "sgesvd.f"
	    }
#line 3414 "sgesvd.f"
	    if (wntun) {
#line 3414 "sgesvd.f"
		nru = 0;
#line 3414 "sgesvd.f"
	    }
#line 3416 "sgesvd.f"
	    if (wntvas || wntvo) {
#line 3416 "sgesvd.f"
		ncvt = *n;
#line 3416 "sgesvd.f"
	    }
#line 3418 "sgesvd.f"
	    if (wntvn) {
#line 3418 "sgesvd.f"
		ncvt = 0;
#line 3418 "sgesvd.f"
	    }
#line 3420 "sgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3427 "sgesvd.f"
		sbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3429 "sgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 3436 "sgesvd.f"
		sbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 3438 "sgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3445 "sgesvd.f"
		sbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3447 "sgesvd.f"
	    }

#line 3449 "sgesvd.f"
	}

#line 3451 "sgesvd.f"
    }

/*     If SBDSQR failed to converge, copy unconverged superdiagonals */
/*     to WORK( 2:MINMN ) */

#line 3456 "sgesvd.f"
    if (*info != 0) {
#line 3457 "sgesvd.f"
	if (ie > 2) {
#line 3458 "sgesvd.f"
	    i__2 = minmn - 1;
#line 3458 "sgesvd.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 3459 "sgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3460 "sgesvd.f"
/* L50: */
#line 3460 "sgesvd.f"
	    }
#line 3461 "sgesvd.f"
	}
#line 3462 "sgesvd.f"
	if (ie < 2) {
#line 3463 "sgesvd.f"
	    for (i__ = minmn - 1; i__ >= 1; --i__) {
#line 3464 "sgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3465 "sgesvd.f"
/* L60: */
#line 3465 "sgesvd.f"
	    }
#line 3466 "sgesvd.f"
	}
#line 3467 "sgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3471 "sgesvd.f"
    if (iscl == 1) {
#line 3472 "sgesvd.f"
	if (anrm > bignum) {
#line 3472 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3472 "sgesvd.f"
	}
#line 3475 "sgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3475 "sgesvd.f"
	    i__2 = minmn - 1;
#line 3475 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3475 "sgesvd.f"
	}
#line 3478 "sgesvd.f"
	if (anrm < smlnum) {
#line 3478 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3478 "sgesvd.f"
	}
#line 3481 "sgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3481 "sgesvd.f"
	    i__2 = minmn - 1;
#line 3481 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3481 "sgesvd.f"
	}
#line 3484 "sgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3488 "sgesvd.f"
    work[1] = (doublereal) maxwrk;

#line 3490 "sgesvd.f"
    return 0;

/*     End of SGESVD */

} /* sgesvd_ */

