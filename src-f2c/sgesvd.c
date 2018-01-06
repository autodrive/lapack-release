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
		if (*n > 1) {
#line 697 "sgesvd.f"
		    i__2 = *n - 1;
#line 697 "sgesvd.f"
		    i__3 = *n - 1;
#line 697 "sgesvd.f"
		    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2],
			     lda, (ftnlen)1);
#line 699 "sgesvd.f"
		}
#line 700 "sgesvd.f"
		ie = 1;
#line 701 "sgesvd.f"
		itauq = ie + *n;
#line 702 "sgesvd.f"
		itaup = itauq + *n;
#line 703 "sgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 708 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 708 "sgesvd.f"
		sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 711 "sgesvd.f"
		ncvt = 0;
#line 712 "sgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 717 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 717 "sgesvd.f"
		    sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 719 "sgesvd.f"
		    ncvt = *n;
#line 720 "sgesvd.f"
		}
#line 721 "sgesvd.f"
		iwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 727 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, dum, &c__1, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 732 "sgesvd.f"
		if (wntvas) {
#line 732 "sgesvd.f"
		    slacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 732 "sgesvd.f"
		}

#line 735 "sgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

/* Computing MAX */
#line 741 "sgesvd.f"
		i__2 = *n << 2;
#line 741 "sgesvd.f"
		if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 745 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 746 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 746 "sgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 750 "sgesvd.f"
			ldwrku = *lda;
#line 751 "sgesvd.f"
			ldwrkr = *lda;
#line 752 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 752 "sgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *n;
#line 752 "sgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 756 "sgesvd.f"
			    ldwrku = *lda;
#line 757 "sgesvd.f"
			    ldwrkr = *n;
#line 758 "sgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 762 "sgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 763 "sgesvd.f"
			    ldwrkr = *n;
#line 764 "sgesvd.f"
			}
#line 764 "sgesvd.f"
		    }
#line 765 "sgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 766 "sgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 771 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 771 "sgesvd.f"
		    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 776 "sgesvd.f"
		    slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 777 "sgesvd.f"
		    i__2 = *n - 1;
#line 777 "sgesvd.f"
		    i__3 = *n - 1;
#line 777 "sgesvd.f"
		    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1], 
			    &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 783 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 783 "sgesvd.f"
		    sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 785 "sgesvd.f"
		    ie = itau;
#line 786 "sgesvd.f"
		    itauq = ie + *n;
#line 787 "sgesvd.f"
		    itaup = itauq + *n;
#line 788 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 793 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 793 "sgesvd.f"
		    sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 800 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 800 "sgesvd.f"
		    sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 803 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (Workspace: need N*N+BDSPAC) */

#line 809 "sgesvd.f"
		    sbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &
			    c__1, &work[ir], &ldwrkr, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 812 "sgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N+2*N, prefer N*N+M*N+N) */

#line 818 "sgesvd.f"
		    i__2 = *m;
#line 818 "sgesvd.f"
		    i__3 = ldwrku;
#line 818 "sgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 819 "sgesvd.f"
			i__4 = *m - i__ + 1;
#line 819 "sgesvd.f"
			chunk = min(i__4,ldwrku);
#line 820 "sgesvd.f"
			sgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 823 "sgesvd.f"
			slacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 825 "sgesvd.f"
/* L10: */
#line 825 "sgesvd.f"
		    }

#line 827 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 831 "sgesvd.f"
		    ie = 1;
#line 832 "sgesvd.f"
		    itauq = ie + *n;
#line 833 "sgesvd.f"
		    itaup = itauq + *n;
#line 834 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 839 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 839 "sgesvd.f"
		    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (Workspace: need 4*N, prefer 3*N+N*NB) */

#line 846 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 846 "sgesvd.f"
		    sorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 848 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 854 "sgesvd.f"
		    sbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &
			    c__1, &a[a_offset], lda, dum, &c__1, &work[iwork],
			     info, (ftnlen)1);

#line 857 "sgesvd.f"
		}

#line 859 "sgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

/* Computing MAX */
#line 865 "sgesvd.f"
		i__3 = *n << 2;
#line 865 "sgesvd.f"
		if (*lwork >= *n * *n + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 869 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 870 "sgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 870 "sgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 874 "sgesvd.f"
			ldwrku = *lda;
#line 875 "sgesvd.f"
			ldwrkr = *lda;
#line 876 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 876 "sgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *n;
#line 876 "sgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 880 "sgesvd.f"
			    ldwrku = *lda;
#line 881 "sgesvd.f"
			    ldwrkr = *n;
#line 882 "sgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 886 "sgesvd.f"
			    ldwrku = (*lwork - *n * *n - *n) / *n;
#line 887 "sgesvd.f"
			    ldwrkr = *n;
#line 888 "sgesvd.f"
			}
#line 888 "sgesvd.f"
		    }
#line 889 "sgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 890 "sgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 895 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 895 "sgesvd.f"
		    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 900 "sgesvd.f"
		    slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 901 "sgesvd.f"
		    if (*n > 1) {
#line 901 "sgesvd.f"
			i__3 = *n - 1;
#line 901 "sgesvd.f"
			i__2 = *n - 1;
#line 901 "sgesvd.f"
			slaset_("L", &i__3, &i__2, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 901 "sgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 908 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 908 "sgesvd.f"
		    sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 910 "sgesvd.f"
		    ie = itau;
#line 911 "sgesvd.f"
		    itauq = ie + *n;
#line 912 "sgesvd.f"
		    itaup = itauq + *n;
#line 913 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 918 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 918 "sgesvd.f"
		    sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 921 "sgesvd.f"
		    slacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 926 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 926 "sgesvd.f"
		    sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need N*N+4*N-1, prefer N*N+3*N+(N-1)*NB) */

#line 933 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 933 "sgesvd.f"
		    sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 935 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (Workspace: need N*N+BDSPAC) */

#line 942 "sgesvd.f"
		    sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, dum, &c__1, 
			    &work[iwork], info, (ftnlen)1);
#line 945 "sgesvd.f"
		    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N+2*N, prefer N*N+M*N+N) */

#line 951 "sgesvd.f"
		    i__3 = *m;
#line 951 "sgesvd.f"
		    i__2 = ldwrku;
#line 951 "sgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 952 "sgesvd.f"
			i__4 = *m - i__ + 1;
#line 952 "sgesvd.f"
			chunk = min(i__4,ldwrku);
#line 953 "sgesvd.f"
			sgemm_("N", "N", &chunk, n, n, &c_b79, &a[i__ + 
				a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 956 "sgesvd.f"
			slacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 958 "sgesvd.f"
/* L20: */
#line 958 "sgesvd.f"
		    }

#line 960 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 964 "sgesvd.f"
		    itau = 1;
#line 965 "sgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need 2*N, prefer N+N*NB) */

#line 970 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 970 "sgesvd.f"
		    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 975 "sgesvd.f"
		    slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 976 "sgesvd.f"
		    if (*n > 1) {
#line 976 "sgesvd.f"
			i__2 = *n - 1;
#line 976 "sgesvd.f"
			i__3 = *n - 1;
#line 976 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				vt_dim1 + 2], ldvt, (ftnlen)1);
#line 976 "sgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (Workspace: need 2*N, prefer N+N*NB) */

#line 983 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 983 "sgesvd.f"
		    sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 985 "sgesvd.f"
		    ie = itau;
#line 986 "sgesvd.f"
		    itauq = ie + *n;
#line 987 "sgesvd.f"
		    itaup = itauq + *n;
#line 988 "sgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 993 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 993 "sgesvd.f"
		    sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1000 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1000 "sgesvd.f"
		    sormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1007 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1007 "sgesvd.f"
		    sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1009 "sgesvd.f"
		    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (Workspace: need BDSPAC) */

#line 1016 "sgesvd.f"
		    sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 1019 "sgesvd.f"
		}

#line 1021 "sgesvd.f"
	    } else if (wntus) {

#line 1023 "sgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1029 "sgesvd.f"
		    i__2 = *n << 2;
#line 1029 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1033 "sgesvd.f"
			ir = 1;
#line 1034 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1038 "sgesvd.f"
			    ldwrkr = *lda;
#line 1039 "sgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1043 "sgesvd.f"
			    ldwrkr = *n;
#line 1044 "sgesvd.f"
			}
#line 1045 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1046 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1051 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1051 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1056 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1058 "sgesvd.f"
			i__2 = *n - 1;
#line 1058 "sgesvd.f"
			i__3 = *n - 1;
#line 1058 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1064 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1064 "sgesvd.f"
			sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1066 "sgesvd.f"
			ie = itau;
#line 1067 "sgesvd.f"
			itauq = ie + *n;
#line 1068 "sgesvd.f"
			itaup = itauq + *n;
#line 1069 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1074 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1074 "sgesvd.f"
			sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1082 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1082 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1085 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1091 "sgesvd.f"
			sbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (Workspace: need N*N) */

#line 1099 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1102 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1106 "sgesvd.f"
			itau = 1;
#line 1107 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1112 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1112 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1114 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1119 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1119 "sgesvd.f"
			sorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1121 "sgesvd.f"
			ie = itau;
#line 1122 "sgesvd.f"
			itauq = ie + *n;
#line 1123 "sgesvd.f"
			itaup = itauq + *n;
#line 1124 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1128 "sgesvd.f"
			if (*n > 1) {
#line 1129 "sgesvd.f"
			    i__2 = *n - 1;
#line 1129 "sgesvd.f"
			    i__3 = *n - 1;
#line 1129 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1131 "sgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1136 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1136 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1143 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1143 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1146 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1152 "sgesvd.f"
			sbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1156 "sgesvd.f"
		    }

#line 1158 "sgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1164 "sgesvd.f"
		    i__2 = *n << 2;
#line 1164 "sgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1168 "sgesvd.f"
			iu = 1;
#line 1169 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1173 "sgesvd.f"
			    ldwrku = *lda;
#line 1174 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1175 "sgesvd.f"
			    ldwrkr = *lda;
#line 1176 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1180 "sgesvd.f"
			    ldwrku = *lda;
#line 1181 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1182 "sgesvd.f"
			    ldwrkr = *n;
#line 1183 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1187 "sgesvd.f"
			    ldwrku = *n;
#line 1188 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1189 "sgesvd.f"
			    ldwrkr = *n;
#line 1190 "sgesvd.f"
			}
#line 1191 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1192 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1197 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1197 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1202 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1204 "sgesvd.f"
			i__2 = *n - 1;
#line 1204 "sgesvd.f"
			i__3 = *n - 1;
#line 1204 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1210 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1210 "sgesvd.f"
			sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1212 "sgesvd.f"
			ie = itau;
#line 1213 "sgesvd.f"
			itauq = ie + *n;
#line 1214 "sgesvd.f"
			itaup = itauq + *n;
#line 1215 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1222 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1222 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1226 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB) */

#line 1232 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1232 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1240 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1240 "sgesvd.f"
			sorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1243 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N+BDSPAC) */

#line 1250 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1258 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (Workspace: need N*N) */

#line 1264 "sgesvd.f"
			slacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1267 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1271 "sgesvd.f"
			itau = 1;
#line 1272 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1277 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1277 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1279 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1284 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1284 "sgesvd.f"
			sorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1286 "sgesvd.f"
			ie = itau;
#line 1287 "sgesvd.f"
			itauq = ie + *n;
#line 1288 "sgesvd.f"
			itaup = itauq + *n;
#line 1289 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1293 "sgesvd.f"
			if (*n > 1) {
#line 1294 "sgesvd.f"
			    i__2 = *n - 1;
#line 1294 "sgesvd.f"
			    i__3 = *n - 1;
#line 1294 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1296 "sgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1301 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1301 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1308 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1308 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1315 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1315 "sgesvd.f"
			sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1317 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1324 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1328 "sgesvd.f"
		    }

#line 1330 "sgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1337 "sgesvd.f"
		    i__2 = *n << 2;
#line 1337 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1341 "sgesvd.f"
			iu = 1;
#line 1342 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1346 "sgesvd.f"
			    ldwrku = *lda;
#line 1347 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1351 "sgesvd.f"
			    ldwrku = *n;
#line 1352 "sgesvd.f"
			}
#line 1353 "sgesvd.f"
			itau = iu + ldwrku * *n;
#line 1354 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1359 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1359 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1364 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1366 "sgesvd.f"
			i__2 = *n - 1;
#line 1366 "sgesvd.f"
			i__3 = *n - 1;
#line 1366 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1372 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1372 "sgesvd.f"
			sorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1374 "sgesvd.f"
			ie = itau;
#line 1375 "sgesvd.f"
			itauq = ie + *n;
#line 1376 "sgesvd.f"
			itaup = itauq + *n;
#line 1377 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1382 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1382 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1386 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1392 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1392 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N+4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1400 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1400 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1402 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1409 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

#line 1417 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b57, &u[u_offset], ldu, 
				(ftnlen)1, (ftnlen)1);

#line 1420 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1424 "sgesvd.f"
			itau = 1;
#line 1425 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1430 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1430 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1432 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1437 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1437 "sgesvd.f"
			sorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1442 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1443 "sgesvd.f"
			if (*n > 1) {
#line 1443 "sgesvd.f"
			    i__2 = *n - 1;
#line 1443 "sgesvd.f"
			    i__3 = *n - 1;
#line 1443 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1443 "sgesvd.f"
			}
#line 1446 "sgesvd.f"
			ie = itau;
#line 1447 "sgesvd.f"
			itauq = ie + *n;
#line 1448 "sgesvd.f"
			itaup = itauq + *n;
#line 1449 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1454 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1454 "sgesvd.f"
			sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1462 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1462 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1469 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1469 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1471 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1478 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1482 "sgesvd.f"
		    }

#line 1484 "sgesvd.f"
		}

#line 1486 "sgesvd.f"
	    } else if (wntua) {

#line 1488 "sgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1494 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1494 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1498 "sgesvd.f"
			ir = 1;
#line 1499 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1503 "sgesvd.f"
			    ldwrkr = *lda;
#line 1504 "sgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1508 "sgesvd.f"
			    ldwrkr = *n;
#line 1509 "sgesvd.f"
			}
#line 1510 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1511 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1516 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1516 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1518 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1522 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1524 "sgesvd.f"
			i__2 = *n - 1;
#line 1524 "sgesvd.f"
			i__3 = *n - 1;
#line 1524 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 1530 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1530 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1532 "sgesvd.f"
			ie = itau;
#line 1533 "sgesvd.f"
			itauq = ie + *n;
#line 1534 "sgesvd.f"
			itaup = itauq + *n;
#line 1535 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1540 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1540 "sgesvd.f"
			sgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1548 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1548 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1551 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1557 "sgesvd.f"
			sbdsqr_("U", n, &c__0, n, &c__0, &s[1], &work[ie], 
				dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (Workspace: need N*N) */

#line 1565 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1570 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1572 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1576 "sgesvd.f"
			itau = 1;
#line 1577 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1582 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1582 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1584 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1589 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1589 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1591 "sgesvd.f"
			ie = itau;
#line 1592 "sgesvd.f"
			itauq = ie + *n;
#line 1593 "sgesvd.f"
			itaup = itauq + *n;
#line 1594 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1598 "sgesvd.f"
			if (*n > 1) {
#line 1599 "sgesvd.f"
			    i__2 = *n - 1;
#line 1599 "sgesvd.f"
			    i__3 = *n - 1;
#line 1599 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1601 "sgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1606 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1606 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1614 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1614 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1617 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

#line 1623 "sgesvd.f"
			sbdsqr_("U", n, &c__0, m, &c__0, &s[1], &work[ie], 
				dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 1627 "sgesvd.f"
		    }

#line 1629 "sgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1635 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1635 "sgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1639 "sgesvd.f"
			iu = 1;
#line 1640 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1644 "sgesvd.f"
			    ldwrku = *lda;
#line 1645 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1646 "sgesvd.f"
			    ldwrkr = *lda;
#line 1647 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1651 "sgesvd.f"
			    ldwrku = *lda;
#line 1652 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1653 "sgesvd.f"
			    ldwrkr = *n;
#line 1654 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1658 "sgesvd.f"
			    ldwrku = *n;
#line 1659 "sgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1660 "sgesvd.f"
			    ldwrkr = *n;
#line 1661 "sgesvd.f"
			}
#line 1662 "sgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1663 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */

#line 1668 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1668 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1670 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */

#line 1675 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1675 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1680 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1682 "sgesvd.f"
			i__2 = *n - 1;
#line 1682 "sgesvd.f"
			i__3 = *n - 1;
#line 1682 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1684 "sgesvd.f"
			ie = itau;
#line 1685 "sgesvd.f"
			itauq = ie + *n;
#line 1686 "sgesvd.f"
			itaup = itauq + *n;
#line 1687 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

#line 1694 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1694 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1698 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB) */

#line 1704 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1704 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N+4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

#line 1712 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1712 "sgesvd.f"
			sorgbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1715 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N+BDSPAC) */

#line 1722 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1730 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1735 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1739 "sgesvd.f"
			slacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1742 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1746 "sgesvd.f"
			itau = 1;
#line 1747 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1752 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1752 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1754 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1759 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1759 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1761 "sgesvd.f"
			ie = itau;
#line 1762 "sgesvd.f"
			itauq = ie + *n;
#line 1763 "sgesvd.f"
			itaup = itauq + *n;
#line 1764 "sgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1768 "sgesvd.f"
			if (*n > 1) {
#line 1769 "sgesvd.f"
			    i__2 = *n - 1;
#line 1769 "sgesvd.f"
			    i__3 = *n - 1;
#line 1769 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1771 "sgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1776 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1776 "sgesvd.f"
			sgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1784 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1784 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1791 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1791 "sgesvd.f"
			sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1793 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

#line 1800 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &a[
				a_offset], lda, &u[u_offset], ldu, dum, &c__1,
				 &work[iwork], info, (ftnlen)1);

#line 1804 "sgesvd.f"
		    }

#line 1806 "sgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1813 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
#line 1813 "sgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1817 "sgesvd.f"
			iu = 1;
#line 1818 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1822 "sgesvd.f"
			    ldwrku = *lda;
#line 1823 "sgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1827 "sgesvd.f"
			    ldwrku = *n;
#line 1828 "sgesvd.f"
			}
#line 1829 "sgesvd.f"
			itau = iu + ldwrku * *n;
#line 1830 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB) */

#line 1835 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1835 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1837 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB) */

#line 1842 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1842 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1847 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1849 "sgesvd.f"
			i__2 = *n - 1;
#line 1849 "sgesvd.f"
			i__3 = *n - 1;
#line 1849 "sgesvd.f"
			slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				1], &ldwrku, (ftnlen)1);
#line 1851 "sgesvd.f"
			ie = itau;
#line 1852 "sgesvd.f"
			itauq = ie + *n;
#line 1853 "sgesvd.f"
			itaup = itauq + *n;
#line 1854 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB) */

#line 1859 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1859 "sgesvd.f"
			sgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1863 "sgesvd.f"
			slacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB) */

#line 1869 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1869 "sgesvd.f"
			sorgbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N+4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

#line 1877 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1877 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1879 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N+BDSPAC) */

#line 1886 "sgesvd.f"
			sbdsqr_("U", n, n, n, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

#line 1894 "sgesvd.f"
			sgemm_("N", "N", m, n, n, &c_b79, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b57, &a[a_offset], lda, 
				(ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1899 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1901 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1905 "sgesvd.f"
			itau = 1;
#line 1906 "sgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N+N*NB) */

#line 1911 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1911 "sgesvd.f"
			sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1913 "sgesvd.f"
			slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N+M, prefer N+M*NB) */

#line 1918 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1918 "sgesvd.f"
			sorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 1923 "sgesvd.f"
			slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1924 "sgesvd.f"
			if (*n > 1) {
#line 1924 "sgesvd.f"
			    i__2 = *n - 1;
#line 1924 "sgesvd.f"
			    i__3 = *n - 1;
#line 1924 "sgesvd.f"
			    slaset_("L", &i__2, &i__3, &c_b57, &c_b57, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1924 "sgesvd.f"
			}
#line 1927 "sgesvd.f"
			ie = itau;
#line 1928 "sgesvd.f"
			itauq = ie + *n;
#line 1929 "sgesvd.f"
			itaup = itauq + *n;
#line 1930 "sgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N+2*N*NB) */

#line 1935 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1935 "sgesvd.f"
			sgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], 
				&work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N+M, prefer 3*N+M*NB) */

#line 1943 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1943 "sgesvd.f"
			sormbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 1950 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1950 "sgesvd.f"
			sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1952 "sgesvd.f"
			iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 1959 "sgesvd.f"
			sbdsqr_("U", n, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 1963 "sgesvd.f"
		    }

#line 1965 "sgesvd.f"
		}

#line 1967 "sgesvd.f"
	    }

#line 1969 "sgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 1976 "sgesvd.f"
	    ie = 1;
#line 1977 "sgesvd.f"
	    itauq = ie + *n;
#line 1978 "sgesvd.f"
	    itaup = itauq + *n;
#line 1979 "sgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB) */

#line 1984 "sgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 1984 "sgesvd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 1987 "sgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 3*N+NCU, prefer 3*N+NCU*NB) */

#line 1993 "sgesvd.f"
		slacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1994 "sgesvd.f"
		if (wntus) {
#line 1994 "sgesvd.f"
		    ncu = *n;
#line 1994 "sgesvd.f"
		}
#line 1996 "sgesvd.f"
		if (wntua) {
#line 1996 "sgesvd.f"
		    ncu = *m;
#line 1996 "sgesvd.f"
		}
#line 1998 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 1998 "sgesvd.f"
		sorgbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2000 "sgesvd.f"
	    }
#line 2001 "sgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 2007 "sgesvd.f"
		slacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2008 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2008 "sgesvd.f"
		sorgbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2010 "sgesvd.f"
	    }
#line 2011 "sgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N, prefer 3*N+N*NB) */

#line 2017 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2017 "sgesvd.f"
		sorgbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2019 "sgesvd.f"
	    }
#line 2020 "sgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 2026 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2026 "sgesvd.f"
		sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2028 "sgesvd.f"
	    }
#line 2029 "sgesvd.f"
	    iwork = ie + *n;
#line 2030 "sgesvd.f"
	    if (wntuas || wntuo) {
#line 2030 "sgesvd.f"
		nru = *m;
#line 2030 "sgesvd.f"
	    }
#line 2032 "sgesvd.f"
	    if (wntun) {
#line 2032 "sgesvd.f"
		nru = 0;
#line 2032 "sgesvd.f"
	    }
#line 2034 "sgesvd.f"
	    if (wntvas || wntvo) {
#line 2034 "sgesvd.f"
		ncvt = *n;
#line 2034 "sgesvd.f"
	    }
#line 2036 "sgesvd.f"
	    if (wntvn) {
#line 2036 "sgesvd.f"
		ncvt = 0;
#line 2036 "sgesvd.f"
	    }
#line 2038 "sgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2045 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2047 "sgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 2054 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 2056 "sgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 2063 "sgesvd.f"
		sbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 2065 "sgesvd.f"
	    }

#line 2067 "sgesvd.f"
	}

#line 2069 "sgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2075 "sgesvd.f"
	if (*n >= mnthr) {

#line 2077 "sgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2082 "sgesvd.f"
		itau = 1;
#line 2083 "sgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need 2*M, prefer M+M*NB) */

#line 2088 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2088 "sgesvd.f"
		sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2093 "sgesvd.f"
		i__2 = *m - 1;
#line 2093 "sgesvd.f"
		i__3 = *m - 1;
#line 2093 "sgesvd.f"
		slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 2094 "sgesvd.f"
		ie = 1;
#line 2095 "sgesvd.f"
		itauq = ie + *m;
#line 2096 "sgesvd.f"
		itaup = itauq + *m;
#line 2097 "sgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2102 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2102 "sgesvd.f"
		sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2105 "sgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2110 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2110 "sgesvd.f"
		    sorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2112 "sgesvd.f"
		}
#line 2113 "sgesvd.f"
		iwork = ie + *m;
#line 2114 "sgesvd.f"
		nru = 0;
#line 2115 "sgesvd.f"
		if (wntuo || wntuas) {
#line 2115 "sgesvd.f"
		    nru = *m;
#line 2115 "sgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

#line 2122 "sgesvd.f"
		sbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &work[ie], dum, &
			c__1, &a[a_offset], lda, dum, &c__1, &work[iwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2127 "sgesvd.f"
		if (wntuas) {
#line 2127 "sgesvd.f"
		    slacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2127 "sgesvd.f"
		}

#line 2130 "sgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

/* Computing MAX */
#line 2136 "sgesvd.f"
		i__2 = *m << 2;
#line 2136 "sgesvd.f"
		if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2140 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2141 "sgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2141 "sgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2145 "sgesvd.f"
			ldwrku = *lda;
#line 2146 "sgesvd.f"
			chunk = *n;
#line 2147 "sgesvd.f"
			ldwrkr = *lda;
#line 2148 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2148 "sgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n + *m;
#line 2148 "sgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2152 "sgesvd.f"
			    ldwrku = *lda;
#line 2153 "sgesvd.f"
			    chunk = *n;
#line 2154 "sgesvd.f"
			    ldwrkr = *m;
#line 2155 "sgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2159 "sgesvd.f"
			    ldwrku = *m;
#line 2160 "sgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2161 "sgesvd.f"
			    ldwrkr = *m;
#line 2162 "sgesvd.f"
			}
#line 2162 "sgesvd.f"
		    }
#line 2163 "sgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2164 "sgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2169 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2169 "sgesvd.f"
		    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2174 "sgesvd.f"
		    slacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2175 "sgesvd.f"
		    i__2 = *m - 1;
#line 2175 "sgesvd.f"
		    i__3 = *m - 1;
#line 2175 "sgesvd.f"
		    slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2181 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2181 "sgesvd.f"
		    sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2183 "sgesvd.f"
		    ie = itau;
#line 2184 "sgesvd.f"
		    itauq = ie + *m;
#line 2185 "sgesvd.f"
		    itaup = itauq + *m;
#line 2186 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2191 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2191 "sgesvd.f"
		    sgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */

#line 2198 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2198 "sgesvd.f"
		    sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2201 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M+BDSPAC) */

#line 2207 "sgesvd.f"
		    sbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[
			    ir], &ldwrkr, dum, &c__1, dum, &c__1, &work[iwork]
			    , info, (ftnlen)1);
#line 2210 "sgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M*N+M) */

#line 2216 "sgesvd.f"
		    i__2 = *n;
#line 2216 "sgesvd.f"
		    i__3 = chunk;
#line 2216 "sgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2217 "sgesvd.f"
			i__4 = *n - i__ + 1;
#line 2217 "sgesvd.f"
			blk = min(i__4,chunk);
#line 2218 "sgesvd.f"
			sgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2221 "sgesvd.f"
			slacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2223 "sgesvd.f"
/* L30: */
#line 2223 "sgesvd.f"
		    }

#line 2225 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2229 "sgesvd.f"
		    ie = 1;
#line 2230 "sgesvd.f"
		    itauq = ie + *m;
#line 2231 "sgesvd.f"
		    itaup = itauq + *m;
#line 2232 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 2237 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2237 "sgesvd.f"
		    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2244 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2244 "sgesvd.f"
		    sorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2246 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2252 "sgesvd.f"
		    sbdsqr_("L", m, n, &c__0, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, dum, &c__1, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);

#line 2255 "sgesvd.f"
		}

#line 2257 "sgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

/* Computing MAX */
#line 2263 "sgesvd.f"
		i__3 = *m << 2;
#line 2263 "sgesvd.f"
		if (*lwork >= *m * *m + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

#line 2267 "sgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2268 "sgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2268 "sgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2272 "sgesvd.f"
			ldwrku = *lda;
#line 2273 "sgesvd.f"
			chunk = *n;
#line 2274 "sgesvd.f"
			ldwrkr = *lda;
#line 2275 "sgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2275 "sgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n + *m;
#line 2275 "sgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2279 "sgesvd.f"
			    ldwrku = *lda;
#line 2280 "sgesvd.f"
			    chunk = *n;
#line 2281 "sgesvd.f"
			    ldwrkr = *m;
#line 2282 "sgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2286 "sgesvd.f"
			    ldwrku = *m;
#line 2287 "sgesvd.f"
			    chunk = (*lwork - *m * *m - *m) / *m;
#line 2288 "sgesvd.f"
			    ldwrkr = *m;
#line 2289 "sgesvd.f"
			}
#line 2289 "sgesvd.f"
		    }
#line 2290 "sgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2291 "sgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2296 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2296 "sgesvd.f"
		    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2301 "sgesvd.f"
		    slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2302 "sgesvd.f"
		    i__3 = *m - 1;
#line 2302 "sgesvd.f"
		    i__2 = *m - 1;
#line 2302 "sgesvd.f"
		    slaset_("U", &i__3, &i__2, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2308 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2308 "sgesvd.f"
		    sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2310 "sgesvd.f"
		    ie = itau;
#line 2311 "sgesvd.f"
		    itauq = ie + *m;
#line 2312 "sgesvd.f"
		    itaup = itauq + *m;
#line 2313 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2318 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2318 "sgesvd.f"
		    sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2321 "sgesvd.f"
		    slacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */

#line 2326 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2326 "sgesvd.f"
		    sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 2333 "sgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2333 "sgesvd.f"
		    sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2335 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M+BDSPAC) */

#line 2342 "sgesvd.f"
		    sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[ir], 
			    &ldwrkr, &u[u_offset], ldu, dum, &c__1, &work[
			    iwork], info, (ftnlen)1);
#line 2345 "sgesvd.f"
		    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M+2*M, prefer M*M+M*N+M)) */

#line 2351 "sgesvd.f"
		    i__3 = *n;
#line 2351 "sgesvd.f"
		    i__2 = chunk;
#line 2351 "sgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2352 "sgesvd.f"
			i__4 = *n - i__ + 1;
#line 2352 "sgesvd.f"
			blk = min(i__4,chunk);
#line 2353 "sgesvd.f"
			sgemm_("N", "N", m, &blk, m, &c_b79, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2356 "sgesvd.f"
			slacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2358 "sgesvd.f"
/* L40: */
#line 2358 "sgesvd.f"
		    }

#line 2360 "sgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2364 "sgesvd.f"
		    itau = 1;
#line 2365 "sgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need 2*M, prefer M+M*NB) */

#line 2370 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2370 "sgesvd.f"
		    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2375 "sgesvd.f"
		    slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2376 "sgesvd.f"
		    i__2 = *m - 1;
#line 2376 "sgesvd.f"
		    i__3 = *m - 1;
#line 2376 "sgesvd.f"
		    slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 << 
			    1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need 2*M, prefer M+M*NB) */

#line 2382 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2382 "sgesvd.f"
		    sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2384 "sgesvd.f"
		    ie = itau;
#line 2385 "sgesvd.f"
		    itauq = ie + *m;
#line 2386 "sgesvd.f"
		    itaup = itauq + *m;
#line 2387 "sgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2392 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2392 "sgesvd.f"
		    sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2399 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2399 "sgesvd.f"
		    sormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2406 "sgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2406 "sgesvd.f"
		    sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2408 "sgesvd.f"
		    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

#line 2415 "sgesvd.f"
		    sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, dum, &c__1, &
			    work[iwork], info, (ftnlen)1);

#line 2418 "sgesvd.f"
		}

#line 2420 "sgesvd.f"
	    } else if (wntvs) {

#line 2422 "sgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2428 "sgesvd.f"
		    i__2 = *m << 2;
#line 2428 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2432 "sgesvd.f"
			ir = 1;
#line 2433 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2437 "sgesvd.f"
			    ldwrkr = *lda;
#line 2438 "sgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2442 "sgesvd.f"
			    ldwrkr = *m;
#line 2443 "sgesvd.f"
			}
#line 2444 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2445 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2450 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2450 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2455 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2457 "sgesvd.f"
			i__2 = *m - 1;
#line 2457 "sgesvd.f"
			i__3 = *m - 1;
#line 2457 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2463 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2463 "sgesvd.f"
			sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2465 "sgesvd.f"
			ie = itau;
#line 2466 "sgesvd.f"
			itauq = ie + *m;
#line 2467 "sgesvd.f"
			itaup = itauq + *m;
#line 2468 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2473 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2473 "sgesvd.f"
			sgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB) */

#line 2482 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2482 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2485 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2491 "sgesvd.f"
			sbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2499 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2502 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2506 "sgesvd.f"
			itau = 1;
#line 2507 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2512 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2512 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2517 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2522 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2522 "sgesvd.f"
			sorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2524 "sgesvd.f"
			ie = itau;
#line 2525 "sgesvd.f"
			itauq = ie + *m;
#line 2526 "sgesvd.f"
			itaup = itauq + *m;
#line 2527 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2531 "sgesvd.f"
			i__2 = *m - 1;
#line 2531 "sgesvd.f"
			i__3 = *m - 1;
#line 2531 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2537 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2537 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2544 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2544 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2547 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2553 "sgesvd.f"
			sbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 2557 "sgesvd.f"
		    }

#line 2559 "sgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 2565 "sgesvd.f"
		    i__2 = *m << 2;
#line 2565 "sgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2569 "sgesvd.f"
			iu = 1;
#line 2570 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2574 "sgesvd.f"
			    ldwrku = *lda;
#line 2575 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2576 "sgesvd.f"
			    ldwrkr = *lda;
#line 2577 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2581 "sgesvd.f"
			    ldwrku = *lda;
#line 2582 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2583 "sgesvd.f"
			    ldwrkr = *m;
#line 2584 "sgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2588 "sgesvd.f"
			    ldwrku = *m;
#line 2589 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2590 "sgesvd.f"
			    ldwrkr = *m;
#line 2591 "sgesvd.f"
			}
#line 2592 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2593 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 2598 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2598 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2603 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2605 "sgesvd.f"
			i__2 = *m - 1;
#line 2605 "sgesvd.f"
			i__3 = *m - 1;
#line 2605 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 2611 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2611 "sgesvd.f"
			sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2613 "sgesvd.f"
			ie = itau;
#line 2614 "sgesvd.f"
			itauq = ie + *m;
#line 2615 "sgesvd.f"
			itaup = itauq + *m;
#line 2616 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 2623 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2623 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2627 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M+4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 2634 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2634 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB) */

#line 2641 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2641 "sgesvd.f"
			sorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2644 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M+BDSPAC) */

#line 2651 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2659 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (Workspace: need M*M) */

#line 2665 "sgesvd.f"
			slacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2668 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2672 "sgesvd.f"
			itau = 1;
#line 2673 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2678 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2678 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2680 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2685 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2685 "sgesvd.f"
			sorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2687 "sgesvd.f"
			ie = itau;
#line 2688 "sgesvd.f"
			itauq = ie + *m;
#line 2689 "sgesvd.f"
			itaup = itauq + *m;
#line 2690 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2694 "sgesvd.f"
			i__2 = *m - 1;
#line 2694 "sgesvd.f"
			i__3 = *m - 1;
#line 2694 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2700 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2700 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2707 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2707 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2714 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2714 "sgesvd.f"
			sorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2716 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, compute left */
/*                    singular vectors of A in A and compute right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2723 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2727 "sgesvd.f"
		    }

#line 2729 "sgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 2736 "sgesvd.f"
		    i__2 = *m << 2;
#line 2736 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2740 "sgesvd.f"
			iu = 1;
#line 2741 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2745 "sgesvd.f"
			    ldwrku = *lda;
#line 2746 "sgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2750 "sgesvd.f"
			    ldwrku = *m;
#line 2751 "sgesvd.f"
			}
#line 2752 "sgesvd.f"
			itau = iu + ldwrku * *m;
#line 2753 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2758 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2758 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2763 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2765 "sgesvd.f"
			i__2 = *m - 1;
#line 2765 "sgesvd.f"
			i__3 = *m - 1;
#line 2765 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2771 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2771 "sgesvd.f"
			sorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2773 "sgesvd.f"
			ie = itau;
#line 2774 "sgesvd.f"
			itauq = ie + *m;
#line 2775 "sgesvd.f"
			itaup = itauq + *m;
#line 2776 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2781 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2781 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2785 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M+4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2792 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2792 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 2799 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2799 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2801 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2808 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

#line 2816 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&a[a_offset], lda, &c_b57, &vt[vt_offset], 
				ldvt, (ftnlen)1, (ftnlen)1);

#line 2819 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2823 "sgesvd.f"
			itau = 1;
#line 2824 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2829 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2829 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2831 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2836 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2836 "sgesvd.f"
			sorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2841 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2842 "sgesvd.f"
			i__2 = *m - 1;
#line 2842 "sgesvd.f"
			i__3 = *m - 1;
#line 2842 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 2844 "sgesvd.f"
			ie = itau;
#line 2845 "sgesvd.f"
			itauq = ie + *m;
#line 2846 "sgesvd.f"
			itaup = itauq + *m;
#line 2847 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 2852 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2852 "sgesvd.f"
			sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 2860 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2860 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 2867 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2867 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2869 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 2876 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 2880 "sgesvd.f"
		    }

#line 2882 "sgesvd.f"
		}

#line 2884 "sgesvd.f"
	    } else if (wntva) {

#line 2886 "sgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 2892 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 2892 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 2896 "sgesvd.f"
			ir = 1;
#line 2897 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2901 "sgesvd.f"
			    ldwrkr = *lda;
#line 2902 "sgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2906 "sgesvd.f"
			    ldwrkr = *m;
#line 2907 "sgesvd.f"
			}
#line 2908 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2909 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 2914 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2914 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2916 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2920 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2922 "sgesvd.f"
			i__2 = *m - 1;
#line 2922 "sgesvd.f"
			i__3 = *m - 1;
#line 2922 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB) */

#line 2928 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2928 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2930 "sgesvd.f"
			ie = itau;
#line 2931 "sgesvd.f"
			itauq = ie + *m;
#line 2932 "sgesvd.f"
			itaup = itauq + *m;
#line 2933 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 2938 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2938 "sgesvd.f"
			sgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need M*M+4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

#line 2947 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2947 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2950 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 2956 "sgesvd.f"
			sbdsqr_("U", m, m, &c__0, &c__0, &s[1], &work[ie], &
				work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 2964 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[ir], &ldwrkr, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 2969 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 2971 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2975 "sgesvd.f"
			itau = 1;
#line 2976 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 2981 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2981 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2983 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 2988 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2988 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2990 "sgesvd.f"
			ie = itau;
#line 2991 "sgesvd.f"
			itauq = ie + *m;
#line 2992 "sgesvd.f"
			itaup = itauq + *m;
#line 2993 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2997 "sgesvd.f"
			i__2 = *m - 1;
#line 2997 "sgesvd.f"
			i__3 = *m - 1;
#line 2997 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 3003 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3003 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3011 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3011 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3014 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3020 "sgesvd.f"
			sbdsqr_("U", m, n, &c__0, &c__0, &s[1], &work[ie], &
				vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

#line 3024 "sgesvd.f"
		    }

#line 3026 "sgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3032 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3032 "sgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3036 "sgesvd.f"
			iu = 1;
#line 3037 "sgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3041 "sgesvd.f"
			    ldwrku = *lda;
#line 3042 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3043 "sgesvd.f"
			    ldwrkr = *lda;
#line 3044 "sgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3048 "sgesvd.f"
			    ldwrku = *lda;
#line 3049 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3050 "sgesvd.f"
			    ldwrkr = *m;
#line 3051 "sgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3055 "sgesvd.f"
			    ldwrku = *m;
#line 3056 "sgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3057 "sgesvd.f"
			    ldwrkr = *m;
#line 3058 "sgesvd.f"
			}
#line 3059 "sgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3060 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */

#line 3065 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3065 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3067 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */

#line 3072 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3072 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3077 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3079 "sgesvd.f"
			i__2 = *m - 1;
#line 3079 "sgesvd.f"
			i__3 = *m - 1;
#line 3079 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3081 "sgesvd.f"
			ie = itau;
#line 3082 "sgesvd.f"
			itauq = ie + *m;
#line 3083 "sgesvd.f"
			itaup = itauq + *m;
#line 3084 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

#line 3091 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3091 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3095 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M+4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

#line 3102 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3102 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB) */

#line 3109 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3109 "sgesvd.f"
			sorgbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3112 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M+BDSPAC) */

#line 3119 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1, 
				&work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3127 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3132 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3136 "sgesvd.f"
			slacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3139 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3143 "sgesvd.f"
			itau = 1;
#line 3144 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 3149 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3149 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3151 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 3156 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3156 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3158 "sgesvd.f"
			ie = itau;
#line 3159 "sgesvd.f"
			itauq = ie + *m;
#line 3160 "sgesvd.f"
			itaup = itauq + *m;
#line 3161 "sgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3165 "sgesvd.f"
			i__2 = *m - 1;
#line 3165 "sgesvd.f"
			i__3 = *m - 1;
#line 3165 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 
				<< 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 3171 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3171 "sgesvd.f"
			sgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3179 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3179 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3186 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3186 "sgesvd.f"
			sorgbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3188 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3195 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3199 "sgesvd.f"
		    }

#line 3201 "sgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3208 "sgesvd.f"
		    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
#line 3208 "sgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3212 "sgesvd.f"
			iu = 1;
#line 3213 "sgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3217 "sgesvd.f"
			    ldwrku = *lda;
#line 3218 "sgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3222 "sgesvd.f"
			    ldwrku = *m;
#line 3223 "sgesvd.f"
			}
#line 3224 "sgesvd.f"
			itau = iu + ldwrku * *m;
#line 3225 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB) */

#line 3230 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3230 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3232 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB) */

#line 3237 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3237 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3242 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3244 "sgesvd.f"
			i__2 = *m - 1;
#line 3244 "sgesvd.f"
			i__3 = *m - 1;
#line 3244 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3246 "sgesvd.f"
			ie = itau;
#line 3247 "sgesvd.f"
			itauq = ie + *m;
#line 3248 "sgesvd.f"
			itaup = itauq + *m;
#line 3249 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */

#line 3254 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3254 "sgesvd.f"
			sgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3258 "sgesvd.f"
			slacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB) */

#line 3264 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3264 "sgesvd.f"
			sorgbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB) */

#line 3271 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3271 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3273 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M+BDSPAC) */

#line 3280 "sgesvd.f"
			sbdsqr_("U", m, m, m, &c__0, &s[1], &work[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
				work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

#line 3288 "sgesvd.f"
			sgemm_("N", "N", m, n, m, &c_b79, &work[iu], &ldwrku, 
				&vt[vt_offset], ldvt, &c_b57, &a[a_offset], 
				lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3293 "sgesvd.f"
			slacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3295 "sgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3299 "sgesvd.f"
			itau = 1;
#line 3300 "sgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M+M*NB) */

#line 3305 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3305 "sgesvd.f"
			sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3307 "sgesvd.f"
			slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M+N, prefer M+N*NB) */

#line 3312 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3312 "sgesvd.f"
			sorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3317 "sgesvd.f"
			slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3318 "sgesvd.f"
			i__2 = *m - 1;
#line 3318 "sgesvd.f"
			i__3 = *m - 1;
#line 3318 "sgesvd.f"
			slaset_("U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 
				<< 1) + 1], ldu, (ftnlen)1);
#line 3320 "sgesvd.f"
			ie = itau;
#line 3321 "sgesvd.f"
			itauq = ie + *m;
#line 3322 "sgesvd.f"
			itaup = itauq + *m;
#line 3323 "sgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M+2*M*NB) */

#line 3328 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3328 "sgesvd.f"
			sgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M+N, prefer 3*M+N*NB) */

#line 3336 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3336 "sgesvd.f"
			sormbr_("P", "L", "T", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3343 "sgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3343 "sgesvd.f"
			sorgbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3345 "sgesvd.f"
			iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

#line 3352 "sgesvd.f"
			sbdsqr_("U", m, n, m, &c__0, &s[1], &work[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, dum, &
				c__1, &work[iwork], info, (ftnlen)1);

#line 3356 "sgesvd.f"
		    }

#line 3358 "sgesvd.f"
		}

#line 3360 "sgesvd.f"
	    }

#line 3362 "sgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3369 "sgesvd.f"
	    ie = 1;
#line 3370 "sgesvd.f"
	    itauq = ie + *m;
#line 3371 "sgesvd.f"
	    itaup = itauq + *m;
#line 3372 "sgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 3377 "sgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3377 "sgesvd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__2, &ierr);
#line 3380 "sgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) */

#line 3386 "sgesvd.f"
		slacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3387 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3387 "sgesvd.f"
		sorgbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3389 "sgesvd.f"
	    }
#line 3390 "sgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB) */

#line 3396 "sgesvd.f"
		slacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3397 "sgesvd.f"
		if (wntva) {
#line 3397 "sgesvd.f"
		    nrvt = *n;
#line 3397 "sgesvd.f"
		}
#line 3399 "sgesvd.f"
		if (wntvs) {
#line 3399 "sgesvd.f"
		    nrvt = *m;
#line 3399 "sgesvd.f"
		}
#line 3401 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3401 "sgesvd.f"
		sorgbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3403 "sgesvd.f"
	    }
#line 3404 "sgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB) */

#line 3410 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3410 "sgesvd.f"
		sorgbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3412 "sgesvd.f"
	    }
#line 3413 "sgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 3419 "sgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3419 "sgesvd.f"
		sorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3421 "sgesvd.f"
	    }
#line 3422 "sgesvd.f"
	    iwork = ie + *m;
#line 3423 "sgesvd.f"
	    if (wntuas || wntuo) {
#line 3423 "sgesvd.f"
		nru = *m;
#line 3423 "sgesvd.f"
	    }
#line 3425 "sgesvd.f"
	    if (wntun) {
#line 3425 "sgesvd.f"
		nru = 0;
#line 3425 "sgesvd.f"
	    }
#line 3427 "sgesvd.f"
	    if (wntvas || wntvo) {
#line 3427 "sgesvd.f"
		ncvt = *n;
#line 3427 "sgesvd.f"
	    }
#line 3429 "sgesvd.f"
	    if (wntvn) {
#line 3429 "sgesvd.f"
		ncvt = 0;
#line 3429 "sgesvd.f"
	    }
#line 3431 "sgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3438 "sgesvd.f"
		sbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3440 "sgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

#line 3447 "sgesvd.f"
		sbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
			a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
			iwork], info, (ftnlen)1);
#line 3449 "sgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

#line 3456 "sgesvd.f"
		sbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
			work[iwork], info, (ftnlen)1);
#line 3458 "sgesvd.f"
	    }

#line 3460 "sgesvd.f"
	}

#line 3462 "sgesvd.f"
    }

/*     If SBDSQR failed to converge, copy unconverged superdiagonals */
/*     to WORK( 2:MINMN ) */

#line 3467 "sgesvd.f"
    if (*info != 0) {
#line 3468 "sgesvd.f"
	if (ie > 2) {
#line 3469 "sgesvd.f"
	    i__2 = minmn - 1;
#line 3469 "sgesvd.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 3470 "sgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3471 "sgesvd.f"
/* L50: */
#line 3471 "sgesvd.f"
	    }
#line 3472 "sgesvd.f"
	}
#line 3473 "sgesvd.f"
	if (ie < 2) {
#line 3474 "sgesvd.f"
	    for (i__ = minmn - 1; i__ >= 1; --i__) {
#line 3475 "sgesvd.f"
		work[i__ + 1] = work[i__ + ie - 1];
#line 3476 "sgesvd.f"
/* L60: */
#line 3476 "sgesvd.f"
	    }
#line 3477 "sgesvd.f"
	}
#line 3478 "sgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3482 "sgesvd.f"
    if (iscl == 1) {
#line 3483 "sgesvd.f"
	if (anrm > bignum) {
#line 3483 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3483 "sgesvd.f"
	}
#line 3486 "sgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3486 "sgesvd.f"
	    i__2 = minmn - 1;
#line 3486 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3486 "sgesvd.f"
	}
#line 3489 "sgesvd.f"
	if (anrm < smlnum) {
#line 3489 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3489 "sgesvd.f"
	}
#line 3492 "sgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3492 "sgesvd.f"
	    i__2 = minmn - 1;
#line 3492 "sgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &work[2],
		     &minmn, &ierr, (ftnlen)1);
#line 3492 "sgesvd.f"
	}
#line 3495 "sgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3499 "sgesvd.f"
    work[1] = (doublereal) maxwrk;

#line 3501 "sgesvd.f"
    return 0;

/*     End of SGESVD */

} /* sgesvd_ */

