#line 1 "cgesvd.f"
/* cgesvd.f -- translated by f2c (version 20100827).
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

#line 1 "cgesvd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__1 = 1;

/* > \brief <b> CGESVD computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU, JOBVT */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), S( * ) */
/*       COMPLEX            A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESVD computes the singular value decomposition (SVD) of a complex */
/* > M-by-N matrix A, optionally computing the left and/or right singular */
/* > vectors. The SVD is written */
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
/* > Note that the routine returns V**H, not V. */
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
/* >          V**H: */
/* >          = 'A':  all N rows of V**H are returned in the array VT; */
/* >          = 'S':  the first min(m,n) rows of V**H (the right singular */
/* >                  vectors) are returned in the array VT; */
/* >          = 'O':  the first min(m,n) rows of V**H (the right singular */
/* >                  vectors) are overwritten on the array A; */
/* >          = 'N':  no rows of V**H (no right singular vectors) are */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          if JOBU = 'O',  A is overwritten with the first min(m,n) */
/* >                          columns of U (the left singular vectors, */
/* >                          stored columnwise); */
/* >          if JOBVT = 'O', A is overwritten with the first min(m,n) */
/* >                          rows of V**H (the right singular vectors, */
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
/* >          U is COMPLEX array, dimension (LDU,UCOL) */
/* >          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'. */
/* >          If JOBU = 'A', U contains the M-by-M unitary matrix U; */
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
/* >          VT is COMPLEX array, dimension (LDVT,N) */
/* >          If JOBVT = 'A', VT contains the N-by-N unitary matrix */
/* >          V**H; */
/* >          if JOBVT = 'S', VT contains the first min(m,n) rows of */
/* >          V**H (the right singular vectors, stored rowwise); */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N)). */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (5*min(M,N)) */
/* >          On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the */
/* >          unconverged superdiagonal elements of an upper bidiagonal */
/* >          matrix B whose diagonal is in S (not necessarily sorted). */
/* >          B satisfies A = U * B * VT, so it has the same singular */
/* >          values as A, and singular vectors related by U and VT. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if CBDSQR did not converge, INFO specifies how many */
/* >                superdiagonals of an intermediate bidiagonal form B */
/* >                did not converge to zero. See the description of RWORK */
/* >                above for details. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup complexGEsing */

/*  ===================================================================== */
/* Subroutine */ int cgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, 
	integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info, ftnlen jobu_len, 
	ftnlen jobvt_len)
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
    static integer nru;
    static doublecomplex cdum[1];
    static integer iscl;
    static doublereal anrm;
    static integer ierr, itau, ncvt, nrvt, lwork_cgebrd__, lwork_cgelqf__, 
	    lwork_cgeqrf__;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn, wrkbl, itaup, itauq, mnthr, iwork;
    static logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    extern /* Subroutine */ int cgebrd_();
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(), clascl_(char *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen), cgeqrf_();
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), cbdsqr_(
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), cungbr_();
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int cunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen), cunglq_(), cungqr_();
    static integer ldwrkr, minwrk, ldwrku, maxwrk;
    static doublereal smlnum;
    static integer irwork;
    static logical lquery, wntuas, wntvas;
    static integer lwork_cungbr_p__, lwork_cungbr_q__, lwork_cunglq_m__, 
	    lwork_cunglq_n__, lwork_cungqr_m__, lwork_cungqr_n__;


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

#line 275 "cgesvd.f"
    /* Parameter adjustments */
#line 275 "cgesvd.f"
    a_dim1 = *lda;
#line 275 "cgesvd.f"
    a_offset = 1 + a_dim1;
#line 275 "cgesvd.f"
    a -= a_offset;
#line 275 "cgesvd.f"
    --s;
#line 275 "cgesvd.f"
    u_dim1 = *ldu;
#line 275 "cgesvd.f"
    u_offset = 1 + u_dim1;
#line 275 "cgesvd.f"
    u -= u_offset;
#line 275 "cgesvd.f"
    vt_dim1 = *ldvt;
#line 275 "cgesvd.f"
    vt_offset = 1 + vt_dim1;
#line 275 "cgesvd.f"
    vt -= vt_offset;
#line 275 "cgesvd.f"
    --work;
#line 275 "cgesvd.f"
    --rwork;
#line 275 "cgesvd.f"

#line 275 "cgesvd.f"
    /* Function Body */
#line 275 "cgesvd.f"
    *info = 0;
#line 276 "cgesvd.f"
    minmn = min(*m,*n);
#line 277 "cgesvd.f"
    wntua = lsame_(jobu, "A", (ftnlen)1, (ftnlen)1);
#line 278 "cgesvd.f"
    wntus = lsame_(jobu, "S", (ftnlen)1, (ftnlen)1);
#line 279 "cgesvd.f"
    wntuas = wntua || wntus;
#line 280 "cgesvd.f"
    wntuo = lsame_(jobu, "O", (ftnlen)1, (ftnlen)1);
#line 281 "cgesvd.f"
    wntun = lsame_(jobu, "N", (ftnlen)1, (ftnlen)1);
#line 282 "cgesvd.f"
    wntva = lsame_(jobvt, "A", (ftnlen)1, (ftnlen)1);
#line 283 "cgesvd.f"
    wntvs = lsame_(jobvt, "S", (ftnlen)1, (ftnlen)1);
#line 284 "cgesvd.f"
    wntvas = wntva || wntvs;
#line 285 "cgesvd.f"
    wntvo = lsame_(jobvt, "O", (ftnlen)1, (ftnlen)1);
#line 286 "cgesvd.f"
    wntvn = lsame_(jobvt, "N", (ftnlen)1, (ftnlen)1);
#line 287 "cgesvd.f"
    lquery = *lwork == -1;

#line 289 "cgesvd.f"
    if (! (wntua || wntus || wntuo || wntun)) {
#line 290 "cgesvd.f"
	*info = -1;
#line 291 "cgesvd.f"
    } else if (! (wntva || wntvs || wntvo || wntvn) || wntvo && wntuo) {
#line 293 "cgesvd.f"
	*info = -2;
#line 294 "cgesvd.f"
    } else if (*m < 0) {
#line 295 "cgesvd.f"
	*info = -3;
#line 296 "cgesvd.f"
    } else if (*n < 0) {
#line 297 "cgesvd.f"
	*info = -4;
#line 298 "cgesvd.f"
    } else if (*lda < max(1,*m)) {
#line 299 "cgesvd.f"
	*info = -6;
#line 300 "cgesvd.f"
    } else if (*ldu < 1 || wntuas && *ldu < *m) {
#line 301 "cgesvd.f"
	*info = -9;
#line 302 "cgesvd.f"
    } else if (*ldvt < 1 || wntva && *ldvt < *n || wntvs && *ldvt < minmn) {
#line 304 "cgesvd.f"
	*info = -11;
#line 305 "cgesvd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to */
/*       real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 315 "cgesvd.f"
    if (*info == 0) {
#line 316 "cgesvd.f"
	minwrk = 1;
#line 317 "cgesvd.f"
	maxwrk = 1;
#line 318 "cgesvd.f"
	if (*m >= *n && minmn > 0) {

/*           Space needed for ZBDSQR is BDSPAC = 5*N */

/* Writing concatenation */
#line 322 "cgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 322 "cgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 322 "cgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 322 "cgesvd.f"
	    mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
/*           Compute space needed for CGEQRF */
#line 324 "cgesvd.f"
	    cgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 325 "cgesvd.f"
	    lwork_cgeqrf__ = (integer) dum[0];
/*           Compute space needed for CUNGQR */
#line 327 "cgesvd.f"
	    cungqr_(m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 328 "cgesvd.f"
	    lwork_cungqr_n__ = (integer) dum[0];
#line 329 "cgesvd.f"
	    cungqr_(m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 330 "cgesvd.f"
	    lwork_cungqr_m__ = (integer) dum[0];
/*           Compute space needed for CGEBRD */
#line 332 "cgesvd.f"
	    cgebrd_(n, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 334 "cgesvd.f"
	    lwork_cgebrd__ = (integer) dum[0];
/*           Compute space needed for CUNGBR */
#line 336 "cgesvd.f"
	    cungbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 338 "cgesvd.f"
	    lwork_cungbr_p__ = (integer) dum[0];
#line 339 "cgesvd.f"
	    cungbr_("Q", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 341 "cgesvd.f"
	    lwork_cungbr_q__ = (integer) dum[0];

/* Writing concatenation */
#line 343 "cgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 343 "cgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 343 "cgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 343 "cgesvd.f"
	    mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
#line 344 "cgesvd.f"
	    if (*m >= mnthr) {
#line 345 "cgesvd.f"
		if (wntun) {

/*                 Path 1 (M much larger than N, JOBU='N') */

#line 349 "cgesvd.f"
		    maxwrk = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 350 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_cgebrd__;
#line 350 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 351 "cgesvd.f"
		    if (wntvo || wntvas) {
/* Computing MAX */
#line 351 "cgesvd.f"
			i__2 = maxwrk, i__3 = (*n << 1) + lwork_cungbr_p__;
#line 351 "cgesvd.f"
			maxwrk = max(i__2,i__3);
#line 351 "cgesvd.f"
		    }
#line 353 "cgesvd.f"
		    minwrk = *n * 3;
#line 354 "cgesvd.f"
		} else if (wntuo && wntvn) {

/*                 Path 2 (M much larger than N, JOBU='O', JOBVT='N') */

#line 358 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 359 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_n__;
#line 359 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 360 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 360 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 361 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 361 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 362 "cgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n;
#line 362 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 363 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 364 "cgesvd.f"
		} else if (wntuo && wntvas) {

/*                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or */
/*                 'A') */

#line 369 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 370 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_n__;
#line 370 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 371 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 371 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 372 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 372 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 373 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_p__;
#line 373 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 374 "cgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n;
#line 374 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 375 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 376 "cgesvd.f"
		} else if (wntus && wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */

#line 380 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 381 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_n__;
#line 381 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 382 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 382 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 383 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 383 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 384 "cgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 385 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 386 "cgesvd.f"
		} else if (wntus && wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */

#line 390 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 391 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_n__;
#line 391 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 392 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 392 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 393 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 393 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 394 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_p__;
#line 394 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 395 "cgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
#line 396 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 397 "cgesvd.f"
		} else if (wntus && wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or */
/*                 'A') */

#line 402 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 403 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_n__;
#line 403 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 404 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 404 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 405 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 405 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 406 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_p__;
#line 406 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 407 "cgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 408 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 409 "cgesvd.f"
		} else if (wntua && wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */

#line 413 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 414 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_m__;
#line 414 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 415 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 415 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 416 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 416 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 417 "cgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 418 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 419 "cgesvd.f"
		} else if (wntua && wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */

#line 423 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 424 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_m__;
#line 424 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 425 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 425 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 426 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 426 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 427 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_p__;
#line 427 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 428 "cgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
#line 429 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 430 "cgesvd.f"
		} else if (wntua && wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or */
/*                 'A') */

#line 435 "cgesvd.f"
		    wrkbl = *n + lwork_cgeqrf__;
/* Computing MAX */
#line 436 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_cungqr_m__;
#line 436 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 437 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cgebrd__;
#line 437 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 438 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 438 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 439 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_cungbr_p__;
#line 439 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 440 "cgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 441 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 442 "cgesvd.f"
		}
#line 443 "cgesvd.f"
	    } else {

/*              Path 10 (M at least N, but not much larger) */

#line 447 "cgesvd.f"
		cgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 449 "cgesvd.f"
		lwork_cgebrd__ = (integer) dum[0];
#line 450 "cgesvd.f"
		maxwrk = (*n << 1) + lwork_cgebrd__;
#line 451 "cgesvd.f"
		if (wntus || wntuo) {
#line 452 "cgesvd.f"
		    cungbr_("Q", m, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 454 "cgesvd.f"
		    lwork_cungbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 455 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 455 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 456 "cgesvd.f"
		}
#line 457 "cgesvd.f"
		if (wntua) {
#line 458 "cgesvd.f"
		    cungbr_("Q", m, m, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 460 "cgesvd.f"
		    lwork_cungbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 461 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_cungbr_q__;
#line 461 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 462 "cgesvd.f"
		}
#line 463 "cgesvd.f"
		if (! wntvn) {
/* Computing MAX */
#line 464 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_cungbr_p__;
#line 464 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 465 "cgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 466 "cgesvd.f"
		}
#line 467 "cgesvd.f"
	    }
#line 468 "cgesvd.f"
	} else if (minmn > 0) {

/*           Space needed for CBDSQR is BDSPAC = 5*M */

/* Writing concatenation */
#line 472 "cgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 472 "cgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 472 "cgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 472 "cgesvd.f"
	    mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
/*           Compute space needed for CGELQF */
#line 474 "cgesvd.f"
	    cgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 475 "cgesvd.f"
	    lwork_cgelqf__ = (integer) dum[0];
/*           Compute space needed for CUNGLQ */
#line 477 "cgesvd.f"
	    cunglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
#line 478 "cgesvd.f"
	    lwork_cunglq_n__ = (integer) dum[0];
#line 479 "cgesvd.f"
	    cunglq_(m, n, m, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 480 "cgesvd.f"
	    lwork_cunglq_m__ = (integer) dum[0];
/*           Compute space needed for CGEBRD */
#line 482 "cgesvd.f"
	    cgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 484 "cgesvd.f"
	    lwork_cgebrd__ = (integer) dum[0];
/*            Compute space needed for CUNGBR P */
#line 486 "cgesvd.f"
	    cungbr_("P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 488 "cgesvd.f"
	    lwork_cungbr_p__ = (integer) dum[0];
/*           Compute space needed for CUNGBR Q */
#line 490 "cgesvd.f"
	    cungbr_("Q", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 492 "cgesvd.f"
	    lwork_cungbr_q__ = (integer) dum[0];
#line 493 "cgesvd.f"
	    if (*n >= mnthr) {
#line 494 "cgesvd.f"
		if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N') */

#line 498 "cgesvd.f"
		    maxwrk = *m + lwork_cgelqf__;
/* Computing MAX */
#line 499 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cgebrd__;
#line 499 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 500 "cgesvd.f"
		    if (wntuo || wntuas) {
/* Computing MAX */
#line 500 "cgesvd.f"
			i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 500 "cgesvd.f"
			maxwrk = max(i__2,i__3);
#line 500 "cgesvd.f"
		    }
#line 502 "cgesvd.f"
		    minwrk = *m * 3;
#line 503 "cgesvd.f"
		} else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O') */

#line 507 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 508 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 508 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 509 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 509 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 510 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 510 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 511 "cgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 511 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 512 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 513 "cgesvd.f"
		} else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='O') */

#line 518 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 519 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 519 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 520 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 520 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 521 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 521 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 522 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 522 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 523 "cgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 523 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 524 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 525 "cgesvd.f"
		} else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */

#line 529 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 530 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 530 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 531 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 531 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 532 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 532 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 533 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 534 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 535 "cgesvd.f"
		} else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */

#line 539 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 540 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 540 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 541 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 541 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 542 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 542 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 543 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 543 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 544 "cgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 545 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 546 "cgesvd.f"
		} else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='S') */

#line 551 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 552 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 552 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 553 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 553 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 554 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 554 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 555 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 555 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 556 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 557 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 558 "cgesvd.f"
		} else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */

#line 562 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 563 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_n__;
#line 563 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 564 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 564 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 565 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 565 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 566 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 567 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 568 "cgesvd.f"
		} else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */

#line 572 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 573 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_n__;
#line 573 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 574 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 574 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 575 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 575 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 576 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 576 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 577 "cgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 578 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 579 "cgesvd.f"
		} else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='A') */

#line 584 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 585 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_n__;
#line 585 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 586 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 586 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 587 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 587 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 588 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 588 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 589 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 590 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 591 "cgesvd.f"
		}
#line 592 "cgesvd.f"
	    } else {

/*              Path 10t(N greater than M, but not much larger) */

#line 596 "cgesvd.f"
		cgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 598 "cgesvd.f"
		lwork_cgebrd__ = (integer) dum[0];
#line 599 "cgesvd.f"
		maxwrk = (*m << 1) + lwork_cgebrd__;
#line 600 "cgesvd.f"
		if (wntvs || wntvo) {
/*                Compute space needed for CUNGBR P */
#line 602 "cgesvd.f"
		    cungbr_("P", m, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 604 "cgesvd.f"
		    lwork_cungbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 605 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 605 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 606 "cgesvd.f"
		}
#line 607 "cgesvd.f"
		if (wntva) {
#line 608 "cgesvd.f"
		    cungbr_("P", n, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 610 "cgesvd.f"
		    lwork_cungbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 611 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 611 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 612 "cgesvd.f"
		}
#line 613 "cgesvd.f"
		if (! wntun) {
/* Computing MAX */
#line 614 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 614 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 615 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 616 "cgesvd.f"
		}
#line 617 "cgesvd.f"
	    }
#line 618 "cgesvd.f"
	}
#line 619 "cgesvd.f"
	maxwrk = max(minwrk,maxwrk);
#line 620 "cgesvd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 622 "cgesvd.f"
	if (*lwork < minwrk && ! lquery) {
#line 623 "cgesvd.f"
	    *info = -13;
#line 624 "cgesvd.f"
	}
#line 625 "cgesvd.f"
    }

#line 627 "cgesvd.f"
    if (*info != 0) {
#line 628 "cgesvd.f"
	i__2 = -(*info);
#line 628 "cgesvd.f"
	xerbla_("CGESVD", &i__2, (ftnlen)6);
#line 629 "cgesvd.f"
	return 0;
#line 630 "cgesvd.f"
    } else if (lquery) {
#line 631 "cgesvd.f"
	return 0;
#line 632 "cgesvd.f"
    }

/*     Quick return if possible */

#line 636 "cgesvd.f"
    if (*m == 0 || *n == 0) {
#line 637 "cgesvd.f"
	return 0;
#line 638 "cgesvd.f"
    }

/*     Get machine constants */

#line 642 "cgesvd.f"
    eps = slamch_("P", (ftnlen)1);
#line 643 "cgesvd.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 644 "cgesvd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 648 "cgesvd.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 649 "cgesvd.f"
    iscl = 0;
#line 650 "cgesvd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 651 "cgesvd.f"
	iscl = 1;
#line 652 "cgesvd.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 653 "cgesvd.f"
    } else if (anrm > bignum) {
#line 654 "cgesvd.f"
	iscl = 1;
#line 655 "cgesvd.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 656 "cgesvd.f"
    }

#line 658 "cgesvd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 664 "cgesvd.f"
	if (*m >= mnthr) {

#line 666 "cgesvd.f"
	    if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N') */
/*              No left singular vectors to be computed */

#line 671 "cgesvd.f"
		itau = 1;
#line 672 "cgesvd.f"
		iwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: need 0) */

#line 678 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 678 "cgesvd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out below R */

#line 683 "cgesvd.f"
		i__2 = *n - 1;
#line 683 "cgesvd.f"
		i__3 = *n - 1;
#line 683 "cgesvd.f"
		claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 685 "cgesvd.f"
		ie = 1;
#line 686 "cgesvd.f"
		itauq = 1;
#line 687 "cgesvd.f"
		itaup = itauq + *n;
#line 688 "cgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 694 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 694 "cgesvd.f"
		cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 697 "cgesvd.f"
		ncvt = 0;
#line 698 "cgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 704 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 704 "cgesvd.f"
		    cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 706 "cgesvd.f"
		    ncvt = *n;
#line 707 "cgesvd.f"
		}
#line 708 "cgesvd.f"
		irwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 715 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			irwork], info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 720 "cgesvd.f"
		if (wntvas) {
#line 720 "cgesvd.f"
		    clacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 720 "cgesvd.f"
		}

#line 723 "cgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

#line 729 "cgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 733 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 734 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 734 "cgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 738 "cgesvd.f"
			ldwrku = *lda;
#line 739 "cgesvd.f"
			ldwrkr = *lda;
#line 740 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 740 "cgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 740 "cgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 744 "cgesvd.f"
			    ldwrku = *lda;
#line 745 "cgesvd.f"
			    ldwrkr = *n;
#line 746 "cgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 750 "cgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 751 "cgesvd.f"
			    ldwrkr = *n;
#line 752 "cgesvd.f"
			}
#line 752 "cgesvd.f"
		    }
#line 753 "cgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 754 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 760 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 760 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 765 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 766 "cgesvd.f"
		    i__2 = *n - 1;
#line 766 "cgesvd.f"
		    i__3 = *n - 1;
#line 766 "cgesvd.f"
		    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1], &
			    ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 773 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 773 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 775 "cgesvd.f"
		    ie = 1;
#line 776 "cgesvd.f"
		    itauq = itau;
#line 777 "cgesvd.f"
		    itaup = itauq + *n;
#line 778 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 784 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 784 "cgesvd.f"
		    cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: need 0) */

#line 792 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 792 "cgesvd.f"
		    cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 795 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 802 "cgesvd.f"
		    cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &work[ir], &ldwrkr, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 805 "cgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 812 "cgesvd.f"
		    i__2 = *m;
#line 812 "cgesvd.f"
		    i__3 = ldwrku;
#line 812 "cgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 813 "cgesvd.f"
			i__4 = *m - i__ + 1;
#line 813 "cgesvd.f"
			chunk = min(i__4,ldwrku);
#line 814 "cgesvd.f"
			cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 817 "cgesvd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 819 "cgesvd.f"
/* L10: */
#line 819 "cgesvd.f"
		    }

#line 821 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 825 "cgesvd.f"
		    ie = 1;
#line 826 "cgesvd.f"
		    itauq = 1;
#line 827 "cgesvd.f"
		    itaup = itauq + *n;
#line 828 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*                 (RWorkspace: N) */

#line 834 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 834 "cgesvd.f"
		    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 842 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 842 "cgesvd.f"
		    cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 844 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (CWorkspace: need 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 851 "cgesvd.f"
		    cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &a[a_offset], lda, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 854 "cgesvd.f"
		}

#line 856 "cgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 862 "cgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 866 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 867 "cgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 867 "cgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 871 "cgesvd.f"
			ldwrku = *lda;
#line 872 "cgesvd.f"
			ldwrkr = *lda;
#line 873 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 873 "cgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 873 "cgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 877 "cgesvd.f"
			    ldwrku = *lda;
#line 878 "cgesvd.f"
			    ldwrkr = *n;
#line 879 "cgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 883 "cgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 884 "cgesvd.f"
			    ldwrkr = *n;
#line 885 "cgesvd.f"
			}
#line 885 "cgesvd.f"
		    }
#line 886 "cgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 887 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 893 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 893 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 898 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 899 "cgesvd.f"
		    if (*n > 1) {
#line 899 "cgesvd.f"
			i__3 = *n - 1;
#line 899 "cgesvd.f"
			i__2 = *n - 1;
#line 899 "cgesvd.f"
			claset_("L", &i__3, &i__2, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 899 "cgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 907 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 907 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 909 "cgesvd.f"
		    ie = 1;
#line 910 "cgesvd.f"
		    itauq = itau;
#line 911 "cgesvd.f"
		    itaup = itauq + *n;
#line 912 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 918 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 918 "cgesvd.f"
		    cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 921 "cgesvd.f"
		    clacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 927 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 927 "cgesvd.f"
		    cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 935 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 935 "cgesvd.f"
		    cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 937 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 945 "cgesvd.f"
		    cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, cdum, &c__1,
			     &rwork[irwork], info, (ftnlen)1);
#line 948 "cgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 955 "cgesvd.f"
		    i__3 = *m;
#line 955 "cgesvd.f"
		    i__2 = ldwrku;
#line 955 "cgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 956 "cgesvd.f"
			i__4 = *m - i__ + 1;
#line 956 "cgesvd.f"
			chunk = min(i__4,ldwrku);
#line 957 "cgesvd.f"
			cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 960 "cgesvd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 962 "cgesvd.f"
/* L20: */
#line 962 "cgesvd.f"
		    }

#line 964 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 968 "cgesvd.f"
		    itau = 1;
#line 969 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 975 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 975 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 980 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 981 "cgesvd.f"
		    if (*n > 1) {
#line 981 "cgesvd.f"
			i__2 = *n - 1;
#line 981 "cgesvd.f"
			i__3 = *n - 1;
#line 981 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 981 "cgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 989 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 989 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 991 "cgesvd.f"
		    ie = 1;
#line 992 "cgesvd.f"
		    itauq = itau;
#line 993 "cgesvd.f"
		    itaup = itauq + *n;
#line 994 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                 (RWorkspace: N) */

#line 1000 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1000 "cgesvd.f"
		    cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                 (RWorkspace: 0) */

#line 1008 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1008 "cgesvd.f"
		    cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 1016 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1016 "cgesvd.f"
		    cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1018 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 1026 "cgesvd.f"
		    cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, 
			    &rwork[irwork], info, (ftnlen)1);

#line 1030 "cgesvd.f"
		}

#line 1032 "cgesvd.f"
	    } else if (wntus) {

#line 1034 "cgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

#line 1040 "cgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1044 "cgesvd.f"
			ir = 1;
#line 1045 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1049 "cgesvd.f"
			    ldwrkr = *lda;
#line 1050 "cgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1054 "cgesvd.f"
			    ldwrkr = *n;
#line 1055 "cgesvd.f"
			}
#line 1056 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1057 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1063 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1063 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1068 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1070 "cgesvd.f"
			i__2 = *n - 1;
#line 1070 "cgesvd.f"
			i__3 = *n - 1;
#line 1070 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1077 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1077 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1079 "cgesvd.f"
			ie = 1;
#line 1080 "cgesvd.f"
			itauq = itau;
#line 1081 "cgesvd.f"
			itaup = itauq + *n;
#line 1082 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1088 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1088 "cgesvd.f"
			cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1097 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1097 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1100 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1107 "cgesvd.f"
			cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1116 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1119 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1123 "cgesvd.f"
			itau = 1;
#line 1124 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1130 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1130 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1132 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1138 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1138 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1140 "cgesvd.f"
			ie = 1;
#line 1141 "cgesvd.f"
			itauq = itau;
#line 1142 "cgesvd.f"
			itaup = itauq + *n;
#line 1143 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1147 "cgesvd.f"
			i__2 = *n - 1;
#line 1147 "cgesvd.f"
			i__3 = *n - 1;
#line 1147 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1154 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1154 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1162 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1162 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1165 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1172 "cgesvd.f"
			cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1176 "cgesvd.f"
		    }

#line 1178 "cgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

#line 1184 "cgesvd.f"
		    if (*lwork >= (*n << 1) * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1188 "cgesvd.f"
			iu = 1;
#line 1189 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1193 "cgesvd.f"
			    ldwrku = *lda;
#line 1194 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1195 "cgesvd.f"
			    ldwrkr = *lda;
#line 1196 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1200 "cgesvd.f"
			    ldwrku = *lda;
#line 1201 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1202 "cgesvd.f"
			    ldwrkr = *n;
#line 1203 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1207 "cgesvd.f"
			    ldwrku = *n;
#line 1208 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1209 "cgesvd.f"
			    ldwrkr = *n;
#line 1210 "cgesvd.f"
			}
#line 1211 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1212 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1218 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1218 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1223 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1225 "cgesvd.f"
			i__2 = *n - 1;
#line 1225 "cgesvd.f"
			i__3 = *n - 1;
#line 1225 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1232 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1232 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1234 "cgesvd.f"
			ie = 1;
#line 1235 "cgesvd.f"
			itauq = itau;
#line 1236 "cgesvd.f"
			itaup = itauq + *n;
#line 1237 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1245 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1245 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1249 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1256 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1256 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1265 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1265 "cgesvd.f"
			cungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1268 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1276 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1286 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1293 "cgesvd.f"
			clacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1296 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1300 "cgesvd.f"
			itau = 1;
#line 1301 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1307 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1307 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1309 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1315 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1315 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1317 "cgesvd.f"
			ie = 1;
#line 1318 "cgesvd.f"
			itauq = itau;
#line 1319 "cgesvd.f"
			itaup = itauq + *n;
#line 1320 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1324 "cgesvd.f"
			i__2 = *n - 1;
#line 1324 "cgesvd.f"
			i__3 = *n - 1;
#line 1324 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1331 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1331 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1339 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1339 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1347 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1347 "cgesvd.f"
			cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1349 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1357 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1361 "cgesvd.f"
		    }

#line 1363 "cgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

#line 1370 "cgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1374 "cgesvd.f"
			iu = 1;
#line 1375 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1379 "cgesvd.f"
			    ldwrku = *lda;
#line 1380 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1384 "cgesvd.f"
			    ldwrku = *n;
#line 1385 "cgesvd.f"
			}
#line 1386 "cgesvd.f"
			itau = iu + ldwrku * *n;
#line 1387 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1393 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1393 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1398 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1400 "cgesvd.f"
			i__2 = *n - 1;
#line 1400 "cgesvd.f"
			i__3 = *n - 1;
#line 1400 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1407 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1407 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1409 "cgesvd.f"
			ie = 1;
#line 1410 "cgesvd.f"
			itauq = itau;
#line 1411 "cgesvd.f"
			itaup = itauq + *n;
#line 1412 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1418 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1418 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1422 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1429 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1429 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1438 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1438 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1440 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1448 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1457 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1460 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1464 "cgesvd.f"
			itau = 1;
#line 1465 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1471 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1471 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1473 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1479 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1479 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1484 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1485 "cgesvd.f"
			if (*n > 1) {
#line 1485 "cgesvd.f"
			    i__2 = *n - 1;
#line 1485 "cgesvd.f"
			    i__3 = *n - 1;
#line 1485 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1485 "cgesvd.f"
			}
#line 1488 "cgesvd.f"
			ie = 1;
#line 1489 "cgesvd.f"
			itauq = itau;
#line 1490 "cgesvd.f"
			itaup = itauq + *n;
#line 1491 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1497 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1497 "cgesvd.f"
			cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1506 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1506 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1514 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1514 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1516 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1524 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1528 "cgesvd.f"
		    }

#line 1530 "cgesvd.f"
		}

#line 1532 "cgesvd.f"
	    } else if (wntua) {

#line 1534 "cgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1540 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1540 "cgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1544 "cgesvd.f"
			ir = 1;
#line 1545 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1549 "cgesvd.f"
			    ldwrkr = *lda;
#line 1550 "cgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1554 "cgesvd.f"
			    ldwrkr = *n;
#line 1555 "cgesvd.f"
			}
#line 1556 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1557 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1563 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1563 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1565 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1569 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1571 "cgesvd.f"
			i__2 = *n - 1;
#line 1571 "cgesvd.f"
			i__3 = *n - 1;
#line 1571 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1578 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1578 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1580 "cgesvd.f"
			ie = 1;
#line 1581 "cgesvd.f"
			itauq = itau;
#line 1582 "cgesvd.f"
			itaup = itauq + *n;
#line 1583 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1589 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1589 "cgesvd.f"
			cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1598 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1598 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1601 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1608 "cgesvd.f"
			cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1617 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1622 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1624 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1628 "cgesvd.f"
			itau = 1;
#line 1629 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1635 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1635 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1637 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1643 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1643 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1645 "cgesvd.f"
			ie = 1;
#line 1646 "cgesvd.f"
			itauq = itau;
#line 1647 "cgesvd.f"
			itaup = itauq + *n;
#line 1648 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1652 "cgesvd.f"
			i__2 = *n - 1;
#line 1652 "cgesvd.f"
			i__3 = *n - 1;
#line 1652 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1659 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1659 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1668 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1668 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1671 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1678 "cgesvd.f"
			cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1682 "cgesvd.f"
		    }

#line 1684 "cgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1690 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1690 "cgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1694 "cgesvd.f"
			iu = 1;
#line 1695 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1699 "cgesvd.f"
			    ldwrku = *lda;
#line 1700 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1701 "cgesvd.f"
			    ldwrkr = *lda;
#line 1702 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1706 "cgesvd.f"
			    ldwrku = *lda;
#line 1707 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1708 "cgesvd.f"
			    ldwrkr = *n;
#line 1709 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1713 "cgesvd.f"
			    ldwrku = *n;
#line 1714 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1715 "cgesvd.f"
			    ldwrkr = *n;
#line 1716 "cgesvd.f"
			}
#line 1717 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1718 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1724 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1724 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1726 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1732 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1732 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1737 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1739 "cgesvd.f"
			i__2 = *n - 1;
#line 1739 "cgesvd.f"
			i__3 = *n - 1;
#line 1739 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1741 "cgesvd.f"
			ie = 1;
#line 1742 "cgesvd.f"
			itauq = itau;
#line 1743 "cgesvd.f"
			itaup = itauq + *n;
#line 1744 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1752 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1752 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1756 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1763 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1763 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1772 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1772 "cgesvd.f"
			cungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1775 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1783 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1793 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1798 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1802 "cgesvd.f"
			clacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1805 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1809 "cgesvd.f"
			itau = 1;
#line 1810 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1816 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1816 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1818 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1824 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1824 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1826 "cgesvd.f"
			ie = 1;
#line 1827 "cgesvd.f"
			itauq = itau;
#line 1828 "cgesvd.f"
			itaup = itauq + *n;
#line 1829 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1833 "cgesvd.f"
			i__2 = *n - 1;
#line 1833 "cgesvd.f"
			i__3 = *n - 1;
#line 1833 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1840 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1840 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1849 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1849 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1857 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1857 "cgesvd.f"
			cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1859 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1867 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1871 "cgesvd.f"
		    }

#line 1873 "cgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1880 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1880 "cgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1884 "cgesvd.f"
			iu = 1;
#line 1885 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1889 "cgesvd.f"
			    ldwrku = *lda;
#line 1890 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1894 "cgesvd.f"
			    ldwrku = *n;
#line 1895 "cgesvd.f"
			}
#line 1896 "cgesvd.f"
			itau = iu + ldwrku * *n;
#line 1897 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1903 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1903 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1905 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1911 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1911 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1916 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1918 "cgesvd.f"
			i__2 = *n - 1;
#line 1918 "cgesvd.f"
			i__3 = *n - 1;
#line 1918 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1920 "cgesvd.f"
			ie = 1;
#line 1921 "cgesvd.f"
			itauq = itau;
#line 1922 "cgesvd.f"
			itaup = itauq + *n;
#line 1923 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1929 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1929 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1933 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1940 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1940 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: need   0) */

#line 1949 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1949 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1951 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1959 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1968 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1973 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1975 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1979 "cgesvd.f"
			itau = 1;
#line 1980 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1986 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1986 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1988 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1994 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1994 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 1999 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 2000 "cgesvd.f"
			if (*n > 1) {
#line 2000 "cgesvd.f"
			    i__2 = *n - 1;
#line 2000 "cgesvd.f"
			    i__3 = *n - 1;
#line 2000 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 2000 "cgesvd.f"
			}
#line 2003 "cgesvd.f"
			ie = 1;
#line 2004 "cgesvd.f"
			itauq = itau;
#line 2005 "cgesvd.f"
			itaup = itauq + *n;
#line 2006 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 2012 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2012 "cgesvd.f"
			cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2021 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2021 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2029 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2029 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 2031 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2039 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2043 "cgesvd.f"
		    }

#line 2045 "cgesvd.f"
		}

#line 2047 "cgesvd.f"
	    }

#line 2049 "cgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 2056 "cgesvd.f"
	    ie = 1;
#line 2057 "cgesvd.f"
	    itauq = 1;
#line 2058 "cgesvd.f"
	    itaup = itauq + *n;
#line 2059 "cgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 2065 "cgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 2065 "cgesvd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 2068 "cgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB) */
/*              (RWorkspace: 0) */

#line 2075 "cgesvd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 2076 "cgesvd.f"
		if (wntus) {
#line 2076 "cgesvd.f"
		    ncu = *n;
#line 2076 "cgesvd.f"
		}
#line 2078 "cgesvd.f"
		if (wntua) {
#line 2078 "cgesvd.f"
		    ncu = *m;
#line 2078 "cgesvd.f"
		}
#line 2080 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2080 "cgesvd.f"
		cungbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2082 "cgesvd.f"
	    }
#line 2083 "cgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2090 "cgesvd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2091 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2091 "cgesvd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2093 "cgesvd.f"
	    }
#line 2094 "cgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 2101 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2101 "cgesvd.f"
		cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2103 "cgesvd.f"
	    }
#line 2104 "cgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2111 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2111 "cgesvd.f"
		cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2113 "cgesvd.f"
	    }
#line 2114 "cgesvd.f"
	    irwork = ie + *n;
#line 2115 "cgesvd.f"
	    if (wntuas || wntuo) {
#line 2115 "cgesvd.f"
		nru = *m;
#line 2115 "cgesvd.f"
	    }
#line 2117 "cgesvd.f"
	    if (wntun) {
#line 2117 "cgesvd.f"
		nru = 0;
#line 2117 "cgesvd.f"
	    }
#line 2119 "cgesvd.f"
	    if (wntvas || wntvo) {
#line 2119 "cgesvd.f"
		ncvt = *n;
#line 2119 "cgesvd.f"
	    }
#line 2121 "cgesvd.f"
	    if (wntvn) {
#line 2121 "cgesvd.f"
		ncvt = 0;
#line 2121 "cgesvd.f"
	    }
#line 2123 "cgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2131 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2134 "cgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2142 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2145 "cgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2153 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2156 "cgesvd.f"
	    }

#line 2158 "cgesvd.f"
	}

#line 2160 "cgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2166 "cgesvd.f"
	if (*n >= mnthr) {

#line 2168 "cgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2173 "cgesvd.f"
		itau = 1;
#line 2174 "cgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 2180 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2180 "cgesvd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2185 "cgesvd.f"
		i__2 = *m - 1;
#line 2185 "cgesvd.f"
		i__3 = *m - 1;
#line 2185 "cgesvd.f"
		claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 2187 "cgesvd.f"
		ie = 1;
#line 2188 "cgesvd.f"
		itauq = 1;
#line 2189 "cgesvd.f"
		itaup = itauq + *m;
#line 2190 "cgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 2196 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2196 "cgesvd.f"
		cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2199 "cgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2205 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2205 "cgesvd.f"
		    cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2207 "cgesvd.f"
		}
#line 2208 "cgesvd.f"
		irwork = ie + *m;
#line 2209 "cgesvd.f"
		nru = 0;
#line 2210 "cgesvd.f"
		if (wntuo || wntuas) {
#line 2210 "cgesvd.f"
		    nru = *m;
#line 2210 "cgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2218 "cgesvd.f"
		cbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &rwork[ie], cdum, &
			c__1, &a[a_offset], lda, cdum, &c__1, &rwork[irwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2223 "cgesvd.f"
		if (wntuas) {
#line 2223 "cgesvd.f"
		    clacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2223 "cgesvd.f"
		}

#line 2226 "cgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

#line 2232 "cgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2236 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2237 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 2237 "cgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2241 "cgesvd.f"
			ldwrku = *lda;
#line 2242 "cgesvd.f"
			chunk = *n;
#line 2243 "cgesvd.f"
			ldwrkr = *lda;
#line 2244 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2244 "cgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 2244 "cgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2248 "cgesvd.f"
			    ldwrku = *lda;
#line 2249 "cgesvd.f"
			    chunk = *n;
#line 2250 "cgesvd.f"
			    ldwrkr = *m;
#line 2251 "cgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2255 "cgesvd.f"
			    ldwrku = *m;
#line 2256 "cgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2257 "cgesvd.f"
			    ldwrkr = *m;
#line 2258 "cgesvd.f"
			}
#line 2258 "cgesvd.f"
		    }
#line 2259 "cgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2260 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2266 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2266 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2271 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2272 "cgesvd.f"
		    i__2 = *m - 1;
#line 2272 "cgesvd.f"
		    i__3 = *m - 1;
#line 2272 "cgesvd.f"
		    claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2279 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2279 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2281 "cgesvd.f"
		    ie = 1;
#line 2282 "cgesvd.f"
		    itauq = itau;
#line 2283 "cgesvd.f"
		    itaup = itauq + *m;
#line 2284 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2290 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2290 "cgesvd.f"
		    cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2298 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2298 "cgesvd.f"
		    cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2301 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2308 "cgesvd.f"
		    cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &work[
			    ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2311 "cgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N) */
/*                 (RWorkspace: 0) */

#line 2318 "cgesvd.f"
		    i__2 = *n;
#line 2318 "cgesvd.f"
		    i__3 = chunk;
#line 2318 "cgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2319 "cgesvd.f"
			i__4 = *n - i__ + 1;
#line 2319 "cgesvd.f"
			blk = min(i__4,chunk);
#line 2320 "cgesvd.f"
			cgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2323 "cgesvd.f"
			clacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2325 "cgesvd.f"
/* L30: */
#line 2325 "cgesvd.f"
		    }

#line 2327 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2331 "cgesvd.f"
		    ie = 1;
#line 2332 "cgesvd.f"
		    itauq = 1;
#line 2333 "cgesvd.f"
		    itaup = itauq + *m;
#line 2334 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*                 (RWorkspace: need M) */

#line 2340 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2340 "cgesvd.f"
		    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2348 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2348 "cgesvd.f"
		    cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2350 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2357 "cgesvd.f"
		    cbdsqr_("L", m, n, &c__0, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 2360 "cgesvd.f"
		}

#line 2362 "cgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 2368 "cgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2372 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2373 "cgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 2373 "cgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2377 "cgesvd.f"
			ldwrku = *lda;
#line 2378 "cgesvd.f"
			chunk = *n;
#line 2379 "cgesvd.f"
			ldwrkr = *lda;
#line 2380 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2380 "cgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 2380 "cgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2384 "cgesvd.f"
			    ldwrku = *lda;
#line 2385 "cgesvd.f"
			    chunk = *n;
#line 2386 "cgesvd.f"
			    ldwrkr = *m;
#line 2387 "cgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2391 "cgesvd.f"
			    ldwrku = *m;
#line 2392 "cgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2393 "cgesvd.f"
			    ldwrkr = *m;
#line 2394 "cgesvd.f"
			}
#line 2394 "cgesvd.f"
		    }
#line 2395 "cgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2396 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2402 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2402 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2407 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2408 "cgesvd.f"
		    i__3 = *m - 1;
#line 2408 "cgesvd.f"
		    i__2 = *m - 1;
#line 2408 "cgesvd.f"
		    claset_("U", &i__3, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2415 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2415 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2417 "cgesvd.f"
		    ie = 1;
#line 2418 "cgesvd.f"
		    itauq = itau;
#line 2419 "cgesvd.f"
		    itaup = itauq + *m;
#line 2420 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2426 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2426 "cgesvd.f"
		    cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2429 "cgesvd.f"
		    clacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2435 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2435 "cgesvd.f"
		    cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2443 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2443 "cgesvd.f"
		    cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2445 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2453 "cgesvd.f"
		    cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[ir],
			     &ldwrkr, &u[u_offset], ldu, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2456 "cgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N)) */
/*                 (RWorkspace: 0) */

#line 2463 "cgesvd.f"
		    i__3 = *n;
#line 2463 "cgesvd.f"
		    i__2 = chunk;
#line 2463 "cgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2464 "cgesvd.f"
			i__4 = *n - i__ + 1;
#line 2464 "cgesvd.f"
			blk = min(i__4,chunk);
#line 2465 "cgesvd.f"
			cgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2468 "cgesvd.f"
			clacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2470 "cgesvd.f"
/* L40: */
#line 2470 "cgesvd.f"
		    }

#line 2472 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2476 "cgesvd.f"
		    itau = 1;
#line 2477 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2483 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2483 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2488 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2489 "cgesvd.f"
		    i__2 = *m - 1;
#line 2489 "cgesvd.f"
		    i__3 = *m - 1;
#line 2489 "cgesvd.f"
		    claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2496 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2496 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2498 "cgesvd.f"
		    ie = 1;
#line 2499 "cgesvd.f"
		    itauq = itau;
#line 2500 "cgesvd.f"
		    itaup = itauq + *m;
#line 2501 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2507 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2507 "cgesvd.f"
		    cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                 (RWorkspace: 0) */

#line 2515 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2515 "cgesvd.f"
		    cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2523 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2523 "cgesvd.f"
		    cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2525 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2533 "cgesvd.f"
		    cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			    rwork[irwork], info, (ftnlen)1);

#line 2536 "cgesvd.f"
		}

#line 2538 "cgesvd.f"
	    } else if (wntvs) {

#line 2540 "cgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

#line 2546 "cgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2550 "cgesvd.f"
			ir = 1;
#line 2551 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2555 "cgesvd.f"
			    ldwrkr = *lda;
#line 2556 "cgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2560 "cgesvd.f"
			    ldwrkr = *m;
#line 2561 "cgesvd.f"
			}
#line 2562 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2563 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2569 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2569 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2574 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2576 "cgesvd.f"
			i__2 = *m - 1;
#line 2576 "cgesvd.f"
			i__3 = *m - 1;
#line 2576 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2583 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2583 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2585 "cgesvd.f"
			ie = 1;
#line 2586 "cgesvd.f"
			itauq = itau;
#line 2587 "cgesvd.f"
			itaup = itauq + *m;
#line 2588 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2594 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2594 "cgesvd.f"
			cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2604 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2604 "cgesvd.f"
			cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2607 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2614 "cgesvd.f"
			cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2623 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2626 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2630 "cgesvd.f"
			itau = 1;
#line 2631 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2637 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2637 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2642 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2648 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2648 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2650 "cgesvd.f"
			ie = 1;
#line 2651 "cgesvd.f"
			itauq = itau;
#line 2652 "cgesvd.f"
			itaup = itauq + *m;
#line 2653 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2657 "cgesvd.f"
			i__2 = *m - 1;
#line 2657 "cgesvd.f"
			i__3 = *m - 1;
#line 2657 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2664 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2664 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2672 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2672 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2675 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2682 "cgesvd.f"
			cbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 2686 "cgesvd.f"
		    }

#line 2688 "cgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

#line 2694 "cgesvd.f"
		    if (*lwork >= (*m << 1) * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2698 "cgesvd.f"
			iu = 1;
#line 2699 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2703 "cgesvd.f"
			    ldwrku = *lda;
#line 2704 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2705 "cgesvd.f"
			    ldwrkr = *lda;
#line 2706 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2710 "cgesvd.f"
			    ldwrku = *lda;
#line 2711 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2712 "cgesvd.f"
			    ldwrkr = *m;
#line 2713 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2717 "cgesvd.f"
			    ldwrku = *m;
#line 2718 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2719 "cgesvd.f"
			    ldwrkr = *m;
#line 2720 "cgesvd.f"
			}
#line 2721 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2722 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2728 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2728 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2733 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2735 "cgesvd.f"
			i__2 = *m - 1;
#line 2735 "cgesvd.f"
			i__3 = *m - 1;
#line 2735 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2742 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2742 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2744 "cgesvd.f"
			ie = 1;
#line 2745 "cgesvd.f"
			itauq = itau;
#line 2746 "cgesvd.f"
			itaup = itauq + *m;
#line 2747 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 2755 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2755 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2759 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2767 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2767 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2775 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2775 "cgesvd.f"
			cungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2778 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2786 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2796 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2803 "cgesvd.f"
			clacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2806 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2810 "cgesvd.f"
			itau = 1;
#line 2811 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2817 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2817 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2819 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2825 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2825 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2827 "cgesvd.f"
			ie = 1;
#line 2828 "cgesvd.f"
			itauq = itau;
#line 2829 "cgesvd.f"
			itaup = itauq + *m;
#line 2830 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2834 "cgesvd.f"
			i__2 = *m - 1;
#line 2834 "cgesvd.f"
			i__3 = *m - 1;
#line 2834 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2841 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2841 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2849 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2849 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2857 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2857 "cgesvd.f"
			cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2859 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2867 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2871 "cgesvd.f"
		    }

#line 2873 "cgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

#line 2880 "cgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2884 "cgesvd.f"
			iu = 1;
#line 2885 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2889 "cgesvd.f"
			    ldwrku = *lda;
#line 2890 "cgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2894 "cgesvd.f"
			    ldwrku = *m;
#line 2895 "cgesvd.f"
			}
#line 2896 "cgesvd.f"
			itau = iu + ldwrku * *m;
#line 2897 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2903 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2903 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2908 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2910 "cgesvd.f"
			i__2 = *m - 1;
#line 2910 "cgesvd.f"
			i__3 = *m - 1;
#line 2910 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2917 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2917 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2919 "cgesvd.f"
			ie = 1;
#line 2920 "cgesvd.f"
			itauq = itau;
#line 2921 "cgesvd.f"
			itaup = itauq + *m;
#line 2922 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2928 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2928 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2932 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2940 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2940 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2948 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2948 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2950 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2958 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2967 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2970 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2974 "cgesvd.f"
			itau = 1;
#line 2975 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2981 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2981 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2983 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2989 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2989 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2994 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2995 "cgesvd.f"
			i__2 = *m - 1;
#line 2995 "cgesvd.f"
			i__3 = *m - 1;
#line 2995 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 2997 "cgesvd.f"
			ie = 1;
#line 2998 "cgesvd.f"
			itauq = itau;
#line 2999 "cgesvd.f"
			itaup = itauq + *m;
#line 3000 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3006 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3006 "cgesvd.f"
			cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3015 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3015 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3023 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3023 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3025 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3033 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3037 "cgesvd.f"
		    }

#line 3039 "cgesvd.f"
		}

#line 3041 "cgesvd.f"
	    } else if (wntva) {

#line 3043 "cgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 3049 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3049 "cgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3053 "cgesvd.f"
			ir = 1;
#line 3054 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 3058 "cgesvd.f"
			    ldwrkr = *lda;
#line 3059 "cgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 3063 "cgesvd.f"
			    ldwrkr = *m;
#line 3064 "cgesvd.f"
			}
#line 3065 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3066 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3072 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3072 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3074 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 3078 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 3080 "cgesvd.f"
			i__2 = *m - 1;
#line 3080 "cgesvd.f"
			i__3 = *m - 1;
#line 3080 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3087 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3087 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3089 "cgesvd.f"
			ie = 1;
#line 3090 "cgesvd.f"
			itauq = itau;
#line 3091 "cgesvd.f"
			itaup = itauq + *m;
#line 3092 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3098 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3098 "cgesvd.f"
			cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3108 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3108 "cgesvd.f"
			cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3111 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3118 "cgesvd.f"
			cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3127 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3132 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3134 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3138 "cgesvd.f"
			itau = 1;
#line 3139 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3145 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3145 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3147 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3153 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3153 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3155 "cgesvd.f"
			ie = 1;
#line 3156 "cgesvd.f"
			itauq = itau;
#line 3157 "cgesvd.f"
			itaup = itauq + *m;
#line 3158 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3162 "cgesvd.f"
			i__2 = *m - 1;
#line 3162 "cgesvd.f"
			i__3 = *m - 1;
#line 3162 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3169 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3169 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3178 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3178 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3181 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3188 "cgesvd.f"
			cbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 3192 "cgesvd.f"
		    }

#line 3194 "cgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3200 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3200 "cgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3204 "cgesvd.f"
			iu = 1;
#line 3205 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3209 "cgesvd.f"
			    ldwrku = *lda;
#line 3210 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3211 "cgesvd.f"
			    ldwrkr = *lda;
#line 3212 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3216 "cgesvd.f"
			    ldwrku = *lda;
#line 3217 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3218 "cgesvd.f"
			    ldwrkr = *m;
#line 3219 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3223 "cgesvd.f"
			    ldwrku = *m;
#line 3224 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3225 "cgesvd.f"
			    ldwrkr = *m;
#line 3226 "cgesvd.f"
			}
#line 3227 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3228 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3234 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3234 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3236 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3242 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3242 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3247 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3249 "cgesvd.f"
			i__2 = *m - 1;
#line 3249 "cgesvd.f"
			i__3 = *m - 1;
#line 3249 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3251 "cgesvd.f"
			ie = 1;
#line 3252 "cgesvd.f"
			itauq = itau;
#line 3253 "cgesvd.f"
			itaup = itauq + *m;
#line 3254 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 3262 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3262 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3266 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3274 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3274 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3282 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3282 "cgesvd.f"
			cungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3285 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3293 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3303 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3308 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3312 "cgesvd.f"
			clacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3315 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3319 "cgesvd.f"
			itau = 1;
#line 3320 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3326 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3326 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3328 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3334 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3334 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3336 "cgesvd.f"
			ie = 1;
#line 3337 "cgesvd.f"
			itauq = itau;
#line 3338 "cgesvd.f"
			itaup = itauq + *m;
#line 3339 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3343 "cgesvd.f"
			i__2 = *m - 1;
#line 3343 "cgesvd.f"
			i__3 = *m - 1;
#line 3343 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3350 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3350 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3359 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3359 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3367 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3367 "cgesvd.f"
			cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3369 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3377 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3381 "cgesvd.f"
		    }

#line 3383 "cgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3390 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3390 "cgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3394 "cgesvd.f"
			iu = 1;
#line 3395 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3399 "cgesvd.f"
			    ldwrku = *lda;
#line 3400 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3404 "cgesvd.f"
			    ldwrku = *m;
#line 3405 "cgesvd.f"
			}
#line 3406 "cgesvd.f"
			itau = iu + ldwrku * *m;
#line 3407 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3413 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3413 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3415 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3421 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3421 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3426 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3428 "cgesvd.f"
			i__2 = *m - 1;
#line 3428 "cgesvd.f"
			i__3 = *m - 1;
#line 3428 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3430 "cgesvd.f"
			ie = 1;
#line 3431 "cgesvd.f"
			itauq = itau;
#line 3432 "cgesvd.f"
			itaup = itauq + *m;
#line 3433 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3439 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3439 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3443 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3450 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3450 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3458 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3458 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3460 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3468 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3477 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3482 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3484 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3488 "cgesvd.f"
			itau = 1;
#line 3489 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3495 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3495 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3497 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3503 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3503 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3508 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3509 "cgesvd.f"
			i__2 = *m - 1;
#line 3509 "cgesvd.f"
			i__3 = *m - 1;
#line 3509 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3511 "cgesvd.f"
			ie = 1;
#line 3512 "cgesvd.f"
			itauq = itau;
#line 3513 "cgesvd.f"
			itaup = itauq + *m;
#line 3514 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3520 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3520 "cgesvd.f"
			cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3529 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3529 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3537 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3537 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3539 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3547 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3551 "cgesvd.f"
		    }

#line 3553 "cgesvd.f"
		}

#line 3555 "cgesvd.f"
	    }

#line 3557 "cgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3564 "cgesvd.f"
	    ie = 1;
#line 3565 "cgesvd.f"
	    itauq = 1;
#line 3566 "cgesvd.f"
	    itaup = itauq + *m;
#line 3567 "cgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 3573 "cgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3573 "cgesvd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 3576 "cgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3583 "cgesvd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3584 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3584 "cgesvd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3586 "cgesvd.f"
	    }
#line 3587 "cgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB) */
/*              (RWorkspace: 0) */

#line 3594 "cgesvd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3595 "cgesvd.f"
		if (wntva) {
#line 3595 "cgesvd.f"
		    nrvt = *n;
#line 3595 "cgesvd.f"
		}
#line 3597 "cgesvd.f"
		if (wntvs) {
#line 3597 "cgesvd.f"
		    nrvt = *m;
#line 3597 "cgesvd.f"
		}
#line 3599 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3599 "cgesvd.f"
		cungbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3601 "cgesvd.f"
	    }
#line 3602 "cgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3609 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3609 "cgesvd.f"
		cungbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3611 "cgesvd.f"
	    }
#line 3612 "cgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 3619 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3619 "cgesvd.f"
		cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3621 "cgesvd.f"
	    }
#line 3622 "cgesvd.f"
	    irwork = ie + *m;
#line 3623 "cgesvd.f"
	    if (wntuas || wntuo) {
#line 3623 "cgesvd.f"
		nru = *m;
#line 3623 "cgesvd.f"
	    }
#line 3625 "cgesvd.f"
	    if (wntun) {
#line 3625 "cgesvd.f"
		nru = 0;
#line 3625 "cgesvd.f"
	    }
#line 3627 "cgesvd.f"
	    if (wntvas || wntvo) {
#line 3627 "cgesvd.f"
		ncvt = *n;
#line 3627 "cgesvd.f"
	    }
#line 3629 "cgesvd.f"
	    if (wntvn) {
#line 3629 "cgesvd.f"
		ncvt = 0;
#line 3629 "cgesvd.f"
	    }
#line 3631 "cgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3639 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3642 "cgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3650 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3653 "cgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3661 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3664 "cgesvd.f"
	    }

#line 3666 "cgesvd.f"
	}

#line 3668 "cgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3672 "cgesvd.f"
    if (iscl == 1) {
#line 3673 "cgesvd.f"
	if (anrm > bignum) {
#line 3673 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3673 "cgesvd.f"
	}
#line 3676 "cgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3676 "cgesvd.f"
	    i__2 = minmn - 1;
#line 3676 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3676 "cgesvd.f"
	}
#line 3679 "cgesvd.f"
	if (anrm < smlnum) {
#line 3679 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3679 "cgesvd.f"
	}
#line 3682 "cgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3682 "cgesvd.f"
	    i__2 = minmn - 1;
#line 3682 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3682 "cgesvd.f"
	}
#line 3685 "cgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3689 "cgesvd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 3691 "cgesvd.f"
    return 0;

/*     End of CGESVD */

} /* cgesvd_ */

