#line 1 "zgesvd.f"
/* zgesvd.f -- translated by f2c (version 20100827).
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

#line 1 "zgesvd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__1 = 1;

/* > \brief <b> ZGESVD computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGESVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU, JOBVT */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGESVD computes the singular value decomposition (SVD) of a complex */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX*16 array, dimension (LDU,UCOL) */
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
/* >          VT is COMPLEX*16 array, dimension (LDVT,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (5*min(M,N)) */
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
/* >          > 0:  if ZBDSQR did not converge, INFO specifies how many */
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

/* > \ingroup complex16GEsing */

/*  ===================================================================== */
/* Subroutine */ int zgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
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
    static integer ierr, itau, ncvt, nrvt, lwork_zgebrd__, lwork_zgelqf__, 
	    lwork_zgeqrf__;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer wrkbl, itaup, itauq, mnthr, iwork;
    static logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen),
	     zgebrd_();
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zgelqf_(), zlascl_(char *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen), zgeqrf_(), 
	    zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen), zlaset_(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen);
    static integer ldwrkr;
    extern /* Subroutine */ int zbdsqr_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen);
    static integer minwrk, ldwrku, maxwrk;
    extern /* Subroutine */ int zungbr_();
    static doublereal smlnum;
    static integer irwork;
    extern /* Subroutine */ int zunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen), zunglq_();
    static logical lquery, wntuas, wntvas;
    extern /* Subroutine */ int zungqr_();
    static integer lwork_zungbr_p__, lwork_zungbr_q__, lwork_zunglq_m__, 
	    lwork_zunglq_n__, lwork_zungqr_m__, lwork_zungqr_n__;


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

#line 275 "zgesvd.f"
    /* Parameter adjustments */
#line 275 "zgesvd.f"
    a_dim1 = *lda;
#line 275 "zgesvd.f"
    a_offset = 1 + a_dim1;
#line 275 "zgesvd.f"
    a -= a_offset;
#line 275 "zgesvd.f"
    --s;
#line 275 "zgesvd.f"
    u_dim1 = *ldu;
#line 275 "zgesvd.f"
    u_offset = 1 + u_dim1;
#line 275 "zgesvd.f"
    u -= u_offset;
#line 275 "zgesvd.f"
    vt_dim1 = *ldvt;
#line 275 "zgesvd.f"
    vt_offset = 1 + vt_dim1;
#line 275 "zgesvd.f"
    vt -= vt_offset;
#line 275 "zgesvd.f"
    --work;
#line 275 "zgesvd.f"
    --rwork;
#line 275 "zgesvd.f"

#line 275 "zgesvd.f"
    /* Function Body */
#line 275 "zgesvd.f"
    *info = 0;
#line 276 "zgesvd.f"
    minmn = min(*m,*n);
#line 277 "zgesvd.f"
    wntua = lsame_(jobu, "A", (ftnlen)1, (ftnlen)1);
#line 278 "zgesvd.f"
    wntus = lsame_(jobu, "S", (ftnlen)1, (ftnlen)1);
#line 279 "zgesvd.f"
    wntuas = wntua || wntus;
#line 280 "zgesvd.f"
    wntuo = lsame_(jobu, "O", (ftnlen)1, (ftnlen)1);
#line 281 "zgesvd.f"
    wntun = lsame_(jobu, "N", (ftnlen)1, (ftnlen)1);
#line 282 "zgesvd.f"
    wntva = lsame_(jobvt, "A", (ftnlen)1, (ftnlen)1);
#line 283 "zgesvd.f"
    wntvs = lsame_(jobvt, "S", (ftnlen)1, (ftnlen)1);
#line 284 "zgesvd.f"
    wntvas = wntva || wntvs;
#line 285 "zgesvd.f"
    wntvo = lsame_(jobvt, "O", (ftnlen)1, (ftnlen)1);
#line 286 "zgesvd.f"
    wntvn = lsame_(jobvt, "N", (ftnlen)1, (ftnlen)1);
#line 287 "zgesvd.f"
    lquery = *lwork == -1;

#line 289 "zgesvd.f"
    if (! (wntua || wntus || wntuo || wntun)) {
#line 290 "zgesvd.f"
	*info = -1;
#line 291 "zgesvd.f"
    } else if (! (wntva || wntvs || wntvo || wntvn) || wntvo && wntuo) {
#line 293 "zgesvd.f"
	*info = -2;
#line 294 "zgesvd.f"
    } else if (*m < 0) {
#line 295 "zgesvd.f"
	*info = -3;
#line 296 "zgesvd.f"
    } else if (*n < 0) {
#line 297 "zgesvd.f"
	*info = -4;
#line 298 "zgesvd.f"
    } else if (*lda < max(1,*m)) {
#line 299 "zgesvd.f"
	*info = -6;
#line 300 "zgesvd.f"
    } else if (*ldu < 1 || wntuas && *ldu < *m) {
#line 301 "zgesvd.f"
	*info = -9;
#line 302 "zgesvd.f"
    } else if (*ldvt < 1 || wntva && *ldvt < *n || wntvs && *ldvt < minmn) {
#line 304 "zgesvd.f"
	*info = -11;
#line 305 "zgesvd.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to */
/*       real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 315 "zgesvd.f"
    if (*info == 0) {
#line 316 "zgesvd.f"
	minwrk = 1;
#line 317 "zgesvd.f"
	maxwrk = 1;
#line 318 "zgesvd.f"
	if (*m >= *n && minmn > 0) {

/*           Space needed for ZBDSQR is BDSPAC = 5*N */

/* Writing concatenation */
#line 322 "zgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 322 "zgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 322 "zgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 322 "zgesvd.f"
	    mnthr = ilaenv_(&c__6, "ZGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
/*           Compute space needed for ZGEQRF */
#line 324 "zgesvd.f"
	    zgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 325 "zgesvd.f"
	    lwork_zgeqrf__ = (integer) dum[0];
/*           Compute space needed for ZUNGQR */
#line 327 "zgesvd.f"
	    zungqr_(m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 328 "zgesvd.f"
	    lwork_zungqr_n__ = (integer) dum[0];
#line 329 "zgesvd.f"
	    zungqr_(m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 330 "zgesvd.f"
	    lwork_zungqr_m__ = (integer) dum[0];
/*           Compute space needed for ZGEBRD */
#line 332 "zgesvd.f"
	    zgebrd_(n, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 334 "zgesvd.f"
	    lwork_zgebrd__ = (integer) dum[0];
/*           Compute space needed for ZUNGBR */
#line 336 "zgesvd.f"
	    zungbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 338 "zgesvd.f"
	    lwork_zungbr_p__ = (integer) dum[0];
#line 339 "zgesvd.f"
	    zungbr_("Q", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 341 "zgesvd.f"
	    lwork_zungbr_q__ = (integer) dum[0];

#line 343 "zgesvd.f"
	    if (*m >= mnthr) {
#line 344 "zgesvd.f"
		if (wntun) {

/*                 Path 1 (M much larger than N, JOBU='N') */

#line 348 "zgesvd.f"
		    maxwrk = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 349 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_zgebrd__;
#line 349 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 350 "zgesvd.f"
		    if (wntvo || wntvas) {
/* Computing MAX */
#line 350 "zgesvd.f"
			i__2 = maxwrk, i__3 = (*n << 1) + lwork_zungbr_p__;
#line 350 "zgesvd.f"
			maxwrk = max(i__2,i__3);
#line 350 "zgesvd.f"
		    }
#line 352 "zgesvd.f"
		    minwrk = *n * 3;
#line 353 "zgesvd.f"
		} else if (wntuo && wntvn) {

/*                 Path 2 (M much larger than N, JOBU='O', JOBVT='N') */

#line 357 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 358 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_n__;
#line 358 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 359 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 359 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 360 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 360 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 361 "zgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n;
#line 361 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 362 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 363 "zgesvd.f"
		} else if (wntuo && wntvas) {

/*                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or */
/*                 'A') */

#line 368 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 369 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_n__;
#line 369 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 370 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 370 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 371 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 371 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 372 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_p__;
#line 372 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 373 "zgesvd.f"
		    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n;
#line 373 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 374 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 375 "zgesvd.f"
		} else if (wntus && wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */

#line 379 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 380 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_n__;
#line 380 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 381 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 381 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 382 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 382 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 383 "zgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 384 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 385 "zgesvd.f"
		} else if (wntus && wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */

#line 389 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 390 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_n__;
#line 390 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 391 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 391 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 392 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 392 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 393 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_p__;
#line 393 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 394 "zgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
#line 395 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 396 "zgesvd.f"
		} else if (wntus && wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or */
/*                 'A') */

#line 401 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 402 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_n__;
#line 402 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 403 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 403 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 404 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 404 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 405 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_p__;
#line 405 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 406 "zgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 407 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 408 "zgesvd.f"
		} else if (wntua && wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */

#line 412 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 413 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_m__;
#line 413 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 414 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 414 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 415 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 415 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 416 "zgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 417 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 418 "zgesvd.f"
		} else if (wntua && wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */

#line 422 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 423 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_m__;
#line 423 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 424 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 424 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 425 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 425 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 426 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_p__;
#line 426 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 427 "zgesvd.f"
		    maxwrk = (*n << 1) * *n + wrkbl;
#line 428 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 429 "zgesvd.f"
		} else if (wntua && wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or */
/*                 'A') */

#line 434 "zgesvd.f"
		    wrkbl = *n + lwork_zgeqrf__;
/* Computing MAX */
#line 435 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *n + lwork_zungqr_m__;
#line 435 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 436 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zgebrd__;
#line 436 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 437 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 437 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 438 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*n << 1) + lwork_zungbr_p__;
#line 438 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 439 "zgesvd.f"
		    maxwrk = *n * *n + wrkbl;
#line 440 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 441 "zgesvd.f"
		}
#line 442 "zgesvd.f"
	    } else {

/*              Path 10 (M at least N, but not much larger) */

#line 446 "zgesvd.f"
		zgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 448 "zgesvd.f"
		lwork_zgebrd__ = (integer) dum[0];
#line 449 "zgesvd.f"
		maxwrk = (*n << 1) + lwork_zgebrd__;
#line 450 "zgesvd.f"
		if (wntus || wntuo) {
#line 451 "zgesvd.f"
		    zungbr_("Q", m, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 453 "zgesvd.f"
		    lwork_zungbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 454 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 454 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 455 "zgesvd.f"
		}
#line 456 "zgesvd.f"
		if (wntua) {
#line 457 "zgesvd.f"
		    zungbr_("Q", m, m, n, &a[a_offset], lda, dum, dum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 459 "zgesvd.f"
		    lwork_zungbr_q__ = (integer) dum[0];
/* Computing MAX */
#line 460 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_zungbr_q__;
#line 460 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 461 "zgesvd.f"
		}
#line 462 "zgesvd.f"
		if (! wntvn) {
/* Computing MAX */
#line 463 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*n << 1) + lwork_zungbr_p__;
#line 463 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 464 "zgesvd.f"
		    minwrk = (*n << 1) + *m;
#line 465 "zgesvd.f"
		}
#line 466 "zgesvd.f"
	    }
#line 467 "zgesvd.f"
	} else if (minmn > 0) {

/*           Space needed for ZBDSQR is BDSPAC = 5*M */

/* Writing concatenation */
#line 471 "zgesvd.f"
	    i__1[0] = 1, a__1[0] = jobu;
#line 471 "zgesvd.f"
	    i__1[1] = 1, a__1[1] = jobvt;
#line 471 "zgesvd.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 471 "zgesvd.f"
	    mnthr = ilaenv_(&c__6, "ZGESVD", ch__1, m, n, &c__0, &c__0, (
		    ftnlen)6, (ftnlen)2);
/*           Compute space needed for ZGELQF */
#line 473 "zgesvd.f"
	    zgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 474 "zgesvd.f"
	    lwork_zgelqf__ = (integer) dum[0];
/*           Compute space needed for ZUNGLQ */
#line 476 "zgesvd.f"
	    zunglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
#line 477 "zgesvd.f"
	    lwork_zunglq_n__ = (integer) dum[0];
#line 478 "zgesvd.f"
	    zunglq_(m, n, m, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
#line 479 "zgesvd.f"
	    lwork_zunglq_m__ = (integer) dum[0];
/*           Compute space needed for ZGEBRD */
#line 481 "zgesvd.f"
	    zgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
		     &ierr);
#line 483 "zgesvd.f"
	    lwork_zgebrd__ = (integer) dum[0];
/*            Compute space needed for ZUNGBR P */
#line 485 "zgesvd.f"
	    zungbr_("P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 487 "zgesvd.f"
	    lwork_zungbr_p__ = (integer) dum[0];
/*           Compute space needed for ZUNGBR Q */
#line 489 "zgesvd.f"
	    zungbr_("Q", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
		    ftnlen)1);
#line 491 "zgesvd.f"
	    lwork_zungbr_q__ = (integer) dum[0];
#line 492 "zgesvd.f"
	    if (*n >= mnthr) {
#line 493 "zgesvd.f"
		if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N') */

#line 497 "zgesvd.f"
		    maxwrk = *m + lwork_zgelqf__;
/* Computing MAX */
#line 498 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zgebrd__;
#line 498 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 499 "zgesvd.f"
		    if (wntuo || wntuas) {
/* Computing MAX */
#line 499 "zgesvd.f"
			i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 499 "zgesvd.f"
			maxwrk = max(i__2,i__3);
#line 499 "zgesvd.f"
		    }
#line 501 "zgesvd.f"
		    minwrk = *m * 3;
#line 502 "zgesvd.f"
		} else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O') */

#line 506 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 507 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 507 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 508 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 508 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 509 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 509 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 510 "zgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 510 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 511 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 512 "zgesvd.f"
		} else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='O') */

#line 517 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 518 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 518 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 519 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 519 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 520 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 520 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 521 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 521 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 522 "zgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 522 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 523 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 524 "zgesvd.f"
		} else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */

#line 528 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 529 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 529 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 530 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 530 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 531 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 531 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 532 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 533 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 534 "zgesvd.f"
		} else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */

#line 538 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 539 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 539 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 540 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 540 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 541 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 541 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 542 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 542 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 543 "zgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 544 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 545 "zgesvd.f"
		} else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='S') */

#line 550 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 551 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 551 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 552 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 552 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 553 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 553 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 554 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 554 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 555 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 556 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 557 "zgesvd.f"
		} else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */

#line 561 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 562 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_n__;
#line 562 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 563 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 563 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 564 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 564 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 565 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 566 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 567 "zgesvd.f"
		} else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */

#line 571 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 572 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_n__;
#line 572 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 573 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 573 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 574 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 574 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 575 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 575 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 576 "zgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 577 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 578 "zgesvd.f"
		} else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='A') */

#line 583 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 584 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_n__;
#line 584 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 585 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 585 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 586 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 586 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 587 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 587 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 588 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 589 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 590 "zgesvd.f"
		}
#line 591 "zgesvd.f"
	    } else {

/*              Path 10t(N greater than M, but not much larger) */

#line 595 "zgesvd.f"
		zgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
			c_n1, &ierr);
#line 597 "zgesvd.f"
		lwork_zgebrd__ = (integer) dum[0];
#line 598 "zgesvd.f"
		maxwrk = (*m << 1) + lwork_zgebrd__;
#line 599 "zgesvd.f"
		if (wntvs || wntvo) {
/*                Compute space needed for ZUNGBR P */
#line 601 "zgesvd.f"
		    zungbr_("P", m, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 603 "zgesvd.f"
		    lwork_zungbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 604 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 604 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 605 "zgesvd.f"
		}
#line 606 "zgesvd.f"
		if (wntva) {
#line 607 "zgesvd.f"
		    zungbr_("P", n, n, m, &a[a_offset], n, dum, dum, &c_n1, &
			    ierr, (ftnlen)1);
#line 609 "zgesvd.f"
		    lwork_zungbr_p__ = (integer) dum[0];
/* Computing MAX */
#line 610 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 610 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 611 "zgesvd.f"
		}
#line 612 "zgesvd.f"
		if (! wntun) {
/* Computing MAX */
#line 613 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 613 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 614 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 615 "zgesvd.f"
		}
#line 616 "zgesvd.f"
	    }
#line 617 "zgesvd.f"
	}
#line 618 "zgesvd.f"
	maxwrk = max(maxwrk,minwrk);
#line 619 "zgesvd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 621 "zgesvd.f"
	if (*lwork < minwrk && ! lquery) {
#line 622 "zgesvd.f"
	    *info = -13;
#line 623 "zgesvd.f"
	}
#line 624 "zgesvd.f"
    }

#line 626 "zgesvd.f"
    if (*info != 0) {
#line 627 "zgesvd.f"
	i__2 = -(*info);
#line 627 "zgesvd.f"
	xerbla_("ZGESVD", &i__2, (ftnlen)6);
#line 628 "zgesvd.f"
	return 0;
#line 629 "zgesvd.f"
    } else if (lquery) {
#line 630 "zgesvd.f"
	return 0;
#line 631 "zgesvd.f"
    }

/*     Quick return if possible */

#line 635 "zgesvd.f"
    if (*m == 0 || *n == 0) {
#line 636 "zgesvd.f"
	return 0;
#line 637 "zgesvd.f"
    }

/*     Get machine constants */

#line 641 "zgesvd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 642 "zgesvd.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 643 "zgesvd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 647 "zgesvd.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 648 "zgesvd.f"
    iscl = 0;
#line 649 "zgesvd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 650 "zgesvd.f"
	iscl = 1;
#line 651 "zgesvd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 652 "zgesvd.f"
    } else if (anrm > bignum) {
#line 653 "zgesvd.f"
	iscl = 1;
#line 654 "zgesvd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 655 "zgesvd.f"
    }

#line 657 "zgesvd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 663 "zgesvd.f"
	if (*m >= mnthr) {

#line 665 "zgesvd.f"
	    if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N') */
/*              No left singular vectors to be computed */

#line 670 "zgesvd.f"
		itau = 1;
#line 671 "zgesvd.f"
		iwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: need 0) */

#line 677 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 677 "zgesvd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out below R */

#line 682 "zgesvd.f"
		i__2 = *n - 1;
#line 682 "zgesvd.f"
		i__3 = *n - 1;
#line 682 "zgesvd.f"
		zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 684 "zgesvd.f"
		ie = 1;
#line 685 "zgesvd.f"
		itauq = 1;
#line 686 "zgesvd.f"
		itaup = itauq + *n;
#line 687 "zgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 693 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 693 "zgesvd.f"
		zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 696 "zgesvd.f"
		ncvt = 0;
#line 697 "zgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 703 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 703 "zgesvd.f"
		    zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 705 "zgesvd.f"
		    ncvt = *n;
#line 706 "zgesvd.f"
		}
#line 707 "zgesvd.f"
		irwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 714 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			irwork], info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 719 "zgesvd.f"
		if (wntvas) {
#line 719 "zgesvd.f"
		    zlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 719 "zgesvd.f"
		}

#line 722 "zgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

#line 728 "zgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 732 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 733 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 733 "zgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 737 "zgesvd.f"
			ldwrku = *lda;
#line 738 "zgesvd.f"
			ldwrkr = *lda;
#line 739 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 739 "zgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 739 "zgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 743 "zgesvd.f"
			    ldwrku = *lda;
#line 744 "zgesvd.f"
			    ldwrkr = *n;
#line 745 "zgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 749 "zgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 750 "zgesvd.f"
			    ldwrkr = *n;
#line 751 "zgesvd.f"
			}
#line 751 "zgesvd.f"
		    }
#line 752 "zgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 753 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 759 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 759 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 764 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 765 "zgesvd.f"
		    i__2 = *n - 1;
#line 765 "zgesvd.f"
		    i__3 = *n - 1;
#line 765 "zgesvd.f"
		    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1], &
			    ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 772 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 772 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 774 "zgesvd.f"
		    ie = 1;
#line 775 "zgesvd.f"
		    itauq = itau;
#line 776 "zgesvd.f"
		    itaup = itauq + *n;
#line 777 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 783 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 783 "zgesvd.f"
		    zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: need 0) */

#line 791 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 791 "zgesvd.f"
		    zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 794 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 801 "zgesvd.f"
		    zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &work[ir], &ldwrkr, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 804 "zgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 811 "zgesvd.f"
		    i__2 = *m;
#line 811 "zgesvd.f"
		    i__3 = ldwrku;
#line 811 "zgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 812 "zgesvd.f"
			i__4 = *m - i__ + 1;
#line 812 "zgesvd.f"
			chunk = min(i__4,ldwrku);
#line 813 "zgesvd.f"
			zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 816 "zgesvd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 818 "zgesvd.f"
/* L10: */
#line 818 "zgesvd.f"
		    }

#line 820 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 824 "zgesvd.f"
		    ie = 1;
#line 825 "zgesvd.f"
		    itauq = 1;
#line 826 "zgesvd.f"
		    itaup = itauq + *n;
#line 827 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*                 (RWorkspace: N) */

#line 833 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 833 "zgesvd.f"
		    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 841 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 841 "zgesvd.f"
		    zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 843 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (CWorkspace: need 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 850 "zgesvd.f"
		    zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &a[a_offset], lda, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 853 "zgesvd.f"
		}

#line 855 "zgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 861 "zgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 865 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 866 "zgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 866 "zgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 870 "zgesvd.f"
			ldwrku = *lda;
#line 871 "zgesvd.f"
			ldwrkr = *lda;
#line 872 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 872 "zgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 872 "zgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 876 "zgesvd.f"
			    ldwrku = *lda;
#line 877 "zgesvd.f"
			    ldwrkr = *n;
#line 878 "zgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 882 "zgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 883 "zgesvd.f"
			    ldwrkr = *n;
#line 884 "zgesvd.f"
			}
#line 884 "zgesvd.f"
		    }
#line 885 "zgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 886 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 892 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 892 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 897 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 898 "zgesvd.f"
		    if (*n > 1) {
#line 898 "zgesvd.f"
			i__3 = *n - 1;
#line 898 "zgesvd.f"
			i__2 = *n - 1;
#line 898 "zgesvd.f"
			zlaset_("L", &i__3, &i__2, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 898 "zgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 906 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 906 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 908 "zgesvd.f"
		    ie = 1;
#line 909 "zgesvd.f"
		    itauq = itau;
#line 910 "zgesvd.f"
		    itaup = itauq + *n;
#line 911 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 917 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 917 "zgesvd.f"
		    zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 920 "zgesvd.f"
		    zlacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 926 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 926 "zgesvd.f"
		    zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 934 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 934 "zgesvd.f"
		    zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 936 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 944 "zgesvd.f"
		    zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, cdum, &c__1,
			     &rwork[irwork], info, (ftnlen)1);
#line 947 "zgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 954 "zgesvd.f"
		    i__3 = *m;
#line 954 "zgesvd.f"
		    i__2 = ldwrku;
#line 954 "zgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 955 "zgesvd.f"
			i__4 = *m - i__ + 1;
#line 955 "zgesvd.f"
			chunk = min(i__4,ldwrku);
#line 956 "zgesvd.f"
			zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 959 "zgesvd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 961 "zgesvd.f"
/* L20: */
#line 961 "zgesvd.f"
		    }

#line 963 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 967 "zgesvd.f"
		    itau = 1;
#line 968 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 974 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 974 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 979 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 980 "zgesvd.f"
		    if (*n > 1) {
#line 980 "zgesvd.f"
			i__2 = *n - 1;
#line 980 "zgesvd.f"
			i__3 = *n - 1;
#line 980 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 980 "zgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 988 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 988 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 990 "zgesvd.f"
		    ie = 1;
#line 991 "zgesvd.f"
		    itauq = itau;
#line 992 "zgesvd.f"
		    itaup = itauq + *n;
#line 993 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                 (RWorkspace: N) */

#line 999 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 999 "zgesvd.f"
		    zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                 (RWorkspace: 0) */

#line 1007 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1007 "zgesvd.f"
		    zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 1015 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1015 "zgesvd.f"
		    zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1017 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 1025 "zgesvd.f"
		    zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, 
			    &rwork[irwork], info, (ftnlen)1);

#line 1029 "zgesvd.f"
		}

#line 1031 "zgesvd.f"
	    } else if (wntus) {

#line 1033 "zgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

#line 1039 "zgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1043 "zgesvd.f"
			ir = 1;
#line 1044 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1048 "zgesvd.f"
			    ldwrkr = *lda;
#line 1049 "zgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1053 "zgesvd.f"
			    ldwrkr = *n;
#line 1054 "zgesvd.f"
			}
#line 1055 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1056 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1062 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1062 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1067 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1069 "zgesvd.f"
			i__2 = *n - 1;
#line 1069 "zgesvd.f"
			i__3 = *n - 1;
#line 1069 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1076 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1076 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1078 "zgesvd.f"
			ie = 1;
#line 1079 "zgesvd.f"
			itauq = itau;
#line 1080 "zgesvd.f"
			itaup = itauq + *n;
#line 1081 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1087 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1087 "zgesvd.f"
			zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1096 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1096 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1099 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1106 "zgesvd.f"
			zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1115 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1118 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1122 "zgesvd.f"
			itau = 1;
#line 1123 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1129 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1129 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1131 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1137 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1137 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1139 "zgesvd.f"
			ie = 1;
#line 1140 "zgesvd.f"
			itauq = itau;
#line 1141 "zgesvd.f"
			itaup = itauq + *n;
#line 1142 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1146 "zgesvd.f"
			i__2 = *n - 1;
#line 1146 "zgesvd.f"
			i__3 = *n - 1;
#line 1146 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1153 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1153 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1161 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1161 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1164 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1171 "zgesvd.f"
			zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1175 "zgesvd.f"
		    }

#line 1177 "zgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

#line 1183 "zgesvd.f"
		    if (*lwork >= (*n << 1) * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1187 "zgesvd.f"
			iu = 1;
#line 1188 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1192 "zgesvd.f"
			    ldwrku = *lda;
#line 1193 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1194 "zgesvd.f"
			    ldwrkr = *lda;
#line 1195 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1199 "zgesvd.f"
			    ldwrku = *lda;
#line 1200 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1201 "zgesvd.f"
			    ldwrkr = *n;
#line 1202 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1206 "zgesvd.f"
			    ldwrku = *n;
#line 1207 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1208 "zgesvd.f"
			    ldwrkr = *n;
#line 1209 "zgesvd.f"
			}
#line 1210 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1211 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1217 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1217 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1222 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1224 "zgesvd.f"
			i__2 = *n - 1;
#line 1224 "zgesvd.f"
			i__3 = *n - 1;
#line 1224 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1231 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1231 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1233 "zgesvd.f"
			ie = 1;
#line 1234 "zgesvd.f"
			itauq = itau;
#line 1235 "zgesvd.f"
			itaup = itauq + *n;
#line 1236 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1244 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1244 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1248 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1255 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1255 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1264 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1264 "zgesvd.f"
			zungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1267 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1275 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1285 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1292 "zgesvd.f"
			zlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1295 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1299 "zgesvd.f"
			itau = 1;
#line 1300 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1306 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1306 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1308 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1314 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1314 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1316 "zgesvd.f"
			ie = 1;
#line 1317 "zgesvd.f"
			itauq = itau;
#line 1318 "zgesvd.f"
			itaup = itauq + *n;
#line 1319 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1323 "zgesvd.f"
			i__2 = *n - 1;
#line 1323 "zgesvd.f"
			i__3 = *n - 1;
#line 1323 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1330 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1330 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1338 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1338 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1346 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1346 "zgesvd.f"
			zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1348 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1356 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1360 "zgesvd.f"
		    }

#line 1362 "zgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

#line 1369 "zgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1373 "zgesvd.f"
			iu = 1;
#line 1374 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1378 "zgesvd.f"
			    ldwrku = *lda;
#line 1379 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1383 "zgesvd.f"
			    ldwrku = *n;
#line 1384 "zgesvd.f"
			}
#line 1385 "zgesvd.f"
			itau = iu + ldwrku * *n;
#line 1386 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1392 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1392 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1397 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1399 "zgesvd.f"
			i__2 = *n - 1;
#line 1399 "zgesvd.f"
			i__3 = *n - 1;
#line 1399 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1406 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1406 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1408 "zgesvd.f"
			ie = 1;
#line 1409 "zgesvd.f"
			itauq = itau;
#line 1410 "zgesvd.f"
			itaup = itauq + *n;
#line 1411 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1417 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1417 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1421 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1428 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1428 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1437 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1437 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1439 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1447 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1456 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1459 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1463 "zgesvd.f"
			itau = 1;
#line 1464 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1470 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1470 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1472 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1478 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1478 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1483 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1484 "zgesvd.f"
			if (*n > 1) {
#line 1484 "zgesvd.f"
			    i__2 = *n - 1;
#line 1484 "zgesvd.f"
			    i__3 = *n - 1;
#line 1484 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1484 "zgesvd.f"
			}
#line 1487 "zgesvd.f"
			ie = 1;
#line 1488 "zgesvd.f"
			itauq = itau;
#line 1489 "zgesvd.f"
			itaup = itauq + *n;
#line 1490 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1496 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1496 "zgesvd.f"
			zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1505 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1505 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1513 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1513 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1515 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1523 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1527 "zgesvd.f"
		    }

#line 1529 "zgesvd.f"
		}

#line 1531 "zgesvd.f"
	    } else if (wntua) {

#line 1533 "zgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1539 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1539 "zgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1543 "zgesvd.f"
			ir = 1;
#line 1544 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1548 "zgesvd.f"
			    ldwrkr = *lda;
#line 1549 "zgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1553 "zgesvd.f"
			    ldwrkr = *n;
#line 1554 "zgesvd.f"
			}
#line 1555 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1556 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1562 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1562 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1564 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1568 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1570 "zgesvd.f"
			i__2 = *n - 1;
#line 1570 "zgesvd.f"
			i__3 = *n - 1;
#line 1570 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1577 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1577 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1579 "zgesvd.f"
			ie = 1;
#line 1580 "zgesvd.f"
			itauq = itau;
#line 1581 "zgesvd.f"
			itaup = itauq + *n;
#line 1582 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1588 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1588 "zgesvd.f"
			zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1597 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1597 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1600 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1607 "zgesvd.f"
			zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1616 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1621 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1623 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1627 "zgesvd.f"
			itau = 1;
#line 1628 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1634 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1634 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1636 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1642 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1642 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1644 "zgesvd.f"
			ie = 1;
#line 1645 "zgesvd.f"
			itauq = itau;
#line 1646 "zgesvd.f"
			itaup = itauq + *n;
#line 1647 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1651 "zgesvd.f"
			i__2 = *n - 1;
#line 1651 "zgesvd.f"
			i__3 = *n - 1;
#line 1651 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1658 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1658 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1667 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1667 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1670 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1677 "zgesvd.f"
			zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1681 "zgesvd.f"
		    }

#line 1683 "zgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1689 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1689 "zgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1693 "zgesvd.f"
			iu = 1;
#line 1694 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1698 "zgesvd.f"
			    ldwrku = *lda;
#line 1699 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1700 "zgesvd.f"
			    ldwrkr = *lda;
#line 1701 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1705 "zgesvd.f"
			    ldwrku = *lda;
#line 1706 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1707 "zgesvd.f"
			    ldwrkr = *n;
#line 1708 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1712 "zgesvd.f"
			    ldwrku = *n;
#line 1713 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1714 "zgesvd.f"
			    ldwrkr = *n;
#line 1715 "zgesvd.f"
			}
#line 1716 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1717 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1723 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1723 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1725 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1731 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1731 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1736 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1738 "zgesvd.f"
			i__2 = *n - 1;
#line 1738 "zgesvd.f"
			i__3 = *n - 1;
#line 1738 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1740 "zgesvd.f"
			ie = 1;
#line 1741 "zgesvd.f"
			itauq = itau;
#line 1742 "zgesvd.f"
			itaup = itauq + *n;
#line 1743 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1751 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1751 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1755 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1762 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1762 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1771 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1771 "zgesvd.f"
			zungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1774 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1782 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1792 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1797 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1801 "zgesvd.f"
			zlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1804 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1808 "zgesvd.f"
			itau = 1;
#line 1809 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1815 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1815 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1817 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1823 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1823 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1825 "zgesvd.f"
			ie = 1;
#line 1826 "zgesvd.f"
			itauq = itau;
#line 1827 "zgesvd.f"
			itaup = itauq + *n;
#line 1828 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1832 "zgesvd.f"
			i__2 = *n - 1;
#line 1832 "zgesvd.f"
			i__3 = *n - 1;
#line 1832 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1839 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1839 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1848 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1848 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1856 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1856 "zgesvd.f"
			zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1858 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1866 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1870 "zgesvd.f"
		    }

#line 1872 "zgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1879 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1879 "zgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1883 "zgesvd.f"
			iu = 1;
#line 1884 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1888 "zgesvd.f"
			    ldwrku = *lda;
#line 1889 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1893 "zgesvd.f"
			    ldwrku = *n;
#line 1894 "zgesvd.f"
			}
#line 1895 "zgesvd.f"
			itau = iu + ldwrku * *n;
#line 1896 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1902 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1902 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1904 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1910 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1910 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1915 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1917 "zgesvd.f"
			i__2 = *n - 1;
#line 1917 "zgesvd.f"
			i__3 = *n - 1;
#line 1917 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1919 "zgesvd.f"
			ie = 1;
#line 1920 "zgesvd.f"
			itauq = itau;
#line 1921 "zgesvd.f"
			itaup = itauq + *n;
#line 1922 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1928 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1928 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1932 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1939 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1939 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: need   0) */

#line 1948 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1948 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1950 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1958 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1967 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1972 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1974 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1978 "zgesvd.f"
			itau = 1;
#line 1979 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1985 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1985 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1987 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1993 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1993 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 1998 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1999 "zgesvd.f"
			if (*n > 1) {
#line 1999 "zgesvd.f"
			    i__2 = *n - 1;
#line 1999 "zgesvd.f"
			    i__3 = *n - 1;
#line 1999 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1999 "zgesvd.f"
			}
#line 2002 "zgesvd.f"
			ie = 1;
#line 2003 "zgesvd.f"
			itauq = itau;
#line 2004 "zgesvd.f"
			itaup = itauq + *n;
#line 2005 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 2011 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2011 "zgesvd.f"
			zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2020 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2020 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2028 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2028 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 2030 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2038 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2042 "zgesvd.f"
		    }

#line 2044 "zgesvd.f"
		}

#line 2046 "zgesvd.f"
	    }

#line 2048 "zgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 2055 "zgesvd.f"
	    ie = 1;
#line 2056 "zgesvd.f"
	    itauq = 1;
#line 2057 "zgesvd.f"
	    itaup = itauq + *n;
#line 2058 "zgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 2064 "zgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 2064 "zgesvd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 2067 "zgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB) */
/*              (RWorkspace: 0) */

#line 2074 "zgesvd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 2075 "zgesvd.f"
		if (wntus) {
#line 2075 "zgesvd.f"
		    ncu = *n;
#line 2075 "zgesvd.f"
		}
#line 2077 "zgesvd.f"
		if (wntua) {
#line 2077 "zgesvd.f"
		    ncu = *m;
#line 2077 "zgesvd.f"
		}
#line 2079 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2079 "zgesvd.f"
		zungbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2081 "zgesvd.f"
	    }
#line 2082 "zgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2089 "zgesvd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2090 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2090 "zgesvd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2092 "zgesvd.f"
	    }
#line 2093 "zgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 2100 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2100 "zgesvd.f"
		zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2102 "zgesvd.f"
	    }
#line 2103 "zgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2110 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2110 "zgesvd.f"
		zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2112 "zgesvd.f"
	    }
#line 2113 "zgesvd.f"
	    irwork = ie + *n;
#line 2114 "zgesvd.f"
	    if (wntuas || wntuo) {
#line 2114 "zgesvd.f"
		nru = *m;
#line 2114 "zgesvd.f"
	    }
#line 2116 "zgesvd.f"
	    if (wntun) {
#line 2116 "zgesvd.f"
		nru = 0;
#line 2116 "zgesvd.f"
	    }
#line 2118 "zgesvd.f"
	    if (wntvas || wntvo) {
#line 2118 "zgesvd.f"
		ncvt = *n;
#line 2118 "zgesvd.f"
	    }
#line 2120 "zgesvd.f"
	    if (wntvn) {
#line 2120 "zgesvd.f"
		ncvt = 0;
#line 2120 "zgesvd.f"
	    }
#line 2122 "zgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2130 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2133 "zgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2141 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2144 "zgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2152 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2155 "zgesvd.f"
	    }

#line 2157 "zgesvd.f"
	}

#line 2159 "zgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2165 "zgesvd.f"
	if (*n >= mnthr) {

#line 2167 "zgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2172 "zgesvd.f"
		itau = 1;
#line 2173 "zgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 2179 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2179 "zgesvd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2184 "zgesvd.f"
		i__2 = *m - 1;
#line 2184 "zgesvd.f"
		i__3 = *m - 1;
#line 2184 "zgesvd.f"
		zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 2186 "zgesvd.f"
		ie = 1;
#line 2187 "zgesvd.f"
		itauq = 1;
#line 2188 "zgesvd.f"
		itaup = itauq + *m;
#line 2189 "zgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 2195 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2195 "zgesvd.f"
		zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2198 "zgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2204 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2204 "zgesvd.f"
		    zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2206 "zgesvd.f"
		}
#line 2207 "zgesvd.f"
		irwork = ie + *m;
#line 2208 "zgesvd.f"
		nru = 0;
#line 2209 "zgesvd.f"
		if (wntuo || wntuas) {
#line 2209 "zgesvd.f"
		    nru = *m;
#line 2209 "zgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2217 "zgesvd.f"
		zbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &rwork[ie], cdum, &
			c__1, &a[a_offset], lda, cdum, &c__1, &rwork[irwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2222 "zgesvd.f"
		if (wntuas) {
#line 2222 "zgesvd.f"
		    zlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2222 "zgesvd.f"
		}

#line 2225 "zgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

#line 2231 "zgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2235 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2236 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 2236 "zgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2240 "zgesvd.f"
			ldwrku = *lda;
#line 2241 "zgesvd.f"
			chunk = *n;
#line 2242 "zgesvd.f"
			ldwrkr = *lda;
#line 2243 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2243 "zgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 2243 "zgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2247 "zgesvd.f"
			    ldwrku = *lda;
#line 2248 "zgesvd.f"
			    chunk = *n;
#line 2249 "zgesvd.f"
			    ldwrkr = *m;
#line 2250 "zgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2254 "zgesvd.f"
			    ldwrku = *m;
#line 2255 "zgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2256 "zgesvd.f"
			    ldwrkr = *m;
#line 2257 "zgesvd.f"
			}
#line 2257 "zgesvd.f"
		    }
#line 2258 "zgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2259 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2265 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2265 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2270 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2271 "zgesvd.f"
		    i__2 = *m - 1;
#line 2271 "zgesvd.f"
		    i__3 = *m - 1;
#line 2271 "zgesvd.f"
		    zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2278 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2278 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2280 "zgesvd.f"
		    ie = 1;
#line 2281 "zgesvd.f"
		    itauq = itau;
#line 2282 "zgesvd.f"
		    itaup = itauq + *m;
#line 2283 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2289 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2289 "zgesvd.f"
		    zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2297 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2297 "zgesvd.f"
		    zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2300 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2307 "zgesvd.f"
		    zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &work[
			    ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2310 "zgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N) */
/*                 (RWorkspace: 0) */

#line 2317 "zgesvd.f"
		    i__2 = *n;
#line 2317 "zgesvd.f"
		    i__3 = chunk;
#line 2317 "zgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2318 "zgesvd.f"
			i__4 = *n - i__ + 1;
#line 2318 "zgesvd.f"
			blk = min(i__4,chunk);
#line 2319 "zgesvd.f"
			zgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2322 "zgesvd.f"
			zlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2324 "zgesvd.f"
/* L30: */
#line 2324 "zgesvd.f"
		    }

#line 2326 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2330 "zgesvd.f"
		    ie = 1;
#line 2331 "zgesvd.f"
		    itauq = 1;
#line 2332 "zgesvd.f"
		    itaup = itauq + *m;
#line 2333 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*                 (RWorkspace: need M) */

#line 2339 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2339 "zgesvd.f"
		    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2347 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2347 "zgesvd.f"
		    zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2349 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2356 "zgesvd.f"
		    zbdsqr_("L", m, n, &c__0, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 2359 "zgesvd.f"
		}

#line 2361 "zgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 2367 "zgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2371 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2372 "zgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 2372 "zgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2376 "zgesvd.f"
			ldwrku = *lda;
#line 2377 "zgesvd.f"
			chunk = *n;
#line 2378 "zgesvd.f"
			ldwrkr = *lda;
#line 2379 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2379 "zgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 2379 "zgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2383 "zgesvd.f"
			    ldwrku = *lda;
#line 2384 "zgesvd.f"
			    chunk = *n;
#line 2385 "zgesvd.f"
			    ldwrkr = *m;
#line 2386 "zgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2390 "zgesvd.f"
			    ldwrku = *m;
#line 2391 "zgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2392 "zgesvd.f"
			    ldwrkr = *m;
#line 2393 "zgesvd.f"
			}
#line 2393 "zgesvd.f"
		    }
#line 2394 "zgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2395 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2401 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2401 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2406 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2407 "zgesvd.f"
		    i__3 = *m - 1;
#line 2407 "zgesvd.f"
		    i__2 = *m - 1;
#line 2407 "zgesvd.f"
		    zlaset_("U", &i__3, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2414 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2414 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2416 "zgesvd.f"
		    ie = 1;
#line 2417 "zgesvd.f"
		    itauq = itau;
#line 2418 "zgesvd.f"
		    itaup = itauq + *m;
#line 2419 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2425 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2425 "zgesvd.f"
		    zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2428 "zgesvd.f"
		    zlacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2434 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2434 "zgesvd.f"
		    zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2442 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2442 "zgesvd.f"
		    zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2444 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2452 "zgesvd.f"
		    zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[ir],
			     &ldwrkr, &u[u_offset], ldu, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2455 "zgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N)) */
/*                 (RWorkspace: 0) */

#line 2462 "zgesvd.f"
		    i__3 = *n;
#line 2462 "zgesvd.f"
		    i__2 = chunk;
#line 2462 "zgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2463 "zgesvd.f"
			i__4 = *n - i__ + 1;
#line 2463 "zgesvd.f"
			blk = min(i__4,chunk);
#line 2464 "zgesvd.f"
			zgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2467 "zgesvd.f"
			zlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2469 "zgesvd.f"
/* L40: */
#line 2469 "zgesvd.f"
		    }

#line 2471 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2475 "zgesvd.f"
		    itau = 1;
#line 2476 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2482 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2482 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2487 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2488 "zgesvd.f"
		    i__2 = *m - 1;
#line 2488 "zgesvd.f"
		    i__3 = *m - 1;
#line 2488 "zgesvd.f"
		    zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2495 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2495 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2497 "zgesvd.f"
		    ie = 1;
#line 2498 "zgesvd.f"
		    itauq = itau;
#line 2499 "zgesvd.f"
		    itaup = itauq + *m;
#line 2500 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2506 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2506 "zgesvd.f"
		    zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                 (RWorkspace: 0) */

#line 2514 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2514 "zgesvd.f"
		    zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2522 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2522 "zgesvd.f"
		    zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2524 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2532 "zgesvd.f"
		    zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			    rwork[irwork], info, (ftnlen)1);

#line 2535 "zgesvd.f"
		}

#line 2537 "zgesvd.f"
	    } else if (wntvs) {

#line 2539 "zgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

#line 2545 "zgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2549 "zgesvd.f"
			ir = 1;
#line 2550 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2554 "zgesvd.f"
			    ldwrkr = *lda;
#line 2555 "zgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2559 "zgesvd.f"
			    ldwrkr = *m;
#line 2560 "zgesvd.f"
			}
#line 2561 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2562 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2568 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2568 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2573 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2575 "zgesvd.f"
			i__2 = *m - 1;
#line 2575 "zgesvd.f"
			i__3 = *m - 1;
#line 2575 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2582 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2582 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2584 "zgesvd.f"
			ie = 1;
#line 2585 "zgesvd.f"
			itauq = itau;
#line 2586 "zgesvd.f"
			itaup = itauq + *m;
#line 2587 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2593 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2593 "zgesvd.f"
			zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2603 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2603 "zgesvd.f"
			zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2606 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2613 "zgesvd.f"
			zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2622 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2625 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2629 "zgesvd.f"
			itau = 1;
#line 2630 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2636 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2636 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2641 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2647 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2647 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2649 "zgesvd.f"
			ie = 1;
#line 2650 "zgesvd.f"
			itauq = itau;
#line 2651 "zgesvd.f"
			itaup = itauq + *m;
#line 2652 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2656 "zgesvd.f"
			i__2 = *m - 1;
#line 2656 "zgesvd.f"
			i__3 = *m - 1;
#line 2656 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2663 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2663 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2671 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2671 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2674 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2681 "zgesvd.f"
			zbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 2685 "zgesvd.f"
		    }

#line 2687 "zgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

#line 2693 "zgesvd.f"
		    if (*lwork >= (*m << 1) * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2697 "zgesvd.f"
			iu = 1;
#line 2698 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2702 "zgesvd.f"
			    ldwrku = *lda;
#line 2703 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2704 "zgesvd.f"
			    ldwrkr = *lda;
#line 2705 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2709 "zgesvd.f"
			    ldwrku = *lda;
#line 2710 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2711 "zgesvd.f"
			    ldwrkr = *m;
#line 2712 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2716 "zgesvd.f"
			    ldwrku = *m;
#line 2717 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2718 "zgesvd.f"
			    ldwrkr = *m;
#line 2719 "zgesvd.f"
			}
#line 2720 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2721 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2727 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2727 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2732 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2734 "zgesvd.f"
			i__2 = *m - 1;
#line 2734 "zgesvd.f"
			i__3 = *m - 1;
#line 2734 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2741 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2741 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2743 "zgesvd.f"
			ie = 1;
#line 2744 "zgesvd.f"
			itauq = itau;
#line 2745 "zgesvd.f"
			itaup = itauq + *m;
#line 2746 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 2754 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2754 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2758 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2766 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2766 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2774 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2774 "zgesvd.f"
			zungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2777 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2785 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2795 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2802 "zgesvd.f"
			zlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2805 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2809 "zgesvd.f"
			itau = 1;
#line 2810 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2816 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2816 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2818 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2824 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2824 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2826 "zgesvd.f"
			ie = 1;
#line 2827 "zgesvd.f"
			itauq = itau;
#line 2828 "zgesvd.f"
			itaup = itauq + *m;
#line 2829 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2833 "zgesvd.f"
			i__2 = *m - 1;
#line 2833 "zgesvd.f"
			i__3 = *m - 1;
#line 2833 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2840 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2840 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2848 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2848 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2856 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2856 "zgesvd.f"
			zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2858 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2866 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2870 "zgesvd.f"
		    }

#line 2872 "zgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

#line 2879 "zgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2883 "zgesvd.f"
			iu = 1;
#line 2884 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2888 "zgesvd.f"
			    ldwrku = *lda;
#line 2889 "zgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2893 "zgesvd.f"
			    ldwrku = *m;
#line 2894 "zgesvd.f"
			}
#line 2895 "zgesvd.f"
			itau = iu + ldwrku * *m;
#line 2896 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2902 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2902 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2907 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2909 "zgesvd.f"
			i__2 = *m - 1;
#line 2909 "zgesvd.f"
			i__3 = *m - 1;
#line 2909 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2916 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2916 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2918 "zgesvd.f"
			ie = 1;
#line 2919 "zgesvd.f"
			itauq = itau;
#line 2920 "zgesvd.f"
			itaup = itauq + *m;
#line 2921 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2927 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2927 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2931 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2939 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2939 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2947 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2947 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2949 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2957 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2966 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2969 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2973 "zgesvd.f"
			itau = 1;
#line 2974 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2980 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2980 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2982 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2988 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2988 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2993 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2994 "zgesvd.f"
			i__2 = *m - 1;
#line 2994 "zgesvd.f"
			i__3 = *m - 1;
#line 2994 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 2996 "zgesvd.f"
			ie = 1;
#line 2997 "zgesvd.f"
			itauq = itau;
#line 2998 "zgesvd.f"
			itaup = itauq + *m;
#line 2999 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3005 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3005 "zgesvd.f"
			zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3014 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3014 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3022 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3022 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3024 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3032 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3036 "zgesvd.f"
		    }

#line 3038 "zgesvd.f"
		}

#line 3040 "zgesvd.f"
	    } else if (wntva) {

#line 3042 "zgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 3048 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3048 "zgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3052 "zgesvd.f"
			ir = 1;
#line 3053 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 3057 "zgesvd.f"
			    ldwrkr = *lda;
#line 3058 "zgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 3062 "zgesvd.f"
			    ldwrkr = *m;
#line 3063 "zgesvd.f"
			}
#line 3064 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3065 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3071 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3071 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3073 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 3077 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 3079 "zgesvd.f"
			i__2 = *m - 1;
#line 3079 "zgesvd.f"
			i__3 = *m - 1;
#line 3079 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3086 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3086 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3088 "zgesvd.f"
			ie = 1;
#line 3089 "zgesvd.f"
			itauq = itau;
#line 3090 "zgesvd.f"
			itaup = itauq + *m;
#line 3091 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3097 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3097 "zgesvd.f"
			zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3107 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3107 "zgesvd.f"
			zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3110 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3117 "zgesvd.f"
			zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3126 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3131 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3133 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3137 "zgesvd.f"
			itau = 1;
#line 3138 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3144 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3144 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3146 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3152 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3152 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3154 "zgesvd.f"
			ie = 1;
#line 3155 "zgesvd.f"
			itauq = itau;
#line 3156 "zgesvd.f"
			itaup = itauq + *m;
#line 3157 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3161 "zgesvd.f"
			i__2 = *m - 1;
#line 3161 "zgesvd.f"
			i__3 = *m - 1;
#line 3161 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3168 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3168 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3177 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3177 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3180 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3187 "zgesvd.f"
			zbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 3191 "zgesvd.f"
		    }

#line 3193 "zgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3199 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3199 "zgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3203 "zgesvd.f"
			iu = 1;
#line 3204 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3208 "zgesvd.f"
			    ldwrku = *lda;
#line 3209 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3210 "zgesvd.f"
			    ldwrkr = *lda;
#line 3211 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3215 "zgesvd.f"
			    ldwrku = *lda;
#line 3216 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3217 "zgesvd.f"
			    ldwrkr = *m;
#line 3218 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3222 "zgesvd.f"
			    ldwrku = *m;
#line 3223 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3224 "zgesvd.f"
			    ldwrkr = *m;
#line 3225 "zgesvd.f"
			}
#line 3226 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3227 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3233 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3233 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3235 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3241 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3241 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3246 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3248 "zgesvd.f"
			i__2 = *m - 1;
#line 3248 "zgesvd.f"
			i__3 = *m - 1;
#line 3248 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3250 "zgesvd.f"
			ie = 1;
#line 3251 "zgesvd.f"
			itauq = itau;
#line 3252 "zgesvd.f"
			itaup = itauq + *m;
#line 3253 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 3261 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3261 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3265 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3273 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3273 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3281 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3281 "zgesvd.f"
			zungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3284 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3292 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3302 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3307 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3311 "zgesvd.f"
			zlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3314 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3318 "zgesvd.f"
			itau = 1;
#line 3319 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3325 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3325 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3327 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3333 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3333 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3335 "zgesvd.f"
			ie = 1;
#line 3336 "zgesvd.f"
			itauq = itau;
#line 3337 "zgesvd.f"
			itaup = itauq + *m;
#line 3338 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3342 "zgesvd.f"
			i__2 = *m - 1;
#line 3342 "zgesvd.f"
			i__3 = *m - 1;
#line 3342 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3349 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3349 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3358 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3358 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3366 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3366 "zgesvd.f"
			zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3368 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3376 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3380 "zgesvd.f"
		    }

#line 3382 "zgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3389 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3389 "zgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3393 "zgesvd.f"
			iu = 1;
#line 3394 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3398 "zgesvd.f"
			    ldwrku = *lda;
#line 3399 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3403 "zgesvd.f"
			    ldwrku = *m;
#line 3404 "zgesvd.f"
			}
#line 3405 "zgesvd.f"
			itau = iu + ldwrku * *m;
#line 3406 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3412 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3412 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3414 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3420 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3420 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3425 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3427 "zgesvd.f"
			i__2 = *m - 1;
#line 3427 "zgesvd.f"
			i__3 = *m - 1;
#line 3427 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3429 "zgesvd.f"
			ie = 1;
#line 3430 "zgesvd.f"
			itauq = itau;
#line 3431 "zgesvd.f"
			itaup = itauq + *m;
#line 3432 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3438 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3438 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3442 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3449 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3449 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3457 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3457 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3459 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3467 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3476 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3481 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3483 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3487 "zgesvd.f"
			itau = 1;
#line 3488 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3494 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3494 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3496 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3502 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3502 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3507 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3508 "zgesvd.f"
			i__2 = *m - 1;
#line 3508 "zgesvd.f"
			i__3 = *m - 1;
#line 3508 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3510 "zgesvd.f"
			ie = 1;
#line 3511 "zgesvd.f"
			itauq = itau;
#line 3512 "zgesvd.f"
			itaup = itauq + *m;
#line 3513 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3519 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3519 "zgesvd.f"
			zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3528 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3528 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3536 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3536 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3538 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3546 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3550 "zgesvd.f"
		    }

#line 3552 "zgesvd.f"
		}

#line 3554 "zgesvd.f"
	    }

#line 3556 "zgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3563 "zgesvd.f"
	    ie = 1;
#line 3564 "zgesvd.f"
	    itauq = 1;
#line 3565 "zgesvd.f"
	    itaup = itauq + *m;
#line 3566 "zgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 3572 "zgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3572 "zgesvd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 3575 "zgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3582 "zgesvd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3583 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3583 "zgesvd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3585 "zgesvd.f"
	    }
#line 3586 "zgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB) */
/*              (RWorkspace: 0) */

#line 3593 "zgesvd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3594 "zgesvd.f"
		if (wntva) {
#line 3594 "zgesvd.f"
		    nrvt = *n;
#line 3594 "zgesvd.f"
		}
#line 3596 "zgesvd.f"
		if (wntvs) {
#line 3596 "zgesvd.f"
		    nrvt = *m;
#line 3596 "zgesvd.f"
		}
#line 3598 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3598 "zgesvd.f"
		zungbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3600 "zgesvd.f"
	    }
#line 3601 "zgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3608 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3608 "zgesvd.f"
		zungbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3610 "zgesvd.f"
	    }
#line 3611 "zgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 3618 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3618 "zgesvd.f"
		zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3620 "zgesvd.f"
	    }
#line 3621 "zgesvd.f"
	    irwork = ie + *m;
#line 3622 "zgesvd.f"
	    if (wntuas || wntuo) {
#line 3622 "zgesvd.f"
		nru = *m;
#line 3622 "zgesvd.f"
	    }
#line 3624 "zgesvd.f"
	    if (wntun) {
#line 3624 "zgesvd.f"
		nru = 0;
#line 3624 "zgesvd.f"
	    }
#line 3626 "zgesvd.f"
	    if (wntvas || wntvo) {
#line 3626 "zgesvd.f"
		ncvt = *n;
#line 3626 "zgesvd.f"
	    }
#line 3628 "zgesvd.f"
	    if (wntvn) {
#line 3628 "zgesvd.f"
		ncvt = 0;
#line 3628 "zgesvd.f"
	    }
#line 3630 "zgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3638 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3641 "zgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3649 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3652 "zgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3660 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3663 "zgesvd.f"
	    }

#line 3665 "zgesvd.f"
	}

#line 3667 "zgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3671 "zgesvd.f"
    if (iscl == 1) {
#line 3672 "zgesvd.f"
	if (anrm > bignum) {
#line 3672 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3672 "zgesvd.f"
	}
#line 3675 "zgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3675 "zgesvd.f"
	    i__2 = minmn - 1;
#line 3675 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3675 "zgesvd.f"
	}
#line 3678 "zgesvd.f"
	if (anrm < smlnum) {
#line 3678 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3678 "zgesvd.f"
	}
#line 3681 "zgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3681 "zgesvd.f"
	    i__2 = minmn - 1;
#line 3681 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3681 "zgesvd.f"
	}
#line 3684 "zgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3688 "zgesvd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 3690 "zgesvd.f"
    return 0;

/*     End of ZGESVD */

} /* zgesvd_ */

