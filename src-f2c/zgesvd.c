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
	     zgebrd_(integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zlascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), zgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *), zlacpy_(
	    char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen), zlaset_(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen);
    static integer ldwrkr;
    extern /* Subroutine */ int zbdsqr_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen);
    static integer minwrk, ldwrku, maxwrk;
    extern /* Subroutine */ int zungbr_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *, ftnlen);
    static doublereal smlnum;
    static integer irwork;
    extern /* Subroutine */ int zunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen), zunglq_(integer *, integer *, integer *
	    , doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static logical lquery, wntuas, wntvas;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer lwork_zungbr_p__, lwork_zungbr_q__, lwork_zunglq_m__, 
	    lwork_zunglq_n__, lwork_zungqr_m__, lwork_zungqr_n__;


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
	    zgeqrf_(m, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 325 "zgesvd.f"
	    lwork_zgeqrf__ = (integer) cdum[0].r;
/*           Compute space needed for ZUNGQR */
#line 327 "zgesvd.f"
	    zungqr_(m, n, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 328 "zgesvd.f"
	    lwork_zungqr_n__ = (integer) cdum[0].r;
#line 329 "zgesvd.f"
	    zungqr_(m, m, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 330 "zgesvd.f"
	    lwork_zungqr_m__ = (integer) cdum[0].r;
/*           Compute space needed for ZGEBRD */
#line 332 "zgesvd.f"
	    zgebrd_(n, n, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum, &
		    c_n1, &ierr);
#line 334 "zgesvd.f"
	    lwork_zgebrd__ = (integer) cdum[0].r;
/*           Compute space needed for ZUNGBR */
#line 336 "zgesvd.f"
	    zungbr_("P", n, n, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr,
		     (ftnlen)1);
#line 338 "zgesvd.f"
	    lwork_zungbr_p__ = (integer) cdum[0].r;
#line 339 "zgesvd.f"
	    zungbr_("Q", n, n, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr,
		     (ftnlen)1);
#line 341 "zgesvd.f"
	    lwork_zungbr_q__ = (integer) cdum[0].r;

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
		zgebrd_(m, n, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum,
			 &c_n1, &ierr);
#line 448 "zgesvd.f"
		lwork_zgebrd__ = (integer) cdum[0].r;
#line 449 "zgesvd.f"
		maxwrk = (*n << 1) + lwork_zgebrd__;
#line 450 "zgesvd.f"
		if (wntus || wntuo) {
#line 451 "zgesvd.f"
		    zungbr_("Q", m, n, n, &a[a_offset], lda, cdum, cdum, &
			    c_n1, &ierr, (ftnlen)1);
#line 453 "zgesvd.f"
		    lwork_zungbr_q__ = (integer) cdum[0].r;
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
		    zungbr_("Q", m, m, n, &a[a_offset], lda, cdum, cdum, &
			    c_n1, &ierr, (ftnlen)1);
#line 459 "zgesvd.f"
		    lwork_zungbr_q__ = (integer) cdum[0].r;
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
		}
#line 465 "zgesvd.f"
		minwrk = (*n << 1) + *m;
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
	    zgelqf_(m, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 474 "zgesvd.f"
	    lwork_zgelqf__ = (integer) cdum[0].r;
/*           Compute space needed for ZUNGLQ */
#line 476 "zgesvd.f"
	    zunglq_(n, n, m, cdum, n, cdum, cdum, &c_n1, &ierr);
#line 478 "zgesvd.f"
	    lwork_zunglq_n__ = (integer) cdum[0].r;
#line 479 "zgesvd.f"
	    zunglq_(m, n, m, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 480 "zgesvd.f"
	    lwork_zunglq_m__ = (integer) cdum[0].r;
/*           Compute space needed for ZGEBRD */
#line 482 "zgesvd.f"
	    zgebrd_(m, m, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum, &
		    c_n1, &ierr);
#line 484 "zgesvd.f"
	    lwork_zgebrd__ = (integer) cdum[0].r;
/*            Compute space needed for ZUNGBR P */
#line 486 "zgesvd.f"
	    zungbr_("P", m, m, m, &a[a_offset], n, cdum, cdum, &c_n1, &ierr, (
		    ftnlen)1);
#line 488 "zgesvd.f"
	    lwork_zungbr_p__ = (integer) cdum[0].r;
/*           Compute space needed for ZUNGBR Q */
#line 490 "zgesvd.f"
	    zungbr_("Q", m, m, m, &a[a_offset], n, cdum, cdum, &c_n1, &ierr, (
		    ftnlen)1);
#line 492 "zgesvd.f"
	    lwork_zungbr_q__ = (integer) cdum[0].r;
#line 493 "zgesvd.f"
	    if (*n >= mnthr) {
#line 494 "zgesvd.f"
		if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N') */

#line 498 "zgesvd.f"
		    maxwrk = *m + lwork_zgelqf__;
/* Computing MAX */
#line 499 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zgebrd__;
#line 499 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 500 "zgesvd.f"
		    if (wntuo || wntuas) {
/* Computing MAX */
#line 500 "zgesvd.f"
			i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 500 "zgesvd.f"
			maxwrk = max(i__2,i__3);
#line 500 "zgesvd.f"
		    }
#line 502 "zgesvd.f"
		    minwrk = *m * 3;
#line 503 "zgesvd.f"
		} else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O') */

#line 507 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 508 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 508 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 509 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 509 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 510 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 510 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 511 "zgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 511 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 512 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 513 "zgesvd.f"
		} else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='O') */

#line 518 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 519 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 519 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 520 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 520 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 521 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 521 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 522 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 522 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 523 "zgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 523 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 524 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 525 "zgesvd.f"
		} else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */

#line 529 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 530 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 530 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 531 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 531 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 532 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 532 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 533 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 534 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 535 "zgesvd.f"
		} else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */

#line 539 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 540 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 540 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 541 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 541 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 542 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 542 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 543 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 543 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 544 "zgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 545 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 546 "zgesvd.f"
		} else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='S') */

#line 551 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 552 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_m__;
#line 552 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 553 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 553 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 554 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 554 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 555 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 555 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 556 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 557 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 558 "zgesvd.f"
		} else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */

#line 562 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 563 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_n__;
#line 563 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 564 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 564 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 565 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 565 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 566 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 567 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 568 "zgesvd.f"
		} else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */

#line 572 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 573 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_n__;
#line 573 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 574 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 574 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 575 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 575 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 576 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 576 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 577 "zgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 578 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 579 "zgesvd.f"
		} else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='A') */

#line 584 "zgesvd.f"
		    wrkbl = *m + lwork_zgelqf__;
/* Computing MAX */
#line 585 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_zunglq_n__;
#line 585 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 586 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zgebrd__;
#line 586 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 587 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 587 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 588 "zgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 588 "zgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 589 "zgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 590 "zgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 591 "zgesvd.f"
		}
#line 592 "zgesvd.f"
	    } else {

/*              Path 10t(N greater than M, but not much larger) */

#line 596 "zgesvd.f"
		zgebrd_(m, n, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum,
			 &c_n1, &ierr);
#line 598 "zgesvd.f"
		lwork_zgebrd__ = (integer) cdum[0].r;
#line 599 "zgesvd.f"
		maxwrk = (*m << 1) + lwork_zgebrd__;
#line 600 "zgesvd.f"
		if (wntvs || wntvo) {
/*                Compute space needed for ZUNGBR P */
#line 602 "zgesvd.f"
		    zungbr_("P", m, n, m, &a[a_offset], n, cdum, cdum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 604 "zgesvd.f"
		    lwork_zungbr_p__ = (integer) cdum[0].r;
/* Computing MAX */
#line 605 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 605 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 606 "zgesvd.f"
		}
#line 607 "zgesvd.f"
		if (wntva) {
#line 608 "zgesvd.f"
		    zungbr_("P", n, n, m, &a[a_offset], n, cdum, cdum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 610 "zgesvd.f"
		    lwork_zungbr_p__ = (integer) cdum[0].r;
/* Computing MAX */
#line 611 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_p__;
#line 611 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 612 "zgesvd.f"
		}
#line 613 "zgesvd.f"
		if (! wntun) {
/* Computing MAX */
#line 614 "zgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_zungbr_q__;
#line 614 "zgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 615 "zgesvd.f"
		}
#line 616 "zgesvd.f"
		minwrk = (*m << 1) + *n;
#line 617 "zgesvd.f"
	    }
#line 618 "zgesvd.f"
	}
#line 619 "zgesvd.f"
	maxwrk = max(maxwrk,minwrk);
#line 620 "zgesvd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 622 "zgesvd.f"
	if (*lwork < minwrk && ! lquery) {
#line 623 "zgesvd.f"
	    *info = -13;
#line 624 "zgesvd.f"
	}
#line 625 "zgesvd.f"
    }

#line 627 "zgesvd.f"
    if (*info != 0) {
#line 628 "zgesvd.f"
	i__2 = -(*info);
#line 628 "zgesvd.f"
	xerbla_("ZGESVD", &i__2, (ftnlen)6);
#line 629 "zgesvd.f"
	return 0;
#line 630 "zgesvd.f"
    } else if (lquery) {
#line 631 "zgesvd.f"
	return 0;
#line 632 "zgesvd.f"
    }

/*     Quick return if possible */

#line 636 "zgesvd.f"
    if (*m == 0 || *n == 0) {
#line 637 "zgesvd.f"
	return 0;
#line 638 "zgesvd.f"
    }

/*     Get machine constants */

#line 642 "zgesvd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 643 "zgesvd.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 644 "zgesvd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 648 "zgesvd.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 649 "zgesvd.f"
    iscl = 0;
#line 650 "zgesvd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 651 "zgesvd.f"
	iscl = 1;
#line 652 "zgesvd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 653 "zgesvd.f"
    } else if (anrm > bignum) {
#line 654 "zgesvd.f"
	iscl = 1;
#line 655 "zgesvd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 656 "zgesvd.f"
    }

#line 658 "zgesvd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 664 "zgesvd.f"
	if (*m >= mnthr) {

#line 666 "zgesvd.f"
	    if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N') */
/*              No left singular vectors to be computed */

#line 671 "zgesvd.f"
		itau = 1;
#line 672 "zgesvd.f"
		iwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: need 0) */

#line 678 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 678 "zgesvd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out below R */

#line 683 "zgesvd.f"
		if (*n > 1) {
#line 684 "zgesvd.f"
		    i__2 = *n - 1;
#line 684 "zgesvd.f"
		    i__3 = *n - 1;
#line 684 "zgesvd.f"
		    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 2], 
			    lda, (ftnlen)1);
#line 686 "zgesvd.f"
		}
#line 687 "zgesvd.f"
		ie = 1;
#line 688 "zgesvd.f"
		itauq = 1;
#line 689 "zgesvd.f"
		itaup = itauq + *n;
#line 690 "zgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 696 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 696 "zgesvd.f"
		zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 699 "zgesvd.f"
		ncvt = 0;
#line 700 "zgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 706 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 706 "zgesvd.f"
		    zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 708 "zgesvd.f"
		    ncvt = *n;
#line 709 "zgesvd.f"
		}
#line 710 "zgesvd.f"
		irwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 717 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			irwork], info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 722 "zgesvd.f"
		if (wntvas) {
#line 722 "zgesvd.f"
		    zlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 722 "zgesvd.f"
		}

#line 725 "zgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

#line 731 "zgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 735 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 736 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 736 "zgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 740 "zgesvd.f"
			ldwrku = *lda;
#line 741 "zgesvd.f"
			ldwrkr = *lda;
#line 742 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 742 "zgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 742 "zgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 746 "zgesvd.f"
			    ldwrku = *lda;
#line 747 "zgesvd.f"
			    ldwrkr = *n;
#line 748 "zgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 752 "zgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 753 "zgesvd.f"
			    ldwrkr = *n;
#line 754 "zgesvd.f"
			}
#line 754 "zgesvd.f"
		    }
#line 755 "zgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 756 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 762 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 762 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 767 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 768 "zgesvd.f"
		    i__2 = *n - 1;
#line 768 "zgesvd.f"
		    i__3 = *n - 1;
#line 768 "zgesvd.f"
		    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1], &
			    ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 775 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 775 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 777 "zgesvd.f"
		    ie = 1;
#line 778 "zgesvd.f"
		    itauq = itau;
#line 779 "zgesvd.f"
		    itaup = itauq + *n;
#line 780 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 786 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 786 "zgesvd.f"
		    zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: need 0) */

#line 794 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 794 "zgesvd.f"
		    zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 797 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 804 "zgesvd.f"
		    zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &work[ir], &ldwrkr, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 807 "zgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 814 "zgesvd.f"
		    i__2 = *m;
#line 814 "zgesvd.f"
		    i__3 = ldwrku;
#line 814 "zgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 815 "zgesvd.f"
			i__4 = *m - i__ + 1;
#line 815 "zgesvd.f"
			chunk = min(i__4,ldwrku);
#line 816 "zgesvd.f"
			zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 819 "zgesvd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 821 "zgesvd.f"
/* L10: */
#line 821 "zgesvd.f"
		    }

#line 823 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 827 "zgesvd.f"
		    ie = 1;
#line 828 "zgesvd.f"
		    itauq = 1;
#line 829 "zgesvd.f"
		    itaup = itauq + *n;
#line 830 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*                 (RWorkspace: N) */

#line 836 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 836 "zgesvd.f"
		    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 844 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 844 "zgesvd.f"
		    zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 846 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (CWorkspace: need 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 853 "zgesvd.f"
		    zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &a[a_offset], lda, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 856 "zgesvd.f"
		}

#line 858 "zgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 864 "zgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 868 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 869 "zgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 869 "zgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 873 "zgesvd.f"
			ldwrku = *lda;
#line 874 "zgesvd.f"
			ldwrkr = *lda;
#line 875 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 875 "zgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 875 "zgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 879 "zgesvd.f"
			    ldwrku = *lda;
#line 880 "zgesvd.f"
			    ldwrkr = *n;
#line 881 "zgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 885 "zgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 886 "zgesvd.f"
			    ldwrkr = *n;
#line 887 "zgesvd.f"
			}
#line 887 "zgesvd.f"
		    }
#line 888 "zgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 889 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 895 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 895 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 900 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 901 "zgesvd.f"
		    if (*n > 1) {
#line 901 "zgesvd.f"
			i__3 = *n - 1;
#line 901 "zgesvd.f"
			i__2 = *n - 1;
#line 901 "zgesvd.f"
			zlaset_("L", &i__3, &i__2, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 901 "zgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 909 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 909 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 911 "zgesvd.f"
		    ie = 1;
#line 912 "zgesvd.f"
		    itauq = itau;
#line 913 "zgesvd.f"
		    itaup = itauq + *n;
#line 914 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 920 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 920 "zgesvd.f"
		    zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 923 "zgesvd.f"
		    zlacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 929 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 929 "zgesvd.f"
		    zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 937 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 937 "zgesvd.f"
		    zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 939 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 947 "zgesvd.f"
		    zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, cdum, &c__1,
			     &rwork[irwork], info, (ftnlen)1);
#line 950 "zgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 957 "zgesvd.f"
		    i__3 = *m;
#line 957 "zgesvd.f"
		    i__2 = ldwrku;
#line 957 "zgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 958 "zgesvd.f"
			i__4 = *m - i__ + 1;
#line 958 "zgesvd.f"
			chunk = min(i__4,ldwrku);
#line 959 "zgesvd.f"
			zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 962 "zgesvd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 964 "zgesvd.f"
/* L20: */
#line 964 "zgesvd.f"
		    }

#line 966 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 970 "zgesvd.f"
		    itau = 1;
#line 971 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 977 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 977 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 982 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 983 "zgesvd.f"
		    if (*n > 1) {
#line 983 "zgesvd.f"
			i__2 = *n - 1;
#line 983 "zgesvd.f"
			i__3 = *n - 1;
#line 983 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 983 "zgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 991 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 991 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 993 "zgesvd.f"
		    ie = 1;
#line 994 "zgesvd.f"
		    itauq = itau;
#line 995 "zgesvd.f"
		    itaup = itauq + *n;
#line 996 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                 (RWorkspace: N) */

#line 1002 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1002 "zgesvd.f"
		    zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                 (RWorkspace: 0) */

#line 1010 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1010 "zgesvd.f"
		    zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 1018 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1018 "zgesvd.f"
		    zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1020 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 1028 "zgesvd.f"
		    zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, 
			    &rwork[irwork], info, (ftnlen)1);

#line 1032 "zgesvd.f"
		}

#line 1034 "zgesvd.f"
	    } else if (wntus) {

#line 1036 "zgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

#line 1042 "zgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1046 "zgesvd.f"
			ir = 1;
#line 1047 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1051 "zgesvd.f"
			    ldwrkr = *lda;
#line 1052 "zgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1056 "zgesvd.f"
			    ldwrkr = *n;
#line 1057 "zgesvd.f"
			}
#line 1058 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1059 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1065 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1065 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1070 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1072 "zgesvd.f"
			i__2 = *n - 1;
#line 1072 "zgesvd.f"
			i__3 = *n - 1;
#line 1072 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1079 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1079 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1081 "zgesvd.f"
			ie = 1;
#line 1082 "zgesvd.f"
			itauq = itau;
#line 1083 "zgesvd.f"
			itaup = itauq + *n;
#line 1084 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1090 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1090 "zgesvd.f"
			zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1099 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1099 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1102 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1109 "zgesvd.f"
			zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1118 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1121 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1125 "zgesvd.f"
			itau = 1;
#line 1126 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1132 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1132 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1134 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1140 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1140 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1142 "zgesvd.f"
			ie = 1;
#line 1143 "zgesvd.f"
			itauq = itau;
#line 1144 "zgesvd.f"
			itaup = itauq + *n;
#line 1145 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1149 "zgesvd.f"
			if (*n > 1) {
#line 1150 "zgesvd.f"
			    i__2 = *n - 1;
#line 1150 "zgesvd.f"
			    i__3 = *n - 1;
#line 1150 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1152 "zgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1158 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1158 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1166 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1166 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1169 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1176 "zgesvd.f"
			zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1180 "zgesvd.f"
		    }

#line 1182 "zgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

#line 1188 "zgesvd.f"
		    if (*lwork >= (*n << 1) * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1192 "zgesvd.f"
			iu = 1;
#line 1193 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1197 "zgesvd.f"
			    ldwrku = *lda;
#line 1198 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1199 "zgesvd.f"
			    ldwrkr = *lda;
#line 1200 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1204 "zgesvd.f"
			    ldwrku = *lda;
#line 1205 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1206 "zgesvd.f"
			    ldwrkr = *n;
#line 1207 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1211 "zgesvd.f"
			    ldwrku = *n;
#line 1212 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1213 "zgesvd.f"
			    ldwrkr = *n;
#line 1214 "zgesvd.f"
			}
#line 1215 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1216 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1222 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1222 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1227 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1229 "zgesvd.f"
			i__2 = *n - 1;
#line 1229 "zgesvd.f"
			i__3 = *n - 1;
#line 1229 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1236 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1236 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1238 "zgesvd.f"
			ie = 1;
#line 1239 "zgesvd.f"
			itauq = itau;
#line 1240 "zgesvd.f"
			itaup = itauq + *n;
#line 1241 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1249 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1249 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1253 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1260 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1260 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1269 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1269 "zgesvd.f"
			zungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1272 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1280 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1290 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1297 "zgesvd.f"
			zlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1300 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1304 "zgesvd.f"
			itau = 1;
#line 1305 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1311 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1311 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1313 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1319 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1319 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1321 "zgesvd.f"
			ie = 1;
#line 1322 "zgesvd.f"
			itauq = itau;
#line 1323 "zgesvd.f"
			itaup = itauq + *n;
#line 1324 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1328 "zgesvd.f"
			if (*n > 1) {
#line 1329 "zgesvd.f"
			    i__2 = *n - 1;
#line 1329 "zgesvd.f"
			    i__3 = *n - 1;
#line 1329 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1331 "zgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1337 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1337 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1345 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1345 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1353 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1353 "zgesvd.f"
			zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1355 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1363 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1367 "zgesvd.f"
		    }

#line 1369 "zgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

#line 1376 "zgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1380 "zgesvd.f"
			iu = 1;
#line 1381 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1385 "zgesvd.f"
			    ldwrku = *lda;
#line 1386 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1390 "zgesvd.f"
			    ldwrku = *n;
#line 1391 "zgesvd.f"
			}
#line 1392 "zgesvd.f"
			itau = iu + ldwrku * *n;
#line 1393 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1399 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1399 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1404 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1406 "zgesvd.f"
			i__2 = *n - 1;
#line 1406 "zgesvd.f"
			i__3 = *n - 1;
#line 1406 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1413 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1413 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1415 "zgesvd.f"
			ie = 1;
#line 1416 "zgesvd.f"
			itauq = itau;
#line 1417 "zgesvd.f"
			itaup = itauq + *n;
#line 1418 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1424 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1424 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1428 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1435 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1435 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1444 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1444 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1446 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1454 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1463 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1466 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1470 "zgesvd.f"
			itau = 1;
#line 1471 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1477 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1477 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1479 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1485 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1485 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1490 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1491 "zgesvd.f"
			if (*n > 1) {
#line 1491 "zgesvd.f"
			    i__2 = *n - 1;
#line 1491 "zgesvd.f"
			    i__3 = *n - 1;
#line 1491 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1491 "zgesvd.f"
			}
#line 1494 "zgesvd.f"
			ie = 1;
#line 1495 "zgesvd.f"
			itauq = itau;
#line 1496 "zgesvd.f"
			itaup = itauq + *n;
#line 1497 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1503 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1503 "zgesvd.f"
			zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1512 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1512 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1520 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1520 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1522 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1530 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1534 "zgesvd.f"
		    }

#line 1536 "zgesvd.f"
		}

#line 1538 "zgesvd.f"
	    } else if (wntua) {

#line 1540 "zgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1546 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1546 "zgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1550 "zgesvd.f"
			ir = 1;
#line 1551 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1555 "zgesvd.f"
			    ldwrkr = *lda;
#line 1556 "zgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1560 "zgesvd.f"
			    ldwrkr = *n;
#line 1561 "zgesvd.f"
			}
#line 1562 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1563 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1569 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1569 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1571 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1575 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1577 "zgesvd.f"
			i__2 = *n - 1;
#line 1577 "zgesvd.f"
			i__3 = *n - 1;
#line 1577 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1584 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1584 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1586 "zgesvd.f"
			ie = 1;
#line 1587 "zgesvd.f"
			itauq = itau;
#line 1588 "zgesvd.f"
			itaup = itauq + *n;
#line 1589 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1595 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1595 "zgesvd.f"
			zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1604 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1604 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1607 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1614 "zgesvd.f"
			zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1623 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1628 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1630 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1634 "zgesvd.f"
			itau = 1;
#line 1635 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1641 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1641 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1643 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1649 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1649 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1651 "zgesvd.f"
			ie = 1;
#line 1652 "zgesvd.f"
			itauq = itau;
#line 1653 "zgesvd.f"
			itaup = itauq + *n;
#line 1654 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1658 "zgesvd.f"
			if (*n > 1) {
#line 1659 "zgesvd.f"
			    i__2 = *n - 1;
#line 1659 "zgesvd.f"
			    i__3 = *n - 1;
#line 1659 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1661 "zgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1667 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1667 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1676 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1676 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1679 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1686 "zgesvd.f"
			zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1690 "zgesvd.f"
		    }

#line 1692 "zgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1698 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1698 "zgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1702 "zgesvd.f"
			iu = 1;
#line 1703 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1707 "zgesvd.f"
			    ldwrku = *lda;
#line 1708 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1709 "zgesvd.f"
			    ldwrkr = *lda;
#line 1710 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1714 "zgesvd.f"
			    ldwrku = *lda;
#line 1715 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1716 "zgesvd.f"
			    ldwrkr = *n;
#line 1717 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1721 "zgesvd.f"
			    ldwrku = *n;
#line 1722 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1723 "zgesvd.f"
			    ldwrkr = *n;
#line 1724 "zgesvd.f"
			}
#line 1725 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1726 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1732 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1732 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1734 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1740 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1740 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1745 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1747 "zgesvd.f"
			i__2 = *n - 1;
#line 1747 "zgesvd.f"
			i__3 = *n - 1;
#line 1747 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1749 "zgesvd.f"
			ie = 1;
#line 1750 "zgesvd.f"
			itauq = itau;
#line 1751 "zgesvd.f"
			itaup = itauq + *n;
#line 1752 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1760 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1760 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1764 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1771 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1771 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1780 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1780 "zgesvd.f"
			zungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1783 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1791 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1801 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1806 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1810 "zgesvd.f"
			zlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1813 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1817 "zgesvd.f"
			itau = 1;
#line 1818 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1824 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1824 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1826 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1832 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1832 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1834 "zgesvd.f"
			ie = 1;
#line 1835 "zgesvd.f"
			itauq = itau;
#line 1836 "zgesvd.f"
			itaup = itauq + *n;
#line 1837 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1841 "zgesvd.f"
			if (*n > 1) {
#line 1842 "zgesvd.f"
			    i__2 = *n - 1;
#line 1842 "zgesvd.f"
			    i__3 = *n - 1;
#line 1842 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1844 "zgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1850 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1850 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1859 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1859 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1867 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1867 "zgesvd.f"
			zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1869 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1877 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1881 "zgesvd.f"
		    }

#line 1883 "zgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1890 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1890 "zgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1894 "zgesvd.f"
			iu = 1;
#line 1895 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1899 "zgesvd.f"
			    ldwrku = *lda;
#line 1900 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1904 "zgesvd.f"
			    ldwrku = *n;
#line 1905 "zgesvd.f"
			}
#line 1906 "zgesvd.f"
			itau = iu + ldwrku * *n;
#line 1907 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1913 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1913 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1915 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1921 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1921 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1926 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1928 "zgesvd.f"
			i__2 = *n - 1;
#line 1928 "zgesvd.f"
			i__3 = *n - 1;
#line 1928 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1930 "zgesvd.f"
			ie = 1;
#line 1931 "zgesvd.f"
			itauq = itau;
#line 1932 "zgesvd.f"
			itaup = itauq + *n;
#line 1933 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1939 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1939 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1943 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1950 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1950 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: need   0) */

#line 1959 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1959 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1961 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1969 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1978 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1983 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1985 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1989 "zgesvd.f"
			itau = 1;
#line 1990 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1996 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1996 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1998 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2004 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2004 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 2009 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 2010 "zgesvd.f"
			if (*n > 1) {
#line 2010 "zgesvd.f"
			    i__2 = *n - 1;
#line 2010 "zgesvd.f"
			    i__3 = *n - 1;
#line 2010 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 2010 "zgesvd.f"
			}
#line 2013 "zgesvd.f"
			ie = 1;
#line 2014 "zgesvd.f"
			itauq = itau;
#line 2015 "zgesvd.f"
			itaup = itauq + *n;
#line 2016 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 2022 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2022 "zgesvd.f"
			zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2031 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2031 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2039 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2039 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 2041 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2049 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2053 "zgesvd.f"
		    }

#line 2055 "zgesvd.f"
		}

#line 2057 "zgesvd.f"
	    }

#line 2059 "zgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 2066 "zgesvd.f"
	    ie = 1;
#line 2067 "zgesvd.f"
	    itauq = 1;
#line 2068 "zgesvd.f"
	    itaup = itauq + *n;
#line 2069 "zgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 2075 "zgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 2075 "zgesvd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 2078 "zgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB) */
/*              (RWorkspace: 0) */

#line 2085 "zgesvd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 2086 "zgesvd.f"
		if (wntus) {
#line 2086 "zgesvd.f"
		    ncu = *n;
#line 2086 "zgesvd.f"
		}
#line 2088 "zgesvd.f"
		if (wntua) {
#line 2088 "zgesvd.f"
		    ncu = *m;
#line 2088 "zgesvd.f"
		}
#line 2090 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2090 "zgesvd.f"
		zungbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2092 "zgesvd.f"
	    }
#line 2093 "zgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2100 "zgesvd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2101 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2101 "zgesvd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2103 "zgesvd.f"
	    }
#line 2104 "zgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 2111 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2111 "zgesvd.f"
		zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2113 "zgesvd.f"
	    }
#line 2114 "zgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2121 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2121 "zgesvd.f"
		zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2123 "zgesvd.f"
	    }
#line 2124 "zgesvd.f"
	    irwork = ie + *n;
#line 2125 "zgesvd.f"
	    if (wntuas || wntuo) {
#line 2125 "zgesvd.f"
		nru = *m;
#line 2125 "zgesvd.f"
	    }
#line 2127 "zgesvd.f"
	    if (wntun) {
#line 2127 "zgesvd.f"
		nru = 0;
#line 2127 "zgesvd.f"
	    }
#line 2129 "zgesvd.f"
	    if (wntvas || wntvo) {
#line 2129 "zgesvd.f"
		ncvt = *n;
#line 2129 "zgesvd.f"
	    }
#line 2131 "zgesvd.f"
	    if (wntvn) {
#line 2131 "zgesvd.f"
		ncvt = 0;
#line 2131 "zgesvd.f"
	    }
#line 2133 "zgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2141 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2144 "zgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2152 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2155 "zgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2163 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2166 "zgesvd.f"
	    }

#line 2168 "zgesvd.f"
	}

#line 2170 "zgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2176 "zgesvd.f"
	if (*n >= mnthr) {

#line 2178 "zgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2183 "zgesvd.f"
		itau = 1;
#line 2184 "zgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 2190 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2190 "zgesvd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2195 "zgesvd.f"
		i__2 = *m - 1;
#line 2195 "zgesvd.f"
		i__3 = *m - 1;
#line 2195 "zgesvd.f"
		zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 2197 "zgesvd.f"
		ie = 1;
#line 2198 "zgesvd.f"
		itauq = 1;
#line 2199 "zgesvd.f"
		itaup = itauq + *m;
#line 2200 "zgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 2206 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2206 "zgesvd.f"
		zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2209 "zgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2215 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2215 "zgesvd.f"
		    zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2217 "zgesvd.f"
		}
#line 2218 "zgesvd.f"
		irwork = ie + *m;
#line 2219 "zgesvd.f"
		nru = 0;
#line 2220 "zgesvd.f"
		if (wntuo || wntuas) {
#line 2220 "zgesvd.f"
		    nru = *m;
#line 2220 "zgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2228 "zgesvd.f"
		zbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &rwork[ie], cdum, &
			c__1, &a[a_offset], lda, cdum, &c__1, &rwork[irwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2233 "zgesvd.f"
		if (wntuas) {
#line 2233 "zgesvd.f"
		    zlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2233 "zgesvd.f"
		}

#line 2236 "zgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

#line 2242 "zgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2246 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2247 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 2247 "zgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2251 "zgesvd.f"
			ldwrku = *lda;
#line 2252 "zgesvd.f"
			chunk = *n;
#line 2253 "zgesvd.f"
			ldwrkr = *lda;
#line 2254 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2254 "zgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 2254 "zgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2258 "zgesvd.f"
			    ldwrku = *lda;
#line 2259 "zgesvd.f"
			    chunk = *n;
#line 2260 "zgesvd.f"
			    ldwrkr = *m;
#line 2261 "zgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2265 "zgesvd.f"
			    ldwrku = *m;
#line 2266 "zgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2267 "zgesvd.f"
			    ldwrkr = *m;
#line 2268 "zgesvd.f"
			}
#line 2268 "zgesvd.f"
		    }
#line 2269 "zgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2270 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2276 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2276 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2281 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2282 "zgesvd.f"
		    i__2 = *m - 1;
#line 2282 "zgesvd.f"
		    i__3 = *m - 1;
#line 2282 "zgesvd.f"
		    zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2289 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2289 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2291 "zgesvd.f"
		    ie = 1;
#line 2292 "zgesvd.f"
		    itauq = itau;
#line 2293 "zgesvd.f"
		    itaup = itauq + *m;
#line 2294 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2300 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2300 "zgesvd.f"
		    zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2308 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2308 "zgesvd.f"
		    zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2311 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2318 "zgesvd.f"
		    zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &work[
			    ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2321 "zgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N) */
/*                 (RWorkspace: 0) */

#line 2328 "zgesvd.f"
		    i__2 = *n;
#line 2328 "zgesvd.f"
		    i__3 = chunk;
#line 2328 "zgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2329 "zgesvd.f"
			i__4 = *n - i__ + 1;
#line 2329 "zgesvd.f"
			blk = min(i__4,chunk);
#line 2330 "zgesvd.f"
			zgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2333 "zgesvd.f"
			zlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2335 "zgesvd.f"
/* L30: */
#line 2335 "zgesvd.f"
		    }

#line 2337 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2341 "zgesvd.f"
		    ie = 1;
#line 2342 "zgesvd.f"
		    itauq = 1;
#line 2343 "zgesvd.f"
		    itaup = itauq + *m;
#line 2344 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*                 (RWorkspace: need M) */

#line 2350 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2350 "zgesvd.f"
		    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2358 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2358 "zgesvd.f"
		    zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2360 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2367 "zgesvd.f"
		    zbdsqr_("L", m, n, &c__0, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 2370 "zgesvd.f"
		}

#line 2372 "zgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 2378 "zgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2382 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2383 "zgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 2383 "zgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2387 "zgesvd.f"
			ldwrku = *lda;
#line 2388 "zgesvd.f"
			chunk = *n;
#line 2389 "zgesvd.f"
			ldwrkr = *lda;
#line 2390 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2390 "zgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 2390 "zgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2394 "zgesvd.f"
			    ldwrku = *lda;
#line 2395 "zgesvd.f"
			    chunk = *n;
#line 2396 "zgesvd.f"
			    ldwrkr = *m;
#line 2397 "zgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2401 "zgesvd.f"
			    ldwrku = *m;
#line 2402 "zgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2403 "zgesvd.f"
			    ldwrkr = *m;
#line 2404 "zgesvd.f"
			}
#line 2404 "zgesvd.f"
		    }
#line 2405 "zgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2406 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2412 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2412 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2417 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2418 "zgesvd.f"
		    i__3 = *m - 1;
#line 2418 "zgesvd.f"
		    i__2 = *m - 1;
#line 2418 "zgesvd.f"
		    zlaset_("U", &i__3, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2425 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2425 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2427 "zgesvd.f"
		    ie = 1;
#line 2428 "zgesvd.f"
		    itauq = itau;
#line 2429 "zgesvd.f"
		    itaup = itauq + *m;
#line 2430 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2436 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2436 "zgesvd.f"
		    zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2439 "zgesvd.f"
		    zlacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2445 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2445 "zgesvd.f"
		    zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2453 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2453 "zgesvd.f"
		    zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2455 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2463 "zgesvd.f"
		    zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[ir],
			     &ldwrkr, &u[u_offset], ldu, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2466 "zgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N)) */
/*                 (RWorkspace: 0) */

#line 2473 "zgesvd.f"
		    i__3 = *n;
#line 2473 "zgesvd.f"
		    i__2 = chunk;
#line 2473 "zgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2474 "zgesvd.f"
			i__4 = *n - i__ + 1;
#line 2474 "zgesvd.f"
			blk = min(i__4,chunk);
#line 2475 "zgesvd.f"
			zgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2478 "zgesvd.f"
			zlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2480 "zgesvd.f"
/* L40: */
#line 2480 "zgesvd.f"
		    }

#line 2482 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2486 "zgesvd.f"
		    itau = 1;
#line 2487 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2493 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2493 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2498 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2499 "zgesvd.f"
		    i__2 = *m - 1;
#line 2499 "zgesvd.f"
		    i__3 = *m - 1;
#line 2499 "zgesvd.f"
		    zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2506 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2506 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2508 "zgesvd.f"
		    ie = 1;
#line 2509 "zgesvd.f"
		    itauq = itau;
#line 2510 "zgesvd.f"
		    itaup = itauq + *m;
#line 2511 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2517 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2517 "zgesvd.f"
		    zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                 (RWorkspace: 0) */

#line 2525 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2525 "zgesvd.f"
		    zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2533 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2533 "zgesvd.f"
		    zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2535 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2543 "zgesvd.f"
		    zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			    rwork[irwork], info, (ftnlen)1);

#line 2546 "zgesvd.f"
		}

#line 2548 "zgesvd.f"
	    } else if (wntvs) {

#line 2550 "zgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

#line 2556 "zgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2560 "zgesvd.f"
			ir = 1;
#line 2561 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2565 "zgesvd.f"
			    ldwrkr = *lda;
#line 2566 "zgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2570 "zgesvd.f"
			    ldwrkr = *m;
#line 2571 "zgesvd.f"
			}
#line 2572 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2573 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2579 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2579 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2584 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2586 "zgesvd.f"
			i__2 = *m - 1;
#line 2586 "zgesvd.f"
			i__3 = *m - 1;
#line 2586 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2593 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2593 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2595 "zgesvd.f"
			ie = 1;
#line 2596 "zgesvd.f"
			itauq = itau;
#line 2597 "zgesvd.f"
			itaup = itauq + *m;
#line 2598 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2604 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2604 "zgesvd.f"
			zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2614 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2614 "zgesvd.f"
			zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2617 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2624 "zgesvd.f"
			zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2633 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2636 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2640 "zgesvd.f"
			itau = 1;
#line 2641 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2647 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2647 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2652 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2658 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2658 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2660 "zgesvd.f"
			ie = 1;
#line 2661 "zgesvd.f"
			itauq = itau;
#line 2662 "zgesvd.f"
			itaup = itauq + *m;
#line 2663 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2667 "zgesvd.f"
			i__2 = *m - 1;
#line 2667 "zgesvd.f"
			i__3 = *m - 1;
#line 2667 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2674 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2674 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2682 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2682 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2685 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2692 "zgesvd.f"
			zbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 2696 "zgesvd.f"
		    }

#line 2698 "zgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

#line 2704 "zgesvd.f"
		    if (*lwork >= (*m << 1) * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2708 "zgesvd.f"
			iu = 1;
#line 2709 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2713 "zgesvd.f"
			    ldwrku = *lda;
#line 2714 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2715 "zgesvd.f"
			    ldwrkr = *lda;
#line 2716 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2720 "zgesvd.f"
			    ldwrku = *lda;
#line 2721 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2722 "zgesvd.f"
			    ldwrkr = *m;
#line 2723 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2727 "zgesvd.f"
			    ldwrku = *m;
#line 2728 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2729 "zgesvd.f"
			    ldwrkr = *m;
#line 2730 "zgesvd.f"
			}
#line 2731 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2732 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2738 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2738 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2743 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2745 "zgesvd.f"
			i__2 = *m - 1;
#line 2745 "zgesvd.f"
			i__3 = *m - 1;
#line 2745 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2752 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2752 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2754 "zgesvd.f"
			ie = 1;
#line 2755 "zgesvd.f"
			itauq = itau;
#line 2756 "zgesvd.f"
			itaup = itauq + *m;
#line 2757 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 2765 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2765 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2769 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2777 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2777 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2785 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2785 "zgesvd.f"
			zungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2788 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2796 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2806 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2813 "zgesvd.f"
			zlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2816 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2820 "zgesvd.f"
			itau = 1;
#line 2821 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2827 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2827 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2829 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2835 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2835 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2837 "zgesvd.f"
			ie = 1;
#line 2838 "zgesvd.f"
			itauq = itau;
#line 2839 "zgesvd.f"
			itaup = itauq + *m;
#line 2840 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2844 "zgesvd.f"
			i__2 = *m - 1;
#line 2844 "zgesvd.f"
			i__3 = *m - 1;
#line 2844 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2851 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2851 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2859 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2859 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2867 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2867 "zgesvd.f"
			zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2869 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2877 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2881 "zgesvd.f"
		    }

#line 2883 "zgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

#line 2890 "zgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2894 "zgesvd.f"
			iu = 1;
#line 2895 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2899 "zgesvd.f"
			    ldwrku = *lda;
#line 2900 "zgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2904 "zgesvd.f"
			    ldwrku = *m;
#line 2905 "zgesvd.f"
			}
#line 2906 "zgesvd.f"
			itau = iu + ldwrku * *m;
#line 2907 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2913 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2913 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2918 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2920 "zgesvd.f"
			i__2 = *m - 1;
#line 2920 "zgesvd.f"
			i__3 = *m - 1;
#line 2920 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2927 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2927 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2929 "zgesvd.f"
			ie = 1;
#line 2930 "zgesvd.f"
			itauq = itau;
#line 2931 "zgesvd.f"
			itaup = itauq + *m;
#line 2932 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2938 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2938 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2942 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2950 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2950 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2958 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2958 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2960 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2968 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2977 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2980 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2984 "zgesvd.f"
			itau = 1;
#line 2985 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2991 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2991 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2993 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2999 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2999 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3004 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3005 "zgesvd.f"
			i__2 = *m - 1;
#line 3005 "zgesvd.f"
			i__3 = *m - 1;
#line 3005 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3007 "zgesvd.f"
			ie = 1;
#line 3008 "zgesvd.f"
			itauq = itau;
#line 3009 "zgesvd.f"
			itaup = itauq + *m;
#line 3010 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3016 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3016 "zgesvd.f"
			zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3025 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3025 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3033 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3033 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3035 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3043 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3047 "zgesvd.f"
		    }

#line 3049 "zgesvd.f"
		}

#line 3051 "zgesvd.f"
	    } else if (wntva) {

#line 3053 "zgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 3059 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3059 "zgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3063 "zgesvd.f"
			ir = 1;
#line 3064 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 3068 "zgesvd.f"
			    ldwrkr = *lda;
#line 3069 "zgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 3073 "zgesvd.f"
			    ldwrkr = *m;
#line 3074 "zgesvd.f"
			}
#line 3075 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3076 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3082 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3082 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3084 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 3088 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 3090 "zgesvd.f"
			i__2 = *m - 1;
#line 3090 "zgesvd.f"
			i__3 = *m - 1;
#line 3090 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3097 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3097 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3099 "zgesvd.f"
			ie = 1;
#line 3100 "zgesvd.f"
			itauq = itau;
#line 3101 "zgesvd.f"
			itaup = itauq + *m;
#line 3102 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3108 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3108 "zgesvd.f"
			zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3118 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3118 "zgesvd.f"
			zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3121 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3128 "zgesvd.f"
			zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3137 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3142 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3144 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3148 "zgesvd.f"
			itau = 1;
#line 3149 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3155 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3155 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3157 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3163 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3163 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3165 "zgesvd.f"
			ie = 1;
#line 3166 "zgesvd.f"
			itauq = itau;
#line 3167 "zgesvd.f"
			itaup = itauq + *m;
#line 3168 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3172 "zgesvd.f"
			i__2 = *m - 1;
#line 3172 "zgesvd.f"
			i__3 = *m - 1;
#line 3172 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3179 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3179 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3188 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3188 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3191 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3198 "zgesvd.f"
			zbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 3202 "zgesvd.f"
		    }

#line 3204 "zgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3210 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3210 "zgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3214 "zgesvd.f"
			iu = 1;
#line 3215 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3219 "zgesvd.f"
			    ldwrku = *lda;
#line 3220 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3221 "zgesvd.f"
			    ldwrkr = *lda;
#line 3222 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3226 "zgesvd.f"
			    ldwrku = *lda;
#line 3227 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3228 "zgesvd.f"
			    ldwrkr = *m;
#line 3229 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3233 "zgesvd.f"
			    ldwrku = *m;
#line 3234 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3235 "zgesvd.f"
			    ldwrkr = *m;
#line 3236 "zgesvd.f"
			}
#line 3237 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3238 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3244 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3244 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3246 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3252 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3252 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3257 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3259 "zgesvd.f"
			i__2 = *m - 1;
#line 3259 "zgesvd.f"
			i__3 = *m - 1;
#line 3259 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3261 "zgesvd.f"
			ie = 1;
#line 3262 "zgesvd.f"
			itauq = itau;
#line 3263 "zgesvd.f"
			itaup = itauq + *m;
#line 3264 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 3272 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3272 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3276 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3284 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3284 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3292 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3292 "zgesvd.f"
			zungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3295 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3303 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3313 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3318 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3322 "zgesvd.f"
			zlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3325 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3329 "zgesvd.f"
			itau = 1;
#line 3330 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3336 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3336 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3338 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3344 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3344 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3346 "zgesvd.f"
			ie = 1;
#line 3347 "zgesvd.f"
			itauq = itau;
#line 3348 "zgesvd.f"
			itaup = itauq + *m;
#line 3349 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3353 "zgesvd.f"
			i__2 = *m - 1;
#line 3353 "zgesvd.f"
			i__3 = *m - 1;
#line 3353 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3360 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3360 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3369 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3369 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3377 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3377 "zgesvd.f"
			zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3379 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3387 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3391 "zgesvd.f"
		    }

#line 3393 "zgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3400 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3400 "zgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3404 "zgesvd.f"
			iu = 1;
#line 3405 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3409 "zgesvd.f"
			    ldwrku = *lda;
#line 3410 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3414 "zgesvd.f"
			    ldwrku = *m;
#line 3415 "zgesvd.f"
			}
#line 3416 "zgesvd.f"
			itau = iu + ldwrku * *m;
#line 3417 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3423 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3423 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3425 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3431 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3431 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3436 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3438 "zgesvd.f"
			i__2 = *m - 1;
#line 3438 "zgesvd.f"
			i__3 = *m - 1;
#line 3438 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3440 "zgesvd.f"
			ie = 1;
#line 3441 "zgesvd.f"
			itauq = itau;
#line 3442 "zgesvd.f"
			itaup = itauq + *m;
#line 3443 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3449 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3449 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3453 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3460 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3460 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3468 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3468 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3470 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3478 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3487 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3492 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3494 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3498 "zgesvd.f"
			itau = 1;
#line 3499 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3505 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3505 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3507 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3513 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3513 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3518 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3519 "zgesvd.f"
			i__2 = *m - 1;
#line 3519 "zgesvd.f"
			i__3 = *m - 1;
#line 3519 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3521 "zgesvd.f"
			ie = 1;
#line 3522 "zgesvd.f"
			itauq = itau;
#line 3523 "zgesvd.f"
			itaup = itauq + *m;
#line 3524 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3530 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3530 "zgesvd.f"
			zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3539 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3539 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3547 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3547 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3549 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3557 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3561 "zgesvd.f"
		    }

#line 3563 "zgesvd.f"
		}

#line 3565 "zgesvd.f"
	    }

#line 3567 "zgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3574 "zgesvd.f"
	    ie = 1;
#line 3575 "zgesvd.f"
	    itauq = 1;
#line 3576 "zgesvd.f"
	    itaup = itauq + *m;
#line 3577 "zgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 3583 "zgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3583 "zgesvd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 3586 "zgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3593 "zgesvd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3594 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3594 "zgesvd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3596 "zgesvd.f"
	    }
#line 3597 "zgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB) */
/*              (RWorkspace: 0) */

#line 3604 "zgesvd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3605 "zgesvd.f"
		if (wntva) {
#line 3605 "zgesvd.f"
		    nrvt = *n;
#line 3605 "zgesvd.f"
		}
#line 3607 "zgesvd.f"
		if (wntvs) {
#line 3607 "zgesvd.f"
		    nrvt = *m;
#line 3607 "zgesvd.f"
		}
#line 3609 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3609 "zgesvd.f"
		zungbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3611 "zgesvd.f"
	    }
#line 3612 "zgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3619 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3619 "zgesvd.f"
		zungbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3621 "zgesvd.f"
	    }
#line 3622 "zgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 3629 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3629 "zgesvd.f"
		zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3631 "zgesvd.f"
	    }
#line 3632 "zgesvd.f"
	    irwork = ie + *m;
#line 3633 "zgesvd.f"
	    if (wntuas || wntuo) {
#line 3633 "zgesvd.f"
		nru = *m;
#line 3633 "zgesvd.f"
	    }
#line 3635 "zgesvd.f"
	    if (wntun) {
#line 3635 "zgesvd.f"
		nru = 0;
#line 3635 "zgesvd.f"
	    }
#line 3637 "zgesvd.f"
	    if (wntvas || wntvo) {
#line 3637 "zgesvd.f"
		ncvt = *n;
#line 3637 "zgesvd.f"
	    }
#line 3639 "zgesvd.f"
	    if (wntvn) {
#line 3639 "zgesvd.f"
		ncvt = 0;
#line 3639 "zgesvd.f"
	    }
#line 3641 "zgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3649 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3652 "zgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3660 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3663 "zgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3671 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3674 "zgesvd.f"
	    }

#line 3676 "zgesvd.f"
	}

#line 3678 "zgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3682 "zgesvd.f"
    if (iscl == 1) {
#line 3683 "zgesvd.f"
	if (anrm > bignum) {
#line 3683 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3683 "zgesvd.f"
	}
#line 3686 "zgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3686 "zgesvd.f"
	    i__2 = minmn - 1;
#line 3686 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3686 "zgesvd.f"
	}
#line 3689 "zgesvd.f"
	if (anrm < smlnum) {
#line 3689 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3689 "zgesvd.f"
	}
#line 3692 "zgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3692 "zgesvd.f"
	    i__2 = minmn - 1;
#line 3692 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3692 "zgesvd.f"
	}
#line 3695 "zgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3699 "zgesvd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 3701 "zgesvd.f"
    return 0;

/*     End of ZGESVD */

} /* zgesvd_ */

