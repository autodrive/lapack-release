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


/*  -- LAPACK driver routine (version 3.6.0) -- */
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
		    minwrk = (*m << 1) + *n;
#line 616 "zgesvd.f"
		}
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
		i__2 = *n - 1;
#line 683 "zgesvd.f"
		i__3 = *n - 1;
#line 683 "zgesvd.f"
		zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 685 "zgesvd.f"
		ie = 1;
#line 686 "zgesvd.f"
		itauq = 1;
#line 687 "zgesvd.f"
		itaup = itauq + *n;
#line 688 "zgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 694 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 694 "zgesvd.f"
		zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 697 "zgesvd.f"
		ncvt = 0;
#line 698 "zgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 704 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 704 "zgesvd.f"
		    zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 706 "zgesvd.f"
		    ncvt = *n;
#line 707 "zgesvd.f"
		}
#line 708 "zgesvd.f"
		irwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 715 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			irwork], info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 720 "zgesvd.f"
		if (wntvas) {
#line 720 "zgesvd.f"
		    zlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 720 "zgesvd.f"
		}

#line 723 "zgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

#line 729 "zgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 733 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 734 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 734 "zgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 738 "zgesvd.f"
			ldwrku = *lda;
#line 739 "zgesvd.f"
			ldwrkr = *lda;
#line 740 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 740 "zgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 740 "zgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 744 "zgesvd.f"
			    ldwrku = *lda;
#line 745 "zgesvd.f"
			    ldwrkr = *n;
#line 746 "zgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 750 "zgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 751 "zgesvd.f"
			    ldwrkr = *n;
#line 752 "zgesvd.f"
			}
#line 752 "zgesvd.f"
		    }
#line 753 "zgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 754 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 760 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 760 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 765 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 766 "zgesvd.f"
		    i__2 = *n - 1;
#line 766 "zgesvd.f"
		    i__3 = *n - 1;
#line 766 "zgesvd.f"
		    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1], &
			    ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 773 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 773 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 775 "zgesvd.f"
		    ie = 1;
#line 776 "zgesvd.f"
		    itauq = itau;
#line 777 "zgesvd.f"
		    itaup = itauq + *n;
#line 778 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 784 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 784 "zgesvd.f"
		    zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: need 0) */

#line 792 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 792 "zgesvd.f"
		    zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 795 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 802 "zgesvd.f"
		    zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &work[ir], &ldwrkr, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 805 "zgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 812 "zgesvd.f"
		    i__2 = *m;
#line 812 "zgesvd.f"
		    i__3 = ldwrku;
#line 812 "zgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 813 "zgesvd.f"
			i__4 = *m - i__ + 1;
#line 813 "zgesvd.f"
			chunk = min(i__4,ldwrku);
#line 814 "zgesvd.f"
			zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 817 "zgesvd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 819 "zgesvd.f"
/* L10: */
#line 819 "zgesvd.f"
		    }

#line 821 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 825 "zgesvd.f"
		    ie = 1;
#line 826 "zgesvd.f"
		    itauq = 1;
#line 827 "zgesvd.f"
		    itaup = itauq + *n;
#line 828 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*                 (RWorkspace: N) */

#line 834 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 834 "zgesvd.f"
		    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 842 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 842 "zgesvd.f"
		    zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 844 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (CWorkspace: need 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 851 "zgesvd.f"
		    zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &a[a_offset], lda, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 854 "zgesvd.f"
		}

#line 856 "zgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 862 "zgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 866 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 867 "zgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 867 "zgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 871 "zgesvd.f"
			ldwrku = *lda;
#line 872 "zgesvd.f"
			ldwrkr = *lda;
#line 873 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 873 "zgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 873 "zgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 877 "zgesvd.f"
			    ldwrku = *lda;
#line 878 "zgesvd.f"
			    ldwrkr = *n;
#line 879 "zgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 883 "zgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 884 "zgesvd.f"
			    ldwrkr = *n;
#line 885 "zgesvd.f"
			}
#line 885 "zgesvd.f"
		    }
#line 886 "zgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 887 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 893 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 893 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 898 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 899 "zgesvd.f"
		    if (*n > 1) {
#line 899 "zgesvd.f"
			i__3 = *n - 1;
#line 899 "zgesvd.f"
			i__2 = *n - 1;
#line 899 "zgesvd.f"
			zlaset_("L", &i__3, &i__2, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 899 "zgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 907 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 907 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 909 "zgesvd.f"
		    ie = 1;
#line 910 "zgesvd.f"
		    itauq = itau;
#line 911 "zgesvd.f"
		    itaup = itauq + *n;
#line 912 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 918 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 918 "zgesvd.f"
		    zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 921 "zgesvd.f"
		    zlacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 927 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 927 "zgesvd.f"
		    zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 935 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 935 "zgesvd.f"
		    zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 937 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 945 "zgesvd.f"
		    zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, cdum, &c__1,
			     &rwork[irwork], info, (ftnlen)1);
#line 948 "zgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 955 "zgesvd.f"
		    i__3 = *m;
#line 955 "zgesvd.f"
		    i__2 = ldwrku;
#line 955 "zgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 956 "zgesvd.f"
			i__4 = *m - i__ + 1;
#line 956 "zgesvd.f"
			chunk = min(i__4,ldwrku);
#line 957 "zgesvd.f"
			zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 960 "zgesvd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 962 "zgesvd.f"
/* L20: */
#line 962 "zgesvd.f"
		    }

#line 964 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 968 "zgesvd.f"
		    itau = 1;
#line 969 "zgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 975 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 975 "zgesvd.f"
		    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 980 "zgesvd.f"
		    zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 981 "zgesvd.f"
		    if (*n > 1) {
#line 981 "zgesvd.f"
			i__2 = *n - 1;
#line 981 "zgesvd.f"
			i__3 = *n - 1;
#line 981 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 981 "zgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 989 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 989 "zgesvd.f"
		    zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 991 "zgesvd.f"
		    ie = 1;
#line 992 "zgesvd.f"
		    itauq = itau;
#line 993 "zgesvd.f"
		    itaup = itauq + *n;
#line 994 "zgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                 (RWorkspace: N) */

#line 1000 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1000 "zgesvd.f"
		    zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                 (RWorkspace: 0) */

#line 1008 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1008 "zgesvd.f"
		    zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 1016 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1016 "zgesvd.f"
		    zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1018 "zgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 1026 "zgesvd.f"
		    zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, 
			    &rwork[irwork], info, (ftnlen)1);

#line 1030 "zgesvd.f"
		}

#line 1032 "zgesvd.f"
	    } else if (wntus) {

#line 1034 "zgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

#line 1040 "zgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1044 "zgesvd.f"
			ir = 1;
#line 1045 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1049 "zgesvd.f"
			    ldwrkr = *lda;
#line 1050 "zgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1054 "zgesvd.f"
			    ldwrkr = *n;
#line 1055 "zgesvd.f"
			}
#line 1056 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1057 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1063 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1063 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1068 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1070 "zgesvd.f"
			i__2 = *n - 1;
#line 1070 "zgesvd.f"
			i__3 = *n - 1;
#line 1070 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1077 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1077 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1079 "zgesvd.f"
			ie = 1;
#line 1080 "zgesvd.f"
			itauq = itau;
#line 1081 "zgesvd.f"
			itaup = itauq + *n;
#line 1082 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1088 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1088 "zgesvd.f"
			zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1097 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1097 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1100 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1107 "zgesvd.f"
			zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1116 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1119 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1123 "zgesvd.f"
			itau = 1;
#line 1124 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1130 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1130 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1132 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1138 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1138 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1140 "zgesvd.f"
			ie = 1;
#line 1141 "zgesvd.f"
			itauq = itau;
#line 1142 "zgesvd.f"
			itaup = itauq + *n;
#line 1143 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1147 "zgesvd.f"
			i__2 = *n - 1;
#line 1147 "zgesvd.f"
			i__3 = *n - 1;
#line 1147 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1154 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1154 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1162 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1162 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1165 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1172 "zgesvd.f"
			zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1176 "zgesvd.f"
		    }

#line 1178 "zgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

#line 1184 "zgesvd.f"
		    if (*lwork >= (*n << 1) * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1188 "zgesvd.f"
			iu = 1;
#line 1189 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1193 "zgesvd.f"
			    ldwrku = *lda;
#line 1194 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1195 "zgesvd.f"
			    ldwrkr = *lda;
#line 1196 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1200 "zgesvd.f"
			    ldwrku = *lda;
#line 1201 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1202 "zgesvd.f"
			    ldwrkr = *n;
#line 1203 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1207 "zgesvd.f"
			    ldwrku = *n;
#line 1208 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1209 "zgesvd.f"
			    ldwrkr = *n;
#line 1210 "zgesvd.f"
			}
#line 1211 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1212 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1218 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1218 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1223 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1225 "zgesvd.f"
			i__2 = *n - 1;
#line 1225 "zgesvd.f"
			i__3 = *n - 1;
#line 1225 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1232 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1232 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1234 "zgesvd.f"
			ie = 1;
#line 1235 "zgesvd.f"
			itauq = itau;
#line 1236 "zgesvd.f"
			itaup = itauq + *n;
#line 1237 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1245 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1245 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1249 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1256 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1256 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1265 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1265 "zgesvd.f"
			zungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1268 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1276 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1286 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1293 "zgesvd.f"
			zlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1296 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1300 "zgesvd.f"
			itau = 1;
#line 1301 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1307 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1307 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1309 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1315 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1315 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1317 "zgesvd.f"
			ie = 1;
#line 1318 "zgesvd.f"
			itauq = itau;
#line 1319 "zgesvd.f"
			itaup = itauq + *n;
#line 1320 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1324 "zgesvd.f"
			i__2 = *n - 1;
#line 1324 "zgesvd.f"
			i__3 = *n - 1;
#line 1324 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1331 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1331 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1339 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1339 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1347 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1347 "zgesvd.f"
			zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1349 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1357 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1361 "zgesvd.f"
		    }

#line 1363 "zgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

#line 1370 "zgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1374 "zgesvd.f"
			iu = 1;
#line 1375 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1379 "zgesvd.f"
			    ldwrku = *lda;
#line 1380 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1384 "zgesvd.f"
			    ldwrku = *n;
#line 1385 "zgesvd.f"
			}
#line 1386 "zgesvd.f"
			itau = iu + ldwrku * *n;
#line 1387 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1393 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1393 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1398 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1400 "zgesvd.f"
			i__2 = *n - 1;
#line 1400 "zgesvd.f"
			i__3 = *n - 1;
#line 1400 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1407 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1407 "zgesvd.f"
			zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1409 "zgesvd.f"
			ie = 1;
#line 1410 "zgesvd.f"
			itauq = itau;
#line 1411 "zgesvd.f"
			itaup = itauq + *n;
#line 1412 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1418 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1418 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1422 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1429 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1429 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1438 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1438 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1440 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1448 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1457 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1460 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1464 "zgesvd.f"
			itau = 1;
#line 1465 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1471 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1471 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1473 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1479 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1479 "zgesvd.f"
			zungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1484 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1485 "zgesvd.f"
			if (*n > 1) {
#line 1485 "zgesvd.f"
			    i__2 = *n - 1;
#line 1485 "zgesvd.f"
			    i__3 = *n - 1;
#line 1485 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1485 "zgesvd.f"
			}
#line 1488 "zgesvd.f"
			ie = 1;
#line 1489 "zgesvd.f"
			itauq = itau;
#line 1490 "zgesvd.f"
			itaup = itauq + *n;
#line 1491 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1497 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1497 "zgesvd.f"
			zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1506 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1506 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1514 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1514 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1516 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1524 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1528 "zgesvd.f"
		    }

#line 1530 "zgesvd.f"
		}

#line 1532 "zgesvd.f"
	    } else if (wntua) {

#line 1534 "zgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1540 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1540 "zgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1544 "zgesvd.f"
			ir = 1;
#line 1545 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1549 "zgesvd.f"
			    ldwrkr = *lda;
#line 1550 "zgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1554 "zgesvd.f"
			    ldwrkr = *n;
#line 1555 "zgesvd.f"
			}
#line 1556 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1557 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1563 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1563 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1565 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1569 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1571 "zgesvd.f"
			i__2 = *n - 1;
#line 1571 "zgesvd.f"
			i__3 = *n - 1;
#line 1571 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1578 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1578 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1580 "zgesvd.f"
			ie = 1;
#line 1581 "zgesvd.f"
			itauq = itau;
#line 1582 "zgesvd.f"
			itaup = itauq + *n;
#line 1583 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1589 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1589 "zgesvd.f"
			zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1598 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1598 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1601 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1608 "zgesvd.f"
			zbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1617 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1622 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1624 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1628 "zgesvd.f"
			itau = 1;
#line 1629 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1635 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1635 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1637 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1643 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1643 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1645 "zgesvd.f"
			ie = 1;
#line 1646 "zgesvd.f"
			itauq = itau;
#line 1647 "zgesvd.f"
			itaup = itauq + *n;
#line 1648 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1652 "zgesvd.f"
			i__2 = *n - 1;
#line 1652 "zgesvd.f"
			i__3 = *n - 1;
#line 1652 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1659 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1659 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1668 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1668 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1671 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1678 "zgesvd.f"
			zbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1682 "zgesvd.f"
		    }

#line 1684 "zgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1690 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1690 "zgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1694 "zgesvd.f"
			iu = 1;
#line 1695 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1699 "zgesvd.f"
			    ldwrku = *lda;
#line 1700 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1701 "zgesvd.f"
			    ldwrkr = *lda;
#line 1702 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1706 "zgesvd.f"
			    ldwrku = *lda;
#line 1707 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1708 "zgesvd.f"
			    ldwrkr = *n;
#line 1709 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1713 "zgesvd.f"
			    ldwrku = *n;
#line 1714 "zgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1715 "zgesvd.f"
			    ldwrkr = *n;
#line 1716 "zgesvd.f"
			}
#line 1717 "zgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1718 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1724 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1724 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1726 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1732 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1732 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1737 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1739 "zgesvd.f"
			i__2 = *n - 1;
#line 1739 "zgesvd.f"
			i__3 = *n - 1;
#line 1739 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1741 "zgesvd.f"
			ie = 1;
#line 1742 "zgesvd.f"
			itauq = itau;
#line 1743 "zgesvd.f"
			itaup = itauq + *n;
#line 1744 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1752 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1752 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1756 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1763 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1763 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1772 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1772 "zgesvd.f"
			zungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1775 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1783 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1793 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1798 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1802 "zgesvd.f"
			zlacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1805 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1809 "zgesvd.f"
			itau = 1;
#line 1810 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1816 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1816 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1818 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1824 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1824 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1826 "zgesvd.f"
			ie = 1;
#line 1827 "zgesvd.f"
			itauq = itau;
#line 1828 "zgesvd.f"
			itaup = itauq + *n;
#line 1829 "zgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1833 "zgesvd.f"
			i__2 = *n - 1;
#line 1833 "zgesvd.f"
			i__3 = *n - 1;
#line 1833 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1840 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1840 "zgesvd.f"
			zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1849 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1849 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1857 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1857 "zgesvd.f"
			zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1859 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1867 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1871 "zgesvd.f"
		    }

#line 1873 "zgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1880 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1880 "zgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1884 "zgesvd.f"
			iu = 1;
#line 1885 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1889 "zgesvd.f"
			    ldwrku = *lda;
#line 1890 "zgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1894 "zgesvd.f"
			    ldwrku = *n;
#line 1895 "zgesvd.f"
			}
#line 1896 "zgesvd.f"
			itau = iu + ldwrku * *n;
#line 1897 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1903 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1903 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1905 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1911 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1911 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1916 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1918 "zgesvd.f"
			i__2 = *n - 1;
#line 1918 "zgesvd.f"
			i__3 = *n - 1;
#line 1918 "zgesvd.f"
			zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1920 "zgesvd.f"
			ie = 1;
#line 1921 "zgesvd.f"
			itauq = itau;
#line 1922 "zgesvd.f"
			itaup = itauq + *n;
#line 1923 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1929 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1929 "zgesvd.f"
			zgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1933 "zgesvd.f"
			zlacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1940 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1940 "zgesvd.f"
			zungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: need   0) */

#line 1949 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1949 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1951 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1959 "zgesvd.f"
			zbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1968 "zgesvd.f"
			zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1973 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1975 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1979 "zgesvd.f"
			itau = 1;
#line 1980 "zgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1986 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1986 "zgesvd.f"
			zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1988 "zgesvd.f"
			zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1994 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1994 "zgesvd.f"
			zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 1999 "zgesvd.f"
			zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 2000 "zgesvd.f"
			if (*n > 1) {
#line 2000 "zgesvd.f"
			    i__2 = *n - 1;
#line 2000 "zgesvd.f"
			    i__3 = *n - 1;
#line 2000 "zgesvd.f"
			    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 2000 "zgesvd.f"
			}
#line 2003 "zgesvd.f"
			ie = 1;
#line 2004 "zgesvd.f"
			itauq = itau;
#line 2005 "zgesvd.f"
			itaup = itauq + *n;
#line 2006 "zgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 2012 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2012 "zgesvd.f"
			zgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2021 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2021 "zgesvd.f"
			zunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2029 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2029 "zgesvd.f"
			zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 2031 "zgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2039 "zgesvd.f"
			zbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2043 "zgesvd.f"
		    }

#line 2045 "zgesvd.f"
		}

#line 2047 "zgesvd.f"
	    }

#line 2049 "zgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 2056 "zgesvd.f"
	    ie = 1;
#line 2057 "zgesvd.f"
	    itauq = 1;
#line 2058 "zgesvd.f"
	    itaup = itauq + *n;
#line 2059 "zgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 2065 "zgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 2065 "zgesvd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 2068 "zgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB) */
/*              (RWorkspace: 0) */

#line 2075 "zgesvd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 2076 "zgesvd.f"
		if (wntus) {
#line 2076 "zgesvd.f"
		    ncu = *n;
#line 2076 "zgesvd.f"
		}
#line 2078 "zgesvd.f"
		if (wntua) {
#line 2078 "zgesvd.f"
		    ncu = *m;
#line 2078 "zgesvd.f"
		}
#line 2080 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2080 "zgesvd.f"
		zungbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2082 "zgesvd.f"
	    }
#line 2083 "zgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2090 "zgesvd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2091 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2091 "zgesvd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2093 "zgesvd.f"
	    }
#line 2094 "zgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 2101 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2101 "zgesvd.f"
		zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2103 "zgesvd.f"
	    }
#line 2104 "zgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2111 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2111 "zgesvd.f"
		zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2113 "zgesvd.f"
	    }
#line 2114 "zgesvd.f"
	    irwork = ie + *n;
#line 2115 "zgesvd.f"
	    if (wntuas || wntuo) {
#line 2115 "zgesvd.f"
		nru = *m;
#line 2115 "zgesvd.f"
	    }
#line 2117 "zgesvd.f"
	    if (wntun) {
#line 2117 "zgesvd.f"
		nru = 0;
#line 2117 "zgesvd.f"
	    }
#line 2119 "zgesvd.f"
	    if (wntvas || wntvo) {
#line 2119 "zgesvd.f"
		ncvt = *n;
#line 2119 "zgesvd.f"
	    }
#line 2121 "zgesvd.f"
	    if (wntvn) {
#line 2121 "zgesvd.f"
		ncvt = 0;
#line 2121 "zgesvd.f"
	    }
#line 2123 "zgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2131 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2134 "zgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2142 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2145 "zgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2153 "zgesvd.f"
		zbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2156 "zgesvd.f"
	    }

#line 2158 "zgesvd.f"
	}

#line 2160 "zgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2166 "zgesvd.f"
	if (*n >= mnthr) {

#line 2168 "zgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2173 "zgesvd.f"
		itau = 1;
#line 2174 "zgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 2180 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2180 "zgesvd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2185 "zgesvd.f"
		i__2 = *m - 1;
#line 2185 "zgesvd.f"
		i__3 = *m - 1;
#line 2185 "zgesvd.f"
		zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 2187 "zgesvd.f"
		ie = 1;
#line 2188 "zgesvd.f"
		itauq = 1;
#line 2189 "zgesvd.f"
		itaup = itauq + *m;
#line 2190 "zgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 2196 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2196 "zgesvd.f"
		zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2199 "zgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2205 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2205 "zgesvd.f"
		    zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2207 "zgesvd.f"
		}
#line 2208 "zgesvd.f"
		irwork = ie + *m;
#line 2209 "zgesvd.f"
		nru = 0;
#line 2210 "zgesvd.f"
		if (wntuo || wntuas) {
#line 2210 "zgesvd.f"
		    nru = *m;
#line 2210 "zgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2218 "zgesvd.f"
		zbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &rwork[ie], cdum, &
			c__1, &a[a_offset], lda, cdum, &c__1, &rwork[irwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2223 "zgesvd.f"
		if (wntuas) {
#line 2223 "zgesvd.f"
		    zlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2223 "zgesvd.f"
		}

#line 2226 "zgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

#line 2232 "zgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2236 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2237 "zgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 2237 "zgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2241 "zgesvd.f"
			ldwrku = *lda;
#line 2242 "zgesvd.f"
			chunk = *n;
#line 2243 "zgesvd.f"
			ldwrkr = *lda;
#line 2244 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2244 "zgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 2244 "zgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2248 "zgesvd.f"
			    ldwrku = *lda;
#line 2249 "zgesvd.f"
			    chunk = *n;
#line 2250 "zgesvd.f"
			    ldwrkr = *m;
#line 2251 "zgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2255 "zgesvd.f"
			    ldwrku = *m;
#line 2256 "zgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2257 "zgesvd.f"
			    ldwrkr = *m;
#line 2258 "zgesvd.f"
			}
#line 2258 "zgesvd.f"
		    }
#line 2259 "zgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2260 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2266 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2266 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2271 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2272 "zgesvd.f"
		    i__2 = *m - 1;
#line 2272 "zgesvd.f"
		    i__3 = *m - 1;
#line 2272 "zgesvd.f"
		    zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2279 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2279 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2281 "zgesvd.f"
		    ie = 1;
#line 2282 "zgesvd.f"
		    itauq = itau;
#line 2283 "zgesvd.f"
		    itaup = itauq + *m;
#line 2284 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2290 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2290 "zgesvd.f"
		    zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2298 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2298 "zgesvd.f"
		    zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2301 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2308 "zgesvd.f"
		    zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &work[
			    ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2311 "zgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N) */
/*                 (RWorkspace: 0) */

#line 2318 "zgesvd.f"
		    i__2 = *n;
#line 2318 "zgesvd.f"
		    i__3 = chunk;
#line 2318 "zgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2319 "zgesvd.f"
			i__4 = *n - i__ + 1;
#line 2319 "zgesvd.f"
			blk = min(i__4,chunk);
#line 2320 "zgesvd.f"
			zgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2323 "zgesvd.f"
			zlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2325 "zgesvd.f"
/* L30: */
#line 2325 "zgesvd.f"
		    }

#line 2327 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2331 "zgesvd.f"
		    ie = 1;
#line 2332 "zgesvd.f"
		    itauq = 1;
#line 2333 "zgesvd.f"
		    itaup = itauq + *m;
#line 2334 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*                 (RWorkspace: need M) */

#line 2340 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2340 "zgesvd.f"
		    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2348 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2348 "zgesvd.f"
		    zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2350 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2357 "zgesvd.f"
		    zbdsqr_("L", m, n, &c__0, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 2360 "zgesvd.f"
		}

#line 2362 "zgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 2368 "zgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2372 "zgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2373 "zgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 2373 "zgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2377 "zgesvd.f"
			ldwrku = *lda;
#line 2378 "zgesvd.f"
			chunk = *n;
#line 2379 "zgesvd.f"
			ldwrkr = *lda;
#line 2380 "zgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2380 "zgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 2380 "zgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2384 "zgesvd.f"
			    ldwrku = *lda;
#line 2385 "zgesvd.f"
			    chunk = *n;
#line 2386 "zgesvd.f"
			    ldwrkr = *m;
#line 2387 "zgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2391 "zgesvd.f"
			    ldwrku = *m;
#line 2392 "zgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2393 "zgesvd.f"
			    ldwrkr = *m;
#line 2394 "zgesvd.f"
			}
#line 2394 "zgesvd.f"
		    }
#line 2395 "zgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2396 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2402 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2402 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2407 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2408 "zgesvd.f"
		    i__3 = *m - 1;
#line 2408 "zgesvd.f"
		    i__2 = *m - 1;
#line 2408 "zgesvd.f"
		    zlaset_("U", &i__3, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2415 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2415 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2417 "zgesvd.f"
		    ie = 1;
#line 2418 "zgesvd.f"
		    itauq = itau;
#line 2419 "zgesvd.f"
		    itaup = itauq + *m;
#line 2420 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2426 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2426 "zgesvd.f"
		    zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2429 "zgesvd.f"
		    zlacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2435 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2435 "zgesvd.f"
		    zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2443 "zgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2443 "zgesvd.f"
		    zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2445 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2453 "zgesvd.f"
		    zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[ir],
			     &ldwrkr, &u[u_offset], ldu, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2456 "zgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N)) */
/*                 (RWorkspace: 0) */

#line 2463 "zgesvd.f"
		    i__3 = *n;
#line 2463 "zgesvd.f"
		    i__2 = chunk;
#line 2463 "zgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2464 "zgesvd.f"
			i__4 = *n - i__ + 1;
#line 2464 "zgesvd.f"
			blk = min(i__4,chunk);
#line 2465 "zgesvd.f"
			zgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2468 "zgesvd.f"
			zlacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2470 "zgesvd.f"
/* L40: */
#line 2470 "zgesvd.f"
		    }

#line 2472 "zgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2476 "zgesvd.f"
		    itau = 1;
#line 2477 "zgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2483 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2483 "zgesvd.f"
		    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2488 "zgesvd.f"
		    zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2489 "zgesvd.f"
		    i__2 = *m - 1;
#line 2489 "zgesvd.f"
		    i__3 = *m - 1;
#line 2489 "zgesvd.f"
		    zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2496 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2496 "zgesvd.f"
		    zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2498 "zgesvd.f"
		    ie = 1;
#line 2499 "zgesvd.f"
		    itauq = itau;
#line 2500 "zgesvd.f"
		    itaup = itauq + *m;
#line 2501 "zgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2507 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2507 "zgesvd.f"
		    zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                 (RWorkspace: 0) */

#line 2515 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2515 "zgesvd.f"
		    zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2523 "zgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2523 "zgesvd.f"
		    zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2525 "zgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2533 "zgesvd.f"
		    zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			    rwork[irwork], info, (ftnlen)1);

#line 2536 "zgesvd.f"
		}

#line 2538 "zgesvd.f"
	    } else if (wntvs) {

#line 2540 "zgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

#line 2546 "zgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2550 "zgesvd.f"
			ir = 1;
#line 2551 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2555 "zgesvd.f"
			    ldwrkr = *lda;
#line 2556 "zgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2560 "zgesvd.f"
			    ldwrkr = *m;
#line 2561 "zgesvd.f"
			}
#line 2562 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2563 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2569 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2569 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2574 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2576 "zgesvd.f"
			i__2 = *m - 1;
#line 2576 "zgesvd.f"
			i__3 = *m - 1;
#line 2576 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2583 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2583 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2585 "zgesvd.f"
			ie = 1;
#line 2586 "zgesvd.f"
			itauq = itau;
#line 2587 "zgesvd.f"
			itaup = itauq + *m;
#line 2588 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2594 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2594 "zgesvd.f"
			zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2604 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2604 "zgesvd.f"
			zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2607 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2614 "zgesvd.f"
			zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2623 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2626 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2630 "zgesvd.f"
			itau = 1;
#line 2631 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2637 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2637 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2642 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2648 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2648 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2650 "zgesvd.f"
			ie = 1;
#line 2651 "zgesvd.f"
			itauq = itau;
#line 2652 "zgesvd.f"
			itaup = itauq + *m;
#line 2653 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2657 "zgesvd.f"
			i__2 = *m - 1;
#line 2657 "zgesvd.f"
			i__3 = *m - 1;
#line 2657 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2664 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2664 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2672 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2672 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2675 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2682 "zgesvd.f"
			zbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 2686 "zgesvd.f"
		    }

#line 2688 "zgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

#line 2694 "zgesvd.f"
		    if (*lwork >= (*m << 1) * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2698 "zgesvd.f"
			iu = 1;
#line 2699 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2703 "zgesvd.f"
			    ldwrku = *lda;
#line 2704 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2705 "zgesvd.f"
			    ldwrkr = *lda;
#line 2706 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2710 "zgesvd.f"
			    ldwrku = *lda;
#line 2711 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2712 "zgesvd.f"
			    ldwrkr = *m;
#line 2713 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2717 "zgesvd.f"
			    ldwrku = *m;
#line 2718 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2719 "zgesvd.f"
			    ldwrkr = *m;
#line 2720 "zgesvd.f"
			}
#line 2721 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2722 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2728 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2728 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2733 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2735 "zgesvd.f"
			i__2 = *m - 1;
#line 2735 "zgesvd.f"
			i__3 = *m - 1;
#line 2735 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2742 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2742 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2744 "zgesvd.f"
			ie = 1;
#line 2745 "zgesvd.f"
			itauq = itau;
#line 2746 "zgesvd.f"
			itaup = itauq + *m;
#line 2747 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 2755 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2755 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2759 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2767 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2767 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2775 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2775 "zgesvd.f"
			zungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2778 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2786 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2796 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2803 "zgesvd.f"
			zlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2806 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2810 "zgesvd.f"
			itau = 1;
#line 2811 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2817 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2817 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2819 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2825 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2825 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2827 "zgesvd.f"
			ie = 1;
#line 2828 "zgesvd.f"
			itauq = itau;
#line 2829 "zgesvd.f"
			itaup = itauq + *m;
#line 2830 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2834 "zgesvd.f"
			i__2 = *m - 1;
#line 2834 "zgesvd.f"
			i__3 = *m - 1;
#line 2834 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2841 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2841 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2849 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2849 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2857 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2857 "zgesvd.f"
			zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2859 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2867 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2871 "zgesvd.f"
		    }

#line 2873 "zgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

#line 2880 "zgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2884 "zgesvd.f"
			iu = 1;
#line 2885 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2889 "zgesvd.f"
			    ldwrku = *lda;
#line 2890 "zgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2894 "zgesvd.f"
			    ldwrku = *m;
#line 2895 "zgesvd.f"
			}
#line 2896 "zgesvd.f"
			itau = iu + ldwrku * *m;
#line 2897 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2903 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2903 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2908 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2910 "zgesvd.f"
			i__2 = *m - 1;
#line 2910 "zgesvd.f"
			i__3 = *m - 1;
#line 2910 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2917 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2917 "zgesvd.f"
			zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2919 "zgesvd.f"
			ie = 1;
#line 2920 "zgesvd.f"
			itauq = itau;
#line 2921 "zgesvd.f"
			itaup = itauq + *m;
#line 2922 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2928 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2928 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2932 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2940 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2940 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2948 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2948 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2950 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2958 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2967 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2970 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2974 "zgesvd.f"
			itau = 1;
#line 2975 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2981 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2981 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2983 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2989 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2989 "zgesvd.f"
			zunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2994 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2995 "zgesvd.f"
			i__2 = *m - 1;
#line 2995 "zgesvd.f"
			i__3 = *m - 1;
#line 2995 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 2997 "zgesvd.f"
			ie = 1;
#line 2998 "zgesvd.f"
			itauq = itau;
#line 2999 "zgesvd.f"
			itaup = itauq + *m;
#line 3000 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3006 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3006 "zgesvd.f"
			zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3015 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3015 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3023 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3023 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3025 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3033 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3037 "zgesvd.f"
		    }

#line 3039 "zgesvd.f"
		}

#line 3041 "zgesvd.f"
	    } else if (wntva) {

#line 3043 "zgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 3049 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3049 "zgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3053 "zgesvd.f"
			ir = 1;
#line 3054 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 3058 "zgesvd.f"
			    ldwrkr = *lda;
#line 3059 "zgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 3063 "zgesvd.f"
			    ldwrkr = *m;
#line 3064 "zgesvd.f"
			}
#line 3065 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3066 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3072 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3072 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3074 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 3078 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 3080 "zgesvd.f"
			i__2 = *m - 1;
#line 3080 "zgesvd.f"
			i__3 = *m - 1;
#line 3080 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3087 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3087 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3089 "zgesvd.f"
			ie = 1;
#line 3090 "zgesvd.f"
			itauq = itau;
#line 3091 "zgesvd.f"
			itaup = itauq + *m;
#line 3092 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3098 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3098 "zgesvd.f"
			zgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3108 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3108 "zgesvd.f"
			zungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3111 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3118 "zgesvd.f"
			zbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3127 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3132 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3134 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3138 "zgesvd.f"
			itau = 1;
#line 3139 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3145 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3145 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3147 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3153 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3153 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3155 "zgesvd.f"
			ie = 1;
#line 3156 "zgesvd.f"
			itauq = itau;
#line 3157 "zgesvd.f"
			itaup = itauq + *m;
#line 3158 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3162 "zgesvd.f"
			i__2 = *m - 1;
#line 3162 "zgesvd.f"
			i__3 = *m - 1;
#line 3162 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3169 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3169 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3178 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3178 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3181 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3188 "zgesvd.f"
			zbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 3192 "zgesvd.f"
		    }

#line 3194 "zgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3200 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3200 "zgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3204 "zgesvd.f"
			iu = 1;
#line 3205 "zgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3209 "zgesvd.f"
			    ldwrku = *lda;
#line 3210 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3211 "zgesvd.f"
			    ldwrkr = *lda;
#line 3212 "zgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3216 "zgesvd.f"
			    ldwrku = *lda;
#line 3217 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3218 "zgesvd.f"
			    ldwrkr = *m;
#line 3219 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3223 "zgesvd.f"
			    ldwrku = *m;
#line 3224 "zgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3225 "zgesvd.f"
			    ldwrkr = *m;
#line 3226 "zgesvd.f"
			}
#line 3227 "zgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3228 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3234 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3234 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3236 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3242 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3242 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3247 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3249 "zgesvd.f"
			i__2 = *m - 1;
#line 3249 "zgesvd.f"
			i__3 = *m - 1;
#line 3249 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3251 "zgesvd.f"
			ie = 1;
#line 3252 "zgesvd.f"
			itauq = itau;
#line 3253 "zgesvd.f"
			itaup = itauq + *m;
#line 3254 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 3262 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3262 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3266 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3274 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3274 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3282 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3282 "zgesvd.f"
			zungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3285 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3293 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3303 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3308 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3312 "zgesvd.f"
			zlacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3315 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3319 "zgesvd.f"
			itau = 1;
#line 3320 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3326 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3326 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3328 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3334 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3334 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3336 "zgesvd.f"
			ie = 1;
#line 3337 "zgesvd.f"
			itauq = itau;
#line 3338 "zgesvd.f"
			itaup = itauq + *m;
#line 3339 "zgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3343 "zgesvd.f"
			i__2 = *m - 1;
#line 3343 "zgesvd.f"
			i__3 = *m - 1;
#line 3343 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3350 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3350 "zgesvd.f"
			zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3359 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3359 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3367 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3367 "zgesvd.f"
			zungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3369 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3377 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3381 "zgesvd.f"
		    }

#line 3383 "zgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3390 "zgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3390 "zgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3394 "zgesvd.f"
			iu = 1;
#line 3395 "zgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3399 "zgesvd.f"
			    ldwrku = *lda;
#line 3400 "zgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3404 "zgesvd.f"
			    ldwrku = *m;
#line 3405 "zgesvd.f"
			}
#line 3406 "zgesvd.f"
			itau = iu + ldwrku * *m;
#line 3407 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3413 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3413 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3415 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3421 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3421 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3426 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3428 "zgesvd.f"
			i__2 = *m - 1;
#line 3428 "zgesvd.f"
			i__3 = *m - 1;
#line 3428 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3430 "zgesvd.f"
			ie = 1;
#line 3431 "zgesvd.f"
			itauq = itau;
#line 3432 "zgesvd.f"
			itaup = itauq + *m;
#line 3433 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3439 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3439 "zgesvd.f"
			zgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3443 "zgesvd.f"
			zlacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3450 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3450 "zgesvd.f"
			zungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3458 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3458 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3460 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3468 "zgesvd.f"
			zbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3477 "zgesvd.f"
			zgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3482 "zgesvd.f"
			zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3484 "zgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3488 "zgesvd.f"
			itau = 1;
#line 3489 "zgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3495 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3495 "zgesvd.f"
			zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3497 "zgesvd.f"
			zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3503 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3503 "zgesvd.f"
			zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3508 "zgesvd.f"
			zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3509 "zgesvd.f"
			i__2 = *m - 1;
#line 3509 "zgesvd.f"
			i__3 = *m - 1;
#line 3509 "zgesvd.f"
			zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3511 "zgesvd.f"
			ie = 1;
#line 3512 "zgesvd.f"
			itauq = itau;
#line 3513 "zgesvd.f"
			itaup = itauq + *m;
#line 3514 "zgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3520 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3520 "zgesvd.f"
			zgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3529 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3529 "zgesvd.f"
			zunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3537 "zgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3537 "zgesvd.f"
			zungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3539 "zgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3547 "zgesvd.f"
			zbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3551 "zgesvd.f"
		    }

#line 3553 "zgesvd.f"
		}

#line 3555 "zgesvd.f"
	    }

#line 3557 "zgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3564 "zgesvd.f"
	    ie = 1;
#line 3565 "zgesvd.f"
	    itauq = 1;
#line 3566 "zgesvd.f"
	    itaup = itauq + *m;
#line 3567 "zgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 3573 "zgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3573 "zgesvd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 3576 "zgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3583 "zgesvd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3584 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3584 "zgesvd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3586 "zgesvd.f"
	    }
#line 3587 "zgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB) */
/*              (RWorkspace: 0) */

#line 3594 "zgesvd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3595 "zgesvd.f"
		if (wntva) {
#line 3595 "zgesvd.f"
		    nrvt = *n;
#line 3595 "zgesvd.f"
		}
#line 3597 "zgesvd.f"
		if (wntvs) {
#line 3597 "zgesvd.f"
		    nrvt = *m;
#line 3597 "zgesvd.f"
		}
#line 3599 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3599 "zgesvd.f"
		zungbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3601 "zgesvd.f"
	    }
#line 3602 "zgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3609 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3609 "zgesvd.f"
		zungbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3611 "zgesvd.f"
	    }
#line 3612 "zgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 3619 "zgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3619 "zgesvd.f"
		zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3621 "zgesvd.f"
	    }
#line 3622 "zgesvd.f"
	    irwork = ie + *m;
#line 3623 "zgesvd.f"
	    if (wntuas || wntuo) {
#line 3623 "zgesvd.f"
		nru = *m;
#line 3623 "zgesvd.f"
	    }
#line 3625 "zgesvd.f"
	    if (wntun) {
#line 3625 "zgesvd.f"
		nru = 0;
#line 3625 "zgesvd.f"
	    }
#line 3627 "zgesvd.f"
	    if (wntvas || wntvo) {
#line 3627 "zgesvd.f"
		ncvt = *n;
#line 3627 "zgesvd.f"
	    }
#line 3629 "zgesvd.f"
	    if (wntvn) {
#line 3629 "zgesvd.f"
		ncvt = 0;
#line 3629 "zgesvd.f"
	    }
#line 3631 "zgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3639 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3642 "zgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3650 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3653 "zgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3661 "zgesvd.f"
		zbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3664 "zgesvd.f"
	    }

#line 3666 "zgesvd.f"
	}

#line 3668 "zgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3672 "zgesvd.f"
    if (iscl == 1) {
#line 3673 "zgesvd.f"
	if (anrm > bignum) {
#line 3673 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3673 "zgesvd.f"
	}
#line 3676 "zgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3676 "zgesvd.f"
	    i__2 = minmn - 1;
#line 3676 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3676 "zgesvd.f"
	}
#line 3679 "zgesvd.f"
	if (anrm < smlnum) {
#line 3679 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3679 "zgesvd.f"
	}
#line 3682 "zgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3682 "zgesvd.f"
	    i__2 = minmn - 1;
#line 3682 "zgesvd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3682 "zgesvd.f"
	}
#line 3685 "zgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3689 "zgesvd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 3691 "zgesvd.f"
    return 0;

/*     End of ZGESVD */

} /* zgesvd_ */

