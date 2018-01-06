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
    extern /* Subroutine */ int cgebrd_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), cgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), cbdsqr_(
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), cungbr_(char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int cunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen), cunglq_(integer *, integer *, integer *
	    , doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer ldwrkr, minwrk, ldwrku, maxwrk;
    static doublereal smlnum;
    static integer irwork;
    static logical lquery, wntuas, wntvas;
    static integer lwork_cungbr_p__, lwork_cungbr_q__, lwork_cunglq_m__, 
	    lwork_cunglq_n__, lwork_cungqr_m__, lwork_cungqr_n__;


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
	    cgeqrf_(m, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 325 "cgesvd.f"
	    lwork_cgeqrf__ = (integer) cdum[0].r;
/*           Compute space needed for CUNGQR */
#line 327 "cgesvd.f"
	    cungqr_(m, n, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 328 "cgesvd.f"
	    lwork_cungqr_n__ = (integer) cdum[0].r;
#line 329 "cgesvd.f"
	    cungqr_(m, m, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 330 "cgesvd.f"
	    lwork_cungqr_m__ = (integer) cdum[0].r;
/*           Compute space needed for CGEBRD */
#line 332 "cgesvd.f"
	    cgebrd_(n, n, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum, &
		    c_n1, &ierr);
#line 334 "cgesvd.f"
	    lwork_cgebrd__ = (integer) cdum[0].r;
/*           Compute space needed for CUNGBR */
#line 336 "cgesvd.f"
	    cungbr_("P", n, n, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr,
		     (ftnlen)1);
#line 338 "cgesvd.f"
	    lwork_cungbr_p__ = (integer) cdum[0].r;
#line 339 "cgesvd.f"
	    cungbr_("Q", n, n, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr,
		     (ftnlen)1);
#line 341 "cgesvd.f"
	    lwork_cungbr_q__ = (integer) cdum[0].r;

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
		cgebrd_(m, n, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum,
			 &c_n1, &ierr);
#line 449 "cgesvd.f"
		lwork_cgebrd__ = (integer) cdum[0].r;
#line 450 "cgesvd.f"
		maxwrk = (*n << 1) + lwork_cgebrd__;
#line 451 "cgesvd.f"
		if (wntus || wntuo) {
#line 452 "cgesvd.f"
		    cungbr_("Q", m, n, n, &a[a_offset], lda, cdum, cdum, &
			    c_n1, &ierr, (ftnlen)1);
#line 454 "cgesvd.f"
		    lwork_cungbr_q__ = (integer) cdum[0].r;
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
		    cungbr_("Q", m, m, n, &a[a_offset], lda, cdum, cdum, &
			    c_n1, &ierr, (ftnlen)1);
#line 460 "cgesvd.f"
		    lwork_cungbr_q__ = (integer) cdum[0].r;
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
		}
#line 466 "cgesvd.f"
		minwrk = (*n << 1) + *m;
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
	    cgelqf_(m, n, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 475 "cgesvd.f"
	    lwork_cgelqf__ = (integer) cdum[0].r;
/*           Compute space needed for CUNGLQ */
#line 477 "cgesvd.f"
	    cunglq_(n, n, m, cdum, n, cdum, cdum, &c_n1, &ierr);
#line 479 "cgesvd.f"
	    lwork_cunglq_n__ = (integer) cdum[0].r;
#line 480 "cgesvd.f"
	    cunglq_(m, n, m, &a[a_offset], lda, cdum, cdum, &c_n1, &ierr);
#line 481 "cgesvd.f"
	    lwork_cunglq_m__ = (integer) cdum[0].r;
/*           Compute space needed for CGEBRD */
#line 483 "cgesvd.f"
	    cgebrd_(m, m, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum, &
		    c_n1, &ierr);
#line 485 "cgesvd.f"
	    lwork_cgebrd__ = (integer) cdum[0].r;
/*            Compute space needed for CUNGBR P */
#line 487 "cgesvd.f"
	    cungbr_("P", m, m, m, &a[a_offset], n, cdum, cdum, &c_n1, &ierr, (
		    ftnlen)1);
#line 489 "cgesvd.f"
	    lwork_cungbr_p__ = (integer) cdum[0].r;
/*           Compute space needed for CUNGBR Q */
#line 491 "cgesvd.f"
	    cungbr_("Q", m, m, m, &a[a_offset], n, cdum, cdum, &c_n1, &ierr, (
		    ftnlen)1);
#line 493 "cgesvd.f"
	    lwork_cungbr_q__ = (integer) cdum[0].r;
#line 494 "cgesvd.f"
	    if (*n >= mnthr) {
#line 495 "cgesvd.f"
		if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N') */

#line 499 "cgesvd.f"
		    maxwrk = *m + lwork_cgelqf__;
/* Computing MAX */
#line 500 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cgebrd__;
#line 500 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 501 "cgesvd.f"
		    if (wntuo || wntuas) {
/* Computing MAX */
#line 501 "cgesvd.f"
			i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 501 "cgesvd.f"
			maxwrk = max(i__2,i__3);
#line 501 "cgesvd.f"
		    }
#line 503 "cgesvd.f"
		    minwrk = *m * 3;
#line 504 "cgesvd.f"
		} else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O') */

#line 508 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 509 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 509 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 510 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 510 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 511 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 511 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 512 "cgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 512 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 513 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 514 "cgesvd.f"
		} else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='O') */

#line 519 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 520 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 520 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 521 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 521 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 522 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 522 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 523 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 523 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 524 "cgesvd.f"
		    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n;
#line 524 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 525 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 526 "cgesvd.f"
		} else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */

#line 530 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 531 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 531 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 532 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 532 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 533 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 533 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 534 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 535 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 536 "cgesvd.f"
		} else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */

#line 540 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 541 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 541 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 542 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 542 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 543 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 543 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 544 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 544 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 545 "cgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 546 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 547 "cgesvd.f"
		} else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='S') */

#line 552 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 553 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_m__;
#line 553 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 554 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 554 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 555 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 555 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 556 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 556 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 557 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 558 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 559 "cgesvd.f"
		} else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */

#line 563 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 564 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_n__;
#line 564 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 565 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 565 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 566 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 566 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 567 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 568 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 569 "cgesvd.f"
		} else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */

#line 573 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 574 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_n__;
#line 574 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 575 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 575 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 576 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 576 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 577 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 577 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 578 "cgesvd.f"
		    maxwrk = (*m << 1) * *m + wrkbl;
#line 579 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 580 "cgesvd.f"
		} else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='A') */

#line 585 "cgesvd.f"
		    wrkbl = *m + lwork_cgelqf__;
/* Computing MAX */
#line 586 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *m + lwork_cunglq_n__;
#line 586 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 587 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cgebrd__;
#line 587 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 588 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 588 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
/* Computing MAX */
#line 589 "cgesvd.f"
		    i__2 = wrkbl, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 589 "cgesvd.f"
		    wrkbl = max(i__2,i__3);
#line 590 "cgesvd.f"
		    maxwrk = *m * *m + wrkbl;
#line 591 "cgesvd.f"
		    minwrk = (*m << 1) + *n;
#line 592 "cgesvd.f"
		}
#line 593 "cgesvd.f"
	    } else {

/*              Path 10t(N greater than M, but not much larger) */

#line 597 "cgesvd.f"
		cgebrd_(m, n, &a[a_offset], lda, &s[1], dum, cdum, cdum, cdum,
			 &c_n1, &ierr);
#line 599 "cgesvd.f"
		lwork_cgebrd__ = (integer) cdum[0].r;
#line 600 "cgesvd.f"
		maxwrk = (*m << 1) + lwork_cgebrd__;
#line 601 "cgesvd.f"
		if (wntvs || wntvo) {
/*                Compute space needed for CUNGBR P */
#line 603 "cgesvd.f"
		    cungbr_("P", m, n, m, &a[a_offset], n, cdum, cdum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 605 "cgesvd.f"
		    lwork_cungbr_p__ = (integer) cdum[0].r;
/* Computing MAX */
#line 606 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 606 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 607 "cgesvd.f"
		}
#line 608 "cgesvd.f"
		if (wntva) {
#line 609 "cgesvd.f"
		    cungbr_("P", n, n, m, &a[a_offset], n, cdum, cdum, &c_n1, 
			    &ierr, (ftnlen)1);
#line 611 "cgesvd.f"
		    lwork_cungbr_p__ = (integer) cdum[0].r;
/* Computing MAX */
#line 612 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_p__;
#line 612 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 613 "cgesvd.f"
		}
#line 614 "cgesvd.f"
		if (! wntun) {
/* Computing MAX */
#line 615 "cgesvd.f"
		    i__2 = maxwrk, i__3 = (*m << 1) + lwork_cungbr_q__;
#line 615 "cgesvd.f"
		    maxwrk = max(i__2,i__3);
#line 616 "cgesvd.f"
		}
#line 617 "cgesvd.f"
		minwrk = (*m << 1) + *n;
#line 618 "cgesvd.f"
	    }
#line 619 "cgesvd.f"
	}
#line 620 "cgesvd.f"
	maxwrk = max(minwrk,maxwrk);
#line 621 "cgesvd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 623 "cgesvd.f"
	if (*lwork < minwrk && ! lquery) {
#line 624 "cgesvd.f"
	    *info = -13;
#line 625 "cgesvd.f"
	}
#line 626 "cgesvd.f"
    }

#line 628 "cgesvd.f"
    if (*info != 0) {
#line 629 "cgesvd.f"
	i__2 = -(*info);
#line 629 "cgesvd.f"
	xerbla_("CGESVD", &i__2, (ftnlen)6);
#line 630 "cgesvd.f"
	return 0;
#line 631 "cgesvd.f"
    } else if (lquery) {
#line 632 "cgesvd.f"
	return 0;
#line 633 "cgesvd.f"
    }

/*     Quick return if possible */

#line 637 "cgesvd.f"
    if (*m == 0 || *n == 0) {
#line 638 "cgesvd.f"
	return 0;
#line 639 "cgesvd.f"
    }

/*     Get machine constants */

#line 643 "cgesvd.f"
    eps = slamch_("P", (ftnlen)1);
#line 644 "cgesvd.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 645 "cgesvd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 649 "cgesvd.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 650 "cgesvd.f"
    iscl = 0;
#line 651 "cgesvd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 652 "cgesvd.f"
	iscl = 1;
#line 653 "cgesvd.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 654 "cgesvd.f"
    } else if (anrm > bignum) {
#line 655 "cgesvd.f"
	iscl = 1;
#line 656 "cgesvd.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 657 "cgesvd.f"
    }

#line 659 "cgesvd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 665 "cgesvd.f"
	if (*m >= mnthr) {

#line 667 "cgesvd.f"
	    if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N') */
/*              No left singular vectors to be computed */

#line 672 "cgesvd.f"
		itau = 1;
#line 673 "cgesvd.f"
		iwork = itau + *n;

/*              Compute A=Q*R */
/*              (CWorkspace: need 2*N, prefer N+N*NB) */
/*              (RWorkspace: need 0) */

#line 679 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 679 "cgesvd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out below R */

#line 684 "cgesvd.f"
		if (*n > 1) {
#line 685 "cgesvd.f"
		    i__2 = *n - 1;
#line 685 "cgesvd.f"
		    i__3 = *n - 1;
#line 685 "cgesvd.f"
		    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 2], 
			    lda, (ftnlen)1);
#line 687 "cgesvd.f"
		}
#line 688 "cgesvd.f"
		ie = 1;
#line 689 "cgesvd.f"
		itauq = 1;
#line 690 "cgesvd.f"
		itaup = itauq + *n;
#line 691 "cgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 697 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 697 "cgesvd.f"
		cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 700 "cgesvd.f"
		ncvt = 0;
#line 701 "cgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 707 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 707 "cgesvd.f"
		    cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 709 "cgesvd.f"
		    ncvt = *n;
#line 710 "cgesvd.f"
		}
#line 711 "cgesvd.f"
		irwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 718 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			irwork], info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 723 "cgesvd.f"
		if (wntvas) {
#line 723 "cgesvd.f"
		    clacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 723 "cgesvd.f"
		}

#line 726 "cgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

#line 732 "cgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 736 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 737 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 737 "cgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 741 "cgesvd.f"
			ldwrku = *lda;
#line 742 "cgesvd.f"
			ldwrkr = *lda;
#line 743 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 743 "cgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 743 "cgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 747 "cgesvd.f"
			    ldwrku = *lda;
#line 748 "cgesvd.f"
			    ldwrkr = *n;
#line 749 "cgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 753 "cgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 754 "cgesvd.f"
			    ldwrkr = *n;
#line 755 "cgesvd.f"
			}
#line 755 "cgesvd.f"
		    }
#line 756 "cgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 757 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 763 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 763 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 768 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 769 "cgesvd.f"
		    i__2 = *n - 1;
#line 769 "cgesvd.f"
		    i__3 = *n - 1;
#line 769 "cgesvd.f"
		    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1], &
			    ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 776 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 776 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 778 "cgesvd.f"
		    ie = 1;
#line 779 "cgesvd.f"
		    itauq = itau;
#line 780 "cgesvd.f"
		    itaup = itauq + *n;
#line 781 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 787 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 787 "cgesvd.f"
		    cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: need 0) */

#line 795 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 795 "cgesvd.f"
		    cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 798 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 805 "cgesvd.f"
		    cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &work[ir], &ldwrkr, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 808 "cgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 815 "cgesvd.f"
		    i__2 = *m;
#line 815 "cgesvd.f"
		    i__3 = ldwrku;
#line 815 "cgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 816 "cgesvd.f"
			i__4 = *m - i__ + 1;
#line 816 "cgesvd.f"
			chunk = min(i__4,ldwrku);
#line 817 "cgesvd.f"
			cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 820 "cgesvd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 822 "cgesvd.f"
/* L10: */
#line 822 "cgesvd.f"
		    }

#line 824 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 828 "cgesvd.f"
		    ie = 1;
#line 829 "cgesvd.f"
		    itauq = 1;
#line 830 "cgesvd.f"
		    itaup = itauq + *n;
#line 831 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*                 (RWorkspace: N) */

#line 837 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 837 "cgesvd.f"
		    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 845 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 845 "cgesvd.f"
		    cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 847 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (CWorkspace: need 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 854 "cgesvd.f"
		    cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &a[a_offset], lda, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 857 "cgesvd.f"
		}

#line 859 "cgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 865 "cgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 869 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 870 "cgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 870 "cgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 874 "cgesvd.f"
			ldwrku = *lda;
#line 875 "cgesvd.f"
			ldwrkr = *lda;
#line 876 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 876 "cgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 876 "cgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 880 "cgesvd.f"
			    ldwrku = *lda;
#line 881 "cgesvd.f"
			    ldwrkr = *n;
#line 882 "cgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 886 "cgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 887 "cgesvd.f"
			    ldwrkr = *n;
#line 888 "cgesvd.f"
			}
#line 888 "cgesvd.f"
		    }
#line 889 "cgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 890 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 896 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 896 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 901 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 902 "cgesvd.f"
		    if (*n > 1) {
#line 902 "cgesvd.f"
			i__3 = *n - 1;
#line 902 "cgesvd.f"
			i__2 = *n - 1;
#line 902 "cgesvd.f"
			claset_("L", &i__3, &i__2, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 902 "cgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 910 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 910 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 912 "cgesvd.f"
		    ie = 1;
#line 913 "cgesvd.f"
		    itauq = itau;
#line 914 "cgesvd.f"
		    itaup = itauq + *n;
#line 915 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 921 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 921 "cgesvd.f"
		    cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 924 "cgesvd.f"
		    clacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 930 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 930 "cgesvd.f"
		    cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 938 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 938 "cgesvd.f"
		    cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 940 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 948 "cgesvd.f"
		    cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, cdum, &c__1,
			     &rwork[irwork], info, (ftnlen)1);
#line 951 "cgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 958 "cgesvd.f"
		    i__3 = *m;
#line 958 "cgesvd.f"
		    i__2 = ldwrku;
#line 958 "cgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 959 "cgesvd.f"
			i__4 = *m - i__ + 1;
#line 959 "cgesvd.f"
			chunk = min(i__4,ldwrku);
#line 960 "cgesvd.f"
			cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 963 "cgesvd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 965 "cgesvd.f"
/* L20: */
#line 965 "cgesvd.f"
		    }

#line 967 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 971 "cgesvd.f"
		    itau = 1;
#line 972 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 978 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 978 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 983 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 984 "cgesvd.f"
		    if (*n > 1) {
#line 984 "cgesvd.f"
			i__2 = *n - 1;
#line 984 "cgesvd.f"
			i__3 = *n - 1;
#line 984 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 984 "cgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 992 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 992 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 994 "cgesvd.f"
		    ie = 1;
#line 995 "cgesvd.f"
		    itauq = itau;
#line 996 "cgesvd.f"
		    itaup = itauq + *n;
#line 997 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                 (RWorkspace: N) */

#line 1003 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1003 "cgesvd.f"
		    cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                 (RWorkspace: 0) */

#line 1011 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1011 "cgesvd.f"
		    cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 1019 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1019 "cgesvd.f"
		    cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1021 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 1029 "cgesvd.f"
		    cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, 
			    &rwork[irwork], info, (ftnlen)1);

#line 1033 "cgesvd.f"
		}

#line 1035 "cgesvd.f"
	    } else if (wntus) {

#line 1037 "cgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

#line 1043 "cgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1047 "cgesvd.f"
			ir = 1;
#line 1048 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1052 "cgesvd.f"
			    ldwrkr = *lda;
#line 1053 "cgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1057 "cgesvd.f"
			    ldwrkr = *n;
#line 1058 "cgesvd.f"
			}
#line 1059 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1060 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1066 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1066 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1071 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1073 "cgesvd.f"
			i__2 = *n - 1;
#line 1073 "cgesvd.f"
			i__3 = *n - 1;
#line 1073 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1080 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1080 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1082 "cgesvd.f"
			ie = 1;
#line 1083 "cgesvd.f"
			itauq = itau;
#line 1084 "cgesvd.f"
			itaup = itauq + *n;
#line 1085 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1091 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1091 "cgesvd.f"
			cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1100 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1100 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1103 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1110 "cgesvd.f"
			cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1119 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1122 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1126 "cgesvd.f"
			itau = 1;
#line 1127 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1133 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1133 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1135 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1141 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1141 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1143 "cgesvd.f"
			ie = 1;
#line 1144 "cgesvd.f"
			itauq = itau;
#line 1145 "cgesvd.f"
			itaup = itauq + *n;
#line 1146 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1150 "cgesvd.f"
			if (*n > 1) {
#line 1151 "cgesvd.f"
			    i__2 = *n - 1;
#line 1151 "cgesvd.f"
			    i__3 = *n - 1;
#line 1151 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1153 "cgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1159 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1159 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1167 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1167 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1170 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1177 "cgesvd.f"
			cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1181 "cgesvd.f"
		    }

#line 1183 "cgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

#line 1189 "cgesvd.f"
		    if (*lwork >= (*n << 1) * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1193 "cgesvd.f"
			iu = 1;
#line 1194 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1198 "cgesvd.f"
			    ldwrku = *lda;
#line 1199 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1200 "cgesvd.f"
			    ldwrkr = *lda;
#line 1201 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1205 "cgesvd.f"
			    ldwrku = *lda;
#line 1206 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1207 "cgesvd.f"
			    ldwrkr = *n;
#line 1208 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1212 "cgesvd.f"
			    ldwrku = *n;
#line 1213 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1214 "cgesvd.f"
			    ldwrkr = *n;
#line 1215 "cgesvd.f"
			}
#line 1216 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1217 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1223 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1223 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1228 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1230 "cgesvd.f"
			i__2 = *n - 1;
#line 1230 "cgesvd.f"
			i__3 = *n - 1;
#line 1230 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1237 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1237 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1239 "cgesvd.f"
			ie = 1;
#line 1240 "cgesvd.f"
			itauq = itau;
#line 1241 "cgesvd.f"
			itaup = itauq + *n;
#line 1242 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1250 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1250 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1254 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1261 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1261 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1270 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1270 "cgesvd.f"
			cungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1273 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1281 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1291 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1298 "cgesvd.f"
			clacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1301 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1305 "cgesvd.f"
			itau = 1;
#line 1306 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1312 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1312 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1314 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1320 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1320 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1322 "cgesvd.f"
			ie = 1;
#line 1323 "cgesvd.f"
			itauq = itau;
#line 1324 "cgesvd.f"
			itaup = itauq + *n;
#line 1325 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1329 "cgesvd.f"
			if (*n > 1) {
#line 1330 "cgesvd.f"
			    i__2 = *n - 1;
#line 1330 "cgesvd.f"
			    i__3 = *n - 1;
#line 1330 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1332 "cgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1338 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1338 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1346 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1346 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1354 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1354 "cgesvd.f"
			cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1356 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1364 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1368 "cgesvd.f"
		    }

#line 1370 "cgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

#line 1377 "cgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1381 "cgesvd.f"
			iu = 1;
#line 1382 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1386 "cgesvd.f"
			    ldwrku = *lda;
#line 1387 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1391 "cgesvd.f"
			    ldwrku = *n;
#line 1392 "cgesvd.f"
			}
#line 1393 "cgesvd.f"
			itau = iu + ldwrku * *n;
#line 1394 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1400 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1400 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1405 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1407 "cgesvd.f"
			i__2 = *n - 1;
#line 1407 "cgesvd.f"
			i__3 = *n - 1;
#line 1407 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1414 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1414 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1416 "cgesvd.f"
			ie = 1;
#line 1417 "cgesvd.f"
			itauq = itau;
#line 1418 "cgesvd.f"
			itaup = itauq + *n;
#line 1419 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1425 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1425 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1429 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1436 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1436 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1445 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1445 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1447 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1455 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1464 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1467 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1471 "cgesvd.f"
			itau = 1;
#line 1472 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1478 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1478 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1480 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1486 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1486 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1491 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1492 "cgesvd.f"
			if (*n > 1) {
#line 1492 "cgesvd.f"
			    i__2 = *n - 1;
#line 1492 "cgesvd.f"
			    i__3 = *n - 1;
#line 1492 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1492 "cgesvd.f"
			}
#line 1495 "cgesvd.f"
			ie = 1;
#line 1496 "cgesvd.f"
			itauq = itau;
#line 1497 "cgesvd.f"
			itaup = itauq + *n;
#line 1498 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1504 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1504 "cgesvd.f"
			cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1513 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1513 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1521 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1521 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1523 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1531 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1535 "cgesvd.f"
		    }

#line 1537 "cgesvd.f"
		}

#line 1539 "cgesvd.f"
	    } else if (wntua) {

#line 1541 "cgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1547 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1547 "cgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1551 "cgesvd.f"
			ir = 1;
#line 1552 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1556 "cgesvd.f"
			    ldwrkr = *lda;
#line 1557 "cgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1561 "cgesvd.f"
			    ldwrkr = *n;
#line 1562 "cgesvd.f"
			}
#line 1563 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1564 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1570 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1570 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1572 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1576 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1578 "cgesvd.f"
			i__2 = *n - 1;
#line 1578 "cgesvd.f"
			i__3 = *n - 1;
#line 1578 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1585 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1585 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1587 "cgesvd.f"
			ie = 1;
#line 1588 "cgesvd.f"
			itauq = itau;
#line 1589 "cgesvd.f"
			itaup = itauq + *n;
#line 1590 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1596 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1596 "cgesvd.f"
			cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1605 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1605 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1608 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1615 "cgesvd.f"
			cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1624 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1629 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1631 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1635 "cgesvd.f"
			itau = 1;
#line 1636 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1642 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1642 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1644 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1650 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1650 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1652 "cgesvd.f"
			ie = 1;
#line 1653 "cgesvd.f"
			itauq = itau;
#line 1654 "cgesvd.f"
			itaup = itauq + *n;
#line 1655 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1659 "cgesvd.f"
			if (*n > 1) {
#line 1660 "cgesvd.f"
			    i__2 = *n - 1;
#line 1660 "cgesvd.f"
			    i__3 = *n - 1;
#line 1660 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1662 "cgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1668 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1668 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1677 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1677 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1680 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1687 "cgesvd.f"
			cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1691 "cgesvd.f"
		    }

#line 1693 "cgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1699 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1699 "cgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1703 "cgesvd.f"
			iu = 1;
#line 1704 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1708 "cgesvd.f"
			    ldwrku = *lda;
#line 1709 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1710 "cgesvd.f"
			    ldwrkr = *lda;
#line 1711 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1715 "cgesvd.f"
			    ldwrku = *lda;
#line 1716 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1717 "cgesvd.f"
			    ldwrkr = *n;
#line 1718 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1722 "cgesvd.f"
			    ldwrku = *n;
#line 1723 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1724 "cgesvd.f"
			    ldwrkr = *n;
#line 1725 "cgesvd.f"
			}
#line 1726 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1727 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1733 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1733 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1735 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1741 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1741 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1746 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1748 "cgesvd.f"
			i__2 = *n - 1;
#line 1748 "cgesvd.f"
			i__3 = *n - 1;
#line 1748 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1750 "cgesvd.f"
			ie = 1;
#line 1751 "cgesvd.f"
			itauq = itau;
#line 1752 "cgesvd.f"
			itaup = itauq + *n;
#line 1753 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1761 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1761 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1765 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1772 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1772 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1781 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1781 "cgesvd.f"
			cungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1784 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1792 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1802 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1807 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1811 "cgesvd.f"
			clacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1814 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1818 "cgesvd.f"
			itau = 1;
#line 1819 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1825 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1825 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1827 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1833 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1833 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1835 "cgesvd.f"
			ie = 1;
#line 1836 "cgesvd.f"
			itauq = itau;
#line 1837 "cgesvd.f"
			itaup = itauq + *n;
#line 1838 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1842 "cgesvd.f"
			if (*n > 1) {
#line 1843 "cgesvd.f"
			    i__2 = *n - 1;
#line 1843 "cgesvd.f"
			    i__3 = *n - 1;
#line 1843 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[
				    a_dim1 + 2], lda, (ftnlen)1);
#line 1845 "cgesvd.f"
			}

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1851 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1851 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1860 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1860 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1868 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1868 "cgesvd.f"
			cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1870 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1878 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1882 "cgesvd.f"
		    }

#line 1884 "cgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1891 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1891 "cgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1895 "cgesvd.f"
			iu = 1;
#line 1896 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1900 "cgesvd.f"
			    ldwrku = *lda;
#line 1901 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1905 "cgesvd.f"
			    ldwrku = *n;
#line 1906 "cgesvd.f"
			}
#line 1907 "cgesvd.f"
			itau = iu + ldwrku * *n;
#line 1908 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1914 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1914 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1916 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1922 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1922 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1927 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1929 "cgesvd.f"
			i__2 = *n - 1;
#line 1929 "cgesvd.f"
			i__3 = *n - 1;
#line 1929 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1931 "cgesvd.f"
			ie = 1;
#line 1932 "cgesvd.f"
			itauq = itau;
#line 1933 "cgesvd.f"
			itaup = itauq + *n;
#line 1934 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1940 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1940 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1944 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1951 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1951 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: need   0) */

#line 1960 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1960 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1962 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1970 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1979 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1984 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1986 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1990 "cgesvd.f"
			itau = 1;
#line 1991 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1997 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1997 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1999 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2005 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2005 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 2010 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 2011 "cgesvd.f"
			if (*n > 1) {
#line 2011 "cgesvd.f"
			    i__2 = *n - 1;
#line 2011 "cgesvd.f"
			    i__3 = *n - 1;
#line 2011 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 2011 "cgesvd.f"
			}
#line 2014 "cgesvd.f"
			ie = 1;
#line 2015 "cgesvd.f"
			itauq = itau;
#line 2016 "cgesvd.f"
			itaup = itauq + *n;
#line 2017 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 2023 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2023 "cgesvd.f"
			cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2032 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2032 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2040 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2040 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 2042 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2050 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2054 "cgesvd.f"
		    }

#line 2056 "cgesvd.f"
		}

#line 2058 "cgesvd.f"
	    }

#line 2060 "cgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 2067 "cgesvd.f"
	    ie = 1;
#line 2068 "cgesvd.f"
	    itauq = 1;
#line 2069 "cgesvd.f"
	    itaup = itauq + *n;
#line 2070 "cgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 2076 "cgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 2076 "cgesvd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 2079 "cgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB) */
/*              (RWorkspace: 0) */

#line 2086 "cgesvd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 2087 "cgesvd.f"
		if (wntus) {
#line 2087 "cgesvd.f"
		    ncu = *n;
#line 2087 "cgesvd.f"
		}
#line 2089 "cgesvd.f"
		if (wntua) {
#line 2089 "cgesvd.f"
		    ncu = *m;
#line 2089 "cgesvd.f"
		}
#line 2091 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2091 "cgesvd.f"
		cungbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2093 "cgesvd.f"
	    }
#line 2094 "cgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2101 "cgesvd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2102 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2102 "cgesvd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2104 "cgesvd.f"
	    }
#line 2105 "cgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 2112 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2112 "cgesvd.f"
		cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2114 "cgesvd.f"
	    }
#line 2115 "cgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2122 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2122 "cgesvd.f"
		cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2124 "cgesvd.f"
	    }
#line 2125 "cgesvd.f"
	    irwork = ie + *n;
#line 2126 "cgesvd.f"
	    if (wntuas || wntuo) {
#line 2126 "cgesvd.f"
		nru = *m;
#line 2126 "cgesvd.f"
	    }
#line 2128 "cgesvd.f"
	    if (wntun) {
#line 2128 "cgesvd.f"
		nru = 0;
#line 2128 "cgesvd.f"
	    }
#line 2130 "cgesvd.f"
	    if (wntvas || wntvo) {
#line 2130 "cgesvd.f"
		ncvt = *n;
#line 2130 "cgesvd.f"
	    }
#line 2132 "cgesvd.f"
	    if (wntvn) {
#line 2132 "cgesvd.f"
		ncvt = 0;
#line 2132 "cgesvd.f"
	    }
#line 2134 "cgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2142 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2145 "cgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2153 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2156 "cgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2164 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2167 "cgesvd.f"
	    }

#line 2169 "cgesvd.f"
	}

#line 2171 "cgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2177 "cgesvd.f"
	if (*n >= mnthr) {

#line 2179 "cgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2184 "cgesvd.f"
		itau = 1;
#line 2185 "cgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 2191 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2191 "cgesvd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2196 "cgesvd.f"
		i__2 = *m - 1;
#line 2196 "cgesvd.f"
		i__3 = *m - 1;
#line 2196 "cgesvd.f"
		claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 2198 "cgesvd.f"
		ie = 1;
#line 2199 "cgesvd.f"
		itauq = 1;
#line 2200 "cgesvd.f"
		itaup = itauq + *m;
#line 2201 "cgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 2207 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2207 "cgesvd.f"
		cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2210 "cgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2216 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2216 "cgesvd.f"
		    cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2218 "cgesvd.f"
		}
#line 2219 "cgesvd.f"
		irwork = ie + *m;
#line 2220 "cgesvd.f"
		nru = 0;
#line 2221 "cgesvd.f"
		if (wntuo || wntuas) {
#line 2221 "cgesvd.f"
		    nru = *m;
#line 2221 "cgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2229 "cgesvd.f"
		cbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &rwork[ie], cdum, &
			c__1, &a[a_offset], lda, cdum, &c__1, &rwork[irwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2234 "cgesvd.f"
		if (wntuas) {
#line 2234 "cgesvd.f"
		    clacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2234 "cgesvd.f"
		}

#line 2237 "cgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

#line 2243 "cgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2247 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2248 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 2248 "cgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2252 "cgesvd.f"
			ldwrku = *lda;
#line 2253 "cgesvd.f"
			chunk = *n;
#line 2254 "cgesvd.f"
			ldwrkr = *lda;
#line 2255 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2255 "cgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 2255 "cgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2259 "cgesvd.f"
			    ldwrku = *lda;
#line 2260 "cgesvd.f"
			    chunk = *n;
#line 2261 "cgesvd.f"
			    ldwrkr = *m;
#line 2262 "cgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2266 "cgesvd.f"
			    ldwrku = *m;
#line 2267 "cgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2268 "cgesvd.f"
			    ldwrkr = *m;
#line 2269 "cgesvd.f"
			}
#line 2269 "cgesvd.f"
		    }
#line 2270 "cgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2271 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2277 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2277 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2282 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2283 "cgesvd.f"
		    i__2 = *m - 1;
#line 2283 "cgesvd.f"
		    i__3 = *m - 1;
#line 2283 "cgesvd.f"
		    claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2290 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2290 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2292 "cgesvd.f"
		    ie = 1;
#line 2293 "cgesvd.f"
		    itauq = itau;
#line 2294 "cgesvd.f"
		    itaup = itauq + *m;
#line 2295 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2301 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2301 "cgesvd.f"
		    cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2309 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2309 "cgesvd.f"
		    cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2312 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2319 "cgesvd.f"
		    cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &work[
			    ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2322 "cgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N) */
/*                 (RWorkspace: 0) */

#line 2329 "cgesvd.f"
		    i__2 = *n;
#line 2329 "cgesvd.f"
		    i__3 = chunk;
#line 2329 "cgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2330 "cgesvd.f"
			i__4 = *n - i__ + 1;
#line 2330 "cgesvd.f"
			blk = min(i__4,chunk);
#line 2331 "cgesvd.f"
			cgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2334 "cgesvd.f"
			clacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2336 "cgesvd.f"
/* L30: */
#line 2336 "cgesvd.f"
		    }

#line 2338 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2342 "cgesvd.f"
		    ie = 1;
#line 2343 "cgesvd.f"
		    itauq = 1;
#line 2344 "cgesvd.f"
		    itaup = itauq + *m;
#line 2345 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*                 (RWorkspace: need M) */

#line 2351 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2351 "cgesvd.f"
		    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2359 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2359 "cgesvd.f"
		    cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2361 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2368 "cgesvd.f"
		    cbdsqr_("L", m, n, &c__0, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 2371 "cgesvd.f"
		}

#line 2373 "cgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 2379 "cgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2383 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2384 "cgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 2384 "cgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2388 "cgesvd.f"
			ldwrku = *lda;
#line 2389 "cgesvd.f"
			chunk = *n;
#line 2390 "cgesvd.f"
			ldwrkr = *lda;
#line 2391 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2391 "cgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 2391 "cgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2395 "cgesvd.f"
			    ldwrku = *lda;
#line 2396 "cgesvd.f"
			    chunk = *n;
#line 2397 "cgesvd.f"
			    ldwrkr = *m;
#line 2398 "cgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2402 "cgesvd.f"
			    ldwrku = *m;
#line 2403 "cgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2404 "cgesvd.f"
			    ldwrkr = *m;
#line 2405 "cgesvd.f"
			}
#line 2405 "cgesvd.f"
		    }
#line 2406 "cgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2407 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2413 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2413 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2418 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2419 "cgesvd.f"
		    i__3 = *m - 1;
#line 2419 "cgesvd.f"
		    i__2 = *m - 1;
#line 2419 "cgesvd.f"
		    claset_("U", &i__3, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2426 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2426 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2428 "cgesvd.f"
		    ie = 1;
#line 2429 "cgesvd.f"
		    itauq = itau;
#line 2430 "cgesvd.f"
		    itaup = itauq + *m;
#line 2431 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2437 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2437 "cgesvd.f"
		    cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2440 "cgesvd.f"
		    clacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2446 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2446 "cgesvd.f"
		    cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2454 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2454 "cgesvd.f"
		    cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2456 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2464 "cgesvd.f"
		    cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[ir],
			     &ldwrkr, &u[u_offset], ldu, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2467 "cgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N)) */
/*                 (RWorkspace: 0) */

#line 2474 "cgesvd.f"
		    i__3 = *n;
#line 2474 "cgesvd.f"
		    i__2 = chunk;
#line 2474 "cgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2475 "cgesvd.f"
			i__4 = *n - i__ + 1;
#line 2475 "cgesvd.f"
			blk = min(i__4,chunk);
#line 2476 "cgesvd.f"
			cgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2479 "cgesvd.f"
			clacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2481 "cgesvd.f"
/* L40: */
#line 2481 "cgesvd.f"
		    }

#line 2483 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2487 "cgesvd.f"
		    itau = 1;
#line 2488 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2494 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2494 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2499 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2500 "cgesvd.f"
		    i__2 = *m - 1;
#line 2500 "cgesvd.f"
		    i__3 = *m - 1;
#line 2500 "cgesvd.f"
		    claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2507 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2507 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2509 "cgesvd.f"
		    ie = 1;
#line 2510 "cgesvd.f"
		    itauq = itau;
#line 2511 "cgesvd.f"
		    itaup = itauq + *m;
#line 2512 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2518 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2518 "cgesvd.f"
		    cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                 (RWorkspace: 0) */

#line 2526 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2526 "cgesvd.f"
		    cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2534 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2534 "cgesvd.f"
		    cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2536 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2544 "cgesvd.f"
		    cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			    rwork[irwork], info, (ftnlen)1);

#line 2547 "cgesvd.f"
		}

#line 2549 "cgesvd.f"
	    } else if (wntvs) {

#line 2551 "cgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

#line 2557 "cgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2561 "cgesvd.f"
			ir = 1;
#line 2562 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2566 "cgesvd.f"
			    ldwrkr = *lda;
#line 2567 "cgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2571 "cgesvd.f"
			    ldwrkr = *m;
#line 2572 "cgesvd.f"
			}
#line 2573 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2574 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2580 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2580 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2585 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2587 "cgesvd.f"
			i__2 = *m - 1;
#line 2587 "cgesvd.f"
			i__3 = *m - 1;
#line 2587 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2594 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2594 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2596 "cgesvd.f"
			ie = 1;
#line 2597 "cgesvd.f"
			itauq = itau;
#line 2598 "cgesvd.f"
			itaup = itauq + *m;
#line 2599 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2605 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2605 "cgesvd.f"
			cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2615 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2615 "cgesvd.f"
			cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2618 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2625 "cgesvd.f"
			cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2634 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2637 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2641 "cgesvd.f"
			itau = 1;
#line 2642 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2648 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2648 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2653 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2659 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2659 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2661 "cgesvd.f"
			ie = 1;
#line 2662 "cgesvd.f"
			itauq = itau;
#line 2663 "cgesvd.f"
			itaup = itauq + *m;
#line 2664 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2668 "cgesvd.f"
			i__2 = *m - 1;
#line 2668 "cgesvd.f"
			i__3 = *m - 1;
#line 2668 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2675 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2675 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2683 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2683 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2686 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2693 "cgesvd.f"
			cbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 2697 "cgesvd.f"
		    }

#line 2699 "cgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

#line 2705 "cgesvd.f"
		    if (*lwork >= (*m << 1) * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2709 "cgesvd.f"
			iu = 1;
#line 2710 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2714 "cgesvd.f"
			    ldwrku = *lda;
#line 2715 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2716 "cgesvd.f"
			    ldwrkr = *lda;
#line 2717 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2721 "cgesvd.f"
			    ldwrku = *lda;
#line 2722 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2723 "cgesvd.f"
			    ldwrkr = *m;
#line 2724 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2728 "cgesvd.f"
			    ldwrku = *m;
#line 2729 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2730 "cgesvd.f"
			    ldwrkr = *m;
#line 2731 "cgesvd.f"
			}
#line 2732 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2733 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2739 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2739 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2744 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2746 "cgesvd.f"
			i__2 = *m - 1;
#line 2746 "cgesvd.f"
			i__3 = *m - 1;
#line 2746 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2753 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2753 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2755 "cgesvd.f"
			ie = 1;
#line 2756 "cgesvd.f"
			itauq = itau;
#line 2757 "cgesvd.f"
			itaup = itauq + *m;
#line 2758 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 2766 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2766 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2770 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2778 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2778 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2786 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2786 "cgesvd.f"
			cungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2789 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2797 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2807 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2814 "cgesvd.f"
			clacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2817 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2821 "cgesvd.f"
			itau = 1;
#line 2822 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2828 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2828 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2830 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2836 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2836 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2838 "cgesvd.f"
			ie = 1;
#line 2839 "cgesvd.f"
			itauq = itau;
#line 2840 "cgesvd.f"
			itaup = itauq + *m;
#line 2841 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2845 "cgesvd.f"
			i__2 = *m - 1;
#line 2845 "cgesvd.f"
			i__3 = *m - 1;
#line 2845 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2852 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2852 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2860 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2860 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2868 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2868 "cgesvd.f"
			cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2870 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2878 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2882 "cgesvd.f"
		    }

#line 2884 "cgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

#line 2891 "cgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2895 "cgesvd.f"
			iu = 1;
#line 2896 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2900 "cgesvd.f"
			    ldwrku = *lda;
#line 2901 "cgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2905 "cgesvd.f"
			    ldwrku = *m;
#line 2906 "cgesvd.f"
			}
#line 2907 "cgesvd.f"
			itau = iu + ldwrku * *m;
#line 2908 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2914 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2914 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2919 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2921 "cgesvd.f"
			i__2 = *m - 1;
#line 2921 "cgesvd.f"
			i__3 = *m - 1;
#line 2921 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2928 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2928 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2930 "cgesvd.f"
			ie = 1;
#line 2931 "cgesvd.f"
			itauq = itau;
#line 2932 "cgesvd.f"
			itaup = itauq + *m;
#line 2933 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2939 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2939 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2943 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2951 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2951 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2959 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2959 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2961 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2969 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2978 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2981 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2985 "cgesvd.f"
			itau = 1;
#line 2986 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2992 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2992 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2994 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3000 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3000 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3005 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3006 "cgesvd.f"
			i__2 = *m - 1;
#line 3006 "cgesvd.f"
			i__3 = *m - 1;
#line 3006 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3008 "cgesvd.f"
			ie = 1;
#line 3009 "cgesvd.f"
			itauq = itau;
#line 3010 "cgesvd.f"
			itaup = itauq + *m;
#line 3011 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3017 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3017 "cgesvd.f"
			cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3026 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3026 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3034 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3034 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3036 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3044 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3048 "cgesvd.f"
		    }

#line 3050 "cgesvd.f"
		}

#line 3052 "cgesvd.f"
	    } else if (wntva) {

#line 3054 "cgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 3060 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3060 "cgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3064 "cgesvd.f"
			ir = 1;
#line 3065 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 3069 "cgesvd.f"
			    ldwrkr = *lda;
#line 3070 "cgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 3074 "cgesvd.f"
			    ldwrkr = *m;
#line 3075 "cgesvd.f"
			}
#line 3076 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3077 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3083 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3083 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3085 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 3089 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 3091 "cgesvd.f"
			i__2 = *m - 1;
#line 3091 "cgesvd.f"
			i__3 = *m - 1;
#line 3091 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3098 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3098 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3100 "cgesvd.f"
			ie = 1;
#line 3101 "cgesvd.f"
			itauq = itau;
#line 3102 "cgesvd.f"
			itaup = itauq + *m;
#line 3103 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3109 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3109 "cgesvd.f"
			cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3119 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3119 "cgesvd.f"
			cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3122 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3129 "cgesvd.f"
			cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3138 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3143 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3145 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3149 "cgesvd.f"
			itau = 1;
#line 3150 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3156 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3156 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3158 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3164 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3164 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3166 "cgesvd.f"
			ie = 1;
#line 3167 "cgesvd.f"
			itauq = itau;
#line 3168 "cgesvd.f"
			itaup = itauq + *m;
#line 3169 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3173 "cgesvd.f"
			i__2 = *m - 1;
#line 3173 "cgesvd.f"
			i__3 = *m - 1;
#line 3173 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3180 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3180 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3189 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3189 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3192 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3199 "cgesvd.f"
			cbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 3203 "cgesvd.f"
		    }

#line 3205 "cgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3211 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3211 "cgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3215 "cgesvd.f"
			iu = 1;
#line 3216 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3220 "cgesvd.f"
			    ldwrku = *lda;
#line 3221 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3222 "cgesvd.f"
			    ldwrkr = *lda;
#line 3223 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3227 "cgesvd.f"
			    ldwrku = *lda;
#line 3228 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3229 "cgesvd.f"
			    ldwrkr = *m;
#line 3230 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3234 "cgesvd.f"
			    ldwrku = *m;
#line 3235 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3236 "cgesvd.f"
			    ldwrkr = *m;
#line 3237 "cgesvd.f"
			}
#line 3238 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3239 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3245 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3245 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3247 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3253 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3253 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3258 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3260 "cgesvd.f"
			i__2 = *m - 1;
#line 3260 "cgesvd.f"
			i__3 = *m - 1;
#line 3260 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3262 "cgesvd.f"
			ie = 1;
#line 3263 "cgesvd.f"
			itauq = itau;
#line 3264 "cgesvd.f"
			itaup = itauq + *m;
#line 3265 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 3273 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3273 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3277 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3285 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3285 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3293 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3293 "cgesvd.f"
			cungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3296 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3304 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3314 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3319 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3323 "cgesvd.f"
			clacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3326 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3330 "cgesvd.f"
			itau = 1;
#line 3331 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3337 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3337 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3339 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3345 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3345 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3347 "cgesvd.f"
			ie = 1;
#line 3348 "cgesvd.f"
			itauq = itau;
#line 3349 "cgesvd.f"
			itaup = itauq + *m;
#line 3350 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3354 "cgesvd.f"
			i__2 = *m - 1;
#line 3354 "cgesvd.f"
			i__3 = *m - 1;
#line 3354 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3361 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3361 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3370 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3370 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3378 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3378 "cgesvd.f"
			cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3380 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3388 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3392 "cgesvd.f"
		    }

#line 3394 "cgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3401 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3401 "cgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3405 "cgesvd.f"
			iu = 1;
#line 3406 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3410 "cgesvd.f"
			    ldwrku = *lda;
#line 3411 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3415 "cgesvd.f"
			    ldwrku = *m;
#line 3416 "cgesvd.f"
			}
#line 3417 "cgesvd.f"
			itau = iu + ldwrku * *m;
#line 3418 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3424 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3424 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3426 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3432 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3432 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3437 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3439 "cgesvd.f"
			i__2 = *m - 1;
#line 3439 "cgesvd.f"
			i__3 = *m - 1;
#line 3439 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3441 "cgesvd.f"
			ie = 1;
#line 3442 "cgesvd.f"
			itauq = itau;
#line 3443 "cgesvd.f"
			itaup = itauq + *m;
#line 3444 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3450 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3450 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3454 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3461 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3461 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3469 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3469 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3471 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3479 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3488 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3493 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3495 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3499 "cgesvd.f"
			itau = 1;
#line 3500 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3506 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3506 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3508 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3514 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3514 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3519 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3520 "cgesvd.f"
			i__2 = *m - 1;
#line 3520 "cgesvd.f"
			i__3 = *m - 1;
#line 3520 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3522 "cgesvd.f"
			ie = 1;
#line 3523 "cgesvd.f"
			itauq = itau;
#line 3524 "cgesvd.f"
			itaup = itauq + *m;
#line 3525 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3531 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3531 "cgesvd.f"
			cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3540 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3540 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3548 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3548 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3550 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3558 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3562 "cgesvd.f"
		    }

#line 3564 "cgesvd.f"
		}

#line 3566 "cgesvd.f"
	    }

#line 3568 "cgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3575 "cgesvd.f"
	    ie = 1;
#line 3576 "cgesvd.f"
	    itauq = 1;
#line 3577 "cgesvd.f"
	    itaup = itauq + *m;
#line 3578 "cgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 3584 "cgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3584 "cgesvd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 3587 "cgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3594 "cgesvd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3595 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3595 "cgesvd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3597 "cgesvd.f"
	    }
#line 3598 "cgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB) */
/*              (RWorkspace: 0) */

#line 3605 "cgesvd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3606 "cgesvd.f"
		if (wntva) {
#line 3606 "cgesvd.f"
		    nrvt = *n;
#line 3606 "cgesvd.f"
		}
#line 3608 "cgesvd.f"
		if (wntvs) {
#line 3608 "cgesvd.f"
		    nrvt = *m;
#line 3608 "cgesvd.f"
		}
#line 3610 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3610 "cgesvd.f"
		cungbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3612 "cgesvd.f"
	    }
#line 3613 "cgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3620 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3620 "cgesvd.f"
		cungbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3622 "cgesvd.f"
	    }
#line 3623 "cgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 3630 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3630 "cgesvd.f"
		cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3632 "cgesvd.f"
	    }
#line 3633 "cgesvd.f"
	    irwork = ie + *m;
#line 3634 "cgesvd.f"
	    if (wntuas || wntuo) {
#line 3634 "cgesvd.f"
		nru = *m;
#line 3634 "cgesvd.f"
	    }
#line 3636 "cgesvd.f"
	    if (wntun) {
#line 3636 "cgesvd.f"
		nru = 0;
#line 3636 "cgesvd.f"
	    }
#line 3638 "cgesvd.f"
	    if (wntvas || wntvo) {
#line 3638 "cgesvd.f"
		ncvt = *n;
#line 3638 "cgesvd.f"
	    }
#line 3640 "cgesvd.f"
	    if (wntvn) {
#line 3640 "cgesvd.f"
		ncvt = 0;
#line 3640 "cgesvd.f"
	    }
#line 3642 "cgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3650 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3653 "cgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3661 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3664 "cgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3672 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3675 "cgesvd.f"
	    }

#line 3677 "cgesvd.f"
	}

#line 3679 "cgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3683 "cgesvd.f"
    if (iscl == 1) {
#line 3684 "cgesvd.f"
	if (anrm > bignum) {
#line 3684 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3684 "cgesvd.f"
	}
#line 3687 "cgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3687 "cgesvd.f"
	    i__2 = minmn - 1;
#line 3687 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3687 "cgesvd.f"
	}
#line 3690 "cgesvd.f"
	if (anrm < smlnum) {
#line 3690 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3690 "cgesvd.f"
	}
#line 3693 "cgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3693 "cgesvd.f"
	    i__2 = minmn - 1;
#line 3693 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3693 "cgesvd.f"
	}
#line 3696 "cgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3700 "cgesvd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 3702 "cgesvd.f"
    return 0;

/*     End of CGESVD */

} /* cgesvd_ */

