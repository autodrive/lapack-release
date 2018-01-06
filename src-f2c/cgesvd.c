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
		    minwrk = (*m << 1) + *n;
#line 617 "cgesvd.f"
		}
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
		i__2 = *n - 1;
#line 684 "cgesvd.f"
		i__3 = *n - 1;
#line 684 "cgesvd.f"
		claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 686 "cgesvd.f"
		ie = 1;
#line 687 "cgesvd.f"
		itauq = 1;
#line 688 "cgesvd.f"
		itaup = itauq + *n;
#line 689 "cgesvd.f"
		iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*              (RWorkspace: need N) */

#line 695 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 695 "cgesvd.f"
		cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 698 "cgesvd.f"
		ncvt = 0;
#line 699 "cgesvd.f"
		if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 705 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 705 "cgesvd.f"
		    cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 707 "cgesvd.f"
		    ncvt = *n;
#line 708 "cgesvd.f"
		}
#line 709 "cgesvd.f"
		irwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 716 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &c__0, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			irwork], info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

#line 721 "cgesvd.f"
		if (wntvas) {
#line 721 "cgesvd.f"
		    clacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 721 "cgesvd.f"
		}

#line 724 "cgesvd.f"
	    } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

#line 730 "cgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 734 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 735 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 735 "cgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

#line 739 "cgesvd.f"
			ldwrku = *lda;
#line 740 "cgesvd.f"
			ldwrkr = *lda;
#line 741 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 741 "cgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 741 "cgesvd.f"
			if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

#line 745 "cgesvd.f"
			    ldwrku = *lda;
#line 746 "cgesvd.f"
			    ldwrkr = *n;
#line 747 "cgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

#line 751 "cgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 752 "cgesvd.f"
			    ldwrkr = *n;
#line 753 "cgesvd.f"
			}
#line 753 "cgesvd.f"
		    }
#line 754 "cgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 755 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 761 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 761 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

#line 766 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 767 "cgesvd.f"
		    i__2 = *n - 1;
#line 767 "cgesvd.f"
		    i__3 = *n - 1;
#line 767 "cgesvd.f"
		    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1], &
			    ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 774 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 774 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 776 "cgesvd.f"
		    ie = 1;
#line 777 "cgesvd.f"
		    itauq = itau;
#line 778 "cgesvd.f"
		    itaup = itauq + *n;
#line 779 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 785 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 785 "cgesvd.f"
		    cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: need 0) */

#line 793 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 793 "cgesvd.f"
		    cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 796 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 803 "cgesvd.f"
		    cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &work[ir], &ldwrkr, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 806 "cgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 813 "cgesvd.f"
		    i__2 = *m;
#line 813 "cgesvd.f"
		    i__3 = ldwrku;
#line 813 "cgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 814 "cgesvd.f"
			i__4 = *m - i__ + 1;
#line 814 "cgesvd.f"
			chunk = min(i__4,ldwrku);
#line 815 "cgesvd.f"
			cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 818 "cgesvd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 820 "cgesvd.f"
/* L10: */
#line 820 "cgesvd.f"
		    }

#line 822 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 826 "cgesvd.f"
		    ie = 1;
#line 827 "cgesvd.f"
		    itauq = 1;
#line 828 "cgesvd.f"
		    itaup = itauq + *n;
#line 829 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*                 (RWorkspace: N) */

#line 835 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 835 "cgesvd.f"
		    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 843 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 843 "cgesvd.f"
		    cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 845 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (CWorkspace: need 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 852 "cgesvd.f"
		    cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], cdum, 
			    &c__1, &a[a_offset], lda, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 855 "cgesvd.f"
		}

#line 857 "cgesvd.f"
	    } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 863 "cgesvd.f"
		if (*lwork >= *n * *n + *n * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 867 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 868 "cgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 868 "cgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 872 "cgesvd.f"
			ldwrku = *lda;
#line 873 "cgesvd.f"
			ldwrkr = *lda;
#line 874 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 874 "cgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 874 "cgesvd.f"
			if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 878 "cgesvd.f"
			    ldwrku = *lda;
#line 879 "cgesvd.f"
			    ldwrkr = *n;
#line 880 "cgesvd.f"
			} else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

#line 884 "cgesvd.f"
			    ldwrku = (*lwork - *n * *n) / *n;
#line 885 "cgesvd.f"
			    ldwrkr = *n;
#line 886 "cgesvd.f"
			}
#line 886 "cgesvd.f"
		    }
#line 887 "cgesvd.f"
		    itau = ir + ldwrkr * *n;
#line 888 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 894 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 894 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 899 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 900 "cgesvd.f"
		    if (*n > 1) {
#line 900 "cgesvd.f"
			i__3 = *n - 1;
#line 900 "cgesvd.f"
			i__2 = *n - 1;
#line 900 "cgesvd.f"
			claset_("L", &i__3, &i__2, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 900 "cgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                 (RWorkspace: 0) */

#line 908 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 908 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 910 "cgesvd.f"
		    ie = 1;
#line 911 "cgesvd.f"
		    itauq = itau;
#line 912 "cgesvd.f"
		    itaup = itauq + *n;
#line 913 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                 (RWorkspace: need N) */

#line 919 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 919 "cgesvd.f"
		    cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__3, &
			    ierr);
#line 922 "cgesvd.f"
		    clacpy_("L", n, n, &vt[vt_offset], ldvt, &work[ir], &
			    ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                 (RWorkspace: 0) */

#line 928 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 928 "cgesvd.f"
		    cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 936 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 936 "cgesvd.f"
		    cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__3, &ierr, (ftnlen)1);
#line 938 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (CWorkspace: need N*N) */
/*                 (RWorkspace: need BDSPAC) */

#line 946 "cgesvd.f"
		    cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &work[ir], &ldwrkr, cdum, &c__1,
			     &rwork[irwork], info, (ftnlen)1);
#line 949 "cgesvd.f"
		    iu = itauq;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need N*N+N, prefer N*N+M*N) */
/*                 (RWorkspace: 0) */

#line 956 "cgesvd.f"
		    i__3 = *m;
#line 956 "cgesvd.f"
		    i__2 = ldwrku;
#line 956 "cgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 957 "cgesvd.f"
			i__4 = *m - i__ + 1;
#line 957 "cgesvd.f"
			chunk = min(i__4,ldwrku);
#line 958 "cgesvd.f"
			cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1]
				, lda, &work[ir], &ldwrkr, &c_b1, &work[iu], &
				ldwrku, (ftnlen)1, (ftnlen)1);
#line 961 "cgesvd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 963 "cgesvd.f"
/* L20: */
#line 963 "cgesvd.f"
		    }

#line 965 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 969 "cgesvd.f"
		    itau = 1;
#line 970 "cgesvd.f"
		    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 976 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 976 "cgesvd.f"
		    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

#line 981 "cgesvd.f"
		    clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
			    ldvt, (ftnlen)1);
#line 982 "cgesvd.f"
		    if (*n > 1) {
#line 982 "cgesvd.f"
			i__2 = *n - 1;
#line 982 "cgesvd.f"
			i__3 = *n - 1;
#line 982 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[vt_dim1 
				+ 2], ldvt, (ftnlen)1);
#line 982 "cgesvd.f"
		    }

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*N, prefer N+N*NB) */
/*                 (RWorkspace: 0) */

#line 990 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 990 "cgesvd.f"
		    cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 992 "cgesvd.f"
		    ie = 1;
#line 993 "cgesvd.f"
		    itauq = itau;
#line 994 "cgesvd.f"
		    itaup = itauq + *n;
#line 995 "cgesvd.f"
		    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                 (RWorkspace: N) */

#line 1001 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1001 "cgesvd.f"
		    cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                 (RWorkspace: 0) */

#line 1009 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1009 "cgesvd.f"
		    cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, &
			    work[itauq], &a[a_offset], lda, &work[iwork], &
			    i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                 (RWorkspace: 0) */

#line 1017 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 1017 "cgesvd.f"
		    cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], 
			    &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1019 "cgesvd.f"
		    irwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 1027 "cgesvd.f"
		    cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
			    vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, 
			    &rwork[irwork], info, (ftnlen)1);

#line 1031 "cgesvd.f"
		}

#line 1033 "cgesvd.f"
	    } else if (wntus) {

#line 1035 "cgesvd.f"
		if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

#line 1041 "cgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1045 "cgesvd.f"
			ir = 1;
#line 1046 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1050 "cgesvd.f"
			    ldwrkr = *lda;
#line 1051 "cgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1055 "cgesvd.f"
			    ldwrkr = *n;
#line 1056 "cgesvd.f"
			}
#line 1057 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1058 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1064 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1064 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1069 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1071 "cgesvd.f"
			i__2 = *n - 1;
#line 1071 "cgesvd.f"
			i__3 = *n - 1;
#line 1071 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1078 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1078 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1080 "cgesvd.f"
			ie = 1;
#line 1081 "cgesvd.f"
			itauq = itau;
#line 1082 "cgesvd.f"
			itaup = itauq + *n;
#line 1083 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1089 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1089 "cgesvd.f"
			cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1098 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1098 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1101 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1108 "cgesvd.f"
			cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1117 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[ir], &ldwrkr, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1120 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1124 "cgesvd.f"
			itau = 1;
#line 1125 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1131 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1131 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1133 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1139 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1139 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1141 "cgesvd.f"
			ie = 1;
#line 1142 "cgesvd.f"
			itauq = itau;
#line 1143 "cgesvd.f"
			itaup = itauq + *n;
#line 1144 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1148 "cgesvd.f"
			i__2 = *n - 1;
#line 1148 "cgesvd.f"
			i__3 = *n - 1;
#line 1148 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1155 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1155 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1163 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1163 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1166 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1173 "cgesvd.f"
			cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1177 "cgesvd.f"
		    }

#line 1179 "cgesvd.f"
		} else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

#line 1185 "cgesvd.f"
		    if (*lwork >= (*n << 1) * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1189 "cgesvd.f"
			iu = 1;
#line 1190 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1194 "cgesvd.f"
			    ldwrku = *lda;
#line 1195 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1196 "cgesvd.f"
			    ldwrkr = *lda;
#line 1197 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1201 "cgesvd.f"
			    ldwrku = *lda;
#line 1202 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1203 "cgesvd.f"
			    ldwrkr = *n;
#line 1204 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1208 "cgesvd.f"
			    ldwrku = *n;
#line 1209 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1210 "cgesvd.f"
			    ldwrkr = *n;
#line 1211 "cgesvd.f"
			}
#line 1212 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1213 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1219 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1219 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1224 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1226 "cgesvd.f"
			i__2 = *n - 1;
#line 1226 "cgesvd.f"
			i__3 = *n - 1;
#line 1226 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1233 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1233 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1235 "cgesvd.f"
			ie = 1;
#line 1236 "cgesvd.f"
			itauq = itau;
#line 1237 "cgesvd.f"
			itaup = itauq + *n;
#line 1238 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1246 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1246 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1250 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1257 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1257 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1266 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1266 "cgesvd.f"
			cungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1269 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1277 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1287 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1294 "cgesvd.f"
			clacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1297 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1301 "cgesvd.f"
			itau = 1;
#line 1302 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1308 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1308 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1310 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1316 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1316 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1318 "cgesvd.f"
			ie = 1;
#line 1319 "cgesvd.f"
			itauq = itau;
#line 1320 "cgesvd.f"
			itaup = itauq + *n;
#line 1321 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1325 "cgesvd.f"
			i__2 = *n - 1;
#line 1325 "cgesvd.f"
			i__3 = *n - 1;
#line 1325 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1332 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1332 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1340 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1340 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1348 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1348 "cgesvd.f"
			cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1350 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1358 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1362 "cgesvd.f"
		    }

#line 1364 "cgesvd.f"
		} else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

#line 1371 "cgesvd.f"
		    if (*lwork >= *n * *n + *n * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 1375 "cgesvd.f"
			iu = 1;
#line 1376 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1380 "cgesvd.f"
			    ldwrku = *lda;
#line 1381 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1385 "cgesvd.f"
			    ldwrku = *n;
#line 1386 "cgesvd.f"
			}
#line 1387 "cgesvd.f"
			itau = iu + ldwrku * *n;
#line 1388 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1394 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1394 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1399 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1401 "cgesvd.f"
			i__2 = *n - 1;
#line 1401 "cgesvd.f"
			i__3 = *n - 1;
#line 1401 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1408 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1408 "cgesvd.f"
			cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1410 "cgesvd.f"
			ie = 1;
#line 1411 "cgesvd.f"
			itauq = itau;
#line 1412 "cgesvd.f"
			itaup = itauq + *n;
#line 1413 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1419 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1419 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1423 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1430 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1430 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1439 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1439 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1441 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1449 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1458 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &
				work[iu], &ldwrku, &c_b1, &u[u_offset], ldu, (
				ftnlen)1, (ftnlen)1);

#line 1461 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1465 "cgesvd.f"
			itau = 1;
#line 1466 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1472 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1472 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1474 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1480 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1480 "cgesvd.f"
			cungqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

#line 1485 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 1486 "cgesvd.f"
			if (*n > 1) {
#line 1486 "cgesvd.f"
			    i__2 = *n - 1;
#line 1486 "cgesvd.f"
			    i__3 = *n - 1;
#line 1486 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 1486 "cgesvd.f"
			}
#line 1489 "cgesvd.f"
			ie = 1;
#line 1490 "cgesvd.f"
			itauq = itau;
#line 1491 "cgesvd.f"
			itaup = itauq + *n;
#line 1492 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1498 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1498 "cgesvd.f"
			cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1507 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1507 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1515 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1515 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1517 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1525 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1529 "cgesvd.f"
		    }

#line 1531 "cgesvd.f"
		}

#line 1533 "cgesvd.f"
	    } else if (wntua) {

#line 1535 "cgesvd.f"
		if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
#line 1541 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1541 "cgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1545 "cgesvd.f"
			ir = 1;
#line 1546 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

#line 1550 "cgesvd.f"
			    ldwrkr = *lda;
#line 1551 "cgesvd.f"
			} else {

/*                       WORK(IR) is N by N */

#line 1555 "cgesvd.f"
			    ldwrkr = *n;
#line 1556 "cgesvd.f"
			}
#line 1557 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1558 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1564 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1564 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1566 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

#line 1570 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 1572 "cgesvd.f"
			i__2 = *n - 1;
#line 1572 "cgesvd.f"
			i__3 = *n - 1;
#line 1572 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 1]
				, &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1579 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1579 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1581 "cgesvd.f"
			ie = 1;
#line 1582 "cgesvd.f"
			itauq = itau;
#line 1583 "cgesvd.f"
			itaup = itauq + *n;
#line 1584 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1590 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1590 "cgesvd.f"
			cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1599 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1599 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1602 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1609 "cgesvd.f"
			cbdsqr_("U", n, &c__0, n, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &work[ir], &ldwrkr, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1618 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[ir], &ldwrkr, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1623 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1625 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1629 "cgesvd.f"
			itau = 1;
#line 1630 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1636 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1636 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1638 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1644 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1644 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1646 "cgesvd.f"
			ie = 1;
#line 1647 "cgesvd.f"
			itauq = itau;
#line 1648 "cgesvd.f"
			itaup = itauq + *n;
#line 1649 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1653 "cgesvd.f"
			i__2 = *n - 1;
#line 1653 "cgesvd.f"
			i__3 = *n - 1;
#line 1653 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1660 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1660 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1669 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1669 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;
#line 1672 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1679 "cgesvd.f"
			cbdsqr_("U", n, &c__0, m, &c__0, &s[1], &rwork[ie], 
				cdum, &c__1, &u[u_offset], ldu, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

#line 1683 "cgesvd.f"
		    }

#line 1685 "cgesvd.f"
		} else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
#line 1691 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1691 "cgesvd.f"
		    if (*lwork >= (*n << 1) * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1695 "cgesvd.f"
			iu = 1;
#line 1696 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

#line 1700 "cgesvd.f"
			    ldwrku = *lda;
#line 1701 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1702 "cgesvd.f"
			    ldwrkr = *lda;
#line 1703 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

#line 1707 "cgesvd.f"
			    ldwrku = *lda;
#line 1708 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1709 "cgesvd.f"
			    ldwrkr = *n;
#line 1710 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

#line 1714 "cgesvd.f"
			    ldwrku = *n;
#line 1715 "cgesvd.f"
			    ir = iu + ldwrku * *n;
#line 1716 "cgesvd.f"
			    ldwrkr = *n;
#line 1717 "cgesvd.f"
			}
#line 1718 "cgesvd.f"
			itau = ir + ldwrkr * *n;
#line 1719 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1725 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1725 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1727 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1733 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1733 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1738 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1740 "cgesvd.f"
			i__2 = *n - 1;
#line 1740 "cgesvd.f"
			i__3 = *n - 1;
#line 1740 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1742 "cgesvd.f"
			ie = 1;
#line 1743 "cgesvd.f"
			itauq = itau;
#line 1744 "cgesvd.f"
			itaup = itauq + *n;
#line 1745 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N, */
/*                                 prefer 2*N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need   N) */

#line 1753 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1753 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1757 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1764 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1764 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   2*N*N+3*N-1, */
/*                                 prefer 2*N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1773 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1773 "cgesvd.f"
			cungbr_("P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1776 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (CWorkspace: need 2*N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1784 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &work[
				ir], &ldwrkr, &work[iu], &ldwrku, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1794 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1799 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

#line 1803 "cgesvd.f"
			clacpy_("F", n, n, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 1806 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1810 "cgesvd.f"
			itau = 1;
#line 1811 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1817 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1817 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1819 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1825 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1825 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 1827 "cgesvd.f"
			ie = 1;
#line 1828 "cgesvd.f"
			itauq = itau;
#line 1829 "cgesvd.f"
			itaup = itauq + *n;
#line 1830 "cgesvd.f"
			iwork = itaup + *n;

/*                    Zero out below R in A */

#line 1834 "cgesvd.f"
			i__2 = *n - 1;
#line 1834 "cgesvd.f"
			i__3 = *n - 1;
#line 1834 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &a[a_dim1 + 
				2], lda, (ftnlen)1);

/*                    Bidiagonalize R in A */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1841 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1841 "cgesvd.f"
			cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1850 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1850 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &a[a_offset], lda, &
				work[itauq], &u[u_offset], ldu, &work[iwork], 
				&i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
				;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 1858 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1858 "cgesvd.f"
			cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 1860 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 1868 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &a[
				a_offset], lda, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 1872 "cgesvd.f"
		    }

#line 1874 "cgesvd.f"
		} else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
#line 1881 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *n * 3;
#line 1881 "cgesvd.f"
		    if (*lwork >= *n * *n + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 1885 "cgesvd.f"
			iu = 1;
#line 1886 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

#line 1890 "cgesvd.f"
			    ldwrku = *lda;
#line 1891 "cgesvd.f"
			} else {

/*                       WORK(IU) is N by N */

#line 1895 "cgesvd.f"
			    ldwrku = *n;
#line 1896 "cgesvd.f"
			}
#line 1897 "cgesvd.f"
			itau = iu + ldwrku * *n;
#line 1898 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1904 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1904 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1906 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1912 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1912 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

#line 1917 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 1919 "cgesvd.f"
			i__2 = *n - 1;
#line 1919 "cgesvd.f"
			i__3 = *n - 1;
#line 1919 "cgesvd.f"
			claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 1]
				, &ldwrku, (ftnlen)1);
#line 1921 "cgesvd.f"
			ie = 1;
#line 1922 "cgesvd.f"
			itauq = itau;
#line 1923 "cgesvd.f"
			itaup = itauq + *n;
#line 1924 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 1930 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1930 "cgesvd.f"
			cgebrd_(n, n, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 1934 "cgesvd.f"
			clacpy_("U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
				 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1941 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1941 "cgesvd.f"
			cungbr_("Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need   N*N+3*N-1, */
/*                                 prefer N*N+2*N+(N-1)*NB) */
/*                    (RWorkspace: need   0) */

#line 1950 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1950 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 1952 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: need BDSPAC) */

#line 1960 "cgesvd.f"
			cbdsqr_("U", n, n, n, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &work[iu], &ldwrku, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (CWorkspace: need N*N) */
/*                    (RWorkspace: 0) */

#line 1969 "cgesvd.f"
			cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &
				work[iu], &ldwrku, &c_b1, &a[a_offset], lda, (
				ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

#line 1974 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

#line 1976 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 1980 "cgesvd.f"
			itau = 1;
#line 1981 "cgesvd.f"
			iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (CWorkspace: need 2*N, prefer N+N*NB) */
/*                    (RWorkspace: 0) */

#line 1987 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1987 "cgesvd.f"
			cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 1989 "cgesvd.f"
			clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (CWorkspace: need N+M, prefer N+M*NB) */
/*                    (RWorkspace: 0) */

#line 1995 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 1995 "cgesvd.f"
			cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

#line 2000 "cgesvd.f"
			clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);
#line 2001 "cgesvd.f"
			if (*n > 1) {
#line 2001 "cgesvd.f"
			    i__2 = *n - 1;
#line 2001 "cgesvd.f"
			    i__3 = *n - 1;
#line 2001 "cgesvd.f"
			    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &vt[
				    vt_dim1 + 2], ldvt, (ftnlen)1);
#line 2001 "cgesvd.f"
			}
#line 2004 "cgesvd.f"
			ie = 1;
#line 2005 "cgesvd.f"
			itauq = itau;
#line 2006 "cgesvd.f"
			itaup = itauq + *n;
#line 2007 "cgesvd.f"
			iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB) */
/*                    (RWorkspace: need N) */

#line 2013 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2013 "cgesvd.f"
			cgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &rwork[ie],
				 &work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB) */
/*                    (RWorkspace: 0) */

#line 2022 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2022 "cgesvd.f"
			cunmbr_("Q", "R", "N", m, n, n, &vt[vt_offset], ldvt, 
				&work[itauq], &u[u_offset], ldu, &work[iwork],
				 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2030 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2030 "cgesvd.f"
			cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[
				itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
				;
#line 2032 "cgesvd.f"
			irwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2040 "cgesvd.f"
			cbdsqr_("U", n, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2044 "cgesvd.f"
		    }

#line 2046 "cgesvd.f"
		}

#line 2048 "cgesvd.f"
	    }

#line 2050 "cgesvd.f"
	} else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

#line 2057 "cgesvd.f"
	    ie = 1;
#line 2058 "cgesvd.f"
	    itauq = 1;
#line 2059 "cgesvd.f"
	    itaup = itauq + *n;
#line 2060 "cgesvd.f"
	    iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
/*           (RWorkspace: need N) */

#line 2066 "cgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 2066 "cgesvd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 2069 "cgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB) */
/*              (RWorkspace: 0) */

#line 2076 "cgesvd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 2077 "cgesvd.f"
		if (wntus) {
#line 2077 "cgesvd.f"
		    ncu = *n;
#line 2077 "cgesvd.f"
		}
#line 2079 "cgesvd.f"
		if (wntua) {
#line 2079 "cgesvd.f"
		    ncu = *m;
#line 2079 "cgesvd.f"
		}
#line 2081 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2081 "cgesvd.f"
		cungbr_("Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2083 "cgesvd.f"
	    }
#line 2084 "cgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2091 "cgesvd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2092 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2092 "cgesvd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2094 "cgesvd.f"
	    }
#line 2095 "cgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N, prefer 2*N+N*NB) */
/*              (RWorkspace: 0) */

#line 2102 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2102 "cgesvd.f"
		cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2104 "cgesvd.f"
	    }
#line 2105 "cgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*              (RWorkspace: 0) */

#line 2112 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2112 "cgesvd.f"
		cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 2114 "cgesvd.f"
	    }
#line 2115 "cgesvd.f"
	    irwork = ie + *n;
#line 2116 "cgesvd.f"
	    if (wntuas || wntuo) {
#line 2116 "cgesvd.f"
		nru = *m;
#line 2116 "cgesvd.f"
	    }
#line 2118 "cgesvd.f"
	    if (wntun) {
#line 2118 "cgesvd.f"
		nru = 0;
#line 2118 "cgesvd.f"
	    }
#line 2120 "cgesvd.f"
	    if (wntvas || wntvo) {
#line 2120 "cgesvd.f"
		ncvt = *n;
#line 2120 "cgesvd.f"
	    }
#line 2122 "cgesvd.f"
	    if (wntvn) {
#line 2122 "cgesvd.f"
		ncvt = 0;
#line 2122 "cgesvd.f"
	    }
#line 2124 "cgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2132 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2135 "cgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2143 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2146 "cgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2154 "cgesvd.f"
		cbdsqr_("U", n, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 2157 "cgesvd.f"
	    }

#line 2159 "cgesvd.f"
	}

#line 2161 "cgesvd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 2167 "cgesvd.f"
	if (*n >= mnthr) {

#line 2169 "cgesvd.f"
	    if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

#line 2174 "cgesvd.f"
		itau = 1;
#line 2175 "cgesvd.f"
		iwork = itau + *m;

/*              Compute A=L*Q */
/*              (CWorkspace: need 2*M, prefer M+M*NB) */
/*              (RWorkspace: 0) */

#line 2181 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2181 "cgesvd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 2186 "cgesvd.f"
		i__2 = *m - 1;
#line 2186 "cgesvd.f"
		i__3 = *m - 1;
#line 2186 "cgesvd.f"
		claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 2188 "cgesvd.f"
		ie = 1;
#line 2189 "cgesvd.f"
		itauq = 1;
#line 2190 "cgesvd.f"
		itaup = itauq + *m;
#line 2191 "cgesvd.f"
		iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*              (RWorkspace: need M) */

#line 2197 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 2197 "cgesvd.f"
		cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[iwork], &i__2, &ierr);
#line 2200 "cgesvd.f"
		if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2206 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2206 "cgesvd.f"
		    cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2208 "cgesvd.f"
		}
#line 2209 "cgesvd.f"
		irwork = ie + *m;
#line 2210 "cgesvd.f"
		nru = 0;
#line 2211 "cgesvd.f"
		if (wntuo || wntuas) {
#line 2211 "cgesvd.f"
		    nru = *m;
#line 2211 "cgesvd.f"
		}

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 2219 "cgesvd.f"
		cbdsqr_("U", m, &c__0, &nru, &c__0, &s[1], &rwork[ie], cdum, &
			c__1, &a[a_offset], lda, cdum, &c__1, &rwork[irwork], 
			info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

#line 2224 "cgesvd.f"
		if (wntuas) {
#line 2224 "cgesvd.f"
		    clacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2224 "cgesvd.f"
		}

#line 2227 "cgesvd.f"
	    } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

#line 2233 "cgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2237 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2238 "cgesvd.f"
		    i__2 = wrkbl, i__3 = *lda * *n;
#line 2238 "cgesvd.f"
		    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2242 "cgesvd.f"
			ldwrku = *lda;
#line 2243 "cgesvd.f"
			chunk = *n;
#line 2244 "cgesvd.f"
			ldwrkr = *lda;
#line 2245 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2245 "cgesvd.f"
			i__2 = wrkbl, i__3 = *lda * *n;
#line 2245 "cgesvd.f"
			if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2249 "cgesvd.f"
			    ldwrku = *lda;
#line 2250 "cgesvd.f"
			    chunk = *n;
#line 2251 "cgesvd.f"
			    ldwrkr = *m;
#line 2252 "cgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2256 "cgesvd.f"
			    ldwrku = *m;
#line 2257 "cgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2258 "cgesvd.f"
			    ldwrkr = *m;
#line 2259 "cgesvd.f"
			}
#line 2259 "cgesvd.f"
		    }
#line 2260 "cgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2261 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2267 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2267 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

#line 2272 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, 
			    (ftnlen)1);
#line 2273 "cgesvd.f"
		    i__2 = *m - 1;
#line 2273 "cgesvd.f"
		    i__3 = *m - 1;
#line 2273 "cgesvd.f"
		    claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
			    ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2280 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2280 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2282 "cgesvd.f"
		    ie = 1;
#line 2283 "cgesvd.f"
		    itauq = itau;
#line 2284 "cgesvd.f"
		    itaup = itauq + *m;
#line 2285 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2291 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2291 "cgesvd.f"
		    cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
			    work[itauq], &work[itaup], &work[iwork], &i__2, &
			    ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2299 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2299 "cgesvd.f"
		    cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2302 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2309 "cgesvd.f"
		    cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &work[
			    ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2312 "cgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N) */
/*                 (RWorkspace: 0) */

#line 2319 "cgesvd.f"
		    i__2 = *n;
#line 2319 "cgesvd.f"
		    i__3 = chunk;
#line 2319 "cgesvd.f"
		    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__3) {
/* Computing MIN */
#line 2320 "cgesvd.f"
			i__4 = *n - i__ + 1;
#line 2320 "cgesvd.f"
			blk = min(i__4,chunk);
#line 2321 "cgesvd.f"
			cgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2324 "cgesvd.f"
			clacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2326 "cgesvd.f"
/* L30: */
#line 2326 "cgesvd.f"
		    }

#line 2328 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2332 "cgesvd.f"
		    ie = 1;
#line 2333 "cgesvd.f"
		    itauq = 1;
#line 2334 "cgesvd.f"
		    itaup = itauq + *m;
#line 2335 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*                 (RWorkspace: need M) */

#line 2341 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2341 "cgesvd.f"
		    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2349 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2349 "cgesvd.f"
		    cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2351 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2358 "cgesvd.f"
		    cbdsqr_("L", m, n, &c__0, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, cdum, &c__1, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);

#line 2361 "cgesvd.f"
		}

#line 2363 "cgesvd.f"
	    } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 2369 "cgesvd.f"
		if (*lwork >= *m * *m + *m * 3) {

/*                 Sufficient workspace for a fast algorithm */

#line 2373 "cgesvd.f"
		    ir = 1;
/* Computing MAX */
#line 2374 "cgesvd.f"
		    i__3 = wrkbl, i__2 = *lda * *n;
#line 2374 "cgesvd.f"
		    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

#line 2378 "cgesvd.f"
			ldwrku = *lda;
#line 2379 "cgesvd.f"
			chunk = *n;
#line 2380 "cgesvd.f"
			ldwrkr = *lda;
#line 2381 "cgesvd.f"
		    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 2381 "cgesvd.f"
			i__3 = wrkbl, i__2 = *lda * *n;
#line 2381 "cgesvd.f"
			if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

#line 2385 "cgesvd.f"
			    ldwrku = *lda;
#line 2386 "cgesvd.f"
			    chunk = *n;
#line 2387 "cgesvd.f"
			    ldwrkr = *m;
#line 2388 "cgesvd.f"
			} else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

#line 2392 "cgesvd.f"
			    ldwrku = *m;
#line 2393 "cgesvd.f"
			    chunk = (*lwork - *m * *m) / *m;
#line 2394 "cgesvd.f"
			    ldwrkr = *m;
#line 2395 "cgesvd.f"
			}
#line 2395 "cgesvd.f"
		    }
#line 2396 "cgesvd.f"
		    itau = ir + ldwrkr * *m;
#line 2397 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2403 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2403 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

#line 2408 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2409 "cgesvd.f"
		    i__3 = *m - 1;
#line 2409 "cgesvd.f"
		    i__2 = *m - 1;
#line 2409 "cgesvd.f"
		    claset_("U", &i__3, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2416 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2416 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__3, &ierr);
#line 2418 "cgesvd.f"
		    ie = 1;
#line 2419 "cgesvd.f"
		    itauq = itau;
#line 2420 "cgesvd.f"
		    itaup = itauq + *m;
#line 2421 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2427 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2427 "cgesvd.f"
		    cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__3, &ierr);
#line 2430 "cgesvd.f"
		    clacpy_("U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, 
			    (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB) */
/*                 (RWorkspace: 0) */

#line 2436 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2436 "cgesvd.f"
		    cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2444 "cgesvd.f"
		    i__3 = *lwork - iwork + 1;
#line 2444 "cgesvd.f"
		    cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__3, &ierr, (ftnlen)1);
#line 2446 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (CWorkspace: need M*M) */
/*                 (RWorkspace: need BDSPAC) */

#line 2454 "cgesvd.f"
		    cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[ir],
			     &ldwrkr, &u[u_offset], ldu, cdum, &c__1, &rwork[
			    irwork], info, (ftnlen)1);
#line 2457 "cgesvd.f"
		    iu = itauq;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (CWorkspace: need M*M+M, prefer M*M+M*N)) */
/*                 (RWorkspace: 0) */

#line 2464 "cgesvd.f"
		    i__3 = *n;
#line 2464 "cgesvd.f"
		    i__2 = chunk;
#line 2464 "cgesvd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
			     i__2) {
/* Computing MIN */
#line 2465 "cgesvd.f"
			i__4 = *n - i__ + 1;
#line 2465 "cgesvd.f"
			blk = min(i__4,chunk);
#line 2466 "cgesvd.f"
			cgemm_("N", "N", m, &blk, m, &c_b2, &work[ir], &
				ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b1, &
				work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
#line 2469 "cgesvd.f"
			clacpy_("F", m, &blk, &work[iu], &ldwrku, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2471 "cgesvd.f"
/* L40: */
#line 2471 "cgesvd.f"
		    }

#line 2473 "cgesvd.f"
		} else {

/*                 Insufficient workspace for a fast algorithm */

#line 2477 "cgesvd.f"
		    itau = 1;
#line 2478 "cgesvd.f"
		    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2484 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2484 "cgesvd.f"
		    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
			    , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

#line 2489 "cgesvd.f"
		    clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			    ftnlen)1);
#line 2490 "cgesvd.f"
		    i__2 = *m - 1;
#line 2490 "cgesvd.f"
		    i__3 = *m - 1;
#line 2490 "cgesvd.f"
		    claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 << 1) 
			    + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (CWorkspace: need 2*M, prefer M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2497 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2497 "cgesvd.f"
		    cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
			    iwork], &i__2, &ierr);
#line 2499 "cgesvd.f"
		    ie = 1;
#line 2500 "cgesvd.f"
		    itauq = itau;
#line 2501 "cgesvd.f"
		    itaup = itauq + *m;
#line 2502 "cgesvd.f"
		    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                 (RWorkspace: need M) */

#line 2508 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2508 "cgesvd.f"
		    cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &work[
			    itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                 (RWorkspace: 0) */

#line 2516 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2516 "cgesvd.f"
		    cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &work[
			    itaup], &a[a_offset], lda, &work[iwork], &i__2, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                 (RWorkspace: 0) */

#line 2524 "cgesvd.f"
		    i__2 = *lwork - iwork + 1;
#line 2524 "cgesvd.f"
		    cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
			    work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2526 "cgesvd.f"
		    irwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (CWorkspace: 0) */
/*                 (RWorkspace: need BDSPAC) */

#line 2534 "cgesvd.f"
		    cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &a[
			    a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			    rwork[irwork], info, (ftnlen)1);

#line 2537 "cgesvd.f"
		}

#line 2539 "cgesvd.f"
	    } else if (wntvs) {

#line 2541 "cgesvd.f"
		if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

#line 2547 "cgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2551 "cgesvd.f"
			ir = 1;
#line 2552 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 2556 "cgesvd.f"
			    ldwrkr = *lda;
#line 2557 "cgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 2561 "cgesvd.f"
			    ldwrkr = *m;
#line 2562 "cgesvd.f"
			}
#line 2563 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2564 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2570 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2570 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 2575 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 2577 "cgesvd.f"
			i__2 = *m - 1;
#line 2577 "cgesvd.f"
			i__3 = *m - 1;
#line 2577 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2584 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2584 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2586 "cgesvd.f"
			ie = 1;
#line 2587 "cgesvd.f"
			itauq = itau;
#line 2588 "cgesvd.f"
			itaup = itauq + *m;
#line 2589 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2595 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2595 "cgesvd.f"
			cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2605 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2605 "cgesvd.f"
			cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2608 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2615 "cgesvd.f"
			cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2624 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2627 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2631 "cgesvd.f"
			itau = 1;
#line 2632 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2638 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2638 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy result to VT */

#line 2643 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2649 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2649 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2651 "cgesvd.f"
			ie = 1;
#line 2652 "cgesvd.f"
			itauq = itau;
#line 2653 "cgesvd.f"
			itaup = itauq + *m;
#line 2654 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2658 "cgesvd.f"
			i__2 = *m - 1;
#line 2658 "cgesvd.f"
			i__3 = *m - 1;
#line 2658 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2665 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2665 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2673 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2673 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 2676 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2683 "cgesvd.f"
			cbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 2687 "cgesvd.f"
		    }

#line 2689 "cgesvd.f"
		} else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

#line 2695 "cgesvd.f"
		    if (*lwork >= (*m << 1) * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2699 "cgesvd.f"
			iu = 1;
#line 2700 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 2704 "cgesvd.f"
			    ldwrku = *lda;
#line 2705 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2706 "cgesvd.f"
			    ldwrkr = *lda;
#line 2707 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 2711 "cgesvd.f"
			    ldwrku = *lda;
#line 2712 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2713 "cgesvd.f"
			    ldwrkr = *m;
#line 2714 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 2718 "cgesvd.f"
			    ldwrku = *m;
#line 2719 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 2720 "cgesvd.f"
			    ldwrkr = *m;
#line 2721 "cgesvd.f"
			}
#line 2722 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 2723 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2729 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2729 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

#line 2734 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2736 "cgesvd.f"
			i__2 = *m - 1;
#line 2736 "cgesvd.f"
			i__3 = *m - 1;
#line 2736 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2743 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2743 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2745 "cgesvd.f"
			ie = 1;
#line 2746 "cgesvd.f"
			itauq = itau;
#line 2747 "cgesvd.f"
			itaup = itauq + *m;
#line 2748 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 2756 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2756 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2760 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2768 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2768 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2776 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2776 "cgesvd.f"
			cungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2779 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2787 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2797 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2804 "cgesvd.f"
			clacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 2807 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2811 "cgesvd.f"
			itau = 1;
#line 2812 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2818 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2818 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2820 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2826 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2826 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2828 "cgesvd.f"
			ie = 1;
#line 2829 "cgesvd.f"
			itauq = itau;
#line 2830 "cgesvd.f"
			itaup = itauq + *m;
#line 2831 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 2835 "cgesvd.f"
			i__2 = *m - 1;
#line 2835 "cgesvd.f"
			i__3 = *m - 1;
#line 2835 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2842 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2842 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 2850 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2850 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2858 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2858 "cgesvd.f"
			cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2860 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 2868 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 2872 "cgesvd.f"
		    }

#line 2874 "cgesvd.f"
		} else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

#line 2881 "cgesvd.f"
		    if (*lwork >= *m * *m + *m * 3) {

/*                    Sufficient workspace for a fast algorithm */

#line 2885 "cgesvd.f"
			iu = 1;
#line 2886 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

#line 2890 "cgesvd.f"
			    ldwrku = *lda;
#line 2891 "cgesvd.f"
			} else {

/*                       WORK(IU) is LDA by M */

#line 2895 "cgesvd.f"
			    ldwrku = *m;
#line 2896 "cgesvd.f"
			}
#line 2897 "cgesvd.f"
			itau = iu + ldwrku * *m;
#line 2898 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2904 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2904 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 2909 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 2911 "cgesvd.f"
			i__2 = *m - 1;
#line 2911 "cgesvd.f"
			i__3 = *m - 1;
#line 2911 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2918 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2918 "cgesvd.f"
			cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 2920 "cgesvd.f"
			ie = 1;
#line 2921 "cgesvd.f"
			itauq = itau;
#line 2922 "cgesvd.f"
			itaup = itauq + *m;
#line 2923 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 2929 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2929 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 2933 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 2941 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2941 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2949 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2949 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 2951 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 2959 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 2968 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				a[a_offset], lda, &c_b1, &vt[vt_offset], ldvt,
				 (ftnlen)1, (ftnlen)1);

#line 2971 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 2975 "cgesvd.f"
			itau = 1;
#line 2976 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2982 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2982 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 2984 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 2990 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 2990 "cgesvd.f"
			cunglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 2995 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 2996 "cgesvd.f"
			i__2 = *m - 1;
#line 2996 "cgesvd.f"
			i__3 = *m - 1;
#line 2996 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 2998 "cgesvd.f"
			ie = 1;
#line 2999 "cgesvd.f"
			itauq = itau;
#line 3000 "cgesvd.f"
			itaup = itauq + *m;
#line 3001 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3007 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3007 "cgesvd.f"
			cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3016 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3016 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3024 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3024 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3026 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3034 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3038 "cgesvd.f"
		    }

#line 3040 "cgesvd.f"
		}

#line 3042 "cgesvd.f"
	    } else if (wntva) {

#line 3044 "cgesvd.f"
		if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
#line 3050 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3050 "cgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3054 "cgesvd.f"
			ir = 1;
#line 3055 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

#line 3059 "cgesvd.f"
			    ldwrkr = *lda;
#line 3060 "cgesvd.f"
			} else {

/*                       WORK(IR) is M by M */

#line 3064 "cgesvd.f"
			    ldwrkr = *m;
#line 3065 "cgesvd.f"
			}
#line 3066 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3067 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3073 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3073 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3075 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

#line 3079 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[ir], &
				ldwrkr, (ftnlen)1);
#line 3081 "cgesvd.f"
			i__2 = *m - 1;
#line 3081 "cgesvd.f"
			i__3 = *m - 1;
#line 3081 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ir + 
				ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3088 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3088 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3090 "cgesvd.f"
			ie = 1;
#line 3091 "cgesvd.f"
			itauq = itau;
#line 3092 "cgesvd.f"
			itaup = itauq + *m;
#line 3093 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3099 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3099 "cgesvd.f"
			cgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need   M*M+3*M-1, */
/*                                 prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3109 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3109 "cgesvd.f"
			cungbr_("P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3112 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3119 "cgesvd.f"
			cbdsqr_("U", m, m, &c__0, &c__0, &s[1], &rwork[ie], &
				work[ir], &ldwrkr, cdum, &c__1, cdum, &c__1, &
				rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3128 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[ir], &ldwrkr, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3133 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3135 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3139 "cgesvd.f"
			itau = 1;
#line 3140 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3146 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3146 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3148 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3154 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3154 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3156 "cgesvd.f"
			ie = 1;
#line 3157 "cgesvd.f"
			itauq = itau;
#line 3158 "cgesvd.f"
			itaup = itauq + *m;
#line 3159 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3163 "cgesvd.f"
			i__2 = *m - 1;
#line 3163 "cgesvd.f"
			i__3 = *m - 1;
#line 3163 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3170 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3170 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3179 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3179 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 3182 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3189 "cgesvd.f"
			cbdsqr_("U", m, n, &c__0, &c__0, &s[1], &rwork[ie], &
				vt[vt_offset], ldvt, cdum, &c__1, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

#line 3193 "cgesvd.f"
		    }

#line 3195 "cgesvd.f"
		} else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
#line 3201 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3201 "cgesvd.f"
		    if (*lwork >= (*m << 1) * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3205 "cgesvd.f"
			iu = 1;
#line 3206 "cgesvd.f"
			if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

#line 3210 "cgesvd.f"
			    ldwrku = *lda;
#line 3211 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3212 "cgesvd.f"
			    ldwrkr = *lda;
#line 3213 "cgesvd.f"
			} else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

#line 3217 "cgesvd.f"
			    ldwrku = *lda;
#line 3218 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3219 "cgesvd.f"
			    ldwrkr = *m;
#line 3220 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

#line 3224 "cgesvd.f"
			    ldwrku = *m;
#line 3225 "cgesvd.f"
			    ir = iu + ldwrku * *m;
#line 3226 "cgesvd.f"
			    ldwrkr = *m;
#line 3227 "cgesvd.f"
			}
#line 3228 "cgesvd.f"
			itau = ir + ldwrkr * *m;
#line 3229 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3235 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3235 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3237 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3243 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3243 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3248 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3250 "cgesvd.f"
			i__2 = *m - 1;
#line 3250 "cgesvd.f"
			i__3 = *m - 1;
#line 3250 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3252 "cgesvd.f"
			ie = 1;
#line 3253 "cgesvd.f"
			itauq = itau;
#line 3254 "cgesvd.f"
			itaup = itauq + *m;
#line 3255 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (CWorkspace: need   2*M*M+3*M, */
/*                                 prefer 2*M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need   M) */

#line 3263 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3263 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3267 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &work[ir], &
				ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need   2*M*M+3*M-1, */
/*                                 prefer 2*M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3275 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3275 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3283 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3283 "cgesvd.f"
			cungbr_("Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3286 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need 2*M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3294 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &work[ir], &ldwrkr, cdum, &c__1,
				 &rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3304 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3309 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

#line 3313 "cgesvd.f"
			clacpy_("F", m, m, &work[ir], &ldwrkr, &a[a_offset], 
				lda, (ftnlen)1);

#line 3316 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3320 "cgesvd.f"
			itau = 1;
#line 3321 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3327 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3327 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3329 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3335 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3335 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);
#line 3337 "cgesvd.f"
			ie = 1;
#line 3338 "cgesvd.f"
			itauq = itau;
#line 3339 "cgesvd.f"
			itaup = itauq + *m;
#line 3340 "cgesvd.f"
			iwork = itaup + *m;

/*                    Zero out above L in A */

#line 3344 "cgesvd.f"
			i__2 = *m - 1;
#line 3344 "cgesvd.f"
			i__3 = *m - 1;
#line 3344 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &a[(a_dim1 <<
				 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3351 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3351 "cgesvd.f"
			cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3360 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3360 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &a[a_offset], lda, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3368 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3368 "cgesvd.f"
			cungbr_("Q", m, m, m, &a[a_offset], lda, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3370 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3378 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &a[a_offset], lda, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3382 "cgesvd.f"
		    }

#line 3384 "cgesvd.f"
		} else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
#line 3391 "cgesvd.f"
		    i__2 = *n + *m, i__3 = *m * 3;
#line 3391 "cgesvd.f"
		    if (*lwork >= *m * *m + max(i__2,i__3)) {

/*                    Sufficient workspace for a fast algorithm */

#line 3395 "cgesvd.f"
			iu = 1;
#line 3396 "cgesvd.f"
			if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

#line 3400 "cgesvd.f"
			    ldwrku = *lda;
#line 3401 "cgesvd.f"
			} else {

/*                       WORK(IU) is M by M */

#line 3405 "cgesvd.f"
			    ldwrku = *m;
#line 3406 "cgesvd.f"
			}
#line 3407 "cgesvd.f"
			itau = iu + ldwrku * *m;
#line 3408 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3414 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3414 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3416 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3422 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3422 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

#line 3427 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &work[iu], &
				ldwrku, (ftnlen)1);
#line 3429 "cgesvd.f"
			i__2 = *m - 1;
#line 3429 "cgesvd.f"
			i__3 = *m - 1;
#line 3429 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[iu + 
				ldwrku], &ldwrku, (ftnlen)1);
#line 3431 "cgesvd.f"
			ie = 1;
#line 3432 "cgesvd.f"
			itauq = itau;
#line 3433 "cgesvd.f"
			itaup = itauq + *m;
#line 3434 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3440 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3440 "cgesvd.f"
			cgebrd_(m, m, &work[iu], &ldwrku, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);
#line 3444 "cgesvd.f"
			clacpy_("L", m, m, &work[iu], &ldwrku, &u[u_offset], 
				ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB) */
/*                    (RWorkspace: 0) */

#line 3451 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3451 "cgesvd.f"
			cungbr_("P", m, m, m, &work[iu], &ldwrku, &work[itaup]
				, &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3459 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3459 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3461 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: need BDSPAC) */

#line 3469 "cgesvd.f"
			cbdsqr_("U", m, m, m, &c__0, &s[1], &rwork[ie], &work[
				iu], &ldwrku, &u[u_offset], ldu, cdum, &c__1, 
				&rwork[irwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (CWorkspace: need M*M) */
/*                    (RWorkspace: 0) */

#line 3478 "cgesvd.f"
			cgemm_("N", "N", m, n, m, &c_b2, &work[iu], &ldwrku, &
				vt[vt_offset], ldvt, &c_b1, &a[a_offset], lda,
				 (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

#line 3483 "cgesvd.f"
			clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

#line 3485 "cgesvd.f"
		    } else {

/*                    Insufficient workspace for a fast algorithm */

#line 3489 "cgesvd.f"
			itau = 1;
#line 3490 "cgesvd.f"
			iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (CWorkspace: need 2*M, prefer M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3496 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3496 "cgesvd.f"
			cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
				iwork], &i__2, &ierr);
#line 3498 "cgesvd.f"
			clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], 
				ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (CWorkspace: need M+N, prefer M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3504 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3504 "cgesvd.f"
			cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
				work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

#line 3509 "cgesvd.f"
			clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], 
				ldu, (ftnlen)1);
#line 3510 "cgesvd.f"
			i__2 = *m - 1;
#line 3510 "cgesvd.f"
			i__3 = *m - 1;
#line 3510 "cgesvd.f"
			claset_("U", &i__2, &i__3, &c_b1, &c_b1, &u[(u_dim1 <<
				 1) + 1], ldu, (ftnlen)1);
#line 3512 "cgesvd.f"
			ie = 1;
#line 3513 "cgesvd.f"
			itauq = itau;
#line 3514 "cgesvd.f"
			itaup = itauq + *m;
#line 3515 "cgesvd.f"
			iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB) */
/*                    (RWorkspace: need M) */

#line 3521 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3521 "cgesvd.f"
			cgebrd_(m, m, &u[u_offset], ldu, &s[1], &rwork[ie], &
				work[itauq], &work[itaup], &work[iwork], &
				i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB) */
/*                    (RWorkspace: 0) */

#line 3530 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3530 "cgesvd.f"
			cunmbr_("P", "L", "C", m, n, m, &u[u_offset], ldu, &
				work[itaup], &vt[vt_offset], ldvt, &work[
				iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*                    (RWorkspace: 0) */

#line 3538 "cgesvd.f"
			i__2 = *lwork - iwork + 1;
#line 3538 "cgesvd.f"
			cungbr_("Q", m, m, m, &u[u_offset], ldu, &work[itauq],
				 &work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3540 "cgesvd.f"
			irwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (CWorkspace: 0) */
/*                    (RWorkspace: need BDSPAC) */

#line 3548 "cgesvd.f"
			cbdsqr_("U", m, n, m, &c__0, &s[1], &rwork[ie], &vt[
				vt_offset], ldvt, &u[u_offset], ldu, cdum, &
				c__1, &rwork[irwork], info, (ftnlen)1);

#line 3552 "cgesvd.f"
		    }

#line 3554 "cgesvd.f"
		}

#line 3556 "cgesvd.f"
	    }

#line 3558 "cgesvd.f"
	} else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

#line 3565 "cgesvd.f"
	    ie = 1;
#line 3566 "cgesvd.f"
	    itauq = 1;
#line 3567 "cgesvd.f"
	    itaup = itauq + *m;
#line 3568 "cgesvd.f"
	    iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
/*           (RWorkspace: M) */

#line 3574 "cgesvd.f"
	    i__2 = *lwork - iwork + 1;
#line 3574 "cgesvd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, &ierr);
#line 3577 "cgesvd.f"
	    if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3584 "cgesvd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 3585 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3585 "cgesvd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3587 "cgesvd.f"
	    }
#line 3588 "cgesvd.f"
	    if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB) */
/*              (RWorkspace: 0) */

#line 3595 "cgesvd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 3596 "cgesvd.f"
		if (wntva) {
#line 3596 "cgesvd.f"
		    nrvt = *n;
#line 3596 "cgesvd.f"
		}
#line 3598 "cgesvd.f"
		if (wntvs) {
#line 3598 "cgesvd.f"
		    nrvt = *m;
#line 3598 "cgesvd.f"
		}
#line 3600 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3600 "cgesvd.f"
		cungbr_("P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], 
			&work[iwork], &i__2, &ierr, (ftnlen)1);
#line 3602 "cgesvd.f"
	    }
#line 3603 "cgesvd.f"
	    if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB) */
/*              (RWorkspace: 0) */

#line 3610 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3610 "cgesvd.f"
		cungbr_("Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3612 "cgesvd.f"
	    }
#line 3613 "cgesvd.f"
	    if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*              (RWorkspace: 0) */

#line 3620 "cgesvd.f"
		i__2 = *lwork - iwork + 1;
#line 3620 "cgesvd.f"
		cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			iwork], &i__2, &ierr, (ftnlen)1);
#line 3622 "cgesvd.f"
	    }
#line 3623 "cgesvd.f"
	    irwork = ie + *m;
#line 3624 "cgesvd.f"
	    if (wntuas || wntuo) {
#line 3624 "cgesvd.f"
		nru = *m;
#line 3624 "cgesvd.f"
	    }
#line 3626 "cgesvd.f"
	    if (wntun) {
#line 3626 "cgesvd.f"
		nru = 0;
#line 3626 "cgesvd.f"
	    }
#line 3628 "cgesvd.f"
	    if (wntvas || wntvo) {
#line 3628 "cgesvd.f"
		ncvt = *n;
#line 3628 "cgesvd.f"
	    }
#line 3630 "cgesvd.f"
	    if (wntvn) {
#line 3630 "cgesvd.f"
		ncvt = 0;
#line 3630 "cgesvd.f"
	    }
#line 3632 "cgesvd.f"
	    if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3640 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3643 "cgesvd.f"
	    } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3651 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &a[
			a_offset], lda, &u[u_offset], ldu, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3654 "cgesvd.f"
	    } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (CWorkspace: 0) */
/*              (RWorkspace: need BDSPAC) */

#line 3662 "cgesvd.f"
		cbdsqr_("L", m, &ncvt, &nru, &c__0, &s[1], &rwork[ie], &vt[
			vt_offset], ldvt, &a[a_offset], lda, cdum, &c__1, &
			rwork[irwork], info, (ftnlen)1);
#line 3665 "cgesvd.f"
	    }

#line 3667 "cgesvd.f"
	}

#line 3669 "cgesvd.f"
    }

/*     Undo scaling if necessary */

#line 3673 "cgesvd.f"
    if (iscl == 1) {
#line 3674 "cgesvd.f"
	if (anrm > bignum) {
#line 3674 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3674 "cgesvd.f"
	}
#line 3677 "cgesvd.f"
	if (*info != 0 && anrm > bignum) {
#line 3677 "cgesvd.f"
	    i__2 = minmn - 1;
#line 3677 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3677 "cgesvd.f"
	}
#line 3680 "cgesvd.f"
	if (anrm < smlnum) {
#line 3680 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 3680 "cgesvd.f"
	}
#line 3683 "cgesvd.f"
	if (*info != 0 && anrm < smlnum) {
#line 3683 "cgesvd.f"
	    i__2 = minmn - 1;
#line 3683 "cgesvd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 3683 "cgesvd.f"
	}
#line 3686 "cgesvd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 3690 "cgesvd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 3692 "cgesvd.f"
    return 0;

/*     End of CGESVD */

} /* cgesvd_ */

