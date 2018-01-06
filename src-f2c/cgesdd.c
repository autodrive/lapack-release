#line 1 "cgesdd.f"
/* cgesdd.f -- translated by f2c (version 20100827).
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

#line 1 "cgesdd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c_n1 = -1;
static integer c__0 = 0;
static integer c__1 = 1;

/* > \brief \b CGESDD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESDD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesdd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesdd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesdd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), S( * ) */
/*       COMPLEX            A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESDD computes the singular value decomposition (SVD) of a complex */
/* > M-by-N matrix A, optionally computing the left and/or right singular */
/* > vectors, by using divide-and-conquer method. The SVD is written */
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
/* > Note that the routine returns VT = V**H, not V. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          Specifies options for computing all or part of the matrix U: */
/* >          = 'A':  all M columns of U and all N rows of V**H are */
/* >                  returned in the arrays U and VT; */
/* >          = 'S':  the first min(M,N) columns of U and the first */
/* >                  min(M,N) rows of V**H are returned in the arrays U */
/* >                  and VT; */
/* >          = 'O':  If M >= N, the first N columns of U are overwritten */
/* >                  in the array A and all rows of V**H are returned in */
/* >                  the array VT; */
/* >                  otherwise, all columns of U are returned in the */
/* >                  array U and the first M rows of V**H are overwritten */
/* >                  in the array A; */
/* >          = 'N':  no columns of U or rows of V**H are computed. */
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
/* >          if JOBZ = 'O',  A is overwritten with the first N columns */
/* >                          of U (the left singular vectors, stored */
/* >                          columnwise) if M >= N; */
/* >                          A is overwritten with the first M rows */
/* >                          of V**H (the right singular vectors, stored */
/* >                          rowwise) otherwise. */
/* >          if JOBZ .ne. 'O', the contents of A are destroyed. */
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
/* >          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N; */
/* >          UCOL = min(M,N) if JOBZ = 'S'. */
/* >          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M */
/* >          unitary matrix U; */
/* >          if JOBZ = 'S', U contains the first min(M,N) columns of U */
/* >          (the left singular vectors, stored columnwise); */
/* >          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U.  LDU >= 1; */
/* >          if JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is COMPLEX array, dimension (LDVT,N) */
/* >          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the */
/* >          N-by-N unitary matrix V**H; */
/* >          if JOBZ = 'S', VT contains the first min(M,N) rows of */
/* >          V**H (the right singular vectors, stored rowwise); */
/* >          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >          The leading dimension of the array VT.  LDVT >= 1; */
/* >          if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N; */
/* >          if JOBZ = 'S', LDVT >= min(M,N). */
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
/* >          The dimension of the array WORK. LWORK >= 1. */
/* >          If LWORK = -1, a workspace query is assumed.  The optimal */
/* >          size for the WORK array is calculated and stored in WORK(1), */
/* >          and no other work except argument checking is performed. */
/* > */
/* >          Let mx = max(M,N) and mn = min(M,N). */
/* >          If JOBZ = 'N', LWORK >= 2*mn + mx. */
/* >          If JOBZ = 'O', LWORK >= 2*mn*mn + 2*mn + mx. */
/* >          If JOBZ = 'S', LWORK >=   mn*mn + 3*mn. */
/* >          If JOBZ = 'A', LWORK >=   mn*mn + 2*mn + mx. */
/* >          These are not tight minimums in all cases; see comments inside code. */
/* >          For good performance, LWORK should generally be larger; */
/* >          a query is recommended. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* >          Let mx = max(M,N) and mn = min(M,N). */
/* >          If JOBZ = 'N',    LRWORK >= 5*mn (LAPACK <= 3.6 needs 7*mn); */
/* >          else if mx >> mn, LRWORK >= 5*mn*mn + 5*mn; */
/* >          else              LRWORK >= max( 5*mn*mn + 5*mn, */
/* >                                           2*mx*mn + 2*mn*mn + mn ). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (8*min(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  The updating process of SBDSDC did not converge. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexGEsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgesdd_(char *jobz, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, 
	integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info, 
	ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer lwork_cunglq_mn__, lwork_cunglq_nn__, lwork_cungqr_mm__, 
	    lwork_cungqr_mn__, i__, ie, lwork_cungbr_p_mn__, il, 
	    lwork_cungbr_p_nn__, lwork_cungbr_q_mn__, lwork_cungbr_q_mm__, ir,
	     iu, blk;
    static doublereal dum[1], eps;
    static integer iru, ivt;
    static doublecomplex cdum[1];
    static integer iscl, lwork_cunmbr_prc_mm__, lwork_cunmbr_prc_mn__, 
	    lwork_cunmbr_prc_nn__;
    static doublereal anrm;
    static integer ierr, itau, lwork_cunmbr_qln_mm__, lwork_cunmbr_qln_mn__, 
	    lwork_cunmbr_qln_nn__, idum[1], irvt;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn, wrkbl, itaup, itauq;
    static logical wntqa;
    static integer nwork;
    extern /* Subroutine */ int clacp2_(char *, integer *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, ftnlen);
    static logical wntqn, wntqo, wntqs;
    static integer mnthr1, mnthr2;
    extern /* Subroutine */ int cgebrd_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clacrm_(integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, doublereal *)
	    , clarcm_(integer *, integer *, doublereal *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *), clascl_(char *, integer *, integer *, doublereal *,
	     doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen), sbdsdc_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), cgeqrf_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen), cungbr_(char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), cunmbr_(char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), cunglq_(integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, integer *);
    static integer ldwrkl;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer ldwrkr, minwrk, ldwrku, maxwrk, ldwkvt;
    static doublereal smlnum;
    static logical wntqas, lquery;
    static integer nrwork, lwork_cgebrd_mm__, lwork_cgebrd_mn__, 
	    lwork_cgebrd_nn__, lwork_cgelqf_mn__, lwork_cgeqrf_mn__;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 295 "cgesdd.f"
    /* Parameter adjustments */
#line 295 "cgesdd.f"
    a_dim1 = *lda;
#line 295 "cgesdd.f"
    a_offset = 1 + a_dim1;
#line 295 "cgesdd.f"
    a -= a_offset;
#line 295 "cgesdd.f"
    --s;
#line 295 "cgesdd.f"
    u_dim1 = *ldu;
#line 295 "cgesdd.f"
    u_offset = 1 + u_dim1;
#line 295 "cgesdd.f"
    u -= u_offset;
#line 295 "cgesdd.f"
    vt_dim1 = *ldvt;
#line 295 "cgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 295 "cgesdd.f"
    vt -= vt_offset;
#line 295 "cgesdd.f"
    --work;
#line 295 "cgesdd.f"
    --rwork;
#line 295 "cgesdd.f"
    --iwork;
#line 295 "cgesdd.f"

#line 295 "cgesdd.f"
    /* Function Body */
#line 295 "cgesdd.f"
    *info = 0;
#line 296 "cgesdd.f"
    minmn = min(*m,*n);
#line 297 "cgesdd.f"
    mnthr1 = (integer) (minmn * 17. / 9.);
#line 298 "cgesdd.f"
    mnthr2 = (integer) (minmn * 5. / 3.);
#line 299 "cgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 300 "cgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 301 "cgesdd.f"
    wntqas = wntqa || wntqs;
#line 302 "cgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 303 "cgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 304 "cgesdd.f"
    lquery = *lwork == -1;
#line 305 "cgesdd.f"
    minwrk = 1;
#line 306 "cgesdd.f"
    maxwrk = 1;

#line 308 "cgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 309 "cgesdd.f"
	*info = -1;
#line 310 "cgesdd.f"
    } else if (*m < 0) {
#line 311 "cgesdd.f"
	*info = -2;
#line 312 "cgesdd.f"
    } else if (*n < 0) {
#line 313 "cgesdd.f"
	*info = -3;
#line 314 "cgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 315 "cgesdd.f"
	*info = -5;
#line 316 "cgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 318 "cgesdd.f"
	*info = -8;
#line 319 "cgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 322 "cgesdd.f"
	*info = -10;
#line 323 "cgesdd.f"
    }

/*     Compute workspace */
/*       Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace allocated at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to */
/*       real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 333 "cgesdd.f"
    if (*info == 0) {
#line 334 "cgesdd.f"
	minwrk = 1;
#line 335 "cgesdd.f"
	maxwrk = 1;
#line 336 "cgesdd.f"
	if (*m >= *n && minmn > 0) {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD (sbdsdc) is */
/*           BDSPAC = 3*N*N + 4*N for singular values and vectors; */
/*           BDSPAC = 4*N         for singular values only; */
/*           not including e, RU, and RVT matrices. */

/*           Compute space preferred for each routine */
#line 345 "cgesdd.f"
	    cgebrd_(m, n, cdum, m, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 347 "cgesdd.f"
	    lwork_cgebrd_mn__ = (integer) cdum[0].r;

#line 349 "cgesdd.f"
	    cgebrd_(n, n, cdum, n, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 351 "cgesdd.f"
	    lwork_cgebrd_nn__ = (integer) cdum[0].r;

#line 353 "cgesdd.f"
	    cgeqrf_(m, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 354 "cgesdd.f"
	    lwork_cgeqrf_mn__ = (integer) cdum[0].r;

#line 356 "cgesdd.f"
	    cungbr_("P", n, n, n, cdum, n, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 358 "cgesdd.f"
	    lwork_cungbr_p_nn__ = (integer) cdum[0].r;

#line 360 "cgesdd.f"
	    cungbr_("Q", m, m, n, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 362 "cgesdd.f"
	    lwork_cungbr_q_mm__ = (integer) cdum[0].r;

#line 364 "cgesdd.f"
	    cungbr_("Q", m, n, n, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 366 "cgesdd.f"
	    lwork_cungbr_q_mn__ = (integer) cdum[0].r;

#line 368 "cgesdd.f"
	    cungqr_(m, m, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 370 "cgesdd.f"
	    lwork_cungqr_mm__ = (integer) cdum[0].r;

#line 372 "cgesdd.f"
	    cungqr_(m, n, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 374 "cgesdd.f"
	    lwork_cungqr_mn__ = (integer) cdum[0].r;

#line 376 "cgesdd.f"
	    cunmbr_("P", "R", "C", n, n, n, cdum, n, cdum, cdum, n, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 378 "cgesdd.f"
	    lwork_cunmbr_prc_nn__ = (integer) cdum[0].r;

#line 380 "cgesdd.f"
	    cunmbr_("Q", "L", "N", m, m, n, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 382 "cgesdd.f"
	    lwork_cunmbr_qln_mm__ = (integer) cdum[0].r;

#line 384 "cgesdd.f"
	    cunmbr_("Q", "L", "N", m, n, n, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 386 "cgesdd.f"
	    lwork_cunmbr_qln_mn__ = (integer) cdum[0].r;

#line 388 "cgesdd.f"
	    cunmbr_("Q", "L", "N", n, n, n, cdum, n, cdum, cdum, n, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 390 "cgesdd.f"
	    lwork_cunmbr_qln_nn__ = (integer) cdum[0].r;

#line 392 "cgesdd.f"
	    if (*m >= mnthr1) {
#line 393 "cgesdd.f"
		if (wntqn) {

/*                 Path 1 (M >> N, JOBZ='N') */

#line 397 "cgesdd.f"
		    maxwrk = *n + lwork_cgeqrf_mn__;
/* Computing MAX */
#line 398 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cgebrd_nn__;
#line 398 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 399 "cgesdd.f"
		    minwrk = *n * 3;
#line 400 "cgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M >> N, JOBZ='O') */

#line 404 "cgesdd.f"
		    wrkbl = *n + lwork_cgeqrf_mn__;
/* Computing MAX */
#line 405 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_cungqr_mn__;
#line 405 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 406 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cgebrd_nn__;
#line 406 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 407 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cunmbr_qln_nn__;
#line 407 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 408 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cunmbr_prc_nn__;
#line 408 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 409 "cgesdd.f"
		    maxwrk = *m * *n + *n * *n + wrkbl;
#line 410 "cgesdd.f"
		    minwrk = (*n << 1) * *n + *n * 3;
#line 411 "cgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M >> N, JOBZ='S') */

#line 415 "cgesdd.f"
		    wrkbl = *n + lwork_cgeqrf_mn__;
/* Computing MAX */
#line 416 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_cungqr_mn__;
#line 416 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 417 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cgebrd_nn__;
#line 417 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 418 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cunmbr_qln_nn__;
#line 418 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 419 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cunmbr_prc_nn__;
#line 419 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 420 "cgesdd.f"
		    maxwrk = *n * *n + wrkbl;
#line 421 "cgesdd.f"
		    minwrk = *n * *n + *n * 3;
#line 422 "cgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M >> N, JOBZ='A') */

#line 426 "cgesdd.f"
		    wrkbl = *n + lwork_cgeqrf_mn__;
/* Computing MAX */
#line 427 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_cungqr_mm__;
#line 427 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 428 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cgebrd_nn__;
#line 428 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 429 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cunmbr_qln_nn__;
#line 429 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 430 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_cunmbr_prc_nn__;
#line 430 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 431 "cgesdd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 432 "cgesdd.f"
		    i__1 = *n * 3, i__2 = *n + *m;
#line 432 "cgesdd.f"
		    minwrk = *n * *n + max(i__1,i__2);
#line 433 "cgesdd.f"
		}
#line 434 "cgesdd.f"
	    } else if (*m >= mnthr2) {

/*              Path 5 (M >> N, but not as much as MNTHR1) */

#line 438 "cgesdd.f"
		maxwrk = (*n << 1) + lwork_cgebrd_mn__;
#line 439 "cgesdd.f"
		minwrk = (*n << 1) + *m;
#line 440 "cgesdd.f"
		if (wntqo) {
/*                 Path 5o (M >> N, JOBZ='O') */
/* Computing MAX */
#line 442 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cungbr_p_nn__;
#line 442 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 443 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cungbr_q_mn__;
#line 443 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 444 "cgesdd.f"
		    maxwrk += *m * *n;
#line 445 "cgesdd.f"
		    minwrk += *n * *n;
#line 446 "cgesdd.f"
		} else if (wntqs) {
/*                 Path 5s (M >> N, JOBZ='S') */
/* Computing MAX */
#line 448 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cungbr_p_nn__;
#line 448 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 449 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cungbr_q_mn__;
#line 449 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 450 "cgesdd.f"
		} else if (wntqa) {
/*                 Path 5a (M >> N, JOBZ='A') */
/* Computing MAX */
#line 452 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cungbr_p_nn__;
#line 452 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 453 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cungbr_q_mm__;
#line 453 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 454 "cgesdd.f"
		}
#line 455 "cgesdd.f"
	    } else {

/*              Path 6 (M >= N, but not much larger) */

#line 459 "cgesdd.f"
		maxwrk = (*n << 1) + lwork_cgebrd_mn__;
#line 460 "cgesdd.f"
		minwrk = (*n << 1) + *m;
#line 461 "cgesdd.f"
		if (wntqo) {
/*                 Path 6o (M >= N, JOBZ='O') */
/* Computing MAX */
#line 463 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cunmbr_prc_nn__;
#line 463 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 464 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cunmbr_qln_mn__;
#line 464 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 465 "cgesdd.f"
		    maxwrk += *m * *n;
#line 466 "cgesdd.f"
		    minwrk += *n * *n;
#line 467 "cgesdd.f"
		} else if (wntqs) {
/*                 Path 6s (M >= N, JOBZ='S') */
/* Computing MAX */
#line 469 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cunmbr_qln_mn__;
#line 469 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 470 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cunmbr_prc_nn__;
#line 470 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 471 "cgesdd.f"
		} else if (wntqa) {
/*                 Path 6a (M >= N, JOBZ='A') */
/* Computing MAX */
#line 473 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cunmbr_qln_mm__;
#line 473 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 474 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_cunmbr_prc_nn__;
#line 474 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 475 "cgesdd.f"
		}
#line 476 "cgesdd.f"
	    }
#line 477 "cgesdd.f"
	} else if (minmn > 0) {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD (sbdsdc) is */
/*           BDSPAC = 3*M*M + 4*M for singular values and vectors; */
/*           BDSPAC = 4*M         for singular values only; */
/*           not including e, RU, and RVT matrices. */

/*           Compute space preferred for each routine */
#line 486 "cgesdd.f"
	    cgebrd_(m, n, cdum, m, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 488 "cgesdd.f"
	    lwork_cgebrd_mn__ = (integer) cdum[0].r;

#line 490 "cgesdd.f"
	    cgebrd_(m, m, cdum, m, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 492 "cgesdd.f"
	    lwork_cgebrd_mm__ = (integer) cdum[0].r;

#line 494 "cgesdd.f"
	    cgelqf_(m, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 495 "cgesdd.f"
	    lwork_cgelqf_mn__ = (integer) cdum[0].r;

#line 497 "cgesdd.f"
	    cungbr_("P", m, n, m, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 499 "cgesdd.f"
	    lwork_cungbr_p_mn__ = (integer) cdum[0].r;

#line 501 "cgesdd.f"
	    cungbr_("P", n, n, m, cdum, n, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 503 "cgesdd.f"
	    lwork_cungbr_p_nn__ = (integer) cdum[0].r;

#line 505 "cgesdd.f"
	    cungbr_("Q", m, m, n, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 507 "cgesdd.f"
	    lwork_cungbr_q_mm__ = (integer) cdum[0].r;

#line 509 "cgesdd.f"
	    cunglq_(m, n, m, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 511 "cgesdd.f"
	    lwork_cunglq_mn__ = (integer) cdum[0].r;

#line 513 "cgesdd.f"
	    cunglq_(n, n, m, cdum, n, cdum, cdum, &c_n1, &ierr);
#line 515 "cgesdd.f"
	    lwork_cunglq_nn__ = (integer) cdum[0].r;

#line 517 "cgesdd.f"
	    cunmbr_("P", "R", "C", m, m, m, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 519 "cgesdd.f"
	    lwork_cunmbr_prc_mm__ = (integer) cdum[0].r;

#line 521 "cgesdd.f"
	    cunmbr_("P", "R", "C", m, n, m, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 523 "cgesdd.f"
	    lwork_cunmbr_prc_mn__ = (integer) cdum[0].r;

#line 525 "cgesdd.f"
	    cunmbr_("P", "R", "C", n, n, m, cdum, n, cdum, cdum, n, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 527 "cgesdd.f"
	    lwork_cunmbr_prc_nn__ = (integer) cdum[0].r;

#line 529 "cgesdd.f"
	    cunmbr_("Q", "L", "N", m, m, m, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 531 "cgesdd.f"
	    lwork_cunmbr_qln_mm__ = (integer) cdum[0].r;

#line 533 "cgesdd.f"
	    if (*n >= mnthr1) {
#line 534 "cgesdd.f"
		if (wntqn) {

/*                 Path 1t (N >> M, JOBZ='N') */

#line 538 "cgesdd.f"
		    maxwrk = *m + lwork_cgelqf_mn__;
/* Computing MAX */
#line 539 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cgebrd_mm__;
#line 539 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 540 "cgesdd.f"
		    minwrk = *m * 3;
#line 541 "cgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N >> M, JOBZ='O') */

#line 545 "cgesdd.f"
		    wrkbl = *m + lwork_cgelqf_mn__;
/* Computing MAX */
#line 546 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_cunglq_mn__;
#line 546 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 547 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cgebrd_mm__;
#line 547 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 548 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cunmbr_qln_mm__;
#line 548 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 549 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cunmbr_prc_mm__;
#line 549 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 550 "cgesdd.f"
		    maxwrk = *m * *n + *m * *m + wrkbl;
#line 551 "cgesdd.f"
		    minwrk = (*m << 1) * *m + *m * 3;
#line 552 "cgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N >> M, JOBZ='S') */

#line 556 "cgesdd.f"
		    wrkbl = *m + lwork_cgelqf_mn__;
/* Computing MAX */
#line 557 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_cunglq_mn__;
#line 557 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 558 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cgebrd_mm__;
#line 558 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 559 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cunmbr_qln_mm__;
#line 559 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 560 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cunmbr_prc_mm__;
#line 560 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 561 "cgesdd.f"
		    maxwrk = *m * *m + wrkbl;
#line 562 "cgesdd.f"
		    minwrk = *m * *m + *m * 3;
#line 563 "cgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N >> M, JOBZ='A') */

#line 567 "cgesdd.f"
		    wrkbl = *m + lwork_cgelqf_mn__;
/* Computing MAX */
#line 568 "cgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_cunglq_nn__;
#line 568 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 569 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cgebrd_mm__;
#line 569 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 570 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cunmbr_qln_mm__;
#line 570 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 571 "cgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_cunmbr_prc_mm__;
#line 571 "cgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 572 "cgesdd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 573 "cgesdd.f"
		    i__1 = *m * 3, i__2 = *m + *n;
#line 573 "cgesdd.f"
		    minwrk = *m * *m + max(i__1,i__2);
#line 574 "cgesdd.f"
		}
#line 575 "cgesdd.f"
	    } else if (*n >= mnthr2) {

/*              Path 5t (N >> M, but not as much as MNTHR1) */

#line 579 "cgesdd.f"
		maxwrk = (*m << 1) + lwork_cgebrd_mn__;
#line 580 "cgesdd.f"
		minwrk = (*m << 1) + *n;
#line 581 "cgesdd.f"
		if (wntqo) {
/*                 Path 5to (N >> M, JOBZ='O') */
/* Computing MAX */
#line 583 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cungbr_q_mm__;
#line 583 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 584 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cungbr_p_mn__;
#line 584 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 585 "cgesdd.f"
		    maxwrk += *m * *n;
#line 586 "cgesdd.f"
		    minwrk += *m * *m;
#line 587 "cgesdd.f"
		} else if (wntqs) {
/*                 Path 5ts (N >> M, JOBZ='S') */
/* Computing MAX */
#line 589 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cungbr_q_mm__;
#line 589 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 590 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cungbr_p_mn__;
#line 590 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 591 "cgesdd.f"
		} else if (wntqa) {
/*                 Path 5ta (N >> M, JOBZ='A') */
/* Computing MAX */
#line 593 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cungbr_q_mm__;
#line 593 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 594 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cungbr_p_nn__;
#line 594 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 595 "cgesdd.f"
		}
#line 596 "cgesdd.f"
	    } else {

/*              Path 6t (N > M, but not much larger) */

#line 600 "cgesdd.f"
		maxwrk = (*m << 1) + lwork_cgebrd_mn__;
#line 601 "cgesdd.f"
		minwrk = (*m << 1) + *n;
#line 602 "cgesdd.f"
		if (wntqo) {
/*                 Path 6to (N > M, JOBZ='O') */
/* Computing MAX */
#line 604 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cunmbr_qln_mm__;
#line 604 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 605 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cunmbr_prc_mn__;
#line 605 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 606 "cgesdd.f"
		    maxwrk += *m * *n;
#line 607 "cgesdd.f"
		    minwrk += *m * *m;
#line 608 "cgesdd.f"
		} else if (wntqs) {
/*                 Path 6ts (N > M, JOBZ='S') */
/* Computing MAX */
#line 610 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cunmbr_qln_mm__;
#line 610 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 611 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cunmbr_prc_mn__;
#line 611 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 612 "cgesdd.f"
		} else if (wntqa) {
/*                 Path 6ta (N > M, JOBZ='A') */
/* Computing MAX */
#line 614 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cunmbr_qln_mm__;
#line 614 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 615 "cgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cunmbr_prc_nn__;
#line 615 "cgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 616 "cgesdd.f"
		}
#line 617 "cgesdd.f"
	    }
#line 618 "cgesdd.f"
	}
#line 619 "cgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 620 "cgesdd.f"
    }
#line 621 "cgesdd.f"
    if (*info == 0) {
#line 622 "cgesdd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 623 "cgesdd.f"
	if (*lwork < minwrk && ! lquery) {
#line 624 "cgesdd.f"
	    *info = -12;
#line 625 "cgesdd.f"
	}
#line 626 "cgesdd.f"
    }

#line 628 "cgesdd.f"
    if (*info != 0) {
#line 629 "cgesdd.f"
	i__1 = -(*info);
#line 629 "cgesdd.f"
	xerbla_("CGESDD", &i__1, (ftnlen)6);
#line 630 "cgesdd.f"
	return 0;
#line 631 "cgesdd.f"
    } else if (lquery) {
#line 632 "cgesdd.f"
	return 0;
#line 633 "cgesdd.f"
    }

/*     Quick return if possible */

#line 637 "cgesdd.f"
    if (*m == 0 || *n == 0) {
#line 638 "cgesdd.f"
	return 0;
#line 639 "cgesdd.f"
    }

/*     Get machine constants */

#line 643 "cgesdd.f"
    eps = slamch_("P", (ftnlen)1);
#line 644 "cgesdd.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 645 "cgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 649 "cgesdd.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 650 "cgesdd.f"
    iscl = 0;
#line 651 "cgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 652 "cgesdd.f"
	iscl = 1;
#line 653 "cgesdd.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 654 "cgesdd.f"
    } else if (anrm > bignum) {
#line 655 "cgesdd.f"
	iscl = 1;
#line 656 "cgesdd.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 657 "cgesdd.f"
    }

#line 659 "cgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 665 "cgesdd.f"
	if (*m >= mnthr1) {

#line 667 "cgesdd.f"
	    if (wntqn) {

/*              Path 1 (M >> N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 672 "cgesdd.f"
		itau = 1;
#line 673 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              CWorkspace: need   N [tau] + N    [work] */
/*              CWorkspace: prefer N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 680 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 680 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 685 "cgesdd.f"
		i__1 = *n - 1;
#line 685 "cgesdd.f"
		i__2 = *n - 1;
#line 685 "cgesdd.f"
		claset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 687 "cgesdd.f"
		ie = 1;
#line 688 "cgesdd.f"
		itauq = 1;
#line 689 "cgesdd.f"
		itaup = itauq + *n;
#line 690 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              CWorkspace: need   2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 697 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 697 "cgesdd.f"
		cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 700 "cgesdd.f"
		nrwork = ie + *n;

/*              Perform bidiagonal SVD, compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + BDSPAC */

#line 706 "cgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 709 "cgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M >> N, JOBZ='O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 715 "cgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 719 "cgesdd.f"
		ldwrku = *n;
#line 720 "cgesdd.f"
		ir = iu + ldwrku * *n;
#line 721 "cgesdd.f"
		if (*lwork >= *m * *n + *n * *n + *n * 3) {

/*                 WORK(IR) is M by N */

#line 725 "cgesdd.f"
		    ldwrkr = *m;
#line 726 "cgesdd.f"
		} else {
#line 727 "cgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 728 "cgesdd.f"
		}
#line 729 "cgesdd.f"
		itau = ir + ldwrkr * *n;
#line 730 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 737 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 737 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK( IR ), zeroing out below it */

#line 742 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 743 "cgesdd.f"
		i__1 = *n - 1;
#line 743 "cgesdd.f"
		i__2 = *n - 1;
#line 743 "cgesdd.f"
		claset_("L", &i__1, &i__2, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 751 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 751 "cgesdd.f"
		cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 753 "cgesdd.f"
		ie = 1;
#line 754 "cgesdd.f"
		itauq = itau;
#line 755 "cgesdd.f"
		itaup = itauq + *n;
#line 756 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 763 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 763 "cgesdd.f"
		cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of R in WORK(IRU) and computing right singular vectors */
/*              of R in WORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 773 "cgesdd.f"
		iru = ie + *n;
#line 774 "cgesdd.f"
		irvt = iru + *n * *n;
#line 775 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 776 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of R */
/*              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 786 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 788 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 788 "cgesdd.f"
		cunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by the right singular vectors of R */
/*              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 798 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 799 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 799 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              CWorkspace: need   N*N [U] + N*N [R] */
/*              CWorkspace: prefer N*N [U] + M*N [R] */
/*              RWorkspace: need   0 */

#line 809 "cgesdd.f"
		i__1 = *m;
#line 809 "cgesdd.f"
		i__2 = ldwrkr;
#line 809 "cgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 810 "cgesdd.f"
		    i__3 = *m - i__ + 1;
#line 810 "cgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 811 "cgesdd.f"
		    cgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1], 
			    lda, &work[iu], &ldwrku, &c_b1, &work[ir], &
			    ldwrkr, (ftnlen)1, (ftnlen)1);
#line 814 "cgesdd.f"
		    clacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 816 "cgesdd.f"
/* L10: */
#line 816 "cgesdd.f"
		}

#line 818 "cgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M >> N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 824 "cgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 828 "cgesdd.f"
		ldwrkr = *n;
#line 829 "cgesdd.f"
		itau = ir + ldwrkr * *n;
#line 830 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              CWorkspace: need   N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 837 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 837 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 842 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 843 "cgesdd.f"
		i__2 = *n - 1;
#line 843 "cgesdd.f"
		i__1 = *n - 1;
#line 843 "cgesdd.f"
		claset_("L", &i__2, &i__1, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 851 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 851 "cgesdd.f"
		cungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 853 "cgesdd.f"
		ie = 1;
#line 854 "cgesdd.f"
		itauq = itau;
#line 855 "cgesdd.f"
		itaup = itauq + *n;
#line 856 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 863 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 863 "cgesdd.f"
		cgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 873 "cgesdd.f"
		iru = ie + *n;
#line 874 "cgesdd.f"
		irvt = iru + *n * *n;
#line 875 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 876 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of R */
/*              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 886 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 887 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 887 "cgesdd.f"
		cunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 897 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 898 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 898 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              CWorkspace: need   N*N [R] */
/*              RWorkspace: need   0 */

#line 907 "cgesdd.f"
		clacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 908 "cgesdd.f"
		cgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &work[ir],
			 &ldwrkr, &c_b1, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 911 "cgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M >> N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 917 "cgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 921 "cgesdd.f"
		ldwrku = *n;
#line 922 "cgesdd.f"
		itau = iu + ldwrku * *n;
#line 923 "cgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              CWorkspace: need   N*N [U] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 930 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 930 "cgesdd.f"
		cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 932 "cgesdd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              CWorkspace: need   N*N [U] + N [tau] + M    [work] */
/*              CWorkspace: prefer N*N [U] + N [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 939 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 939 "cgesdd.f"
		cungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out below it */

#line 944 "cgesdd.f"
		i__2 = *n - 1;
#line 944 "cgesdd.f"
		i__1 = *n - 1;
#line 944 "cgesdd.f"
		claset_("L", &i__2, &i__1, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 946 "cgesdd.f"
		ie = 1;
#line 947 "cgesdd.f"
		itauq = itau;
#line 948 "cgesdd.f"
		itaup = itauq + *n;
#line 949 "cgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 956 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 956 "cgesdd.f"
		cgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 959 "cgesdd.f"
		iru = ie + *n;
#line 960 "cgesdd.f"
		irvt = iru + *n * *n;
#line 961 "cgesdd.f"
		nrwork = irvt + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 969 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by left singular vectors of R */
/*              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 979 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 981 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 981 "cgesdd.f"
		cunmbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 991 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 992 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 992 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              CWorkspace: need   N*N [U] */
/*              RWorkspace: need   0 */

#line 1001 "cgesdd.f"
		cgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &work[iu],
			 &ldwrku, &c_b1, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 1006 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 1008 "cgesdd.f"
	    }

#line 1010 "cgesdd.f"
	} else if (*m >= mnthr2) {

/*           MNTHR2 <= M < MNTHR1 */

/*           Path 5 (M >> N, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           CUNGBR and matrix multiplication to compute singular vectors */

#line 1018 "cgesdd.f"
	    ie = 1;
#line 1019 "cgesdd.f"
	    nrwork = ie + *n;
#line 1020 "cgesdd.f"
	    itauq = 1;
#line 1021 "cgesdd.f"
	    itaup = itauq + *n;
#line 1022 "cgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*N [tauq, taup] + M        [work] */
/*           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   N [e] */

#line 1029 "cgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1029 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 1032 "cgesdd.f"
	    if (wntqn) {

/*              Path 5n (M >> N, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + BDSPAC */

#line 1039 "cgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1041 "cgesdd.f"
	    } else if (wntqo) {
#line 1042 "cgesdd.f"
		iu = nwork;
#line 1043 "cgesdd.f"
		iru = nrwork;
#line 1044 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1045 "cgesdd.f"
		nrwork = irvt + *n * *n;

/*              Path 5o (M >> N, JOBZ='O') */
/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1053 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1054 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1054 "cgesdd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1062 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1062 "cgesdd.f"
		cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

#line 1065 "cgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 1069 "cgesdd.f"
		    ldwrku = *m;
#line 1070 "cgesdd.f"
		} else {

/*                 WORK(IU) is LDWRKU by N */

#line 1074 "cgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 1075 "cgesdd.f"
		}
#line 1076 "cgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1084 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in WORK(IU), copying to VT */
/*              CWorkspace: need   2*N [tauq, taup] + N*N [U] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork] */

#line 1093 "cgesdd.f"
		clarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &work[iu]
			, &ldwrku, &rwork[nrwork]);
#line 1095 "cgesdd.f"
		clacpy_("F", n, n, &work[iu], &ldwrku, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in A by real matrix RWORK(IRU), storing the */
/*              result in WORK(IU), copying to A */
/*              CWorkspace: need   2*N [tauq, taup] + N*N [U] */
/*              CWorkspace: prefer 2*N [tauq, taup] + M*N [U] */
/*              RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork] */
/*              RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1104 "cgesdd.f"
		nrwork = irvt;
#line 1105 "cgesdd.f"
		i__2 = *m;
#line 1105 "cgesdd.f"
		i__1 = ldwrku;
#line 1105 "cgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1106 "cgesdd.f"
		    i__3 = *m - i__ + 1;
#line 1106 "cgesdd.f"
		    chunk = min(i__3,ldwrku);
#line 1107 "cgesdd.f"
		    clacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru], n, 
			    &work[iu], &ldwrku, &rwork[nrwork]);
#line 1109 "cgesdd.f"
		    clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 1111 "cgesdd.f"
/* L20: */
#line 1111 "cgesdd.f"
		}

#line 1113 "cgesdd.f"
	    } else if (wntqs) {

/*              Path 5s (M >> N, JOBZ='S') */
/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1121 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1122 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1122 "cgesdd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1130 "cgesdd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1131 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1131 "cgesdd.f"
		cungbr_("Q", m, n, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1140 "cgesdd.f"
		iru = nrwork;
#line 1141 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1142 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1143 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork] */

#line 1152 "cgesdd.f"
		clarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1154 "cgesdd.f"
		clacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1161 "cgesdd.f"
		nrwork = irvt;
#line 1162 "cgesdd.f"
		clacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1164 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1165 "cgesdd.f"
	    } else {

/*              Path 5a (M >> N, JOBZ='A') */
/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1173 "cgesdd.f"
		clacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1174 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1174 "cgesdd.f"
		cungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*N [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1182 "cgesdd.f"
		clacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1183 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1183 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1192 "cgesdd.f"
		iru = nrwork;
#line 1193 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1194 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1195 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork] */

#line 1204 "cgesdd.f"
		clarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1206 "cgesdd.f"
		clacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1213 "cgesdd.f"
		nrwork = irvt;
#line 1214 "cgesdd.f"
		clacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1216 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1217 "cgesdd.f"
	    }

#line 1219 "cgesdd.f"
	} else {

/*           M .LT. MNTHR2 */

/*           Path 6 (M >= N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */
/*           Use CUNMBR to compute singular vectors */

#line 1227 "cgesdd.f"
	    ie = 1;
#line 1228 "cgesdd.f"
	    nrwork = ie + *n;
#line 1229 "cgesdd.f"
	    itauq = 1;
#line 1230 "cgesdd.f"
	    itaup = itauq + *n;
#line 1231 "cgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*N [tauq, taup] + M        [work] */
/*           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   N [e] */

#line 1238 "cgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1238 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);
#line 1241 "cgesdd.f"
	    if (wntqn) {

/*              Path 6n (M >= N, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + BDSPAC */

#line 1248 "cgesdd.f"
		sbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1250 "cgesdd.f"
	    } else if (wntqo) {
#line 1251 "cgesdd.f"
		iu = nwork;
#line 1252 "cgesdd.f"
		iru = nrwork;
#line 1253 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1254 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1255 "cgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 1259 "cgesdd.f"
		    ldwrku = *m;
#line 1260 "cgesdd.f"
		} else {

/*                 WORK( IU ) is LDWRKU by N */

#line 1264 "cgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 1265 "cgesdd.f"
		}
#line 1266 "cgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Path 6o (M >= N, JOBZ='O') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1275 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1285 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1286 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1286 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1290 "cgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 Path 6o-fast */
/*                 Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*                 Overwrite WORK(IU) by left singular vectors of A, copying */
/*                 to A */
/*                 CWorkspace: need   2*N [tauq, taup] + M*N [U] + N    [work] */
/*                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U] + N*NB [work] */
/*                 RWorkspace: need   N [e] + N*N [RU] */

#line 1300 "cgesdd.f"
		    claset_("F", m, n, &c_b1, &c_b1, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1302 "cgesdd.f"
		    clacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1304 "cgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1304 "cgesdd.f"
		    cunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1307 "cgesdd.f"
		    clacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 1308 "cgesdd.f"
		} else {

/*                 Path 6o-slow */
/*                 Generate Q in A */
/*                 CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work] */
/*                 CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work] */
/*                 RWorkspace: need   0 */

#line 1316 "cgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1316 "cgesdd.f"
		    cungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__1, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 CWorkspace: need   2*N [tauq, taup] + N*N [U] */
/*                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U] */
/*                 RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork] */
/*                 RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1326 "cgesdd.f"
		    nrwork = irvt;
#line 1327 "cgesdd.f"
		    i__1 = *m;
#line 1327 "cgesdd.f"
		    i__2 = ldwrku;
#line 1327 "cgesdd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ +=
			     i__2) {
/* Computing MIN */
#line 1328 "cgesdd.f"
			i__3 = *m - i__ + 1;
#line 1328 "cgesdd.f"
			chunk = min(i__3,ldwrku);
#line 1329 "cgesdd.f"
			clacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru],
				 n, &work[iu], &ldwrku, &rwork[nrwork]);
#line 1332 "cgesdd.f"
			clacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 1334 "cgesdd.f"
/* L30: */
#line 1334 "cgesdd.f"
		    }
#line 1335 "cgesdd.f"
		}

#line 1337 "cgesdd.f"
	    } else if (wntqs) {

/*              Path 6s (M >= N, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1346 "cgesdd.f"
		iru = nrwork;
#line 1347 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1348 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1349 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1359 "cgesdd.f"
		claset_("F", m, n, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1360 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1361 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1361 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1371 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1372 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1372 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1375 "cgesdd.f"
	    } else {

/*              Path 6a (M >= N, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1384 "cgesdd.f"
		iru = nrwork;
#line 1385 "cgesdd.f"
		irvt = iru + *n * *n;
#line 1386 "cgesdd.f"
		nrwork = irvt + *n * *n;
#line 1387 "cgesdd.f"
		sbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 1393 "cgesdd.f"
		claset_("F", m, m, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1394 "cgesdd.f"
		if (*m > *n) {
#line 1395 "cgesdd.f"
		    i__2 = *m - *n;
#line 1395 "cgesdd.f"
		    i__1 = *m - *n;
#line 1395 "cgesdd.f"
		    claset_("F", &i__2, &i__1, &c_b1, &c_b2, &u[*n + 1 + (*n 
			    + 1) * u_dim1], ldu, (ftnlen)1);
#line 1397 "cgesdd.f"
		}

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1405 "cgesdd.f"
		clacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1406 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1406 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1416 "cgesdd.f"
		clacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1417 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1417 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1420 "cgesdd.f"
	    }

#line 1422 "cgesdd.f"
	}

#line 1424 "cgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 1430 "cgesdd.f"
	if (*n >= mnthr1) {

#line 1432 "cgesdd.f"
	    if (wntqn) {

/*              Path 1t (N >> M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 1437 "cgesdd.f"
		itau = 1;
#line 1438 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              CWorkspace: need   M [tau] + M    [work] */
/*              CWorkspace: prefer M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1445 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1445 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 1450 "cgesdd.f"
		i__2 = *m - 1;
#line 1450 "cgesdd.f"
		i__1 = *m - 1;
#line 1450 "cgesdd.f"
		claset_("U", &i__2, &i__1, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1452 "cgesdd.f"
		ie = 1;
#line 1453 "cgesdd.f"
		itauq = 1;
#line 1454 "cgesdd.f"
		itaup = itauq + *m;
#line 1455 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              CWorkspace: need   2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1462 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1462 "cgesdd.f"
		cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 1465 "cgesdd.f"
		nrwork = ie + *m;

/*              Perform bidiagonal SVD, compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + BDSPAC */

#line 1471 "cgesdd.f"
		sbdsdc_("U", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 1474 "cgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N >> M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1480 "cgesdd.f"
		ivt = 1;
#line 1481 "cgesdd.f"
		ldwkvt = *m;

/*              WORK(IVT) is M by M */

#line 1485 "cgesdd.f"
		il = ivt + ldwkvt * *m;
#line 1486 "cgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3) {

/*                 WORK(IL) M by N */

#line 1490 "cgesdd.f"
		    ldwrkl = *m;
#line 1491 "cgesdd.f"
		    chunk = *n;
#line 1492 "cgesdd.f"
		} else {

/*                 WORK(IL) is M by CHUNK */

#line 1496 "cgesdd.f"
		    ldwrkl = *m;
#line 1497 "cgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1498 "cgesdd.f"
		}
#line 1499 "cgesdd.f"
		itau = il + ldwrkl * chunk;
#line 1500 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1507 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1507 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1512 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1513 "cgesdd.f"
		i__2 = *m - 1;
#line 1513 "cgesdd.f"
		i__1 = *m - 1;
#line 1513 "cgesdd.f"
		claset_("U", &i__2, &i__1, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1521 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1521 "cgesdd.f"
		cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1523 "cgesdd.f"
		ie = 1;
#line 1524 "cgesdd.f"
		itauq = itau;
#line 1525 "cgesdd.f"
		itaup = itauq + *m;
#line 1526 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1533 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1533 "cgesdd.f"
		cgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC */

#line 1543 "cgesdd.f"
		iru = ie + *m;
#line 1544 "cgesdd.f"
		irvt = iru + *m * *m;
#line 1545 "cgesdd.f"
		nrwork = irvt + *m * *m;
#line 1546 "cgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of L */
/*              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1556 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1557 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1557 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by the right singular vectors of L */
/*              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1567 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1569 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1569 "cgesdd.f"
		cunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              CWorkspace: need   M*M [VT] + M*M [L] */
/*              CWorkspace: prefer M*M [VT] + M*N [L] */
/*              RWorkspace: need   0 */

#line 1579 "cgesdd.f"
		i__2 = *n;
#line 1579 "cgesdd.f"
		i__1 = chunk;
#line 1579 "cgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1580 "cgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1580 "cgesdd.f"
		    blk = min(i__3,chunk);
#line 1581 "cgesdd.f"
		    cgemm_("N", "N", m, &blk, m, &c_b2, &work[ivt], m, &a[i__ 
			    * a_dim1 + 1], lda, &c_b1, &work[il], &ldwrkl, (
			    ftnlen)1, (ftnlen)1);
#line 1584 "cgesdd.f"
		    clacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1586 "cgesdd.f"
/* L40: */
#line 1586 "cgesdd.f"
		}

#line 1588 "cgesdd.f"
	    } else if (wntqs) {

/*              Path 3t (N >> M, JOBZ='S') */
/*              M right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1594 "cgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1598 "cgesdd.f"
		ldwrkl = *m;
#line 1599 "cgesdd.f"
		itau = il + ldwrkl * *m;
#line 1600 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              CWorkspace: need   M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1607 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1607 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1612 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1613 "cgesdd.f"
		i__1 = *m - 1;
#line 1613 "cgesdd.f"
		i__2 = *m - 1;
#line 1613 "cgesdd.f"
		claset_("U", &i__1, &i__2, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1621 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1621 "cgesdd.f"
		cunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1623 "cgesdd.f"
		ie = 1;
#line 1624 "cgesdd.f"
		itauq = itau;
#line 1625 "cgesdd.f"
		itaup = itauq + *m;
#line 1626 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1633 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1633 "cgesdd.f"
		cgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC */

#line 1643 "cgesdd.f"
		iru = ie + *m;
#line 1644 "cgesdd.f"
		irvt = iru + *m * *m;
#line 1645 "cgesdd.f"
		nrwork = irvt + *m * *m;
#line 1646 "cgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1656 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1657 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1657 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by left singular vectors of L */
/*              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1667 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1668 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1668 "cgesdd.f"
		cunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy VT to WORK(IL), multiply right singular vectors of L */
/*              in WORK(IL) by Q in A, storing result in VT */
/*              CWorkspace: need   M*M [L] */
/*              RWorkspace: need   0 */

#line 1677 "cgesdd.f"
		clacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1678 "cgesdd.f"
		cgemm_("N", "N", m, n, m, &c_b2, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b1, &vt[vt_offset], ldvt, (ftnlen)
			1, (ftnlen)1);

#line 1681 "cgesdd.f"
	    } else if (wntqa) {

/*              Path 4t (N >> M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1687 "cgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1691 "cgesdd.f"
		ldwkvt = *m;
#line 1692 "cgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1693 "cgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              CWorkspace: need   M*M [VT] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1700 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1700 "cgesdd.f"
		cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
#line 1702 "cgesdd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              CWorkspace: need   M*M [VT] + M [tau] + N    [work] */
/*              CWorkspace: prefer M*M [VT] + M [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1709 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1709 "cgesdd.f"
		cunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

/*              Produce L in A, zeroing out above it */

#line 1714 "cgesdd.f"
		i__1 = *m - 1;
#line 1714 "cgesdd.f"
		i__2 = *m - 1;
#line 1714 "cgesdd.f"
		claset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1716 "cgesdd.f"
		ie = 1;
#line 1717 "cgesdd.f"
		itauq = itau;
#line 1718 "cgesdd.f"
		itaup = itauq + *m;
#line 1719 "cgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1726 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1726 "cgesdd.f"
		cgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC */

#line 1736 "cgesdd.f"
		iru = ie + *m;
#line 1737 "cgesdd.f"
		irvt = iru + *m * *m;
#line 1738 "cgesdd.f"
		nrwork = irvt + *m * *m;
#line 1739 "cgesdd.f"
		sbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1749 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1750 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1750 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by right singular vectors of L */
/*              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1760 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1762 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1762 "cgesdd.f"
		cunmbr_("P", "R", "C", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              CWorkspace: need   M*M [VT] */
/*              RWorkspace: need   0 */

#line 1771 "cgesdd.f"
		cgemm_("N", "N", m, n, m, &c_b2, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b1, &a[a_offset], lda, (ftnlen)1,
			 (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1776 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1778 "cgesdd.f"
	    }

#line 1780 "cgesdd.f"
	} else if (*n >= mnthr2) {

/*           MNTHR2 <= N < MNTHR1 */

/*           Path 5t (N >> M, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           CUNGBR and matrix multiplication to compute singular vectors */

#line 1788 "cgesdd.f"
	    ie = 1;
#line 1789 "cgesdd.f"
	    nrwork = ie + *m;
#line 1790 "cgesdd.f"
	    itauq = 1;
#line 1791 "cgesdd.f"
	    itaup = itauq + *m;
#line 1792 "cgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*M [tauq, taup] + N        [work] */
/*           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   M [e] */

#line 1799 "cgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1799 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);

#line 1803 "cgesdd.f"
	    if (wntqn) {

/*              Path 5tn (N >> M, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + BDSPAC */

#line 1810 "cgesdd.f"
		sbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1812 "cgesdd.f"
	    } else if (wntqo) {
#line 1813 "cgesdd.f"
		irvt = nrwork;
#line 1814 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1815 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1816 "cgesdd.f"
		ivt = nwork;

/*              Path 5to (N >> M, JOBZ='O') */
/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1824 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1825 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1825 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Generate P**H in A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1833 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1833 "cgesdd.f"
		cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

#line 1836 "cgesdd.f"
		ldwkvt = *m;
#line 1837 "cgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 1841 "cgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1842 "cgesdd.f"
		    chunk = *n;
#line 1843 "cgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 1847 "cgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 1848 "cgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 1849 "cgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 1857 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRVT) */
/*              storing the result in WORK(IVT), copying to U */
/*              CWorkspace: need   2*M [tauq, taup] + M*M [VT] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork] */

#line 1866 "cgesdd.f"
		clacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &work[ivt], &
			ldwkvt, &rwork[nrwork]);
#line 1868 "cgesdd.f"
		clacpy_("F", m, m, &work[ivt], &ldwkvt, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply RWORK(IRVT) by P**H in A, storing the */
/*              result in WORK(IVT), copying to A */
/*              CWorkspace: need   2*M [tauq, taup] + M*M [VT] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] */
/*              RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork] */
/*              RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 1877 "cgesdd.f"
		nrwork = iru;
#line 1878 "cgesdd.f"
		i__1 = *n;
#line 1878 "cgesdd.f"
		i__2 = chunk;
#line 1878 "cgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1879 "cgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1879 "cgesdd.f"
		    blk = min(i__3,chunk);
#line 1880 "cgesdd.f"
		    clarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1], 
			    lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 1882 "cgesdd.f"
		    clacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
			    a_dim1 + 1], lda, (ftnlen)1);
#line 1884 "cgesdd.f"
/* L50: */
#line 1884 "cgesdd.f"
		}
#line 1885 "cgesdd.f"
	    } else if (wntqs) {

/*              Path 5ts (N >> M, JOBZ='S') */
/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1893 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1894 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1894 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1902 "cgesdd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1903 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1903 "cgesdd.f"
		cungbr_("P", m, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 1912 "cgesdd.f"
		irvt = nrwork;
#line 1913 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1914 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1915 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork] */

#line 1924 "cgesdd.f"
		clacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1926 "cgesdd.f"
		clacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 1933 "cgesdd.f"
		nrwork = iru;
#line 1934 "cgesdd.f"
		clarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1936 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1937 "cgesdd.f"
	    } else {

/*              Path 5ta (N >> M, JOBZ='A') */
/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1945 "cgesdd.f"
		clacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1946 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1946 "cgesdd.f"
		cungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*M [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1954 "cgesdd.f"
		clacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1955 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1955 "cgesdd.f"
		cungbr_("P", n, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 1964 "cgesdd.f"
		irvt = nrwork;
#line 1965 "cgesdd.f"
		iru = irvt + *m * *m;
#line 1966 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 1967 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork] */

#line 1976 "cgesdd.f"
		clacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1978 "cgesdd.f"
		clacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 1985 "cgesdd.f"
		nrwork = iru;
#line 1986 "cgesdd.f"
		clarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1988 "cgesdd.f"
		clacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1989 "cgesdd.f"
	    }

#line 1991 "cgesdd.f"
	} else {

/*           N .LT. MNTHR2 */

/*           Path 6t (N > M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           Use CUNMBR to compute singular vectors */

#line 1999 "cgesdd.f"
	    ie = 1;
#line 2000 "cgesdd.f"
	    nrwork = ie + *m;
#line 2001 "cgesdd.f"
	    itauq = 1;
#line 2002 "cgesdd.f"
	    itaup = itauq + *m;
#line 2003 "cgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*M [tauq, taup] + N        [work] */
/*           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   M [e] */

#line 2010 "cgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 2010 "cgesdd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 2013 "cgesdd.f"
	    if (wntqn) {

/*              Path 6tn (N > M, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + BDSPAC */

#line 2020 "cgesdd.f"
		sbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 2022 "cgesdd.f"
	    } else if (wntqo) {
/*              Path 6to (N > M, JOBZ='O') */
#line 2024 "cgesdd.f"
		ldwkvt = *m;
#line 2025 "cgesdd.f"
		ivt = nwork;
#line 2026 "cgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 2030 "cgesdd.f"
		    claset_("F", m, n, &c_b1, &c_b1, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 2032 "cgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 2033 "cgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 2037 "cgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 2038 "cgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 2039 "cgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 2047 "cgesdd.f"
		irvt = nrwork;
#line 2048 "cgesdd.f"
		iru = irvt + *m * *m;
#line 2049 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 2050 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] */

#line 2060 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2061 "cgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 2061 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 2065 "cgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 Path 6to-fast */
/*                 Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*                 Overwrite WORK(IVT) by right singular vectors of A, */
/*                 copying to A */
/*                 CWorkspace: need   2*M [tauq, taup] + M*N [VT] + M    [work] */
/*                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] + M*NB [work] */
/*                 RWorkspace: need   M [e] + M*M [RVT] */

#line 2075 "cgesdd.f"
		    clacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 2077 "cgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 2077 "cgesdd.f"
		    cunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2080 "cgesdd.f"
		    clacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 2081 "cgesdd.f"
		} else {

/*                 Path 6to-slow */
/*                 Generate P**H in A */
/*                 CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work] */
/*                 CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work] */
/*                 RWorkspace: need   0 */

#line 2089 "cgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 2089 "cgesdd.f"
		    cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 CWorkspace: need   2*M [tauq, taup] + M*M [VT] */
/*                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] */
/*                 RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork] */
/*                 RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 2099 "cgesdd.f"
		    nrwork = iru;
#line 2100 "cgesdd.f"
		    i__2 = *n;
#line 2100 "cgesdd.f"
		    i__1 = chunk;
#line 2100 "cgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 2101 "cgesdd.f"
			i__3 = *n - i__ + 1;
#line 2101 "cgesdd.f"
			blk = min(i__3,chunk);
#line 2102 "cgesdd.f"
			clarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1]
				, lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 2105 "cgesdd.f"
			clacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2107 "cgesdd.f"
/* L60: */
#line 2107 "cgesdd.f"
		    }
#line 2108 "cgesdd.f"
		}
#line 2109 "cgesdd.f"
	    } else if (wntqs) {

/*              Path 6ts (N > M, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 2118 "cgesdd.f"
		irvt = nrwork;
#line 2119 "cgesdd.f"
		iru = irvt + *m * *m;
#line 2120 "cgesdd.f"
		nrwork = iru + *m * *m;
#line 2121 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] */

#line 2131 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2132 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2132 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] */

#line 2142 "cgesdd.f"
		claset_("F", m, n, &c_b1, &c_b1, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2143 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2144 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2144 "cgesdd.f"
		cunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2147 "cgesdd.f"
	    } else {

/*              Path 6ta (N > M, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 2156 "cgesdd.f"
		irvt = nrwork;
#line 2157 "cgesdd.f"
		iru = irvt + *m * *m;
#line 2158 "cgesdd.f"
		nrwork = iru + *m * *m;

#line 2160 "cgesdd.f"
		sbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] */

#line 2170 "cgesdd.f"
		clacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2171 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2171 "cgesdd.f"
		cunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Set all of VT to identity matrix */

#line 2177 "cgesdd.f"
		claset_("F", n, n, &c_b1, &c_b2, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] */

#line 2185 "cgesdd.f"
		clacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2186 "cgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2186 "cgesdd.f"
		cunmbr_("P", "R", "C", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2189 "cgesdd.f"
	    }

#line 2191 "cgesdd.f"
	}

#line 2193 "cgesdd.f"
    }

/*     Undo scaling if necessary */

#line 2197 "cgesdd.f"
    if (iscl == 1) {
#line 2198 "cgesdd.f"
	if (anrm > bignum) {
#line 2198 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2198 "cgesdd.f"
	}
#line 2201 "cgesdd.f"
	if (*info != 0 && anrm > bignum) {
#line 2201 "cgesdd.f"
	    i__1 = minmn - 1;
#line 2201 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2201 "cgesdd.f"
	}
#line 2204 "cgesdd.f"
	if (anrm < smlnum) {
#line 2204 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2204 "cgesdd.f"
	}
#line 2207 "cgesdd.f"
	if (*info != 0 && anrm < smlnum) {
#line 2207 "cgesdd.f"
	    i__1 = minmn - 1;
#line 2207 "cgesdd.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2207 "cgesdd.f"
	}
#line 2210 "cgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 2214 "cgesdd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 2216 "cgesdd.f"
    return 0;

/*     End of CGESDD */

} /* cgesdd_ */

