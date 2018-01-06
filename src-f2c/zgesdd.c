#line 1 "zgesdd.f"
/* zgesdd.f -- translated by f2c (version 20100827).
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

#line 1 "zgesdd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c_n1 = -1;
static integer c__0 = 0;
static integer c__1 = 1;

/* > \brief \b ZGESDD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGESDD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesdd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesdd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesdd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGESDD computes the singular value decomposition (SVD) of a complex */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX*16 array, dimension (LDU,UCOL) */
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
/* >          VT is COMPLEX*16 array, dimension (LDVT,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
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
/* >          > 0:  The updating process of DBDSDC did not converge. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complex16GEsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgesdd_(char *jobz, integer *m, integer *n, 
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
    static integer lwork_zgebrd_mm__, lwork_zgebrd_mn__, lwork_zgebrd_nn__, 
	    lwork_zgelqf_mn__, lwork_zgeqrf_mn__, lwork_zunglq_mn__, 
	    lwork_zunglq_nn__, lwork_zungqr_mm__, lwork_zungqr_mn__, i__, ie, 
	    il, ir, iu, lwork_zungbr_p_mn__, lwork_zungbr_p_nn__, 
	    lwork_zungbr_q_mn__, lwork_zungbr_q_mm__, blk;
    static doublereal dum[1], eps;
    static integer iru, ivt;
    static doublecomplex cdum[1];
    static integer iscl;
    static doublereal anrm;
    static integer idum[1], ierr, itau, irvt, lwork_zunmbr_prc_mm__, 
	    lwork_zunmbr_prc_mn__, lwork_zunmbr_prc_nn__, 
	    lwork_zunmbr_qln_mm__, lwork_zunmbr_qln_mn__, 
	    lwork_zunmbr_qln_nn__;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk, minmn;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer wrkbl, itaup, itauq;
    static logical wntqa;
    static integer nwork;
    static logical wntqn, wntqo, wntqs;
    extern /* Subroutine */ int zlacp2_(char *, integer *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, ftnlen);
    static integer mnthr1, mnthr2;
    extern /* Subroutine */ int dbdsdc_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen),
	     zgebrd_(integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *);
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zlacrm_(integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, doublereal *)
	    , zlarcm_(integer *, integer *, doublereal *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *), zlascl_(char *, integer *, integer *, doublereal *,
	     doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen), zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static integer ldwrkl;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer ldwrkr, minwrk, ldwrku, maxwrk;
    extern /* Subroutine */ int zungbr_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *, ftnlen);
    static integer ldwkvt;
    static doublereal smlnum;
    static logical wntqas;
    extern /* Subroutine */ int zunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen), zunglq_(integer *, integer *, integer *
	    , doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static logical lquery;
    static integer nrwork;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);


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

#line 295 "zgesdd.f"
    /* Parameter adjustments */
#line 295 "zgesdd.f"
    a_dim1 = *lda;
#line 295 "zgesdd.f"
    a_offset = 1 + a_dim1;
#line 295 "zgesdd.f"
    a -= a_offset;
#line 295 "zgesdd.f"
    --s;
#line 295 "zgesdd.f"
    u_dim1 = *ldu;
#line 295 "zgesdd.f"
    u_offset = 1 + u_dim1;
#line 295 "zgesdd.f"
    u -= u_offset;
#line 295 "zgesdd.f"
    vt_dim1 = *ldvt;
#line 295 "zgesdd.f"
    vt_offset = 1 + vt_dim1;
#line 295 "zgesdd.f"
    vt -= vt_offset;
#line 295 "zgesdd.f"
    --work;
#line 295 "zgesdd.f"
    --rwork;
#line 295 "zgesdd.f"
    --iwork;
#line 295 "zgesdd.f"

#line 295 "zgesdd.f"
    /* Function Body */
#line 295 "zgesdd.f"
    *info = 0;
#line 296 "zgesdd.f"
    minmn = min(*m,*n);
#line 297 "zgesdd.f"
    mnthr1 = (integer) (minmn * 17. / 9.);
#line 298 "zgesdd.f"
    mnthr2 = (integer) (minmn * 5. / 3.);
#line 299 "zgesdd.f"
    wntqa = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
#line 300 "zgesdd.f"
    wntqs = lsame_(jobz, "S", (ftnlen)1, (ftnlen)1);
#line 301 "zgesdd.f"
    wntqas = wntqa || wntqs;
#line 302 "zgesdd.f"
    wntqo = lsame_(jobz, "O", (ftnlen)1, (ftnlen)1);
#line 303 "zgesdd.f"
    wntqn = lsame_(jobz, "N", (ftnlen)1, (ftnlen)1);
#line 304 "zgesdd.f"
    lquery = *lwork == -1;
#line 305 "zgesdd.f"
    minwrk = 1;
#line 306 "zgesdd.f"
    maxwrk = 1;

#line 308 "zgesdd.f"
    if (! (wntqa || wntqs || wntqo || wntqn)) {
#line 309 "zgesdd.f"
	*info = -1;
#line 310 "zgesdd.f"
    } else if (*m < 0) {
#line 311 "zgesdd.f"
	*info = -2;
#line 312 "zgesdd.f"
    } else if (*n < 0) {
#line 313 "zgesdd.f"
	*info = -3;
#line 314 "zgesdd.f"
    } else if (*lda < max(1,*m)) {
#line 315 "zgesdd.f"
	*info = -5;
#line 316 "zgesdd.f"
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *
	    m) {
#line 318 "zgesdd.f"
	*info = -8;
#line 319 "zgesdd.f"
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn || 
	    wntqo && *m >= *n && *ldvt < *n) {
#line 322 "zgesdd.f"
	*info = -10;
#line 323 "zgesdd.f"
    }

/*     Compute workspace */
/*       Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace allocated at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to */
/*       real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 333 "zgesdd.f"
    if (*info == 0) {
#line 334 "zgesdd.f"
	minwrk = 1;
#line 335 "zgesdd.f"
	maxwrk = 1;
#line 336 "zgesdd.f"
	if (*m >= *n && minmn > 0) {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD (dbdsdc) is */
/*           BDSPAC = 3*N*N + 4*N for singular values and vectors; */
/*           BDSPAC = 4*N         for singular values only; */
/*           not including e, RU, and RVT matrices. */

/*           Compute space preferred for each routine */
#line 345 "zgesdd.f"
	    zgebrd_(m, n, cdum, m, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 347 "zgesdd.f"
	    lwork_zgebrd_mn__ = (integer) cdum[0].r;

#line 349 "zgesdd.f"
	    zgebrd_(n, n, cdum, n, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 351 "zgesdd.f"
	    lwork_zgebrd_nn__ = (integer) cdum[0].r;

#line 353 "zgesdd.f"
	    zgeqrf_(m, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 354 "zgesdd.f"
	    lwork_zgeqrf_mn__ = (integer) cdum[0].r;

#line 356 "zgesdd.f"
	    zungbr_("P", n, n, n, cdum, n, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 358 "zgesdd.f"
	    lwork_zungbr_p_nn__ = (integer) cdum[0].r;

#line 360 "zgesdd.f"
	    zungbr_("Q", m, m, n, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 362 "zgesdd.f"
	    lwork_zungbr_q_mm__ = (integer) cdum[0].r;

#line 364 "zgesdd.f"
	    zungbr_("Q", m, n, n, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 366 "zgesdd.f"
	    lwork_zungbr_q_mn__ = (integer) cdum[0].r;

#line 368 "zgesdd.f"
	    zungqr_(m, m, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 370 "zgesdd.f"
	    lwork_zungqr_mm__ = (integer) cdum[0].r;

#line 372 "zgesdd.f"
	    zungqr_(m, n, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 374 "zgesdd.f"
	    lwork_zungqr_mn__ = (integer) cdum[0].r;

#line 376 "zgesdd.f"
	    zunmbr_("P", "R", "C", n, n, n, cdum, n, cdum, cdum, n, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 378 "zgesdd.f"
	    lwork_zunmbr_prc_nn__ = (integer) cdum[0].r;

#line 380 "zgesdd.f"
	    zunmbr_("Q", "L", "N", m, m, n, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 382 "zgesdd.f"
	    lwork_zunmbr_qln_mm__ = (integer) cdum[0].r;

#line 384 "zgesdd.f"
	    zunmbr_("Q", "L", "N", m, n, n, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 386 "zgesdd.f"
	    lwork_zunmbr_qln_mn__ = (integer) cdum[0].r;

#line 388 "zgesdd.f"
	    zunmbr_("Q", "L", "N", n, n, n, cdum, n, cdum, cdum, n, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 390 "zgesdd.f"
	    lwork_zunmbr_qln_nn__ = (integer) cdum[0].r;

#line 392 "zgesdd.f"
	    if (*m >= mnthr1) {
#line 393 "zgesdd.f"
		if (wntqn) {

/*                 Path 1 (M >> N, JOBZ='N') */

#line 397 "zgesdd.f"
		    maxwrk = *n + lwork_zgeqrf_mn__;
/* Computing MAX */
#line 398 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zgebrd_nn__;
#line 398 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 399 "zgesdd.f"
		    minwrk = *n * 3;
#line 400 "zgesdd.f"
		} else if (wntqo) {

/*                 Path 2 (M >> N, JOBZ='O') */

#line 404 "zgesdd.f"
		    wrkbl = *n + lwork_zgeqrf_mn__;
/* Computing MAX */
#line 405 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_zungqr_mn__;
#line 405 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 406 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zgebrd_nn__;
#line 406 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 407 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zunmbr_qln_nn__;
#line 407 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 408 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zunmbr_prc_nn__;
#line 408 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 409 "zgesdd.f"
		    maxwrk = *m * *n + *n * *n + wrkbl;
#line 410 "zgesdd.f"
		    minwrk = (*n << 1) * *n + *n * 3;
#line 411 "zgesdd.f"
		} else if (wntqs) {

/*                 Path 3 (M >> N, JOBZ='S') */

#line 415 "zgesdd.f"
		    wrkbl = *n + lwork_zgeqrf_mn__;
/* Computing MAX */
#line 416 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_zungqr_mn__;
#line 416 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 417 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zgebrd_nn__;
#line 417 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 418 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zunmbr_qln_nn__;
#line 418 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 419 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zunmbr_prc_nn__;
#line 419 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 420 "zgesdd.f"
		    maxwrk = *n * *n + wrkbl;
#line 421 "zgesdd.f"
		    minwrk = *n * *n + *n * 3;
#line 422 "zgesdd.f"
		} else if (wntqa) {

/*                 Path 4 (M >> N, JOBZ='A') */

#line 426 "zgesdd.f"
		    wrkbl = *n + lwork_zgeqrf_mn__;
/* Computing MAX */
#line 427 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *n + lwork_zungqr_mm__;
#line 427 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 428 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zgebrd_nn__;
#line 428 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 429 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zunmbr_qln_nn__;
#line 429 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 430 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*n << 1) + lwork_zunmbr_prc_nn__;
#line 430 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 431 "zgesdd.f"
		    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
#line 432 "zgesdd.f"
		    i__1 = *n * 3, i__2 = *n + *m;
#line 432 "zgesdd.f"
		    minwrk = *n * *n + max(i__1,i__2);
#line 433 "zgesdd.f"
		}
#line 434 "zgesdd.f"
	    } else if (*m >= mnthr2) {

/*              Path 5 (M >> N, but not as much as MNTHR1) */

#line 438 "zgesdd.f"
		maxwrk = (*n << 1) + lwork_zgebrd_mn__;
#line 439 "zgesdd.f"
		minwrk = (*n << 1) + *m;
#line 440 "zgesdd.f"
		if (wntqo) {
/*                 Path 5o (M >> N, JOBZ='O') */
/* Computing MAX */
#line 442 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zungbr_p_nn__;
#line 442 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 443 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zungbr_q_mn__;
#line 443 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 444 "zgesdd.f"
		    maxwrk += *m * *n;
#line 445 "zgesdd.f"
		    minwrk += *n * *n;
#line 446 "zgesdd.f"
		} else if (wntqs) {
/*                 Path 5s (M >> N, JOBZ='S') */
/* Computing MAX */
#line 448 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zungbr_p_nn__;
#line 448 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 449 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zungbr_q_mn__;
#line 449 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 450 "zgesdd.f"
		} else if (wntqa) {
/*                 Path 5a (M >> N, JOBZ='A') */
/* Computing MAX */
#line 452 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zungbr_p_nn__;
#line 452 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 453 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zungbr_q_mm__;
#line 453 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 454 "zgesdd.f"
		}
#line 455 "zgesdd.f"
	    } else {

/*              Path 6 (M >= N, but not much larger) */

#line 459 "zgesdd.f"
		maxwrk = (*n << 1) + lwork_zgebrd_mn__;
#line 460 "zgesdd.f"
		minwrk = (*n << 1) + *m;
#line 461 "zgesdd.f"
		if (wntqo) {
/*                 Path 6o (M >= N, JOBZ='O') */
/* Computing MAX */
#line 463 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zunmbr_prc_nn__;
#line 463 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 464 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zunmbr_qln_mn__;
#line 464 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 465 "zgesdd.f"
		    maxwrk += *m * *n;
#line 466 "zgesdd.f"
		    minwrk += *n * *n;
#line 467 "zgesdd.f"
		} else if (wntqs) {
/*                 Path 6s (M >= N, JOBZ='S') */
/* Computing MAX */
#line 469 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zunmbr_qln_mn__;
#line 469 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 470 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zunmbr_prc_nn__;
#line 470 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 471 "zgesdd.f"
		} else if (wntqa) {
/*                 Path 6a (M >= N, JOBZ='A') */
/* Computing MAX */
#line 473 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zunmbr_qln_mm__;
#line 473 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 474 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*n << 1) + lwork_zunmbr_prc_nn__;
#line 474 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 475 "zgesdd.f"
		}
#line 476 "zgesdd.f"
	    }
#line 477 "zgesdd.f"
	} else if (minmn > 0) {

/*           There is no complex work space needed for bidiagonal SVD */
/*           The real work space needed for bidiagonal SVD (dbdsdc) is */
/*           BDSPAC = 3*M*M + 4*M for singular values and vectors; */
/*           BDSPAC = 4*M         for singular values only; */
/*           not including e, RU, and RVT matrices. */

/*           Compute space preferred for each routine */
#line 486 "zgesdd.f"
	    zgebrd_(m, n, cdum, m, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 488 "zgesdd.f"
	    lwork_zgebrd_mn__ = (integer) cdum[0].r;

#line 490 "zgesdd.f"
	    zgebrd_(m, m, cdum, m, dum, dum, cdum, cdum, cdum, &c_n1, &ierr);
#line 492 "zgesdd.f"
	    lwork_zgebrd_mm__ = (integer) cdum[0].r;

#line 494 "zgesdd.f"
	    zgelqf_(m, n, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 495 "zgesdd.f"
	    lwork_zgelqf_mn__ = (integer) cdum[0].r;

#line 497 "zgesdd.f"
	    zungbr_("P", m, n, m, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 499 "zgesdd.f"
	    lwork_zungbr_p_mn__ = (integer) cdum[0].r;

#line 501 "zgesdd.f"
	    zungbr_("P", n, n, m, cdum, n, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 503 "zgesdd.f"
	    lwork_zungbr_p_nn__ = (integer) cdum[0].r;

#line 505 "zgesdd.f"
	    zungbr_("Q", m, m, n, cdum, m, cdum, cdum, &c_n1, &ierr, (ftnlen)
		    1);
#line 507 "zgesdd.f"
	    lwork_zungbr_q_mm__ = (integer) cdum[0].r;

#line 509 "zgesdd.f"
	    zunglq_(m, n, m, cdum, m, cdum, cdum, &c_n1, &ierr);
#line 511 "zgesdd.f"
	    lwork_zunglq_mn__ = (integer) cdum[0].r;

#line 513 "zgesdd.f"
	    zunglq_(n, n, m, cdum, n, cdum, cdum, &c_n1, &ierr);
#line 515 "zgesdd.f"
	    lwork_zunglq_nn__ = (integer) cdum[0].r;

#line 517 "zgesdd.f"
	    zunmbr_("P", "R", "C", m, m, m, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 519 "zgesdd.f"
	    lwork_zunmbr_prc_mm__ = (integer) cdum[0].r;

#line 521 "zgesdd.f"
	    zunmbr_("P", "R", "C", m, n, m, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 523 "zgesdd.f"
	    lwork_zunmbr_prc_mn__ = (integer) cdum[0].r;

#line 525 "zgesdd.f"
	    zunmbr_("P", "R", "C", n, n, m, cdum, n, cdum, cdum, n, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 527 "zgesdd.f"
	    lwork_zunmbr_prc_nn__ = (integer) cdum[0].r;

#line 529 "zgesdd.f"
	    zunmbr_("Q", "L", "N", m, m, m, cdum, m, cdum, cdum, m, cdum, &
		    c_n1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 531 "zgesdd.f"
	    lwork_zunmbr_qln_mm__ = (integer) cdum[0].r;

#line 533 "zgesdd.f"
	    if (*n >= mnthr1) {
#line 534 "zgesdd.f"
		if (wntqn) {

/*                 Path 1t (N >> M, JOBZ='N') */

#line 538 "zgesdd.f"
		    maxwrk = *m + lwork_zgelqf_mn__;
/* Computing MAX */
#line 539 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zgebrd_mm__;
#line 539 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 540 "zgesdd.f"
		    minwrk = *m * 3;
#line 541 "zgesdd.f"
		} else if (wntqo) {

/*                 Path 2t (N >> M, JOBZ='O') */

#line 545 "zgesdd.f"
		    wrkbl = *m + lwork_zgelqf_mn__;
/* Computing MAX */
#line 546 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_zunglq_mn__;
#line 546 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 547 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zgebrd_mm__;
#line 547 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 548 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zunmbr_qln_mm__;
#line 548 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 549 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zunmbr_prc_mm__;
#line 549 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 550 "zgesdd.f"
		    maxwrk = *m * *n + *m * *m + wrkbl;
#line 551 "zgesdd.f"
		    minwrk = (*m << 1) * *m + *m * 3;
#line 552 "zgesdd.f"
		} else if (wntqs) {

/*                 Path 3t (N >> M, JOBZ='S') */

#line 556 "zgesdd.f"
		    wrkbl = *m + lwork_zgelqf_mn__;
/* Computing MAX */
#line 557 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_zunglq_mn__;
#line 557 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 558 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zgebrd_mm__;
#line 558 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 559 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zunmbr_qln_mm__;
#line 559 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 560 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zunmbr_prc_mm__;
#line 560 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 561 "zgesdd.f"
		    maxwrk = *m * *m + wrkbl;
#line 562 "zgesdd.f"
		    minwrk = *m * *m + *m * 3;
#line 563 "zgesdd.f"
		} else if (wntqa) {

/*                 Path 4t (N >> M, JOBZ='A') */

#line 567 "zgesdd.f"
		    wrkbl = *m + lwork_zgelqf_mn__;
/* Computing MAX */
#line 568 "zgesdd.f"
		    i__1 = wrkbl, i__2 = *m + lwork_zunglq_nn__;
#line 568 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 569 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zgebrd_mm__;
#line 569 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 570 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zunmbr_qln_mm__;
#line 570 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
/* Computing MAX */
#line 571 "zgesdd.f"
		    i__1 = wrkbl, i__2 = (*m << 1) + lwork_zunmbr_prc_mm__;
#line 571 "zgesdd.f"
		    wrkbl = max(i__1,i__2);
#line 572 "zgesdd.f"
		    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
#line 573 "zgesdd.f"
		    i__1 = *m * 3, i__2 = *m + *n;
#line 573 "zgesdd.f"
		    minwrk = *m * *m + max(i__1,i__2);
#line 574 "zgesdd.f"
		}
#line 575 "zgesdd.f"
	    } else if (*n >= mnthr2) {

/*              Path 5t (N >> M, but not as much as MNTHR1) */

#line 579 "zgesdd.f"
		maxwrk = (*m << 1) + lwork_zgebrd_mn__;
#line 580 "zgesdd.f"
		minwrk = (*m << 1) + *n;
#line 581 "zgesdd.f"
		if (wntqo) {
/*                 Path 5to (N >> M, JOBZ='O') */
/* Computing MAX */
#line 583 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zungbr_q_mm__;
#line 583 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 584 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zungbr_p_mn__;
#line 584 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 585 "zgesdd.f"
		    maxwrk += *m * *n;
#line 586 "zgesdd.f"
		    minwrk += *m * *m;
#line 587 "zgesdd.f"
		} else if (wntqs) {
/*                 Path 5ts (N >> M, JOBZ='S') */
/* Computing MAX */
#line 589 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zungbr_q_mm__;
#line 589 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 590 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zungbr_p_mn__;
#line 590 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 591 "zgesdd.f"
		} else if (wntqa) {
/*                 Path 5ta (N >> M, JOBZ='A') */
/* Computing MAX */
#line 593 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zungbr_q_mm__;
#line 593 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 594 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zungbr_p_nn__;
#line 594 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 595 "zgesdd.f"
		}
#line 596 "zgesdd.f"
	    } else {

/*              Path 6t (N > M, but not much larger) */

#line 600 "zgesdd.f"
		maxwrk = (*m << 1) + lwork_zgebrd_mn__;
#line 601 "zgesdd.f"
		minwrk = (*m << 1) + *n;
#line 602 "zgesdd.f"
		if (wntqo) {
/*                 Path 6to (N > M, JOBZ='O') */
/* Computing MAX */
#line 604 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zunmbr_qln_mm__;
#line 604 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 605 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zunmbr_prc_mn__;
#line 605 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 606 "zgesdd.f"
		    maxwrk += *m * *n;
#line 607 "zgesdd.f"
		    minwrk += *m * *m;
#line 608 "zgesdd.f"
		} else if (wntqs) {
/*                 Path 6ts (N > M, JOBZ='S') */
/* Computing MAX */
#line 610 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zunmbr_qln_mm__;
#line 610 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 611 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zunmbr_prc_mn__;
#line 611 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 612 "zgesdd.f"
		} else if (wntqa) {
/*                 Path 6ta (N > M, JOBZ='A') */
/* Computing MAX */
#line 614 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zunmbr_qln_mm__;
#line 614 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 615 "zgesdd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zunmbr_prc_nn__;
#line 615 "zgesdd.f"
		    maxwrk = max(i__1,i__2);
#line 616 "zgesdd.f"
		}
#line 617 "zgesdd.f"
	    }
#line 618 "zgesdd.f"
	}
#line 619 "zgesdd.f"
	maxwrk = max(maxwrk,minwrk);
#line 620 "zgesdd.f"
    }
#line 621 "zgesdd.f"
    if (*info == 0) {
#line 622 "zgesdd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 623 "zgesdd.f"
	if (*lwork < minwrk && ! lquery) {
#line 624 "zgesdd.f"
	    *info = -12;
#line 625 "zgesdd.f"
	}
#line 626 "zgesdd.f"
    }

#line 628 "zgesdd.f"
    if (*info != 0) {
#line 629 "zgesdd.f"
	i__1 = -(*info);
#line 629 "zgesdd.f"
	xerbla_("ZGESDD", &i__1, (ftnlen)6);
#line 630 "zgesdd.f"
	return 0;
#line 631 "zgesdd.f"
    } else if (lquery) {
#line 632 "zgesdd.f"
	return 0;
#line 633 "zgesdd.f"
    }

/*     Quick return if possible */

#line 637 "zgesdd.f"
    if (*m == 0 || *n == 0) {
#line 638 "zgesdd.f"
	return 0;
#line 639 "zgesdd.f"
    }

/*     Get machine constants */

#line 643 "zgesdd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 644 "zgesdd.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 645 "zgesdd.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 649 "zgesdd.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 650 "zgesdd.f"
    iscl = 0;
#line 651 "zgesdd.f"
    if (anrm > 0. && anrm < smlnum) {
#line 652 "zgesdd.f"
	iscl = 1;
#line 653 "zgesdd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 654 "zgesdd.f"
    } else if (anrm > bignum) {
#line 655 "zgesdd.f"
	iscl = 1;
#line 656 "zgesdd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 657 "zgesdd.f"
    }

#line 659 "zgesdd.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

#line 665 "zgesdd.f"
	if (*m >= mnthr1) {

#line 667 "zgesdd.f"
	    if (wntqn) {

/*              Path 1 (M >> N, JOBZ='N') */
/*              No singular vectors to be computed */

#line 672 "zgesdd.f"
		itau = 1;
#line 673 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              CWorkspace: need   N [tau] + N    [work] */
/*              CWorkspace: prefer N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 680 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 680 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Zero out below R */

#line 685 "zgesdd.f"
		i__1 = *n - 1;
#line 685 "zgesdd.f"
		i__2 = *n - 1;
#line 685 "zgesdd.f"
		zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 687 "zgesdd.f"
		ie = 1;
#line 688 "zgesdd.f"
		itauq = 1;
#line 689 "zgesdd.f"
		itaup = itauq + *n;
#line 690 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              CWorkspace: need   2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 697 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 697 "zgesdd.f"
		zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
#line 700 "zgesdd.f"
		nrwork = ie + *n;

/*              Perform bidiagonal SVD, compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + BDSPAC */

#line 706 "zgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 709 "zgesdd.f"
	    } else if (wntqo) {

/*              Path 2 (M >> N, JOBZ='O') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

#line 715 "zgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 719 "zgesdd.f"
		ldwrku = *n;
#line 720 "zgesdd.f"
		ir = iu + ldwrku * *n;
#line 721 "zgesdd.f"
		if (*lwork >= *m * *n + *n * *n + *n * 3) {

/*                 WORK(IR) is M by N */

#line 725 "zgesdd.f"
		    ldwrkr = *m;
#line 726 "zgesdd.f"
		} else {
#line 727 "zgesdd.f"
		    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
#line 728 "zgesdd.f"
		}
#line 729 "zgesdd.f"
		itau = ir + ldwrkr * *n;
#line 730 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 737 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 737 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy R to WORK( IR ), zeroing out below it */

#line 742 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 743 "zgesdd.f"
		i__1 = *n - 1;
#line 743 "zgesdd.f"
		i__2 = *n - 1;
#line 743 "zgesdd.f"
		zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 751 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 751 "zgesdd.f"
		zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 753 "zgesdd.f"
		ie = 1;
#line 754 "zgesdd.f"
		itauq = itau;
#line 755 "zgesdd.f"
		itaup = itauq + *n;
#line 756 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 763 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 763 "zgesdd.f"
		zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of R in WORK(IRU) and computing right singular vectors */
/*              of R in WORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 773 "zgesdd.f"
		iru = ie + *n;
#line 774 "zgesdd.f"
		irvt = iru + *n * *n;
#line 775 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 776 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of R */
/*              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 786 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 788 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 788 "zgesdd.f"
		zunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by the right singular vectors of R */
/*              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 798 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 799 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 799 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IU), storing result in WORK(IR) and copying to A */
/*              CWorkspace: need   N*N [U] + N*N [R] */
/*              CWorkspace: prefer N*N [U] + M*N [R] */
/*              RWorkspace: need   0 */

#line 809 "zgesdd.f"
		i__1 = *m;
#line 809 "zgesdd.f"
		i__2 = ldwrkr;
#line 809 "zgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 810 "zgesdd.f"
		    i__3 = *m - i__ + 1;
#line 810 "zgesdd.f"
		    chunk = min(i__3,ldwrkr);
#line 811 "zgesdd.f"
		    zgemm_("N", "N", &chunk, n, n, &c_b2, &a[i__ + a_dim1], 
			    lda, &work[iu], &ldwrku, &c_b1, &work[ir], &
			    ldwrkr, (ftnlen)1, (ftnlen)1);
#line 814 "zgesdd.f"
		    zlacpy_("F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 816 "zgesdd.f"
/* L10: */
#line 816 "zgesdd.f"
		}

#line 818 "zgesdd.f"
	    } else if (wntqs) {

/*              Path 3 (M >> N, JOBZ='S') */
/*              N left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 824 "zgesdd.f"
		ir = 1;

/*              WORK(IR) is N by N */

#line 828 "zgesdd.f"
		ldwrkr = *n;
#line 829 "zgesdd.f"
		itau = ir + ldwrkr * *n;
#line 830 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R */
/*              CWorkspace: need   N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 837 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 837 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy R to WORK(IR), zeroing out below it */

#line 842 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 843 "zgesdd.f"
		i__2 = *n - 1;
#line 843 "zgesdd.f"
		i__1 = *n - 1;
#line 843 "zgesdd.f"
		zlaset_("L", &i__2, &i__1, &c_b1, &c_b1, &work[ir + 1], &
			ldwrkr, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   N*N [R] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 851 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 851 "zgesdd.f"
		zungqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 853 "zgesdd.f"
		ie = 1;
#line 854 "zgesdd.f"
		itauq = itau;
#line 855 "zgesdd.f"
		itaup = itauq + *n;
#line 856 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in WORK(IR) */
/*              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 863 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 863 "zgesdd.f"
		zgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 873 "zgesdd.f"
		iru = ie + *n;
#line 874 "zgesdd.f"
		irvt = iru + *n * *n;
#line 875 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 876 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of R */
/*              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 886 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 887 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 887 "zgesdd.f"
		zunmbr_("Q", "L", "N", n, n, n, &work[ir], &ldwrkr, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 897 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 898 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 898 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &work[ir], &ldwrkr, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in A by left singular vectors of R in */
/*              WORK(IR), storing result in U */
/*              CWorkspace: need   N*N [R] */
/*              RWorkspace: need   0 */

#line 907 "zgesdd.f"
		zlacpy_("F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (
			ftnlen)1);
#line 908 "zgesdd.f"
		zgemm_("N", "N", m, n, n, &c_b2, &a[a_offset], lda, &work[ir],
			 &ldwrkr, &c_b1, &u[u_offset], ldu, (ftnlen)1, (
			ftnlen)1);

#line 911 "zgesdd.f"
	    } else if (wntqa) {

/*              Path 4 (M >> N, JOBZ='A') */
/*              M left singular vectors to be computed in U and */
/*              N right singular vectors to be computed in VT */

#line 917 "zgesdd.f"
		iu = 1;

/*              WORK(IU) is N by N */

#line 921 "zgesdd.f"
		ldwrku = *n;
#line 922 "zgesdd.f"
		itau = iu + ldwrku * *n;
#line 923 "zgesdd.f"
		nwork = itau + *n;

/*              Compute A=Q*R, copying result to U */
/*              CWorkspace: need   N*N [U] + N [tau] + N    [work] */
/*              CWorkspace: prefer N*N [U] + N [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 930 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 930 "zgesdd.f"
		zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);
#line 932 "zgesdd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Generate Q in U */
/*              CWorkspace: need   N*N [U] + N [tau] + M    [work] */
/*              CWorkspace: prefer N*N [U] + N [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 939 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 939 "zgesdd.f"
		zungqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__2, &ierr);

/*              Produce R in A, zeroing out below it */

#line 944 "zgesdd.f"
		i__2 = *n - 1;
#line 944 "zgesdd.f"
		i__1 = *n - 1;
#line 944 "zgesdd.f"
		zlaset_("L", &i__2, &i__1, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 946 "zgesdd.f"
		ie = 1;
#line 947 "zgesdd.f"
		itauq = itau;
#line 948 "zgesdd.f"
		itaup = itauq + *n;
#line 949 "zgesdd.f"
		nwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N      [work] */
/*              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + 2*N*NB [work] */
/*              RWorkspace: need   N [e] */

#line 956 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 956 "zgesdd.f"
		zgebrd_(n, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 959 "zgesdd.f"
		iru = ie + *n;
#line 960 "zgesdd.f"
		irvt = iru + *n * *n;
#line 961 "zgesdd.f"
		nrwork = irvt + *n * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 969 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by left singular vectors of R */
/*              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 979 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			ftnlen)1);
#line 981 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 981 "zgesdd.f"
		zunmbr_("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of R */
/*              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 991 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 992 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 992 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by left singular vectors of R in */
/*              WORK(IU), storing result in A */
/*              CWorkspace: need   N*N [U] */
/*              RWorkspace: need   0 */

#line 1001 "zgesdd.f"
		zgemm_("N", "N", m, n, n, &c_b2, &u[u_offset], ldu, &work[iu],
			 &ldwrku, &c_b1, &a[a_offset], lda, (ftnlen)1, (
			ftnlen)1);

/*              Copy left singular vectors of A from A to U */

#line 1006 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

#line 1008 "zgesdd.f"
	    }

#line 1010 "zgesdd.f"
	} else if (*m >= mnthr2) {

/*           MNTHR2 <= M < MNTHR1 */

/*           Path 5 (M >> N, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           ZUNGBR and matrix multiplication to compute singular vectors */

#line 1018 "zgesdd.f"
	    ie = 1;
#line 1019 "zgesdd.f"
	    nrwork = ie + *n;
#line 1020 "zgesdd.f"
	    itauq = 1;
#line 1021 "zgesdd.f"
	    itaup = itauq + *n;
#line 1022 "zgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*N [tauq, taup] + M        [work] */
/*           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   N [e] */

#line 1029 "zgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 1029 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 1032 "zgesdd.f"
	    if (wntqn) {

/*              Path 5n (M >> N, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + BDSPAC */

#line 1039 "zgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1041 "zgesdd.f"
	    } else if (wntqo) {
#line 1042 "zgesdd.f"
		iu = nwork;
#line 1043 "zgesdd.f"
		iru = nrwork;
#line 1044 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1045 "zgesdd.f"
		nrwork = irvt + *n * *n;

/*              Path 5o (M >> N, JOBZ='O') */
/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1053 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1054 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1054 "zgesdd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1062 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1062 "zgesdd.f"
		zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

#line 1065 "zgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 1069 "zgesdd.f"
		    ldwrku = *m;
#line 1070 "zgesdd.f"
		} else {

/*                 WORK(IU) is LDWRKU by N */

#line 1074 "zgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 1075 "zgesdd.f"
		}
#line 1076 "zgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1084 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in WORK(IU), copying to VT */
/*              CWorkspace: need   2*N [tauq, taup] + N*N [U] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork] */

#line 1093 "zgesdd.f"
		zlarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &work[iu]
			, &ldwrku, &rwork[nrwork]);
#line 1095 "zgesdd.f"
		zlacpy_("F", n, n, &work[iu], &ldwrku, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in A by real matrix RWORK(IRU), storing the */
/*              result in WORK(IU), copying to A */
/*              CWorkspace: need   2*N [tauq, taup] + N*N [U] */
/*              CWorkspace: prefer 2*N [tauq, taup] + M*N [U] */
/*              RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork] */
/*              RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1104 "zgesdd.f"
		nrwork = irvt;
#line 1105 "zgesdd.f"
		i__2 = *m;
#line 1105 "zgesdd.f"
		i__1 = ldwrku;
#line 1105 "zgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1106 "zgesdd.f"
		    i__3 = *m - i__ + 1;
#line 1106 "zgesdd.f"
		    chunk = min(i__3,ldwrku);
#line 1107 "zgesdd.f"
		    zlacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru], n, 
			    &work[iu], &ldwrku, &rwork[nrwork]);
#line 1109 "zgesdd.f"
		    zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
			    a_dim1], lda, (ftnlen)1);
#line 1111 "zgesdd.f"
/* L20: */
#line 1111 "zgesdd.f"
		}

#line 1113 "zgesdd.f"
	    } else if (wntqs) {

/*              Path 5s (M >> N, JOBZ='S') */
/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1121 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1122 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1122 "zgesdd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1130 "zgesdd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1131 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1131 "zgesdd.f"
		zungbr_("Q", m, n, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1140 "zgesdd.f"
		iru = nrwork;
#line 1141 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1142 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1143 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork] */

#line 1152 "zgesdd.f"
		zlarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1154 "zgesdd.f"
		zlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1161 "zgesdd.f"
		nrwork = irvt;
#line 1162 "zgesdd.f"
		zlacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1164 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1165 "zgesdd.f"
	    } else {

/*              Path 5a (M >> N, JOBZ='A') */
/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1173 "zgesdd.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1174 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1174 "zgesdd.f"
		zungbr_("P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__1, &ierr, (ftnlen)1);

/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*N [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1182 "zgesdd.f"
		zlacpy_("L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1183 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1183 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1192 "zgesdd.f"
		iru = nrwork;
#line 1193 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1194 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1195 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork] */

#line 1204 "zgesdd.f"
		zlarcm_(n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1206 "zgesdd.f"
		zlacpy_("F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1213 "zgesdd.f"
		nrwork = irvt;
#line 1214 "zgesdd.f"
		zlacrm_(m, n, &u[u_offset], ldu, &rwork[iru], n, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1216 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1217 "zgesdd.f"
	    }

#line 1219 "zgesdd.f"
	} else {

/*           M .LT. MNTHR2 */

/*           Path 6 (M >= N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */
/*           Use ZUNMBR to compute singular vectors */

#line 1227 "zgesdd.f"
	    ie = 1;
#line 1228 "zgesdd.f"
	    nrwork = ie + *n;
#line 1229 "zgesdd.f"
	    itauq = 1;
#line 1230 "zgesdd.f"
	    itaup = itauq + *n;
#line 1231 "zgesdd.f"
	    nwork = itaup + *n;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*N [tauq, taup] + M        [work] */
/*           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   N [e] */

#line 1238 "zgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1238 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);
#line 1241 "zgesdd.f"
	    if (wntqn) {

/*              Path 6n (M >= N, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + BDSPAC */

#line 1248 "zgesdd.f"
		dbdsdc_("U", "N", n, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1250 "zgesdd.f"
	    } else if (wntqo) {
#line 1251 "zgesdd.f"
		iu = nwork;
#line 1252 "zgesdd.f"
		iru = nrwork;
#line 1253 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1254 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1255 "zgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 WORK( IU ) is M by N */

#line 1259 "zgesdd.f"
		    ldwrku = *m;
#line 1260 "zgesdd.f"
		} else {

/*                 WORK( IU ) is LDWRKU by N */

#line 1264 "zgesdd.f"
		    ldwrku = (*lwork - *n * 3) / *n;
#line 1265 "zgesdd.f"
		}
#line 1266 "zgesdd.f"
		nwork = iu + ldwrku * *n;

/*              Path 6o (M >= N, JOBZ='O') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1275 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1285 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1286 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1286 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1290 "zgesdd.f"
		if (*lwork >= *m * *n + *n * 3) {

/*                 Path 6o-fast */
/*                 Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*                 Overwrite WORK(IU) by left singular vectors of A, copying */
/*                 to A */
/*                 CWorkspace: need   2*N [tauq, taup] + M*N [U] + N    [work] */
/*                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U] + N*NB [work] */
/*                 RWorkspace: need   N [e] + N*N [RU] */

#line 1300 "zgesdd.f"
		    zlaset_("F", m, n, &c_b1, &c_b1, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1302 "zgesdd.f"
		    zlacp2_("F", n, n, &rwork[iru], n, &work[iu], &ldwrku, (
			    ftnlen)1);
#line 1304 "zgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1304 "zgesdd.f"
		    zunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			    itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1307 "zgesdd.f"
		    zlacpy_("F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, 
			    (ftnlen)1);
#line 1308 "zgesdd.f"
		} else {

/*                 Path 6o-slow */
/*                 Generate Q in A */
/*                 CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work] */
/*                 CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work] */
/*                 RWorkspace: need   0 */

#line 1316 "zgesdd.f"
		    i__1 = *lwork - nwork + 1;
#line 1316 "zgesdd.f"
		    zungbr_("Q", m, n, n, &a[a_offset], lda, &work[itauq], &
			    work[nwork], &i__1, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 CWorkspace: need   2*N [tauq, taup] + N*N [U] */
/*                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U] */
/*                 RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork] */
/*                 RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here */

#line 1326 "zgesdd.f"
		    nrwork = irvt;
#line 1327 "zgesdd.f"
		    i__1 = *m;
#line 1327 "zgesdd.f"
		    i__2 = ldwrku;
#line 1327 "zgesdd.f"
		    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ +=
			     i__2) {
/* Computing MIN */
#line 1328 "zgesdd.f"
			i__3 = *m - i__ + 1;
#line 1328 "zgesdd.f"
			chunk = min(i__3,ldwrku);
#line 1329 "zgesdd.f"
			zlacrm_(&chunk, n, &a[i__ + a_dim1], lda, &rwork[iru],
				 n, &work[iu], &ldwrku, &rwork[nrwork]);
#line 1332 "zgesdd.f"
			zlacpy_("F", &chunk, n, &work[iu], &ldwrku, &a[i__ + 
				a_dim1], lda, (ftnlen)1);
#line 1334 "zgesdd.f"
/* L30: */
#line 1334 "zgesdd.f"
		    }
#line 1335 "zgesdd.f"
		}

#line 1337 "zgesdd.f"
	    } else if (wntqs) {

/*              Path 6s (M >= N, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1346 "zgesdd.f"
		iru = nrwork;
#line 1347 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1348 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1349 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1359 "zgesdd.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1360 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1361 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1361 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, n, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1371 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1372 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1372 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1375 "zgesdd.f"
	    } else {

/*              Path 6a (M >= N, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC */

#line 1384 "zgesdd.f"
		iru = nrwork;
#line 1385 "zgesdd.f"
		irvt = iru + *n * *n;
#line 1386 "zgesdd.f"
		nrwork = irvt + *n * *n;
#line 1387 "zgesdd.f"
		dbdsdc_("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &
			rwork[irvt], n, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Set the right corner of U to identity matrix */

#line 1393 "zgesdd.f"
		zlaset_("F", m, m, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)1)
			;
#line 1394 "zgesdd.f"
		if (*m > *n) {
#line 1395 "zgesdd.f"
		    i__2 = *m - *n;
#line 1395 "zgesdd.f"
		    i__1 = *m - *n;
#line 1395 "zgesdd.f"
		    zlaset_("F", &i__2, &i__1, &c_b1, &c_b2, &u[*n + 1 + (*n 
			    + 1) * u_dim1], ldu, (ftnlen)1);
#line 1397 "zgesdd.f"
		}

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1405 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[iru], n, &u[u_offset], ldu, (ftnlen)
			1);
#line 1406 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1406 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*N [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] */

#line 1416 "zgesdd.f"
		zlacp2_("F", n, n, &rwork[irvt], n, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1417 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1417 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1420 "zgesdd.f"
	    }

#line 1422 "zgesdd.f"
	}

#line 1424 "zgesdd.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

#line 1430 "zgesdd.f"
	if (*n >= mnthr1) {

#line 1432 "zgesdd.f"
	    if (wntqn) {

/*              Path 1t (N >> M, JOBZ='N') */
/*              No singular vectors to be computed */

#line 1437 "zgesdd.f"
		itau = 1;
#line 1438 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              CWorkspace: need   M [tau] + M    [work] */
/*              CWorkspace: prefer M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1445 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1445 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Zero out above L */

#line 1450 "zgesdd.f"
		i__2 = *m - 1;
#line 1450 "zgesdd.f"
		i__1 = *m - 1;
#line 1450 "zgesdd.f"
		zlaset_("U", &i__2, &i__1, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1452 "zgesdd.f"
		ie = 1;
#line 1453 "zgesdd.f"
		itauq = 1;
#line 1454 "zgesdd.f"
		itaup = itauq + *m;
#line 1455 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              CWorkspace: need   2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1462 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1462 "zgesdd.f"
		zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);
#line 1465 "zgesdd.f"
		nrwork = ie + *m;

/*              Perform bidiagonal SVD, compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + BDSPAC */

#line 1471 "zgesdd.f"
		dbdsdc_("U", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);

#line 1474 "zgesdd.f"
	    } else if (wntqo) {

/*              Path 2t (N >> M, JOBZ='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

#line 1480 "zgesdd.f"
		ivt = 1;
#line 1481 "zgesdd.f"
		ldwkvt = *m;

/*              WORK(IVT) is M by M */

#line 1485 "zgesdd.f"
		il = ivt + ldwkvt * *m;
#line 1486 "zgesdd.f"
		if (*lwork >= *m * *n + *m * *m + *m * 3) {

/*                 WORK(IL) M by N */

#line 1490 "zgesdd.f"
		    ldwrkl = *m;
#line 1491 "zgesdd.f"
		    chunk = *n;
#line 1492 "zgesdd.f"
		} else {

/*                 WORK(IL) is M by CHUNK */

#line 1496 "zgesdd.f"
		    ldwrkl = *m;
#line 1497 "zgesdd.f"
		    chunk = (*lwork - *m * *m - *m * 3) / *m;
#line 1498 "zgesdd.f"
		}
#line 1499 "zgesdd.f"
		itau = il + ldwrkl * chunk;
#line 1500 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1507 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1507 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__2, &ierr);

/*              Copy L to WORK(IL), zeroing about above it */

#line 1512 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1513 "zgesdd.f"
		i__2 = *m - 1;
#line 1513 "zgesdd.f"
		i__1 = *m - 1;
#line 1513 "zgesdd.f"
		zlaset_("U", &i__2, &i__1, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1521 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1521 "zgesdd.f"
		zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__2, &ierr);
#line 1523 "zgesdd.f"
		ie = 1;
#line 1524 "zgesdd.f"
		itauq = itau;
#line 1525 "zgesdd.f"
		itaup = itauq + *m;
#line 1526 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1533 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1533 "zgesdd.f"
		zgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__2, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC */

#line 1543 "zgesdd.f"
		iru = ie + *m;
#line 1544 "zgesdd.f"
		irvt = iru + *m * *m;
#line 1545 "zgesdd.f"
		nrwork = irvt + *m * *m;
#line 1546 "zgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix WORK(IU) */
/*              Overwrite WORK(IU) by the left singular vectors of L */
/*              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1556 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1557 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1557 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by the right singular vectors of L */
/*              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1567 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1569 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1569 "zgesdd.f"
		zunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IL) by Q */
/*              in A, storing result in WORK(IL) and copying to A */
/*              CWorkspace: need   M*M [VT] + M*M [L] */
/*              CWorkspace: prefer M*M [VT] + M*N [L] */
/*              RWorkspace: need   0 */

#line 1579 "zgesdd.f"
		i__2 = *n;
#line 1579 "zgesdd.f"
		i__1 = chunk;
#line 1579 "zgesdd.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 1580 "zgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1580 "zgesdd.f"
		    blk = min(i__3,chunk);
#line 1581 "zgesdd.f"
		    zgemm_("N", "N", m, &blk, m, &c_b2, &work[ivt], m, &a[i__ 
			    * a_dim1 + 1], lda, &c_b1, &work[il], &ldwrkl, (
			    ftnlen)1, (ftnlen)1);
#line 1584 "zgesdd.f"
		    zlacpy_("F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 
			    + 1], lda, (ftnlen)1);
#line 1586 "zgesdd.f"
/* L40: */
#line 1586 "zgesdd.f"
		}

#line 1588 "zgesdd.f"
	    } else if (wntqs) {

/*              Path 3t (N >> M, JOBZ='S') */
/*              M right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1594 "zgesdd.f"
		il = 1;

/*              WORK(IL) is M by M */

#line 1598 "zgesdd.f"
		ldwrkl = *m;
#line 1599 "zgesdd.f"
		itau = il + ldwrkl * *m;
#line 1600 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q */
/*              CWorkspace: need   M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1607 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1607 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

/*              Copy L to WORK(IL), zeroing out above it */

#line 1612 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1613 "zgesdd.f"
		i__1 = *m - 1;
#line 1613 "zgesdd.f"
		i__2 = *m - 1;
#line 1613 "zgesdd.f"
		zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &work[il + ldwrkl], &
			ldwrkl, (ftnlen)1);

/*              Generate Q in A */
/*              CWorkspace: need   M*M [L] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1621 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1621 "zgesdd.f"
		zunglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork],
			 &i__1, &ierr);
#line 1623 "zgesdd.f"
		ie = 1;
#line 1624 "zgesdd.f"
		itauq = itau;
#line 1625 "zgesdd.f"
		itaup = itauq + *m;
#line 1626 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in WORK(IL) */
/*              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1633 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1633 "zgesdd.f"
		zgebrd_(m, m, &work[il], &ldwrkl, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC */

#line 1643 "zgesdd.f"
		iru = ie + *m;
#line 1644 "zgesdd.f"
		irvt = iru + *m * *m;
#line 1645 "zgesdd.f"
		nrwork = irvt + *m * *m;
#line 1646 "zgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1656 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1657 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1657 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, m, &work[il], &ldwrkl, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by left singular vectors of L */
/*              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1667 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1668 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1668 "zgesdd.f"
		zunmbr_("P", "R", "C", m, m, m, &work[il], &ldwrkl, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy VT to WORK(IL), multiply right singular vectors of L */
/*              in WORK(IL) by Q in A, storing result in VT */
/*              CWorkspace: need   M*M [L] */
/*              RWorkspace: need   0 */

#line 1677 "zgesdd.f"
		zlacpy_("F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (
			ftnlen)1);
#line 1678 "zgesdd.f"
		zgemm_("N", "N", m, n, m, &c_b2, &work[il], &ldwrkl, &a[
			a_offset], lda, &c_b1, &vt[vt_offset], ldvt, (ftnlen)
			1, (ftnlen)1);

#line 1681 "zgesdd.f"
	    } else if (wntqa) {

/*              Path 4t (N >> M, JOBZ='A') */
/*              N right singular vectors to be computed in VT and */
/*              M left singular vectors to be computed in U */

#line 1687 "zgesdd.f"
		ivt = 1;

/*              WORK(IVT) is M by M */

#line 1691 "zgesdd.f"
		ldwkvt = *m;
#line 1692 "zgesdd.f"
		itau = ivt + ldwkvt * *m;
#line 1693 "zgesdd.f"
		nwork = itau + *m;

/*              Compute A=L*Q, copying result to VT */
/*              CWorkspace: need   M*M [VT] + M [tau] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + M [tau] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1700 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1700 "zgesdd.f"
		zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
#line 1702 "zgesdd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Generate Q in VT */
/*              CWorkspace: need   M*M [VT] + M [tau] + N    [work] */
/*              CWorkspace: prefer M*M [VT] + M [tau] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1709 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1709 "zgesdd.f"
		zunglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

/*              Produce L in A, zeroing out above it */

#line 1714 "zgesdd.f"
		i__1 = *m - 1;
#line 1714 "zgesdd.f"
		i__2 = *m - 1;
#line 1714 "zgesdd.f"
		zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1716 "zgesdd.f"
		ie = 1;
#line 1717 "zgesdd.f"
		itauq = itau;
#line 1718 "zgesdd.f"
		itaup = itauq + *m;
#line 1719 "zgesdd.f"
		nwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M      [work] */
/*              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + 2*M*NB [work] */
/*              RWorkspace: need   M [e] */

#line 1726 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1726 "zgesdd.f"
		zgebrd_(m, m, &a[a_offset], lda, &s[1], &rwork[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC */

#line 1736 "zgesdd.f"
		iru = ie + *m;
#line 1737 "zgesdd.f"
		irvt = iru + *m * *m;
#line 1738 "zgesdd.f"
		nrwork = irvt + *m * *m;
#line 1739 "zgesdd.f"
		dbdsdc_("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of L */
/*              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1749 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 1750 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1750 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*              Overwrite WORK(IVT) by right singular vectors of L */
/*              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1760 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			ftnlen)1);
#line 1762 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1762 "zgesdd.f"
		zunmbr_("P", "R", "C", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Multiply right singular vectors of L in WORK(IVT) by */
/*              Q in VT, storing result in A */
/*              CWorkspace: need   M*M [VT] */
/*              RWorkspace: need   0 */

#line 1771 "zgesdd.f"
		zgemm_("N", "N", m, n, m, &c_b2, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &c_b1, &a[a_offset], lda, (ftnlen)1,
			 (ftnlen)1);

/*              Copy right singular vectors of A from A to VT */

#line 1776 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);

#line 1778 "zgesdd.f"
	    }

#line 1780 "zgesdd.f"
	} else if (*n >= mnthr2) {

/*           MNTHR2 <= N < MNTHR1 */

/*           Path 5t (N >> M, but not as much as MNTHR1) */
/*           Reduce to bidiagonal form without QR decomposition, use */
/*           ZUNGBR and matrix multiplication to compute singular vectors */

#line 1788 "zgesdd.f"
	    ie = 1;
#line 1789 "zgesdd.f"
	    nrwork = ie + *m;
#line 1790 "zgesdd.f"
	    itauq = 1;
#line 1791 "zgesdd.f"
	    itaup = itauq + *m;
#line 1792 "zgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*M [tauq, taup] + N        [work] */
/*           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   M [e] */

#line 1799 "zgesdd.f"
	    i__1 = *lwork - nwork + 1;
#line 1799 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, &ierr);

#line 1803 "zgesdd.f"
	    if (wntqn) {

/*              Path 5tn (N >> M, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + BDSPAC */

#line 1810 "zgesdd.f"
		dbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 1812 "zgesdd.f"
	    } else if (wntqo) {
#line 1813 "zgesdd.f"
		irvt = nrwork;
#line 1814 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1815 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1816 "zgesdd.f"
		ivt = nwork;

/*              Path 5to (N >> M, JOBZ='O') */
/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1824 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1825 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1825 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

/*              Generate P**H in A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1833 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 1833 "zgesdd.f"
		zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
			nwork], &i__1, &ierr, (ftnlen)1);

#line 1836 "zgesdd.f"
		ldwkvt = *m;
#line 1837 "zgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 1841 "zgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 1842 "zgesdd.f"
		    chunk = *n;
#line 1843 "zgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 1847 "zgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 1848 "zgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 1849 "zgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 1857 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRVT) */
/*              storing the result in WORK(IVT), copying to U */
/*              CWorkspace: need   2*M [tauq, taup] + M*M [VT] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork] */

#line 1866 "zgesdd.f"
		zlacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &work[ivt], &
			ldwkvt, &rwork[nrwork]);
#line 1868 "zgesdd.f"
		zlacpy_("F", m, m, &work[ivt], &ldwkvt, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply RWORK(IRVT) by P**H in A, storing the */
/*              result in WORK(IVT), copying to A */
/*              CWorkspace: need   2*M [tauq, taup] + M*M [VT] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] */
/*              RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork] */
/*              RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 1877 "zgesdd.f"
		nrwork = iru;
#line 1878 "zgesdd.f"
		i__1 = *n;
#line 1878 "zgesdd.f"
		i__2 = chunk;
#line 1878 "zgesdd.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 1879 "zgesdd.f"
		    i__3 = *n - i__ + 1;
#line 1879 "zgesdd.f"
		    blk = min(i__3,chunk);
#line 1880 "zgesdd.f"
		    zlarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1], 
			    lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 1882 "zgesdd.f"
		    zlacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
			    a_dim1 + 1], lda, (ftnlen)1);
#line 1884 "zgesdd.f"
/* L50: */
#line 1884 "zgesdd.f"
		}
#line 1885 "zgesdd.f"
	    } else if (wntqs) {

/*              Path 5ts (N >> M, JOBZ='S') */
/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1893 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1894 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1894 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1902 "zgesdd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1903 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1903 "zgesdd.f"
		zungbr_("P", m, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 1912 "zgesdd.f"
		irvt = nrwork;
#line 1913 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1914 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1915 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork] */

#line 1924 "zgesdd.f"
		zlacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1926 "zgesdd.f"
		zlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 1933 "zgesdd.f"
		nrwork = iru;
#line 1934 "zgesdd.f"
		zlarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1936 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1937 "zgesdd.f"
	    } else {

/*              Path 5ta (N >> M, JOBZ='A') */
/*              Copy A to U, generate Q */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   0 */

#line 1945 "zgesdd.f"
		zlacpy_("L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1946 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1946 "zgesdd.f"
		zungbr_("Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
			nwork], &i__2, &ierr, (ftnlen)1);

/*              Copy A to VT, generate P**H */
/*              CWorkspace: need   2*M [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   0 */

#line 1954 "zgesdd.f"
		zlacpy_("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1955 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 1955 "zgesdd.f"
		zungbr_("P", n, n, m, &vt[vt_offset], ldvt, &work[itaup], &
			work[nwork], &i__2, &ierr, (ftnlen)1);

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 1964 "zgesdd.f"
		irvt = nrwork;
#line 1965 "zgesdd.f"
		iru = irvt + *m * *m;
#line 1966 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 1967 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Multiply Q in U by real matrix RWORK(IRU), storing the */
/*              result in A, copying to U */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork] */

#line 1976 "zgesdd.f"
		zlacrm_(m, m, &u[u_offset], ldu, &rwork[iru], m, &a[a_offset],
			 lda, &rwork[nrwork]);
#line 1978 "zgesdd.f"
		zlacpy_("F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);

/*              Multiply real matrix RWORK(IRVT) by P**H in VT, */
/*              storing the result in A, copying to VT */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 1985 "zgesdd.f"
		nrwork = iru;
#line 1986 "zgesdd.f"
		zlarcm_(m, n, &rwork[irvt], m, &vt[vt_offset], ldvt, &a[
			a_offset], lda, &rwork[nrwork]);
#line 1988 "zgesdd.f"
		zlacpy_("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 1989 "zgesdd.f"
	    }

#line 1991 "zgesdd.f"
	} else {

/*           N .LT. MNTHR2 */

/*           Path 6t (N > M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           Use ZUNMBR to compute singular vectors */

#line 1999 "zgesdd.f"
	    ie = 1;
#line 2000 "zgesdd.f"
	    nrwork = ie + *m;
#line 2001 "zgesdd.f"
	    itauq = 1;
#line 2002 "zgesdd.f"
	    itaup = itauq + *m;
#line 2003 "zgesdd.f"
	    nwork = itaup + *m;

/*           Bidiagonalize A */
/*           CWorkspace: need   2*M [tauq, taup] + N        [work] */
/*           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work] */
/*           RWorkspace: need   M [e] */

#line 2010 "zgesdd.f"
	    i__2 = *lwork - nwork + 1;
#line 2010 "zgesdd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__2, &ierr);
#line 2013 "zgesdd.f"
	    if (wntqn) {

/*              Path 6tn (N > M, JOBZ='N') */
/*              Compute singular values only */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + BDSPAC */

#line 2020 "zgesdd.f"
		dbdsdc_("L", "N", m, &s[1], &rwork[ie], dum, &c__1, dum, &
			c__1, dum, idum, &rwork[nrwork], &iwork[1], info, (
			ftnlen)1, (ftnlen)1);
#line 2022 "zgesdd.f"
	    } else if (wntqo) {
/*              Path 6to (N > M, JOBZ='O') */
#line 2024 "zgesdd.f"
		ldwkvt = *m;
#line 2025 "zgesdd.f"
		ivt = nwork;
#line 2026 "zgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 WORK( IVT ) is M by N */

#line 2030 "zgesdd.f"
		    zlaset_("F", m, n, &c_b1, &c_b1, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 2032 "zgesdd.f"
		    nwork = ivt + ldwkvt * *n;
#line 2033 "zgesdd.f"
		} else {

/*                 WORK( IVT ) is M by CHUNK */

#line 2037 "zgesdd.f"
		    chunk = (*lwork - *m * 3) / *m;
#line 2038 "zgesdd.f"
		    nwork = ivt + ldwkvt * chunk;
#line 2039 "zgesdd.f"
		}

/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 2047 "zgesdd.f"
		irvt = nrwork;
#line 2048 "zgesdd.f"
		iru = irvt + *m * *m;
#line 2049 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 2050 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] */

#line 2060 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2061 "zgesdd.f"
		i__2 = *lwork - nwork + 1;
#line 2061 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 2065 "zgesdd.f"
		if (*lwork >= *m * *n + *m * 3) {

/*                 Path 6to-fast */
/*                 Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT) */
/*                 Overwrite WORK(IVT) by right singular vectors of A, */
/*                 copying to A */
/*                 CWorkspace: need   2*M [tauq, taup] + M*N [VT] + M    [work] */
/*                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] + M*NB [work] */
/*                 RWorkspace: need   M [e] + M*M [RVT] */

#line 2075 "zgesdd.f"
		    zlacp2_("F", m, m, &rwork[irvt], m, &work[ivt], &ldwkvt, (
			    ftnlen)1);
#line 2077 "zgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 2077 "zgesdd.f"
		    zunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			    itaup], &work[ivt], &ldwkvt, &work[nwork], &i__2, 
			    &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2080 "zgesdd.f"
		    zlacpy_("F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda,
			     (ftnlen)1);
#line 2081 "zgesdd.f"
		} else {

/*                 Path 6to-slow */
/*                 Generate P**H in A */
/*                 CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work] */
/*                 CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work] */
/*                 RWorkspace: need   0 */

#line 2089 "zgesdd.f"
		    i__2 = *lwork - nwork + 1;
#line 2089 "zgesdd.f"
		    zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &
			    work[nwork], &i__2, &ierr, (ftnlen)1);

/*                 Multiply Q in A by real matrix RWORK(IRU), storing the */
/*                 result in WORK(IU), copying to A */
/*                 CWorkspace: need   2*M [tauq, taup] + M*M [VT] */
/*                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] */
/*                 RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork] */
/*                 RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here */

#line 2099 "zgesdd.f"
		    nrwork = iru;
#line 2100 "zgesdd.f"
		    i__2 = *n;
#line 2100 "zgesdd.f"
		    i__1 = chunk;
#line 2100 "zgesdd.f"
		    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
			     i__1) {
/* Computing MIN */
#line 2101 "zgesdd.f"
			i__3 = *n - i__ + 1;
#line 2101 "zgesdd.f"
			blk = min(i__3,chunk);
#line 2102 "zgesdd.f"
			zlarcm_(m, &blk, &rwork[irvt], m, &a[i__ * a_dim1 + 1]
				, lda, &work[ivt], &ldwkvt, &rwork[nrwork]);
#line 2105 "zgesdd.f"
			zlacpy_("F", m, &blk, &work[ivt], &ldwkvt, &a[i__ * 
				a_dim1 + 1], lda, (ftnlen)1);
#line 2107 "zgesdd.f"
/* L60: */
#line 2107 "zgesdd.f"
		    }
#line 2108 "zgesdd.f"
		}
#line 2109 "zgesdd.f"
	    } else if (wntqs) {

/*              Path 6ts (N > M, JOBZ='S') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 2118 "zgesdd.f"
		irvt = nrwork;
#line 2119 "zgesdd.f"
		iru = irvt + *m * *m;
#line 2120 "zgesdd.f"
		nrwork = iru + *m * *m;
#line 2121 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] */

#line 2131 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2132 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2132 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] */

#line 2142 "zgesdd.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2143 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2144 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2144 "zgesdd.f"
		zunmbr_("P", "R", "C", m, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2147 "zgesdd.f"
	    } else {

/*              Path 6ta (N > M, JOBZ='A') */
/*              Perform bidiagonal SVD, computing left singular vectors */
/*              of bidiagonal matrix in RWORK(IRU) and computing right */
/*              singular vectors of bidiagonal matrix in RWORK(IRVT) */
/*              CWorkspace: need   0 */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC */

#line 2156 "zgesdd.f"
		irvt = nrwork;
#line 2157 "zgesdd.f"
		iru = irvt + *m * *m;
#line 2158 "zgesdd.f"
		nrwork = iru + *m * *m;

#line 2160 "zgesdd.f"
		dbdsdc_("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &
			rwork[irvt], m, dum, idum, &rwork[nrwork], &iwork[1], 
			info, (ftnlen)1, (ftnlen)1);

/*              Copy real matrix RWORK(IRU) to complex matrix U */
/*              Overwrite U by left singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + M    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] */

#line 2170 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[iru], m, &u[u_offset], ldu, (ftnlen)
			1);
#line 2171 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2171 "zgesdd.f"
		zunmbr_("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Set all of VT to identity matrix */

#line 2177 "zgesdd.f"
		zlaset_("F", n, n, &c_b1, &c_b2, &vt[vt_offset], ldvt, (
			ftnlen)1);

/*              Copy real matrix RWORK(IRVT) to complex matrix VT */
/*              Overwrite VT by right singular vectors of A */
/*              CWorkspace: need   2*M [tauq, taup] + N    [work] */
/*              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work] */
/*              RWorkspace: need   M [e] + M*M [RVT] */

#line 2185 "zgesdd.f"
		zlacp2_("F", m, m, &rwork[irvt], m, &vt[vt_offset], ldvt, (
			ftnlen)1);
#line 2186 "zgesdd.f"
		i__1 = *lwork - nwork + 1;
#line 2186 "zgesdd.f"
		zunmbr_("P", "R", "C", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 2189 "zgesdd.f"
	    }

#line 2191 "zgesdd.f"
	}

#line 2193 "zgesdd.f"
    }

/*     Undo scaling if necessary */

#line 2197 "zgesdd.f"
    if (iscl == 1) {
#line 2198 "zgesdd.f"
	if (anrm > bignum) {
#line 2198 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2198 "zgesdd.f"
	}
#line 2201 "zgesdd.f"
	if (*info != 0 && anrm > bignum) {
#line 2201 "zgesdd.f"
	    i__1 = minmn - 1;
#line 2201 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2201 "zgesdd.f"
	}
#line 2204 "zgesdd.f"
	if (anrm < smlnum) {
#line 2204 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr, (ftnlen)1);
#line 2204 "zgesdd.f"
	}
#line 2207 "zgesdd.f"
	if (*info != 0 && anrm < smlnum) {
#line 2207 "zgesdd.f"
	    i__1 = minmn - 1;
#line 2207 "zgesdd.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &i__1, &c__1, &rwork[
		    ie], &minmn, &ierr, (ftnlen)1);
#line 2207 "zgesdd.f"
	}
#line 2210 "zgesdd.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 2214 "zgesdd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 2216 "zgesdd.f"
    return 0;

/*     End of ZGESDD */

} /* zgesdd_ */

