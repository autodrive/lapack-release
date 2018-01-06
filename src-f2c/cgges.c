#line 1 "cgges.f"
/* cgges.f -- translated by f2c (version 20100827).
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

#line 1 "cgges.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> CGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f
or GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgges.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgges.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgges.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, */
/*                         SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, */
/*                         LWORK, RWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR, SORT */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), */
/*      $                   WORK( * ) */
/*       .. */
/*       .. Function Arguments .. */
/*       LOGICAL            SELCTG */
/*       EXTERNAL           SELCTG */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGES computes for a pair of N-by-N complex nonsymmetric matrices */
/* > (A,B), the generalized eigenvalues, the generalized complex Schur */
/* > form (S, T), and optionally left and/or right Schur vectors (VSL */
/* > and VSR). This gives the generalized Schur factorization */
/* > */
/* >         (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H ) */
/* > */
/* > where (VSR)**H is the conjugate-transpose of VSR. */
/* > */
/* > Optionally, it also orders the eigenvalues so that a selected cluster */
/* > of eigenvalues appears in the leading diagonal blocks of the upper */
/* > triangular matrix S and the upper triangular matrix T. The leading */
/* > columns of VSL and VSR then form an unitary basis for the */
/* > corresponding left and right eigenspaces (deflating subspaces). */
/* > */
/* > (If only the generalized eigenvalues are needed, use the driver */
/* > CGGEV instead, which is faster.) */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar w */
/* > or a ratio alpha/beta = w, such that  A - w*B is singular.  It is */
/* > usually represented as the pair (alpha,beta), as there is a */
/* > reasonable interpretation for beta=0, and even for both being zero. */
/* > */
/* > A pair of matrices (S,T) is in generalized complex Schur form if S */
/* > and T are upper triangular and, in addition, the diagonal elements */
/* > of T are non-negative real numbers. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBVSL */
/* > \verbatim */
/* >          JOBVSL is CHARACTER*1 */
/* >          = 'N':  do not compute the left Schur vectors; */
/* >          = 'V':  compute the left Schur vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVSR */
/* > \verbatim */
/* >          JOBVSR is CHARACTER*1 */
/* >          = 'N':  do not compute the right Schur vectors; */
/* >          = 'V':  compute the right Schur vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] SORT */
/* > \verbatim */
/* >          SORT is CHARACTER*1 */
/* >          Specifies whether or not to order the eigenvalues on the */
/* >          diagonal of the generalized Schur form. */
/* >          = 'N':  Eigenvalues are not ordered; */
/* >          = 'S':  Eigenvalues are ordered (see SELCTG). */
/* > \endverbatim */
/* > */
/* > \param[in] SELCTG */
/* > \verbatim */
/* >          SELCTG is a LOGICAL FUNCTION of two COMPLEX arguments */
/* >          SELCTG must be declared EXTERNAL in the calling subroutine. */
/* >          If SORT = 'N', SELCTG is not referenced. */
/* >          If SORT = 'S', SELCTG is used to select eigenvalues to sort */
/* >          to the top left of the Schur form. */
/* >          An eigenvalue ALPHA(j)/BETA(j) is selected if */
/* >          SELCTG(ALPHA(j),BETA(j)) is true. */
/* > */
/* >          Note that a selected complex eigenvalue may no longer satisfy */
/* >          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since */
/* >          ordering may change the value of complex eigenvalues */
/* >          (especially if the eigenvalue is ill-conditioned), in this */
/* >          case INFO is set to N+2 (See INFO below). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A, B, VSL, and VSR.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          On entry, the first of the pair of matrices. */
/* >          On exit, A has been overwritten by its generalized Schur */
/* >          form S. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB, N) */
/* >          On entry, the second of the pair of matrices. */
/* >          On exit, B has been overwritten by its generalized Schur */
/* >          form T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SDIM */
/* > \verbatim */
/* >          SDIM is INTEGER */
/* >          If SORT = 'N', SDIM = 0. */
/* >          If SORT = 'S', SDIM = number of eigenvalues (after sorting) */
/* >          for which SELCTG is true. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX array, dimension (N) */
/* >          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the */
/* >          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j), */
/* >          j=1,...,N  are the diagonals of the complex Schur form (A,B) */
/* >          output by CGGES. The  BETA(j) will be non-negative real. */
/* > */
/* >          Note: the quotients ALPHA(j)/BETA(j) may easily over- or */
/* >          underflow, and BETA(j) may even be zero.  Thus, the user */
/* >          should avoid naively computing the ratio alpha/beta. */
/* >          However, ALPHA will be always less than and usually */
/* >          comparable with norm(A) in magnitude, and BETA always less */
/* >          than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* >          VSL is COMPLEX array, dimension (LDVSL,N) */
/* >          If JOBVSL = 'V', VSL will contain the left Schur vectors. */
/* >          Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* >          LDVSL is INTEGER */
/* >          The leading dimension of the matrix VSL. LDVSL >= 1, and */
/* >          if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* >          VSR is COMPLEX array, dimension (LDVSR,N) */
/* >          If JOBVSR = 'V', VSR will contain the right Schur vectors. */
/* >          Not referenced if JOBVSR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSR */
/* > \verbatim */
/* >          LDVSR is INTEGER */
/* >          The leading dimension of the matrix VSR. LDVSR >= 1, and */
/* >          if JOBVSR = 'V', LDVSR >= N. */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,2*N). */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (8*N) */
/* > \endverbatim */
/* > */
/* > \param[out] BWORK */
/* > \verbatim */
/* >          BWORK is LOGICAL array, dimension (N) */
/* >          Not referenced if SORT = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          =1,...,N: */
/* >                The QZ iteration failed.  (A,B) are not in Schur */
/* >                form, but ALPHA(j) and BETA(j) should be correct for */
/* >                j=INFO+1,...,N. */
/* >          > N:  =N+1: other than QZ iteration failed in CHGEQZ */
/* >                =N+2: after reordering, roundoff changed values of */
/* >                      some complex eigenvalues so that leading */
/* >                      eigenvalues in the Generalized Schur form no */
/* >                      longer satisfy SELCTG=.TRUE.  This could also */
/* >                      be caused due to scaling. */
/* >                =N+3: reordering failed in CTGSEN. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexGEeigen */

/*  ===================================================================== */
/* Subroutine */ int cgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer 
	*ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, 
	logical *bwork, integer *info, ftnlen jobvsl_len, ftnlen jobvsr_len, 
	ftnlen sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal dif[2];
    static integer ihi, ilo;
    static doublereal eps, anrm, bnrm;
    static integer idum[1], ierr, itau, iwrk;
    static doublereal pvsl, pvsr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ileft, icols;
    static logical cursl, ilvsl, ilvsr;
    static integer irwrk, irows;
    extern /* Subroutine */ int cggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), cggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), slabad_(doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int cgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int chgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen, ftnlen), ctgsen_(
	    integer *, logical *, logical *, logical *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublecomplex *, integer *, integer *, integer *, 
	    integer *);
    static integer ijobvl, iright, ijobvr;
    static doublereal anrmto;
    static integer lwkmin;
    static logical lastsl;
    static doublereal bnrmto;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), cunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static doublereal smlnum;
    static logical wantst, lquery;
    static integer lwkopt;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Function Arguments .. */
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

/*     Decode the input arguments */

#line 334 "cgges.f"
    /* Parameter adjustments */
#line 334 "cgges.f"
    a_dim1 = *lda;
#line 334 "cgges.f"
    a_offset = 1 + a_dim1;
#line 334 "cgges.f"
    a -= a_offset;
#line 334 "cgges.f"
    b_dim1 = *ldb;
#line 334 "cgges.f"
    b_offset = 1 + b_dim1;
#line 334 "cgges.f"
    b -= b_offset;
#line 334 "cgges.f"
    --alpha;
#line 334 "cgges.f"
    --beta;
#line 334 "cgges.f"
    vsl_dim1 = *ldvsl;
#line 334 "cgges.f"
    vsl_offset = 1 + vsl_dim1;
#line 334 "cgges.f"
    vsl -= vsl_offset;
#line 334 "cgges.f"
    vsr_dim1 = *ldvsr;
#line 334 "cgges.f"
    vsr_offset = 1 + vsr_dim1;
#line 334 "cgges.f"
    vsr -= vsr_offset;
#line 334 "cgges.f"
    --work;
#line 334 "cgges.f"
    --rwork;
#line 334 "cgges.f"
    --bwork;
#line 334 "cgges.f"

#line 334 "cgges.f"
    /* Function Body */
#line 334 "cgges.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 335 "cgges.f"
	ijobvl = 1;
#line 336 "cgges.f"
	ilvsl = FALSE_;
#line 337 "cgges.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 338 "cgges.f"
	ijobvl = 2;
#line 339 "cgges.f"
	ilvsl = TRUE_;
#line 340 "cgges.f"
    } else {
#line 341 "cgges.f"
	ijobvl = -1;
#line 342 "cgges.f"
	ilvsl = FALSE_;
#line 343 "cgges.f"
    }

#line 345 "cgges.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 346 "cgges.f"
	ijobvr = 1;
#line 347 "cgges.f"
	ilvsr = FALSE_;
#line 348 "cgges.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 349 "cgges.f"
	ijobvr = 2;
#line 350 "cgges.f"
	ilvsr = TRUE_;
#line 351 "cgges.f"
    } else {
#line 352 "cgges.f"
	ijobvr = -1;
#line 353 "cgges.f"
	ilvsr = FALSE_;
#line 354 "cgges.f"
    }

#line 356 "cgges.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 360 "cgges.f"
    *info = 0;
#line 361 "cgges.f"
    lquery = *lwork == -1;
#line 362 "cgges.f"
    if (ijobvl <= 0) {
#line 363 "cgges.f"
	*info = -1;
#line 364 "cgges.f"
    } else if (ijobvr <= 0) {
#line 365 "cgges.f"
	*info = -2;
#line 366 "cgges.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 367 "cgges.f"
	*info = -3;
#line 368 "cgges.f"
    } else if (*n < 0) {
#line 369 "cgges.f"
	*info = -5;
#line 370 "cgges.f"
    } else if (*lda < max(1,*n)) {
#line 371 "cgges.f"
	*info = -7;
#line 372 "cgges.f"
    } else if (*ldb < max(1,*n)) {
#line 373 "cgges.f"
	*info = -9;
#line 374 "cgges.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 375 "cgges.f"
	*info = -14;
#line 376 "cgges.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 377 "cgges.f"
	*info = -16;
#line 378 "cgges.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 387 "cgges.f"
    if (*info == 0) {
/* Computing MAX */
#line 388 "cgges.f"
	i__1 = 1, i__2 = *n << 1;
#line 388 "cgges.f"
	lwkmin = max(i__1,i__2);
/* Computing MAX */
#line 389 "cgges.f"
	i__1 = 1, i__2 = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", n, &c__1, n, 
		&c__0, (ftnlen)6, (ftnlen)1);
#line 389 "cgges.f"
	lwkopt = max(i__1,i__2);
/* Computing MAX */
#line 390 "cgges.f"
	i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "CUNMQR", " ", n, &
		c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 390 "cgges.f"
	lwkopt = max(i__1,i__2);
#line 392 "cgges.f"
	if (ilvsl) {
/* Computing MAX */
#line 393 "cgges.f"
	    i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "CUNGQR", " ", n, &
		    c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 393 "cgges.f"
	    lwkopt = max(i__1,i__2);
#line 395 "cgges.f"
	}
#line 396 "cgges.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 398 "cgges.f"
	if (*lwork < lwkmin && ! lquery) {
#line 398 "cgges.f"
	    *info = -18;
#line 398 "cgges.f"
	}
#line 400 "cgges.f"
    }

#line 402 "cgges.f"
    if (*info != 0) {
#line 403 "cgges.f"
	i__1 = -(*info);
#line 403 "cgges.f"
	xerbla_("CGGES ", &i__1, (ftnlen)6);
#line 404 "cgges.f"
	return 0;
#line 405 "cgges.f"
    } else if (lquery) {
#line 406 "cgges.f"
	return 0;
#line 407 "cgges.f"
    }

/*     Quick return if possible */

#line 411 "cgges.f"
    if (*n == 0) {
#line 412 "cgges.f"
	*sdim = 0;
#line 413 "cgges.f"
	return 0;
#line 414 "cgges.f"
    }

/*     Get machine constants */

#line 418 "cgges.f"
    eps = slamch_("P", (ftnlen)1);
#line 419 "cgges.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 420 "cgges.f"
    bignum = 1. / smlnum;
#line 421 "cgges.f"
    slabad_(&smlnum, &bignum);
#line 422 "cgges.f"
    smlnum = sqrt(smlnum) / eps;
#line 423 "cgges.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 427 "cgges.f"
    anrm = clange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 428 "cgges.f"
    ilascl = FALSE_;
#line 429 "cgges.f"
    if (anrm > 0. && anrm < smlnum) {
#line 430 "cgges.f"
	anrmto = smlnum;
#line 431 "cgges.f"
	ilascl = TRUE_;
#line 432 "cgges.f"
    } else if (anrm > bignum) {
#line 433 "cgges.f"
	anrmto = bignum;
#line 434 "cgges.f"
	ilascl = TRUE_;
#line 435 "cgges.f"
    }

#line 437 "cgges.f"
    if (ilascl) {
#line 437 "cgges.f"
	clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 437 "cgges.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 442 "cgges.f"
    bnrm = clange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 443 "cgges.f"
    ilbscl = FALSE_;
#line 444 "cgges.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 445 "cgges.f"
	bnrmto = smlnum;
#line 446 "cgges.f"
	ilbscl = TRUE_;
#line 447 "cgges.f"
    } else if (bnrm > bignum) {
#line 448 "cgges.f"
	bnrmto = bignum;
#line 449 "cgges.f"
	ilbscl = TRUE_;
#line 450 "cgges.f"
    }

#line 452 "cgges.f"
    if (ilbscl) {
#line 452 "cgges.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 452 "cgges.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Real Workspace: need 6*N) */

#line 458 "cgges.f"
    ileft = 1;
#line 459 "cgges.f"
    iright = *n + 1;
#line 460 "cgges.f"
    irwrk = iright + *n;
#line 461 "cgges.f"
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 467 "cgges.f"
    irows = ihi + 1 - ilo;
#line 468 "cgges.f"
    icols = *n + 1 - ilo;
#line 469 "cgges.f"
    itau = 1;
#line 470 "cgges.f"
    iwrk = itau + irows;
#line 471 "cgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 471 "cgges.f"
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 477 "cgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 477 "cgges.f"
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 484 "cgges.f"
    if (ilvsl) {
#line 485 "cgges.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 486 "cgges.f"
	if (irows > 1) {
#line 487 "cgges.f"
	    i__1 = irows - 1;
#line 487 "cgges.f"
	    i__2 = irows - 1;
#line 487 "cgges.f"
	    clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 489 "cgges.f"
	}
#line 490 "cgges.f"
	i__1 = *lwork + 1 - iwrk;
#line 490 "cgges.f"
	cungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 492 "cgges.f"
    }

/*     Initialize VSR */

#line 496 "cgges.f"
    if (ilvsr) {
#line 496 "cgges.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 496 "cgges.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 502 "cgges.f"
    cgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 505 "cgges.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Complex Workspace: need N) */
/*     (Real Workspace: need N) */

#line 511 "cgges.f"
    iwrk = itau;
#line 512 "cgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 512 "cgges.f"
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 515 "cgges.f"
    if (ierr != 0) {
#line 516 "cgges.f"
	if (ierr > 0 && ierr <= *n) {
#line 517 "cgges.f"
	    *info = ierr;
#line 518 "cgges.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 519 "cgges.f"
	    *info = ierr - *n;
#line 520 "cgges.f"
	} else {
#line 521 "cgges.f"
	    *info = *n + 1;
#line 522 "cgges.f"
	}
#line 523 "cgges.f"
	goto L30;
#line 524 "cgges.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */
/*     (Workspace: none needed) */

#line 529 "cgges.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before selecting */

#line 533 "cgges.f"
	if (ilascl) {
#line 533 "cgges.f"
	    clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, &c__1, &alpha[1], n,
		     &ierr, (ftnlen)1);
#line 533 "cgges.f"
	}
#line 535 "cgges.f"
	if (ilbscl) {
#line 535 "cgges.f"
	    clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 535 "cgges.f"
	}

/*        Select eigenvalues */

#line 540 "cgges.f"
	i__1 = *n;
#line 540 "cgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 541 "cgges.f"
	    bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
#line 542 "cgges.f"
/* L10: */
#line 542 "cgges.f"
	}

#line 544 "cgges.f"
	i__1 = *lwork - iwrk + 1;
#line 544 "cgges.f"
	ctgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, sdim, &pvsl, &pvsr, dif, &work[iwrk],
		 &i__1, idum, &c__1, &ierr);
#line 547 "cgges.f"
	if (ierr == 1) {
#line 547 "cgges.f"
	    *info = *n + 3;
#line 547 "cgges.f"
	}

#line 550 "cgges.f"
    }

/*     Apply back-permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 555 "cgges.f"
    if (ilvsl) {
#line 555 "cgges.f"
	cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 555 "cgges.f"
    }
#line 558 "cgges.f"
    if (ilvsr) {
#line 558 "cgges.f"
	cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 558 "cgges.f"
    }

/*     Undo scaling */

#line 564 "cgges.f"
    if (ilascl) {
#line 565 "cgges.f"
	clascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 566 "cgges.f"
	clascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 567 "cgges.f"
    }

#line 569 "cgges.f"
    if (ilbscl) {
#line 570 "cgges.f"
	clascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 571 "cgges.f"
	clascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 572 "cgges.f"
    }

#line 574 "cgges.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 578 "cgges.f"
	lastsl = TRUE_;
#line 579 "cgges.f"
	*sdim = 0;
#line 580 "cgges.f"
	i__1 = *n;
#line 580 "cgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 581 "cgges.f"
	    cursl = (*selctg)(&alpha[i__], &beta[i__]);
#line 582 "cgges.f"
	    if (cursl) {
#line 582 "cgges.f"
		++(*sdim);
#line 582 "cgges.f"
	    }
#line 584 "cgges.f"
	    if (cursl && ! lastsl) {
#line 584 "cgges.f"
		*info = *n + 2;
#line 584 "cgges.f"
	    }
#line 586 "cgges.f"
	    lastsl = cursl;
#line 587 "cgges.f"
/* L20: */
#line 587 "cgges.f"
	}

#line 589 "cgges.f"
    }

#line 591 "cgges.f"
L30:

#line 593 "cgges.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 595 "cgges.f"
    return 0;

/*     End of CGGES */

} /* cgges_ */

