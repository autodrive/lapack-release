#line 1 "zgges.f"
/* zgges.f -- translated by f2c (version 20100827).
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

#line 1 "zgges.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> ZGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f
or GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgges.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgges.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgges.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, */
/*                         SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, */
/*                         LWORK, RWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR, SORT */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), */
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
/* > ZGGES computes for a pair of N-by-N complex nonsymmetric matrices */
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
/* > ZGGEV instead, which is faster.) */
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
/* >          SELCTG is a LOGICAL FUNCTION of two COMPLEX*16 arguments */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB, N) */
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
/* >          ALPHA is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 array, dimension (N) */
/* >          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the */
/* >          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j), */
/* >          j=1,...,N  are the diagonals of the complex Schur form (A,B) */
/* >          output by ZGGES. The  BETA(j) will be non-negative real. */
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
/* >          VSL is COMPLEX*16 array, dimension (LDVSL,N) */
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
/* >          VSR is COMPLEX*16 array, dimension (LDVSR,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (8*N) */
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
/* >          > N:  =N+1: other than QZ iteration failed in ZHGEQZ */
/* >                =N+2: after reordering, roundoff changed values of */
/* >                      some complex eigenvalues so that leading */
/* >                      eigenvalues in the Generalized Schur form no */
/* >                      longer satisfy SELCTG=.TRUE.  This could also */
/* >                      be caused due to scaling. */
/* >                =N+3: reordering failed in ZTGSEN. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), zggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    static integer ijobvl, iright;
    extern /* Subroutine */ int zgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static doublereal anrmto;
    static integer lwkmin;
    static logical lastsl;
    static doublereal bnrmto;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), zhgeqz_(
	    char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen), ztgsen_(integer 
	    *, logical *, logical *, logical *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, integer *, integer *, integer *);
    static doublereal smlnum;
    static logical wantst, lquery;
    static integer lwkopt;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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

#line 334 "zgges.f"
    /* Parameter adjustments */
#line 334 "zgges.f"
    a_dim1 = *lda;
#line 334 "zgges.f"
    a_offset = 1 + a_dim1;
#line 334 "zgges.f"
    a -= a_offset;
#line 334 "zgges.f"
    b_dim1 = *ldb;
#line 334 "zgges.f"
    b_offset = 1 + b_dim1;
#line 334 "zgges.f"
    b -= b_offset;
#line 334 "zgges.f"
    --alpha;
#line 334 "zgges.f"
    --beta;
#line 334 "zgges.f"
    vsl_dim1 = *ldvsl;
#line 334 "zgges.f"
    vsl_offset = 1 + vsl_dim1;
#line 334 "zgges.f"
    vsl -= vsl_offset;
#line 334 "zgges.f"
    vsr_dim1 = *ldvsr;
#line 334 "zgges.f"
    vsr_offset = 1 + vsr_dim1;
#line 334 "zgges.f"
    vsr -= vsr_offset;
#line 334 "zgges.f"
    --work;
#line 334 "zgges.f"
    --rwork;
#line 334 "zgges.f"
    --bwork;
#line 334 "zgges.f"

#line 334 "zgges.f"
    /* Function Body */
#line 334 "zgges.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 335 "zgges.f"
	ijobvl = 1;
#line 336 "zgges.f"
	ilvsl = FALSE_;
#line 337 "zgges.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 338 "zgges.f"
	ijobvl = 2;
#line 339 "zgges.f"
	ilvsl = TRUE_;
#line 340 "zgges.f"
    } else {
#line 341 "zgges.f"
	ijobvl = -1;
#line 342 "zgges.f"
	ilvsl = FALSE_;
#line 343 "zgges.f"
    }

#line 345 "zgges.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 346 "zgges.f"
	ijobvr = 1;
#line 347 "zgges.f"
	ilvsr = FALSE_;
#line 348 "zgges.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 349 "zgges.f"
	ijobvr = 2;
#line 350 "zgges.f"
	ilvsr = TRUE_;
#line 351 "zgges.f"
    } else {
#line 352 "zgges.f"
	ijobvr = -1;
#line 353 "zgges.f"
	ilvsr = FALSE_;
#line 354 "zgges.f"
    }

#line 356 "zgges.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 360 "zgges.f"
    *info = 0;
#line 361 "zgges.f"
    lquery = *lwork == -1;
#line 362 "zgges.f"
    if (ijobvl <= 0) {
#line 363 "zgges.f"
	*info = -1;
#line 364 "zgges.f"
    } else if (ijobvr <= 0) {
#line 365 "zgges.f"
	*info = -2;
#line 366 "zgges.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 367 "zgges.f"
	*info = -3;
#line 368 "zgges.f"
    } else if (*n < 0) {
#line 369 "zgges.f"
	*info = -5;
#line 370 "zgges.f"
    } else if (*lda < max(1,*n)) {
#line 371 "zgges.f"
	*info = -7;
#line 372 "zgges.f"
    } else if (*ldb < max(1,*n)) {
#line 373 "zgges.f"
	*info = -9;
#line 374 "zgges.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 375 "zgges.f"
	*info = -14;
#line 376 "zgges.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 377 "zgges.f"
	*info = -16;
#line 378 "zgges.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 387 "zgges.f"
    if (*info == 0) {
/* Computing MAX */
#line 388 "zgges.f"
	i__1 = 1, i__2 = *n << 1;
#line 388 "zgges.f"
	lwkmin = max(i__1,i__2);
/* Computing MAX */
#line 389 "zgges.f"
	i__1 = 1, i__2 = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", n, &c__1, n, 
		&c__0, (ftnlen)6, (ftnlen)1);
#line 389 "zgges.f"
	lwkopt = max(i__1,i__2);
/* Computing MAX */
#line 390 "zgges.f"
	i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "ZUNMQR", " ", n, &
		c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 390 "zgges.f"
	lwkopt = max(i__1,i__2);
#line 392 "zgges.f"
	if (ilvsl) {
/* Computing MAX */
#line 393 "zgges.f"
	    i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", " ", n, &
		    c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 393 "zgges.f"
	    lwkopt = max(i__1,i__2);
#line 395 "zgges.f"
	}
#line 396 "zgges.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 398 "zgges.f"
	if (*lwork < lwkmin && ! lquery) {
#line 398 "zgges.f"
	    *info = -18;
#line 398 "zgges.f"
	}
#line 400 "zgges.f"
    }

#line 402 "zgges.f"
    if (*info != 0) {
#line 403 "zgges.f"
	i__1 = -(*info);
#line 403 "zgges.f"
	xerbla_("ZGGES ", &i__1, (ftnlen)6);
#line 404 "zgges.f"
	return 0;
#line 405 "zgges.f"
    } else if (lquery) {
#line 406 "zgges.f"
	return 0;
#line 407 "zgges.f"
    }

/*     Quick return if possible */

#line 411 "zgges.f"
    if (*n == 0) {
#line 412 "zgges.f"
	*sdim = 0;
#line 413 "zgges.f"
	return 0;
#line 414 "zgges.f"
    }

/*     Get machine constants */

#line 418 "zgges.f"
    eps = dlamch_("P", (ftnlen)1);
#line 419 "zgges.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 420 "zgges.f"
    bignum = 1. / smlnum;
#line 421 "zgges.f"
    dlabad_(&smlnum, &bignum);
#line 422 "zgges.f"
    smlnum = sqrt(smlnum) / eps;
#line 423 "zgges.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 427 "zgges.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 428 "zgges.f"
    ilascl = FALSE_;
#line 429 "zgges.f"
    if (anrm > 0. && anrm < smlnum) {
#line 430 "zgges.f"
	anrmto = smlnum;
#line 431 "zgges.f"
	ilascl = TRUE_;
#line 432 "zgges.f"
    } else if (anrm > bignum) {
#line 433 "zgges.f"
	anrmto = bignum;
#line 434 "zgges.f"
	ilascl = TRUE_;
#line 435 "zgges.f"
    }

#line 437 "zgges.f"
    if (ilascl) {
#line 437 "zgges.f"
	zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 437 "zgges.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 442 "zgges.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 443 "zgges.f"
    ilbscl = FALSE_;
#line 444 "zgges.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 445 "zgges.f"
	bnrmto = smlnum;
#line 446 "zgges.f"
	ilbscl = TRUE_;
#line 447 "zgges.f"
    } else if (bnrm > bignum) {
#line 448 "zgges.f"
	bnrmto = bignum;
#line 449 "zgges.f"
	ilbscl = TRUE_;
#line 450 "zgges.f"
    }

#line 452 "zgges.f"
    if (ilbscl) {
#line 452 "zgges.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 452 "zgges.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Real Workspace: need 6*N) */

#line 458 "zgges.f"
    ileft = 1;
#line 459 "zgges.f"
    iright = *n + 1;
#line 460 "zgges.f"
    irwrk = iright + *n;
#line 461 "zgges.f"
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 467 "zgges.f"
    irows = ihi + 1 - ilo;
#line 468 "zgges.f"
    icols = *n + 1 - ilo;
#line 469 "zgges.f"
    itau = 1;
#line 470 "zgges.f"
    iwrk = itau + irows;
#line 471 "zgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 471 "zgges.f"
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 477 "zgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 477 "zgges.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 484 "zgges.f"
    if (ilvsl) {
#line 485 "zgges.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 486 "zgges.f"
	if (irows > 1) {
#line 487 "zgges.f"
	    i__1 = irows - 1;
#line 487 "zgges.f"
	    i__2 = irows - 1;
#line 487 "zgges.f"
	    zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 489 "zgges.f"
	}
#line 490 "zgges.f"
	i__1 = *lwork + 1 - iwrk;
#line 490 "zgges.f"
	zungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 492 "zgges.f"
    }

/*     Initialize VSR */

#line 496 "zgges.f"
    if (ilvsr) {
#line 496 "zgges.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 496 "zgges.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 502 "zgges.f"
    zgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 505 "zgges.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Complex Workspace: need N) */
/*     (Real Workspace: need N) */

#line 511 "zgges.f"
    iwrk = itau;
#line 512 "zgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 512 "zgges.f"
    zhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 515 "zgges.f"
    if (ierr != 0) {
#line 516 "zgges.f"
	if (ierr > 0 && ierr <= *n) {
#line 517 "zgges.f"
	    *info = ierr;
#line 518 "zgges.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 519 "zgges.f"
	    *info = ierr - *n;
#line 520 "zgges.f"
	} else {
#line 521 "zgges.f"
	    *info = *n + 1;
#line 522 "zgges.f"
	}
#line 523 "zgges.f"
	goto L30;
#line 524 "zgges.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */
/*     (Workspace: none needed) */

#line 529 "zgges.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before selecting */

#line 533 "zgges.f"
	if (ilascl) {
#line 533 "zgges.f"
	    zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, &c__1, &alpha[1], n,
		     &ierr, (ftnlen)1);
#line 533 "zgges.f"
	}
#line 535 "zgges.f"
	if (ilbscl) {
#line 535 "zgges.f"
	    zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 535 "zgges.f"
	}

/*        Select eigenvalues */

#line 540 "zgges.f"
	i__1 = *n;
#line 540 "zgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 541 "zgges.f"
	    bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
#line 542 "zgges.f"
/* L10: */
#line 542 "zgges.f"
	}

#line 544 "zgges.f"
	i__1 = *lwork - iwrk + 1;
#line 544 "zgges.f"
	ztgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, sdim, &pvsl, &pvsr, dif, &work[iwrk],
		 &i__1, idum, &c__1, &ierr);
#line 547 "zgges.f"
	if (ierr == 1) {
#line 547 "zgges.f"
	    *info = *n + 3;
#line 547 "zgges.f"
	}

#line 550 "zgges.f"
    }

/*     Apply back-permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 555 "zgges.f"
    if (ilvsl) {
#line 555 "zgges.f"
	zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 555 "zgges.f"
    }
#line 558 "zgges.f"
    if (ilvsr) {
#line 558 "zgges.f"
	zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 558 "zgges.f"
    }

/*     Undo scaling */

#line 564 "zgges.f"
    if (ilascl) {
#line 565 "zgges.f"
	zlascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 566 "zgges.f"
	zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 567 "zgges.f"
    }

#line 569 "zgges.f"
    if (ilbscl) {
#line 570 "zgges.f"
	zlascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 571 "zgges.f"
	zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 572 "zgges.f"
    }

#line 574 "zgges.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 578 "zgges.f"
	lastsl = TRUE_;
#line 579 "zgges.f"
	*sdim = 0;
#line 580 "zgges.f"
	i__1 = *n;
#line 580 "zgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 581 "zgges.f"
	    cursl = (*selctg)(&alpha[i__], &beta[i__]);
#line 582 "zgges.f"
	    if (cursl) {
#line 582 "zgges.f"
		++(*sdim);
#line 582 "zgges.f"
	    }
#line 584 "zgges.f"
	    if (cursl && ! lastsl) {
#line 584 "zgges.f"
		*info = *n + 2;
#line 584 "zgges.f"
	    }
#line 586 "zgges.f"
	    lastsl = cursl;
#line 587 "zgges.f"
/* L20: */
#line 587 "zgges.f"
	}

#line 589 "zgges.f"
    }

#line 591 "zgges.f"
L30:

#line 593 "zgges.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 595 "zgges.f"
    return 0;

/*     End of ZGGES */

} /* zgges_ */

