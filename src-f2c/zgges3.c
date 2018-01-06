#line 1 "zgges3.f"
/* zgges3.f -- translated by f2c (version 20100827).
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

#line 1 "zgges3.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;

/* > \brief <b> ZGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices (blocked algorithm)</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGES3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgges3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgges3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgges3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, */
/*      $                   LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, */
/*      $                   WORK, LWORK, RWORK, BWORK, INFO ) */

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
/* > ZGGES3 computes for a pair of N-by-N complex nonsymmetric matrices */
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
/* >          output by ZGGES3. The  BETA(j) will be non-negative real. */
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
/* >          The dimension of the array WORK. */
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

/* > \date January 2015 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zgges3_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
    doublecomplex z__1;

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
    extern /* Subroutine */ int zgghd3_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), zggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    static integer ijobvl, iright;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static doublereal anrmto, bnrmto;
    static logical lastsl;
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


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     January 2015 */

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

#line 331 "zgges3.f"
    /* Parameter adjustments */
#line 331 "zgges3.f"
    a_dim1 = *lda;
#line 331 "zgges3.f"
    a_offset = 1 + a_dim1;
#line 331 "zgges3.f"
    a -= a_offset;
#line 331 "zgges3.f"
    b_dim1 = *ldb;
#line 331 "zgges3.f"
    b_offset = 1 + b_dim1;
#line 331 "zgges3.f"
    b -= b_offset;
#line 331 "zgges3.f"
    --alpha;
#line 331 "zgges3.f"
    --beta;
#line 331 "zgges3.f"
    vsl_dim1 = *ldvsl;
#line 331 "zgges3.f"
    vsl_offset = 1 + vsl_dim1;
#line 331 "zgges3.f"
    vsl -= vsl_offset;
#line 331 "zgges3.f"
    vsr_dim1 = *ldvsr;
#line 331 "zgges3.f"
    vsr_offset = 1 + vsr_dim1;
#line 331 "zgges3.f"
    vsr -= vsr_offset;
#line 331 "zgges3.f"
    --work;
#line 331 "zgges3.f"
    --rwork;
#line 331 "zgges3.f"
    --bwork;
#line 331 "zgges3.f"

#line 331 "zgges3.f"
    /* Function Body */
#line 331 "zgges3.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 332 "zgges3.f"
	ijobvl = 1;
#line 333 "zgges3.f"
	ilvsl = FALSE_;
#line 334 "zgges3.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 335 "zgges3.f"
	ijobvl = 2;
#line 336 "zgges3.f"
	ilvsl = TRUE_;
#line 337 "zgges3.f"
    } else {
#line 338 "zgges3.f"
	ijobvl = -1;
#line 339 "zgges3.f"
	ilvsl = FALSE_;
#line 340 "zgges3.f"
    }

#line 342 "zgges3.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 343 "zgges3.f"
	ijobvr = 1;
#line 344 "zgges3.f"
	ilvsr = FALSE_;
#line 345 "zgges3.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 346 "zgges3.f"
	ijobvr = 2;
#line 347 "zgges3.f"
	ilvsr = TRUE_;
#line 348 "zgges3.f"
    } else {
#line 349 "zgges3.f"
	ijobvr = -1;
#line 350 "zgges3.f"
	ilvsr = FALSE_;
#line 351 "zgges3.f"
    }

#line 353 "zgges3.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 357 "zgges3.f"
    *info = 0;
#line 358 "zgges3.f"
    lquery = *lwork == -1;
#line 359 "zgges3.f"
    if (ijobvl <= 0) {
#line 360 "zgges3.f"
	*info = -1;
#line 361 "zgges3.f"
    } else if (ijobvr <= 0) {
#line 362 "zgges3.f"
	*info = -2;
#line 363 "zgges3.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 364 "zgges3.f"
	*info = -3;
#line 365 "zgges3.f"
    } else if (*n < 0) {
#line 366 "zgges3.f"
	*info = -5;
#line 367 "zgges3.f"
    } else if (*lda < max(1,*n)) {
#line 368 "zgges3.f"
	*info = -7;
#line 369 "zgges3.f"
    } else if (*ldb < max(1,*n)) {
#line 370 "zgges3.f"
	*info = -9;
#line 371 "zgges3.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 372 "zgges3.f"
	*info = -14;
#line 373 "zgges3.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 374 "zgges3.f"
	*info = -16;
#line 375 "zgges3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 375 "zgges3.f"
	i__1 = 1, i__2 = *n << 1;
#line 375 "zgges3.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 376 "zgges3.f"
	    *info = -18;
#line 377 "zgges3.f"
	}
#line 377 "zgges3.f"
    }

/*     Compute workspace */

#line 381 "zgges3.f"
    if (*info == 0) {
#line 382 "zgges3.f"
	zgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
/* Computing MAX */
#line 383 "zgges3.f"
	i__1 = 1, i__2 = *n + (integer) work[1].r;
#line 383 "zgges3.f"
	lwkopt = max(i__1,i__2);
#line 384 "zgges3.f"
	zunmqr_("L", "C", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], 
		lda, &work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 386 "zgges3.f"
	i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 386 "zgges3.f"
	lwkopt = max(i__1,i__2);
#line 387 "zgges3.f"
	if (ilvsl) {
#line 388 "zgges3.f"
	    zungqr_(n, n, n, &vsl[vsl_offset], ldvsl, &work[1], &work[1], &
		    c_n1, &ierr);
/* Computing MAX */
#line 389 "zgges3.f"
	    i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 389 "zgges3.f"
	    lwkopt = max(i__1,i__2);
#line 390 "zgges3.f"
	}
#line 391 "zgges3.f"
	zgghd3_(jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[
		1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 393 "zgges3.f"
	i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 393 "zgges3.f"
	lwkopt = max(i__1,i__2);
#line 394 "zgges3.f"
	zhgeqz_("S", jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, &work[1], &c_n1, &rwork[1], &ierr, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 397 "zgges3.f"
	i__1 = lwkopt, i__2 = (integer) work[1].r;
#line 397 "zgges3.f"
	lwkopt = max(i__1,i__2);
#line 398 "zgges3.f"
	if (wantst) {
#line 399 "zgges3.f"
	    ztgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &
		    b[b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], 
		    ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, &pvsr, dif, &
		    work[1], &c_n1, idum, &c__1, &ierr);
/* Computing MAX */
#line 402 "zgges3.f"
	    i__1 = lwkopt, i__2 = (integer) work[1].r;
#line 402 "zgges3.f"
	    lwkopt = max(i__1,i__2);
#line 403 "zgges3.f"
	}
#line 404 "zgges3.f"
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 404 "zgges3.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 405 "zgges3.f"
    }

#line 407 "zgges3.f"
    if (*info != 0) {
#line 408 "zgges3.f"
	i__1 = -(*info);
#line 408 "zgges3.f"
	xerbla_("ZGGES3 ", &i__1, (ftnlen)7);
#line 409 "zgges3.f"
	return 0;
#line 410 "zgges3.f"
    } else if (lquery) {
#line 411 "zgges3.f"
	return 0;
#line 412 "zgges3.f"
    }

/*     Quick return if possible */

#line 416 "zgges3.f"
    if (*n == 0) {
#line 417 "zgges3.f"
	*sdim = 0;
#line 418 "zgges3.f"
	return 0;
#line 419 "zgges3.f"
    }

/*     Get machine constants */

#line 423 "zgges3.f"
    eps = dlamch_("P", (ftnlen)1);
#line 424 "zgges3.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 425 "zgges3.f"
    bignum = 1. / smlnum;
#line 426 "zgges3.f"
    dlabad_(&smlnum, &bignum);
#line 427 "zgges3.f"
    smlnum = sqrt(smlnum) / eps;
#line 428 "zgges3.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 432 "zgges3.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 433 "zgges3.f"
    ilascl = FALSE_;
#line 434 "zgges3.f"
    if (anrm > 0. && anrm < smlnum) {
#line 435 "zgges3.f"
	anrmto = smlnum;
#line 436 "zgges3.f"
	ilascl = TRUE_;
#line 437 "zgges3.f"
    } else if (anrm > bignum) {
#line 438 "zgges3.f"
	anrmto = bignum;
#line 439 "zgges3.f"
	ilascl = TRUE_;
#line 440 "zgges3.f"
    }

#line 442 "zgges3.f"
    if (ilascl) {
#line 442 "zgges3.f"
	zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 442 "zgges3.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 447 "zgges3.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 448 "zgges3.f"
    ilbscl = FALSE_;
#line 449 "zgges3.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 450 "zgges3.f"
	bnrmto = smlnum;
#line 451 "zgges3.f"
	ilbscl = TRUE_;
#line 452 "zgges3.f"
    } else if (bnrm > bignum) {
#line 453 "zgges3.f"
	bnrmto = bignum;
#line 454 "zgges3.f"
	ilbscl = TRUE_;
#line 455 "zgges3.f"
    }

#line 457 "zgges3.f"
    if (ilbscl) {
#line 457 "zgges3.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 457 "zgges3.f"
    }

/*     Permute the matrix to make it more nearly triangular */

#line 462 "zgges3.f"
    ileft = 1;
#line 463 "zgges3.f"
    iright = *n + 1;
#line 464 "zgges3.f"
    irwrk = iright + *n;
#line 465 "zgges3.f"
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */

#line 470 "zgges3.f"
    irows = ihi + 1 - ilo;
#line 471 "zgges3.f"
    icols = *n + 1 - ilo;
#line 472 "zgges3.f"
    itau = 1;
#line 473 "zgges3.f"
    iwrk = itau + irows;
#line 474 "zgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 474 "zgges3.f"
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */

#line 479 "zgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 479 "zgges3.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */

#line 485 "zgges3.f"
    if (ilvsl) {
#line 486 "zgges3.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 487 "zgges3.f"
	if (irows > 1) {
#line 488 "zgges3.f"
	    i__1 = irows - 1;
#line 488 "zgges3.f"
	    i__2 = irows - 1;
#line 488 "zgges3.f"
	    zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 490 "zgges3.f"
	}
#line 491 "zgges3.f"
	i__1 = *lwork + 1 - iwrk;
#line 491 "zgges3.f"
	zungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 493 "zgges3.f"
    }

/*     Initialize VSR */

#line 497 "zgges3.f"
    if (ilvsr) {
#line 497 "zgges3.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 497 "zgges3.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 502 "zgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 502 "zgges3.f"
    zgghd3_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk]
	    , &i__1, &ierr, (ftnlen)1, (ftnlen)1);

#line 505 "zgges3.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */

#line 509 "zgges3.f"
    iwrk = itau;
#line 510 "zgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 510 "zgges3.f"
    zhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 513 "zgges3.f"
    if (ierr != 0) {
#line 514 "zgges3.f"
	if (ierr > 0 && ierr <= *n) {
#line 515 "zgges3.f"
	    *info = ierr;
#line 516 "zgges3.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 517 "zgges3.f"
	    *info = ierr - *n;
#line 518 "zgges3.f"
	} else {
#line 519 "zgges3.f"
	    *info = *n + 1;
#line 520 "zgges3.f"
	}
#line 521 "zgges3.f"
	goto L30;
#line 522 "zgges3.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */

#line 526 "zgges3.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before selecting */

#line 530 "zgges3.f"
	if (ilascl) {
#line 530 "zgges3.f"
	    zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, &c__1, &alpha[1], n,
		     &ierr, (ftnlen)1);
#line 530 "zgges3.f"
	}
#line 532 "zgges3.f"
	if (ilbscl) {
#line 532 "zgges3.f"
	    zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 532 "zgges3.f"
	}

/*        Select eigenvalues */

#line 537 "zgges3.f"
	i__1 = *n;
#line 537 "zgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 538 "zgges3.f"
	    bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
#line 539 "zgges3.f"
/* L10: */
#line 539 "zgges3.f"
	}

#line 541 "zgges3.f"
	i__1 = *lwork - iwrk + 1;
#line 541 "zgges3.f"
	ztgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, sdim, &pvsl, &pvsr, dif, &work[iwrk],
		 &i__1, idum, &c__1, &ierr);
#line 544 "zgges3.f"
	if (ierr == 1) {
#line 544 "zgges3.f"
	    *info = *n + 3;
#line 544 "zgges3.f"
	}

#line 547 "zgges3.f"
    }

/*     Apply back-permutation to VSL and VSR */

#line 551 "zgges3.f"
    if (ilvsl) {
#line 551 "zgges3.f"
	zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 551 "zgges3.f"
    }
#line 554 "zgges3.f"
    if (ilvsr) {
#line 554 "zgges3.f"
	zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 554 "zgges3.f"
    }

/*     Undo scaling */

#line 560 "zgges3.f"
    if (ilascl) {
#line 561 "zgges3.f"
	zlascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 562 "zgges3.f"
	zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 563 "zgges3.f"
    }

#line 565 "zgges3.f"
    if (ilbscl) {
#line 566 "zgges3.f"
	zlascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 567 "zgges3.f"
	zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 568 "zgges3.f"
    }

#line 570 "zgges3.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 574 "zgges3.f"
	lastsl = TRUE_;
#line 575 "zgges3.f"
	*sdim = 0;
#line 576 "zgges3.f"
	i__1 = *n;
#line 576 "zgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 577 "zgges3.f"
	    cursl = (*selctg)(&alpha[i__], &beta[i__]);
#line 578 "zgges3.f"
	    if (cursl) {
#line 578 "zgges3.f"
		++(*sdim);
#line 578 "zgges3.f"
	    }
#line 580 "zgges3.f"
	    if (cursl && ! lastsl) {
#line 580 "zgges3.f"
		*info = *n + 2;
#line 580 "zgges3.f"
	    }
#line 582 "zgges3.f"
	    lastsl = cursl;
#line 583 "zgges3.f"
/* L20: */
#line 583 "zgges3.f"
	}

#line 585 "zgges3.f"
    }

#line 587 "zgges3.f"
L30:

#line 589 "zgges3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 589 "zgges3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

#line 591 "zgges3.f"
    return 0;

/*     End of ZGGES3 */

} /* zgges3_ */

