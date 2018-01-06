#line 1 "cgges3.f"
/* cgges3.f -- translated by f2c (version 20100827).
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

#line 1 "cgges3.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;

/* > \brief <b> CGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices (blocked algorithm)</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGES3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgges3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgges3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgges3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, */
/*      $                   LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, */
/*      $                   WORK, LWORK, RWORK, BWORK, INFO ) */

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
/* > CGGES3 computes for a pair of N-by-N complex nonsymmetric matrices */
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
/* >          output by CGGES3. The  BETA(j) will be non-negative real. */
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

/* > \date January 2015 */

/* > \ingroup complexGEeigen */

/*  ===================================================================== */
/* Subroutine */ int cgges3_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
    static integer irwrk;
    extern /* Subroutine */ int cgghd3_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);
    static integer irows;
    extern /* Subroutine */ int cggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), cggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), slabad_(doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
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
    static doublereal bignum;
    extern /* Subroutine */ int chgeqz_(), ctgsen_(integer *, logical *, 
	    logical *, logical *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, doublecomplex *, 
	    integer *, integer *, integer *, integer *);
    static integer ijobvl, iright, ijobvr;
    static doublereal anrmto, bnrmto;
    static logical lastsl;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), cunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static doublereal smlnum;
    static logical wantst, lquery;
    static integer lwkopt;


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

#line 331 "cgges3.f"
    /* Parameter adjustments */
#line 331 "cgges3.f"
    a_dim1 = *lda;
#line 331 "cgges3.f"
    a_offset = 1 + a_dim1;
#line 331 "cgges3.f"
    a -= a_offset;
#line 331 "cgges3.f"
    b_dim1 = *ldb;
#line 331 "cgges3.f"
    b_offset = 1 + b_dim1;
#line 331 "cgges3.f"
    b -= b_offset;
#line 331 "cgges3.f"
    --alpha;
#line 331 "cgges3.f"
    --beta;
#line 331 "cgges3.f"
    vsl_dim1 = *ldvsl;
#line 331 "cgges3.f"
    vsl_offset = 1 + vsl_dim1;
#line 331 "cgges3.f"
    vsl -= vsl_offset;
#line 331 "cgges3.f"
    vsr_dim1 = *ldvsr;
#line 331 "cgges3.f"
    vsr_offset = 1 + vsr_dim1;
#line 331 "cgges3.f"
    vsr -= vsr_offset;
#line 331 "cgges3.f"
    --work;
#line 331 "cgges3.f"
    --rwork;
#line 331 "cgges3.f"
    --bwork;
#line 331 "cgges3.f"

#line 331 "cgges3.f"
    /* Function Body */
#line 331 "cgges3.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 332 "cgges3.f"
	ijobvl = 1;
#line 333 "cgges3.f"
	ilvsl = FALSE_;
#line 334 "cgges3.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 335 "cgges3.f"
	ijobvl = 2;
#line 336 "cgges3.f"
	ilvsl = TRUE_;
#line 337 "cgges3.f"
    } else {
#line 338 "cgges3.f"
	ijobvl = -1;
#line 339 "cgges3.f"
	ilvsl = FALSE_;
#line 340 "cgges3.f"
    }

#line 342 "cgges3.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 343 "cgges3.f"
	ijobvr = 1;
#line 344 "cgges3.f"
	ilvsr = FALSE_;
#line 345 "cgges3.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 346 "cgges3.f"
	ijobvr = 2;
#line 347 "cgges3.f"
	ilvsr = TRUE_;
#line 348 "cgges3.f"
    } else {
#line 349 "cgges3.f"
	ijobvr = -1;
#line 350 "cgges3.f"
	ilvsr = FALSE_;
#line 351 "cgges3.f"
    }

#line 353 "cgges3.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 357 "cgges3.f"
    *info = 0;
#line 358 "cgges3.f"
    lquery = *lwork == -1;
#line 359 "cgges3.f"
    if (ijobvl <= 0) {
#line 360 "cgges3.f"
	*info = -1;
#line 361 "cgges3.f"
    } else if (ijobvr <= 0) {
#line 362 "cgges3.f"
	*info = -2;
#line 363 "cgges3.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 364 "cgges3.f"
	*info = -3;
#line 365 "cgges3.f"
    } else if (*n < 0) {
#line 366 "cgges3.f"
	*info = -5;
#line 367 "cgges3.f"
    } else if (*lda < max(1,*n)) {
#line 368 "cgges3.f"
	*info = -7;
#line 369 "cgges3.f"
    } else if (*ldb < max(1,*n)) {
#line 370 "cgges3.f"
	*info = -9;
#line 371 "cgges3.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 372 "cgges3.f"
	*info = -14;
#line 373 "cgges3.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 374 "cgges3.f"
	*info = -16;
#line 375 "cgges3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 375 "cgges3.f"
	i__1 = 1, i__2 = *n << 1;
#line 375 "cgges3.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 376 "cgges3.f"
	    *info = -18;
#line 377 "cgges3.f"
	}
#line 377 "cgges3.f"
    }

/*     Compute workspace */

#line 381 "cgges3.f"
    if (*info == 0) {
#line 382 "cgges3.f"
	cgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
/* Computing MAX */
#line 383 "cgges3.f"
	i__1 = 1, i__2 = *n + (integer) work[1].r;
#line 383 "cgges3.f"
	lwkopt = max(i__1,i__2);
#line 384 "cgges3.f"
	cunmqr_("L", "C", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], 
		lda, &work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 386 "cgges3.f"
	i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 386 "cgges3.f"
	lwkopt = max(i__1,i__2);
#line 387 "cgges3.f"
	if (ilvsl) {
#line 388 "cgges3.f"
	    cungqr_(n, n, n, &vsl[vsl_offset], ldvsl, &work[1], &work[1], &
		    c_n1, &ierr);
/* Computing MAX */
#line 390 "cgges3.f"
	    i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 390 "cgges3.f"
	    lwkopt = max(i__1,i__2);
#line 391 "cgges3.f"
	}
#line 392 "cgges3.f"
	cgghd3_(jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[
		1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 394 "cgges3.f"
	i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 394 "cgges3.f"
	lwkopt = max(i__1,i__2);
#line 395 "cgges3.f"
	chgeqz_("S", jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, &work[1], &c_n1, &work[1], &ierr, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 398 "cgges3.f"
	i__1 = lwkopt, i__2 = (integer) work[1].r;
#line 398 "cgges3.f"
	lwkopt = max(i__1,i__2);
#line 399 "cgges3.f"
	if (wantst) {
#line 400 "cgges3.f"
	    ctgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &
		    b[b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], 
		    ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, &pvsr, dif, &
		    work[1], &c_n1, idum, &c__1, &ierr);
/* Computing MAX */
#line 403 "cgges3.f"
	    i__1 = lwkopt, i__2 = (integer) work[1].r;
#line 403 "cgges3.f"
	    lwkopt = max(i__1,i__2);
#line 404 "cgges3.f"
	}
#line 405 "cgges3.f"
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 405 "cgges3.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 406 "cgges3.f"
    }

#line 409 "cgges3.f"
    if (*info != 0) {
#line 410 "cgges3.f"
	i__1 = -(*info);
#line 410 "cgges3.f"
	xerbla_("CGGES3 ", &i__1, (ftnlen)7);
#line 411 "cgges3.f"
	return 0;
#line 412 "cgges3.f"
    } else if (lquery) {
#line 413 "cgges3.f"
	return 0;
#line 414 "cgges3.f"
    }

/*     Quick return if possible */

#line 418 "cgges3.f"
    if (*n == 0) {
#line 419 "cgges3.f"
	*sdim = 0;
#line 420 "cgges3.f"
	return 0;
#line 421 "cgges3.f"
    }

/*     Get machine constants */

#line 425 "cgges3.f"
    eps = slamch_("P", (ftnlen)1);
#line 426 "cgges3.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 427 "cgges3.f"
    bignum = 1. / smlnum;
#line 428 "cgges3.f"
    slabad_(&smlnum, &bignum);
#line 429 "cgges3.f"
    smlnum = sqrt(smlnum) / eps;
#line 430 "cgges3.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 434 "cgges3.f"
    anrm = clange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 435 "cgges3.f"
    ilascl = FALSE_;
#line 436 "cgges3.f"
    if (anrm > 0. && anrm < smlnum) {
#line 437 "cgges3.f"
	anrmto = smlnum;
#line 438 "cgges3.f"
	ilascl = TRUE_;
#line 439 "cgges3.f"
    } else if (anrm > bignum) {
#line 440 "cgges3.f"
	anrmto = bignum;
#line 441 "cgges3.f"
	ilascl = TRUE_;
#line 442 "cgges3.f"
    }

#line 444 "cgges3.f"
    if (ilascl) {
#line 444 "cgges3.f"
	clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 444 "cgges3.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 449 "cgges3.f"
    bnrm = clange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 450 "cgges3.f"
    ilbscl = FALSE_;
#line 451 "cgges3.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 452 "cgges3.f"
	bnrmto = smlnum;
#line 453 "cgges3.f"
	ilbscl = TRUE_;
#line 454 "cgges3.f"
    } else if (bnrm > bignum) {
#line 455 "cgges3.f"
	bnrmto = bignum;
#line 456 "cgges3.f"
	ilbscl = TRUE_;
#line 457 "cgges3.f"
    }

#line 459 "cgges3.f"
    if (ilbscl) {
#line 459 "cgges3.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 459 "cgges3.f"
    }

/*     Permute the matrix to make it more nearly triangular */

#line 464 "cgges3.f"
    ileft = 1;
#line 465 "cgges3.f"
    iright = *n + 1;
#line 466 "cgges3.f"
    irwrk = iright + *n;
#line 467 "cgges3.f"
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */

#line 472 "cgges3.f"
    irows = ihi + 1 - ilo;
#line 473 "cgges3.f"
    icols = *n + 1 - ilo;
#line 474 "cgges3.f"
    itau = 1;
#line 475 "cgges3.f"
    iwrk = itau + irows;
#line 476 "cgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 476 "cgges3.f"
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */

#line 481 "cgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 481 "cgges3.f"
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */

#line 487 "cgges3.f"
    if (ilvsl) {
#line 488 "cgges3.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 489 "cgges3.f"
	if (irows > 1) {
#line 490 "cgges3.f"
	    i__1 = irows - 1;
#line 490 "cgges3.f"
	    i__2 = irows - 1;
#line 490 "cgges3.f"
	    clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 492 "cgges3.f"
	}
#line 493 "cgges3.f"
	i__1 = *lwork + 1 - iwrk;
#line 493 "cgges3.f"
	cungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 495 "cgges3.f"
    }

/*     Initialize VSR */

#line 499 "cgges3.f"
    if (ilvsr) {
#line 499 "cgges3.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 499 "cgges3.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 504 "cgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 504 "cgges3.f"
    cgghd3_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk]
	    , &i__1, &ierr, (ftnlen)1, (ftnlen)1);

#line 507 "cgges3.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */

#line 511 "cgges3.f"
    iwrk = itau;
#line 512 "cgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 512 "cgges3.f"
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 515 "cgges3.f"
    if (ierr != 0) {
#line 516 "cgges3.f"
	if (ierr > 0 && ierr <= *n) {
#line 517 "cgges3.f"
	    *info = ierr;
#line 518 "cgges3.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 519 "cgges3.f"
	    *info = ierr - *n;
#line 520 "cgges3.f"
	} else {
#line 521 "cgges3.f"
	    *info = *n + 1;
#line 522 "cgges3.f"
	}
#line 523 "cgges3.f"
	goto L30;
#line 524 "cgges3.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */

#line 528 "cgges3.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before selecting */

#line 532 "cgges3.f"
	if (ilascl) {
#line 532 "cgges3.f"
	    clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, &c__1, &alpha[1], n,
		     &ierr, (ftnlen)1);
#line 532 "cgges3.f"
	}
#line 534 "cgges3.f"
	if (ilbscl) {
#line 534 "cgges3.f"
	    clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 534 "cgges3.f"
	}

/*        Select eigenvalues */

#line 539 "cgges3.f"
	i__1 = *n;
#line 539 "cgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 540 "cgges3.f"
	    bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
#line 541 "cgges3.f"
/* L10: */
#line 541 "cgges3.f"
	}

#line 543 "cgges3.f"
	i__1 = *lwork - iwrk + 1;
#line 543 "cgges3.f"
	ctgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, sdim, &pvsl, &pvsr, dif, &work[iwrk],
		 &i__1, idum, &c__1, &ierr);
#line 546 "cgges3.f"
	if (ierr == 1) {
#line 546 "cgges3.f"
	    *info = *n + 3;
#line 546 "cgges3.f"
	}

#line 549 "cgges3.f"
    }

/*     Apply back-permutation to VSL and VSR */

#line 553 "cgges3.f"
    if (ilvsl) {
#line 553 "cgges3.f"
	cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 553 "cgges3.f"
    }
#line 556 "cgges3.f"
    if (ilvsr) {
#line 556 "cgges3.f"
	cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 556 "cgges3.f"
    }

/*     Undo scaling */

#line 562 "cgges3.f"
    if (ilascl) {
#line 563 "cgges3.f"
	clascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 564 "cgges3.f"
	clascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 565 "cgges3.f"
    }

#line 567 "cgges3.f"
    if (ilbscl) {
#line 568 "cgges3.f"
	clascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 569 "cgges3.f"
	clascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 570 "cgges3.f"
    }

#line 572 "cgges3.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 576 "cgges3.f"
	lastsl = TRUE_;
#line 577 "cgges3.f"
	*sdim = 0;
#line 578 "cgges3.f"
	i__1 = *n;
#line 578 "cgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 579 "cgges3.f"
	    cursl = (*selctg)(&alpha[i__], &beta[i__]);
#line 580 "cgges3.f"
	    if (cursl) {
#line 580 "cgges3.f"
		++(*sdim);
#line 580 "cgges3.f"
	    }
#line 582 "cgges3.f"
	    if (cursl && ! lastsl) {
#line 582 "cgges3.f"
		*info = *n + 2;
#line 582 "cgges3.f"
	    }
#line 584 "cgges3.f"
	    lastsl = cursl;
#line 585 "cgges3.f"
/* L20: */
#line 585 "cgges3.f"
	}

#line 587 "cgges3.f"
    }

#line 589 "cgges3.f"
L30:

#line 591 "cgges3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 591 "cgges3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

#line 593 "cgges3.f"
    return 0;

/*     End of CGGES3 */

} /* cgges3_ */

