#line 1 "sgges3.f"
/* sgges3.f -- translated by f2c (version 20100827).
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

#line 1 "sgges3.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b36 = 0.;
static doublereal c_b37 = 1.;

/* > \brief <b> SGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices (blocked algorithm)</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGES3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgges3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgges3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgges3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, */
/*      $                   LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, */
/*      $                   VSR, LDVSR, WORK, LWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR, SORT */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ), */
/*      $                   VSR( LDVSR, * ), WORK( * ) */
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
/* > SGGES3 computes for a pair of N-by-N real nonsymmetric matrices (A,B), */
/* > the generalized eigenvalues, the generalized real Schur form (S,T), */
/* > optionally, the left and/or right matrices of Schur vectors (VSL and */
/* > VSR). This gives the generalized Schur factorization */
/* > */
/* >          (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T ) */
/* > */
/* > Optionally, it also orders the eigenvalues so that a selected cluster */
/* > of eigenvalues appears in the leading diagonal blocks of the upper */
/* > quasi-triangular matrix S and the upper triangular matrix T.The */
/* > leading columns of VSL and VSR then form an orthonormal basis for the */
/* > corresponding left and right eigenspaces (deflating subspaces). */
/* > */
/* > (If only the generalized eigenvalues are needed, use the driver */
/* > SGGEV instead, which is faster.) */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar w */
/* > or a ratio alpha/beta = w, such that  A - w*B is singular.  It is */
/* > usually represented as the pair (alpha,beta), as there is a */
/* > reasonable interpretation for beta=0 or both being zero. */
/* > */
/* > A pair of matrices (S,T) is in generalized real Schur form if T is */
/* > upper triangular with non-negative diagonal and S is block upper */
/* > triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond */
/* > to real generalized eigenvalues, while 2-by-2 blocks of S will be */
/* > "standardized" by making the corresponding elements of T have the */
/* > form: */
/* >         [  a  0  ] */
/* >         [  0  b  ] */
/* > */
/* > and the pair of corresponding 2-by-2 blocks in S and T will have a */
/* > complex conjugate pair of generalized eigenvalues. */
/* > */
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
/* >          = 'S':  Eigenvalues are ordered (see SELCTG); */
/* > \endverbatim */
/* > */
/* > \param[in] SELCTG */
/* > \verbatim */
/* >          SELCTG is a LOGICAL FUNCTION of three REAL arguments */
/* >          SELCTG must be declared EXTERNAL in the calling subroutine. */
/* >          If SORT = 'N', SELCTG is not referenced. */
/* >          If SORT = 'S', SELCTG is used to select eigenvalues to sort */
/* >          to the top left of the Schur form. */
/* >          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if */
/* >          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either */
/* >          one of a complex conjugate pair of eigenvalues is selected, */
/* >          then both complex eigenvalues are selected. */
/* > */
/* >          Note that in the ill-conditioned case, a selected complex */
/* >          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j), */
/* >          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2 */
/* >          in this case. */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          B is REAL array, dimension (LDB, N) */
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
/* >          for which SELCTG is true.  (Complex conjugate pairs for which */
/* >          SELCTG is true for either eigenvalue count as 2.) */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* >          ALPHAR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is REAL array, dimension (N) */
/* >          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will */
/* >          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i, */
/* >          and  BETA(j),j=1,...,N are the diagonals of the complex Schur */
/* >          form (S,T) that would result if the 2-by-2 diagonal blocks of */
/* >          the real Schur form of (A,B) were further reduced to */
/* >          triangular form using 2-by-2 complex unitary transformations. */
/* >          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if */
/* >          positive, then the j-th and (j+1)-st eigenvalues are a */
/* >          complex conjugate pair, with ALPHAI(j+1) negative. */
/* > */
/* >          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) */
/* >          may easily over- or underflow, and BETA(j) may even be zero. */
/* >          Thus, the user should avoid naively computing the ratio. */
/* >          However, ALPHAR and ALPHAI will be always less than and */
/* >          usually comparable with norm(A) in magnitude, and BETA always */
/* >          less than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* >          VSL is REAL array, dimension (LDVSL,N) */
/* >          If JOBVSL = 'V', VSL will contain the left Schur vectors. */
/* >          Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* >          LDVSL is INTEGER */
/* >          The leading dimension of the matrix VSL. LDVSL >=1, and */
/* >          if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* >          VSR is REAL array, dimension (LDVSR,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* >          = 1,...,N: */
/* >                The QZ iteration failed.  (A,B) are not in Schur */
/* >                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should */
/* >                be correct for j=INFO+1,...,N. */
/* >          > N:  =N+1: other than QZ iteration failed in SHGEQZ. */
/* >                =N+2: after reordering, roundoff changed values of */
/* >                      some complex eigenvalues so that leading */
/* >                      eigenvalues in the Generalized Schur form no */
/* >                      longer satisfy SELCTG=.TRUE.  This could also */
/* >                      be caused due to scaling. */
/* >                =N+3: reordering failed in STGSEN. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2015 */

/* > \ingroup realGEeigen */

/*  ===================================================================== */
/* Subroutine */ int sgges3_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai, 
	doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, 
	integer *ldvsr, doublereal *work, integer *lwork, logical *bwork, 
	integer *info, ftnlen jobvsl_len, ftnlen jobvsr_len, ftnlen sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ip;
    static doublereal dif[2];
    static integer ihi, ilo;
    static doublereal eps, anrm, bnrm;
    static integer idum[1], ierr, itau, iwrk;
    static doublereal pvsl, pvsr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ileft, icols;
    static logical cursl, ilvsl, ilvsr;
    static integer irows;
    extern /* Subroutine */ int sgghd3_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical lst2sl;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *), sggbak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), sggbal_(char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical ilascl, ilbscl;
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal safmax, bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer ijobvl, iright;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer ijobvr;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal anrmto, bnrmto;
    static logical lastsl;
    extern /* Subroutine */ int shgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), stgsen_(integer *, logical *, 
	    logical *, logical *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *);
    static doublereal smlnum;
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical wantst, lquery;
    static integer lwkopt;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 340 "sgges3.f"
    /* Parameter adjustments */
#line 340 "sgges3.f"
    a_dim1 = *lda;
#line 340 "sgges3.f"
    a_offset = 1 + a_dim1;
#line 340 "sgges3.f"
    a -= a_offset;
#line 340 "sgges3.f"
    b_dim1 = *ldb;
#line 340 "sgges3.f"
    b_offset = 1 + b_dim1;
#line 340 "sgges3.f"
    b -= b_offset;
#line 340 "sgges3.f"
    --alphar;
#line 340 "sgges3.f"
    --alphai;
#line 340 "sgges3.f"
    --beta;
#line 340 "sgges3.f"
    vsl_dim1 = *ldvsl;
#line 340 "sgges3.f"
    vsl_offset = 1 + vsl_dim1;
#line 340 "sgges3.f"
    vsl -= vsl_offset;
#line 340 "sgges3.f"
    vsr_dim1 = *ldvsr;
#line 340 "sgges3.f"
    vsr_offset = 1 + vsr_dim1;
#line 340 "sgges3.f"
    vsr -= vsr_offset;
#line 340 "sgges3.f"
    --work;
#line 340 "sgges3.f"
    --bwork;
#line 340 "sgges3.f"

#line 340 "sgges3.f"
    /* Function Body */
#line 340 "sgges3.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 341 "sgges3.f"
	ijobvl = 1;
#line 342 "sgges3.f"
	ilvsl = FALSE_;
#line 343 "sgges3.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 344 "sgges3.f"
	ijobvl = 2;
#line 345 "sgges3.f"
	ilvsl = TRUE_;
#line 346 "sgges3.f"
    } else {
#line 347 "sgges3.f"
	ijobvl = -1;
#line 348 "sgges3.f"
	ilvsl = FALSE_;
#line 349 "sgges3.f"
    }

#line 351 "sgges3.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 352 "sgges3.f"
	ijobvr = 1;
#line 353 "sgges3.f"
	ilvsr = FALSE_;
#line 354 "sgges3.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 355 "sgges3.f"
	ijobvr = 2;
#line 356 "sgges3.f"
	ilvsr = TRUE_;
#line 357 "sgges3.f"
    } else {
#line 358 "sgges3.f"
	ijobvr = -1;
#line 359 "sgges3.f"
	ilvsr = FALSE_;
#line 360 "sgges3.f"
    }

#line 362 "sgges3.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 366 "sgges3.f"
    *info = 0;
#line 367 "sgges3.f"
    lquery = *lwork == -1;
#line 368 "sgges3.f"
    if (ijobvl <= 0) {
#line 369 "sgges3.f"
	*info = -1;
#line 370 "sgges3.f"
    } else if (ijobvr <= 0) {
#line 371 "sgges3.f"
	*info = -2;
#line 372 "sgges3.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 373 "sgges3.f"
	*info = -3;
#line 374 "sgges3.f"
    } else if (*n < 0) {
#line 375 "sgges3.f"
	*info = -5;
#line 376 "sgges3.f"
    } else if (*lda < max(1,*n)) {
#line 377 "sgges3.f"
	*info = -7;
#line 378 "sgges3.f"
    } else if (*ldb < max(1,*n)) {
#line 379 "sgges3.f"
	*info = -9;
#line 380 "sgges3.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 381 "sgges3.f"
	*info = -15;
#line 382 "sgges3.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 383 "sgges3.f"
	*info = -17;
#line 384 "sgges3.f"
    } else if (*lwork < *n * 6 + 16 && ! lquery) {
#line 385 "sgges3.f"
	*info = -19;
#line 386 "sgges3.f"
    }

/*     Compute workspace */

#line 390 "sgges3.f"
    if (*info == 0) {
#line 391 "sgges3.f"
	sgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
/* Computing MAX */
#line 392 "sgges3.f"
	i__1 = *n * 6 + 16, i__2 = *n * 3 + (integer) work[1];
#line 392 "sgges3.f"
	lwkopt = max(i__1,i__2);
#line 393 "sgges3.f"
	sormqr_("L", "T", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], 
		lda, &work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 395 "sgges3.f"
	i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 395 "sgges3.f"
	lwkopt = max(i__1,i__2);
#line 396 "sgges3.f"
	if (ilvsl) {
#line 397 "sgges3.f"
	    sorgqr_(n, n, n, &vsl[vsl_offset], ldvsl, &work[1], &work[1], &
		    c_n1, &ierr);
/* Computing MAX */
#line 398 "sgges3.f"
	    i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 398 "sgges3.f"
	    lwkopt = max(i__1,i__2);
#line 399 "sgges3.f"
	}
#line 400 "sgges3.f"
	sgghd3_(jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[
		1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 402 "sgges3.f"
	i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 402 "sgges3.f"
	lwkopt = max(i__1,i__2);
#line 403 "sgges3.f"
	shgeqz_("S", jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[1], &c_n1, 
		&ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 406 "sgges3.f"
	i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 406 "sgges3.f"
	lwkopt = max(i__1,i__2);
#line 407 "sgges3.f"
	if (wantst) {
#line 408 "sgges3.f"
	    stgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &
		    b[b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		    vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, 
		    &pvsr, dif, &work[1], &c_n1, idum, &c__1, &ierr);
/* Computing MAX */
#line 412 "sgges3.f"
	    i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 412 "sgges3.f"
	    lwkopt = max(i__1,i__2);
#line 413 "sgges3.f"
	}
#line 414 "sgges3.f"
	work[1] = (doublereal) lwkopt;
#line 415 "sgges3.f"
    }

#line 417 "sgges3.f"
    if (*info != 0) {
#line 418 "sgges3.f"
	i__1 = -(*info);
#line 418 "sgges3.f"
	xerbla_("SGGES3 ", &i__1, (ftnlen)7);
#line 419 "sgges3.f"
	return 0;
#line 420 "sgges3.f"
    } else if (lquery) {
#line 421 "sgges3.f"
	return 0;
#line 422 "sgges3.f"
    }

/*     Quick return if possible */

#line 426 "sgges3.f"
    if (*n == 0) {
#line 427 "sgges3.f"
	*sdim = 0;
#line 428 "sgges3.f"
	return 0;
#line 429 "sgges3.f"
    }

/*     Get machine constants */

#line 433 "sgges3.f"
    eps = slamch_("P", (ftnlen)1);
#line 434 "sgges3.f"
    safmin = slamch_("S", (ftnlen)1);
#line 435 "sgges3.f"
    safmax = 1. / safmin;
#line 436 "sgges3.f"
    slabad_(&safmin, &safmax);
#line 437 "sgges3.f"
    smlnum = sqrt(safmin) / eps;
#line 438 "sgges3.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 442 "sgges3.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 443 "sgges3.f"
    ilascl = FALSE_;
#line 444 "sgges3.f"
    if (anrm > 0. && anrm < smlnum) {
#line 445 "sgges3.f"
	anrmto = smlnum;
#line 446 "sgges3.f"
	ilascl = TRUE_;
#line 447 "sgges3.f"
    } else if (anrm > bignum) {
#line 448 "sgges3.f"
	anrmto = bignum;
#line 449 "sgges3.f"
	ilascl = TRUE_;
#line 450 "sgges3.f"
    }
#line 451 "sgges3.f"
    if (ilascl) {
#line 451 "sgges3.f"
	slascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 451 "sgges3.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 456 "sgges3.f"
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 457 "sgges3.f"
    ilbscl = FALSE_;
#line 458 "sgges3.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 459 "sgges3.f"
	bnrmto = smlnum;
#line 460 "sgges3.f"
	ilbscl = TRUE_;
#line 461 "sgges3.f"
    } else if (bnrm > bignum) {
#line 462 "sgges3.f"
	bnrmto = bignum;
#line 463 "sgges3.f"
	ilbscl = TRUE_;
#line 464 "sgges3.f"
    }
#line 465 "sgges3.f"
    if (ilbscl) {
#line 465 "sgges3.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 465 "sgges3.f"
    }

/*     Permute the matrix to make it more nearly triangular */

#line 470 "sgges3.f"
    ileft = 1;
#line 471 "sgges3.f"
    iright = *n + 1;
#line 472 "sgges3.f"
    iwrk = iright + *n;
#line 473 "sgges3.f"
    sggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */

#line 478 "sgges3.f"
    irows = ihi + 1 - ilo;
#line 479 "sgges3.f"
    icols = *n + 1 - ilo;
#line 480 "sgges3.f"
    itau = iwrk;
#line 481 "sgges3.f"
    iwrk = itau + irows;
#line 482 "sgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 482 "sgges3.f"
    sgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */

#line 487 "sgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 487 "sgges3.f"
    sormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */

#line 493 "sgges3.f"
    if (ilvsl) {
#line 494 "sgges3.f"
	slaset_("Full", n, n, &c_b36, &c_b37, &vsl[vsl_offset], ldvsl, (
		ftnlen)4);
#line 495 "sgges3.f"
	if (irows > 1) {
#line 496 "sgges3.f"
	    i__1 = irows - 1;
#line 496 "sgges3.f"
	    i__2 = irows - 1;
#line 496 "sgges3.f"
	    slacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 498 "sgges3.f"
	}
#line 499 "sgges3.f"
	i__1 = *lwork + 1 - iwrk;
#line 499 "sgges3.f"
	sorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 501 "sgges3.f"
    }

/*     Initialize VSR */

#line 505 "sgges3.f"
    if (ilvsr) {
#line 505 "sgges3.f"
	slaset_("Full", n, n, &c_b36, &c_b37, &vsr[vsr_offset], ldvsr, (
		ftnlen)4);
#line 505 "sgges3.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 510 "sgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 510 "sgges3.f"
    sgghd3_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk]
	    , &i__1, &ierr, (ftnlen)1, (ftnlen)1);

/*     Perform QZ algorithm, computing Schur vectors if desired */

#line 515 "sgges3.f"
    iwrk = itau;
#line 516 "sgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 516 "sgges3.f"
    shgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
	    , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 519 "sgges3.f"
    if (ierr != 0) {
#line 520 "sgges3.f"
	if (ierr > 0 && ierr <= *n) {
#line 521 "sgges3.f"
	    *info = ierr;
#line 522 "sgges3.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 523 "sgges3.f"
	    *info = ierr - *n;
#line 524 "sgges3.f"
	} else {
#line 525 "sgges3.f"
	    *info = *n + 1;
#line 526 "sgges3.f"
	}
#line 527 "sgges3.f"
	goto L40;
#line 528 "sgges3.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */

#line 532 "sgges3.f"
    *sdim = 0;
#line 533 "sgges3.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 537 "sgges3.f"
	if (ilascl) {
#line 538 "sgges3.f"
	    slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], 
		    n, &ierr, (ftnlen)1);
#line 540 "sgges3.f"
	    slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], 
		    n, &ierr, (ftnlen)1);
#line 542 "sgges3.f"
	}
#line 543 "sgges3.f"
	if (ilbscl) {
#line 543 "sgges3.f"
	    slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 543 "sgges3.f"
	}

/*        Select eigenvalues */

#line 548 "sgges3.f"
	i__1 = *n;
#line 548 "sgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 549 "sgges3.f"
	    bwork[i__] = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 550 "sgges3.f"
/* L10: */
#line 550 "sgges3.f"
	}

#line 552 "sgges3.f"
	i__1 = *lwork - iwrk + 1;
#line 552 "sgges3.f"
	stgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, &
		pvsr, dif, &work[iwrk], &i__1, idum, &c__1, &ierr);
#line 556 "sgges3.f"
	if (ierr == 1) {
#line 556 "sgges3.f"
	    *info = *n + 3;
#line 556 "sgges3.f"
	}

#line 559 "sgges3.f"
    }

/*     Apply back-permutation to VSL and VSR */

#line 563 "sgges3.f"
    if (ilvsl) {
#line 563 "sgges3.f"
	sggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 563 "sgges3.f"
    }

#line 567 "sgges3.f"
    if (ilvsr) {
#line 567 "sgges3.f"
	sggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 567 "sgges3.f"
    }

/*     Check if unscaling would cause over/underflow, if so, rescale */
/*     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of */
/*     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I) */

#line 575 "sgges3.f"
    if (ilascl) {
#line 576 "sgges3.f"
	i__1 = *n;
#line 576 "sgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 577 "sgges3.f"
	    if (alphai[i__] != 0.) {
#line 578 "sgges3.f"
		if (alphar[i__] / safmax > anrmto / anrm || safmin / alphar[
			i__] > anrm / anrmto) {
#line 580 "sgges3.f"
		    work[1] = (d__1 = a[i__ + i__ * a_dim1] / alphar[i__], 
			    abs(d__1));
#line 581 "sgges3.f"
		    beta[i__] *= work[1];
#line 582 "sgges3.f"
		    alphar[i__] *= work[1];
#line 583 "sgges3.f"
		    alphai[i__] *= work[1];
#line 584 "sgges3.f"
		} else if (alphai[i__] / safmax > anrmto / anrm || safmin / 
			alphai[i__] > anrm / anrmto) {
#line 586 "sgges3.f"
		    work[1] = (d__1 = a[i__ + (i__ + 1) * a_dim1] / alphai[
			    i__], abs(d__1));
#line 587 "sgges3.f"
		    beta[i__] *= work[1];
#line 588 "sgges3.f"
		    alphar[i__] *= work[1];
#line 589 "sgges3.f"
		    alphai[i__] *= work[1];
#line 590 "sgges3.f"
		}
#line 591 "sgges3.f"
	    }
#line 592 "sgges3.f"
/* L50: */
#line 592 "sgges3.f"
	}
#line 593 "sgges3.f"
    }

#line 595 "sgges3.f"
    if (ilbscl) {
#line 596 "sgges3.f"
	i__1 = *n;
#line 596 "sgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 597 "sgges3.f"
	    if (alphai[i__] != 0.) {
#line 598 "sgges3.f"
		if (beta[i__] / safmax > bnrmto / bnrm || safmin / beta[i__] 
			> bnrm / bnrmto) {
#line 600 "sgges3.f"
		    work[1] = (d__1 = b[i__ + i__ * b_dim1] / beta[i__], abs(
			    d__1));
#line 601 "sgges3.f"
		    beta[i__] *= work[1];
#line 602 "sgges3.f"
		    alphar[i__] *= work[1];
#line 603 "sgges3.f"
		    alphai[i__] *= work[1];
#line 604 "sgges3.f"
		}
#line 605 "sgges3.f"
	    }
#line 606 "sgges3.f"
/* L60: */
#line 606 "sgges3.f"
	}
#line 607 "sgges3.f"
    }

/*     Undo scaling */

#line 611 "sgges3.f"
    if (ilascl) {
#line 612 "sgges3.f"
	slascl_("H", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 613 "sgges3.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 614 "sgges3.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 615 "sgges3.f"
    }

#line 617 "sgges3.f"
    if (ilbscl) {
#line 618 "sgges3.f"
	slascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 619 "sgges3.f"
	slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 620 "sgges3.f"
    }

#line 622 "sgges3.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 626 "sgges3.f"
	lastsl = TRUE_;
#line 627 "sgges3.f"
	lst2sl = TRUE_;
#line 628 "sgges3.f"
	*sdim = 0;
#line 629 "sgges3.f"
	ip = 0;
#line 630 "sgges3.f"
	i__1 = *n;
#line 630 "sgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 631 "sgges3.f"
	    cursl = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 632 "sgges3.f"
	    if (alphai[i__] == 0.) {
#line 633 "sgges3.f"
		if (cursl) {
#line 633 "sgges3.f"
		    ++(*sdim);
#line 633 "sgges3.f"
		}
#line 635 "sgges3.f"
		ip = 0;
#line 636 "sgges3.f"
		if (cursl && ! lastsl) {
#line 636 "sgges3.f"
		    *info = *n + 2;
#line 636 "sgges3.f"
		}
#line 638 "sgges3.f"
	    } else {
#line 639 "sgges3.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 643 "sgges3.f"
		    cursl = cursl || lastsl;
#line 644 "sgges3.f"
		    lastsl = cursl;
#line 645 "sgges3.f"
		    if (cursl) {
#line 645 "sgges3.f"
			*sdim += 2;
#line 645 "sgges3.f"
		    }
#line 647 "sgges3.f"
		    ip = -1;
#line 648 "sgges3.f"
		    if (cursl && ! lst2sl) {
#line 648 "sgges3.f"
			*info = *n + 2;
#line 648 "sgges3.f"
		    }
#line 650 "sgges3.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 654 "sgges3.f"
		    ip = 1;
#line 655 "sgges3.f"
		}
#line 656 "sgges3.f"
	    }
#line 657 "sgges3.f"
	    lst2sl = lastsl;
#line 658 "sgges3.f"
	    lastsl = cursl;
#line 659 "sgges3.f"
/* L30: */
#line 659 "sgges3.f"
	}

#line 661 "sgges3.f"
    }

#line 663 "sgges3.f"
L40:

#line 665 "sgges3.f"
    work[1] = (doublereal) lwkopt;

#line 667 "sgges3.f"
    return 0;

/*     End of SGGES3 */

} /* sgges3_ */

