#line 1 "dgges3.f"
/* dgges3.f -- translated by f2c (version 20100827).
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

#line 1 "dgges3.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b36 = 0.;
static doublereal c_b37 = 1.;

/* > \brief <b> DGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices (blocked algorithm)</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGES3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgges3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgges3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgges3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, */
/*                          SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, */
/*                          LDVSR, WORK, LWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR, SORT */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
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
/* > DGGES3 computes for a pair of N-by-N real nonsymmetric matrices (A,B), */
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
/* > DGGEV instead, which is faster.) */
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
/* >          SELCTG is a LOGICAL FUNCTION of three DOUBLE PRECISION arguments */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
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
/* >          ALPHAR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION array, dimension (N) */
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
/* >          VSL is DOUBLE PRECISION array, dimension (LDVSL,N) */
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
/* >          VSR is DOUBLE PRECISION array, dimension (LDVSR,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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
/* >          > N:  =N+1: other than QZ iteration failed in DHGEQZ. */
/* >                =N+2: after reordering, roundoff changed values of */
/* >                      some complex eigenvalues so that leading */
/* >                      eigenvalues in the Generalized Schur form no */
/* >                      longer satisfy SELCTG=.TRUE.  This could also */
/* >                      be caused due to scaling. */
/* >                =N+3: reordering failed in DTGSEN. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2015 */

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dgges3_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
    extern /* Subroutine */ int dgghd3_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer irows;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dggbak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dggbal_(char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical lst2sl;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dtgsen_(integer *, logical *, 
	    logical *, logical *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *);
    static integer ijobvl, iright, ijobvr;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal anrmto, bnrmto;
    static logical lastsl;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
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

#line 340 "dgges3.f"
    /* Parameter adjustments */
#line 340 "dgges3.f"
    a_dim1 = *lda;
#line 340 "dgges3.f"
    a_offset = 1 + a_dim1;
#line 340 "dgges3.f"
    a -= a_offset;
#line 340 "dgges3.f"
    b_dim1 = *ldb;
#line 340 "dgges3.f"
    b_offset = 1 + b_dim1;
#line 340 "dgges3.f"
    b -= b_offset;
#line 340 "dgges3.f"
    --alphar;
#line 340 "dgges3.f"
    --alphai;
#line 340 "dgges3.f"
    --beta;
#line 340 "dgges3.f"
    vsl_dim1 = *ldvsl;
#line 340 "dgges3.f"
    vsl_offset = 1 + vsl_dim1;
#line 340 "dgges3.f"
    vsl -= vsl_offset;
#line 340 "dgges3.f"
    vsr_dim1 = *ldvsr;
#line 340 "dgges3.f"
    vsr_offset = 1 + vsr_dim1;
#line 340 "dgges3.f"
    vsr -= vsr_offset;
#line 340 "dgges3.f"
    --work;
#line 340 "dgges3.f"
    --bwork;
#line 340 "dgges3.f"

#line 340 "dgges3.f"
    /* Function Body */
#line 340 "dgges3.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 341 "dgges3.f"
	ijobvl = 1;
#line 342 "dgges3.f"
	ilvsl = FALSE_;
#line 343 "dgges3.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 344 "dgges3.f"
	ijobvl = 2;
#line 345 "dgges3.f"
	ilvsl = TRUE_;
#line 346 "dgges3.f"
    } else {
#line 347 "dgges3.f"
	ijobvl = -1;
#line 348 "dgges3.f"
	ilvsl = FALSE_;
#line 349 "dgges3.f"
    }

#line 351 "dgges3.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 352 "dgges3.f"
	ijobvr = 1;
#line 353 "dgges3.f"
	ilvsr = FALSE_;
#line 354 "dgges3.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 355 "dgges3.f"
	ijobvr = 2;
#line 356 "dgges3.f"
	ilvsr = TRUE_;
#line 357 "dgges3.f"
    } else {
#line 358 "dgges3.f"
	ijobvr = -1;
#line 359 "dgges3.f"
	ilvsr = FALSE_;
#line 360 "dgges3.f"
    }

#line 362 "dgges3.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 366 "dgges3.f"
    *info = 0;
#line 367 "dgges3.f"
    lquery = *lwork == -1;
#line 368 "dgges3.f"
    if (ijobvl <= 0) {
#line 369 "dgges3.f"
	*info = -1;
#line 370 "dgges3.f"
    } else if (ijobvr <= 0) {
#line 371 "dgges3.f"
	*info = -2;
#line 372 "dgges3.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 373 "dgges3.f"
	*info = -3;
#line 374 "dgges3.f"
    } else if (*n < 0) {
#line 375 "dgges3.f"
	*info = -5;
#line 376 "dgges3.f"
    } else if (*lda < max(1,*n)) {
#line 377 "dgges3.f"
	*info = -7;
#line 378 "dgges3.f"
    } else if (*ldb < max(1,*n)) {
#line 379 "dgges3.f"
	*info = -9;
#line 380 "dgges3.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 381 "dgges3.f"
	*info = -15;
#line 382 "dgges3.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 383 "dgges3.f"
	*info = -17;
#line 384 "dgges3.f"
    } else if (*lwork < *n * 6 + 16 && ! lquery) {
#line 385 "dgges3.f"
	*info = -19;
#line 386 "dgges3.f"
    }

/*     Compute workspace */

#line 390 "dgges3.f"
    if (*info == 0) {
#line 391 "dgges3.f"
	dgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
/* Computing MAX */
#line 392 "dgges3.f"
	i__1 = *n * 6 + 16, i__2 = *n * 3 + (integer) work[1];
#line 392 "dgges3.f"
	lwkopt = max(i__1,i__2);
#line 393 "dgges3.f"
	dormqr_("L", "T", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], 
		lda, &work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 395 "dgges3.f"
	i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 395 "dgges3.f"
	lwkopt = max(i__1,i__2);
#line 396 "dgges3.f"
	if (ilvsl) {
#line 397 "dgges3.f"
	    dorgqr_(n, n, n, &vsl[vsl_offset], ldvsl, &work[1], &work[1], &
		    c_n1, &ierr);
/* Computing MAX */
#line 398 "dgges3.f"
	    i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 398 "dgges3.f"
	    lwkopt = max(i__1,i__2);
#line 399 "dgges3.f"
	}
#line 400 "dgges3.f"
	dgghd3_(jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[
		1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 402 "dgges3.f"
	i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 402 "dgges3.f"
	lwkopt = max(i__1,i__2);
#line 403 "dgges3.f"
	dhgeqz_("S", jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[1], &c_n1, 
		&ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 406 "dgges3.f"
	i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 406 "dgges3.f"
	lwkopt = max(i__1,i__2);
#line 407 "dgges3.f"
	if (wantst) {
#line 408 "dgges3.f"
	    dtgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &
		    b[b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		    vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, 
		    &pvsr, dif, &work[1], &c_n1, idum, &c__1, &ierr);
/* Computing MAX */
#line 412 "dgges3.f"
	    i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 412 "dgges3.f"
	    lwkopt = max(i__1,i__2);
#line 413 "dgges3.f"
	}
#line 414 "dgges3.f"
	work[1] = (doublereal) lwkopt;
#line 415 "dgges3.f"
    }

#line 417 "dgges3.f"
    if (*info != 0) {
#line 418 "dgges3.f"
	i__1 = -(*info);
#line 418 "dgges3.f"
	xerbla_("DGGES3 ", &i__1, (ftnlen)7);
#line 419 "dgges3.f"
	return 0;
#line 420 "dgges3.f"
    } else if (lquery) {
#line 421 "dgges3.f"
	return 0;
#line 422 "dgges3.f"
    }

/*     Quick return if possible */

#line 426 "dgges3.f"
    if (*n == 0) {
#line 427 "dgges3.f"
	*sdim = 0;
#line 428 "dgges3.f"
	return 0;
#line 429 "dgges3.f"
    }

/*     Get machine constants */

#line 433 "dgges3.f"
    eps = dlamch_("P", (ftnlen)1);
#line 434 "dgges3.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 435 "dgges3.f"
    safmax = 1. / safmin;
#line 436 "dgges3.f"
    dlabad_(&safmin, &safmax);
#line 437 "dgges3.f"
    smlnum = sqrt(safmin) / eps;
#line 438 "dgges3.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 442 "dgges3.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 443 "dgges3.f"
    ilascl = FALSE_;
#line 444 "dgges3.f"
    if (anrm > 0. && anrm < smlnum) {
#line 445 "dgges3.f"
	anrmto = smlnum;
#line 446 "dgges3.f"
	ilascl = TRUE_;
#line 447 "dgges3.f"
    } else if (anrm > bignum) {
#line 448 "dgges3.f"
	anrmto = bignum;
#line 449 "dgges3.f"
	ilascl = TRUE_;
#line 450 "dgges3.f"
    }
#line 451 "dgges3.f"
    if (ilascl) {
#line 451 "dgges3.f"
	dlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 451 "dgges3.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 456 "dgges3.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 457 "dgges3.f"
    ilbscl = FALSE_;
#line 458 "dgges3.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 459 "dgges3.f"
	bnrmto = smlnum;
#line 460 "dgges3.f"
	ilbscl = TRUE_;
#line 461 "dgges3.f"
    } else if (bnrm > bignum) {
#line 462 "dgges3.f"
	bnrmto = bignum;
#line 463 "dgges3.f"
	ilbscl = TRUE_;
#line 464 "dgges3.f"
    }
#line 465 "dgges3.f"
    if (ilbscl) {
#line 465 "dgges3.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 465 "dgges3.f"
    }

/*     Permute the matrix to make it more nearly triangular */

#line 470 "dgges3.f"
    ileft = 1;
#line 471 "dgges3.f"
    iright = *n + 1;
#line 472 "dgges3.f"
    iwrk = iright + *n;
#line 473 "dgges3.f"
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */

#line 478 "dgges3.f"
    irows = ihi + 1 - ilo;
#line 479 "dgges3.f"
    icols = *n + 1 - ilo;
#line 480 "dgges3.f"
    itau = iwrk;
#line 481 "dgges3.f"
    iwrk = itau + irows;
#line 482 "dgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 482 "dgges3.f"
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */

#line 487 "dgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 487 "dgges3.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */

#line 493 "dgges3.f"
    if (ilvsl) {
#line 494 "dgges3.f"
	dlaset_("Full", n, n, &c_b36, &c_b37, &vsl[vsl_offset], ldvsl, (
		ftnlen)4);
#line 495 "dgges3.f"
	if (irows > 1) {
#line 496 "dgges3.f"
	    i__1 = irows - 1;
#line 496 "dgges3.f"
	    i__2 = irows - 1;
#line 496 "dgges3.f"
	    dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 498 "dgges3.f"
	}
#line 499 "dgges3.f"
	i__1 = *lwork + 1 - iwrk;
#line 499 "dgges3.f"
	dorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 501 "dgges3.f"
    }

/*     Initialize VSR */

#line 505 "dgges3.f"
    if (ilvsr) {
#line 505 "dgges3.f"
	dlaset_("Full", n, n, &c_b36, &c_b37, &vsr[vsr_offset], ldvsr, (
		ftnlen)4);
#line 505 "dgges3.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 510 "dgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 510 "dgges3.f"
    dgghd3_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk]
	    , &i__1, &ierr, (ftnlen)1, (ftnlen)1);

/*     Perform QZ algorithm, computing Schur vectors if desired */

#line 516 "dgges3.f"
    iwrk = itau;
#line 517 "dgges3.f"
    i__1 = *lwork + 1 - iwrk;
#line 517 "dgges3.f"
    dhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
	    , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 520 "dgges3.f"
    if (ierr != 0) {
#line 521 "dgges3.f"
	if (ierr > 0 && ierr <= *n) {
#line 522 "dgges3.f"
	    *info = ierr;
#line 523 "dgges3.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 524 "dgges3.f"
	    *info = ierr - *n;
#line 525 "dgges3.f"
	} else {
#line 526 "dgges3.f"
	    *info = *n + 1;
#line 527 "dgges3.f"
	}
#line 528 "dgges3.f"
	goto L50;
#line 529 "dgges3.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */

#line 533 "dgges3.f"
    *sdim = 0;
#line 534 "dgges3.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 538 "dgges3.f"
	if (ilascl) {
#line 539 "dgges3.f"
	    dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], 
		    n, &ierr, (ftnlen)1);
#line 541 "dgges3.f"
	    dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], 
		    n, &ierr, (ftnlen)1);
#line 543 "dgges3.f"
	}
#line 544 "dgges3.f"
	if (ilbscl) {
#line 544 "dgges3.f"
	    dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 544 "dgges3.f"
	}

/*        Select eigenvalues */

#line 549 "dgges3.f"
	i__1 = *n;
#line 549 "dgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 550 "dgges3.f"
	    bwork[i__] = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 551 "dgges3.f"
/* L10: */
#line 551 "dgges3.f"
	}

#line 553 "dgges3.f"
	i__1 = *lwork - iwrk + 1;
#line 553 "dgges3.f"
	dtgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, &
		pvsr, dif, &work[iwrk], &i__1, idum, &c__1, &ierr);
#line 557 "dgges3.f"
	if (ierr == 1) {
#line 557 "dgges3.f"
	    *info = *n + 3;
#line 557 "dgges3.f"
	}

#line 560 "dgges3.f"
    }

/*     Apply back-permutation to VSL and VSR */

#line 564 "dgges3.f"
    if (ilvsl) {
#line 564 "dgges3.f"
	dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 564 "dgges3.f"
    }

#line 568 "dgges3.f"
    if (ilvsr) {
#line 568 "dgges3.f"
	dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 568 "dgges3.f"
    }

/*     Check if unscaling would cause over/underflow, if so, rescale */
/*     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of */
/*     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I) */

#line 576 "dgges3.f"
    if (ilascl) {
#line 577 "dgges3.f"
	i__1 = *n;
#line 577 "dgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 578 "dgges3.f"
	    if (alphai[i__] != 0.) {
#line 579 "dgges3.f"
		if (alphar[i__] / safmax > anrmto / anrm || safmin / alphar[
			i__] > anrm / anrmto) {
#line 581 "dgges3.f"
		    work[1] = (d__1 = a[i__ + i__ * a_dim1] / alphar[i__], 
			    abs(d__1));
#line 582 "dgges3.f"
		    beta[i__] *= work[1];
#line 583 "dgges3.f"
		    alphar[i__] *= work[1];
#line 584 "dgges3.f"
		    alphai[i__] *= work[1];
#line 585 "dgges3.f"
		} else if (alphai[i__] / safmax > anrmto / anrm || safmin / 
			alphai[i__] > anrm / anrmto) {
#line 589 "dgges3.f"
		    work[1] = (d__1 = a[i__ + (i__ + 1) * a_dim1] / alphai[
			    i__], abs(d__1));
#line 590 "dgges3.f"
		    beta[i__] *= work[1];
#line 591 "dgges3.f"
		    alphar[i__] *= work[1];
#line 592 "dgges3.f"
		    alphai[i__] *= work[1];
#line 593 "dgges3.f"
		}
#line 594 "dgges3.f"
	    }
#line 595 "dgges3.f"
/* L20: */
#line 595 "dgges3.f"
	}
#line 596 "dgges3.f"
    }

#line 598 "dgges3.f"
    if (ilbscl) {
#line 599 "dgges3.f"
	i__1 = *n;
#line 599 "dgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 600 "dgges3.f"
	    if (alphai[i__] != 0.) {
#line 601 "dgges3.f"
		if (beta[i__] / safmax > bnrmto / bnrm || safmin / beta[i__] 
			> bnrm / bnrmto) {
#line 603 "dgges3.f"
		    work[1] = (d__1 = b[i__ + i__ * b_dim1] / beta[i__], abs(
			    d__1));
#line 604 "dgges3.f"
		    beta[i__] *= work[1];
#line 605 "dgges3.f"
		    alphar[i__] *= work[1];
#line 606 "dgges3.f"
		    alphai[i__] *= work[1];
#line 607 "dgges3.f"
		}
#line 608 "dgges3.f"
	    }
#line 609 "dgges3.f"
/* L30: */
#line 609 "dgges3.f"
	}
#line 610 "dgges3.f"
    }

/*     Undo scaling */

#line 614 "dgges3.f"
    if (ilascl) {
#line 615 "dgges3.f"
	dlascl_("H", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 616 "dgges3.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 617 "dgges3.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 618 "dgges3.f"
    }

#line 620 "dgges3.f"
    if (ilbscl) {
#line 621 "dgges3.f"
	dlascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 622 "dgges3.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 623 "dgges3.f"
    }

#line 625 "dgges3.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 629 "dgges3.f"
	lastsl = TRUE_;
#line 630 "dgges3.f"
	lst2sl = TRUE_;
#line 631 "dgges3.f"
	*sdim = 0;
#line 632 "dgges3.f"
	ip = 0;
#line 633 "dgges3.f"
	i__1 = *n;
#line 633 "dgges3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 634 "dgges3.f"
	    cursl = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 635 "dgges3.f"
	    if (alphai[i__] == 0.) {
#line 636 "dgges3.f"
		if (cursl) {
#line 636 "dgges3.f"
		    ++(*sdim);
#line 636 "dgges3.f"
		}
#line 638 "dgges3.f"
		ip = 0;
#line 639 "dgges3.f"
		if (cursl && ! lastsl) {
#line 639 "dgges3.f"
		    *info = *n + 2;
#line 639 "dgges3.f"
		}
#line 641 "dgges3.f"
	    } else {
#line 642 "dgges3.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 646 "dgges3.f"
		    cursl = cursl || lastsl;
#line 647 "dgges3.f"
		    lastsl = cursl;
#line 648 "dgges3.f"
		    if (cursl) {
#line 648 "dgges3.f"
			*sdim += 2;
#line 648 "dgges3.f"
		    }
#line 650 "dgges3.f"
		    ip = -1;
#line 651 "dgges3.f"
		    if (cursl && ! lst2sl) {
#line 651 "dgges3.f"
			*info = *n + 2;
#line 651 "dgges3.f"
		    }
#line 653 "dgges3.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 657 "dgges3.f"
		    ip = 1;
#line 658 "dgges3.f"
		}
#line 659 "dgges3.f"
	    }
#line 660 "dgges3.f"
	    lst2sl = lastsl;
#line 661 "dgges3.f"
	    lastsl = cursl;
#line 662 "dgges3.f"
/* L40: */
#line 662 "dgges3.f"
	}

#line 664 "dgges3.f"
    }

#line 666 "dgges3.f"
L50:

#line 668 "dgges3.f"
    work[1] = (doublereal) lwkopt;

#line 670 "dgges3.f"
    return 0;

/*     End of DGGES3 */

} /* dgges3_ */

