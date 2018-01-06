#line 1 "dgges.f"
/* dgges.f -- translated by f2c (version 20100827).
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

#line 1 "dgges.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
static doublereal c_b38 = 0.;
static doublereal c_b39 = 1.;

/* > \brief <b> DGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f
or GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgges.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgges.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgges.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, */
/*                         SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, */
/*                         LDVSR, WORK, LWORK, BWORK, INFO ) */

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
/* > DGGES computes for a pair of N-by-N real nonsymmetric matrices (A,B), */
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
/* >          If N = 0, LWORK >= 1, else LWORK >= 8*N+16. */
/* >          For good performance , LWORK must generally be larger. */
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

/* > \date November 2011 */

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dggbak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dggbal_(char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical lst2sl;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlascl_(char *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
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
    static integer ijobvl, iright;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal anrmto, bnrmto;
    static logical lastsl;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static logical wantst, lquery;


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 344 "dgges.f"
    /* Parameter adjustments */
#line 344 "dgges.f"
    a_dim1 = *lda;
#line 344 "dgges.f"
    a_offset = 1 + a_dim1;
#line 344 "dgges.f"
    a -= a_offset;
#line 344 "dgges.f"
    b_dim1 = *ldb;
#line 344 "dgges.f"
    b_offset = 1 + b_dim1;
#line 344 "dgges.f"
    b -= b_offset;
#line 344 "dgges.f"
    --alphar;
#line 344 "dgges.f"
    --alphai;
#line 344 "dgges.f"
    --beta;
#line 344 "dgges.f"
    vsl_dim1 = *ldvsl;
#line 344 "dgges.f"
    vsl_offset = 1 + vsl_dim1;
#line 344 "dgges.f"
    vsl -= vsl_offset;
#line 344 "dgges.f"
    vsr_dim1 = *ldvsr;
#line 344 "dgges.f"
    vsr_offset = 1 + vsr_dim1;
#line 344 "dgges.f"
    vsr -= vsr_offset;
#line 344 "dgges.f"
    --work;
#line 344 "dgges.f"
    --bwork;
#line 344 "dgges.f"

#line 344 "dgges.f"
    /* Function Body */
#line 344 "dgges.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 345 "dgges.f"
	ijobvl = 1;
#line 346 "dgges.f"
	ilvsl = FALSE_;
#line 347 "dgges.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 348 "dgges.f"
	ijobvl = 2;
#line 349 "dgges.f"
	ilvsl = TRUE_;
#line 350 "dgges.f"
    } else {
#line 351 "dgges.f"
	ijobvl = -1;
#line 352 "dgges.f"
	ilvsl = FALSE_;
#line 353 "dgges.f"
    }

#line 355 "dgges.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 356 "dgges.f"
	ijobvr = 1;
#line 357 "dgges.f"
	ilvsr = FALSE_;
#line 358 "dgges.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 359 "dgges.f"
	ijobvr = 2;
#line 360 "dgges.f"
	ilvsr = TRUE_;
#line 361 "dgges.f"
    } else {
#line 362 "dgges.f"
	ijobvr = -1;
#line 363 "dgges.f"
	ilvsr = FALSE_;
#line 364 "dgges.f"
    }

#line 366 "dgges.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 370 "dgges.f"
    *info = 0;
#line 371 "dgges.f"
    lquery = *lwork == -1;
#line 372 "dgges.f"
    if (ijobvl <= 0) {
#line 373 "dgges.f"
	*info = -1;
#line 374 "dgges.f"
    } else if (ijobvr <= 0) {
#line 375 "dgges.f"
	*info = -2;
#line 376 "dgges.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 377 "dgges.f"
	*info = -3;
#line 378 "dgges.f"
    } else if (*n < 0) {
#line 379 "dgges.f"
	*info = -5;
#line 380 "dgges.f"
    } else if (*lda < max(1,*n)) {
#line 381 "dgges.f"
	*info = -7;
#line 382 "dgges.f"
    } else if (*ldb < max(1,*n)) {
#line 383 "dgges.f"
	*info = -9;
#line 384 "dgges.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 385 "dgges.f"
	*info = -15;
#line 386 "dgges.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 387 "dgges.f"
	*info = -17;
#line 388 "dgges.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 397 "dgges.f"
    if (*info == 0) {
#line 398 "dgges.f"
	if (*n > 0) {
/* Computing MAX */
#line 399 "dgges.f"
	    i__1 = *n << 3, i__2 = *n * 6 + 16;
#line 399 "dgges.f"
	    minwrk = max(i__1,i__2);
#line 400 "dgges.f"
	    maxwrk = minwrk - *n + *n * ilaenv_(&c__1, "DGEQRF", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 402 "dgges.f"
	    i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "DORMQR", 
		    " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 402 "dgges.f"
	    maxwrk = max(i__1,i__2);
#line 404 "dgges.f"
	    if (ilvsl) {
/* Computing MAX */
#line 405 "dgges.f"
		i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "DOR"\
			"GQR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 405 "dgges.f"
		maxwrk = max(i__1,i__2);
#line 407 "dgges.f"
	    }
#line 408 "dgges.f"
	} else {
#line 409 "dgges.f"
	    minwrk = 1;
#line 410 "dgges.f"
	    maxwrk = 1;
#line 411 "dgges.f"
	}
#line 412 "dgges.f"
	work[1] = (doublereal) maxwrk;

#line 414 "dgges.f"
	if (*lwork < minwrk && ! lquery) {
#line 414 "dgges.f"
	    *info = -19;
#line 414 "dgges.f"
	}
#line 416 "dgges.f"
    }

#line 418 "dgges.f"
    if (*info != 0) {
#line 419 "dgges.f"
	i__1 = -(*info);
#line 419 "dgges.f"
	xerbla_("DGGES ", &i__1, (ftnlen)6);
#line 420 "dgges.f"
	return 0;
#line 421 "dgges.f"
    } else if (lquery) {
#line 422 "dgges.f"
	return 0;
#line 423 "dgges.f"
    }

/*     Quick return if possible */

#line 427 "dgges.f"
    if (*n == 0) {
#line 428 "dgges.f"
	*sdim = 0;
#line 429 "dgges.f"
	return 0;
#line 430 "dgges.f"
    }

/*     Get machine constants */

#line 434 "dgges.f"
    eps = dlamch_("P", (ftnlen)1);
#line 435 "dgges.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 436 "dgges.f"
    safmax = 1. / safmin;
#line 437 "dgges.f"
    dlabad_(&safmin, &safmax);
#line 438 "dgges.f"
    smlnum = sqrt(safmin) / eps;
#line 439 "dgges.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 443 "dgges.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 444 "dgges.f"
    ilascl = FALSE_;
#line 445 "dgges.f"
    if (anrm > 0. && anrm < smlnum) {
#line 446 "dgges.f"
	anrmto = smlnum;
#line 447 "dgges.f"
	ilascl = TRUE_;
#line 448 "dgges.f"
    } else if (anrm > bignum) {
#line 449 "dgges.f"
	anrmto = bignum;
#line 450 "dgges.f"
	ilascl = TRUE_;
#line 451 "dgges.f"
    }
#line 452 "dgges.f"
    if (ilascl) {
#line 452 "dgges.f"
	dlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 452 "dgges.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 457 "dgges.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 458 "dgges.f"
    ilbscl = FALSE_;
#line 459 "dgges.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 460 "dgges.f"
	bnrmto = smlnum;
#line 461 "dgges.f"
	ilbscl = TRUE_;
#line 462 "dgges.f"
    } else if (bnrm > bignum) {
#line 463 "dgges.f"
	bnrmto = bignum;
#line 464 "dgges.f"
	ilbscl = TRUE_;
#line 465 "dgges.f"
    }
#line 466 "dgges.f"
    if (ilbscl) {
#line 466 "dgges.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 466 "dgges.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Workspace: need 6*N + 2*N space for storing balancing factors) */

#line 472 "dgges.f"
    ileft = 1;
#line 473 "dgges.f"
    iright = *n + 1;
#line 474 "dgges.f"
    iwrk = iright + *n;
#line 475 "dgges.f"
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Workspace: need N, prefer N*NB) */

#line 481 "dgges.f"
    irows = ihi + 1 - ilo;
#line 482 "dgges.f"
    icols = *n + 1 - ilo;
#line 483 "dgges.f"
    itau = iwrk;
#line 484 "dgges.f"
    iwrk = itau + irows;
#line 485 "dgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 485 "dgges.f"
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Workspace: need N, prefer N*NB) */

#line 491 "dgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 491 "dgges.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Workspace: need N, prefer N*NB) */

#line 498 "dgges.f"
    if (ilvsl) {
#line 499 "dgges.f"
	dlaset_("Full", n, n, &c_b38, &c_b39, &vsl[vsl_offset], ldvsl, (
		ftnlen)4);
#line 500 "dgges.f"
	if (irows > 1) {
#line 501 "dgges.f"
	    i__1 = irows - 1;
#line 501 "dgges.f"
	    i__2 = irows - 1;
#line 501 "dgges.f"
	    dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 503 "dgges.f"
	}
#line 504 "dgges.f"
	i__1 = *lwork + 1 - iwrk;
#line 504 "dgges.f"
	dorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 506 "dgges.f"
    }

/*     Initialize VSR */

#line 510 "dgges.f"
    if (ilvsr) {
#line 510 "dgges.f"
	dlaset_("Full", n, n, &c_b38, &c_b39, &vsr[vsr_offset], ldvsr, (
		ftnlen)4);
#line 510 "dgges.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 516 "dgges.f"
    dgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Workspace: need N) */

#line 522 "dgges.f"
    iwrk = itau;
#line 523 "dgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 523 "dgges.f"
    dhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
	    , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 526 "dgges.f"
    if (ierr != 0) {
#line 527 "dgges.f"
	if (ierr > 0 && ierr <= *n) {
#line 528 "dgges.f"
	    *info = ierr;
#line 529 "dgges.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 530 "dgges.f"
	    *info = ierr - *n;
#line 531 "dgges.f"
	} else {
#line 532 "dgges.f"
	    *info = *n + 1;
#line 533 "dgges.f"
	}
#line 534 "dgges.f"
	goto L50;
#line 535 "dgges.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */
/*     (Workspace: need 4*N+16 ) */

#line 540 "dgges.f"
    *sdim = 0;
#line 541 "dgges.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 545 "dgges.f"
	if (ilascl) {
#line 546 "dgges.f"
	    dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], 
		    n, &ierr, (ftnlen)1);
#line 548 "dgges.f"
	    dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], 
		    n, &ierr, (ftnlen)1);
#line 550 "dgges.f"
	}
#line 551 "dgges.f"
	if (ilbscl) {
#line 551 "dgges.f"
	    dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 551 "dgges.f"
	}

/*        Select eigenvalues */

#line 556 "dgges.f"
	i__1 = *n;
#line 556 "dgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 557 "dgges.f"
	    bwork[i__] = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 558 "dgges.f"
/* L10: */
#line 558 "dgges.f"
	}

#line 560 "dgges.f"
	i__1 = *lwork - iwrk + 1;
#line 560 "dgges.f"
	dtgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, &
		pvsr, dif, &work[iwrk], &i__1, idum, &c__1, &ierr);
#line 564 "dgges.f"
	if (ierr == 1) {
#line 564 "dgges.f"
	    *info = *n + 3;
#line 564 "dgges.f"
	}

#line 567 "dgges.f"
    }

/*     Apply back-permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 572 "dgges.f"
    if (ilvsl) {
#line 572 "dgges.f"
	dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 572 "dgges.f"
    }

#line 576 "dgges.f"
    if (ilvsr) {
#line 576 "dgges.f"
	dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 576 "dgges.f"
    }

/*     Check if unscaling would cause over/underflow, if so, rescale */
/*     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of */
/*     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I) */

#line 584 "dgges.f"
    if (ilascl) {
#line 585 "dgges.f"
	i__1 = *n;
#line 585 "dgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 586 "dgges.f"
	    if (alphai[i__] != 0.) {
#line 587 "dgges.f"
		if (alphar[i__] / safmax > anrmto / anrm || safmin / alphar[
			i__] > anrm / anrmto) {
#line 589 "dgges.f"
		    work[1] = (d__1 = a[i__ + i__ * a_dim1] / alphar[i__], 
			    abs(d__1));
#line 590 "dgges.f"
		    beta[i__] *= work[1];
#line 591 "dgges.f"
		    alphar[i__] *= work[1];
#line 592 "dgges.f"
		    alphai[i__] *= work[1];
#line 593 "dgges.f"
		} else if (alphai[i__] / safmax > anrmto / anrm || safmin / 
			alphai[i__] > anrm / anrmto) {
#line 597 "dgges.f"
		    work[1] = (d__1 = a[i__ + (i__ + 1) * a_dim1] / alphai[
			    i__], abs(d__1));
#line 598 "dgges.f"
		    beta[i__] *= work[1];
#line 599 "dgges.f"
		    alphar[i__] *= work[1];
#line 600 "dgges.f"
		    alphai[i__] *= work[1];
#line 601 "dgges.f"
		}
#line 602 "dgges.f"
	    }
#line 603 "dgges.f"
/* L20: */
#line 603 "dgges.f"
	}
#line 604 "dgges.f"
    }

#line 606 "dgges.f"
    if (ilbscl) {
#line 607 "dgges.f"
	i__1 = *n;
#line 607 "dgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 608 "dgges.f"
	    if (alphai[i__] != 0.) {
#line 609 "dgges.f"
		if (beta[i__] / safmax > bnrmto / bnrm || safmin / beta[i__] 
			> bnrm / bnrmto) {
#line 611 "dgges.f"
		    work[1] = (d__1 = b[i__ + i__ * b_dim1] / beta[i__], abs(
			    d__1));
#line 612 "dgges.f"
		    beta[i__] *= work[1];
#line 613 "dgges.f"
		    alphar[i__] *= work[1];
#line 614 "dgges.f"
		    alphai[i__] *= work[1];
#line 615 "dgges.f"
		}
#line 616 "dgges.f"
	    }
#line 617 "dgges.f"
/* L30: */
#line 617 "dgges.f"
	}
#line 618 "dgges.f"
    }

/*     Undo scaling */

#line 622 "dgges.f"
    if (ilascl) {
#line 623 "dgges.f"
	dlascl_("H", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 624 "dgges.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 625 "dgges.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 626 "dgges.f"
    }

#line 628 "dgges.f"
    if (ilbscl) {
#line 629 "dgges.f"
	dlascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 630 "dgges.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 631 "dgges.f"
    }

#line 633 "dgges.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 637 "dgges.f"
	lastsl = TRUE_;
#line 638 "dgges.f"
	lst2sl = TRUE_;
#line 639 "dgges.f"
	*sdim = 0;
#line 640 "dgges.f"
	ip = 0;
#line 641 "dgges.f"
	i__1 = *n;
#line 641 "dgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 642 "dgges.f"
	    cursl = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 643 "dgges.f"
	    if (alphai[i__] == 0.) {
#line 644 "dgges.f"
		if (cursl) {
#line 644 "dgges.f"
		    ++(*sdim);
#line 644 "dgges.f"
		}
#line 646 "dgges.f"
		ip = 0;
#line 647 "dgges.f"
		if (cursl && ! lastsl) {
#line 647 "dgges.f"
		    *info = *n + 2;
#line 647 "dgges.f"
		}
#line 649 "dgges.f"
	    } else {
#line 650 "dgges.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 654 "dgges.f"
		    cursl = cursl || lastsl;
#line 655 "dgges.f"
		    lastsl = cursl;
#line 656 "dgges.f"
		    if (cursl) {
#line 656 "dgges.f"
			*sdim += 2;
#line 656 "dgges.f"
		    }
#line 658 "dgges.f"
		    ip = -1;
#line 659 "dgges.f"
		    if (cursl && ! lst2sl) {
#line 659 "dgges.f"
			*info = *n + 2;
#line 659 "dgges.f"
		    }
#line 661 "dgges.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 665 "dgges.f"
		    ip = 1;
#line 666 "dgges.f"
		}
#line 667 "dgges.f"
	    }
#line 668 "dgges.f"
	    lst2sl = lastsl;
#line 669 "dgges.f"
	    lastsl = cursl;
#line 670 "dgges.f"
/* L40: */
#line 670 "dgges.f"
	}

#line 672 "dgges.f"
    }

#line 674 "dgges.f"
L50:

#line 676 "dgges.f"
    work[1] = (doublereal) maxwrk;

#line 678 "dgges.f"
    return 0;

/*     End of DGGES */

} /* dgges_ */

