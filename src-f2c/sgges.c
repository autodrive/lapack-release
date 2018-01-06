#line 1 "sgges.f"
/* sgges.f -- translated by f2c (version 20100827).
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

#line 1 "sgges.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
static doublereal c_b38 = 0.;
static doublereal c_b39 = 1.;

/* > \brief <b> SGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f
or GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgges.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgges.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgges.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, */
/*                         SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, */
/*                         LDVSR, WORK, LWORK, BWORK, INFO ) */

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
/* > SGGES computes for a pair of N-by-N real nonsymmetric matrices (A,B), */
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
/* >          If N = 0, LWORK >= 1, else LWORK >= max(8*N,6*N+16). */
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

/* > \date November 2011 */

/* > \ingroup realGEeigen */

/*  ===================================================================== */
/* Subroutine */ int sgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
    extern /* Subroutine */ int sgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical wantst, lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 344 "sgges.f"
    /* Parameter adjustments */
#line 344 "sgges.f"
    a_dim1 = *lda;
#line 344 "sgges.f"
    a_offset = 1 + a_dim1;
#line 344 "sgges.f"
    a -= a_offset;
#line 344 "sgges.f"
    b_dim1 = *ldb;
#line 344 "sgges.f"
    b_offset = 1 + b_dim1;
#line 344 "sgges.f"
    b -= b_offset;
#line 344 "sgges.f"
    --alphar;
#line 344 "sgges.f"
    --alphai;
#line 344 "sgges.f"
    --beta;
#line 344 "sgges.f"
    vsl_dim1 = *ldvsl;
#line 344 "sgges.f"
    vsl_offset = 1 + vsl_dim1;
#line 344 "sgges.f"
    vsl -= vsl_offset;
#line 344 "sgges.f"
    vsr_dim1 = *ldvsr;
#line 344 "sgges.f"
    vsr_offset = 1 + vsr_dim1;
#line 344 "sgges.f"
    vsr -= vsr_offset;
#line 344 "sgges.f"
    --work;
#line 344 "sgges.f"
    --bwork;
#line 344 "sgges.f"

#line 344 "sgges.f"
    /* Function Body */
#line 344 "sgges.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 345 "sgges.f"
	ijobvl = 1;
#line 346 "sgges.f"
	ilvsl = FALSE_;
#line 347 "sgges.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 348 "sgges.f"
	ijobvl = 2;
#line 349 "sgges.f"
	ilvsl = TRUE_;
#line 350 "sgges.f"
    } else {
#line 351 "sgges.f"
	ijobvl = -1;
#line 352 "sgges.f"
	ilvsl = FALSE_;
#line 353 "sgges.f"
    }

#line 355 "sgges.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 356 "sgges.f"
	ijobvr = 1;
#line 357 "sgges.f"
	ilvsr = FALSE_;
#line 358 "sgges.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 359 "sgges.f"
	ijobvr = 2;
#line 360 "sgges.f"
	ilvsr = TRUE_;
#line 361 "sgges.f"
    } else {
#line 362 "sgges.f"
	ijobvr = -1;
#line 363 "sgges.f"
	ilvsr = FALSE_;
#line 364 "sgges.f"
    }

#line 366 "sgges.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 370 "sgges.f"
    *info = 0;
#line 371 "sgges.f"
    lquery = *lwork == -1;
#line 372 "sgges.f"
    if (ijobvl <= 0) {
#line 373 "sgges.f"
	*info = -1;
#line 374 "sgges.f"
    } else if (ijobvr <= 0) {
#line 375 "sgges.f"
	*info = -2;
#line 376 "sgges.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 377 "sgges.f"
	*info = -3;
#line 378 "sgges.f"
    } else if (*n < 0) {
#line 379 "sgges.f"
	*info = -5;
#line 380 "sgges.f"
    } else if (*lda < max(1,*n)) {
#line 381 "sgges.f"
	*info = -7;
#line 382 "sgges.f"
    } else if (*ldb < max(1,*n)) {
#line 383 "sgges.f"
	*info = -9;
#line 384 "sgges.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 385 "sgges.f"
	*info = -15;
#line 386 "sgges.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 387 "sgges.f"
	*info = -17;
#line 388 "sgges.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 397 "sgges.f"
    if (*info == 0) {
#line 398 "sgges.f"
	if (*n > 0) {
/* Computing MAX */
#line 399 "sgges.f"
	    i__1 = *n << 3, i__2 = *n * 6 + 16;
#line 399 "sgges.f"
	    minwrk = max(i__1,i__2);
#line 400 "sgges.f"
	    maxwrk = minwrk - *n + *n * ilaenv_(&c__1, "SGEQRF", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 402 "sgges.f"
	    i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "SORMQR", 
		    " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 402 "sgges.f"
	    maxwrk = max(i__1,i__2);
#line 404 "sgges.f"
	    if (ilvsl) {
/* Computing MAX */
#line 405 "sgges.f"
		i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "SOR"\
			"GQR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 405 "sgges.f"
		maxwrk = max(i__1,i__2);
#line 407 "sgges.f"
	    }
#line 408 "sgges.f"
	} else {
#line 409 "sgges.f"
	    minwrk = 1;
#line 410 "sgges.f"
	    maxwrk = 1;
#line 411 "sgges.f"
	}
#line 412 "sgges.f"
	work[1] = (doublereal) maxwrk;

#line 414 "sgges.f"
	if (*lwork < minwrk && ! lquery) {
#line 414 "sgges.f"
	    *info = -19;
#line 414 "sgges.f"
	}
#line 416 "sgges.f"
    }

#line 418 "sgges.f"
    if (*info != 0) {
#line 419 "sgges.f"
	i__1 = -(*info);
#line 419 "sgges.f"
	xerbla_("SGGES ", &i__1, (ftnlen)6);
#line 420 "sgges.f"
	return 0;
#line 421 "sgges.f"
    } else if (lquery) {
#line 422 "sgges.f"
	return 0;
#line 423 "sgges.f"
    }

/*     Quick return if possible */

#line 427 "sgges.f"
    if (*n == 0) {
#line 428 "sgges.f"
	*sdim = 0;
#line 429 "sgges.f"
	return 0;
#line 430 "sgges.f"
    }

/*     Get machine constants */

#line 434 "sgges.f"
    eps = slamch_("P", (ftnlen)1);
#line 435 "sgges.f"
    safmin = slamch_("S", (ftnlen)1);
#line 436 "sgges.f"
    safmax = 1. / safmin;
#line 437 "sgges.f"
    slabad_(&safmin, &safmax);
#line 438 "sgges.f"
    smlnum = sqrt(safmin) / eps;
#line 439 "sgges.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 443 "sgges.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 444 "sgges.f"
    ilascl = FALSE_;
#line 445 "sgges.f"
    if (anrm > 0. && anrm < smlnum) {
#line 446 "sgges.f"
	anrmto = smlnum;
#line 447 "sgges.f"
	ilascl = TRUE_;
#line 448 "sgges.f"
    } else if (anrm > bignum) {
#line 449 "sgges.f"
	anrmto = bignum;
#line 450 "sgges.f"
	ilascl = TRUE_;
#line 451 "sgges.f"
    }
#line 452 "sgges.f"
    if (ilascl) {
#line 452 "sgges.f"
	slascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 452 "sgges.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 457 "sgges.f"
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 458 "sgges.f"
    ilbscl = FALSE_;
#line 459 "sgges.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 460 "sgges.f"
	bnrmto = smlnum;
#line 461 "sgges.f"
	ilbscl = TRUE_;
#line 462 "sgges.f"
    } else if (bnrm > bignum) {
#line 463 "sgges.f"
	bnrmto = bignum;
#line 464 "sgges.f"
	ilbscl = TRUE_;
#line 465 "sgges.f"
    }
#line 466 "sgges.f"
    if (ilbscl) {
#line 466 "sgges.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 466 "sgges.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Workspace: need 6*N + 2*N space for storing balancing factors) */

#line 472 "sgges.f"
    ileft = 1;
#line 473 "sgges.f"
    iright = *n + 1;
#line 474 "sgges.f"
    iwrk = iright + *n;
#line 475 "sgges.f"
    sggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Workspace: need N, prefer N*NB) */

#line 481 "sgges.f"
    irows = ihi + 1 - ilo;
#line 482 "sgges.f"
    icols = *n + 1 - ilo;
#line 483 "sgges.f"
    itau = iwrk;
#line 484 "sgges.f"
    iwrk = itau + irows;
#line 485 "sgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 485 "sgges.f"
    sgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Workspace: need N, prefer N*NB) */

#line 491 "sgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 491 "sgges.f"
    sormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Workspace: need N, prefer N*NB) */

#line 498 "sgges.f"
    if (ilvsl) {
#line 499 "sgges.f"
	slaset_("Full", n, n, &c_b38, &c_b39, &vsl[vsl_offset], ldvsl, (
		ftnlen)4);
#line 500 "sgges.f"
	if (irows > 1) {
#line 501 "sgges.f"
	    i__1 = irows - 1;
#line 501 "sgges.f"
	    i__2 = irows - 1;
#line 501 "sgges.f"
	    slacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 503 "sgges.f"
	}
#line 504 "sgges.f"
	i__1 = *lwork + 1 - iwrk;
#line 504 "sgges.f"
	sorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 506 "sgges.f"
    }

/*     Initialize VSR */

#line 510 "sgges.f"
    if (ilvsr) {
#line 510 "sgges.f"
	slaset_("Full", n, n, &c_b38, &c_b39, &vsr[vsr_offset], ldvsr, (
		ftnlen)4);
#line 510 "sgges.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 516 "sgges.f"
    sgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Workspace: need N) */

#line 522 "sgges.f"
    iwrk = itau;
#line 523 "sgges.f"
    i__1 = *lwork + 1 - iwrk;
#line 523 "sgges.f"
    shgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
	    , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 526 "sgges.f"
    if (ierr != 0) {
#line 527 "sgges.f"
	if (ierr > 0 && ierr <= *n) {
#line 528 "sgges.f"
	    *info = ierr;
#line 529 "sgges.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 530 "sgges.f"
	    *info = ierr - *n;
#line 531 "sgges.f"
	} else {
#line 532 "sgges.f"
	    *info = *n + 1;
#line 533 "sgges.f"
	}
#line 534 "sgges.f"
	goto L40;
#line 535 "sgges.f"
    }

/*     Sort eigenvalues ALPHA/BETA if desired */
/*     (Workspace: need 4*N+16 ) */

#line 540 "sgges.f"
    *sdim = 0;
#line 541 "sgges.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 545 "sgges.f"
	if (ilascl) {
#line 546 "sgges.f"
	    slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], 
		    n, &ierr, (ftnlen)1);
#line 548 "sgges.f"
	    slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], 
		    n, &ierr, (ftnlen)1);
#line 550 "sgges.f"
	}
#line 551 "sgges.f"
	if (ilbscl) {
#line 551 "sgges.f"
	    slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 551 "sgges.f"
	}

/*        Select eigenvalues */

#line 556 "sgges.f"
	i__1 = *n;
#line 556 "sgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 557 "sgges.f"
	    bwork[i__] = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 558 "sgges.f"
/* L10: */
#line 558 "sgges.f"
	}

#line 560 "sgges.f"
	i__1 = *lwork - iwrk + 1;
#line 560 "sgges.f"
	stgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl, &
		pvsr, dif, &work[iwrk], &i__1, idum, &c__1, &ierr);
#line 564 "sgges.f"
	if (ierr == 1) {
#line 564 "sgges.f"
	    *info = *n + 3;
#line 564 "sgges.f"
	}

#line 567 "sgges.f"
    }

/*     Apply back-permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 572 "sgges.f"
    if (ilvsl) {
#line 572 "sgges.f"
	sggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 572 "sgges.f"
    }

#line 576 "sgges.f"
    if (ilvsr) {
#line 576 "sgges.f"
	sggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 576 "sgges.f"
    }

/*     Check if unscaling would cause over/underflow, if so, rescale */
/*     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of */
/*     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I) */

#line 584 "sgges.f"
    if (ilascl) {
#line 585 "sgges.f"
	i__1 = *n;
#line 585 "sgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 586 "sgges.f"
	    if (alphai[i__] != 0.) {
#line 587 "sgges.f"
		if (alphar[i__] / safmax > anrmto / anrm || safmin / alphar[
			i__] > anrm / anrmto) {
#line 589 "sgges.f"
		    work[1] = (d__1 = a[i__ + i__ * a_dim1] / alphar[i__], 
			    abs(d__1));
#line 590 "sgges.f"
		    beta[i__] *= work[1];
#line 591 "sgges.f"
		    alphar[i__] *= work[1];
#line 592 "sgges.f"
		    alphai[i__] *= work[1];
#line 593 "sgges.f"
		} else if (alphai[i__] / safmax > anrmto / anrm || safmin / 
			alphai[i__] > anrm / anrmto) {
#line 595 "sgges.f"
		    work[1] = (d__1 = a[i__ + (i__ + 1) * a_dim1] / alphai[
			    i__], abs(d__1));
#line 596 "sgges.f"
		    beta[i__] *= work[1];
#line 597 "sgges.f"
		    alphar[i__] *= work[1];
#line 598 "sgges.f"
		    alphai[i__] *= work[1];
#line 599 "sgges.f"
		}
#line 600 "sgges.f"
	    }
#line 601 "sgges.f"
/* L50: */
#line 601 "sgges.f"
	}
#line 602 "sgges.f"
    }

#line 604 "sgges.f"
    if (ilbscl) {
#line 605 "sgges.f"
	i__1 = *n;
#line 605 "sgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 606 "sgges.f"
	    if (alphai[i__] != 0.) {
#line 607 "sgges.f"
		if (beta[i__] / safmax > bnrmto / bnrm || safmin / beta[i__] 
			> bnrm / bnrmto) {
#line 609 "sgges.f"
		    work[1] = (d__1 = b[i__ + i__ * b_dim1] / beta[i__], abs(
			    d__1));
#line 610 "sgges.f"
		    beta[i__] *= work[1];
#line 611 "sgges.f"
		    alphar[i__] *= work[1];
#line 612 "sgges.f"
		    alphai[i__] *= work[1];
#line 613 "sgges.f"
		}
#line 614 "sgges.f"
	    }
#line 615 "sgges.f"
/* L60: */
#line 615 "sgges.f"
	}
#line 616 "sgges.f"
    }

/*     Undo scaling */

#line 620 "sgges.f"
    if (ilascl) {
#line 621 "sgges.f"
	slascl_("H", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 622 "sgges.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 623 "sgges.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 624 "sgges.f"
    }

#line 626 "sgges.f"
    if (ilbscl) {
#line 627 "sgges.f"
	slascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 628 "sgges.f"
	slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 629 "sgges.f"
    }

#line 631 "sgges.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 635 "sgges.f"
	lastsl = TRUE_;
#line 636 "sgges.f"
	lst2sl = TRUE_;
#line 637 "sgges.f"
	*sdim = 0;
#line 638 "sgges.f"
	ip = 0;
#line 639 "sgges.f"
	i__1 = *n;
#line 639 "sgges.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 640 "sgges.f"
	    cursl = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 641 "sgges.f"
	    if (alphai[i__] == 0.) {
#line 642 "sgges.f"
		if (cursl) {
#line 642 "sgges.f"
		    ++(*sdim);
#line 642 "sgges.f"
		}
#line 644 "sgges.f"
		ip = 0;
#line 645 "sgges.f"
		if (cursl && ! lastsl) {
#line 645 "sgges.f"
		    *info = *n + 2;
#line 645 "sgges.f"
		}
#line 647 "sgges.f"
	    } else {
#line 648 "sgges.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 652 "sgges.f"
		    cursl = cursl || lastsl;
#line 653 "sgges.f"
		    lastsl = cursl;
#line 654 "sgges.f"
		    if (cursl) {
#line 654 "sgges.f"
			*sdim += 2;
#line 654 "sgges.f"
		    }
#line 656 "sgges.f"
		    ip = -1;
#line 657 "sgges.f"
		    if (cursl && ! lst2sl) {
#line 657 "sgges.f"
			*info = *n + 2;
#line 657 "sgges.f"
		    }
#line 659 "sgges.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 663 "sgges.f"
		    ip = 1;
#line 664 "sgges.f"
		}
#line 665 "sgges.f"
	    }
#line 666 "sgges.f"
	    lst2sl = lastsl;
#line 667 "sgges.f"
	    lastsl = cursl;
#line 668 "sgges.f"
/* L30: */
#line 668 "sgges.f"
	}

#line 670 "sgges.f"
    }

#line 672 "sgges.f"
L40:

#line 674 "sgges.f"
    work[1] = (doublereal) maxwrk;

#line 676 "sgges.f"
    return 0;

/*     End of SGGES */

} /* sgges_ */

