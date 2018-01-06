#line 1 "zggesx.f"
/* zggesx.f -- translated by f2c (version 20100827).
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

#line 1 "zggesx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> ZGGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGESX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggesx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggesx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggesx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, */
/*                          B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, */
/*                          LDVSR, RCONDE, RCONDV, WORK, LWORK, RWORK, */
/*                          IWORK, LIWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR, SENSE, SORT */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N, */
/*      $                   SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RCONDE( 2 ), RCONDV( 2 ), RWORK( * ) */
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
/* > ZGGESX computes for a pair of N-by-N complex nonsymmetric matrices */
/* > (A,B), the generalized eigenvalues, the complex Schur form (S,T), */
/* > and, optionally, the left and/or right matrices of Schur vectors (VSL */
/* > and VSR).  This gives the generalized Schur factorization */
/* > */
/* >      (A,B) = ( (VSL) S (VSR)**H, (VSL) T (VSR)**H ) */
/* > */
/* > where (VSR)**H is the conjugate-transpose of VSR. */
/* > */
/* > Optionally, it also orders the eigenvalues so that a selected cluster */
/* > of eigenvalues appears in the leading diagonal blocks of the upper */
/* > triangular matrix S and the upper triangular matrix T; computes */
/* > a reciprocal condition number for the average of the selected */
/* > eigenvalues (RCONDE); and computes a reciprocal condition number for */
/* > the right and left deflating subspaces corresponding to the selected */
/* > eigenvalues (RCONDV). The leading columns of VSL and VSR then form */
/* > an orthonormal basis for the corresponding left and right eigenspaces */
/* > (deflating subspaces). */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar w */
/* > or a ratio alpha/beta = w, such that  A - w*B is singular.  It is */
/* > usually represented as the pair (alpha,beta), as there is a */
/* > reasonable interpretation for beta=0 or for both being zero. */
/* > */
/* > A pair of matrices (S,T) is in generalized complex Schur form if T is */
/* > upper triangular with non-negative diagonal and S is upper */
/* > triangular. */
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
/* >          SELCTG is procedure) LOGICAL FUNCTION of two COMPLEX*16 arguments */
/* >          SELCTG must be declared EXTERNAL in the calling subroutine. */
/* >          If SORT = 'N', SELCTG is not referenced. */
/* >          If SORT = 'S', SELCTG is used to select eigenvalues to sort */
/* >          to the top left of the Schur form. */
/* >          Note that a selected complex eigenvalue may no longer satisfy */
/* >          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since */
/* >          ordering may change the value of complex eigenvalues */
/* >          (especially if the eigenvalue is ill-conditioned), in this */
/* >          case INFO is set to N+3 see INFO below). */
/* > \endverbatim */
/* > */
/* > \param[in] SENSE */
/* > \verbatim */
/* >          SENSE is CHARACTER*1 */
/* >          Determines which reciprocal condition numbers are computed. */
/* >          = 'N' : None are computed; */
/* >          = 'E' : Computed for average of selected eigenvalues only; */
/* >          = 'V' : Computed for selected deflating subspaces only; */
/* >          = 'B' : Computed for both. */
/* >          If SENSE = 'E', 'V', or 'B', SORT must equal 'S'. */
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
/* >          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the */
/* >          generalized eigenvalues.  ALPHA(j) and BETA(j),j=1,...,N  are */
/* >          the diagonals of the complex Schur form (S,T).  BETA(j) will */
/* >          be non-negative real. */
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
/* >          The leading dimension of the matrix VSL. LDVSL >=1, and */
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
/* > \param[out] RCONDE */
/* > \verbatim */
/* >          RCONDE is DOUBLE PRECISION array, dimension ( 2 ) */
/* >          If SENSE = 'E' or 'B', RCONDE(1) and RCONDE(2) contain the */
/* >          reciprocal condition numbers for the average of the selected */
/* >          eigenvalues. */
/* >          Not referenced if SENSE = 'N' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is DOUBLE PRECISION array, dimension ( 2 ) */
/* >          If SENSE = 'V' or 'B', RCONDV(1) and RCONDV(2) contain the */
/* >          reciprocal condition number for the selected deflating */
/* >          subspaces. */
/* >          Not referenced if SENSE = 'N' or 'E'. */
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
/* >          If N = 0, LWORK >= 1, else if SENSE = 'E', 'V', or 'B', */
/* >          LWORK >= MAX(1,2*N,2*SDIM*(N-SDIM)), else */
/* >          LWORK >= MAX(1,2*N).  Note that 2*SDIM*(N-SDIM) <= N*N/2. */
/* >          Note also that an error is only returned if */
/* >          LWORK < MAX(1,2*N), but if SENSE = 'E' or 'V' or 'B' this may */
/* >          not be large enough. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the bound on the optimal size of the WORK */
/* >          array and the minimum size of the IWORK array, returns these */
/* >          values as the first entries of the WORK and IWORK arrays, and */
/* >          no error message related to LWORK or LIWORK is issued by */
/* >          XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension ( 8*N ) */
/* >          Real workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. */
/* >          If SENSE = 'N' or N = 0, LIWORK >= 1, otherwise */
/* >          LIWORK >= N+2. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the bound on the optimal size of the */
/* >          WORK array and the minimum size of the IWORK array, returns */
/* >          these values as the first entries of the WORK and IWORK */
/* >          arrays, and no error message related to LWORK or LIWORK is */
/* >          issued by XERBLA. */
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

/* > \date November 2011 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha, 
	doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, 
	doublecomplex *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, integer *liwork, logical *bwork, integer *info, 
	ftnlen jobvsl_len, ftnlen jobvsr_len, ftnlen sort_len, ftnlen 
	sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal pl, pr, dif[2];
    static integer ihi, ilo;
    static doublereal eps;
    static integer ijob;
    static doublereal anrm, bnrm;
    static integer ierr, itau, iwrk, lwrk;
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
    static logical wantsb;
    static integer liwmin;
    static logical wantse, lastsl;
    static doublereal anrmto, bnrmto;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static integer maxwrk;
    static logical wantsn;
    static integer minwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int zhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen, ftnlen), zlacpy_(char *,
	     integer *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, ftnlen), zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static logical wantst, lquery, wantsv;
    extern /* Subroutine */ int ztgsen_(integer *, logical *, logical *, 
	    logical *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, integer *,
	     integer *, integer *), zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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

#line 395 "zggesx.f"
    /* Parameter adjustments */
#line 395 "zggesx.f"
    a_dim1 = *lda;
#line 395 "zggesx.f"
    a_offset = 1 + a_dim1;
#line 395 "zggesx.f"
    a -= a_offset;
#line 395 "zggesx.f"
    b_dim1 = *ldb;
#line 395 "zggesx.f"
    b_offset = 1 + b_dim1;
#line 395 "zggesx.f"
    b -= b_offset;
#line 395 "zggesx.f"
    --alpha;
#line 395 "zggesx.f"
    --beta;
#line 395 "zggesx.f"
    vsl_dim1 = *ldvsl;
#line 395 "zggesx.f"
    vsl_offset = 1 + vsl_dim1;
#line 395 "zggesx.f"
    vsl -= vsl_offset;
#line 395 "zggesx.f"
    vsr_dim1 = *ldvsr;
#line 395 "zggesx.f"
    vsr_offset = 1 + vsr_dim1;
#line 395 "zggesx.f"
    vsr -= vsr_offset;
#line 395 "zggesx.f"
    --rconde;
#line 395 "zggesx.f"
    --rcondv;
#line 395 "zggesx.f"
    --work;
#line 395 "zggesx.f"
    --rwork;
#line 395 "zggesx.f"
    --iwork;
#line 395 "zggesx.f"
    --bwork;
#line 395 "zggesx.f"

#line 395 "zggesx.f"
    /* Function Body */
#line 395 "zggesx.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 396 "zggesx.f"
	ijobvl = 1;
#line 397 "zggesx.f"
	ilvsl = FALSE_;
#line 398 "zggesx.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 399 "zggesx.f"
	ijobvl = 2;
#line 400 "zggesx.f"
	ilvsl = TRUE_;
#line 401 "zggesx.f"
    } else {
#line 402 "zggesx.f"
	ijobvl = -1;
#line 403 "zggesx.f"
	ilvsl = FALSE_;
#line 404 "zggesx.f"
    }

#line 406 "zggesx.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 407 "zggesx.f"
	ijobvr = 1;
#line 408 "zggesx.f"
	ilvsr = FALSE_;
#line 409 "zggesx.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 410 "zggesx.f"
	ijobvr = 2;
#line 411 "zggesx.f"
	ilvsr = TRUE_;
#line 412 "zggesx.f"
    } else {
#line 413 "zggesx.f"
	ijobvr = -1;
#line 414 "zggesx.f"
	ilvsr = FALSE_;
#line 415 "zggesx.f"
    }

#line 417 "zggesx.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 418 "zggesx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 419 "zggesx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 420 "zggesx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 421 "zggesx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 422 "zggesx.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 423 "zggesx.f"
    if (wantsn) {
#line 424 "zggesx.f"
	ijob = 0;
#line 425 "zggesx.f"
    } else if (wantse) {
#line 426 "zggesx.f"
	ijob = 1;
#line 427 "zggesx.f"
    } else if (wantsv) {
#line 428 "zggesx.f"
	ijob = 2;
#line 429 "zggesx.f"
    } else if (wantsb) {
#line 430 "zggesx.f"
	ijob = 4;
#line 431 "zggesx.f"
    }

/*     Test the input arguments */

#line 435 "zggesx.f"
    *info = 0;
#line 436 "zggesx.f"
    if (ijobvl <= 0) {
#line 437 "zggesx.f"
	*info = -1;
#line 438 "zggesx.f"
    } else if (ijobvr <= 0) {
#line 439 "zggesx.f"
	*info = -2;
#line 440 "zggesx.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 441 "zggesx.f"
	*info = -3;
#line 442 "zggesx.f"
    } else if (! (wantsn || wantse || wantsv || wantsb) || ! wantst && ! 
	    wantsn) {
#line 444 "zggesx.f"
	*info = -5;
#line 445 "zggesx.f"
    } else if (*n < 0) {
#line 446 "zggesx.f"
	*info = -6;
#line 447 "zggesx.f"
    } else if (*lda < max(1,*n)) {
#line 448 "zggesx.f"
	*info = -8;
#line 449 "zggesx.f"
    } else if (*ldb < max(1,*n)) {
#line 450 "zggesx.f"
	*info = -10;
#line 451 "zggesx.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 452 "zggesx.f"
	*info = -15;
#line 453 "zggesx.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 454 "zggesx.f"
	*info = -17;
#line 455 "zggesx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 464 "zggesx.f"
    if (*info == 0) {
#line 465 "zggesx.f"
	if (*n > 0) {
#line 466 "zggesx.f"
	    minwrk = *n << 1;
#line 467 "zggesx.f"
	    maxwrk = *n * (ilaenv_(&c__1, "ZGEQRF", " ", n, &c__1, n, &c__0, (
		    ftnlen)6, (ftnlen)1) + 1);
/* Computing MAX */
#line 468 "zggesx.f"
	    i__1 = maxwrk, i__2 = *n * (ilaenv_(&c__1, "ZUNMQR", " ", n, &
		    c__1, n, &c_n1, (ftnlen)6, (ftnlen)1) + 1);
#line 468 "zggesx.f"
	    maxwrk = max(i__1,i__2);
#line 470 "zggesx.f"
	    if (ilvsl) {
/* Computing MAX */
#line 471 "zggesx.f"
		i__1 = maxwrk, i__2 = *n * (ilaenv_(&c__1, "ZUNGQR", " ", n, &
			c__1, n, &c_n1, (ftnlen)6, (ftnlen)1) + 1);
#line 471 "zggesx.f"
		maxwrk = max(i__1,i__2);
#line 473 "zggesx.f"
	    }
#line 474 "zggesx.f"
	    lwrk = maxwrk;
#line 475 "zggesx.f"
	    if (ijob >= 1) {
/* Computing MAX */
#line 475 "zggesx.f"
		i__1 = lwrk, i__2 = *n * *n / 2;
#line 475 "zggesx.f"
		lwrk = max(i__1,i__2);
#line 475 "zggesx.f"
	    }
#line 477 "zggesx.f"
	} else {
#line 478 "zggesx.f"
	    minwrk = 1;
#line 479 "zggesx.f"
	    maxwrk = 1;
#line 480 "zggesx.f"
	    lwrk = 1;
#line 481 "zggesx.f"
	}
#line 482 "zggesx.f"
	work[1].r = (doublereal) lwrk, work[1].i = 0.;
#line 483 "zggesx.f"
	if (wantsn || *n == 0) {
#line 484 "zggesx.f"
	    liwmin = 1;
#line 485 "zggesx.f"
	} else {
#line 486 "zggesx.f"
	    liwmin = *n + 2;
#line 487 "zggesx.f"
	}
#line 488 "zggesx.f"
	iwork[1] = liwmin;

#line 490 "zggesx.f"
	if (*lwork < minwrk && ! lquery) {
#line 491 "zggesx.f"
	    *info = -21;
#line 492 "zggesx.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 493 "zggesx.f"
	    *info = -24;
#line 494 "zggesx.f"
	}
#line 495 "zggesx.f"
    }

#line 497 "zggesx.f"
    if (*info != 0) {
#line 498 "zggesx.f"
	i__1 = -(*info);
#line 498 "zggesx.f"
	xerbla_("ZGGESX", &i__1, (ftnlen)6);
#line 499 "zggesx.f"
	return 0;
#line 500 "zggesx.f"
    } else if (lquery) {
#line 501 "zggesx.f"
	return 0;
#line 502 "zggesx.f"
    }

/*     Quick return if possible */

#line 506 "zggesx.f"
    if (*n == 0) {
#line 507 "zggesx.f"
	*sdim = 0;
#line 508 "zggesx.f"
	return 0;
#line 509 "zggesx.f"
    }

/*     Get machine constants */

#line 513 "zggesx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 514 "zggesx.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 515 "zggesx.f"
    bignum = 1. / smlnum;
#line 516 "zggesx.f"
    dlabad_(&smlnum, &bignum);
#line 517 "zggesx.f"
    smlnum = sqrt(smlnum) / eps;
#line 518 "zggesx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 522 "zggesx.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 523 "zggesx.f"
    ilascl = FALSE_;
#line 524 "zggesx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 525 "zggesx.f"
	anrmto = smlnum;
#line 526 "zggesx.f"
	ilascl = TRUE_;
#line 527 "zggesx.f"
    } else if (anrm > bignum) {
#line 528 "zggesx.f"
	anrmto = bignum;
#line 529 "zggesx.f"
	ilascl = TRUE_;
#line 530 "zggesx.f"
    }
#line 531 "zggesx.f"
    if (ilascl) {
#line 531 "zggesx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 531 "zggesx.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 536 "zggesx.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 537 "zggesx.f"
    ilbscl = FALSE_;
#line 538 "zggesx.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 539 "zggesx.f"
	bnrmto = smlnum;
#line 540 "zggesx.f"
	ilbscl = TRUE_;
#line 541 "zggesx.f"
    } else if (bnrm > bignum) {
#line 542 "zggesx.f"
	bnrmto = bignum;
#line 543 "zggesx.f"
	ilbscl = TRUE_;
#line 544 "zggesx.f"
    }
#line 545 "zggesx.f"
    if (ilbscl) {
#line 545 "zggesx.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 545 "zggesx.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Real Workspace: need 6*N) */

#line 551 "zggesx.f"
    ileft = 1;
#line 552 "zggesx.f"
    iright = *n + 1;
#line 553 "zggesx.f"
    irwrk = iright + *n;
#line 554 "zggesx.f"
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 560 "zggesx.f"
    irows = ihi + 1 - ilo;
#line 561 "zggesx.f"
    icols = *n + 1 - ilo;
#line 562 "zggesx.f"
    itau = 1;
#line 563 "zggesx.f"
    iwrk = itau + irows;
#line 564 "zggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 564 "zggesx.f"
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the unitary transformation to matrix A */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 570 "zggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 570 "zggesx.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 577 "zggesx.f"
    if (ilvsl) {
#line 578 "zggesx.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 579 "zggesx.f"
	if (irows > 1) {
#line 580 "zggesx.f"
	    i__1 = irows - 1;
#line 580 "zggesx.f"
	    i__2 = irows - 1;
#line 580 "zggesx.f"
	    zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 582 "zggesx.f"
	}
#line 583 "zggesx.f"
	i__1 = *lwork + 1 - iwrk;
#line 583 "zggesx.f"
	zungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 585 "zggesx.f"
    }

/*     Initialize VSR */

#line 589 "zggesx.f"
    if (ilvsr) {
#line 589 "zggesx.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 589 "zggesx.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 595 "zggesx.f"
    zgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 598 "zggesx.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Complex Workspace: need N) */
/*     (Real Workspace:    need N) */

#line 604 "zggesx.f"
    iwrk = itau;
#line 605 "zggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 605 "zggesx.f"
    zhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 608 "zggesx.f"
    if (ierr != 0) {
#line 609 "zggesx.f"
	if (ierr > 0 && ierr <= *n) {
#line 610 "zggesx.f"
	    *info = ierr;
#line 611 "zggesx.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 612 "zggesx.f"
	    *info = ierr - *n;
#line 613 "zggesx.f"
	} else {
#line 614 "zggesx.f"
	    *info = *n + 1;
#line 615 "zggesx.f"
	}
#line 616 "zggesx.f"
	goto L40;
#line 617 "zggesx.f"
    }

/*     Sort eigenvalues ALPHA/BETA and compute the reciprocal of */
/*     condition number(s) */

#line 622 "zggesx.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 626 "zggesx.f"
	if (ilascl) {
#line 626 "zggesx.f"
	    zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n,
		     &ierr, (ftnlen)1);
#line 626 "zggesx.f"
	}
#line 628 "zggesx.f"
	if (ilbscl) {
#line 628 "zggesx.f"
	    zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 628 "zggesx.f"
	}

/*        Select eigenvalues */

#line 633 "zggesx.f"
	i__1 = *n;
#line 633 "zggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 634 "zggesx.f"
	    bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
#line 635 "zggesx.f"
/* L10: */
#line 635 "zggesx.f"
	}

/*        Reorder eigenvalues, transform Generalized Schur vectors, and */
/*        compute reciprocal condition numbers */
/*        (Complex Workspace: If IJOB >= 1, need MAX(1, 2*SDIM*(N-SDIM)) */
/*                            otherwise, need 1 ) */

#line 642 "zggesx.f"
	i__1 = *lwork - iwrk + 1;
#line 642 "zggesx.f"
	ztgsen_(&ijob, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, sdim, &pl, &pr, dif, &work[iwrk], &
		i__1, &iwork[1], liwork, &ierr);

#line 647 "zggesx.f"
	if (ijob >= 1) {
/* Computing MAX */
#line 647 "zggesx.f"
	    i__1 = maxwrk, i__2 = (*sdim << 1) * (*n - *sdim);
#line 647 "zggesx.f"
	    maxwrk = max(i__1,i__2);
#line 647 "zggesx.f"
	}
#line 649 "zggesx.f"
	if (ierr == -21) {

/*            not enough complex workspace */

#line 653 "zggesx.f"
	    *info = -21;
#line 654 "zggesx.f"
	} else {
#line 655 "zggesx.f"
	    if (ijob == 1 || ijob == 4) {
#line 656 "zggesx.f"
		rconde[1] = pl;
#line 657 "zggesx.f"
		rconde[2] = pr;
#line 658 "zggesx.f"
	    }
#line 659 "zggesx.f"
	    if (ijob == 2 || ijob == 4) {
#line 660 "zggesx.f"
		rcondv[1] = dif[0];
#line 661 "zggesx.f"
		rcondv[2] = dif[1];
#line 662 "zggesx.f"
	    }
#line 663 "zggesx.f"
	    if (ierr == 1) {
#line 663 "zggesx.f"
		*info = *n + 3;
#line 663 "zggesx.f"
	    }
#line 665 "zggesx.f"
	}

#line 667 "zggesx.f"
    }

/*     Apply permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 672 "zggesx.f"
    if (ilvsl) {
#line 672 "zggesx.f"
	zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 672 "zggesx.f"
    }

#line 676 "zggesx.f"
    if (ilvsr) {
#line 676 "zggesx.f"
	zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 676 "zggesx.f"
    }

/*     Undo scaling */

#line 682 "zggesx.f"
    if (ilascl) {
#line 683 "zggesx.f"
	zlascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 684 "zggesx.f"
	zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 685 "zggesx.f"
    }

#line 687 "zggesx.f"
    if (ilbscl) {
#line 688 "zggesx.f"
	zlascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 689 "zggesx.f"
	zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 690 "zggesx.f"
    }

#line 692 "zggesx.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 696 "zggesx.f"
	lastsl = TRUE_;
#line 697 "zggesx.f"
	*sdim = 0;
#line 698 "zggesx.f"
	i__1 = *n;
#line 698 "zggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 699 "zggesx.f"
	    cursl = (*selctg)(&alpha[i__], &beta[i__]);
#line 700 "zggesx.f"
	    if (cursl) {
#line 700 "zggesx.f"
		++(*sdim);
#line 700 "zggesx.f"
	    }
#line 702 "zggesx.f"
	    if (cursl && ! lastsl) {
#line 702 "zggesx.f"
		*info = *n + 2;
#line 702 "zggesx.f"
	    }
#line 704 "zggesx.f"
	    lastsl = cursl;
#line 705 "zggesx.f"
/* L30: */
#line 705 "zggesx.f"
	}

#line 707 "zggesx.f"
    }

#line 709 "zggesx.f"
L40:

#line 711 "zggesx.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 712 "zggesx.f"
    iwork[1] = liwmin;

#line 714 "zggesx.f"
    return 0;

/*     End of ZGGESX */

} /* zggesx_ */

