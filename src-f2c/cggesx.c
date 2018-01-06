#line 1 "cggesx.f"
/* cggesx.f -- translated by f2c (version 20100827).
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

#line 1 "cggesx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> CGGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGESX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggesx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggesx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggesx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, */
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
/*       REAL               RCONDE( 2 ), RCONDV( 2 ), RWORK( * ) */
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
/* > CGGESX computes for a pair of N-by-N complex nonsymmetric matrices */
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
/* >          SELCTG is procedure) LOGICAL FUNCTION of two COMPLEX arguments */
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
/* >          VSL is COMPLEX array, dimension (LDVSL,N) */
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
/* > \param[out] RCONDE */
/* > \verbatim */
/* >          RCONDE is REAL array, dimension ( 2 ) */
/* >          If SENSE = 'E' or 'B', RCONDE(1) and RCONDE(2) contain the */
/* >          reciprocal condition numbers for the average of the selected */
/* >          eigenvalues. */
/* >          Not referenced if SENSE = 'N' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is REAL array, dimension ( 2 ) */
/* >          If SENSE = 'V' or 'B', RCONDV(1) and RCONDV(2) contain the */
/* >          reciprocal condition number for the selected deflating */
/* >          subspaces. */
/* >          Not referenced if SENSE = 'N' or 'E'. */
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
/* >          RWORK is REAL array, dimension ( 8*N ) */
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
/* >          The dimension of the array WORK. */
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
/* Subroutine */ int cggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
	    ), clacpy_(char *, integer *, integer *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen), claset_(char *, integer *,
	     integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
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
    static logical wantsb;
    static integer liwmin;
    static logical wantse, lastsl;
    static doublereal anrmto, bnrmto;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static logical wantst, lquery, wantsv;


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

#line 395 "cggesx.f"
    /* Parameter adjustments */
#line 395 "cggesx.f"
    a_dim1 = *lda;
#line 395 "cggesx.f"
    a_offset = 1 + a_dim1;
#line 395 "cggesx.f"
    a -= a_offset;
#line 395 "cggesx.f"
    b_dim1 = *ldb;
#line 395 "cggesx.f"
    b_offset = 1 + b_dim1;
#line 395 "cggesx.f"
    b -= b_offset;
#line 395 "cggesx.f"
    --alpha;
#line 395 "cggesx.f"
    --beta;
#line 395 "cggesx.f"
    vsl_dim1 = *ldvsl;
#line 395 "cggesx.f"
    vsl_offset = 1 + vsl_dim1;
#line 395 "cggesx.f"
    vsl -= vsl_offset;
#line 395 "cggesx.f"
    vsr_dim1 = *ldvsr;
#line 395 "cggesx.f"
    vsr_offset = 1 + vsr_dim1;
#line 395 "cggesx.f"
    vsr -= vsr_offset;
#line 395 "cggesx.f"
    --rconde;
#line 395 "cggesx.f"
    --rcondv;
#line 395 "cggesx.f"
    --work;
#line 395 "cggesx.f"
    --rwork;
#line 395 "cggesx.f"
    --iwork;
#line 395 "cggesx.f"
    --bwork;
#line 395 "cggesx.f"

#line 395 "cggesx.f"
    /* Function Body */
#line 395 "cggesx.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 396 "cggesx.f"
	ijobvl = 1;
#line 397 "cggesx.f"
	ilvsl = FALSE_;
#line 398 "cggesx.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 399 "cggesx.f"
	ijobvl = 2;
#line 400 "cggesx.f"
	ilvsl = TRUE_;
#line 401 "cggesx.f"
    } else {
#line 402 "cggesx.f"
	ijobvl = -1;
#line 403 "cggesx.f"
	ilvsl = FALSE_;
#line 404 "cggesx.f"
    }

#line 406 "cggesx.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 407 "cggesx.f"
	ijobvr = 1;
#line 408 "cggesx.f"
	ilvsr = FALSE_;
#line 409 "cggesx.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 410 "cggesx.f"
	ijobvr = 2;
#line 411 "cggesx.f"
	ilvsr = TRUE_;
#line 412 "cggesx.f"
    } else {
#line 413 "cggesx.f"
	ijobvr = -1;
#line 414 "cggesx.f"
	ilvsr = FALSE_;
#line 415 "cggesx.f"
    }

#line 417 "cggesx.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 418 "cggesx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 419 "cggesx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 420 "cggesx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 421 "cggesx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 422 "cggesx.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 423 "cggesx.f"
    if (wantsn) {
#line 424 "cggesx.f"
	ijob = 0;
#line 425 "cggesx.f"
    } else if (wantse) {
#line 426 "cggesx.f"
	ijob = 1;
#line 427 "cggesx.f"
    } else if (wantsv) {
#line 428 "cggesx.f"
	ijob = 2;
#line 429 "cggesx.f"
    } else if (wantsb) {
#line 430 "cggesx.f"
	ijob = 4;
#line 431 "cggesx.f"
    }

/*     Test the input arguments */

#line 435 "cggesx.f"
    *info = 0;
#line 436 "cggesx.f"
    if (ijobvl <= 0) {
#line 437 "cggesx.f"
	*info = -1;
#line 438 "cggesx.f"
    } else if (ijobvr <= 0) {
#line 439 "cggesx.f"
	*info = -2;
#line 440 "cggesx.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 441 "cggesx.f"
	*info = -3;
#line 442 "cggesx.f"
    } else if (! (wantsn || wantse || wantsv || wantsb) || ! wantst && ! 
	    wantsn) {
#line 444 "cggesx.f"
	*info = -5;
#line 445 "cggesx.f"
    } else if (*n < 0) {
#line 446 "cggesx.f"
	*info = -6;
#line 447 "cggesx.f"
    } else if (*lda < max(1,*n)) {
#line 448 "cggesx.f"
	*info = -8;
#line 449 "cggesx.f"
    } else if (*ldb < max(1,*n)) {
#line 450 "cggesx.f"
	*info = -10;
#line 451 "cggesx.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 452 "cggesx.f"
	*info = -15;
#line 453 "cggesx.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 454 "cggesx.f"
	*info = -17;
#line 455 "cggesx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 464 "cggesx.f"
    if (*info == 0) {
#line 465 "cggesx.f"
	if (*n > 0) {
#line 466 "cggesx.f"
	    minwrk = *n << 1;
#line 467 "cggesx.f"
	    maxwrk = *n * (ilaenv_(&c__1, "CGEQRF", " ", n, &c__1, n, &c__0, (
		    ftnlen)6, (ftnlen)1) + 1);
/* Computing MAX */
#line 468 "cggesx.f"
	    i__1 = maxwrk, i__2 = *n * (ilaenv_(&c__1, "CUNMQR", " ", n, &
		    c__1, n, &c_n1, (ftnlen)6, (ftnlen)1) + 1);
#line 468 "cggesx.f"
	    maxwrk = max(i__1,i__2);
#line 470 "cggesx.f"
	    if (ilvsl) {
/* Computing MAX */
#line 471 "cggesx.f"
		i__1 = maxwrk, i__2 = *n * (ilaenv_(&c__1, "CUNGQR", " ", n, &
			c__1, n, &c_n1, (ftnlen)6, (ftnlen)1) + 1);
#line 471 "cggesx.f"
		maxwrk = max(i__1,i__2);
#line 473 "cggesx.f"
	    }
#line 474 "cggesx.f"
	    lwrk = maxwrk;
#line 475 "cggesx.f"
	    if (ijob >= 1) {
/* Computing MAX */
#line 475 "cggesx.f"
		i__1 = lwrk, i__2 = *n * *n / 2;
#line 475 "cggesx.f"
		lwrk = max(i__1,i__2);
#line 475 "cggesx.f"
	    }
#line 477 "cggesx.f"
	} else {
#line 478 "cggesx.f"
	    minwrk = 1;
#line 479 "cggesx.f"
	    maxwrk = 1;
#line 480 "cggesx.f"
	    lwrk = 1;
#line 481 "cggesx.f"
	}
#line 482 "cggesx.f"
	work[1].r = (doublereal) lwrk, work[1].i = 0.;
#line 483 "cggesx.f"
	if (wantsn || *n == 0) {
#line 484 "cggesx.f"
	    liwmin = 1;
#line 485 "cggesx.f"
	} else {
#line 486 "cggesx.f"
	    liwmin = *n + 2;
#line 487 "cggesx.f"
	}
#line 488 "cggesx.f"
	iwork[1] = liwmin;

#line 490 "cggesx.f"
	if (*lwork < minwrk && ! lquery) {
#line 491 "cggesx.f"
	    *info = -21;
#line 492 "cggesx.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 493 "cggesx.f"
	    *info = -24;
#line 494 "cggesx.f"
	}
#line 495 "cggesx.f"
    }

#line 497 "cggesx.f"
    if (*info != 0) {
#line 498 "cggesx.f"
	i__1 = -(*info);
#line 498 "cggesx.f"
	xerbla_("CGGESX", &i__1, (ftnlen)6);
#line 499 "cggesx.f"
	return 0;
#line 500 "cggesx.f"
    } else if (lquery) {
#line 501 "cggesx.f"
	return 0;
#line 502 "cggesx.f"
    }

/*     Quick return if possible */

#line 506 "cggesx.f"
    if (*n == 0) {
#line 507 "cggesx.f"
	*sdim = 0;
#line 508 "cggesx.f"
	return 0;
#line 509 "cggesx.f"
    }

/*     Get machine constants */

#line 513 "cggesx.f"
    eps = slamch_("P", (ftnlen)1);
#line 514 "cggesx.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 515 "cggesx.f"
    bignum = 1. / smlnum;
#line 516 "cggesx.f"
    slabad_(&smlnum, &bignum);
#line 517 "cggesx.f"
    smlnum = sqrt(smlnum) / eps;
#line 518 "cggesx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 522 "cggesx.f"
    anrm = clange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 523 "cggesx.f"
    ilascl = FALSE_;
#line 524 "cggesx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 525 "cggesx.f"
	anrmto = smlnum;
#line 526 "cggesx.f"
	ilascl = TRUE_;
#line 527 "cggesx.f"
    } else if (anrm > bignum) {
#line 528 "cggesx.f"
	anrmto = bignum;
#line 529 "cggesx.f"
	ilascl = TRUE_;
#line 530 "cggesx.f"
    }
#line 531 "cggesx.f"
    if (ilascl) {
#line 531 "cggesx.f"
	clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 531 "cggesx.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 536 "cggesx.f"
    bnrm = clange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 537 "cggesx.f"
    ilbscl = FALSE_;
#line 538 "cggesx.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 539 "cggesx.f"
	bnrmto = smlnum;
#line 540 "cggesx.f"
	ilbscl = TRUE_;
#line 541 "cggesx.f"
    } else if (bnrm > bignum) {
#line 542 "cggesx.f"
	bnrmto = bignum;
#line 543 "cggesx.f"
	ilbscl = TRUE_;
#line 544 "cggesx.f"
    }
#line 545 "cggesx.f"
    if (ilbscl) {
#line 545 "cggesx.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 545 "cggesx.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Real Workspace: need 6*N) */

#line 551 "cggesx.f"
    ileft = 1;
#line 552 "cggesx.f"
    iright = *n + 1;
#line 553 "cggesx.f"
    irwrk = iright + *n;
#line 554 "cggesx.f"
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 560 "cggesx.f"
    irows = ihi + 1 - ilo;
#line 561 "cggesx.f"
    icols = *n + 1 - ilo;
#line 562 "cggesx.f"
    itau = 1;
#line 563 "cggesx.f"
    iwrk = itau + irows;
#line 564 "cggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 564 "cggesx.f"
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the unitary transformation to matrix A */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 570 "cggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 570 "cggesx.f"
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 577 "cggesx.f"
    if (ilvsl) {
#line 578 "cggesx.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 579 "cggesx.f"
	if (irows > 1) {
#line 580 "cggesx.f"
	    i__1 = irows - 1;
#line 580 "cggesx.f"
	    i__2 = irows - 1;
#line 580 "cggesx.f"
	    clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 582 "cggesx.f"
	}
#line 583 "cggesx.f"
	i__1 = *lwork + 1 - iwrk;
#line 583 "cggesx.f"
	cungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 585 "cggesx.f"
    }

/*     Initialize VSR */

#line 589 "cggesx.f"
    if (ilvsr) {
#line 589 "cggesx.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 589 "cggesx.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 595 "cggesx.f"
    cgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 598 "cggesx.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Complex Workspace: need N) */
/*     (Real Workspace:    need N) */

#line 604 "cggesx.f"
    iwrk = itau;
#line 605 "cggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 605 "cggesx.f"
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 608 "cggesx.f"
    if (ierr != 0) {
#line 609 "cggesx.f"
	if (ierr > 0 && ierr <= *n) {
#line 610 "cggesx.f"
	    *info = ierr;
#line 611 "cggesx.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 612 "cggesx.f"
	    *info = ierr - *n;
#line 613 "cggesx.f"
	} else {
#line 614 "cggesx.f"
	    *info = *n + 1;
#line 615 "cggesx.f"
	}
#line 616 "cggesx.f"
	goto L40;
#line 617 "cggesx.f"
    }

/*     Sort eigenvalues ALPHA/BETA and compute the reciprocal of */
/*     condition number(s) */

#line 622 "cggesx.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 626 "cggesx.f"
	if (ilascl) {
#line 626 "cggesx.f"
	    clascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n,
		     &ierr, (ftnlen)1);
#line 626 "cggesx.f"
	}
#line 628 "cggesx.f"
	if (ilbscl) {
#line 628 "cggesx.f"
	    clascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 628 "cggesx.f"
	}

/*        Select eigenvalues */

#line 633 "cggesx.f"
	i__1 = *n;
#line 633 "cggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 634 "cggesx.f"
	    bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
#line 635 "cggesx.f"
/* L10: */
#line 635 "cggesx.f"
	}

/*        Reorder eigenvalues, transform Generalized Schur vectors, and */
/*        compute reciprocal condition numbers */
/*        (Complex Workspace: If IJOB >= 1, need MAX(1, 2*SDIM*(N-SDIM)) */
/*                            otherwise, need 1 ) */

#line 642 "cggesx.f"
	i__1 = *lwork - iwrk + 1;
#line 642 "cggesx.f"
	ctgsen_(&ijob, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, sdim, &pl, &pr, dif, &work[iwrk], &
		i__1, &iwork[1], liwork, &ierr);

#line 647 "cggesx.f"
	if (ijob >= 1) {
/* Computing MAX */
#line 647 "cggesx.f"
	    i__1 = maxwrk, i__2 = (*sdim << 1) * (*n - *sdim);
#line 647 "cggesx.f"
	    maxwrk = max(i__1,i__2);
#line 647 "cggesx.f"
	}
#line 649 "cggesx.f"
	if (ierr == -21) {

/*            not enough complex workspace */

#line 653 "cggesx.f"
	    *info = -21;
#line 654 "cggesx.f"
	} else {
#line 655 "cggesx.f"
	    if (ijob == 1 || ijob == 4) {
#line 656 "cggesx.f"
		rconde[1] = pl;
#line 657 "cggesx.f"
		rconde[2] = pr;
#line 658 "cggesx.f"
	    }
#line 659 "cggesx.f"
	    if (ijob == 2 || ijob == 4) {
#line 660 "cggesx.f"
		rcondv[1] = dif[0];
#line 661 "cggesx.f"
		rcondv[2] = dif[1];
#line 662 "cggesx.f"
	    }
#line 663 "cggesx.f"
	    if (ierr == 1) {
#line 663 "cggesx.f"
		*info = *n + 3;
#line 663 "cggesx.f"
	    }
#line 665 "cggesx.f"
	}

#line 667 "cggesx.f"
    }

/*     Apply permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 672 "cggesx.f"
    if (ilvsl) {
#line 672 "cggesx.f"
	cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 672 "cggesx.f"
    }

#line 676 "cggesx.f"
    if (ilvsr) {
#line 676 "cggesx.f"
	cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 676 "cggesx.f"
    }

/*     Undo scaling */

#line 682 "cggesx.f"
    if (ilascl) {
#line 683 "cggesx.f"
	clascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 684 "cggesx.f"
	clascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 685 "cggesx.f"
    }

#line 687 "cggesx.f"
    if (ilbscl) {
#line 688 "cggesx.f"
	clascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 689 "cggesx.f"
	clascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 690 "cggesx.f"
    }

#line 692 "cggesx.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 696 "cggesx.f"
	lastsl = TRUE_;
#line 697 "cggesx.f"
	*sdim = 0;
#line 698 "cggesx.f"
	i__1 = *n;
#line 698 "cggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 699 "cggesx.f"
	    cursl = (*selctg)(&alpha[i__], &beta[i__]);
#line 700 "cggesx.f"
	    if (cursl) {
#line 700 "cggesx.f"
		++(*sdim);
#line 700 "cggesx.f"
	    }
#line 702 "cggesx.f"
	    if (cursl && ! lastsl) {
#line 702 "cggesx.f"
		*info = *n + 2;
#line 702 "cggesx.f"
	    }
#line 704 "cggesx.f"
	    lastsl = cursl;
#line 705 "cggesx.f"
/* L30: */
#line 705 "cggesx.f"
	}

#line 707 "cggesx.f"
    }

#line 709 "cggesx.f"
L40:

#line 711 "cggesx.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 712 "cggesx.f"
    iwork[1] = liwmin;

#line 714 "cggesx.f"
    return 0;

/*     End of CGGESX */

} /* cggesx_ */

