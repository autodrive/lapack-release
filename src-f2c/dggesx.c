#line 1 "dggesx.f"
/* dggesx.f -- translated by f2c (version 20100827).
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

#line 1 "dggesx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
static doublereal c_b42 = 0.;
static doublereal c_b43 = 1.;

/* > \brief <b> DGGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGESX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggesx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggesx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggesx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, */
/*                          B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, */
/*                          VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK, IWORK, */
/*                          LIWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR, SENSE, SORT */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N, */
/*      $                   SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), RCONDE( 2 ), */
/*      $                   RCONDV( 2 ), VSL( LDVSL, * ), VSR( LDVSR, * ), */
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
/* > DGGESX computes for a pair of N-by-N real nonsymmetric matrices */
/* > (A,B), the generalized eigenvalues, the real Schur form (S,T), and, */
/* > optionally, the left and/or right matrices of Schur vectors (VSL and */
/* > VSR).  This gives the generalized Schur factorization */
/* > */
/* >      (A,B) = ( (VSL) S (VSR)**T, (VSL) T (VSR)**T ) */
/* > */
/* > Optionally, it also orders the eigenvalues so that a selected cluster */
/* > of eigenvalues appears in the leading diagonal blocks of the upper */
/* > quasi-triangular matrix S and the upper triangular matrix T; computes */
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
/* >          = 'S':  Eigenvalues are ordered (see SELCTG). */
/* > \endverbatim */
/* > */
/* > \param[in] SELCTG */
/* > \verbatim */
/* >          SELCTG is procedure) LOGICAL FUNCTION of three DOUBLE PRECISION arguments */
/* >          SELCTG must be declared EXTERNAL in the calling subroutine. */
/* >          If SORT = 'N', SELCTG is not referenced. */
/* >          If SORT = 'S', SELCTG is used to select eigenvalues to sort */
/* >          to the top left of the Schur form. */
/* >          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if */
/* >          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either */
/* >          one of a complex conjugate pair of eigenvalues is selected, */
/* >          then both complex eigenvalues are selected. */
/* >          Note that a selected complex eigenvalue may no longer satisfy */
/* >          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) = .TRUE. after ordering, */
/* >          since ordering may change the value of complex eigenvalues */
/* >          (especially if the eigenvalue is ill-conditioned), in this */
/* >          case INFO is set to N+3. */
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
/* >          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i */
/* >          and BETA(j),j=1,...,N  are the diagonals of the complex Schur */
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
/* >          reciprocal condition numbers for the selected deflating */
/* >          subspaces. */
/* >          Not referenced if SENSE = 'N' or 'E'. */
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
/* >          If N = 0, LWORK >= 1, else if SENSE = 'E', 'V', or 'B', */
/* >          LWORK >= max( 8*N, 6*N+16, 2*SDIM*(N-SDIM) ), else */
/* >          LWORK >= max( 8*N, 6*N+16 ). */
/* >          Note that 2*SDIM*(N-SDIM) <= N*N/2. */
/* >          Note also that an error is only returned if */
/* >          LWORK < max( 8*N, 6*N+16), but if SENSE = 'E' or 'V' or 'B' */
/* >          this may not be large enough. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the bound on the optimal size of the WORK */
/* >          array and the minimum size of the IWORK array, returns these */
/* >          values as the first entries of the WORK and IWORK arrays, and */
/* >          no error message related to LWORK or LIWORK is issued by */
/* >          XERBLA. */
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
/* >          LIWORK >= N+6. */
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
/* >                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should */
/* >                be correct for j=INFO+1,...,N. */
/* >          > N:  =N+1: other than QZ iteration failed in DHGEQZ */
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

/* > \date December 2016 */

/* > \ingroup doubleGEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  An approximate (asymptotic) bound on the average absolute error of */
/* >  the selected eigenvalues is */
/* > */
/* >       EPS * norm((A, B)) / RCONDE( 1 ). */
/* > */
/* >  An approximate (asymptotic) bound on the maximum angular error in */
/* >  the computed deflating subspaces is */
/* > */
/* >       EPS * norm((A, B)) / RCONDV( 2 ). */
/* > */
/* >  See LAPACK User's Guide, section 4.11 for more information. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, integer *sdim, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl,
	 doublereal *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, logical *bwork, integer *info, ftnlen jobvsl_len, ftnlen 
	jobvsr_len, ftnlen sort_len, ftnlen sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ip;
    static doublereal pl, pr, dif[2];
    static integer ihi, ilo;
    static doublereal eps;
    static integer ijob;
    static doublereal anrm, bnrm;
    static integer ierr, itau, iwrk, lwrk;
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
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer ijobvl, iright;
    extern /* Subroutine */ int dtgsen_(integer *, logical *, logical *, 
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ijobvr;
    static logical wantsb;
    static integer liwmin;
    static logical wantse, lastsl;
    static doublereal anrmto, bnrmto;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
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

#line 428 "dggesx.f"
    /* Parameter adjustments */
#line 428 "dggesx.f"
    a_dim1 = *lda;
#line 428 "dggesx.f"
    a_offset = 1 + a_dim1;
#line 428 "dggesx.f"
    a -= a_offset;
#line 428 "dggesx.f"
    b_dim1 = *ldb;
#line 428 "dggesx.f"
    b_offset = 1 + b_dim1;
#line 428 "dggesx.f"
    b -= b_offset;
#line 428 "dggesx.f"
    --alphar;
#line 428 "dggesx.f"
    --alphai;
#line 428 "dggesx.f"
    --beta;
#line 428 "dggesx.f"
    vsl_dim1 = *ldvsl;
#line 428 "dggesx.f"
    vsl_offset = 1 + vsl_dim1;
#line 428 "dggesx.f"
    vsl -= vsl_offset;
#line 428 "dggesx.f"
    vsr_dim1 = *ldvsr;
#line 428 "dggesx.f"
    vsr_offset = 1 + vsr_dim1;
#line 428 "dggesx.f"
    vsr -= vsr_offset;
#line 428 "dggesx.f"
    --rconde;
#line 428 "dggesx.f"
    --rcondv;
#line 428 "dggesx.f"
    --work;
#line 428 "dggesx.f"
    --iwork;
#line 428 "dggesx.f"
    --bwork;
#line 428 "dggesx.f"

#line 428 "dggesx.f"
    /* Function Body */
#line 428 "dggesx.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 429 "dggesx.f"
	ijobvl = 1;
#line 430 "dggesx.f"
	ilvsl = FALSE_;
#line 431 "dggesx.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 432 "dggesx.f"
	ijobvl = 2;
#line 433 "dggesx.f"
	ilvsl = TRUE_;
#line 434 "dggesx.f"
    } else {
#line 435 "dggesx.f"
	ijobvl = -1;
#line 436 "dggesx.f"
	ilvsl = FALSE_;
#line 437 "dggesx.f"
    }

#line 439 "dggesx.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 440 "dggesx.f"
	ijobvr = 1;
#line 441 "dggesx.f"
	ilvsr = FALSE_;
#line 442 "dggesx.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 443 "dggesx.f"
	ijobvr = 2;
#line 444 "dggesx.f"
	ilvsr = TRUE_;
#line 445 "dggesx.f"
    } else {
#line 446 "dggesx.f"
	ijobvr = -1;
#line 447 "dggesx.f"
	ilvsr = FALSE_;
#line 448 "dggesx.f"
    }

#line 450 "dggesx.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 451 "dggesx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 452 "dggesx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 453 "dggesx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 454 "dggesx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 455 "dggesx.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 456 "dggesx.f"
    if (wantsn) {
#line 457 "dggesx.f"
	ijob = 0;
#line 458 "dggesx.f"
    } else if (wantse) {
#line 459 "dggesx.f"
	ijob = 1;
#line 460 "dggesx.f"
    } else if (wantsv) {
#line 461 "dggesx.f"
	ijob = 2;
#line 462 "dggesx.f"
    } else if (wantsb) {
#line 463 "dggesx.f"
	ijob = 4;
#line 464 "dggesx.f"
    }

/*     Test the input arguments */

#line 468 "dggesx.f"
    *info = 0;
#line 469 "dggesx.f"
    if (ijobvl <= 0) {
#line 470 "dggesx.f"
	*info = -1;
#line 471 "dggesx.f"
    } else if (ijobvr <= 0) {
#line 472 "dggesx.f"
	*info = -2;
#line 473 "dggesx.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 474 "dggesx.f"
	*info = -3;
#line 475 "dggesx.f"
    } else if (! (wantsn || wantse || wantsv || wantsb) || ! wantst && ! 
	    wantsn) {
#line 477 "dggesx.f"
	*info = -5;
#line 478 "dggesx.f"
    } else if (*n < 0) {
#line 479 "dggesx.f"
	*info = -6;
#line 480 "dggesx.f"
    } else if (*lda < max(1,*n)) {
#line 481 "dggesx.f"
	*info = -8;
#line 482 "dggesx.f"
    } else if (*ldb < max(1,*n)) {
#line 483 "dggesx.f"
	*info = -10;
#line 484 "dggesx.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 485 "dggesx.f"
	*info = -16;
#line 486 "dggesx.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 487 "dggesx.f"
	*info = -18;
#line 488 "dggesx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 497 "dggesx.f"
    if (*info == 0) {
#line 498 "dggesx.f"
	if (*n > 0) {
/* Computing MAX */
#line 499 "dggesx.f"
	    i__1 = *n << 3, i__2 = *n * 6 + 16;
#line 499 "dggesx.f"
	    minwrk = max(i__1,i__2);
#line 500 "dggesx.f"
	    maxwrk = minwrk - *n + *n * ilaenv_(&c__1, "DGEQRF", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 502 "dggesx.f"
	    i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "DORMQR", 
		    " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 502 "dggesx.f"
	    maxwrk = max(i__1,i__2);
#line 504 "dggesx.f"
	    if (ilvsl) {
/* Computing MAX */
#line 505 "dggesx.f"
		i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "DOR"\
			"GQR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 505 "dggesx.f"
		maxwrk = max(i__1,i__2);
#line 507 "dggesx.f"
	    }
#line 508 "dggesx.f"
	    lwrk = maxwrk;
#line 509 "dggesx.f"
	    if (ijob >= 1) {
/* Computing MAX */
#line 509 "dggesx.f"
		i__1 = lwrk, i__2 = *n * *n / 2;
#line 509 "dggesx.f"
		lwrk = max(i__1,i__2);
#line 509 "dggesx.f"
	    }
#line 511 "dggesx.f"
	} else {
#line 512 "dggesx.f"
	    minwrk = 1;
#line 513 "dggesx.f"
	    maxwrk = 1;
#line 514 "dggesx.f"
	    lwrk = 1;
#line 515 "dggesx.f"
	}
#line 516 "dggesx.f"
	work[1] = (doublereal) lwrk;
#line 517 "dggesx.f"
	if (wantsn || *n == 0) {
#line 518 "dggesx.f"
	    liwmin = 1;
#line 519 "dggesx.f"
	} else {
#line 520 "dggesx.f"
	    liwmin = *n + 6;
#line 521 "dggesx.f"
	}
#line 522 "dggesx.f"
	iwork[1] = liwmin;

#line 524 "dggesx.f"
	if (*lwork < minwrk && ! lquery) {
#line 525 "dggesx.f"
	    *info = -22;
#line 526 "dggesx.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 527 "dggesx.f"
	    *info = -24;
#line 528 "dggesx.f"
	}
#line 529 "dggesx.f"
    }

#line 531 "dggesx.f"
    if (*info != 0) {
#line 532 "dggesx.f"
	i__1 = -(*info);
#line 532 "dggesx.f"
	xerbla_("DGGESX", &i__1, (ftnlen)6);
#line 533 "dggesx.f"
	return 0;
#line 534 "dggesx.f"
    } else if (lquery) {
#line 535 "dggesx.f"
	return 0;
#line 536 "dggesx.f"
    }

/*     Quick return if possible */

#line 540 "dggesx.f"
    if (*n == 0) {
#line 541 "dggesx.f"
	*sdim = 0;
#line 542 "dggesx.f"
	return 0;
#line 543 "dggesx.f"
    }

/*     Get machine constants */

#line 547 "dggesx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 548 "dggesx.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 549 "dggesx.f"
    safmax = 1. / safmin;
#line 550 "dggesx.f"
    dlabad_(&safmin, &safmax);
#line 551 "dggesx.f"
    smlnum = sqrt(safmin) / eps;
#line 552 "dggesx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 556 "dggesx.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 557 "dggesx.f"
    ilascl = FALSE_;
#line 558 "dggesx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 559 "dggesx.f"
	anrmto = smlnum;
#line 560 "dggesx.f"
	ilascl = TRUE_;
#line 561 "dggesx.f"
    } else if (anrm > bignum) {
#line 562 "dggesx.f"
	anrmto = bignum;
#line 563 "dggesx.f"
	ilascl = TRUE_;
#line 564 "dggesx.f"
    }
#line 565 "dggesx.f"
    if (ilascl) {
#line 565 "dggesx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 565 "dggesx.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 570 "dggesx.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 571 "dggesx.f"
    ilbscl = FALSE_;
#line 572 "dggesx.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 573 "dggesx.f"
	bnrmto = smlnum;
#line 574 "dggesx.f"
	ilbscl = TRUE_;
#line 575 "dggesx.f"
    } else if (bnrm > bignum) {
#line 576 "dggesx.f"
	bnrmto = bignum;
#line 577 "dggesx.f"
	ilbscl = TRUE_;
#line 578 "dggesx.f"
    }
#line 579 "dggesx.f"
    if (ilbscl) {
#line 579 "dggesx.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 579 "dggesx.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Workspace: need 6*N + 2*N for permutation parameters) */

#line 585 "dggesx.f"
    ileft = 1;
#line 586 "dggesx.f"
    iright = *n + 1;
#line 587 "dggesx.f"
    iwrk = iright + *n;
#line 588 "dggesx.f"
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Workspace: need N, prefer N*NB) */

#line 594 "dggesx.f"
    irows = ihi + 1 - ilo;
#line 595 "dggesx.f"
    icols = *n + 1 - ilo;
#line 596 "dggesx.f"
    itau = iwrk;
#line 597 "dggesx.f"
    iwrk = itau + irows;
#line 598 "dggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 598 "dggesx.f"
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Workspace: need N, prefer N*NB) */

#line 604 "dggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 604 "dggesx.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Workspace: need N, prefer N*NB) */

#line 611 "dggesx.f"
    if (ilvsl) {
#line 612 "dggesx.f"
	dlaset_("Full", n, n, &c_b42, &c_b43, &vsl[vsl_offset], ldvsl, (
		ftnlen)4);
#line 613 "dggesx.f"
	if (irows > 1) {
#line 614 "dggesx.f"
	    i__1 = irows - 1;
#line 614 "dggesx.f"
	    i__2 = irows - 1;
#line 614 "dggesx.f"
	    dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 616 "dggesx.f"
	}
#line 617 "dggesx.f"
	i__1 = *lwork + 1 - iwrk;
#line 617 "dggesx.f"
	dorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 619 "dggesx.f"
    }

/*     Initialize VSR */

#line 623 "dggesx.f"
    if (ilvsr) {
#line 623 "dggesx.f"
	dlaset_("Full", n, n, &c_b42, &c_b43, &vsr[vsr_offset], ldvsr, (
		ftnlen)4);
#line 623 "dggesx.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 629 "dggesx.f"
    dgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 632 "dggesx.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Workspace: need N) */

#line 637 "dggesx.f"
    iwrk = itau;
#line 638 "dggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 638 "dggesx.f"
    dhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
	    , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 641 "dggesx.f"
    if (ierr != 0) {
#line 642 "dggesx.f"
	if (ierr > 0 && ierr <= *n) {
#line 643 "dggesx.f"
	    *info = ierr;
#line 644 "dggesx.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 645 "dggesx.f"
	    *info = ierr - *n;
#line 646 "dggesx.f"
	} else {
#line 647 "dggesx.f"
	    *info = *n + 1;
#line 648 "dggesx.f"
	}
#line 649 "dggesx.f"
	goto L60;
#line 650 "dggesx.f"
    }

/*     Sort eigenvalues ALPHA/BETA and compute the reciprocal of */
/*     condition number(s) */
/*     (Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) ) */
/*                 otherwise, need 8*(N+1) ) */

#line 657 "dggesx.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 661 "dggesx.f"
	if (ilascl) {
#line 662 "dggesx.f"
	    dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], 
		    n, &ierr, (ftnlen)1);
#line 664 "dggesx.f"
	    dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], 
		    n, &ierr, (ftnlen)1);
#line 666 "dggesx.f"
	}
#line 667 "dggesx.f"
	if (ilbscl) {
#line 667 "dggesx.f"
	    dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 667 "dggesx.f"
	}

/*        Select eigenvalues */

#line 672 "dggesx.f"
	i__1 = *n;
#line 672 "dggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 673 "dggesx.f"
	    bwork[i__] = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 674 "dggesx.f"
/* L10: */
#line 674 "dggesx.f"
	}

/*        Reorder eigenvalues, transform Generalized Schur vectors, and */
/*        compute reciprocal condition numbers */

#line 679 "dggesx.f"
	i__1 = *lwork - iwrk + 1;
#line 679 "dggesx.f"
	dtgsen_(&ijob, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pl, &pr, 
		dif, &work[iwrk], &i__1, &iwork[1], liwork, &ierr);

#line 684 "dggesx.f"
	if (ijob >= 1) {
/* Computing MAX */
#line 684 "dggesx.f"
	    i__1 = maxwrk, i__2 = (*sdim << 1) * (*n - *sdim);
#line 684 "dggesx.f"
	    maxwrk = max(i__1,i__2);
#line 684 "dggesx.f"
	}
#line 686 "dggesx.f"
	if (ierr == -22) {

/*            not enough real workspace */

#line 690 "dggesx.f"
	    *info = -22;
#line 691 "dggesx.f"
	} else {
#line 692 "dggesx.f"
	    if (ijob == 1 || ijob == 4) {
#line 693 "dggesx.f"
		rconde[1] = pl;
#line 694 "dggesx.f"
		rconde[2] = pr;
#line 695 "dggesx.f"
	    }
#line 696 "dggesx.f"
	    if (ijob == 2 || ijob == 4) {
#line 697 "dggesx.f"
		rcondv[1] = dif[0];
#line 698 "dggesx.f"
		rcondv[2] = dif[1];
#line 699 "dggesx.f"
	    }
#line 700 "dggesx.f"
	    if (ierr == 1) {
#line 700 "dggesx.f"
		*info = *n + 3;
#line 700 "dggesx.f"
	    }
#line 702 "dggesx.f"
	}

#line 704 "dggesx.f"
    }

/*     Apply permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 709 "dggesx.f"
    if (ilvsl) {
#line 709 "dggesx.f"
	dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 709 "dggesx.f"
    }

#line 713 "dggesx.f"
    if (ilvsr) {
#line 713 "dggesx.f"
	dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 713 "dggesx.f"
    }

/*     Check if unscaling would cause over/underflow, if so, rescale */
/*     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of */
/*     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I) */

#line 721 "dggesx.f"
    if (ilascl) {
#line 722 "dggesx.f"
	i__1 = *n;
#line 722 "dggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 723 "dggesx.f"
	    if (alphai[i__] != 0.) {
#line 724 "dggesx.f"
		if (alphar[i__] / safmax > anrmto / anrm || safmin / alphar[
			i__] > anrm / anrmto) {
#line 726 "dggesx.f"
		    work[1] = (d__1 = a[i__ + i__ * a_dim1] / alphar[i__], 
			    abs(d__1));
#line 727 "dggesx.f"
		    beta[i__] *= work[1];
#line 728 "dggesx.f"
		    alphar[i__] *= work[1];
#line 729 "dggesx.f"
		    alphai[i__] *= work[1];
#line 730 "dggesx.f"
		} else if (alphai[i__] / safmax > anrmto / anrm || safmin / 
			alphai[i__] > anrm / anrmto) {
#line 734 "dggesx.f"
		    work[1] = (d__1 = a[i__ + (i__ + 1) * a_dim1] / alphai[
			    i__], abs(d__1));
#line 735 "dggesx.f"
		    beta[i__] *= work[1];
#line 736 "dggesx.f"
		    alphar[i__] *= work[1];
#line 737 "dggesx.f"
		    alphai[i__] *= work[1];
#line 738 "dggesx.f"
		}
#line 739 "dggesx.f"
	    }
#line 740 "dggesx.f"
/* L20: */
#line 740 "dggesx.f"
	}
#line 741 "dggesx.f"
    }

#line 743 "dggesx.f"
    if (ilbscl) {
#line 744 "dggesx.f"
	i__1 = *n;
#line 744 "dggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 745 "dggesx.f"
	    if (alphai[i__] != 0.) {
#line 746 "dggesx.f"
		if (beta[i__] / safmax > bnrmto / bnrm || safmin / beta[i__] 
			> bnrm / bnrmto) {
#line 748 "dggesx.f"
		    work[1] = (d__1 = b[i__ + i__ * b_dim1] / beta[i__], abs(
			    d__1));
#line 749 "dggesx.f"
		    beta[i__] *= work[1];
#line 750 "dggesx.f"
		    alphar[i__] *= work[1];
#line 751 "dggesx.f"
		    alphai[i__] *= work[1];
#line 752 "dggesx.f"
		}
#line 753 "dggesx.f"
	    }
#line 754 "dggesx.f"
/* L30: */
#line 754 "dggesx.f"
	}
#line 755 "dggesx.f"
    }

/*     Undo scaling */

#line 759 "dggesx.f"
    if (ilascl) {
#line 760 "dggesx.f"
	dlascl_("H", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 761 "dggesx.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 762 "dggesx.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 763 "dggesx.f"
    }

#line 765 "dggesx.f"
    if (ilbscl) {
#line 766 "dggesx.f"
	dlascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 767 "dggesx.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 768 "dggesx.f"
    }

#line 770 "dggesx.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 774 "dggesx.f"
	lastsl = TRUE_;
#line 775 "dggesx.f"
	lst2sl = TRUE_;
#line 776 "dggesx.f"
	*sdim = 0;
#line 777 "dggesx.f"
	ip = 0;
#line 778 "dggesx.f"
	i__1 = *n;
#line 778 "dggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 779 "dggesx.f"
	    cursl = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 780 "dggesx.f"
	    if (alphai[i__] == 0.) {
#line 781 "dggesx.f"
		if (cursl) {
#line 781 "dggesx.f"
		    ++(*sdim);
#line 781 "dggesx.f"
		}
#line 783 "dggesx.f"
		ip = 0;
#line 784 "dggesx.f"
		if (cursl && ! lastsl) {
#line 784 "dggesx.f"
		    *info = *n + 2;
#line 784 "dggesx.f"
		}
#line 786 "dggesx.f"
	    } else {
#line 787 "dggesx.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 791 "dggesx.f"
		    cursl = cursl || lastsl;
#line 792 "dggesx.f"
		    lastsl = cursl;
#line 793 "dggesx.f"
		    if (cursl) {
#line 793 "dggesx.f"
			*sdim += 2;
#line 793 "dggesx.f"
		    }
#line 795 "dggesx.f"
		    ip = -1;
#line 796 "dggesx.f"
		    if (cursl && ! lst2sl) {
#line 796 "dggesx.f"
			*info = *n + 2;
#line 796 "dggesx.f"
		    }
#line 798 "dggesx.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 802 "dggesx.f"
		    ip = 1;
#line 803 "dggesx.f"
		}
#line 804 "dggesx.f"
	    }
#line 805 "dggesx.f"
	    lst2sl = lastsl;
#line 806 "dggesx.f"
	    lastsl = cursl;
#line 807 "dggesx.f"
/* L50: */
#line 807 "dggesx.f"
	}

#line 809 "dggesx.f"
    }

#line 811 "dggesx.f"
L60:

#line 813 "dggesx.f"
    work[1] = (doublereal) maxwrk;
#line 814 "dggesx.f"
    iwork[1] = liwmin;

#line 816 "dggesx.f"
    return 0;

/*     End of DGGESX */

} /* dggesx_ */

