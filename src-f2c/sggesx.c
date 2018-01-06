#line 1 "sggesx.f"
/* sggesx.f -- translated by f2c (version 20100827).
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

#line 1 "sggesx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
static doublereal c_b42 = 0.;
static doublereal c_b43 = 1.;

/* > \brief <b> SGGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGESX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggesx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggesx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggesx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, */
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
/*       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
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
/* > SGGESX computes for a pair of N-by-N real nonsymmetric matrices */
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
/* >          SELCTG is procedure) LOGICAL FUNCTION of three REAL arguments */
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
/* >          reciprocal condition numbers for the selected deflating */
/* >          subspaces. */
/* >          Not referenced if SENSE = 'N' or 'E'. */
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
/* >          > N:  =N+1: other than QZ iteration failed in SHGEQZ */
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
/* Subroutine */ int sggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
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
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static logical wantsb, wantse, lastsl;
    static integer liwmin;
    static doublereal anrmto, bnrmto;
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    extern /* Subroutine */ int shgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), slaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), sorgqr_(integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    stgsen_(integer *, logical *, logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    static logical wantst, lquery, wantsv;
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

#line 428 "sggesx.f"
    /* Parameter adjustments */
#line 428 "sggesx.f"
    a_dim1 = *lda;
#line 428 "sggesx.f"
    a_offset = 1 + a_dim1;
#line 428 "sggesx.f"
    a -= a_offset;
#line 428 "sggesx.f"
    b_dim1 = *ldb;
#line 428 "sggesx.f"
    b_offset = 1 + b_dim1;
#line 428 "sggesx.f"
    b -= b_offset;
#line 428 "sggesx.f"
    --alphar;
#line 428 "sggesx.f"
    --alphai;
#line 428 "sggesx.f"
    --beta;
#line 428 "sggesx.f"
    vsl_dim1 = *ldvsl;
#line 428 "sggesx.f"
    vsl_offset = 1 + vsl_dim1;
#line 428 "sggesx.f"
    vsl -= vsl_offset;
#line 428 "sggesx.f"
    vsr_dim1 = *ldvsr;
#line 428 "sggesx.f"
    vsr_offset = 1 + vsr_dim1;
#line 428 "sggesx.f"
    vsr -= vsr_offset;
#line 428 "sggesx.f"
    --rconde;
#line 428 "sggesx.f"
    --rcondv;
#line 428 "sggesx.f"
    --work;
#line 428 "sggesx.f"
    --iwork;
#line 428 "sggesx.f"
    --bwork;
#line 428 "sggesx.f"

#line 428 "sggesx.f"
    /* Function Body */
#line 428 "sggesx.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 429 "sggesx.f"
	ijobvl = 1;
#line 430 "sggesx.f"
	ilvsl = FALSE_;
#line 431 "sggesx.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 432 "sggesx.f"
	ijobvl = 2;
#line 433 "sggesx.f"
	ilvsl = TRUE_;
#line 434 "sggesx.f"
    } else {
#line 435 "sggesx.f"
	ijobvl = -1;
#line 436 "sggesx.f"
	ilvsl = FALSE_;
#line 437 "sggesx.f"
    }

#line 439 "sggesx.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 440 "sggesx.f"
	ijobvr = 1;
#line 441 "sggesx.f"
	ilvsr = FALSE_;
#line 442 "sggesx.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 443 "sggesx.f"
	ijobvr = 2;
#line 444 "sggesx.f"
	ilvsr = TRUE_;
#line 445 "sggesx.f"
    } else {
#line 446 "sggesx.f"
	ijobvr = -1;
#line 447 "sggesx.f"
	ilvsr = FALSE_;
#line 448 "sggesx.f"
    }

#line 450 "sggesx.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 451 "sggesx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 452 "sggesx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 453 "sggesx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 454 "sggesx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 455 "sggesx.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 456 "sggesx.f"
    if (wantsn) {
#line 457 "sggesx.f"
	ijob = 0;
#line 458 "sggesx.f"
    } else if (wantse) {
#line 459 "sggesx.f"
	ijob = 1;
#line 460 "sggesx.f"
    } else if (wantsv) {
#line 461 "sggesx.f"
	ijob = 2;
#line 462 "sggesx.f"
    } else if (wantsb) {
#line 463 "sggesx.f"
	ijob = 4;
#line 464 "sggesx.f"
    }

/*     Test the input arguments */

#line 468 "sggesx.f"
    *info = 0;
#line 469 "sggesx.f"
    if (ijobvl <= 0) {
#line 470 "sggesx.f"
	*info = -1;
#line 471 "sggesx.f"
    } else if (ijobvr <= 0) {
#line 472 "sggesx.f"
	*info = -2;
#line 473 "sggesx.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 474 "sggesx.f"
	*info = -3;
#line 475 "sggesx.f"
    } else if (! (wantsn || wantse || wantsv || wantsb) || ! wantst && ! 
	    wantsn) {
#line 477 "sggesx.f"
	*info = -5;
#line 478 "sggesx.f"
    } else if (*n < 0) {
#line 479 "sggesx.f"
	*info = -6;
#line 480 "sggesx.f"
    } else if (*lda < max(1,*n)) {
#line 481 "sggesx.f"
	*info = -8;
#line 482 "sggesx.f"
    } else if (*ldb < max(1,*n)) {
#line 483 "sggesx.f"
	*info = -10;
#line 484 "sggesx.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 485 "sggesx.f"
	*info = -16;
#line 486 "sggesx.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 487 "sggesx.f"
	*info = -18;
#line 488 "sggesx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 497 "sggesx.f"
    if (*info == 0) {
#line 498 "sggesx.f"
	if (*n > 0) {
/* Computing MAX */
#line 499 "sggesx.f"
	    i__1 = *n << 3, i__2 = *n * 6 + 16;
#line 499 "sggesx.f"
	    minwrk = max(i__1,i__2);
#line 500 "sggesx.f"
	    maxwrk = minwrk - *n + *n * ilaenv_(&c__1, "SGEQRF", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 502 "sggesx.f"
	    i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "SORMQR", 
		    " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 502 "sggesx.f"
	    maxwrk = max(i__1,i__2);
#line 504 "sggesx.f"
	    if (ilvsl) {
/* Computing MAX */
#line 505 "sggesx.f"
		i__1 = maxwrk, i__2 = minwrk - *n + *n * ilaenv_(&c__1, "SOR"\
			"GQR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 505 "sggesx.f"
		maxwrk = max(i__1,i__2);
#line 507 "sggesx.f"
	    }
#line 508 "sggesx.f"
	    lwrk = maxwrk;
#line 509 "sggesx.f"
	    if (ijob >= 1) {
/* Computing MAX */
#line 509 "sggesx.f"
		i__1 = lwrk, i__2 = *n * *n / 2;
#line 509 "sggesx.f"
		lwrk = max(i__1,i__2);
#line 509 "sggesx.f"
	    }
#line 511 "sggesx.f"
	} else {
#line 512 "sggesx.f"
	    minwrk = 1;
#line 513 "sggesx.f"
	    maxwrk = 1;
#line 514 "sggesx.f"
	    lwrk = 1;
#line 515 "sggesx.f"
	}
#line 516 "sggesx.f"
	work[1] = (doublereal) lwrk;
#line 517 "sggesx.f"
	if (wantsn || *n == 0) {
#line 518 "sggesx.f"
	    liwmin = 1;
#line 519 "sggesx.f"
	} else {
#line 520 "sggesx.f"
	    liwmin = *n + 6;
#line 521 "sggesx.f"
	}
#line 522 "sggesx.f"
	iwork[1] = liwmin;

#line 524 "sggesx.f"
	if (*lwork < minwrk && ! lquery) {
#line 525 "sggesx.f"
	    *info = -22;
#line 526 "sggesx.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 527 "sggesx.f"
	    *info = -24;
#line 528 "sggesx.f"
	}
#line 529 "sggesx.f"
    }

#line 531 "sggesx.f"
    if (*info != 0) {
#line 532 "sggesx.f"
	i__1 = -(*info);
#line 532 "sggesx.f"
	xerbla_("SGGESX", &i__1, (ftnlen)6);
#line 533 "sggesx.f"
	return 0;
#line 534 "sggesx.f"
    } else if (lquery) {
#line 535 "sggesx.f"
	return 0;
#line 536 "sggesx.f"
    }

/*     Quick return if possible */

#line 540 "sggesx.f"
    if (*n == 0) {
#line 541 "sggesx.f"
	*sdim = 0;
#line 542 "sggesx.f"
	return 0;
#line 543 "sggesx.f"
    }

/*     Get machine constants */

#line 547 "sggesx.f"
    eps = slamch_("P", (ftnlen)1);
#line 548 "sggesx.f"
    safmin = slamch_("S", (ftnlen)1);
#line 549 "sggesx.f"
    safmax = 1. / safmin;
#line 550 "sggesx.f"
    slabad_(&safmin, &safmax);
#line 551 "sggesx.f"
    smlnum = sqrt(safmin) / eps;
#line 552 "sggesx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 556 "sggesx.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 557 "sggesx.f"
    ilascl = FALSE_;
#line 558 "sggesx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 559 "sggesx.f"
	anrmto = smlnum;
#line 560 "sggesx.f"
	ilascl = TRUE_;
#line 561 "sggesx.f"
    } else if (anrm > bignum) {
#line 562 "sggesx.f"
	anrmto = bignum;
#line 563 "sggesx.f"
	ilascl = TRUE_;
#line 564 "sggesx.f"
    }
#line 565 "sggesx.f"
    if (ilascl) {
#line 565 "sggesx.f"
	slascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 565 "sggesx.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 570 "sggesx.f"
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 571 "sggesx.f"
    ilbscl = FALSE_;
#line 572 "sggesx.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 573 "sggesx.f"
	bnrmto = smlnum;
#line 574 "sggesx.f"
	ilbscl = TRUE_;
#line 575 "sggesx.f"
    } else if (bnrm > bignum) {
#line 576 "sggesx.f"
	bnrmto = bignum;
#line 577 "sggesx.f"
	ilbscl = TRUE_;
#line 578 "sggesx.f"
    }
#line 579 "sggesx.f"
    if (ilbscl) {
#line 579 "sggesx.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 579 "sggesx.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Workspace: need 6*N + 2*N for permutation parameters) */

#line 585 "sggesx.f"
    ileft = 1;
#line 586 "sggesx.f"
    iright = *n + 1;
#line 587 "sggesx.f"
    iwrk = iright + *n;
#line 588 "sggesx.f"
    sggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Workspace: need N, prefer N*NB) */

#line 594 "sggesx.f"
    irows = ihi + 1 - ilo;
#line 595 "sggesx.f"
    icols = *n + 1 - ilo;
#line 596 "sggesx.f"
    itau = iwrk;
#line 597 "sggesx.f"
    iwrk = itau + irows;
#line 598 "sggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 598 "sggesx.f"
    sgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Workspace: need N, prefer N*NB) */

#line 604 "sggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 604 "sggesx.f"
    sormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VSL */
/*     (Workspace: need N, prefer N*NB) */

#line 611 "sggesx.f"
    if (ilvsl) {
#line 612 "sggesx.f"
	slaset_("Full", n, n, &c_b42, &c_b43, &vsl[vsl_offset], ldvsl, (
		ftnlen)4);
#line 613 "sggesx.f"
	if (irows > 1) {
#line 614 "sggesx.f"
	    i__1 = irows - 1;
#line 614 "sggesx.f"
	    i__2 = irows - 1;
#line 614 "sggesx.f"
	    slacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 616 "sggesx.f"
	}
#line 617 "sggesx.f"
	i__1 = *lwork + 1 - iwrk;
#line 617 "sggesx.f"
	sorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 619 "sggesx.f"
    }

/*     Initialize VSR */

#line 623 "sggesx.f"
    if (ilvsr) {
#line 623 "sggesx.f"
	slaset_("Full", n, n, &c_b42, &c_b43, &vsr[vsr_offset], ldvsr, (
		ftnlen)4);
#line 623 "sggesx.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 629 "sggesx.f"
    sgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 632 "sggesx.f"
    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Workspace: need N) */

#line 637 "sggesx.f"
    iwrk = itau;
#line 638 "sggesx.f"
    i__1 = *lwork + 1 - iwrk;
#line 638 "sggesx.f"
    shgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
	    , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 641 "sggesx.f"
    if (ierr != 0) {
#line 642 "sggesx.f"
	if (ierr > 0 && ierr <= *n) {
#line 643 "sggesx.f"
	    *info = ierr;
#line 644 "sggesx.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 645 "sggesx.f"
	    *info = ierr - *n;
#line 646 "sggesx.f"
	} else {
#line 647 "sggesx.f"
	    *info = *n + 1;
#line 648 "sggesx.f"
	}
#line 649 "sggesx.f"
	goto L50;
#line 650 "sggesx.f"
    }

/*     Sort eigenvalues ALPHA/BETA and compute the reciprocal of */
/*     condition number(s) */
/*     (Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) ) */
/*                 otherwise, need 8*(N+1) ) */

#line 657 "sggesx.f"
    if (wantst) {

/*        Undo scaling on eigenvalues before SELCTGing */

#line 661 "sggesx.f"
	if (ilascl) {
#line 662 "sggesx.f"
	    slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], 
		    n, &ierr, (ftnlen)1);
#line 664 "sggesx.f"
	    slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], 
		    n, &ierr, (ftnlen)1);
#line 666 "sggesx.f"
	}
#line 667 "sggesx.f"
	if (ilbscl) {
#line 667 "sggesx.f"
	    slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, 
		    &ierr, (ftnlen)1);
#line 667 "sggesx.f"
	}

/*        Select eigenvalues */

#line 672 "sggesx.f"
	i__1 = *n;
#line 672 "sggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 673 "sggesx.f"
	    bwork[i__] = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 674 "sggesx.f"
/* L10: */
#line 674 "sggesx.f"
	}

/*        Reorder eigenvalues, transform Generalized Schur vectors, and */
/*        compute reciprocal condition numbers */

#line 679 "sggesx.f"
	i__1 = *lwork - iwrk + 1;
#line 679 "sggesx.f"
	stgsen_(&ijob, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[
		vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pl, &pr, 
		dif, &work[iwrk], &i__1, &iwork[1], liwork, &ierr);

#line 684 "sggesx.f"
	if (ijob >= 1) {
/* Computing MAX */
#line 684 "sggesx.f"
	    i__1 = maxwrk, i__2 = (*sdim << 1) * (*n - *sdim);
#line 684 "sggesx.f"
	    maxwrk = max(i__1,i__2);
#line 684 "sggesx.f"
	}
#line 686 "sggesx.f"
	if (ierr == -22) {

/*            not enough real workspace */

#line 690 "sggesx.f"
	    *info = -22;
#line 691 "sggesx.f"
	} else {
#line 692 "sggesx.f"
	    if (ijob == 1 || ijob == 4) {
#line 693 "sggesx.f"
		rconde[1] = pl;
#line 694 "sggesx.f"
		rconde[2] = pr;
#line 695 "sggesx.f"
	    }
#line 696 "sggesx.f"
	    if (ijob == 2 || ijob == 4) {
#line 697 "sggesx.f"
		rcondv[1] = dif[0];
#line 698 "sggesx.f"
		rcondv[2] = dif[1];
#line 699 "sggesx.f"
	    }
#line 700 "sggesx.f"
	    if (ierr == 1) {
#line 700 "sggesx.f"
		*info = *n + 3;
#line 700 "sggesx.f"
	    }
#line 702 "sggesx.f"
	}

#line 704 "sggesx.f"
    }

/*     Apply permutation to VSL and VSR */
/*     (Workspace: none needed) */

#line 709 "sggesx.f"
    if (ilvsl) {
#line 709 "sggesx.f"
	sggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &ierr, (ftnlen)1, (ftnlen)1);
#line 709 "sggesx.f"
    }

#line 713 "sggesx.f"
    if (ilvsr) {
#line 713 "sggesx.f"
	sggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &ierr, (ftnlen)1, (ftnlen)1);
#line 713 "sggesx.f"
    }

/*     Check if unscaling would cause over/underflow, if so, rescale */
/*     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of */
/*     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I) */

#line 721 "sggesx.f"
    if (ilascl) {
#line 722 "sggesx.f"
	i__1 = *n;
#line 722 "sggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 723 "sggesx.f"
	    if (alphai[i__] != 0.) {
#line 724 "sggesx.f"
		if (alphar[i__] / safmax > anrmto / anrm || safmin / alphar[
			i__] > anrm / anrmto) {
#line 727 "sggesx.f"
		    work[1] = (d__1 = a[i__ + i__ * a_dim1] / alphar[i__], 
			    abs(d__1));
#line 728 "sggesx.f"
		    beta[i__] *= work[1];
#line 729 "sggesx.f"
		    alphar[i__] *= work[1];
#line 730 "sggesx.f"
		    alphai[i__] *= work[1];
#line 731 "sggesx.f"
		} else if (alphai[i__] / safmax > anrmto / anrm || safmin / 
			alphai[i__] > anrm / anrmto) {
#line 734 "sggesx.f"
		    work[1] = (d__1 = a[i__ + (i__ + 1) * a_dim1] / alphai[
			    i__], abs(d__1));
#line 735 "sggesx.f"
		    beta[i__] *= work[1];
#line 736 "sggesx.f"
		    alphar[i__] *= work[1];
#line 737 "sggesx.f"
		    alphai[i__] *= work[1];
#line 738 "sggesx.f"
		}
#line 739 "sggesx.f"
	    }
#line 740 "sggesx.f"
/* L20: */
#line 740 "sggesx.f"
	}
#line 741 "sggesx.f"
    }

#line 743 "sggesx.f"
    if (ilbscl) {
#line 744 "sggesx.f"
	i__1 = *n;
#line 744 "sggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 745 "sggesx.f"
	    if (alphai[i__] != 0.) {
#line 746 "sggesx.f"
		if (beta[i__] / safmax > bnrmto / bnrm || safmin / beta[i__] 
			> bnrm / bnrmto) {
#line 748 "sggesx.f"
		    work[1] = (d__1 = b[i__ + i__ * b_dim1] / beta[i__], abs(
			    d__1));
#line 749 "sggesx.f"
		    beta[i__] *= work[1];
#line 750 "sggesx.f"
		    alphar[i__] *= work[1];
#line 751 "sggesx.f"
		    alphai[i__] *= work[1];
#line 752 "sggesx.f"
		}
#line 753 "sggesx.f"
	    }
#line 754 "sggesx.f"
/* L25: */
#line 754 "sggesx.f"
	}
#line 755 "sggesx.f"
    }

/*     Undo scaling */

#line 759 "sggesx.f"
    if (ilascl) {
#line 760 "sggesx.f"
	slascl_("H", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 761 "sggesx.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 762 "sggesx.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 763 "sggesx.f"
    }

#line 765 "sggesx.f"
    if (ilbscl) {
#line 766 "sggesx.f"
	slascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 767 "sggesx.f"
	slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 768 "sggesx.f"
    }

#line 770 "sggesx.f"
    if (wantst) {

/*        Check if reordering is correct */

#line 774 "sggesx.f"
	lastsl = TRUE_;
#line 775 "sggesx.f"
	lst2sl = TRUE_;
#line 776 "sggesx.f"
	*sdim = 0;
#line 777 "sggesx.f"
	ip = 0;
#line 778 "sggesx.f"
	i__1 = *n;
#line 778 "sggesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 779 "sggesx.f"
	    cursl = (*selctg)(&alphar[i__], &alphai[i__], &beta[i__]);
#line 780 "sggesx.f"
	    if (alphai[i__] == 0.) {
#line 781 "sggesx.f"
		if (cursl) {
#line 781 "sggesx.f"
		    ++(*sdim);
#line 781 "sggesx.f"
		}
#line 783 "sggesx.f"
		ip = 0;
#line 784 "sggesx.f"
		if (cursl && ! lastsl) {
#line 784 "sggesx.f"
		    *info = *n + 2;
#line 784 "sggesx.f"
		}
#line 786 "sggesx.f"
	    } else {
#line 787 "sggesx.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 791 "sggesx.f"
		    cursl = cursl || lastsl;
#line 792 "sggesx.f"
		    lastsl = cursl;
#line 793 "sggesx.f"
		    if (cursl) {
#line 793 "sggesx.f"
			*sdim += 2;
#line 793 "sggesx.f"
		    }
#line 795 "sggesx.f"
		    ip = -1;
#line 796 "sggesx.f"
		    if (cursl && ! lst2sl) {
#line 796 "sggesx.f"
			*info = *n + 2;
#line 796 "sggesx.f"
		    }
#line 798 "sggesx.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 802 "sggesx.f"
		    ip = 1;
#line 803 "sggesx.f"
		}
#line 804 "sggesx.f"
	    }
#line 805 "sggesx.f"
	    lst2sl = lastsl;
#line 806 "sggesx.f"
	    lastsl = cursl;
#line 807 "sggesx.f"
/* L40: */
#line 807 "sggesx.f"
	}

#line 809 "sggesx.f"
    }

#line 811 "sggesx.f"
L50:

#line 813 "sggesx.f"
    work[1] = (doublereal) maxwrk;
#line 814 "sggesx.f"
    iwork[1] = liwmin;

#line 816 "sggesx.f"
    return 0;

/*     End of SGGESX */

} /* sggesx_ */

