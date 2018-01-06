#line 1 "zggevx.f"
/* zggevx.f -- translated by f2c (version 20100827).
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

#line 1 "zggevx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;

/* > \brief <b> ZGGEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, */
/*                          ALPHA, BETA, VL, LDVL, VR, LDVR, ILO, IHI, */
/*                          LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, */
/*                          WORK, LWORK, RWORK, IWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       DOUBLE PRECISION   ABNRM, BBNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   LSCALE( * ), RCONDE( * ), RCONDV( * ), */
/*      $                   RSCALE( * ), RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGEVX computes for a pair of N-by-N complex nonsymmetric matrices */
/* > (A,B) the generalized eigenvalues, and optionally, the left and/or */
/* > right generalized eigenvectors. */
/* > */
/* > Optionally, it also computes a balancing transformation to improve */
/* > the conditioning of the eigenvalues and eigenvectors (ILO, IHI, */
/* > LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for */
/* > the eigenvalues (RCONDE), and reciprocal condition numbers for the */
/* > right eigenvectors (RCONDV). */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar */
/* > lambda or a ratio alpha/beta = lambda, such that A - lambda*B is */
/* > singular. It is usually represented as the pair (alpha,beta), as */
/* > there is a reasonable interpretation for beta=0, and even for both */
/* > being zero. */
/* > */
/* > The right eigenvector v(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* >                  A * v(j) = lambda(j) * B * v(j) . */
/* > The left eigenvector u(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* >                  u(j)**H * A  = lambda(j) * u(j)**H * B. */
/* > where u(j)**H is the conjugate-transpose of u(j). */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] BALANC */
/* > \verbatim */
/* >          BALANC is CHARACTER*1 */
/* >          Specifies the balance option to be performed: */
/* >          = 'N':  do not diagonally scale or permute; */
/* >          = 'P':  permute only; */
/* >          = 'S':  scale only; */
/* >          = 'B':  both permute and scale. */
/* >          Computed reciprocal condition numbers will be for the */
/* >          matrices after permuting and/or balancing. Permuting does */
/* >          not change condition numbers (in exact arithmetic), but */
/* >          balancing does. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVL */
/* > \verbatim */
/* >          JOBVL is CHARACTER*1 */
/* >          = 'N':  do not compute the left generalized eigenvectors; */
/* >          = 'V':  compute the left generalized eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* >          JOBVR is CHARACTER*1 */
/* >          = 'N':  do not compute the right generalized eigenvectors; */
/* >          = 'V':  compute the right generalized eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] SENSE */
/* > \verbatim */
/* >          SENSE is CHARACTER*1 */
/* >          Determines which reciprocal condition numbers are computed. */
/* >          = 'N': none are computed; */
/* >          = 'E': computed for eigenvalues only; */
/* >          = 'V': computed for eigenvectors only; */
/* >          = 'B': computed for eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A, B, VL, and VR.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the matrix A in the pair (A,B). */
/* >          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V' */
/* >          or both, then A contains the first part of the complex Schur */
/* >          form of the "balanced" versions of the input A and B. */
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
/* >          On entry, the matrix B in the pair (A,B). */
/* >          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V' */
/* >          or both, then B contains the second part of the complex */
/* >          Schur form of the "balanced" versions of the input A and B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of B.  LDB >= max(1,N). */
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
/* >          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the generalized */
/* >          eigenvalues. */
/* > */
/* >          Note: the quotient ALPHA(j)/BETA(j) ) may easily over- or */
/* >          underflow, and BETA(j) may even be zero.  Thus, the user */
/* >          should avoid naively computing the ratio ALPHA/BETA. */
/* >          However, ALPHA will be always less than and usually */
/* >          comparable with norm(A) in magnitude, and BETA always less */
/* >          than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is COMPLEX*16 array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left generalized eigenvectors u(j) are */
/* >          stored one after another in the columns of VL, in the same */
/* >          order as their eigenvalues. */
/* >          Each eigenvector will be scaled so the largest component */
/* >          will have abs(real part) + abs(imag. part) = 1. */
/* >          Not referenced if JOBVL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the matrix VL. LDVL >= 1, and */
/* >          if JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* >          VR is COMPLEX*16 array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right generalized eigenvectors v(j) are */
/* >          stored one after another in the columns of VR, in the same */
/* >          order as their eigenvalues. */
/* >          Each eigenvector will be scaled so the largest component */
/* >          will have abs(real part) + abs(imag. part) = 1. */
/* >          Not referenced if JOBVR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the matrix VR. LDVR >= 1, and */
/* >          if JOBVR = 'V', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          ILO and IHI are integer values such that on exit */
/* >          A(i,j) = 0 and B(i,j) = 0 if i > j and */
/* >          j = 1,...,ILO-1 or i = IHI+1,...,N. */
/* >          If BALANC = 'N' or 'S', ILO = 1 and IHI = N. */
/* > \endverbatim */
/* > */
/* > \param[out] LSCALE */
/* > \verbatim */
/* >          LSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and scaling factors applied */
/* >          to the left side of A and B.  If PL(j) is the index of the */
/* >          row interchanged with row j, and DL(j) is the scaling */
/* >          factor applied to row j, then */
/* >            LSCALE(j) = PL(j)  for j = 1,...,ILO-1 */
/* >                      = DL(j)  for j = ILO,...,IHI */
/* >                      = PL(j)  for j = IHI+1,...,N. */
/* >          The order in which the interchanges are made is N to IHI+1, */
/* >          then 1 to ILO-1. */
/* > \endverbatim */
/* > */
/* > \param[out] RSCALE */
/* > \verbatim */
/* >          RSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and scaling factors applied */
/* >          to the right side of A and B.  If PR(j) is the index of the */
/* >          column interchanged with column j, and DR(j) is the scaling */
/* >          factor applied to column j, then */
/* >            RSCALE(j) = PR(j)  for j = 1,...,ILO-1 */
/* >                      = DR(j)  for j = ILO,...,IHI */
/* >                      = PR(j)  for j = IHI+1,...,N */
/* >          The order in which the interchanges are made is N to IHI+1, */
/* >          then 1 to ILO-1. */
/* > \endverbatim */
/* > */
/* > \param[out] ABNRM */
/* > \verbatim */
/* >          ABNRM is DOUBLE PRECISION */
/* >          The one-norm of the balanced matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] BBNRM */
/* > \verbatim */
/* >          BBNRM is DOUBLE PRECISION */
/* >          The one-norm of the balanced matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDE */
/* > \verbatim */
/* >          RCONDE is DOUBLE PRECISION array, dimension (N) */
/* >          If SENSE = 'E' or 'B', the reciprocal condition numbers of */
/* >          the eigenvalues, stored in consecutive elements of the array. */
/* >          If SENSE = 'N' or 'V', RCONDE is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is DOUBLE PRECISION array, dimension (N) */
/* >          If JOB = 'V' or 'B', the estimated reciprocal condition */
/* >          numbers of the eigenvectors, stored in consecutive elements */
/* >          of the array. If the eigenvalues cannot be reordered to */
/* >          compute RCONDV(j), RCONDV(j) is set to 0; this can only occur */
/* >          when the true value would be very small anyway. */
/* >          If SENSE = 'N' or 'E', RCONDV is not referenced. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,2*N). */
/* >          If SENSE = 'E', LWORK >= max(1,4*N). */
/* >          If SENSE = 'V' or 'B', LWORK >= max(1,2*N*N+2*N). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (lrwork) */
/* >          lrwork must be at least max(1,6*N) if BALANC = 'S' or 'B', */
/* >          and at least max(1,2*N) otherwise. */
/* >          Real workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N+2) */
/* >          If SENSE = 'E', IWORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] BWORK */
/* > \verbatim */
/* >          BWORK is LOGICAL array, dimension (N) */
/* >          If SENSE = 'N', BWORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          = 1,...,N: */
/* >                The QZ iteration failed.  No eigenvectors have been */
/* >                calculated, but ALPHA(j) and BETA(j) should be correct */
/* >                for j=INFO+1,...,N. */
/* >          > N:  =N+1: other than QZ iteration failed in ZHGEQZ. */
/* >                =N+2: error return from ZTGEVC. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup complex16GEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Balancing a matrix pair (A,B) includes, first, permuting rows and */
/* >  columns to isolate eigenvalues, second, applying diagonal similarity */
/* >  transformation to the rows and columns to make the rows and columns */
/* >  as close in norm as possible. The computed reciprocal condition */
/* >  numbers correspond to the balanced matrix. Permuting rows and columns */
/* >  will not change the condition numbers (in exact arithmetic) but */
/* >  diagonal scaling will.  For further explanation of balancing, see */
/* >  section 4.11.1.2 of LAPACK Users' Guide. */
/* > */
/* >  An approximate error bound on the chordal distance between the i-th */
/* >  computed generalized eigenvalue w and the corresponding exact */
/* >  eigenvalue lambda is */
/* > */
/* >       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I) */
/* > */
/* >  An approximate error bound for the angle between the i-th computed */
/* >  eigenvector VL(i) or VR(i) is given by */
/* > */
/* >       EPS * norm(ABNRM, BBNRM) / DIF(i). */
/* > */
/* >  For further explanation of the reciprocal condition numbers RCONDE */
/* >  and RCONDV, see section 4.11 of LAPACK User's Guide. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *alpha, doublecomplex *beta, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, logical *bwork, integer *info, ftnlen balanc_len, 
	ftnlen jobvl_len, ftnlen jobvr_len, ftnlen sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, m, jc, in, jr;
    static doublereal eps;
    static logical ilv;
    static doublereal anrm, bnrm;
    static integer ierr, itau;
    static doublereal temp;
    static logical ilvl, ilvr;
    static integer iwrk, iwrk1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer icols;
    static logical noscl;
    static integer irows;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), zggbak_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen), zggbal_(
	    char *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical ldumma[1];
    static char chtemp[1];
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static integer ijobvl;
    extern /* Subroutine */ int zgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer ijobvr;
    static logical wantsb;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static doublereal anrmto;
    static logical wantse;
    static doublereal bnrmto;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), ztgevc_(
	    char *, char *, logical *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     doublereal *, integer *, ftnlen, ftnlen), ztgsna_(char *, char *,
	     logical *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, doublecomplex *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int zhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static integer maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    static logical lquery, wantsv;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.4.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the input arguments */

#line 443 "zggevx.f"
    /* Parameter adjustments */
#line 443 "zggevx.f"
    a_dim1 = *lda;
#line 443 "zggevx.f"
    a_offset = 1 + a_dim1;
#line 443 "zggevx.f"
    a -= a_offset;
#line 443 "zggevx.f"
    b_dim1 = *ldb;
#line 443 "zggevx.f"
    b_offset = 1 + b_dim1;
#line 443 "zggevx.f"
    b -= b_offset;
#line 443 "zggevx.f"
    --alpha;
#line 443 "zggevx.f"
    --beta;
#line 443 "zggevx.f"
    vl_dim1 = *ldvl;
#line 443 "zggevx.f"
    vl_offset = 1 + vl_dim1;
#line 443 "zggevx.f"
    vl -= vl_offset;
#line 443 "zggevx.f"
    vr_dim1 = *ldvr;
#line 443 "zggevx.f"
    vr_offset = 1 + vr_dim1;
#line 443 "zggevx.f"
    vr -= vr_offset;
#line 443 "zggevx.f"
    --lscale;
#line 443 "zggevx.f"
    --rscale;
#line 443 "zggevx.f"
    --rconde;
#line 443 "zggevx.f"
    --rcondv;
#line 443 "zggevx.f"
    --work;
#line 443 "zggevx.f"
    --rwork;
#line 443 "zggevx.f"
    --iwork;
#line 443 "zggevx.f"
    --bwork;
#line 443 "zggevx.f"

#line 443 "zggevx.f"
    /* Function Body */
#line 443 "zggevx.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 444 "zggevx.f"
	ijobvl = 1;
#line 445 "zggevx.f"
	ilvl = FALSE_;
#line 446 "zggevx.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 447 "zggevx.f"
	ijobvl = 2;
#line 448 "zggevx.f"
	ilvl = TRUE_;
#line 449 "zggevx.f"
    } else {
#line 450 "zggevx.f"
	ijobvl = -1;
#line 451 "zggevx.f"
	ilvl = FALSE_;
#line 452 "zggevx.f"
    }

#line 454 "zggevx.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 455 "zggevx.f"
	ijobvr = 1;
#line 456 "zggevx.f"
	ilvr = FALSE_;
#line 457 "zggevx.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 458 "zggevx.f"
	ijobvr = 2;
#line 459 "zggevx.f"
	ilvr = TRUE_;
#line 460 "zggevx.f"
    } else {
#line 461 "zggevx.f"
	ijobvr = -1;
#line 462 "zggevx.f"
	ilvr = FALSE_;
#line 463 "zggevx.f"
    }
#line 464 "zggevx.f"
    ilv = ilvl || ilvr;

#line 466 "zggevx.f"
    noscl = lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (
	    ftnlen)1, (ftnlen)1);
#line 467 "zggevx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 468 "zggevx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 469 "zggevx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 470 "zggevx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 474 "zggevx.f"
    *info = 0;
#line 475 "zggevx.f"
    lquery = *lwork == -1;
#line 476 "zggevx.f"
    if (! (noscl || lsame_(balanc, "S", (ftnlen)1, (ftnlen)1) || lsame_(
	    balanc, "B", (ftnlen)1, (ftnlen)1))) {
#line 478 "zggevx.f"
	*info = -1;
#line 479 "zggevx.f"
    } else if (ijobvl <= 0) {
#line 480 "zggevx.f"
	*info = -2;
#line 481 "zggevx.f"
    } else if (ijobvr <= 0) {
#line 482 "zggevx.f"
	*info = -3;
#line 483 "zggevx.f"
    } else if (! (wantsn || wantse || wantsb || wantsv)) {
#line 485 "zggevx.f"
	*info = -4;
#line 486 "zggevx.f"
    } else if (*n < 0) {
#line 487 "zggevx.f"
	*info = -5;
#line 488 "zggevx.f"
    } else if (*lda < max(1,*n)) {
#line 489 "zggevx.f"
	*info = -7;
#line 490 "zggevx.f"
    } else if (*ldb < max(1,*n)) {
#line 491 "zggevx.f"
	*info = -9;
#line 492 "zggevx.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 493 "zggevx.f"
	*info = -13;
#line 494 "zggevx.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 495 "zggevx.f"
	*info = -15;
#line 496 "zggevx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. The workspace is */
/*       computed assuming ILO = 1 and IHI = N, the worst case.) */

#line 506 "zggevx.f"
    if (*info == 0) {
#line 507 "zggevx.f"
	if (*n == 0) {
#line 508 "zggevx.f"
	    minwrk = 1;
#line 509 "zggevx.f"
	    maxwrk = 1;
#line 510 "zggevx.f"
	} else {
#line 511 "zggevx.f"
	    minwrk = *n << 1;
#line 512 "zggevx.f"
	    if (wantse) {
#line 513 "zggevx.f"
		minwrk = *n << 2;
#line 514 "zggevx.f"
	    } else if (wantsv || wantsb) {
#line 515 "zggevx.f"
		minwrk = (*n << 1) * (*n + 1);
#line 516 "zggevx.f"
	    }
#line 517 "zggevx.f"
	    maxwrk = minwrk;
/* Computing MAX */
#line 518 "zggevx.f"
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 518 "zggevx.f"
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 520 "zggevx.f"
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "ZUNMQR", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 520 "zggevx.f"
	    maxwrk = max(i__1,i__2);
#line 522 "zggevx.f"
	    if (ilvl) {
/* Computing MAX */
#line 523 "zggevx.f"
		i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", 
			" ", n, &c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 523 "zggevx.f"
		maxwrk = max(i__1,i__2);
#line 525 "zggevx.f"
	    }
#line 526 "zggevx.f"
	}
#line 527 "zggevx.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 529 "zggevx.f"
	if (*lwork < minwrk && ! lquery) {
#line 530 "zggevx.f"
	    *info = -25;
#line 531 "zggevx.f"
	}
#line 532 "zggevx.f"
    }

#line 534 "zggevx.f"
    if (*info != 0) {
#line 535 "zggevx.f"
	i__1 = -(*info);
#line 535 "zggevx.f"
	xerbla_("ZGGEVX", &i__1, (ftnlen)6);
#line 536 "zggevx.f"
	return 0;
#line 537 "zggevx.f"
    } else if (lquery) {
#line 538 "zggevx.f"
	return 0;
#line 539 "zggevx.f"
    }

/*     Quick return if possible */

#line 543 "zggevx.f"
    if (*n == 0) {
#line 543 "zggevx.f"
	return 0;
#line 543 "zggevx.f"
    }

/*     Get machine constants */

#line 548 "zggevx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 549 "zggevx.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 550 "zggevx.f"
    bignum = 1. / smlnum;
#line 551 "zggevx.f"
    dlabad_(&smlnum, &bignum);
#line 552 "zggevx.f"
    smlnum = sqrt(smlnum) / eps;
#line 553 "zggevx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 557 "zggevx.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 558 "zggevx.f"
    ilascl = FALSE_;
#line 559 "zggevx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 560 "zggevx.f"
	anrmto = smlnum;
#line 561 "zggevx.f"
	ilascl = TRUE_;
#line 562 "zggevx.f"
    } else if (anrm > bignum) {
#line 563 "zggevx.f"
	anrmto = bignum;
#line 564 "zggevx.f"
	ilascl = TRUE_;
#line 565 "zggevx.f"
    }
#line 566 "zggevx.f"
    if (ilascl) {
#line 566 "zggevx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 566 "zggevx.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 571 "zggevx.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 572 "zggevx.f"
    ilbscl = FALSE_;
#line 573 "zggevx.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 574 "zggevx.f"
	bnrmto = smlnum;
#line 575 "zggevx.f"
	ilbscl = TRUE_;
#line 576 "zggevx.f"
    } else if (bnrm > bignum) {
#line 577 "zggevx.f"
	bnrmto = bignum;
#line 578 "zggevx.f"
	ilbscl = TRUE_;
#line 579 "zggevx.f"
    }
#line 580 "zggevx.f"
    if (ilbscl) {
#line 580 "zggevx.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 580 "zggevx.f"
    }

/*     Permute and/or balance the matrix pair (A,B) */
/*     (Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise) */

#line 586 "zggevx.f"
    zggbal_(balanc, n, &a[a_offset], lda, &b[b_offset], ldb, ilo, ihi, &
	    lscale[1], &rscale[1], &rwork[1], &ierr, (ftnlen)1);

/*     Compute ABNRM and BBNRM */

#line 591 "zggevx.f"
    *abnrm = zlange_("1", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 592 "zggevx.f"
    if (ilascl) {
#line 593 "zggevx.f"
	rwork[1] = *abnrm;
#line 594 "zggevx.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, &c__1, &c__1, &rwork[1], &
		c__1, &ierr, (ftnlen)1);
#line 596 "zggevx.f"
	*abnrm = rwork[1];
#line 597 "zggevx.f"
    }

#line 599 "zggevx.f"
    *bbnrm = zlange_("1", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 600 "zggevx.f"
    if (ilbscl) {
#line 601 "zggevx.f"
	rwork[1] = *bbnrm;
#line 602 "zggevx.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, &c__1, &c__1, &rwork[1], &
		c__1, &ierr, (ftnlen)1);
#line 604 "zggevx.f"
	*bbnrm = rwork[1];
#line 605 "zggevx.f"
    }

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Complex Workspace: need N, prefer N*NB ) */

#line 610 "zggevx.f"
    irows = *ihi + 1 - *ilo;
#line 611 "zggevx.f"
    if (ilv || ! wantsn) {
#line 612 "zggevx.f"
	icols = *n + 1 - *ilo;
#line 613 "zggevx.f"
    } else {
#line 614 "zggevx.f"
	icols = irows;
#line 615 "zggevx.f"
    }
#line 616 "zggevx.f"
    itau = 1;
#line 617 "zggevx.f"
    iwrk = itau + irows;
#line 618 "zggevx.f"
    i__1 = *lwork + 1 - iwrk;
#line 618 "zggevx.f"
    zgeqrf_(&irows, &icols, &b[*ilo + *ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the unitary transformation to A */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 624 "zggevx.f"
    i__1 = *lwork + 1 - iwrk;
#line 624 "zggevx.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[*ilo + *ilo * b_dim1], ldb, &
	    work[itau], &a[*ilo + *ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL and/or VR */
/*     (Workspace: need N, prefer N*NB) */

#line 631 "zggevx.f"
    if (ilvl) {
#line 632 "zggevx.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vl[vl_offset], ldvl, (ftnlen)4);
#line 633 "zggevx.f"
	if (irows > 1) {
#line 634 "zggevx.f"
	    i__1 = irows - 1;
#line 634 "zggevx.f"
	    i__2 = irows - 1;
#line 634 "zggevx.f"
	    zlacpy_("L", &i__1, &i__2, &b[*ilo + 1 + *ilo * b_dim1], ldb, &vl[
		    *ilo + 1 + *ilo * vl_dim1], ldvl, (ftnlen)1);
#line 636 "zggevx.f"
	}
#line 637 "zggevx.f"
	i__1 = *lwork + 1 - iwrk;
#line 637 "zggevx.f"
	zungqr_(&irows, &irows, &irows, &vl[*ilo + *ilo * vl_dim1], ldvl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 639 "zggevx.f"
    }

#line 641 "zggevx.f"
    if (ilvr) {
#line 641 "zggevx.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vr[vr_offset], ldvr, (ftnlen)4);
#line 641 "zggevx.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 647 "zggevx.f"
    if (ilv || ! wantsn) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 651 "zggevx.f"
	zgghrd_(jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr, (
		ftnlen)1, (ftnlen)1);
#line 653 "zggevx.f"
    } else {
#line 654 "zggevx.f"
	zgghrd_("N", "N", &irows, &c__1, &irows, &a[*ilo + *ilo * a_dim1], 
		lda, &b[*ilo + *ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 656 "zggevx.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur forms and Schur vectors) */
/*     (Complex Workspace: need N) */
/*     (Real Workspace: need N) */

#line 663 "zggevx.f"
    iwrk = itau;
#line 664 "zggevx.f"
    if (ilv || ! wantsn) {
#line 665 "zggevx.f"
	*(unsigned char *)chtemp = 'S';
#line 666 "zggevx.f"
    } else {
#line 667 "zggevx.f"
	*(unsigned char *)chtemp = 'E';
#line 668 "zggevx.f"
    }

#line 670 "zggevx.f"
    i__1 = *lwork + 1 - iwrk;
#line 670 "zggevx.f"
    zhgeqz_(chtemp, jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[b_offset]
	    , ldb, &alpha[1], &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], 
	    ldvr, &work[iwrk], &i__1, &rwork[1], &ierr, (ftnlen)1, (ftnlen)1, 
	    (ftnlen)1);
#line 673 "zggevx.f"
    if (ierr != 0) {
#line 674 "zggevx.f"
	if (ierr > 0 && ierr <= *n) {
#line 675 "zggevx.f"
	    *info = ierr;
#line 676 "zggevx.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 677 "zggevx.f"
	    *info = ierr - *n;
#line 678 "zggevx.f"
	} else {
#line 679 "zggevx.f"
	    *info = *n + 1;
#line 680 "zggevx.f"
	}
#line 681 "zggevx.f"
	goto L90;
#line 682 "zggevx.f"
    }

/*     Compute Eigenvectors and estimate condition numbers if desired */
/*     ZTGEVC: (Complex Workspace: need 2*N ) */
/*             (Real Workspace:    need 2*N ) */
/*     ZTGSNA: (Complex Workspace: need 2*N*N if SENSE='V' or 'B') */
/*             (Integer Workspace: need N+2 ) */

#line 690 "zggevx.f"
    if (ilv || ! wantsn) {
#line 691 "zggevx.f"
	if (ilv) {
#line 692 "zggevx.f"
	    if (ilvl) {
#line 693 "zggevx.f"
		if (ilvr) {
#line 694 "zggevx.f"
		    *(unsigned char *)chtemp = 'B';
#line 695 "zggevx.f"
		} else {
#line 696 "zggevx.f"
		    *(unsigned char *)chtemp = 'L';
#line 697 "zggevx.f"
		}
#line 698 "zggevx.f"
	    } else {
#line 699 "zggevx.f"
		*(unsigned char *)chtemp = 'R';
#line 700 "zggevx.f"
	    }

#line 702 "zggevx.f"
	    ztgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], 
		    ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &
		    work[iwrk], &rwork[1], &ierr, (ftnlen)1, (ftnlen)1);
#line 705 "zggevx.f"
	    if (ierr != 0) {
#line 706 "zggevx.f"
		*info = *n + 2;
#line 707 "zggevx.f"
		goto L90;
#line 708 "zggevx.f"
	    }
#line 709 "zggevx.f"
	}

#line 711 "zggevx.f"
	if (! wantsn) {

/*           compute eigenvectors (DTGEVC) and estimate condition */
/*           numbers (DTGSNA). Note that the definition of the condition */
/*           number is not invariant under transformation (u,v) to */
/*           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized */
/*           Schur form (S,T), Q and Z are orthogonal matrices. In order */
/*           to avoid using extra 2*N*N workspace, we have to */
/*           re-calculate eigenvectors and estimate the condition numbers */
/*           one at a time. */

#line 722 "zggevx.f"
	    i__1 = *n;
#line 722 "zggevx.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

#line 724 "zggevx.f"
		i__2 = *n;
#line 724 "zggevx.f"
		for (j = 1; j <= i__2; ++j) {
#line 725 "zggevx.f"
		    bwork[j] = FALSE_;
#line 726 "zggevx.f"
/* L10: */
#line 726 "zggevx.f"
		}
#line 727 "zggevx.f"
		bwork[i__] = TRUE_;

#line 729 "zggevx.f"
		iwrk = *n + 1;
#line 730 "zggevx.f"
		iwrk1 = iwrk + *n;

#line 732 "zggevx.f"
		if (wantse || wantsb) {
#line 733 "zggevx.f"
		    ztgevc_("B", "S", &bwork[1], n, &a[a_offset], lda, &b[
			    b_offset], ldb, &work[1], n, &work[iwrk], n, &
			    c__1, &m, &work[iwrk1], &rwork[1], &ierr, (ftnlen)
			    1, (ftnlen)1);
#line 736 "zggevx.f"
		    if (ierr != 0) {
#line 737 "zggevx.f"
			*info = *n + 2;
#line 738 "zggevx.f"
			goto L90;
#line 739 "zggevx.f"
		    }
#line 740 "zggevx.f"
		}

#line 742 "zggevx.f"
		i__2 = *lwork - iwrk1 + 1;
#line 742 "zggevx.f"
		ztgsna_(sense, "S", &bwork[1], n, &a[a_offset], lda, &b[
			b_offset], ldb, &work[1], n, &work[iwrk], n, &rconde[
			i__], &rcondv[i__], &c__1, &m, &work[iwrk1], &i__2, &
			iwork[1], &ierr, (ftnlen)1, (ftnlen)1);

#line 747 "zggevx.f"
/* L20: */
#line 747 "zggevx.f"
	    }
#line 748 "zggevx.f"
	}
#line 749 "zggevx.f"
    }

/*     Undo balancing on VL and VR and normalization */
/*     (Workspace: none needed) */

#line 754 "zggevx.f"
    if (ilvl) {
#line 755 "zggevx.f"
	zggbak_(balanc, "L", n, ilo, ihi, &lscale[1], &rscale[1], n, &vl[
		vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);

#line 758 "zggevx.f"
	i__1 = *n;
#line 758 "zggevx.f"
	for (jc = 1; jc <= i__1; ++jc) {
#line 759 "zggevx.f"
	    temp = 0.;
#line 760 "zggevx.f"
	    i__2 = *n;
#line 760 "zggevx.f"
	    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 761 "zggevx.f"
		i__3 = jr + jc * vl_dim1;
#line 761 "zggevx.f"
		d__3 = temp, d__4 = (d__1 = vl[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&vl[jr + jc * vl_dim1]), abs(d__2));
#line 761 "zggevx.f"
		temp = max(d__3,d__4);
#line 762 "zggevx.f"
/* L30: */
#line 762 "zggevx.f"
	    }
#line 763 "zggevx.f"
	    if (temp < smlnum) {
#line 763 "zggevx.f"
		goto L50;
#line 763 "zggevx.f"
	    }
#line 765 "zggevx.f"
	    temp = 1. / temp;
#line 766 "zggevx.f"
	    i__2 = *n;
#line 766 "zggevx.f"
	    for (jr = 1; jr <= i__2; ++jr) {
#line 767 "zggevx.f"
		i__3 = jr + jc * vl_dim1;
#line 767 "zggevx.f"
		i__4 = jr + jc * vl_dim1;
#line 767 "zggevx.f"
		z__1.r = temp * vl[i__4].r, z__1.i = temp * vl[i__4].i;
#line 767 "zggevx.f"
		vl[i__3].r = z__1.r, vl[i__3].i = z__1.i;
#line 768 "zggevx.f"
/* L40: */
#line 768 "zggevx.f"
	    }
#line 769 "zggevx.f"
L50:
#line 769 "zggevx.f"
	    ;
#line 769 "zggevx.f"
	}
#line 770 "zggevx.f"
    }

#line 772 "zggevx.f"
    if (ilvr) {
#line 773 "zggevx.f"
	zggbak_(balanc, "R", n, ilo, ihi, &lscale[1], &rscale[1], n, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 775 "zggevx.f"
	i__1 = *n;
#line 775 "zggevx.f"
	for (jc = 1; jc <= i__1; ++jc) {
#line 776 "zggevx.f"
	    temp = 0.;
#line 777 "zggevx.f"
	    i__2 = *n;
#line 777 "zggevx.f"
	    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 778 "zggevx.f"
		i__3 = jr + jc * vr_dim1;
#line 778 "zggevx.f"
		d__3 = temp, d__4 = (d__1 = vr[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&vr[jr + jc * vr_dim1]), abs(d__2));
#line 778 "zggevx.f"
		temp = max(d__3,d__4);
#line 779 "zggevx.f"
/* L60: */
#line 779 "zggevx.f"
	    }
#line 780 "zggevx.f"
	    if (temp < smlnum) {
#line 780 "zggevx.f"
		goto L80;
#line 780 "zggevx.f"
	    }
#line 782 "zggevx.f"
	    temp = 1. / temp;
#line 783 "zggevx.f"
	    i__2 = *n;
#line 783 "zggevx.f"
	    for (jr = 1; jr <= i__2; ++jr) {
#line 784 "zggevx.f"
		i__3 = jr + jc * vr_dim1;
#line 784 "zggevx.f"
		i__4 = jr + jc * vr_dim1;
#line 784 "zggevx.f"
		z__1.r = temp * vr[i__4].r, z__1.i = temp * vr[i__4].i;
#line 784 "zggevx.f"
		vr[i__3].r = z__1.r, vr[i__3].i = z__1.i;
#line 785 "zggevx.f"
/* L70: */
#line 785 "zggevx.f"
	    }
#line 786 "zggevx.f"
L80:
#line 786 "zggevx.f"
	    ;
#line 786 "zggevx.f"
	}
#line 787 "zggevx.f"
    }

/*     Undo scaling if necessary */

#line 791 "zggevx.f"
L90:

#line 793 "zggevx.f"
    if (ilascl) {
#line 793 "zggevx.f"
	zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 793 "zggevx.f"
    }

#line 796 "zggevx.f"
    if (ilbscl) {
#line 796 "zggevx.f"
	zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 796 "zggevx.f"
    }

#line 799 "zggevx.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 800 "zggevx.f"
    return 0;

/*     End of ZGGEVX */

} /* zggevx_ */

