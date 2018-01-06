#line 1 "dggevx.f"
/* dggevx.f -- translated by f2c (version 20100827).
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

#line 1 "dggevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b59 = 0.;
static doublereal c_b60 = 1.;

/* > \brief <b> DGGEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, */
/*                          ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, */
/*                          IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, */
/*                          RCONDV, WORK, LWORK, IWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       DOUBLE PRECISION   ABNRM, BBNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), LSCALE( * ), */
/*      $                   RCONDE( * ), RCONDV( * ), RSCALE( * ), */
/*      $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B) */
/* > the generalized eigenvalues, and optionally, the left and/or right */
/* > generalized eigenvectors. */
/* > */
/* > Optionally also, it computes a balancing transformation to improve */
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
/* > */
/* >                  A * v(j) = lambda(j) * B * v(j) . */
/* > */
/* > The left eigenvector u(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* > */
/* >                  u(j)**H * A  = lambda(j) * u(j)**H * B. */
/* > */
/* > where u(j)**H is the conjugate-transpose of u(j). */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] BALANC */
/* > \verbatim */
/* >          BALANC is CHARACTER*1 */
/* >          Specifies the balance option to be performed. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the matrix A in the pair (A,B). */
/* >          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V' */
/* >          or both, then A contains the first part of the real Schur */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
/* >          On entry, the matrix B in the pair (A,B). */
/* >          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V' */
/* >          or both, then B contains the second part of the real Schur */
/* >          form of the "balanced" versions of the input A and B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of B.  LDB >= max(1,N). */
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
/* >          be the generalized eigenvalues.  If ALPHAI(j) is zero, then */
/* >          the j-th eigenvalue is real; if positive, then the j-th and */
/* >          (j+1)-st eigenvalues are a complex conjugate pair, with */
/* >          ALPHAI(j+1) negative. */
/* > */
/* >          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) */
/* >          may easily over- or underflow, and BETA(j) may even be zero. */
/* >          Thus, the user should avoid naively computing the ratio */
/* >          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less */
/* >          than and usually comparable with norm(A) in magnitude, and */
/* >          BETA always less than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* >          after another in the columns of VL, in the same order as */
/* >          their eigenvalues. If the j-th eigenvalue is real, then */
/* >          u(j) = VL(:,j), the j-th column of VL. If the j-th and */
/* >          (j+1)-th eigenvalues form a complex conjugate pair, then */
/* >          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1). */
/* >          Each eigenvector will be scaled so the largest component have */
/* >          abs(real part) + abs(imag. part) = 1. */
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
/* >          VR is DOUBLE PRECISION array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* >          after another in the columns of VR, in the same order as */
/* >          their eigenvalues. If the j-th eigenvalue is real, then */
/* >          v(j) = VR(:,j), the j-th column of VR. If the j-th and */
/* >          (j+1)-th eigenvalues form a complex conjugate pair, then */
/* >          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1). */
/* >          Each eigenvector will be scaled so the largest component have */
/* >          abs(real part) + abs(imag. part) = 1. */
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
/* >          For a complex conjugate pair of eigenvalues two consecutive */
/* >          elements of RCONDE are set to the same value. Thus RCONDE(j), */
/* >          RCONDV(j), and the j-th columns of VL and VR all correspond */
/* >          to the j-th eigenpair. */
/* >          If SENSE = 'N or 'V', RCONDE is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is DOUBLE PRECISION array, dimension (N) */
/* >          If SENSE = 'V' or 'B', the estimated reciprocal condition */
/* >          numbers of the eigenvectors, stored in consecutive elements */
/* >          of the array. For a complex eigenvector two consecutive */
/* >          elements of RCONDV are set to the same value. If the */
/* >          eigenvalues cannot be reordered to compute RCONDV(j), */
/* >          RCONDV(j) is set to 0; this can only occur when the true */
/* >          value would be very small anyway. */
/* >          If SENSE = 'N' or 'E', RCONDV is not referenced. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,2*N). */
/* >          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V', */
/* >          LWORK >= max(1,6*N). */
/* >          If SENSE = 'E' or 'B', LWORK >= max(1,10*N). */
/* >          If SENSE = 'V' or 'B', LWORK >= 2*N*N+8*N+16. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N+6) */
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
/* >                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) */
/* >                should be correct for j=INFO+1,...,N. */
/* >          > N:  =N+1: other than QZ iteration failed in DHGEQZ. */
/* >                =N+2: error return from DTGEVC. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup doubleGEeigen */

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
/* Subroutine */ int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, logical *
	bwork, integer *info, ftnlen balanc_len, ftnlen jobvl_len, ftnlen 
	jobvr_len, ftnlen sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, m, jc, in, mm, jr;
    static doublereal eps;
    static logical ilv, pair;
    static doublereal anrm, bnrm;
    static integer ierr, itau;
    static doublereal temp;
    static logical ilvl, ilvr;
    static integer iwrk, iwrk1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer icols;
    static logical noscl;
    static integer irows;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dggbak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dggbal_(char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);
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
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ldumma[1];
    static char chtemp[1];
    static doublereal bignum;
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dtgevc_(char *, char *, 
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static integer ijobvl;
    extern /* Subroutine */ int dtgsna_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ijobvr;
    static logical wantsb;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal anrmto;
    static logical wantse;
    static doublereal bnrmto;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    static logical lquery, wantsv;


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
/*     .. Executable Statements .. */

/*     Decode the input arguments */

#line 450 "dggevx.f"
    /* Parameter adjustments */
#line 450 "dggevx.f"
    a_dim1 = *lda;
#line 450 "dggevx.f"
    a_offset = 1 + a_dim1;
#line 450 "dggevx.f"
    a -= a_offset;
#line 450 "dggevx.f"
    b_dim1 = *ldb;
#line 450 "dggevx.f"
    b_offset = 1 + b_dim1;
#line 450 "dggevx.f"
    b -= b_offset;
#line 450 "dggevx.f"
    --alphar;
#line 450 "dggevx.f"
    --alphai;
#line 450 "dggevx.f"
    --beta;
#line 450 "dggevx.f"
    vl_dim1 = *ldvl;
#line 450 "dggevx.f"
    vl_offset = 1 + vl_dim1;
#line 450 "dggevx.f"
    vl -= vl_offset;
#line 450 "dggevx.f"
    vr_dim1 = *ldvr;
#line 450 "dggevx.f"
    vr_offset = 1 + vr_dim1;
#line 450 "dggevx.f"
    vr -= vr_offset;
#line 450 "dggevx.f"
    --lscale;
#line 450 "dggevx.f"
    --rscale;
#line 450 "dggevx.f"
    --rconde;
#line 450 "dggevx.f"
    --rcondv;
#line 450 "dggevx.f"
    --work;
#line 450 "dggevx.f"
    --iwork;
#line 450 "dggevx.f"
    --bwork;
#line 450 "dggevx.f"

#line 450 "dggevx.f"
    /* Function Body */
#line 450 "dggevx.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 451 "dggevx.f"
	ijobvl = 1;
#line 452 "dggevx.f"
	ilvl = FALSE_;
#line 453 "dggevx.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 454 "dggevx.f"
	ijobvl = 2;
#line 455 "dggevx.f"
	ilvl = TRUE_;
#line 456 "dggevx.f"
    } else {
#line 457 "dggevx.f"
	ijobvl = -1;
#line 458 "dggevx.f"
	ilvl = FALSE_;
#line 459 "dggevx.f"
    }

#line 461 "dggevx.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 462 "dggevx.f"
	ijobvr = 1;
#line 463 "dggevx.f"
	ilvr = FALSE_;
#line 464 "dggevx.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 465 "dggevx.f"
	ijobvr = 2;
#line 466 "dggevx.f"
	ilvr = TRUE_;
#line 467 "dggevx.f"
    } else {
#line 468 "dggevx.f"
	ijobvr = -1;
#line 469 "dggevx.f"
	ilvr = FALSE_;
#line 470 "dggevx.f"
    }
#line 471 "dggevx.f"
    ilv = ilvl || ilvr;

#line 473 "dggevx.f"
    noscl = lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (
	    ftnlen)1, (ftnlen)1);
#line 474 "dggevx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 475 "dggevx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 476 "dggevx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 477 "dggevx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 481 "dggevx.f"
    *info = 0;
#line 482 "dggevx.f"
    lquery = *lwork == -1;
#line 483 "dggevx.f"
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "S", (
	    ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) 
	    || lsame_(balanc, "B", (ftnlen)1, (ftnlen)1))) {
#line 486 "dggevx.f"
	*info = -1;
#line 487 "dggevx.f"
    } else if (ijobvl <= 0) {
#line 488 "dggevx.f"
	*info = -2;
#line 489 "dggevx.f"
    } else if (ijobvr <= 0) {
#line 490 "dggevx.f"
	*info = -3;
#line 491 "dggevx.f"
    } else if (! (wantsn || wantse || wantsb || wantsv)) {
#line 493 "dggevx.f"
	*info = -4;
#line 494 "dggevx.f"
    } else if (*n < 0) {
#line 495 "dggevx.f"
	*info = -5;
#line 496 "dggevx.f"
    } else if (*lda < max(1,*n)) {
#line 497 "dggevx.f"
	*info = -7;
#line 498 "dggevx.f"
    } else if (*ldb < max(1,*n)) {
#line 499 "dggevx.f"
	*info = -9;
#line 500 "dggevx.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 501 "dggevx.f"
	*info = -14;
#line 502 "dggevx.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 503 "dggevx.f"
	*info = -16;
#line 504 "dggevx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. The workspace is */
/*       computed assuming ILO = 1 and IHI = N, the worst case.) */

#line 514 "dggevx.f"
    if (*info == 0) {
#line 515 "dggevx.f"
	if (*n == 0) {
#line 516 "dggevx.f"
	    minwrk = 1;
#line 517 "dggevx.f"
	    maxwrk = 1;
#line 518 "dggevx.f"
	} else {
#line 519 "dggevx.f"
	    if (noscl && ! ilv) {
#line 520 "dggevx.f"
		minwrk = *n << 1;
#line 521 "dggevx.f"
	    } else {
#line 522 "dggevx.f"
		minwrk = *n * 6;
#line 523 "dggevx.f"
	    }
#line 524 "dggevx.f"
	    if (wantse || wantsb) {
#line 525 "dggevx.f"
		minwrk = *n * 10;
#line 526 "dggevx.f"
	    }
#line 527 "dggevx.f"
	    if (wantsv || wantsb) {
/* Computing MAX */
#line 528 "dggevx.f"
		i__1 = minwrk, i__2 = (*n << 1) * (*n + 4) + 16;
#line 528 "dggevx.f"
		minwrk = max(i__1,i__2);
#line 529 "dggevx.f"
	    }
#line 530 "dggevx.f"
	    maxwrk = minwrk;
/* Computing MAX */
#line 531 "dggevx.f"
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 531 "dggevx.f"
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 533 "dggevx.f"
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "DORMQR", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 533 "dggevx.f"
	    maxwrk = max(i__1,i__2);
#line 535 "dggevx.f"
	    if (ilvl) {
/* Computing MAX */
#line 536 "dggevx.f"
		i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "DORGQR", 
			" ", n, &c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 536 "dggevx.f"
		maxwrk = max(i__1,i__2);
#line 538 "dggevx.f"
	    }
#line 539 "dggevx.f"
	}
#line 540 "dggevx.f"
	work[1] = (doublereal) maxwrk;

#line 542 "dggevx.f"
	if (*lwork < minwrk && ! lquery) {
#line 543 "dggevx.f"
	    *info = -26;
#line 544 "dggevx.f"
	}
#line 545 "dggevx.f"
    }

#line 547 "dggevx.f"
    if (*info != 0) {
#line 548 "dggevx.f"
	i__1 = -(*info);
#line 548 "dggevx.f"
	xerbla_("DGGEVX", &i__1, (ftnlen)6);
#line 549 "dggevx.f"
	return 0;
#line 550 "dggevx.f"
    } else if (lquery) {
#line 551 "dggevx.f"
	return 0;
#line 552 "dggevx.f"
    }

/*     Quick return if possible */

#line 556 "dggevx.f"
    if (*n == 0) {
#line 556 "dggevx.f"
	return 0;
#line 556 "dggevx.f"
    }


/*     Get machine constants */

#line 562 "dggevx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 563 "dggevx.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 564 "dggevx.f"
    bignum = 1. / smlnum;
#line 565 "dggevx.f"
    dlabad_(&smlnum, &bignum);
#line 566 "dggevx.f"
    smlnum = sqrt(smlnum) / eps;
#line 567 "dggevx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 571 "dggevx.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 572 "dggevx.f"
    ilascl = FALSE_;
#line 573 "dggevx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 574 "dggevx.f"
	anrmto = smlnum;
#line 575 "dggevx.f"
	ilascl = TRUE_;
#line 576 "dggevx.f"
    } else if (anrm > bignum) {
#line 577 "dggevx.f"
	anrmto = bignum;
#line 578 "dggevx.f"
	ilascl = TRUE_;
#line 579 "dggevx.f"
    }
#line 580 "dggevx.f"
    if (ilascl) {
#line 580 "dggevx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 580 "dggevx.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 585 "dggevx.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 586 "dggevx.f"
    ilbscl = FALSE_;
#line 587 "dggevx.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 588 "dggevx.f"
	bnrmto = smlnum;
#line 589 "dggevx.f"
	ilbscl = TRUE_;
#line 590 "dggevx.f"
    } else if (bnrm > bignum) {
#line 591 "dggevx.f"
	bnrmto = bignum;
#line 592 "dggevx.f"
	ilbscl = TRUE_;
#line 593 "dggevx.f"
    }
#line 594 "dggevx.f"
    if (ilbscl) {
#line 594 "dggevx.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 594 "dggevx.f"
    }

/*     Permute and/or balance the matrix pair (A,B) */
/*     (Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise) */

#line 600 "dggevx.f"
    dggbal_(balanc, n, &a[a_offset], lda, &b[b_offset], ldb, ilo, ihi, &
	    lscale[1], &rscale[1], &work[1], &ierr, (ftnlen)1);

/*     Compute ABNRM and BBNRM */

#line 605 "dggevx.f"
    *abnrm = dlange_("1", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 606 "dggevx.f"
    if (ilascl) {
#line 607 "dggevx.f"
	work[1] = *abnrm;
#line 608 "dggevx.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, &c__1, &c__1, &work[1], &
		c__1, &ierr, (ftnlen)1);
#line 610 "dggevx.f"
	*abnrm = work[1];
#line 611 "dggevx.f"
    }

#line 613 "dggevx.f"
    *bbnrm = dlange_("1", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 614 "dggevx.f"
    if (ilbscl) {
#line 615 "dggevx.f"
	work[1] = *bbnrm;
#line 616 "dggevx.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, &c__1, &c__1, &work[1], &
		c__1, &ierr, (ftnlen)1);
#line 618 "dggevx.f"
	*bbnrm = work[1];
#line 619 "dggevx.f"
    }

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Workspace: need N, prefer N*NB ) */

#line 624 "dggevx.f"
    irows = *ihi + 1 - *ilo;
#line 625 "dggevx.f"
    if (ilv || ! wantsn) {
#line 626 "dggevx.f"
	icols = *n + 1 - *ilo;
#line 627 "dggevx.f"
    } else {
#line 628 "dggevx.f"
	icols = irows;
#line 629 "dggevx.f"
    }
#line 630 "dggevx.f"
    itau = 1;
#line 631 "dggevx.f"
    iwrk = itau + irows;
#line 632 "dggevx.f"
    i__1 = *lwork + 1 - iwrk;
#line 632 "dggevx.f"
    dgeqrf_(&irows, &icols, &b[*ilo + *ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to A */
/*     (Workspace: need N, prefer N*NB) */

#line 638 "dggevx.f"
    i__1 = *lwork + 1 - iwrk;
#line 638 "dggevx.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[*ilo + *ilo * b_dim1], ldb, &
	    work[itau], &a[*ilo + *ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL and/or VR */
/*     (Workspace: need N, prefer N*NB) */

#line 645 "dggevx.f"
    if (ilvl) {
#line 646 "dggevx.f"
	dlaset_("Full", n, n, &c_b59, &c_b60, &vl[vl_offset], ldvl, (ftnlen)4)
		;
#line 647 "dggevx.f"
	if (irows > 1) {
#line 648 "dggevx.f"
	    i__1 = irows - 1;
#line 648 "dggevx.f"
	    i__2 = irows - 1;
#line 648 "dggevx.f"
	    dlacpy_("L", &i__1, &i__2, &b[*ilo + 1 + *ilo * b_dim1], ldb, &vl[
		    *ilo + 1 + *ilo * vl_dim1], ldvl, (ftnlen)1);
#line 650 "dggevx.f"
	}
#line 651 "dggevx.f"
	i__1 = *lwork + 1 - iwrk;
#line 651 "dggevx.f"
	dorgqr_(&irows, &irows, &irows, &vl[*ilo + *ilo * vl_dim1], ldvl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 653 "dggevx.f"
    }

#line 655 "dggevx.f"
    if (ilvr) {
#line 655 "dggevx.f"
	dlaset_("Full", n, n, &c_b59, &c_b60, &vr[vr_offset], ldvr, (ftnlen)4)
		;
#line 655 "dggevx.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 661 "dggevx.f"
    if (ilv || ! wantsn) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 665 "dggevx.f"
	dgghrd_(jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr, (
		ftnlen)1, (ftnlen)1);
#line 667 "dggevx.f"
    } else {
#line 668 "dggevx.f"
	dgghrd_("N", "N", &irows, &c__1, &irows, &a[*ilo + *ilo * a_dim1], 
		lda, &b[*ilo + *ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 670 "dggevx.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur forms and Schur vectors) */
/*     (Workspace: need N) */

#line 676 "dggevx.f"
    if (ilv || ! wantsn) {
#line 677 "dggevx.f"
	*(unsigned char *)chtemp = 'S';
#line 678 "dggevx.f"
    } else {
#line 679 "dggevx.f"
	*(unsigned char *)chtemp = 'E';
#line 680 "dggevx.f"
    }

#line 682 "dggevx.f"
    dhgeqz_(chtemp, jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[b_offset]
	    , ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], ldvl, &
	    vr[vr_offset], ldvr, &work[1], lwork, &ierr, (ftnlen)1, (ftnlen)1,
	     (ftnlen)1);
#line 685 "dggevx.f"
    if (ierr != 0) {
#line 686 "dggevx.f"
	if (ierr > 0 && ierr <= *n) {
#line 687 "dggevx.f"
	    *info = ierr;
#line 688 "dggevx.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 689 "dggevx.f"
	    *info = ierr - *n;
#line 690 "dggevx.f"
	} else {
#line 691 "dggevx.f"
	    *info = *n + 1;
#line 692 "dggevx.f"
	}
#line 693 "dggevx.f"
	goto L130;
#line 694 "dggevx.f"
    }

/*     Compute Eigenvectors and estimate condition numbers if desired */
/*     (Workspace: DTGEVC: need 6*N */
/*                 DTGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B', */
/*                         need N otherwise ) */

#line 701 "dggevx.f"
    if (ilv || ! wantsn) {
#line 702 "dggevx.f"
	if (ilv) {
#line 703 "dggevx.f"
	    if (ilvl) {
#line 704 "dggevx.f"
		if (ilvr) {
#line 705 "dggevx.f"
		    *(unsigned char *)chtemp = 'B';
#line 706 "dggevx.f"
		} else {
#line 707 "dggevx.f"
		    *(unsigned char *)chtemp = 'L';
#line 708 "dggevx.f"
		}
#line 709 "dggevx.f"
	    } else {
#line 710 "dggevx.f"
		*(unsigned char *)chtemp = 'R';
#line 711 "dggevx.f"
	    }

#line 713 "dggevx.f"
	    dtgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], 
		    ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &
		    work[1], &ierr, (ftnlen)1, (ftnlen)1);
#line 715 "dggevx.f"
	    if (ierr != 0) {
#line 716 "dggevx.f"
		*info = *n + 2;
#line 717 "dggevx.f"
		goto L130;
#line 718 "dggevx.f"
	    }
#line 719 "dggevx.f"
	}

#line 721 "dggevx.f"
	if (! wantsn) {

/*           compute eigenvectors (DTGEVC) and estimate condition */
/*           numbers (DTGSNA). Note that the definition of the condition */
/*           number is not invariant under transformation (u,v) to */
/*           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized */
/*           Schur form (S,T), Q and Z are orthogonal matrices. In order */
/*           to avoid using extra 2*N*N workspace, we have to recalculate */
/*           eigenvectors and estimate one condition numbers at a time. */

#line 731 "dggevx.f"
	    pair = FALSE_;
#line 732 "dggevx.f"
	    i__1 = *n;
#line 732 "dggevx.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

#line 734 "dggevx.f"
		if (pair) {
#line 735 "dggevx.f"
		    pair = FALSE_;
#line 736 "dggevx.f"
		    goto L20;
#line 737 "dggevx.f"
		}
#line 738 "dggevx.f"
		mm = 1;
#line 739 "dggevx.f"
		if (i__ < *n) {
#line 740 "dggevx.f"
		    if (a[i__ + 1 + i__ * a_dim1] != 0.) {
#line 741 "dggevx.f"
			pair = TRUE_;
#line 742 "dggevx.f"
			mm = 2;
#line 743 "dggevx.f"
		    }
#line 744 "dggevx.f"
		}

#line 746 "dggevx.f"
		i__2 = *n;
#line 746 "dggevx.f"
		for (j = 1; j <= i__2; ++j) {
#line 747 "dggevx.f"
		    bwork[j] = FALSE_;
#line 748 "dggevx.f"
/* L10: */
#line 748 "dggevx.f"
		}
#line 749 "dggevx.f"
		if (mm == 1) {
#line 750 "dggevx.f"
		    bwork[i__] = TRUE_;
#line 751 "dggevx.f"
		} else if (mm == 2) {
#line 752 "dggevx.f"
		    bwork[i__] = TRUE_;
#line 753 "dggevx.f"
		    bwork[i__ + 1] = TRUE_;
#line 754 "dggevx.f"
		}

#line 756 "dggevx.f"
		iwrk = mm * *n + 1;
#line 757 "dggevx.f"
		iwrk1 = iwrk + mm * *n;

/*              Compute a pair of left and right eigenvectors. */
/*              (compute workspace: need up to 4*N + 6*N) */

#line 762 "dggevx.f"
		if (wantse || wantsb) {
#line 763 "dggevx.f"
		    dtgevc_("B", "S", &bwork[1], n, &a[a_offset], lda, &b[
			    b_offset], ldb, &work[1], n, &work[iwrk], n, &mm, 
			    &m, &work[iwrk1], &ierr, (ftnlen)1, (ftnlen)1);
#line 766 "dggevx.f"
		    if (ierr != 0) {
#line 767 "dggevx.f"
			*info = *n + 2;
#line 768 "dggevx.f"
			goto L130;
#line 769 "dggevx.f"
		    }
#line 770 "dggevx.f"
		}

#line 772 "dggevx.f"
		i__2 = *lwork - iwrk1 + 1;
#line 772 "dggevx.f"
		dtgsna_(sense, "S", &bwork[1], n, &a[a_offset], lda, &b[
			b_offset], ldb, &work[1], n, &work[iwrk], n, &rconde[
			i__], &rcondv[i__], &mm, &m, &work[iwrk1], &i__2, &
			iwork[1], &ierr, (ftnlen)1, (ftnlen)1);

#line 777 "dggevx.f"
L20:
#line 777 "dggevx.f"
		;
#line 777 "dggevx.f"
	    }
#line 778 "dggevx.f"
	}
#line 779 "dggevx.f"
    }

/*     Undo balancing on VL and VR and normalization */
/*     (Workspace: none needed) */

#line 784 "dggevx.f"
    if (ilvl) {
#line 785 "dggevx.f"
	dggbak_(balanc, "L", n, ilo, ihi, &lscale[1], &rscale[1], n, &vl[
		vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);

#line 788 "dggevx.f"
	i__1 = *n;
#line 788 "dggevx.f"
	for (jc = 1; jc <= i__1; ++jc) {
#line 789 "dggevx.f"
	    if (alphai[jc] < 0.) {
#line 789 "dggevx.f"
		goto L70;
#line 789 "dggevx.f"
	    }
#line 791 "dggevx.f"
	    temp = 0.;
#line 792 "dggevx.f"
	    if (alphai[jc] == 0.) {
#line 793 "dggevx.f"
		i__2 = *n;
#line 793 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 794 "dggevx.f"
		    d__2 = temp, d__3 = (d__1 = vl[jr + jc * vl_dim1], abs(
			    d__1));
#line 794 "dggevx.f"
		    temp = max(d__2,d__3);
#line 795 "dggevx.f"
/* L30: */
#line 795 "dggevx.f"
		}
#line 796 "dggevx.f"
	    } else {
#line 797 "dggevx.f"
		i__2 = *n;
#line 797 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 798 "dggevx.f"
		    d__3 = temp, d__4 = (d__1 = vl[jr + jc * vl_dim1], abs(
			    d__1)) + (d__2 = vl[jr + (jc + 1) * vl_dim1], abs(
			    d__2));
#line 798 "dggevx.f"
		    temp = max(d__3,d__4);
#line 800 "dggevx.f"
/* L40: */
#line 800 "dggevx.f"
		}
#line 801 "dggevx.f"
	    }
#line 802 "dggevx.f"
	    if (temp < smlnum) {
#line 802 "dggevx.f"
		goto L70;
#line 802 "dggevx.f"
	    }
#line 804 "dggevx.f"
	    temp = 1. / temp;
#line 805 "dggevx.f"
	    if (alphai[jc] == 0.) {
#line 806 "dggevx.f"
		i__2 = *n;
#line 806 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 807 "dggevx.f"
		    vl[jr + jc * vl_dim1] *= temp;
#line 808 "dggevx.f"
/* L50: */
#line 808 "dggevx.f"
		}
#line 809 "dggevx.f"
	    } else {
#line 810 "dggevx.f"
		i__2 = *n;
#line 810 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 811 "dggevx.f"
		    vl[jr + jc * vl_dim1] *= temp;
#line 812 "dggevx.f"
		    vl[jr + (jc + 1) * vl_dim1] *= temp;
#line 813 "dggevx.f"
/* L60: */
#line 813 "dggevx.f"
		}
#line 814 "dggevx.f"
	    }
#line 815 "dggevx.f"
L70:
#line 815 "dggevx.f"
	    ;
#line 815 "dggevx.f"
	}
#line 816 "dggevx.f"
    }
#line 817 "dggevx.f"
    if (ilvr) {
#line 818 "dggevx.f"
	dggbak_(balanc, "R", n, ilo, ihi, &lscale[1], &rscale[1], n, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 820 "dggevx.f"
	i__1 = *n;
#line 820 "dggevx.f"
	for (jc = 1; jc <= i__1; ++jc) {
#line 821 "dggevx.f"
	    if (alphai[jc] < 0.) {
#line 821 "dggevx.f"
		goto L120;
#line 821 "dggevx.f"
	    }
#line 823 "dggevx.f"
	    temp = 0.;
#line 824 "dggevx.f"
	    if (alphai[jc] == 0.) {
#line 825 "dggevx.f"
		i__2 = *n;
#line 825 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 826 "dggevx.f"
		    d__2 = temp, d__3 = (d__1 = vr[jr + jc * vr_dim1], abs(
			    d__1));
#line 826 "dggevx.f"
		    temp = max(d__2,d__3);
#line 827 "dggevx.f"
/* L80: */
#line 827 "dggevx.f"
		}
#line 828 "dggevx.f"
	    } else {
#line 829 "dggevx.f"
		i__2 = *n;
#line 829 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 830 "dggevx.f"
		    d__3 = temp, d__4 = (d__1 = vr[jr + jc * vr_dim1], abs(
			    d__1)) + (d__2 = vr[jr + (jc + 1) * vr_dim1], abs(
			    d__2));
#line 830 "dggevx.f"
		    temp = max(d__3,d__4);
#line 832 "dggevx.f"
/* L90: */
#line 832 "dggevx.f"
		}
#line 833 "dggevx.f"
	    }
#line 834 "dggevx.f"
	    if (temp < smlnum) {
#line 834 "dggevx.f"
		goto L120;
#line 834 "dggevx.f"
	    }
#line 836 "dggevx.f"
	    temp = 1. / temp;
#line 837 "dggevx.f"
	    if (alphai[jc] == 0.) {
#line 838 "dggevx.f"
		i__2 = *n;
#line 838 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 839 "dggevx.f"
		    vr[jr + jc * vr_dim1] *= temp;
#line 840 "dggevx.f"
/* L100: */
#line 840 "dggevx.f"
		}
#line 841 "dggevx.f"
	    } else {
#line 842 "dggevx.f"
		i__2 = *n;
#line 842 "dggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 843 "dggevx.f"
		    vr[jr + jc * vr_dim1] *= temp;
#line 844 "dggevx.f"
		    vr[jr + (jc + 1) * vr_dim1] *= temp;
#line 845 "dggevx.f"
/* L110: */
#line 845 "dggevx.f"
		}
#line 846 "dggevx.f"
	    }
#line 847 "dggevx.f"
L120:
#line 847 "dggevx.f"
	    ;
#line 847 "dggevx.f"
	}
#line 848 "dggevx.f"
    }

/*     Undo scaling if necessary */

#line 852 "dggevx.f"
L130:

#line 854 "dggevx.f"
    if (ilascl) {
#line 855 "dggevx.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 856 "dggevx.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 857 "dggevx.f"
    }

#line 859 "dggevx.f"
    if (ilbscl) {
#line 860 "dggevx.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 861 "dggevx.f"
    }

#line 863 "dggevx.f"
    work[1] = (doublereal) maxwrk;
#line 864 "dggevx.f"
    return 0;

/*     End of DGGEVX */

} /* dggevx_ */

