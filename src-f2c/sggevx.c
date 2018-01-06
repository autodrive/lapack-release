#line 1 "sggevx.f"
/* sggevx.f -- translated by f2c (version 20100827).
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

#line 1 "sggevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b57 = 0.;
static doublereal c_b58 = 1.;

/* > \brief <b> SGGEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, */
/*                          ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, */
/*                          IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, */
/*                          RCONDV, WORK, LWORK, IWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       REAL               ABNRM, BBNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), LSCALE( * ), */
/*      $                   RCONDE( * ), RCONDV( * ), RSCALE( * ), */
/*      $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B) */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          B is REAL array, dimension (LDB, N) */
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
/* >          VL is REAL array, dimension (LDVL,N) */
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
/* >          VR is REAL array, dimension (LDVR,N) */
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
/* >          LSCALE is REAL array, dimension (N) */
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
/* >          RSCALE is REAL array, dimension (N) */
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
/* >          ABNRM is REAL */
/* >          The one-norm of the balanced matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] BBNRM */
/* > \verbatim */
/* >          BBNRM is REAL */
/* >          The one-norm of the balanced matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDE */
/* > \verbatim */
/* >          RCONDE is REAL array, dimension (N) */
/* >          If SENSE = 'E' or 'B', the reciprocal condition numbers of */
/* >          the eigenvalues, stored in consecutive elements of the array. */
/* >          For a complex conjugate pair of eigenvalues two consecutive */
/* >          elements of RCONDE are set to the same value. Thus RCONDE(j), */
/* >          RCONDV(j), and the j-th columns of VL and VR all correspond */
/* >          to the j-th eigenpair. */
/* >          If SENSE = 'N' or 'V', RCONDE is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is REAL array, dimension (N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= max(1,2*N). */
/* >          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V', */
/* >          LWORK >= max(1,6*N). */
/* >          If SENSE = 'E', LWORK >= max(1,10*N). */
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
/* >          > N:  =N+1: other than QZ iteration failed in SHGEQZ. */
/* >                =N+2: error return from STGEVC. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup realGEeigen */

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
/* Subroutine */ int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
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
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *), sggbak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), sggbal_(char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical ilascl, ilbscl;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), sgghrd_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical ldumma[1];
    static char chtemp[1];
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal slange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    static integer ijobvl;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer ijobvr;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static logical wantsb;
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal anrmto;
    static logical wantse;
    static doublereal bnrmto;
    extern /* Subroutine */ int shgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), stgevc_(char *, char *, 
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    stgsna_(char *, char *, logical *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen, 
	    ftnlen);
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical lquery, wantsv;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
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

#line 450 "sggevx.f"
    /* Parameter adjustments */
#line 450 "sggevx.f"
    a_dim1 = *lda;
#line 450 "sggevx.f"
    a_offset = 1 + a_dim1;
#line 450 "sggevx.f"
    a -= a_offset;
#line 450 "sggevx.f"
    b_dim1 = *ldb;
#line 450 "sggevx.f"
    b_offset = 1 + b_dim1;
#line 450 "sggevx.f"
    b -= b_offset;
#line 450 "sggevx.f"
    --alphar;
#line 450 "sggevx.f"
    --alphai;
#line 450 "sggevx.f"
    --beta;
#line 450 "sggevx.f"
    vl_dim1 = *ldvl;
#line 450 "sggevx.f"
    vl_offset = 1 + vl_dim1;
#line 450 "sggevx.f"
    vl -= vl_offset;
#line 450 "sggevx.f"
    vr_dim1 = *ldvr;
#line 450 "sggevx.f"
    vr_offset = 1 + vr_dim1;
#line 450 "sggevx.f"
    vr -= vr_offset;
#line 450 "sggevx.f"
    --lscale;
#line 450 "sggevx.f"
    --rscale;
#line 450 "sggevx.f"
    --rconde;
#line 450 "sggevx.f"
    --rcondv;
#line 450 "sggevx.f"
    --work;
#line 450 "sggevx.f"
    --iwork;
#line 450 "sggevx.f"
    --bwork;
#line 450 "sggevx.f"

#line 450 "sggevx.f"
    /* Function Body */
#line 450 "sggevx.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 451 "sggevx.f"
	ijobvl = 1;
#line 452 "sggevx.f"
	ilvl = FALSE_;
#line 453 "sggevx.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 454 "sggevx.f"
	ijobvl = 2;
#line 455 "sggevx.f"
	ilvl = TRUE_;
#line 456 "sggevx.f"
    } else {
#line 457 "sggevx.f"
	ijobvl = -1;
#line 458 "sggevx.f"
	ilvl = FALSE_;
#line 459 "sggevx.f"
    }

#line 461 "sggevx.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 462 "sggevx.f"
	ijobvr = 1;
#line 463 "sggevx.f"
	ilvr = FALSE_;
#line 464 "sggevx.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 465 "sggevx.f"
	ijobvr = 2;
#line 466 "sggevx.f"
	ilvr = TRUE_;
#line 467 "sggevx.f"
    } else {
#line 468 "sggevx.f"
	ijobvr = -1;
#line 469 "sggevx.f"
	ilvr = FALSE_;
#line 470 "sggevx.f"
    }
#line 471 "sggevx.f"
    ilv = ilvl || ilvr;

#line 473 "sggevx.f"
    noscl = lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (
	    ftnlen)1, (ftnlen)1);
#line 474 "sggevx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 475 "sggevx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 476 "sggevx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 477 "sggevx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);

/*     Test the input arguments */

#line 481 "sggevx.f"
    *info = 0;
#line 482 "sggevx.f"
    lquery = *lwork == -1;
#line 483 "sggevx.f"
    if (! (noscl || lsame_(balanc, "S", (ftnlen)1, (ftnlen)1) || lsame_(
	    balanc, "B", (ftnlen)1, (ftnlen)1))) {
#line 485 "sggevx.f"
	*info = -1;
#line 486 "sggevx.f"
    } else if (ijobvl <= 0) {
#line 487 "sggevx.f"
	*info = -2;
#line 488 "sggevx.f"
    } else if (ijobvr <= 0) {
#line 489 "sggevx.f"
	*info = -3;
#line 490 "sggevx.f"
    } else if (! (wantsn || wantse || wantsb || wantsv)) {
#line 492 "sggevx.f"
	*info = -4;
#line 493 "sggevx.f"
    } else if (*n < 0) {
#line 494 "sggevx.f"
	*info = -5;
#line 495 "sggevx.f"
    } else if (*lda < max(1,*n)) {
#line 496 "sggevx.f"
	*info = -7;
#line 497 "sggevx.f"
    } else if (*ldb < max(1,*n)) {
#line 498 "sggevx.f"
	*info = -9;
#line 499 "sggevx.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 500 "sggevx.f"
	*info = -14;
#line 501 "sggevx.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 502 "sggevx.f"
	*info = -16;
#line 503 "sggevx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. The workspace is */
/*       computed assuming ILO = 1 and IHI = N, the worst case.) */

#line 513 "sggevx.f"
    if (*info == 0) {
#line 514 "sggevx.f"
	if (*n == 0) {
#line 515 "sggevx.f"
	    minwrk = 1;
#line 516 "sggevx.f"
	    maxwrk = 1;
#line 517 "sggevx.f"
	} else {
#line 518 "sggevx.f"
	    if (noscl && ! ilv) {
#line 519 "sggevx.f"
		minwrk = *n << 1;
#line 520 "sggevx.f"
	    } else {
#line 521 "sggevx.f"
		minwrk = *n * 6;
#line 522 "sggevx.f"
	    }
#line 523 "sggevx.f"
	    if (wantse) {
#line 524 "sggevx.f"
		minwrk = *n * 10;
#line 525 "sggevx.f"
	    } else if (wantsv || wantsb) {
#line 526 "sggevx.f"
		minwrk = (*n << 1) * (*n + 4) + 16;
#line 527 "sggevx.f"
	    }
#line 528 "sggevx.f"
	    maxwrk = minwrk;
/* Computing MAX */
#line 529 "sggevx.f"
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 529 "sggevx.f"
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 531 "sggevx.f"
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "SORMQR", " ", n, &
		    c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 531 "sggevx.f"
	    maxwrk = max(i__1,i__2);
#line 533 "sggevx.f"
	    if (ilvl) {
/* Computing MAX */
#line 534 "sggevx.f"
		i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "SORGQR", 
			" ", n, &c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 534 "sggevx.f"
		maxwrk = max(i__1,i__2);
#line 536 "sggevx.f"
	    }
#line 537 "sggevx.f"
	}
#line 538 "sggevx.f"
	work[1] = (doublereal) maxwrk;

#line 540 "sggevx.f"
	if (*lwork < minwrk && ! lquery) {
#line 541 "sggevx.f"
	    *info = -26;
#line 542 "sggevx.f"
	}
#line 543 "sggevx.f"
    }

#line 545 "sggevx.f"
    if (*info != 0) {
#line 546 "sggevx.f"
	i__1 = -(*info);
#line 546 "sggevx.f"
	xerbla_("SGGEVX", &i__1, (ftnlen)6);
#line 547 "sggevx.f"
	return 0;
#line 548 "sggevx.f"
    } else if (lquery) {
#line 549 "sggevx.f"
	return 0;
#line 550 "sggevx.f"
    }

/*     Quick return if possible */

#line 554 "sggevx.f"
    if (*n == 0) {
#line 554 "sggevx.f"
	return 0;
#line 554 "sggevx.f"
    }


/*     Get machine constants */

#line 560 "sggevx.f"
    eps = slamch_("P", (ftnlen)1);
#line 561 "sggevx.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 562 "sggevx.f"
    bignum = 1. / smlnum;
#line 563 "sggevx.f"
    slabad_(&smlnum, &bignum);
#line 564 "sggevx.f"
    smlnum = sqrt(smlnum) / eps;
#line 565 "sggevx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 569 "sggevx.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 570 "sggevx.f"
    ilascl = FALSE_;
#line 571 "sggevx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 572 "sggevx.f"
	anrmto = smlnum;
#line 573 "sggevx.f"
	ilascl = TRUE_;
#line 574 "sggevx.f"
    } else if (anrm > bignum) {
#line 575 "sggevx.f"
	anrmto = bignum;
#line 576 "sggevx.f"
	ilascl = TRUE_;
#line 577 "sggevx.f"
    }
#line 578 "sggevx.f"
    if (ilascl) {
#line 578 "sggevx.f"
	slascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 578 "sggevx.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 583 "sggevx.f"
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 584 "sggevx.f"
    ilbscl = FALSE_;
#line 585 "sggevx.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 586 "sggevx.f"
	bnrmto = smlnum;
#line 587 "sggevx.f"
	ilbscl = TRUE_;
#line 588 "sggevx.f"
    } else if (bnrm > bignum) {
#line 589 "sggevx.f"
	bnrmto = bignum;
#line 590 "sggevx.f"
	ilbscl = TRUE_;
#line 591 "sggevx.f"
    }
#line 592 "sggevx.f"
    if (ilbscl) {
#line 592 "sggevx.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 592 "sggevx.f"
    }

/*     Permute and/or balance the matrix pair (A,B) */
/*     (Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise) */

#line 598 "sggevx.f"
    sggbal_(balanc, n, &a[a_offset], lda, &b[b_offset], ldb, ilo, ihi, &
	    lscale[1], &rscale[1], &work[1], &ierr, (ftnlen)1);

/*     Compute ABNRM and BBNRM */

#line 603 "sggevx.f"
    *abnrm = slange_("1", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 604 "sggevx.f"
    if (ilascl) {
#line 605 "sggevx.f"
	work[1] = *abnrm;
#line 606 "sggevx.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, &c__1, &c__1, &work[1], &
		c__1, &ierr, (ftnlen)1);
#line 608 "sggevx.f"
	*abnrm = work[1];
#line 609 "sggevx.f"
    }

#line 611 "sggevx.f"
    *bbnrm = slange_("1", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 612 "sggevx.f"
    if (ilbscl) {
#line 613 "sggevx.f"
	work[1] = *bbnrm;
#line 614 "sggevx.f"
	slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, &c__1, &c__1, &work[1], &
		c__1, &ierr, (ftnlen)1);
#line 616 "sggevx.f"
	*bbnrm = work[1];
#line 617 "sggevx.f"
    }

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Workspace: need N, prefer N*NB ) */

#line 622 "sggevx.f"
    irows = *ihi + 1 - *ilo;
#line 623 "sggevx.f"
    if (ilv || ! wantsn) {
#line 624 "sggevx.f"
	icols = *n + 1 - *ilo;
#line 625 "sggevx.f"
    } else {
#line 626 "sggevx.f"
	icols = irows;
#line 627 "sggevx.f"
    }
#line 628 "sggevx.f"
    itau = 1;
#line 629 "sggevx.f"
    iwrk = itau + irows;
#line 630 "sggevx.f"
    i__1 = *lwork + 1 - iwrk;
#line 630 "sggevx.f"
    sgeqrf_(&irows, &icols, &b[*ilo + *ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to A */
/*     (Workspace: need N, prefer N*NB) */

#line 636 "sggevx.f"
    i__1 = *lwork + 1 - iwrk;
#line 636 "sggevx.f"
    sormqr_("L", "T", &irows, &icols, &irows, &b[*ilo + *ilo * b_dim1], ldb, &
	    work[itau], &a[*ilo + *ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL and/or VR */
/*     (Workspace: need N, prefer N*NB) */

#line 643 "sggevx.f"
    if (ilvl) {
#line 644 "sggevx.f"
	slaset_("Full", n, n, &c_b57, &c_b58, &vl[vl_offset], ldvl, (ftnlen)4)
		;
#line 645 "sggevx.f"
	if (irows > 1) {
#line 646 "sggevx.f"
	    i__1 = irows - 1;
#line 646 "sggevx.f"
	    i__2 = irows - 1;
#line 646 "sggevx.f"
	    slacpy_("L", &i__1, &i__2, &b[*ilo + 1 + *ilo * b_dim1], ldb, &vl[
		    *ilo + 1 + *ilo * vl_dim1], ldvl, (ftnlen)1);
#line 648 "sggevx.f"
	}
#line 649 "sggevx.f"
	i__1 = *lwork + 1 - iwrk;
#line 649 "sggevx.f"
	sorgqr_(&irows, &irows, &irows, &vl[*ilo + *ilo * vl_dim1], ldvl, &
		work[itau], &work[iwrk], &i__1, &ierr);
#line 651 "sggevx.f"
    }

#line 653 "sggevx.f"
    if (ilvr) {
#line 653 "sggevx.f"
	slaset_("Full", n, n, &c_b57, &c_b58, &vr[vr_offset], ldvr, (ftnlen)4)
		;
#line 653 "sggevx.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 659 "sggevx.f"
    if (ilv || ! wantsn) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 663 "sggevx.f"
	sgghrd_(jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr, (
		ftnlen)1, (ftnlen)1);
#line 665 "sggevx.f"
    } else {
#line 666 "sggevx.f"
	sgghrd_("N", "N", &irows, &c__1, &irows, &a[*ilo + *ilo * a_dim1], 
		lda, &b[*ilo + *ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 668 "sggevx.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur forms and Schur vectors) */
/*     (Workspace: need N) */

#line 674 "sggevx.f"
    if (ilv || ! wantsn) {
#line 675 "sggevx.f"
	*(unsigned char *)chtemp = 'S';
#line 676 "sggevx.f"
    } else {
#line 677 "sggevx.f"
	*(unsigned char *)chtemp = 'E';
#line 678 "sggevx.f"
    }

#line 680 "sggevx.f"
    shgeqz_(chtemp, jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[b_offset]
	    , ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], ldvl, &
	    vr[vr_offset], ldvr, &work[1], lwork, &ierr, (ftnlen)1, (ftnlen)1,
	     (ftnlen)1);
#line 683 "sggevx.f"
    if (ierr != 0) {
#line 684 "sggevx.f"
	if (ierr > 0 && ierr <= *n) {
#line 685 "sggevx.f"
	    *info = ierr;
#line 686 "sggevx.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 687 "sggevx.f"
	    *info = ierr - *n;
#line 688 "sggevx.f"
	} else {
#line 689 "sggevx.f"
	    *info = *n + 1;
#line 690 "sggevx.f"
	}
#line 691 "sggevx.f"
	goto L130;
#line 692 "sggevx.f"
    }

/*     Compute Eigenvectors and estimate condition numbers if desired */
/*     (Workspace: STGEVC: need 6*N */
/*                 STGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B', */
/*                         need N otherwise ) */

#line 699 "sggevx.f"
    if (ilv || ! wantsn) {
#line 700 "sggevx.f"
	if (ilv) {
#line 701 "sggevx.f"
	    if (ilvl) {
#line 702 "sggevx.f"
		if (ilvr) {
#line 703 "sggevx.f"
		    *(unsigned char *)chtemp = 'B';
#line 704 "sggevx.f"
		} else {
#line 705 "sggevx.f"
		    *(unsigned char *)chtemp = 'L';
#line 706 "sggevx.f"
		}
#line 707 "sggevx.f"
	    } else {
#line 708 "sggevx.f"
		*(unsigned char *)chtemp = 'R';
#line 709 "sggevx.f"
	    }

#line 711 "sggevx.f"
	    stgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], 
		    ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &
		    work[1], &ierr, (ftnlen)1, (ftnlen)1);
#line 713 "sggevx.f"
	    if (ierr != 0) {
#line 714 "sggevx.f"
		*info = *n + 2;
#line 715 "sggevx.f"
		goto L130;
#line 716 "sggevx.f"
	    }
#line 717 "sggevx.f"
	}

#line 719 "sggevx.f"
	if (! wantsn) {

/*           compute eigenvectors (STGEVC) and estimate condition */
/*           numbers (STGSNA). Note that the definition of the condition */
/*           number is not invariant under transformation (u,v) to */
/*           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized */
/*           Schur form (S,T), Q and Z are orthogonal matrices. In order */
/*           to avoid using extra 2*N*N workspace, we have to recalculate */
/*           eigenvectors and estimate one condition numbers at a time. */

#line 729 "sggevx.f"
	    pair = FALSE_;
#line 730 "sggevx.f"
	    i__1 = *n;
#line 730 "sggevx.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

#line 732 "sggevx.f"
		if (pair) {
#line 733 "sggevx.f"
		    pair = FALSE_;
#line 734 "sggevx.f"
		    goto L20;
#line 735 "sggevx.f"
		}
#line 736 "sggevx.f"
		mm = 1;
#line 737 "sggevx.f"
		if (i__ < *n) {
#line 738 "sggevx.f"
		    if (a[i__ + 1 + i__ * a_dim1] != 0.) {
#line 739 "sggevx.f"
			pair = TRUE_;
#line 740 "sggevx.f"
			mm = 2;
#line 741 "sggevx.f"
		    }
#line 742 "sggevx.f"
		}

#line 744 "sggevx.f"
		i__2 = *n;
#line 744 "sggevx.f"
		for (j = 1; j <= i__2; ++j) {
#line 745 "sggevx.f"
		    bwork[j] = FALSE_;
#line 746 "sggevx.f"
/* L10: */
#line 746 "sggevx.f"
		}
#line 747 "sggevx.f"
		if (mm == 1) {
#line 748 "sggevx.f"
		    bwork[i__] = TRUE_;
#line 749 "sggevx.f"
		} else if (mm == 2) {
#line 750 "sggevx.f"
		    bwork[i__] = TRUE_;
#line 751 "sggevx.f"
		    bwork[i__ + 1] = TRUE_;
#line 752 "sggevx.f"
		}

#line 754 "sggevx.f"
		iwrk = mm * *n + 1;
#line 755 "sggevx.f"
		iwrk1 = iwrk + mm * *n;

/*              Compute a pair of left and right eigenvectors. */
/*              (compute workspace: need up to 4*N + 6*N) */

#line 760 "sggevx.f"
		if (wantse || wantsb) {
#line 761 "sggevx.f"
		    stgevc_("B", "S", &bwork[1], n, &a[a_offset], lda, &b[
			    b_offset], ldb, &work[1], n, &work[iwrk], n, &mm, 
			    &m, &work[iwrk1], &ierr, (ftnlen)1, (ftnlen)1);
#line 764 "sggevx.f"
		    if (ierr != 0) {
#line 765 "sggevx.f"
			*info = *n + 2;
#line 766 "sggevx.f"
			goto L130;
#line 767 "sggevx.f"
		    }
#line 768 "sggevx.f"
		}

#line 770 "sggevx.f"
		i__2 = *lwork - iwrk1 + 1;
#line 770 "sggevx.f"
		stgsna_(sense, "S", &bwork[1], n, &a[a_offset], lda, &b[
			b_offset], ldb, &work[1], n, &work[iwrk], n, &rconde[
			i__], &rcondv[i__], &mm, &m, &work[iwrk1], &i__2, &
			iwork[1], &ierr, (ftnlen)1, (ftnlen)1);

#line 775 "sggevx.f"
L20:
#line 775 "sggevx.f"
		;
#line 775 "sggevx.f"
	    }
#line 776 "sggevx.f"
	}
#line 777 "sggevx.f"
    }

/*     Undo balancing on VL and VR and normalization */
/*     (Workspace: none needed) */

#line 782 "sggevx.f"
    if (ilvl) {
#line 783 "sggevx.f"
	sggbak_(balanc, "L", n, ilo, ihi, &lscale[1], &rscale[1], n, &vl[
		vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);

#line 786 "sggevx.f"
	i__1 = *n;
#line 786 "sggevx.f"
	for (jc = 1; jc <= i__1; ++jc) {
#line 787 "sggevx.f"
	    if (alphai[jc] < 0.) {
#line 787 "sggevx.f"
		goto L70;
#line 787 "sggevx.f"
	    }
#line 789 "sggevx.f"
	    temp = 0.;
#line 790 "sggevx.f"
	    if (alphai[jc] == 0.) {
#line 791 "sggevx.f"
		i__2 = *n;
#line 791 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 792 "sggevx.f"
		    d__2 = temp, d__3 = (d__1 = vl[jr + jc * vl_dim1], abs(
			    d__1));
#line 792 "sggevx.f"
		    temp = max(d__2,d__3);
#line 793 "sggevx.f"
/* L30: */
#line 793 "sggevx.f"
		}
#line 794 "sggevx.f"
	    } else {
#line 795 "sggevx.f"
		i__2 = *n;
#line 795 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 796 "sggevx.f"
		    d__3 = temp, d__4 = (d__1 = vl[jr + jc * vl_dim1], abs(
			    d__1)) + (d__2 = vl[jr + (jc + 1) * vl_dim1], abs(
			    d__2));
#line 796 "sggevx.f"
		    temp = max(d__3,d__4);
#line 798 "sggevx.f"
/* L40: */
#line 798 "sggevx.f"
		}
#line 799 "sggevx.f"
	    }
#line 800 "sggevx.f"
	    if (temp < smlnum) {
#line 800 "sggevx.f"
		goto L70;
#line 800 "sggevx.f"
	    }
#line 802 "sggevx.f"
	    temp = 1. / temp;
#line 803 "sggevx.f"
	    if (alphai[jc] == 0.) {
#line 804 "sggevx.f"
		i__2 = *n;
#line 804 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 805 "sggevx.f"
		    vl[jr + jc * vl_dim1] *= temp;
#line 806 "sggevx.f"
/* L50: */
#line 806 "sggevx.f"
		}
#line 807 "sggevx.f"
	    } else {
#line 808 "sggevx.f"
		i__2 = *n;
#line 808 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 809 "sggevx.f"
		    vl[jr + jc * vl_dim1] *= temp;
#line 810 "sggevx.f"
		    vl[jr + (jc + 1) * vl_dim1] *= temp;
#line 811 "sggevx.f"
/* L60: */
#line 811 "sggevx.f"
		}
#line 812 "sggevx.f"
	    }
#line 813 "sggevx.f"
L70:
#line 813 "sggevx.f"
	    ;
#line 813 "sggevx.f"
	}
#line 814 "sggevx.f"
    }
#line 815 "sggevx.f"
    if (ilvr) {
#line 816 "sggevx.f"
	sggbak_(balanc, "R", n, ilo, ihi, &lscale[1], &rscale[1], n, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 818 "sggevx.f"
	i__1 = *n;
#line 818 "sggevx.f"
	for (jc = 1; jc <= i__1; ++jc) {
#line 819 "sggevx.f"
	    if (alphai[jc] < 0.) {
#line 819 "sggevx.f"
		goto L120;
#line 819 "sggevx.f"
	    }
#line 821 "sggevx.f"
	    temp = 0.;
#line 822 "sggevx.f"
	    if (alphai[jc] == 0.) {
#line 823 "sggevx.f"
		i__2 = *n;
#line 823 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 824 "sggevx.f"
		    d__2 = temp, d__3 = (d__1 = vr[jr + jc * vr_dim1], abs(
			    d__1));
#line 824 "sggevx.f"
		    temp = max(d__2,d__3);
#line 825 "sggevx.f"
/* L80: */
#line 825 "sggevx.f"
		}
#line 826 "sggevx.f"
	    } else {
#line 827 "sggevx.f"
		i__2 = *n;
#line 827 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 828 "sggevx.f"
		    d__3 = temp, d__4 = (d__1 = vr[jr + jc * vr_dim1], abs(
			    d__1)) + (d__2 = vr[jr + (jc + 1) * vr_dim1], abs(
			    d__2));
#line 828 "sggevx.f"
		    temp = max(d__3,d__4);
#line 830 "sggevx.f"
/* L90: */
#line 830 "sggevx.f"
		}
#line 831 "sggevx.f"
	    }
#line 832 "sggevx.f"
	    if (temp < smlnum) {
#line 832 "sggevx.f"
		goto L120;
#line 832 "sggevx.f"
	    }
#line 834 "sggevx.f"
	    temp = 1. / temp;
#line 835 "sggevx.f"
	    if (alphai[jc] == 0.) {
#line 836 "sggevx.f"
		i__2 = *n;
#line 836 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 837 "sggevx.f"
		    vr[jr + jc * vr_dim1] *= temp;
#line 838 "sggevx.f"
/* L100: */
#line 838 "sggevx.f"
		}
#line 839 "sggevx.f"
	    } else {
#line 840 "sggevx.f"
		i__2 = *n;
#line 840 "sggevx.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 841 "sggevx.f"
		    vr[jr + jc * vr_dim1] *= temp;
#line 842 "sggevx.f"
		    vr[jr + (jc + 1) * vr_dim1] *= temp;
#line 843 "sggevx.f"
/* L110: */
#line 843 "sggevx.f"
		}
#line 844 "sggevx.f"
	    }
#line 845 "sggevx.f"
L120:
#line 845 "sggevx.f"
	    ;
#line 845 "sggevx.f"
	}
#line 846 "sggevx.f"
    }

/*     Undo scaling if necessary */

#line 850 "sggevx.f"
L130:

#line 852 "sggevx.f"
    if (ilascl) {
#line 853 "sggevx.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 854 "sggevx.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 855 "sggevx.f"
    }

#line 857 "sggevx.f"
    if (ilbscl) {
#line 858 "sggevx.f"
	slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 859 "sggevx.f"
    }

#line 861 "sggevx.f"
    work[1] = (doublereal) maxwrk;
#line 862 "sggevx.f"
    return 0;

/*     End of SGGEVX */

} /* sggevx_ */

