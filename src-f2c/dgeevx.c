#line 1 "dgeevx.f"
/* dgeevx.f -- translated by f2c (version 20100827).
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

#line 1 "dgeevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, */
/*                          VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, */
/*                          RCONDE, RCONDV, WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE */
/*       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N */
/*       DOUBLE PRECISION   ABNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ), */
/*      $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WI( * ), WORK( * ), WR( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEEVX computes for an N-by-N real nonsymmetric matrix A, the */
/* > eigenvalues and, optionally, the left and/or right eigenvectors. */
/* > */
/* > Optionally also, it computes a balancing transformation to improve */
/* > the conditioning of the eigenvalues and eigenvectors (ILO, IHI, */
/* > SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues */
/* > (RCONDE), and reciprocal condition numbers for the right */
/* > eigenvectors (RCONDV). */
/* > */
/* > The right eigenvector v(j) of A satisfies */
/* >                  A * v(j) = lambda(j) * v(j) */
/* > where lambda(j) is its eigenvalue. */
/* > The left eigenvector u(j) of A satisfies */
/* >               u(j)**H * A = lambda(j) * u(j)**H */
/* > where u(j)**H denotes the conjugate-transpose of u(j). */
/* > */
/* > The computed eigenvectors are normalized to have Euclidean norm */
/* > equal to 1 and largest component real. */
/* > */
/* > Balancing a matrix means permuting the rows and columns to make it */
/* > more nearly upper triangular, and applying a diagonal similarity */
/* > transformation D * A * D**(-1), where D is a diagonal matrix, to */
/* > make its rows and columns closer in norm and the condition numbers */
/* > of its eigenvalues and eigenvectors smaller.  The computed */
/* > reciprocal condition numbers correspond to the balanced matrix. */
/* > Permuting rows and columns will not change the condition numbers */
/* > (in exact arithmetic) but diagonal scaling will.  For further */
/* > explanation of balancing, see section 4.10.2 of the LAPACK */
/* > Users' Guide. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] BALANC */
/* > \verbatim */
/* >          BALANC is CHARACTER*1 */
/* >          Indicates how the input matrix should be diagonally scaled */
/* >          and/or permuted to improve the conditioning of its */
/* >          eigenvalues. */
/* >          = 'N': Do not diagonally scale or permute; */
/* >          = 'P': Perform permutations to make the matrix more nearly */
/* >                 upper triangular. Do not diagonally scale; */
/* >          = 'S': Diagonally scale the matrix, i.e. replace A by */
/* >                 D*A*D**(-1), where D is a diagonal matrix chosen */
/* >                 to make the rows and columns of A more equal in */
/* >                 norm. Do not permute; */
/* >          = 'B': Both diagonally scale and permute A. */
/* > */
/* >          Computed reciprocal condition numbers will be for the matrix */
/* >          after balancing and/or permuting. Permuting does not change */
/* >          condition numbers (in exact arithmetic), but balancing does. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVL */
/* > \verbatim */
/* >          JOBVL is CHARACTER*1 */
/* >          = 'N': left eigenvectors of A are not computed; */
/* >          = 'V': left eigenvectors of A are computed. */
/* >          If SENSE = 'E' or 'B', JOBVL must = 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* >          JOBVR is CHARACTER*1 */
/* >          = 'N': right eigenvectors of A are not computed; */
/* >          = 'V': right eigenvectors of A are computed. */
/* >          If SENSE = 'E' or 'B', JOBVR must = 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] SENSE */
/* > \verbatim */
/* >          SENSE is CHARACTER*1 */
/* >          Determines which reciprocal condition numbers are computed. */
/* >          = 'N': None are computed; */
/* >          = 'E': Computed for eigenvalues only; */
/* >          = 'V': Computed for right eigenvectors only; */
/* >          = 'B': Computed for eigenvalues and right eigenvectors. */
/* > */
/* >          If SENSE = 'E' or 'B', both left and right eigenvectors */
/* >          must also be computed (JOBVL = 'V' and JOBVR = 'V'). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the N-by-N matrix A. */
/* >          On exit, A has been overwritten.  If JOBVL = 'V' or */
/* >          JOBVR = 'V', A contains the real Schur form of the balanced */
/* >          version of the input matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WR */
/* > \verbatim */
/* >          WR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is DOUBLE PRECISION array, dimension (N) */
/* >          WR and WI contain the real and imaginary parts, */
/* >          respectively, of the computed eigenvalues.  Complex */
/* >          conjugate pairs of eigenvalues will appear consecutively */
/* >          with the eigenvalue having the positive imaginary part */
/* >          first. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* >          after another in the columns of VL, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVL = 'N', VL is not referenced. */
/* >          If the j-th eigenvalue is real, then u(j) = VL(:,j), */
/* >          the j-th column of VL. */
/* >          If the j-th and (j+1)-st eigenvalues form a complex */
/* >          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and */
/* >          u(j+1) = VL(:,j) - i*VL(:,j+1). */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL.  LDVL >= 1; if */
/* >          JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* >          VR is DOUBLE PRECISION array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* >          after another in the columns of VR, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVR = 'N', VR is not referenced. */
/* >          If the j-th eigenvalue is real, then v(j) = VR(:,j), */
/* >          the j-th column of VR. */
/* >          If the j-th and (j+1)-st eigenvalues form a complex */
/* >          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and */
/* >          v(j+1) = VR(:,j) - i*VR(:,j+1). */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR.  LDVR >= 1, and if */
/* >          JOBVR = 'V', LDVR >= N. */
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
/* >          ILO and IHI are integer values determined when A was */
/* >          balanced.  The balanced A(i,j) = 0 if I > J and */
/* >          J = 1,...,ILO-1 or I = IHI+1,...,N. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and scaling factors applied */
/* >          when balancing A.  If P(j) is the index of the row and column */
/* >          interchanged with row and column j, and D(j) is the scaling */
/* >          factor applied to row and column j, then */
/* >          SCALE(J) = P(J),    for J = 1,...,ILO-1 */
/* >                   = D(J),    for J = ILO,...,IHI */
/* >                   = P(J)     for J = IHI+1,...,N. */
/* >          The order in which the interchanges are made is N to IHI+1, */
/* >          then 1 to ILO-1. */
/* > \endverbatim */
/* > */
/* > \param[out] ABNRM */
/* > \verbatim */
/* >          ABNRM is DOUBLE PRECISION */
/* >          The one-norm of the balanced matrix (the maximum */
/* >          of the sum of absolute values of elements of any column). */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDE */
/* > \verbatim */
/* >          RCONDE is DOUBLE PRECISION array, dimension (N) */
/* >          RCONDE(j) is the reciprocal condition number of the j-th */
/* >          eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is DOUBLE PRECISION array, dimension (N) */
/* >          RCONDV(j) is the reciprocal condition number of the j-th */
/* >          right eigenvector. */
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
/* >          The dimension of the array WORK.   If SENSE = 'N' or 'E', */
/* >          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V', */
/* >          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6). */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (2*N-2) */
/* >          If SENSE = 'N' or 'E', not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the QR algorithm failed to compute all the */
/* >                eigenvalues, and no eigenvectors or condition numbers */
/* >                have been computed; elements 1:ILO-1 and i+1:N of WR */
/* >                and WI contain eigenvalues which have converged. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *wr, 
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	integer *ldvr, integer *ilo, integer *ihi, doublereal *scale, 
	doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal 
	*work, integer *lwork, integer *iwork, integer *info, ftnlen 
	balanc_len, ftnlen jobvl_len, ftnlen jobvr_len, ftnlen sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k;
    static doublereal r__, cs, sn;
    static char job[1];
    static doublereal scl, dum[1], eps;
    static char side[1];
    static doublereal anrm;
    static integer ierr, itau;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer iwrk, nout;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer icond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dgebal_(char *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static logical scalea;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen);
    static logical select[1];
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dorghr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dhseqr_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dtrevc_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen), dtrsna_(char *, char *, logical *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static logical wantvl, wntsnb;
    static integer hswork;
    static logical wntsne;
    static doublereal smlnum;
    static logical lquery, wantvr, wntsnn, wntsnv;


/*  -- LAPACK driver routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

/*     Test the input arguments */

#line 361 "dgeevx.f"
    /* Parameter adjustments */
#line 361 "dgeevx.f"
    a_dim1 = *lda;
#line 361 "dgeevx.f"
    a_offset = 1 + a_dim1;
#line 361 "dgeevx.f"
    a -= a_offset;
#line 361 "dgeevx.f"
    --wr;
#line 361 "dgeevx.f"
    --wi;
#line 361 "dgeevx.f"
    vl_dim1 = *ldvl;
#line 361 "dgeevx.f"
    vl_offset = 1 + vl_dim1;
#line 361 "dgeevx.f"
    vl -= vl_offset;
#line 361 "dgeevx.f"
    vr_dim1 = *ldvr;
#line 361 "dgeevx.f"
    vr_offset = 1 + vr_dim1;
#line 361 "dgeevx.f"
    vr -= vr_offset;
#line 361 "dgeevx.f"
    --scale;
#line 361 "dgeevx.f"
    --rconde;
#line 361 "dgeevx.f"
    --rcondv;
#line 361 "dgeevx.f"
    --work;
#line 361 "dgeevx.f"
    --iwork;
#line 361 "dgeevx.f"

#line 361 "dgeevx.f"
    /* Function Body */
#line 361 "dgeevx.f"
    *info = 0;
#line 362 "dgeevx.f"
    lquery = *lwork == -1;
#line 363 "dgeevx.f"
    wantvl = lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1);
#line 364 "dgeevx.f"
    wantvr = lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1);
#line 365 "dgeevx.f"
    wntsnn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 366 "dgeevx.f"
    wntsne = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 367 "dgeevx.f"
    wntsnv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 368 "dgeevx.f"
    wntsnb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 369 "dgeevx.f"
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "S", (
	    ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) 
	    || lsame_(balanc, "B", (ftnlen)1, (ftnlen)1))) {
#line 372 "dgeevx.f"
	*info = -1;
#line 373 "dgeevx.f"
    } else if (! wantvl && ! lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 374 "dgeevx.f"
	*info = -2;
#line 375 "dgeevx.f"
    } else if (! wantvr && ! lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 376 "dgeevx.f"
	*info = -3;
#line 377 "dgeevx.f"
    } else if (! (wntsnn || wntsne || wntsnb || wntsnv) || (wntsne || wntsnb) 
	    && ! (wantvl && wantvr)) {
#line 380 "dgeevx.f"
	*info = -4;
#line 381 "dgeevx.f"
    } else if (*n < 0) {
#line 382 "dgeevx.f"
	*info = -5;
#line 383 "dgeevx.f"
    } else if (*lda < max(1,*n)) {
#line 384 "dgeevx.f"
	*info = -7;
#line 385 "dgeevx.f"
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
#line 386 "dgeevx.f"
	*info = -11;
#line 387 "dgeevx.f"
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
#line 388 "dgeevx.f"
	*info = -13;
#line 389 "dgeevx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by DHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 401 "dgeevx.f"
    if (*info == 0) {
#line 402 "dgeevx.f"
	if (*n == 0) {
#line 403 "dgeevx.f"
	    minwrk = 1;
#line 404 "dgeevx.f"
	    maxwrk = 1;
#line 405 "dgeevx.f"
	} else {
#line 406 "dgeevx.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "DGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);

#line 408 "dgeevx.f"
	    if (wantvl) {
#line 409 "dgeevx.f"
		dhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vl[vl_offset], ldvl, &work[1], &c_n1, info, (
			ftnlen)1, (ftnlen)1);
#line 411 "dgeevx.f"
	    } else if (wantvr) {
#line 412 "dgeevx.f"
		dhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			ftnlen)1, (ftnlen)1);
#line 414 "dgeevx.f"
	    } else {
#line 415 "dgeevx.f"
		if (wntsnn) {
#line 416 "dgeevx.f"
		    dhseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &wr[1], 
			    &wi[1], &vr[vr_offset], ldvr, &work[1], &c_n1, 
			    info, (ftnlen)1, (ftnlen)1);
#line 418 "dgeevx.f"
		} else {
#line 419 "dgeevx.f"
		    dhseqr_("S", "N", n, &c__1, n, &a[a_offset], lda, &wr[1], 
			    &wi[1], &vr[vr_offset], ldvr, &work[1], &c_n1, 
			    info, (ftnlen)1, (ftnlen)1);
#line 421 "dgeevx.f"
		}
#line 422 "dgeevx.f"
	    }
#line 423 "dgeevx.f"
	    hswork = (integer) work[1];

#line 425 "dgeevx.f"
	    if (! wantvl && ! wantvr) {
#line 426 "dgeevx.f"
		minwrk = *n << 1;
#line 427 "dgeevx.f"
		if (! wntsnn) {
/* Computing MAX */
#line 427 "dgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + *n * 6;
#line 427 "dgeevx.f"
		    minwrk = max(i__1,i__2);
#line 427 "dgeevx.f"
		}
#line 429 "dgeevx.f"
		maxwrk = max(maxwrk,hswork);
#line 430 "dgeevx.f"
		if (! wntsnn) {
/* Computing MAX */
#line 430 "dgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + *n * 6;
#line 430 "dgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 430 "dgeevx.f"
		}
#line 432 "dgeevx.f"
	    } else {
#line 433 "dgeevx.f"
		minwrk = *n * 3;
#line 434 "dgeevx.f"
		if (! wntsnn && ! wntsne) {
/* Computing MAX */
#line 434 "dgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + *n * 6;
#line 434 "dgeevx.f"
		    minwrk = max(i__1,i__2);
#line 434 "dgeevx.f"
		}
#line 436 "dgeevx.f"
		maxwrk = max(maxwrk,hswork);
/* Computing MAX */
#line 437 "dgeevx.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "DORGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 437 "dgeevx.f"
		maxwrk = max(i__1,i__2);
#line 439 "dgeevx.f"
		if (! wntsnn && ! wntsne) {
/* Computing MAX */
#line 439 "dgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + *n * 6;
#line 439 "dgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 439 "dgeevx.f"
		}
/* Computing MAX */
#line 441 "dgeevx.f"
		i__1 = maxwrk, i__2 = *n * 3;
#line 441 "dgeevx.f"
		maxwrk = max(i__1,i__2);
#line 442 "dgeevx.f"
	    }
#line 443 "dgeevx.f"
	    maxwrk = max(maxwrk,minwrk);
#line 444 "dgeevx.f"
	}
#line 445 "dgeevx.f"
	work[1] = (doublereal) maxwrk;

#line 447 "dgeevx.f"
	if (*lwork < minwrk && ! lquery) {
#line 448 "dgeevx.f"
	    *info = -21;
#line 449 "dgeevx.f"
	}
#line 450 "dgeevx.f"
    }

#line 452 "dgeevx.f"
    if (*info != 0) {
#line 453 "dgeevx.f"
	i__1 = -(*info);
#line 453 "dgeevx.f"
	xerbla_("DGEEVX", &i__1, (ftnlen)6);
#line 454 "dgeevx.f"
	return 0;
#line 455 "dgeevx.f"
    } else if (lquery) {
#line 456 "dgeevx.f"
	return 0;
#line 457 "dgeevx.f"
    }

/*     Quick return if possible */

#line 461 "dgeevx.f"
    if (*n == 0) {
#line 461 "dgeevx.f"
	return 0;
#line 461 "dgeevx.f"
    }

/*     Get machine constants */

#line 466 "dgeevx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 467 "dgeevx.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 468 "dgeevx.f"
    bignum = 1. / smlnum;
#line 469 "dgeevx.f"
    dlabad_(&smlnum, &bignum);
#line 470 "dgeevx.f"
    smlnum = sqrt(smlnum) / eps;
#line 471 "dgeevx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 475 "dgeevx.f"
    icond = 0;
#line 476 "dgeevx.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 477 "dgeevx.f"
    scalea = FALSE_;
#line 478 "dgeevx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 479 "dgeevx.f"
	scalea = TRUE_;
#line 480 "dgeevx.f"
	cscale = smlnum;
#line 481 "dgeevx.f"
    } else if (anrm > bignum) {
#line 482 "dgeevx.f"
	scalea = TRUE_;
#line 483 "dgeevx.f"
	cscale = bignum;
#line 484 "dgeevx.f"
    }
#line 485 "dgeevx.f"
    if (scalea) {
#line 485 "dgeevx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 485 "dgeevx.f"
    }

/*     Balance the matrix and compute ABNRM */

#line 490 "dgeevx.f"
    dgebal_(balanc, n, &a[a_offset], lda, ilo, ihi, &scale[1], &ierr, (ftnlen)
	    1);
#line 491 "dgeevx.f"
    *abnrm = dlange_("1", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 492 "dgeevx.f"
    if (scalea) {
#line 493 "dgeevx.f"
	dum[0] = *abnrm;
#line 494 "dgeevx.f"
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &c__1, &
		ierr, (ftnlen)1);
#line 495 "dgeevx.f"
	*abnrm = dum[0];
#line 496 "dgeevx.f"
    }

/*     Reduce to upper Hessenberg form */
/*     (Workspace: need 2*N, prefer N+N*NB) */

#line 501 "dgeevx.f"
    itau = 1;
#line 502 "dgeevx.f"
    iwrk = itau + *n;
#line 503 "dgeevx.f"
    i__1 = *lwork - iwrk + 1;
#line 503 "dgeevx.f"
    dgehrd_(n, ilo, ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &
	    ierr);

#line 506 "dgeevx.f"
    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

#line 511 "dgeevx.f"
	*(unsigned char *)side = 'L';
#line 512 "dgeevx.f"
	dlacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VL */
/*        (Workspace: need 2*N-1, prefer N+(N-1)*NB) */

#line 517 "dgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 517 "dgeevx.f"
	dorghr_(n, ilo, ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

#line 523 "dgeevx.f"
	iwrk = itau;
#line 524 "dgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 524 "dgeevx.f"
	dhseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vl[
		vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 527 "dgeevx.f"
	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

#line 532 "dgeevx.f"
	    *(unsigned char *)side = 'B';
#line 533 "dgeevx.f"
	    dlacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (
		    ftnlen)1);
#line 534 "dgeevx.f"
	}

#line 536 "dgeevx.f"
    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

#line 541 "dgeevx.f"
	*(unsigned char *)side = 'R';
#line 542 "dgeevx.f"
	dlacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VR */
/*        (Workspace: need 2*N-1, prefer N+(N-1)*NB) */

#line 547 "dgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 547 "dgeevx.f"
	dorghr_(n, ilo, ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

#line 553 "dgeevx.f"
	iwrk = itau;
#line 554 "dgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 554 "dgeevx.f"
	dhseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 557 "dgeevx.f"
    } else {

/*        Compute eigenvalues only */
/*        If condition numbers desired, compute Schur form */

#line 562 "dgeevx.f"
	if (wntsnn) {
#line 563 "dgeevx.f"
	    *(unsigned char *)job = 'E';
#line 564 "dgeevx.f"
	} else {
#line 565 "dgeevx.f"
	    *(unsigned char *)job = 'S';
#line 566 "dgeevx.f"
	}

/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

#line 570 "dgeevx.f"
	iwrk = itau;
#line 571 "dgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 571 "dgeevx.f"
	dhseqr_(job, "N", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);
#line 573 "dgeevx.f"
    }

/*     If INFO > 0 from DHSEQR, then quit */

#line 577 "dgeevx.f"
    if (*info > 0) {
#line 577 "dgeevx.f"
	goto L50;
#line 577 "dgeevx.f"
    }

#line 580 "dgeevx.f"
    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (Workspace: need 3*N) */

#line 585 "dgeevx.f"
	dtrevc_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
		 &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &ierr, (ftnlen)
		1, (ftnlen)1);
#line 587 "dgeevx.f"
    }

/*     Compute condition numbers if desired */
/*     (Workspace: need N*N+6*N unless SENSE = 'E') */

#line 592 "dgeevx.f"
    if (! wntsnn) {
#line 593 "dgeevx.f"
	dtrsna_(sense, "A", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, &rconde[1], &rcondv[1], n, &nout, 
		&work[iwrk], n, &iwork[1], &icond, (ftnlen)1, (ftnlen)1);
#line 596 "dgeevx.f"
    }

#line 598 "dgeevx.f"
    if (wantvl) {

/*        Undo balancing of left eigenvectors */

#line 602 "dgeevx.f"
	dgebak_(balanc, "L", n, ilo, ihi, &scale[1], n, &vl[vl_offset], ldvl, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize left eigenvectors and make largest component real */

#line 607 "dgeevx.f"
	i__1 = *n;
#line 607 "dgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 608 "dgeevx.f"
	    if (wi[i__] == 0.) {
#line 609 "dgeevx.f"
		scl = 1. / dnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 610 "dgeevx.f"
		dscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 611 "dgeevx.f"
	    } else if (wi[i__] > 0.) {
#line 612 "dgeevx.f"
		d__1 = dnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 612 "dgeevx.f"
		d__2 = dnrm2_(n, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
#line 612 "dgeevx.f"
		scl = 1. / dlapy2_(&d__1, &d__2);
#line 614 "dgeevx.f"
		dscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 615 "dgeevx.f"
		dscal_(n, &scl, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
#line 616 "dgeevx.f"
		i__2 = *n;
#line 616 "dgeevx.f"
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
#line 617 "dgeevx.f"
		    d__1 = vl[k + i__ * vl_dim1];
/* Computing 2nd power */
#line 617 "dgeevx.f"
		    d__2 = vl[k + (i__ + 1) * vl_dim1];
#line 617 "dgeevx.f"
		    work[k] = d__1 * d__1 + d__2 * d__2;
#line 618 "dgeevx.f"
/* L10: */
#line 618 "dgeevx.f"
		}
#line 619 "dgeevx.f"
		k = idamax_(n, &work[1], &c__1);
#line 620 "dgeevx.f"
		dlartg_(&vl[k + i__ * vl_dim1], &vl[k + (i__ + 1) * vl_dim1], 
			&cs, &sn, &r__);
#line 621 "dgeevx.f"
		drot_(n, &vl[i__ * vl_dim1 + 1], &c__1, &vl[(i__ + 1) * 
			vl_dim1 + 1], &c__1, &cs, &sn);
#line 622 "dgeevx.f"
		vl[k + (i__ + 1) * vl_dim1] = 0.;
#line 623 "dgeevx.f"
	    }
#line 624 "dgeevx.f"
/* L20: */
#line 624 "dgeevx.f"
	}
#line 625 "dgeevx.f"
    }

#line 627 "dgeevx.f"
    if (wantvr) {

/*        Undo balancing of right eigenvectors */

#line 631 "dgeevx.f"
	dgebak_(balanc, "R", n, ilo, ihi, &scale[1], n, &vr[vr_offset], ldvr, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize right eigenvectors and make largest component real */

#line 636 "dgeevx.f"
	i__1 = *n;
#line 636 "dgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 637 "dgeevx.f"
	    if (wi[i__] == 0.) {
#line 638 "dgeevx.f"
		scl = 1. / dnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 639 "dgeevx.f"
		dscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 640 "dgeevx.f"
	    } else if (wi[i__] > 0.) {
#line 641 "dgeevx.f"
		d__1 = dnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 641 "dgeevx.f"
		d__2 = dnrm2_(n, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
#line 641 "dgeevx.f"
		scl = 1. / dlapy2_(&d__1, &d__2);
#line 643 "dgeevx.f"
		dscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 644 "dgeevx.f"
		dscal_(n, &scl, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
#line 645 "dgeevx.f"
		i__2 = *n;
#line 645 "dgeevx.f"
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
#line 646 "dgeevx.f"
		    d__1 = vr[k + i__ * vr_dim1];
/* Computing 2nd power */
#line 646 "dgeevx.f"
		    d__2 = vr[k + (i__ + 1) * vr_dim1];
#line 646 "dgeevx.f"
		    work[k] = d__1 * d__1 + d__2 * d__2;
#line 647 "dgeevx.f"
/* L30: */
#line 647 "dgeevx.f"
		}
#line 648 "dgeevx.f"
		k = idamax_(n, &work[1], &c__1);
#line 649 "dgeevx.f"
		dlartg_(&vr[k + i__ * vr_dim1], &vr[k + (i__ + 1) * vr_dim1], 
			&cs, &sn, &r__);
#line 650 "dgeevx.f"
		drot_(n, &vr[i__ * vr_dim1 + 1], &c__1, &vr[(i__ + 1) * 
			vr_dim1 + 1], &c__1, &cs, &sn);
#line 651 "dgeevx.f"
		vr[k + (i__ + 1) * vr_dim1] = 0.;
#line 652 "dgeevx.f"
	    }
#line 653 "dgeevx.f"
/* L40: */
#line 653 "dgeevx.f"
	}
#line 654 "dgeevx.f"
    }

/*     Undo scaling if necessary */

#line 658 "dgeevx.f"
L50:
#line 659 "dgeevx.f"
    if (scalea) {
#line 660 "dgeevx.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 660 "dgeevx.f"
	i__3 = *n - *info;
#line 660 "dgeevx.f"
	i__2 = max(i__3,1);
#line 660 "dgeevx.f"
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 
		1], &i__2, &ierr, (ftnlen)1);
#line 662 "dgeevx.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 662 "dgeevx.f"
	i__3 = *n - *info;
#line 662 "dgeevx.f"
	i__2 = max(i__3,1);
#line 662 "dgeevx.f"
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 
		1], &i__2, &ierr, (ftnlen)1);
#line 664 "dgeevx.f"
	if (*info == 0) {
#line 665 "dgeevx.f"
	    if ((wntsnv || wntsnb) && icond == 0) {
#line 665 "dgeevx.f"
		dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &rcondv[
			1], n, &ierr, (ftnlen)1);
#line 665 "dgeevx.f"
	    }
#line 668 "dgeevx.f"
	} else {
#line 669 "dgeevx.f"
	    i__1 = *ilo - 1;
#line 669 "dgeevx.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], 
		    n, &ierr, (ftnlen)1);
#line 671 "dgeevx.f"
	    i__1 = *ilo - 1;
#line 671 "dgeevx.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], 
		    n, &ierr, (ftnlen)1);
#line 673 "dgeevx.f"
	}
#line 674 "dgeevx.f"
    }

#line 676 "dgeevx.f"
    work[1] = (doublereal) maxwrk;
#line 677 "dgeevx.f"
    return 0;

/*     End of DGEEVX */

} /* dgeevx_ */

