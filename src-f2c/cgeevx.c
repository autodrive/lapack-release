#line 1 "cgeevx.f"
/* cgeevx.f -- translated by f2c (version 20100827).
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

#line 1 "cgeevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> CGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, */
/*                          LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, */
/*                          RCONDV, WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE */
/*       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N */
/*       REAL               ABNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RCONDE( * ), RCONDV( * ), RWORK( * ), */
/*      $                   SCALE( * ) */
/*       COMPLEX            A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEEVX computes for an N-by-N complex nonsymmetric matrix A, the */
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
/* > where u(j)**H denotes the conjugate transpose of u(j). */
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
/* >          = 'S': Diagonally scale the matrix, ie. replace A by */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the N-by-N matrix A. */
/* >          On exit, A has been overwritten.  If JOBVL = 'V' or */
/* >          JOBVR = 'V', A contains the Schur form of the balanced */
/* >          version of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX array, dimension (N) */
/* >          W contains the computed eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is COMPLEX array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* >          after another in the columns of VL, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVL = 'N', VL is not referenced. */
/* >          u(j) = VL(:,j), the j-th column of VL. */
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
/* >          VR is COMPLEX array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* >          after another in the columns of VR, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVR = 'N', VR is not referenced. */
/* >          v(j) = VR(:,j), the j-th column of VR. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR.  LDVR >= 1; if */
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
/* >          SCALE is REAL array, dimension (N) */
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
/* >          ABNRM is REAL */
/* >          The one-norm of the balanced matrix (the maximum */
/* >          of the sum of absolute values of elements of any column). */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDE */
/* > \verbatim */
/* >          RCONDE is REAL array, dimension (N) */
/* >          RCONDE(j) is the reciprocal condition number of the j-th */
/* >          eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is REAL array, dimension (N) */
/* >          RCONDV(j) is the reciprocal condition number of the j-th */
/* >          right eigenvector. */
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
/* >          The dimension of the array WORK.  If SENSE = 'N' or 'E', */
/* >          LWORK >= max(1,2*N), and if SENSE = 'V' or 'B', */
/* >          LWORK >= N*N+2*N. */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the QR algorithm failed to compute all the */
/* >                eigenvalues, and no eigenvectors or condition numbers */
/* >                have been computed; elements 1:ILO-1 and i+1:N of W */
/* >                contain eigenvalues which have converged. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexGEeigen */

/*  ===================================================================== */
/* Subroutine */ int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *w, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *scale, doublereal *abnrm, 
	doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *info, ftnlen balanc_len, ftnlen 
	jobvl_len, ftnlen jobvr_len, ftnlen sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, k;
    static char job[1];
    static doublereal scl, dum[1], eps;
    static doublecomplex tmp;
    static char side[1];
    static doublereal anrm;
    static integer ierr, itau, iwrk, nout;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer icond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int cgebak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen, ftnlen), cgebal_(char *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen), slabad_(doublereal *, doublereal *);
    static logical scalea;
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal cscale;
    extern /* Subroutine */ int cgehrd_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), clacpy_(char *, integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int chseqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), ctrevc_(char *, char *, logical *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     doublereal *, integer *, ftnlen, ftnlen), cunghr_(integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, integer *), ctrsna_(char *, char *, 
	    logical *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublereal *, doublereal *
	    , integer *, integer *, doublecomplex *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static logical wantvl, wntsnb;
    static integer hswork;
    static logical wntsne;
    static doublereal smlnum;
    static logical lquery, wantvr, wntsnn, wntsnv;


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 342 "cgeevx.f"
    /* Parameter adjustments */
#line 342 "cgeevx.f"
    a_dim1 = *lda;
#line 342 "cgeevx.f"
    a_offset = 1 + a_dim1;
#line 342 "cgeevx.f"
    a -= a_offset;
#line 342 "cgeevx.f"
    --w;
#line 342 "cgeevx.f"
    vl_dim1 = *ldvl;
#line 342 "cgeevx.f"
    vl_offset = 1 + vl_dim1;
#line 342 "cgeevx.f"
    vl -= vl_offset;
#line 342 "cgeevx.f"
    vr_dim1 = *ldvr;
#line 342 "cgeevx.f"
    vr_offset = 1 + vr_dim1;
#line 342 "cgeevx.f"
    vr -= vr_offset;
#line 342 "cgeevx.f"
    --scale;
#line 342 "cgeevx.f"
    --rconde;
#line 342 "cgeevx.f"
    --rcondv;
#line 342 "cgeevx.f"
    --work;
#line 342 "cgeevx.f"
    --rwork;
#line 342 "cgeevx.f"

#line 342 "cgeevx.f"
    /* Function Body */
#line 342 "cgeevx.f"
    *info = 0;
#line 343 "cgeevx.f"
    lquery = *lwork == -1;
#line 344 "cgeevx.f"
    wantvl = lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1);
#line 345 "cgeevx.f"
    wantvr = lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1);
#line 346 "cgeevx.f"
    wntsnn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 347 "cgeevx.f"
    wntsne = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 348 "cgeevx.f"
    wntsnv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 349 "cgeevx.f"
    wntsnb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 350 "cgeevx.f"
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "S", (
	    ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) 
	    || lsame_(balanc, "B", (ftnlen)1, (ftnlen)1))) {
#line 352 "cgeevx.f"
	*info = -1;
#line 353 "cgeevx.f"
    } else if (! wantvl && ! lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 354 "cgeevx.f"
	*info = -2;
#line 355 "cgeevx.f"
    } else if (! wantvr && ! lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 356 "cgeevx.f"
	*info = -3;
#line 357 "cgeevx.f"
    } else if (! (wntsnn || wntsne || wntsnb || wntsnv) || (wntsne || wntsnb) 
	    && ! (wantvl && wantvr)) {
#line 360 "cgeevx.f"
	*info = -4;
#line 361 "cgeevx.f"
    } else if (*n < 0) {
#line 362 "cgeevx.f"
	*info = -5;
#line 363 "cgeevx.f"
    } else if (*lda < max(1,*n)) {
#line 364 "cgeevx.f"
	*info = -7;
#line 365 "cgeevx.f"
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
#line 366 "cgeevx.f"
	*info = -10;
#line 367 "cgeevx.f"
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
#line 368 "cgeevx.f"
	*info = -12;
#line 369 "cgeevx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to real */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by CHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 382 "cgeevx.f"
    if (*info == 0) {
#line 383 "cgeevx.f"
	if (*n == 0) {
#line 384 "cgeevx.f"
	    minwrk = 1;
#line 385 "cgeevx.f"
	    maxwrk = 1;
#line 386 "cgeevx.f"
	} else {
#line 387 "cgeevx.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "CGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);

#line 389 "cgeevx.f"
	    if (wantvl) {
#line 390 "cgeevx.f"
		chseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vl[
			vl_offset], ldvl, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 392 "cgeevx.f"
	    } else if (wantvr) {
#line 393 "cgeevx.f"
		chseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 395 "cgeevx.f"
	    } else {
#line 396 "cgeevx.f"
		if (wntsnn) {
#line 397 "cgeevx.f"
		    chseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &
			    vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			    ftnlen)1, (ftnlen)1);
#line 399 "cgeevx.f"
		} else {
#line 400 "cgeevx.f"
		    chseqr_("S", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &
			    vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			    ftnlen)1, (ftnlen)1);
#line 402 "cgeevx.f"
		}
#line 403 "cgeevx.f"
	    }
#line 404 "cgeevx.f"
	    hswork = (integer) work[1].r;

#line 406 "cgeevx.f"
	    if (! wantvl && ! wantvr) {
#line 407 "cgeevx.f"
		minwrk = *n << 1;
#line 408 "cgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 408 "cgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + (*n << 1);
#line 408 "cgeevx.f"
		    minwrk = max(i__1,i__2);
#line 408 "cgeevx.f"
		}
#line 410 "cgeevx.f"
		maxwrk = max(maxwrk,hswork);
#line 411 "cgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 411 "cgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + (*n << 1);
#line 411 "cgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 411 "cgeevx.f"
		}
#line 413 "cgeevx.f"
	    } else {
#line 414 "cgeevx.f"
		minwrk = *n << 1;
#line 415 "cgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 415 "cgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + (*n << 1);
#line 415 "cgeevx.f"
		    minwrk = max(i__1,i__2);
#line 415 "cgeevx.f"
		}
#line 417 "cgeevx.f"
		maxwrk = max(maxwrk,hswork);
/* Computing MAX */
#line 418 "cgeevx.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 418 "cgeevx.f"
		maxwrk = max(i__1,i__2);
#line 420 "cgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 420 "cgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + (*n << 1);
#line 420 "cgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 420 "cgeevx.f"
		}
/* Computing MAX */
#line 422 "cgeevx.f"
		i__1 = maxwrk, i__2 = *n << 1;
#line 422 "cgeevx.f"
		maxwrk = max(i__1,i__2);
#line 423 "cgeevx.f"
	    }
#line 424 "cgeevx.f"
	    maxwrk = max(maxwrk,minwrk);
#line 425 "cgeevx.f"
	}
#line 426 "cgeevx.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 428 "cgeevx.f"
	if (*lwork < minwrk && ! lquery) {
#line 429 "cgeevx.f"
	    *info = -20;
#line 430 "cgeevx.f"
	}
#line 431 "cgeevx.f"
    }

#line 433 "cgeevx.f"
    if (*info != 0) {
#line 434 "cgeevx.f"
	i__1 = -(*info);
#line 434 "cgeevx.f"
	xerbla_("CGEEVX", &i__1, (ftnlen)6);
#line 435 "cgeevx.f"
	return 0;
#line 436 "cgeevx.f"
    } else if (lquery) {
#line 437 "cgeevx.f"
	return 0;
#line 438 "cgeevx.f"
    }

/*     Quick return if possible */

#line 442 "cgeevx.f"
    if (*n == 0) {
#line 442 "cgeevx.f"
	return 0;
#line 442 "cgeevx.f"
    }

/*     Get machine constants */

#line 447 "cgeevx.f"
    eps = slamch_("P", (ftnlen)1);
#line 448 "cgeevx.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 449 "cgeevx.f"
    bignum = 1. / smlnum;
#line 450 "cgeevx.f"
    slabad_(&smlnum, &bignum);
#line 451 "cgeevx.f"
    smlnum = sqrt(smlnum) / eps;
#line 452 "cgeevx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 456 "cgeevx.f"
    icond = 0;
#line 457 "cgeevx.f"
    anrm = clange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 458 "cgeevx.f"
    scalea = FALSE_;
#line 459 "cgeevx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 460 "cgeevx.f"
	scalea = TRUE_;
#line 461 "cgeevx.f"
	cscale = smlnum;
#line 462 "cgeevx.f"
    } else if (anrm > bignum) {
#line 463 "cgeevx.f"
	scalea = TRUE_;
#line 464 "cgeevx.f"
	cscale = bignum;
#line 465 "cgeevx.f"
    }
#line 466 "cgeevx.f"
    if (scalea) {
#line 466 "cgeevx.f"
	clascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 466 "cgeevx.f"
    }

/*     Balance the matrix and compute ABNRM */

#line 471 "cgeevx.f"
    cgebal_(balanc, n, &a[a_offset], lda, ilo, ihi, &scale[1], &ierr, (ftnlen)
	    1);
#line 472 "cgeevx.f"
    *abnrm = clange_("1", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 473 "cgeevx.f"
    if (scalea) {
#line 474 "cgeevx.f"
	dum[0] = *abnrm;
#line 475 "cgeevx.f"
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &c__1, &
		ierr, (ftnlen)1);
#line 476 "cgeevx.f"
	*abnrm = dum[0];
#line 477 "cgeevx.f"
    }

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

#line 483 "cgeevx.f"
    itau = 1;
#line 484 "cgeevx.f"
    iwrk = itau + *n;
#line 485 "cgeevx.f"
    i__1 = *lwork - iwrk + 1;
#line 485 "cgeevx.f"
    cgehrd_(n, ilo, ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &
	    ierr);

#line 488 "cgeevx.f"
    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

#line 493 "cgeevx.f"
	*(unsigned char *)side = 'L';
#line 494 "cgeevx.f"
	clacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1)
		;

/*        Generate unitary matrix in VL */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 500 "cgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 500 "cgeevx.f"
	cunghr_(n, ilo, ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 507 "cgeevx.f"
	iwrk = itau;
#line 508 "cgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 508 "cgeevx.f"
	chseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &w[1], &vl[
		vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 511 "cgeevx.f"
	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

#line 516 "cgeevx.f"
	    *(unsigned char *)side = 'B';
#line 517 "cgeevx.f"
	    clacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (
		    ftnlen)1);
#line 518 "cgeevx.f"
	}

#line 520 "cgeevx.f"
    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

#line 525 "cgeevx.f"
	*(unsigned char *)side = 'R';
#line 526 "cgeevx.f"
	clacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1)
		;

/*        Generate unitary matrix in VR */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 532 "cgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 532 "cgeevx.f"
	cunghr_(n, ilo, ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 539 "cgeevx.f"
	iwrk = itau;
#line 540 "cgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 540 "cgeevx.f"
	chseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 543 "cgeevx.f"
    } else {

/*        Compute eigenvalues only */
/*        If condition numbers desired, compute Schur form */

#line 548 "cgeevx.f"
	if (wntsnn) {
#line 549 "cgeevx.f"
	    *(unsigned char *)job = 'E';
#line 550 "cgeevx.f"
	} else {
#line 551 "cgeevx.f"
	    *(unsigned char *)job = 'S';
#line 552 "cgeevx.f"
	}

/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 557 "cgeevx.f"
	iwrk = itau;
#line 558 "cgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 558 "cgeevx.f"
	chseqr_(job, "N", n, ilo, ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);
#line 560 "cgeevx.f"
    }

/*     If INFO > 0 from CHSEQR, then quit */

#line 564 "cgeevx.f"
    if (*info > 0) {
#line 564 "cgeevx.f"
	goto L50;
#line 564 "cgeevx.f"
    }

#line 567 "cgeevx.f"
    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (CWorkspace: need 2*N) */
/*        (RWorkspace: need N) */

#line 573 "cgeevx.f"
	ctrevc_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
		 &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &rwork[1], &
		ierr, (ftnlen)1, (ftnlen)1);
#line 575 "cgeevx.f"
    }

/*     Compute condition numbers if desired */
/*     (CWorkspace: need N*N+2*N unless SENSE = 'E') */
/*     (RWorkspace: need 2*N unless SENSE = 'E') */

#line 581 "cgeevx.f"
    if (! wntsnn) {
#line 582 "cgeevx.f"
	ctrsna_(sense, "A", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, &rconde[1], &rcondv[1], n, &nout, 
		&work[iwrk], n, &rwork[1], &icond, (ftnlen)1, (ftnlen)1);
#line 585 "cgeevx.f"
    }

#line 587 "cgeevx.f"
    if (wantvl) {

/*        Undo balancing of left eigenvectors */

#line 591 "cgeevx.f"
	cgebak_(balanc, "L", n, ilo, ihi, &scale[1], n, &vl[vl_offset], ldvl, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize left eigenvectors and make largest component real */

#line 596 "cgeevx.f"
	i__1 = *n;
#line 596 "cgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 597 "cgeevx.f"
	    scl = 1. / scnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 598 "cgeevx.f"
	    csscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 599 "cgeevx.f"
	    i__2 = *n;
#line 599 "cgeevx.f"
	    for (k = 1; k <= i__2; ++k) {
#line 600 "cgeevx.f"
		i__3 = k + i__ * vl_dim1;
/* Computing 2nd power */
#line 600 "cgeevx.f"
		d__1 = vl[i__3].r;
/* Computing 2nd power */
#line 600 "cgeevx.f"
		d__2 = d_imag(&vl[k + i__ * vl_dim1]);
#line 600 "cgeevx.f"
		rwork[k] = d__1 * d__1 + d__2 * d__2;
#line 602 "cgeevx.f"
/* L10: */
#line 602 "cgeevx.f"
	    }
#line 603 "cgeevx.f"
	    k = isamax_(n, &rwork[1], &c__1);
#line 604 "cgeevx.f"
	    d_cnjg(&z__2, &vl[k + i__ * vl_dim1]);
#line 604 "cgeevx.f"
	    d__1 = sqrt(rwork[k]);
#line 604 "cgeevx.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 604 "cgeevx.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 605 "cgeevx.f"
	    cscal_(n, &tmp, &vl[i__ * vl_dim1 + 1], &c__1);
#line 606 "cgeevx.f"
	    i__2 = k + i__ * vl_dim1;
#line 606 "cgeevx.f"
	    i__3 = k + i__ * vl_dim1;
#line 606 "cgeevx.f"
	    d__1 = vl[i__3].r;
#line 606 "cgeevx.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 606 "cgeevx.f"
	    vl[i__2].r = z__1.r, vl[i__2].i = z__1.i;
#line 607 "cgeevx.f"
/* L20: */
#line 607 "cgeevx.f"
	}
#line 608 "cgeevx.f"
    }

#line 610 "cgeevx.f"
    if (wantvr) {

/*        Undo balancing of right eigenvectors */

#line 614 "cgeevx.f"
	cgebak_(balanc, "R", n, ilo, ihi, &scale[1], n, &vr[vr_offset], ldvr, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize right eigenvectors and make largest component real */

#line 619 "cgeevx.f"
	i__1 = *n;
#line 619 "cgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 620 "cgeevx.f"
	    scl = 1. / scnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 621 "cgeevx.f"
	    csscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 622 "cgeevx.f"
	    i__2 = *n;
#line 622 "cgeevx.f"
	    for (k = 1; k <= i__2; ++k) {
#line 623 "cgeevx.f"
		i__3 = k + i__ * vr_dim1;
/* Computing 2nd power */
#line 623 "cgeevx.f"
		d__1 = vr[i__3].r;
/* Computing 2nd power */
#line 623 "cgeevx.f"
		d__2 = d_imag(&vr[k + i__ * vr_dim1]);
#line 623 "cgeevx.f"
		rwork[k] = d__1 * d__1 + d__2 * d__2;
#line 625 "cgeevx.f"
/* L30: */
#line 625 "cgeevx.f"
	    }
#line 626 "cgeevx.f"
	    k = isamax_(n, &rwork[1], &c__1);
#line 627 "cgeevx.f"
	    d_cnjg(&z__2, &vr[k + i__ * vr_dim1]);
#line 627 "cgeevx.f"
	    d__1 = sqrt(rwork[k]);
#line 627 "cgeevx.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 627 "cgeevx.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 628 "cgeevx.f"
	    cscal_(n, &tmp, &vr[i__ * vr_dim1 + 1], &c__1);
#line 629 "cgeevx.f"
	    i__2 = k + i__ * vr_dim1;
#line 629 "cgeevx.f"
	    i__3 = k + i__ * vr_dim1;
#line 629 "cgeevx.f"
	    d__1 = vr[i__3].r;
#line 629 "cgeevx.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 629 "cgeevx.f"
	    vr[i__2].r = z__1.r, vr[i__2].i = z__1.i;
#line 630 "cgeevx.f"
/* L40: */
#line 630 "cgeevx.f"
	}
#line 631 "cgeevx.f"
    }

/*     Undo scaling if necessary */

#line 635 "cgeevx.f"
L50:
#line 636 "cgeevx.f"
    if (scalea) {
#line 637 "cgeevx.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 637 "cgeevx.f"
	i__3 = *n - *info;
#line 637 "cgeevx.f"
	i__2 = max(i__3,1);
#line 637 "cgeevx.f"
	clascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[*info + 1]
		, &i__2, &ierr, (ftnlen)1);
#line 639 "cgeevx.f"
	if (*info == 0) {
#line 640 "cgeevx.f"
	    if ((wntsnv || wntsnb) && icond == 0) {
#line 640 "cgeevx.f"
		slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &rcondv[
			1], n, &ierr, (ftnlen)1);
#line 640 "cgeevx.f"
	    }
#line 643 "cgeevx.f"
	} else {
#line 644 "cgeevx.f"
	    i__1 = *ilo - 1;
#line 644 "cgeevx.f"
	    clascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[1], n,
		     &ierr, (ftnlen)1);
#line 645 "cgeevx.f"
	}
#line 646 "cgeevx.f"
    }

#line 648 "cgeevx.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 649 "cgeevx.f"
    return 0;

/*     End of CGEEVX */

} /* cgeevx_ */

