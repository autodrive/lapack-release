#line 1 "zgeevx.f"
/* zgeevx.f -- translated by f2c (version 20100827).
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

#line 1 "zgeevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> ZGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, */
/*                          LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, */
/*                          RCONDV, WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE */
/*       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N */
/*       DOUBLE PRECISION   ABNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RCONDE( * ), RCONDV( * ), RWORK( * ), */
/*      $                   SCALE( * ) */
/*       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEEVX computes for an N-by-N complex nonsymmetric matrix A, the */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          W is COMPLEX*16 array, dimension (N) */
/* >          W contains the computed eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is COMPLEX*16 array, dimension (LDVL,N) */
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
/* >          VR is COMPLEX*16 array, dimension (LDVR,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \date June 2016 */

/*  @precisions fortran z -> c */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
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
    static integer lwork_trevc__;
    static char side[1];
    static doublereal anrm;
    static integer ierr, itau, iwrk, nout, icond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    static logical scalea;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), zgebak_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen), zgebal_(char *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int zgehrd_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), zlacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen);
    static integer minwrk, maxwrk;
    static logical wantvl, wntsnb;
    static integer hswork;
    static logical wntsne;
    static doublereal smlnum;
    extern /* Subroutine */ int zhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static logical lquery, wantvr;
    extern /* Subroutine */ int ztrsna_(char *, char *, logical *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, integer *,
	     integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), zunghr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static logical wntsnn, wntsnv;
    extern /* Subroutine */ int ztrevc3_(char *, char *, logical *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 345 "zgeevx.f"
    /* Parameter adjustments */
#line 345 "zgeevx.f"
    a_dim1 = *lda;
#line 345 "zgeevx.f"
    a_offset = 1 + a_dim1;
#line 345 "zgeevx.f"
    a -= a_offset;
#line 345 "zgeevx.f"
    --w;
#line 345 "zgeevx.f"
    vl_dim1 = *ldvl;
#line 345 "zgeevx.f"
    vl_offset = 1 + vl_dim1;
#line 345 "zgeevx.f"
    vl -= vl_offset;
#line 345 "zgeevx.f"
    vr_dim1 = *ldvr;
#line 345 "zgeevx.f"
    vr_offset = 1 + vr_dim1;
#line 345 "zgeevx.f"
    vr -= vr_offset;
#line 345 "zgeevx.f"
    --scale;
#line 345 "zgeevx.f"
    --rconde;
#line 345 "zgeevx.f"
    --rcondv;
#line 345 "zgeevx.f"
    --work;
#line 345 "zgeevx.f"
    --rwork;
#line 345 "zgeevx.f"

#line 345 "zgeevx.f"
    /* Function Body */
#line 345 "zgeevx.f"
    *info = 0;
#line 346 "zgeevx.f"
    lquery = *lwork == -1;
#line 347 "zgeevx.f"
    wantvl = lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1);
#line 348 "zgeevx.f"
    wantvr = lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1);
#line 349 "zgeevx.f"
    wntsnn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 350 "zgeevx.f"
    wntsne = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 351 "zgeevx.f"
    wntsnv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 352 "zgeevx.f"
    wntsnb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 353 "zgeevx.f"
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "S", (
	    ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) 
	    || lsame_(balanc, "B", (ftnlen)1, (ftnlen)1))) {
#line 355 "zgeevx.f"
	*info = -1;
#line 356 "zgeevx.f"
    } else if (! wantvl && ! lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 357 "zgeevx.f"
	*info = -2;
#line 358 "zgeevx.f"
    } else if (! wantvr && ! lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 359 "zgeevx.f"
	*info = -3;
#line 360 "zgeevx.f"
    } else if (! (wntsnn || wntsne || wntsnb || wntsnv) || (wntsne || wntsnb) 
	    && ! (wantvl && wantvr)) {
#line 363 "zgeevx.f"
	*info = -4;
#line 364 "zgeevx.f"
    } else if (*n < 0) {
#line 365 "zgeevx.f"
	*info = -5;
#line 366 "zgeevx.f"
    } else if (*lda < max(1,*n)) {
#line 367 "zgeevx.f"
	*info = -7;
#line 368 "zgeevx.f"
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
#line 369 "zgeevx.f"
	*info = -10;
#line 370 "zgeevx.f"
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
#line 371 "zgeevx.f"
	*info = -12;
#line 372 "zgeevx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to real */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by ZHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 385 "zgeevx.f"
    if (*info == 0) {
#line 386 "zgeevx.f"
	if (*n == 0) {
#line 387 "zgeevx.f"
	    minwrk = 1;
#line 388 "zgeevx.f"
	    maxwrk = 1;
#line 389 "zgeevx.f"
	} else {
#line 390 "zgeevx.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);

#line 392 "zgeevx.f"
	    if (wantvl) {
#line 393 "zgeevx.f"
		ztrevc3_("L", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &rwork[1], &c_n1, &ierr, (ftnlen)1, (
			ftnlen)1);
#line 396 "zgeevx.f"
		lwork_trevc__ = (integer) work[1].r;
#line 397 "zgeevx.f"
		maxwrk = max(maxwrk,lwork_trevc__);
#line 398 "zgeevx.f"
		zhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vl[
			vl_offset], ldvl, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 400 "zgeevx.f"
	    } else if (wantvr) {
#line 401 "zgeevx.f"
		ztrevc3_("R", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &rwork[1], &c_n1, &ierr, (ftnlen)1, (
			ftnlen)1);
#line 404 "zgeevx.f"
		lwork_trevc__ = (integer) work[1].r;
#line 405 "zgeevx.f"
		maxwrk = max(maxwrk,lwork_trevc__);
#line 406 "zgeevx.f"
		zhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 408 "zgeevx.f"
	    } else {
#line 409 "zgeevx.f"
		if (wntsnn) {
#line 410 "zgeevx.f"
		    zhseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &
			    vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			    ftnlen)1, (ftnlen)1);
#line 412 "zgeevx.f"
		} else {
#line 413 "zgeevx.f"
		    zhseqr_("S", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &
			    vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			    ftnlen)1, (ftnlen)1);
#line 415 "zgeevx.f"
		}
#line 416 "zgeevx.f"
	    }
#line 417 "zgeevx.f"
	    hswork = (integer) work[1].r;

#line 419 "zgeevx.f"
	    if (! wantvl && ! wantvr) {
#line 420 "zgeevx.f"
		minwrk = *n << 1;
#line 421 "zgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 421 "zgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + (*n << 1);
#line 421 "zgeevx.f"
		    minwrk = max(i__1,i__2);
#line 421 "zgeevx.f"
		}
#line 423 "zgeevx.f"
		maxwrk = max(maxwrk,hswork);
#line 424 "zgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 424 "zgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + (*n << 1);
#line 424 "zgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 424 "zgeevx.f"
		}
#line 426 "zgeevx.f"
	    } else {
#line 427 "zgeevx.f"
		minwrk = *n << 1;
#line 428 "zgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 428 "zgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + (*n << 1);
#line 428 "zgeevx.f"
		    minwrk = max(i__1,i__2);
#line 428 "zgeevx.f"
		}
#line 430 "zgeevx.f"
		maxwrk = max(maxwrk,hswork);
/* Computing MAX */
#line 431 "zgeevx.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "ZUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 431 "zgeevx.f"
		maxwrk = max(i__1,i__2);
#line 433 "zgeevx.f"
		if (! (wntsnn || wntsne)) {
/* Computing MAX */
#line 433 "zgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + (*n << 1);
#line 433 "zgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 433 "zgeevx.f"
		}
/* Computing MAX */
#line 435 "zgeevx.f"
		i__1 = maxwrk, i__2 = *n << 1;
#line 435 "zgeevx.f"
		maxwrk = max(i__1,i__2);
#line 436 "zgeevx.f"
	    }
#line 437 "zgeevx.f"
	    maxwrk = max(maxwrk,minwrk);
#line 438 "zgeevx.f"
	}
#line 439 "zgeevx.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 441 "zgeevx.f"
	if (*lwork < minwrk && ! lquery) {
#line 442 "zgeevx.f"
	    *info = -20;
#line 443 "zgeevx.f"
	}
#line 444 "zgeevx.f"
    }

#line 446 "zgeevx.f"
    if (*info != 0) {
#line 447 "zgeevx.f"
	i__1 = -(*info);
#line 447 "zgeevx.f"
	xerbla_("ZGEEVX", &i__1, (ftnlen)6);
#line 448 "zgeevx.f"
	return 0;
#line 449 "zgeevx.f"
    } else if (lquery) {
#line 450 "zgeevx.f"
	return 0;
#line 451 "zgeevx.f"
    }

/*     Quick return if possible */

#line 455 "zgeevx.f"
    if (*n == 0) {
#line 455 "zgeevx.f"
	return 0;
#line 455 "zgeevx.f"
    }

/*     Get machine constants */

#line 460 "zgeevx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 461 "zgeevx.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 462 "zgeevx.f"
    bignum = 1. / smlnum;
#line 463 "zgeevx.f"
    dlabad_(&smlnum, &bignum);
#line 464 "zgeevx.f"
    smlnum = sqrt(smlnum) / eps;
#line 465 "zgeevx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 469 "zgeevx.f"
    icond = 0;
#line 470 "zgeevx.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 471 "zgeevx.f"
    scalea = FALSE_;
#line 472 "zgeevx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 473 "zgeevx.f"
	scalea = TRUE_;
#line 474 "zgeevx.f"
	cscale = smlnum;
#line 475 "zgeevx.f"
    } else if (anrm > bignum) {
#line 476 "zgeevx.f"
	scalea = TRUE_;
#line 477 "zgeevx.f"
	cscale = bignum;
#line 478 "zgeevx.f"
    }
#line 479 "zgeevx.f"
    if (scalea) {
#line 479 "zgeevx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 479 "zgeevx.f"
    }

/*     Balance the matrix and compute ABNRM */

#line 484 "zgeevx.f"
    zgebal_(balanc, n, &a[a_offset], lda, ilo, ihi, &scale[1], &ierr, (ftnlen)
	    1);
#line 485 "zgeevx.f"
    *abnrm = zlange_("1", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 486 "zgeevx.f"
    if (scalea) {
#line 487 "zgeevx.f"
	dum[0] = *abnrm;
#line 488 "zgeevx.f"
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &c__1, &
		ierr, (ftnlen)1);
#line 489 "zgeevx.f"
	*abnrm = dum[0];
#line 490 "zgeevx.f"
    }

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

#line 496 "zgeevx.f"
    itau = 1;
#line 497 "zgeevx.f"
    iwrk = itau + *n;
#line 498 "zgeevx.f"
    i__1 = *lwork - iwrk + 1;
#line 498 "zgeevx.f"
    zgehrd_(n, ilo, ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &
	    ierr);

#line 501 "zgeevx.f"
    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

#line 506 "zgeevx.f"
	*(unsigned char *)side = 'L';
#line 507 "zgeevx.f"
	zlacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1)
		;

/*        Generate unitary matrix in VL */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 513 "zgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 513 "zgeevx.f"
	zunghr_(n, ilo, ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 520 "zgeevx.f"
	iwrk = itau;
#line 521 "zgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 521 "zgeevx.f"
	zhseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &w[1], &vl[
		vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 524 "zgeevx.f"
	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

#line 529 "zgeevx.f"
	    *(unsigned char *)side = 'B';
#line 530 "zgeevx.f"
	    zlacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (
		    ftnlen)1);
#line 531 "zgeevx.f"
	}

#line 533 "zgeevx.f"
    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

#line 538 "zgeevx.f"
	*(unsigned char *)side = 'R';
#line 539 "zgeevx.f"
	zlacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1)
		;

/*        Generate unitary matrix in VR */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 545 "zgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 545 "zgeevx.f"
	zunghr_(n, ilo, ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 552 "zgeevx.f"
	iwrk = itau;
#line 553 "zgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 553 "zgeevx.f"
	zhseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 556 "zgeevx.f"
    } else {

/*        Compute eigenvalues only */
/*        If condition numbers desired, compute Schur form */

#line 561 "zgeevx.f"
	if (wntsnn) {
#line 562 "zgeevx.f"
	    *(unsigned char *)job = 'E';
#line 563 "zgeevx.f"
	} else {
#line 564 "zgeevx.f"
	    *(unsigned char *)job = 'S';
#line 565 "zgeevx.f"
	}

/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 570 "zgeevx.f"
	iwrk = itau;
#line 571 "zgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 571 "zgeevx.f"
	zhseqr_(job, "N", n, ilo, ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);
#line 573 "zgeevx.f"
    }

/*     If INFO .NE. 0 from ZHSEQR, then quit */

#line 577 "zgeevx.f"
    if (*info != 0) {
#line 577 "zgeevx.f"
	goto L50;
#line 577 "zgeevx.f"
    }

#line 580 "zgeevx.f"
    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (CWorkspace: need 2*N, prefer N + 2*N*NB) */
/*        (RWorkspace: need N) */

#line 586 "zgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 586 "zgeevx.f"
	ztrevc3_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &i__1, &
		rwork[1], n, &ierr, (ftnlen)1, (ftnlen)1);
#line 589 "zgeevx.f"
    }

/*     Compute condition numbers if desired */
/*     (CWorkspace: need N*N+2*N unless SENSE = 'E') */
/*     (RWorkspace: need 2*N unless SENSE = 'E') */

#line 595 "zgeevx.f"
    if (! wntsnn) {
#line 596 "zgeevx.f"
	ztrsna_(sense, "A", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, &rconde[1], &rcondv[1], n, &nout, 
		&work[iwrk], n, &rwork[1], &icond, (ftnlen)1, (ftnlen)1);
#line 599 "zgeevx.f"
    }

#line 601 "zgeevx.f"
    if (wantvl) {

/*        Undo balancing of left eigenvectors */

#line 605 "zgeevx.f"
	zgebak_(balanc, "L", n, ilo, ihi, &scale[1], n, &vl[vl_offset], ldvl, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize left eigenvectors and make largest component real */

#line 610 "zgeevx.f"
	i__1 = *n;
#line 610 "zgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 611 "zgeevx.f"
	    scl = 1. / dznrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 612 "zgeevx.f"
	    zdscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 613 "zgeevx.f"
	    i__2 = *n;
#line 613 "zgeevx.f"
	    for (k = 1; k <= i__2; ++k) {
#line 614 "zgeevx.f"
		i__3 = k + i__ * vl_dim1;
/* Computing 2nd power */
#line 614 "zgeevx.f"
		d__1 = vl[i__3].r;
/* Computing 2nd power */
#line 614 "zgeevx.f"
		d__2 = d_imag(&vl[k + i__ * vl_dim1]);
#line 614 "zgeevx.f"
		rwork[k] = d__1 * d__1 + d__2 * d__2;
#line 616 "zgeevx.f"
/* L10: */
#line 616 "zgeevx.f"
	    }
#line 617 "zgeevx.f"
	    k = idamax_(n, &rwork[1], &c__1);
#line 618 "zgeevx.f"
	    d_cnjg(&z__2, &vl[k + i__ * vl_dim1]);
#line 618 "zgeevx.f"
	    d__1 = sqrt(rwork[k]);
#line 618 "zgeevx.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 618 "zgeevx.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 619 "zgeevx.f"
	    zscal_(n, &tmp, &vl[i__ * vl_dim1 + 1], &c__1);
#line 620 "zgeevx.f"
	    i__2 = k + i__ * vl_dim1;
#line 620 "zgeevx.f"
	    i__3 = k + i__ * vl_dim1;
#line 620 "zgeevx.f"
	    d__1 = vl[i__3].r;
#line 620 "zgeevx.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 620 "zgeevx.f"
	    vl[i__2].r = z__1.r, vl[i__2].i = z__1.i;
#line 621 "zgeevx.f"
/* L20: */
#line 621 "zgeevx.f"
	}
#line 622 "zgeevx.f"
    }

#line 624 "zgeevx.f"
    if (wantvr) {

/*        Undo balancing of right eigenvectors */

#line 628 "zgeevx.f"
	zgebak_(balanc, "R", n, ilo, ihi, &scale[1], n, &vr[vr_offset], ldvr, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize right eigenvectors and make largest component real */

#line 633 "zgeevx.f"
	i__1 = *n;
#line 633 "zgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 634 "zgeevx.f"
	    scl = 1. / dznrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 635 "zgeevx.f"
	    zdscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 636 "zgeevx.f"
	    i__2 = *n;
#line 636 "zgeevx.f"
	    for (k = 1; k <= i__2; ++k) {
#line 637 "zgeevx.f"
		i__3 = k + i__ * vr_dim1;
/* Computing 2nd power */
#line 637 "zgeevx.f"
		d__1 = vr[i__3].r;
/* Computing 2nd power */
#line 637 "zgeevx.f"
		d__2 = d_imag(&vr[k + i__ * vr_dim1]);
#line 637 "zgeevx.f"
		rwork[k] = d__1 * d__1 + d__2 * d__2;
#line 639 "zgeevx.f"
/* L30: */
#line 639 "zgeevx.f"
	    }
#line 640 "zgeevx.f"
	    k = idamax_(n, &rwork[1], &c__1);
#line 641 "zgeevx.f"
	    d_cnjg(&z__2, &vr[k + i__ * vr_dim1]);
#line 641 "zgeevx.f"
	    d__1 = sqrt(rwork[k]);
#line 641 "zgeevx.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 641 "zgeevx.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 642 "zgeevx.f"
	    zscal_(n, &tmp, &vr[i__ * vr_dim1 + 1], &c__1);
#line 643 "zgeevx.f"
	    i__2 = k + i__ * vr_dim1;
#line 643 "zgeevx.f"
	    i__3 = k + i__ * vr_dim1;
#line 643 "zgeevx.f"
	    d__1 = vr[i__3].r;
#line 643 "zgeevx.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 643 "zgeevx.f"
	    vr[i__2].r = z__1.r, vr[i__2].i = z__1.i;
#line 644 "zgeevx.f"
/* L40: */
#line 644 "zgeevx.f"
	}
#line 645 "zgeevx.f"
    }

/*     Undo scaling if necessary */

#line 649 "zgeevx.f"
L50:
#line 650 "zgeevx.f"
    if (scalea) {
#line 651 "zgeevx.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 651 "zgeevx.f"
	i__3 = *n - *info;
#line 651 "zgeevx.f"
	i__2 = max(i__3,1);
#line 651 "zgeevx.f"
	zlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[*info + 1]
		, &i__2, &ierr, (ftnlen)1);
#line 653 "zgeevx.f"
	if (*info == 0) {
#line 654 "zgeevx.f"
	    if ((wntsnv || wntsnb) && icond == 0) {
#line 654 "zgeevx.f"
		dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &rcondv[
			1], n, &ierr, (ftnlen)1);
#line 654 "zgeevx.f"
	    }
#line 657 "zgeevx.f"
	} else {
#line 658 "zgeevx.f"
	    i__1 = *ilo - 1;
#line 658 "zgeevx.f"
	    zlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[1], n,
		     &ierr, (ftnlen)1);
#line 659 "zgeevx.f"
	}
#line 660 "zgeevx.f"
    }

#line 662 "zgeevx.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 663 "zgeevx.f"
    return 0;

/*     End of ZGEEVX */

} /* zgeevx_ */

