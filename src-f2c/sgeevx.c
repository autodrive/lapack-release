#line 1 "sgeevx.f"
/* sgeevx.f -- translated by f2c (version 20100827).
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

#line 1 "sgeevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> SGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, */
/*                          VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, */
/*                          RCONDE, RCONDV, WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE */
/*       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N */
/*       REAL               ABNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), RCONDE( * ), RCONDV( * ), */
/*      $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WI( * ), WORK( * ), WR( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEEVX computes for an N-by-N real nonsymmetric matrix A, the */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          WR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is REAL array, dimension (N) */
/* >          WR and WI contain the real and imaginary parts, */
/* >          respectively, of the computed eigenvalues.  Complex */
/* >          conjugate pairs of eigenvalues will appear consecutively */
/* >          with the eigenvalue having the positive imaginary part */
/* >          first. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is REAL array, dimension (LDVL,N) */
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
/* >          VR is REAL array, dimension (LDVR,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \date June 2016 */

/*  @generated from dgeevx.f, fortran d -> s, Tue Apr 19 01:47:44 2016 */

/* > \ingroup realGEeigen */

/*  ===================================================================== */
/* Subroutine */ int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
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
    static integer lwork_trevc__;
    static char side[1];
    static doublereal anrm;
    static integer ierr, itau, iwrk, nout;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static integer icond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal slapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    static logical scalea;
    static doublereal cscale;
    extern /* Subroutine */ int sgebak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sgebal_(char *, integer *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int sgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), sorghr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), shseqr_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    extern /* Subroutine */ int strsna_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen);
    static logical wantvl, wntsnb;
    static integer hswork;
    static logical wntsne;
    static doublereal smlnum;
    static logical lquery, wantvr, wntsnn, wntsnv;
    extern /* Subroutine */ int strevc3_(char *, char *, logical *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.1) -- */
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

#line 364 "sgeevx.f"
    /* Parameter adjustments */
#line 364 "sgeevx.f"
    a_dim1 = *lda;
#line 364 "sgeevx.f"
    a_offset = 1 + a_dim1;
#line 364 "sgeevx.f"
    a -= a_offset;
#line 364 "sgeevx.f"
    --wr;
#line 364 "sgeevx.f"
    --wi;
#line 364 "sgeevx.f"
    vl_dim1 = *ldvl;
#line 364 "sgeevx.f"
    vl_offset = 1 + vl_dim1;
#line 364 "sgeevx.f"
    vl -= vl_offset;
#line 364 "sgeevx.f"
    vr_dim1 = *ldvr;
#line 364 "sgeevx.f"
    vr_offset = 1 + vr_dim1;
#line 364 "sgeevx.f"
    vr -= vr_offset;
#line 364 "sgeevx.f"
    --scale;
#line 364 "sgeevx.f"
    --rconde;
#line 364 "sgeevx.f"
    --rcondv;
#line 364 "sgeevx.f"
    --work;
#line 364 "sgeevx.f"
    --iwork;
#line 364 "sgeevx.f"

#line 364 "sgeevx.f"
    /* Function Body */
#line 364 "sgeevx.f"
    *info = 0;
#line 365 "sgeevx.f"
    lquery = *lwork == -1;
#line 366 "sgeevx.f"
    wantvl = lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1);
#line 367 "sgeevx.f"
    wantvr = lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1);
#line 368 "sgeevx.f"
    wntsnn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 369 "sgeevx.f"
    wntsne = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 370 "sgeevx.f"
    wntsnv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 371 "sgeevx.f"
    wntsnb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 372 "sgeevx.f"
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "S", (
	    ftnlen)1, (ftnlen)1) || lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) 
	    || lsame_(balanc, "B", (ftnlen)1, (ftnlen)1))) {
#line 375 "sgeevx.f"
	*info = -1;
#line 376 "sgeevx.f"
    } else if (! wantvl && ! lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 377 "sgeevx.f"
	*info = -2;
#line 378 "sgeevx.f"
    } else if (! wantvr && ! lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 379 "sgeevx.f"
	*info = -3;
#line 380 "sgeevx.f"
    } else if (! (wntsnn || wntsne || wntsnb || wntsnv) || (wntsne || wntsnb) 
	    && ! (wantvl && wantvr)) {
#line 383 "sgeevx.f"
	*info = -4;
#line 384 "sgeevx.f"
    } else if (*n < 0) {
#line 385 "sgeevx.f"
	*info = -5;
#line 386 "sgeevx.f"
    } else if (*lda < max(1,*n)) {
#line 387 "sgeevx.f"
	*info = -7;
#line 388 "sgeevx.f"
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
#line 389 "sgeevx.f"
	*info = -11;
#line 390 "sgeevx.f"
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
#line 391 "sgeevx.f"
	*info = -13;
#line 392 "sgeevx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by SHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 404 "sgeevx.f"
    if (*info == 0) {
#line 405 "sgeevx.f"
	if (*n == 0) {
#line 406 "sgeevx.f"
	    minwrk = 1;
#line 407 "sgeevx.f"
	    maxwrk = 1;
#line 408 "sgeevx.f"
	} else {
#line 409 "sgeevx.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "SGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);

#line 411 "sgeevx.f"
	    if (wantvl) {
#line 412 "sgeevx.f"
		strevc3_("L", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
#line 415 "sgeevx.f"
		lwork_trevc__ = (integer) work[1];
/* Computing MAX */
#line 416 "sgeevx.f"
		i__1 = maxwrk, i__2 = *n + lwork_trevc__;
#line 416 "sgeevx.f"
		maxwrk = max(i__1,i__2);
#line 417 "sgeevx.f"
		shseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vl[vl_offset], ldvl, &work[1], &c_n1, info, (
			ftnlen)1, (ftnlen)1);
#line 419 "sgeevx.f"
	    } else if (wantvr) {
#line 420 "sgeevx.f"
		strevc3_("R", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
#line 423 "sgeevx.f"
		lwork_trevc__ = (integer) work[1];
/* Computing MAX */
#line 424 "sgeevx.f"
		i__1 = maxwrk, i__2 = *n + lwork_trevc__;
#line 424 "sgeevx.f"
		maxwrk = max(i__1,i__2);
#line 425 "sgeevx.f"
		shseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			ftnlen)1, (ftnlen)1);
#line 427 "sgeevx.f"
	    } else {
#line 428 "sgeevx.f"
		if (wntsnn) {
#line 429 "sgeevx.f"
		    shseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &wr[1], 
			    &wi[1], &vr[vr_offset], ldvr, &work[1], &c_n1, 
			    info, (ftnlen)1, (ftnlen)1);
#line 431 "sgeevx.f"
		} else {
#line 432 "sgeevx.f"
		    shseqr_("S", "N", n, &c__1, n, &a[a_offset], lda, &wr[1], 
			    &wi[1], &vr[vr_offset], ldvr, &work[1], &c_n1, 
			    info, (ftnlen)1, (ftnlen)1);
#line 434 "sgeevx.f"
		}
#line 435 "sgeevx.f"
	    }
#line 436 "sgeevx.f"
	    hswork = (integer) work[1];

#line 438 "sgeevx.f"
	    if (! wantvl && ! wantvr) {
#line 439 "sgeevx.f"
		minwrk = *n << 1;
#line 440 "sgeevx.f"
		if (! wntsnn) {
/* Computing MAX */
#line 440 "sgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + *n * 6;
#line 440 "sgeevx.f"
		    minwrk = max(i__1,i__2);
#line 440 "sgeevx.f"
		}
#line 442 "sgeevx.f"
		maxwrk = max(maxwrk,hswork);
#line 443 "sgeevx.f"
		if (! wntsnn) {
/* Computing MAX */
#line 443 "sgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + *n * 6;
#line 443 "sgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 443 "sgeevx.f"
		}
#line 445 "sgeevx.f"
	    } else {
#line 446 "sgeevx.f"
		minwrk = *n * 3;
#line 447 "sgeevx.f"
		if (! wntsnn && ! wntsne) {
/* Computing MAX */
#line 447 "sgeevx.f"
		    i__1 = minwrk, i__2 = *n * *n + *n * 6;
#line 447 "sgeevx.f"
		    minwrk = max(i__1,i__2);
#line 447 "sgeevx.f"
		}
#line 449 "sgeevx.f"
		maxwrk = max(maxwrk,hswork);
/* Computing MAX */
#line 450 "sgeevx.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "SORGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 450 "sgeevx.f"
		maxwrk = max(i__1,i__2);
#line 452 "sgeevx.f"
		if (! wntsnn && ! wntsne) {
/* Computing MAX */
#line 452 "sgeevx.f"
		    i__1 = maxwrk, i__2 = *n * *n + *n * 6;
#line 452 "sgeevx.f"
		    maxwrk = max(i__1,i__2);
#line 452 "sgeevx.f"
		}
/* Computing MAX */
#line 454 "sgeevx.f"
		i__1 = maxwrk, i__2 = *n * 3;
#line 454 "sgeevx.f"
		maxwrk = max(i__1,i__2);
#line 455 "sgeevx.f"
	    }
#line 456 "sgeevx.f"
	    maxwrk = max(maxwrk,minwrk);
#line 457 "sgeevx.f"
	}
#line 458 "sgeevx.f"
	work[1] = (doublereal) maxwrk;

#line 460 "sgeevx.f"
	if (*lwork < minwrk && ! lquery) {
#line 461 "sgeevx.f"
	    *info = -21;
#line 462 "sgeevx.f"
	}
#line 463 "sgeevx.f"
    }

#line 465 "sgeevx.f"
    if (*info != 0) {
#line 466 "sgeevx.f"
	i__1 = -(*info);
#line 466 "sgeevx.f"
	xerbla_("SGEEVX", &i__1, (ftnlen)6);
#line 467 "sgeevx.f"
	return 0;
#line 468 "sgeevx.f"
    } else if (lquery) {
#line 469 "sgeevx.f"
	return 0;
#line 470 "sgeevx.f"
    }

/*     Quick return if possible */

#line 474 "sgeevx.f"
    if (*n == 0) {
#line 474 "sgeevx.f"
	return 0;
#line 474 "sgeevx.f"
    }

/*     Get machine constants */

#line 479 "sgeevx.f"
    eps = slamch_("P", (ftnlen)1);
#line 480 "sgeevx.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 481 "sgeevx.f"
    bignum = 1. / smlnum;
#line 482 "sgeevx.f"
    slabad_(&smlnum, &bignum);
#line 483 "sgeevx.f"
    smlnum = sqrt(smlnum) / eps;
#line 484 "sgeevx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 488 "sgeevx.f"
    icond = 0;
#line 489 "sgeevx.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 490 "sgeevx.f"
    scalea = FALSE_;
#line 491 "sgeevx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 492 "sgeevx.f"
	scalea = TRUE_;
#line 493 "sgeevx.f"
	cscale = smlnum;
#line 494 "sgeevx.f"
    } else if (anrm > bignum) {
#line 495 "sgeevx.f"
	scalea = TRUE_;
#line 496 "sgeevx.f"
	cscale = bignum;
#line 497 "sgeevx.f"
    }
#line 498 "sgeevx.f"
    if (scalea) {
#line 498 "sgeevx.f"
	slascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 498 "sgeevx.f"
    }

/*     Balance the matrix and compute ABNRM */

#line 503 "sgeevx.f"
    sgebal_(balanc, n, &a[a_offset], lda, ilo, ihi, &scale[1], &ierr, (ftnlen)
	    1);
#line 504 "sgeevx.f"
    *abnrm = slange_("1", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 505 "sgeevx.f"
    if (scalea) {
#line 506 "sgeevx.f"
	dum[0] = *abnrm;
#line 507 "sgeevx.f"
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &c__1, &
		ierr, (ftnlen)1);
#line 508 "sgeevx.f"
	*abnrm = dum[0];
#line 509 "sgeevx.f"
    }

/*     Reduce to upper Hessenberg form */
/*     (Workspace: need 2*N, prefer N+N*NB) */

#line 514 "sgeevx.f"
    itau = 1;
#line 515 "sgeevx.f"
    iwrk = itau + *n;
#line 516 "sgeevx.f"
    i__1 = *lwork - iwrk + 1;
#line 516 "sgeevx.f"
    sgehrd_(n, ilo, ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &
	    ierr);

#line 519 "sgeevx.f"
    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

#line 524 "sgeevx.f"
	*(unsigned char *)side = 'L';
#line 525 "sgeevx.f"
	slacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VL */
/*        (Workspace: need 2*N-1, prefer N+(N-1)*NB) */

#line 530 "sgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 530 "sgeevx.f"
	sorghr_(n, ilo, ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

#line 536 "sgeevx.f"
	iwrk = itau;
#line 537 "sgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 537 "sgeevx.f"
	shseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vl[
		vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 540 "sgeevx.f"
	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

#line 545 "sgeevx.f"
	    *(unsigned char *)side = 'B';
#line 546 "sgeevx.f"
	    slacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (
		    ftnlen)1);
#line 547 "sgeevx.f"
	}

#line 549 "sgeevx.f"
    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

#line 554 "sgeevx.f"
	*(unsigned char *)side = 'R';
#line 555 "sgeevx.f"
	slacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VR */
/*        (Workspace: need 2*N-1, prefer N+(N-1)*NB) */

#line 560 "sgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 560 "sgeevx.f"
	sorghr_(n, ilo, ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &
		i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

#line 566 "sgeevx.f"
	iwrk = itau;
#line 567 "sgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 567 "sgeevx.f"
	shseqr_("S", "V", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 570 "sgeevx.f"
    } else {

/*        Compute eigenvalues only */
/*        If condition numbers desired, compute Schur form */

#line 575 "sgeevx.f"
	if (wntsnn) {
#line 576 "sgeevx.f"
	    *(unsigned char *)job = 'E';
#line 577 "sgeevx.f"
	} else {
#line 578 "sgeevx.f"
	    *(unsigned char *)job = 'S';
#line 579 "sgeevx.f"
	}

/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

#line 583 "sgeevx.f"
	iwrk = itau;
#line 584 "sgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 584 "sgeevx.f"
	shseqr_(job, "N", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);
#line 586 "sgeevx.f"
    }

/*     If INFO .NE. 0 from SHSEQR, then quit */

#line 590 "sgeevx.f"
    if (*info != 0) {
#line 590 "sgeevx.f"
	goto L50;
#line 590 "sgeevx.f"
    }

#line 593 "sgeevx.f"
    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (Workspace: need 3*N, prefer N + 2*N*NB) */

#line 598 "sgeevx.f"
	i__1 = *lwork - iwrk + 1;
#line 598 "sgeevx.f"
	strevc3_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &i__1, &
		ierr, (ftnlen)1, (ftnlen)1);
#line 600 "sgeevx.f"
    }

/*     Compute condition numbers if desired */
/*     (Workspace: need N*N+6*N unless SENSE = 'E') */

#line 605 "sgeevx.f"
    if (! wntsnn) {
#line 606 "sgeevx.f"
	strsna_(sense, "A", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, &rconde[1], &rcondv[1], n, &nout, 
		&work[iwrk], n, &iwork[1], &icond, (ftnlen)1, (ftnlen)1);
#line 609 "sgeevx.f"
    }

#line 611 "sgeevx.f"
    if (wantvl) {

/*        Undo balancing of left eigenvectors */

#line 615 "sgeevx.f"
	sgebak_(balanc, "L", n, ilo, ihi, &scale[1], n, &vl[vl_offset], ldvl, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize left eigenvectors and make largest component real */

#line 620 "sgeevx.f"
	i__1 = *n;
#line 620 "sgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 621 "sgeevx.f"
	    if (wi[i__] == 0.) {
#line 622 "sgeevx.f"
		scl = 1. / snrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 623 "sgeevx.f"
		sscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 624 "sgeevx.f"
	    } else if (wi[i__] > 0.) {
#line 625 "sgeevx.f"
		d__1 = snrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 625 "sgeevx.f"
		d__2 = snrm2_(n, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
#line 625 "sgeevx.f"
		scl = 1. / slapy2_(&d__1, &d__2);
#line 627 "sgeevx.f"
		sscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 628 "sgeevx.f"
		sscal_(n, &scl, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
#line 629 "sgeevx.f"
		i__2 = *n;
#line 629 "sgeevx.f"
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
#line 630 "sgeevx.f"
		    d__1 = vl[k + i__ * vl_dim1];
/* Computing 2nd power */
#line 630 "sgeevx.f"
		    d__2 = vl[k + (i__ + 1) * vl_dim1];
#line 630 "sgeevx.f"
		    work[k] = d__1 * d__1 + d__2 * d__2;
#line 631 "sgeevx.f"
/* L10: */
#line 631 "sgeevx.f"
		}
#line 632 "sgeevx.f"
		k = isamax_(n, &work[1], &c__1);
#line 633 "sgeevx.f"
		slartg_(&vl[k + i__ * vl_dim1], &vl[k + (i__ + 1) * vl_dim1], 
			&cs, &sn, &r__);
#line 634 "sgeevx.f"
		srot_(n, &vl[i__ * vl_dim1 + 1], &c__1, &vl[(i__ + 1) * 
			vl_dim1 + 1], &c__1, &cs, &sn);
#line 635 "sgeevx.f"
		vl[k + (i__ + 1) * vl_dim1] = 0.;
#line 636 "sgeevx.f"
	    }
#line 637 "sgeevx.f"
/* L20: */
#line 637 "sgeevx.f"
	}
#line 638 "sgeevx.f"
    }

#line 640 "sgeevx.f"
    if (wantvr) {

/*        Undo balancing of right eigenvectors */

#line 644 "sgeevx.f"
	sgebak_(balanc, "R", n, ilo, ihi, &scale[1], n, &vr[vr_offset], ldvr, 
		&ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize right eigenvectors and make largest component real */

#line 649 "sgeevx.f"
	i__1 = *n;
#line 649 "sgeevx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 650 "sgeevx.f"
	    if (wi[i__] == 0.) {
#line 651 "sgeevx.f"
		scl = 1. / snrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 652 "sgeevx.f"
		sscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 653 "sgeevx.f"
	    } else if (wi[i__] > 0.) {
#line 654 "sgeevx.f"
		d__1 = snrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 654 "sgeevx.f"
		d__2 = snrm2_(n, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
#line 654 "sgeevx.f"
		scl = 1. / slapy2_(&d__1, &d__2);
#line 656 "sgeevx.f"
		sscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 657 "sgeevx.f"
		sscal_(n, &scl, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
#line 658 "sgeevx.f"
		i__2 = *n;
#line 658 "sgeevx.f"
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
#line 659 "sgeevx.f"
		    d__1 = vr[k + i__ * vr_dim1];
/* Computing 2nd power */
#line 659 "sgeevx.f"
		    d__2 = vr[k + (i__ + 1) * vr_dim1];
#line 659 "sgeevx.f"
		    work[k] = d__1 * d__1 + d__2 * d__2;
#line 660 "sgeevx.f"
/* L30: */
#line 660 "sgeevx.f"
		}
#line 661 "sgeevx.f"
		k = isamax_(n, &work[1], &c__1);
#line 662 "sgeevx.f"
		slartg_(&vr[k + i__ * vr_dim1], &vr[k + (i__ + 1) * vr_dim1], 
			&cs, &sn, &r__);
#line 663 "sgeevx.f"
		srot_(n, &vr[i__ * vr_dim1 + 1], &c__1, &vr[(i__ + 1) * 
			vr_dim1 + 1], &c__1, &cs, &sn);
#line 664 "sgeevx.f"
		vr[k + (i__ + 1) * vr_dim1] = 0.;
#line 665 "sgeevx.f"
	    }
#line 666 "sgeevx.f"
/* L40: */
#line 666 "sgeevx.f"
	}
#line 667 "sgeevx.f"
    }

/*     Undo scaling if necessary */

#line 671 "sgeevx.f"
L50:
#line 672 "sgeevx.f"
    if (scalea) {
#line 673 "sgeevx.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 673 "sgeevx.f"
	i__3 = *n - *info;
#line 673 "sgeevx.f"
	i__2 = max(i__3,1);
#line 673 "sgeevx.f"
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 
		1], &i__2, &ierr, (ftnlen)1);
#line 675 "sgeevx.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 675 "sgeevx.f"
	i__3 = *n - *info;
#line 675 "sgeevx.f"
	i__2 = max(i__3,1);
#line 675 "sgeevx.f"
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 
		1], &i__2, &ierr, (ftnlen)1);
#line 677 "sgeevx.f"
	if (*info == 0) {
#line 678 "sgeevx.f"
	    if ((wntsnv || wntsnb) && icond == 0) {
#line 678 "sgeevx.f"
		slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &rcondv[
			1], n, &ierr, (ftnlen)1);
#line 678 "sgeevx.f"
	    }
#line 681 "sgeevx.f"
	} else {
#line 682 "sgeevx.f"
	    i__1 = *ilo - 1;
#line 682 "sgeevx.f"
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], 
		    n, &ierr, (ftnlen)1);
#line 684 "sgeevx.f"
	    i__1 = *ilo - 1;
#line 684 "sgeevx.f"
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], 
		    n, &ierr, (ftnlen)1);
#line 686 "sgeevx.f"
	}
#line 687 "sgeevx.f"
    }

#line 689 "sgeevx.f"
    work[1] = (doublereal) maxwrk;
#line 690 "sgeevx.f"
    return 0;

/*     End of SGEEVX */

} /* sgeevx_ */

