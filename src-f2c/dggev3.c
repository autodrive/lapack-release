#line 1 "dggev3.f"
/* dggev3.f -- translated by f2c (version 20100827).
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

#line 1 "dggev3.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b38 = 0.;
static doublereal c_b39 = 1.;

/* > \brief <b> DGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices (blocked algorithm)</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGEV3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggev3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggev3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggev3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, */
/*      $                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, */
/*      $                   INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGEV3 computes for a pair of N-by-N real nonsymmetric matrices (A,B) */
/* > the generalized eigenvalues, and optionally, the left and/or right */
/* > generalized eigenvectors. */
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
/* >                  A * v(j) = lambda(j) * B * v(j). */
/* > */
/* > The left eigenvector u(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* > */
/* >                  u(j)**H * A  = lambda(j) * u(j)**H * B . */
/* > */
/* > where u(j)**H is the conjugate-transpose of u(j). */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* >          On exit, A has been overwritten. */
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
/* >          On exit, B has been overwritten. */
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
/* >          alpha/beta.  However, ALPHAR and ALPHAI will be always less */
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
/* >          Each eigenvector is scaled so the largest component has */
/* >          abs(real part)+abs(imag. part)=1. */
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
/* >          Each eigenvector is scaled so the largest component has */
/* >          abs(real part)+abs(imag. part)=1. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \date January 2015 */

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dggev3_(char *jobvl, char *jobvr, integer *n, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info, ftnlen jobvl_len, ftnlen jobvr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer jc, in, jr, ihi, ilo;
    static doublereal eps;
    static logical ilv;
    static doublereal anrm, bnrm;
    static integer ierr, itau;
    static doublereal temp;
    static logical ilvl, ilvr;
    static integer iwrk;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ileft, icols;
    extern /* Subroutine */ int dgghd3_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer irows;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dggbak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dggbal_(char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), dtgevc_(char *, char *, logical *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static logical ldumma[1];
    static char chtemp[1];
    static doublereal bignum;
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer ijobvl, iright, ijobvr;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal anrmto, bnrmto;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal smlnum;
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     January 2015 */

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

#line 278 "dggev3.f"
    /* Parameter adjustments */
#line 278 "dggev3.f"
    a_dim1 = *lda;
#line 278 "dggev3.f"
    a_offset = 1 + a_dim1;
#line 278 "dggev3.f"
    a -= a_offset;
#line 278 "dggev3.f"
    b_dim1 = *ldb;
#line 278 "dggev3.f"
    b_offset = 1 + b_dim1;
#line 278 "dggev3.f"
    b -= b_offset;
#line 278 "dggev3.f"
    --alphar;
#line 278 "dggev3.f"
    --alphai;
#line 278 "dggev3.f"
    --beta;
#line 278 "dggev3.f"
    vl_dim1 = *ldvl;
#line 278 "dggev3.f"
    vl_offset = 1 + vl_dim1;
#line 278 "dggev3.f"
    vl -= vl_offset;
#line 278 "dggev3.f"
    vr_dim1 = *ldvr;
#line 278 "dggev3.f"
    vr_offset = 1 + vr_dim1;
#line 278 "dggev3.f"
    vr -= vr_offset;
#line 278 "dggev3.f"
    --work;
#line 278 "dggev3.f"

#line 278 "dggev3.f"
    /* Function Body */
#line 278 "dggev3.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 279 "dggev3.f"
	ijobvl = 1;
#line 280 "dggev3.f"
	ilvl = FALSE_;
#line 281 "dggev3.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 282 "dggev3.f"
	ijobvl = 2;
#line 283 "dggev3.f"
	ilvl = TRUE_;
#line 284 "dggev3.f"
    } else {
#line 285 "dggev3.f"
	ijobvl = -1;
#line 286 "dggev3.f"
	ilvl = FALSE_;
#line 287 "dggev3.f"
    }

#line 289 "dggev3.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 290 "dggev3.f"
	ijobvr = 1;
#line 291 "dggev3.f"
	ilvr = FALSE_;
#line 292 "dggev3.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 293 "dggev3.f"
	ijobvr = 2;
#line 294 "dggev3.f"
	ilvr = TRUE_;
#line 295 "dggev3.f"
    } else {
#line 296 "dggev3.f"
	ijobvr = -1;
#line 297 "dggev3.f"
	ilvr = FALSE_;
#line 298 "dggev3.f"
    }
#line 299 "dggev3.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

#line 303 "dggev3.f"
    *info = 0;
#line 304 "dggev3.f"
    lquery = *lwork == -1;
#line 305 "dggev3.f"
    if (ijobvl <= 0) {
#line 306 "dggev3.f"
	*info = -1;
#line 307 "dggev3.f"
    } else if (ijobvr <= 0) {
#line 308 "dggev3.f"
	*info = -2;
#line 309 "dggev3.f"
    } else if (*n < 0) {
#line 310 "dggev3.f"
	*info = -3;
#line 311 "dggev3.f"
    } else if (*lda < max(1,*n)) {
#line 312 "dggev3.f"
	*info = -5;
#line 313 "dggev3.f"
    } else if (*ldb < max(1,*n)) {
#line 314 "dggev3.f"
	*info = -7;
#line 315 "dggev3.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 316 "dggev3.f"
	*info = -12;
#line 317 "dggev3.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 318 "dggev3.f"
	*info = -14;
#line 319 "dggev3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 319 "dggev3.f"
	i__1 = 1, i__2 = *n << 3;
#line 319 "dggev3.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 320 "dggev3.f"
	    *info = -16;
#line 321 "dggev3.f"
	}
#line 321 "dggev3.f"
    }

/*     Compute workspace */

#line 325 "dggev3.f"
    if (*info == 0) {
#line 326 "dggev3.f"
	dgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
/* Computing MAX */
#line 327 "dggev3.f"
	i__1 = 1, i__2 = *n << 3, i__1 = max(i__1,i__2), i__2 = *n * 3 + (
		integer) work[1];
#line 327 "dggev3.f"
	lwkopt = max(i__1,i__2);
#line 328 "dggev3.f"
	dormqr_("L", "T", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], 
		lda, &work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 330 "dggev3.f"
	i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 330 "dggev3.f"
	lwkopt = max(i__1,i__2);
#line 331 "dggev3.f"
	if (ilvl) {
#line 332 "dggev3.f"
	    dorgqr_(n, n, n, &vl[vl_offset], ldvl, &work[1], &work[1], &c_n1, 
		    &ierr);
/* Computing MAX */
#line 333 "dggev3.f"
	    i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 333 "dggev3.f"
	    lwkopt = max(i__1,i__2);
#line 334 "dggev3.f"
	}
#line 335 "dggev3.f"
	if (ilv) {
#line 336 "dggev3.f"
	    dgghd3_(jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[b_offset]
		    , ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[
		    1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 338 "dggev3.f"
	    i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 338 "dggev3.f"
	    lwkopt = max(i__1,i__2);
#line 339 "dggev3.f"
	    dhgeqz_("S", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[
		    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[
		    vl_offset], ldvl, &vr[vr_offset], ldvr, &work[1], &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 342 "dggev3.f"
	    i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 342 "dggev3.f"
	    lwkopt = max(i__1,i__2);
#line 343 "dggev3.f"
	} else {
#line 344 "dggev3.f"
	    dgghd3_("N", "N", n, &c__1, n, &a[a_offset], lda, &b[b_offset], 
		    ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[1],
		     &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 346 "dggev3.f"
	    i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 346 "dggev3.f"
	    lwkopt = max(i__1,i__2);
#line 347 "dggev3.f"
	    dhgeqz_("E", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[
		    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[
		    vl_offset], ldvl, &vr[vr_offset], ldvr, &work[1], &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 350 "dggev3.f"
	    i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 350 "dggev3.f"
	    lwkopt = max(i__1,i__2);
#line 351 "dggev3.f"
	}
#line 353 "dggev3.f"
	work[1] = (doublereal) lwkopt;
#line 354 "dggev3.f"
    }

#line 356 "dggev3.f"
    if (*info != 0) {
#line 357 "dggev3.f"
	i__1 = -(*info);
#line 357 "dggev3.f"
	xerbla_("DGGEV3 ", &i__1, (ftnlen)7);
#line 358 "dggev3.f"
	return 0;
#line 359 "dggev3.f"
    } else if (lquery) {
#line 360 "dggev3.f"
	return 0;
#line 361 "dggev3.f"
    }

/*     Quick return if possible */

#line 365 "dggev3.f"
    if (*n == 0) {
#line 365 "dggev3.f"
	return 0;
#line 365 "dggev3.f"
    }

/*     Get machine constants */

#line 370 "dggev3.f"
    eps = dlamch_("P", (ftnlen)1);
#line 371 "dggev3.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 372 "dggev3.f"
    bignum = 1. / smlnum;
#line 373 "dggev3.f"
    dlabad_(&smlnum, &bignum);
#line 374 "dggev3.f"
    smlnum = sqrt(smlnum) / eps;
#line 375 "dggev3.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 379 "dggev3.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 380 "dggev3.f"
    ilascl = FALSE_;
#line 381 "dggev3.f"
    if (anrm > 0. && anrm < smlnum) {
#line 382 "dggev3.f"
	anrmto = smlnum;
#line 383 "dggev3.f"
	ilascl = TRUE_;
#line 384 "dggev3.f"
    } else if (anrm > bignum) {
#line 385 "dggev3.f"
	anrmto = bignum;
#line 386 "dggev3.f"
	ilascl = TRUE_;
#line 387 "dggev3.f"
    }
#line 388 "dggev3.f"
    if (ilascl) {
#line 388 "dggev3.f"
	dlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 388 "dggev3.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 393 "dggev3.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 394 "dggev3.f"
    ilbscl = FALSE_;
#line 395 "dggev3.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 396 "dggev3.f"
	bnrmto = smlnum;
#line 397 "dggev3.f"
	ilbscl = TRUE_;
#line 398 "dggev3.f"
    } else if (bnrm > bignum) {
#line 399 "dggev3.f"
	bnrmto = bignum;
#line 400 "dggev3.f"
	ilbscl = TRUE_;
#line 401 "dggev3.f"
    }
#line 402 "dggev3.f"
    if (ilbscl) {
#line 402 "dggev3.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 402 "dggev3.f"
    }

/*     Permute the matrices A, B to isolate eigenvalues if possible */

#line 407 "dggev3.f"
    ileft = 1;
#line 408 "dggev3.f"
    iright = *n + 1;
#line 409 "dggev3.f"
    iwrk = iright + *n;
#line 410 "dggev3.f"
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */

#line 415 "dggev3.f"
    irows = ihi + 1 - ilo;
#line 416 "dggev3.f"
    if (ilv) {
#line 417 "dggev3.f"
	icols = *n + 1 - ilo;
#line 418 "dggev3.f"
    } else {
#line 419 "dggev3.f"
	icols = irows;
#line 420 "dggev3.f"
    }
#line 421 "dggev3.f"
    itau = iwrk;
#line 422 "dggev3.f"
    iwrk = itau + irows;
#line 423 "dggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 423 "dggev3.f"
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */

#line 428 "dggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 428 "dggev3.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL */

#line 434 "dggev3.f"
    if (ilvl) {
#line 435 "dggev3.f"
	dlaset_("Full", n, n, &c_b38, &c_b39, &vl[vl_offset], ldvl, (ftnlen)4)
		;
#line 436 "dggev3.f"
	if (irows > 1) {
#line 437 "dggev3.f"
	    i__1 = irows - 1;
#line 437 "dggev3.f"
	    i__2 = irows - 1;
#line 437 "dggev3.f"
	    dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[
		    ilo + 1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 439 "dggev3.f"
	}
#line 440 "dggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 440 "dggev3.f"
	dorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwrk], &i__1, &ierr);
#line 442 "dggev3.f"
    }

/*     Initialize VR */

#line 446 "dggev3.f"
    if (ilvr) {
#line 446 "dggev3.f"
	dlaset_("Full", n, n, &c_b38, &c_b39, &vr[vr_offset], ldvr, (ftnlen)4)
		;
#line 446 "dggev3.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 451 "dggev3.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 455 "dggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 455 "dggev3.f"
	dgghd3_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[iwrk], 
		&i__1, &ierr, (ftnlen)1, (ftnlen)1);
#line 457 "dggev3.f"
    } else {
#line 458 "dggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 458 "dggev3.f"
	dgghd3_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, &ierr, (ftnlen)1, (
		ftnlen)1);
#line 461 "dggev3.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur forms and Schur vectors) */

#line 466 "dggev3.f"
    iwrk = itau;
#line 467 "dggev3.f"
    if (ilv) {
#line 468 "dggev3.f"
	*(unsigned char *)chtemp = 'S';
#line 469 "dggev3.f"
    } else {
#line 470 "dggev3.f"
	*(unsigned char *)chtemp = 'E';
#line 471 "dggev3.f"
    }
#line 472 "dggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 472 "dggev3.f"
    dhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], 
	    ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1, &ierr, (ftnlen)1, 
	    (ftnlen)1, (ftnlen)1);
#line 475 "dggev3.f"
    if (ierr != 0) {
#line 476 "dggev3.f"
	if (ierr > 0 && ierr <= *n) {
#line 477 "dggev3.f"
	    *info = ierr;
#line 478 "dggev3.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 479 "dggev3.f"
	    *info = ierr - *n;
#line 480 "dggev3.f"
	} else {
#line 481 "dggev3.f"
	    *info = *n + 1;
#line 482 "dggev3.f"
	}
#line 483 "dggev3.f"
	goto L110;
#line 484 "dggev3.f"
    }

/*     Compute Eigenvectors */

#line 488 "dggev3.f"
    if (ilv) {
#line 489 "dggev3.f"
	if (ilvl) {
#line 490 "dggev3.f"
	    if (ilvr) {
#line 491 "dggev3.f"
		*(unsigned char *)chtemp = 'B';
#line 492 "dggev3.f"
	    } else {
#line 493 "dggev3.f"
		*(unsigned char *)chtemp = 'L';
#line 494 "dggev3.f"
	    }
#line 495 "dggev3.f"
	} else {
#line 496 "dggev3.f"
	    *(unsigned char *)chtemp = 'R';
#line 497 "dggev3.f"
	}
#line 498 "dggev3.f"
	dtgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwrk], &ierr, (ftnlen)1, (ftnlen)1);
#line 500 "dggev3.f"
	if (ierr != 0) {
#line 501 "dggev3.f"
	    *info = *n + 2;
#line 502 "dggev3.f"
	    goto L110;
#line 503 "dggev3.f"
	}

/*        Undo balancing on VL and VR and normalization */

#line 507 "dggev3.f"
	if (ilvl) {
#line 508 "dggev3.f"
	    dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vl[vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);
#line 510 "dggev3.f"
	    i__1 = *n;
#line 510 "dggev3.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 511 "dggev3.f"
		if (alphai[jc] < 0.) {
#line 511 "dggev3.f"
		    goto L50;
#line 511 "dggev3.f"
		}
#line 513 "dggev3.f"
		temp = 0.;
#line 514 "dggev3.f"
		if (alphai[jc] == 0.) {
#line 515 "dggev3.f"
		    i__2 = *n;
#line 515 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 516 "dggev3.f"
			d__2 = temp, d__3 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1));
#line 516 "dggev3.f"
			temp = max(d__2,d__3);
#line 517 "dggev3.f"
/* L10: */
#line 517 "dggev3.f"
		    }
#line 518 "dggev3.f"
		} else {
#line 519 "dggev3.f"
		    i__2 = *n;
#line 519 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 520 "dggev3.f"
			d__3 = temp, d__4 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1)) + (d__2 = vl[jr + (jc + 1) * 
				vl_dim1], abs(d__2));
#line 520 "dggev3.f"
			temp = max(d__3,d__4);
#line 522 "dggev3.f"
/* L20: */
#line 522 "dggev3.f"
		    }
#line 523 "dggev3.f"
		}
#line 524 "dggev3.f"
		if (temp < smlnum) {
#line 524 "dggev3.f"
		    goto L50;
#line 524 "dggev3.f"
		}
#line 526 "dggev3.f"
		temp = 1. / temp;
#line 527 "dggev3.f"
		if (alphai[jc] == 0.) {
#line 528 "dggev3.f"
		    i__2 = *n;
#line 528 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 529 "dggev3.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 530 "dggev3.f"
/* L30: */
#line 530 "dggev3.f"
		    }
#line 531 "dggev3.f"
		} else {
#line 532 "dggev3.f"
		    i__2 = *n;
#line 532 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 533 "dggev3.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 534 "dggev3.f"
			vl[jr + (jc + 1) * vl_dim1] *= temp;
#line 535 "dggev3.f"
/* L40: */
#line 535 "dggev3.f"
		    }
#line 536 "dggev3.f"
		}
#line 537 "dggev3.f"
L50:
#line 537 "dggev3.f"
		;
#line 537 "dggev3.f"
	    }
#line 538 "dggev3.f"
	}
#line 539 "dggev3.f"
	if (ilvr) {
#line 540 "dggev3.f"
	    dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vr[vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 542 "dggev3.f"
	    i__1 = *n;
#line 542 "dggev3.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 543 "dggev3.f"
		if (alphai[jc] < 0.) {
#line 543 "dggev3.f"
		    goto L100;
#line 543 "dggev3.f"
		}
#line 545 "dggev3.f"
		temp = 0.;
#line 546 "dggev3.f"
		if (alphai[jc] == 0.) {
#line 547 "dggev3.f"
		    i__2 = *n;
#line 547 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 548 "dggev3.f"
			d__2 = temp, d__3 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1));
#line 548 "dggev3.f"
			temp = max(d__2,d__3);
#line 549 "dggev3.f"
/* L60: */
#line 549 "dggev3.f"
		    }
#line 550 "dggev3.f"
		} else {
#line 551 "dggev3.f"
		    i__2 = *n;
#line 551 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 552 "dggev3.f"
			d__3 = temp, d__4 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1)) + (d__2 = vr[jr + (jc + 1) * 
				vr_dim1], abs(d__2));
#line 552 "dggev3.f"
			temp = max(d__3,d__4);
#line 554 "dggev3.f"
/* L70: */
#line 554 "dggev3.f"
		    }
#line 555 "dggev3.f"
		}
#line 556 "dggev3.f"
		if (temp < smlnum) {
#line 556 "dggev3.f"
		    goto L100;
#line 556 "dggev3.f"
		}
#line 558 "dggev3.f"
		temp = 1. / temp;
#line 559 "dggev3.f"
		if (alphai[jc] == 0.) {
#line 560 "dggev3.f"
		    i__2 = *n;
#line 560 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 561 "dggev3.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 562 "dggev3.f"
/* L80: */
#line 562 "dggev3.f"
		    }
#line 563 "dggev3.f"
		} else {
#line 564 "dggev3.f"
		    i__2 = *n;
#line 564 "dggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 565 "dggev3.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 566 "dggev3.f"
			vr[jr + (jc + 1) * vr_dim1] *= temp;
#line 567 "dggev3.f"
/* L90: */
#line 567 "dggev3.f"
		    }
#line 568 "dggev3.f"
		}
#line 569 "dggev3.f"
L100:
#line 569 "dggev3.f"
		;
#line 569 "dggev3.f"
	    }
#line 570 "dggev3.f"
	}

/*        End of eigenvector calculation */

#line 574 "dggev3.f"
    }

/*     Undo scaling if necessary */

#line 578 "dggev3.f"
L110:

#line 580 "dggev3.f"
    if (ilascl) {
#line 581 "dggev3.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 582 "dggev3.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 583 "dggev3.f"
    }

#line 585 "dggev3.f"
    if (ilbscl) {
#line 586 "dggev3.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 587 "dggev3.f"
    }

#line 589 "dggev3.f"
    work[1] = (doublereal) lwkopt;
#line 590 "dggev3.f"
    return 0;

/*     End of DGGEV3 */

} /* dggev3_ */

