#line 1 "dggev.f"
/* dggev.f -- translated by f2c (version 20100827).
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

#line 1 "dggev.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
static doublereal c_b36 = 0.;
static doublereal c_b37 = 1.;

/* > \brief <b> DGGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, */
/*                         BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) */

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
/* > DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B) */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,8*N). */
/* >          For good performance, LWORK must generally be larger. */
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

/* > \date April 2012 */

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dggev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
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
    static integer ileft, icols, irows;
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
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ijobvl, iright, ijobvr;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal anrmto, bnrmto;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static logical lquery;


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

#line 280 "dggev.f"
    /* Parameter adjustments */
#line 280 "dggev.f"
    a_dim1 = *lda;
#line 280 "dggev.f"
    a_offset = 1 + a_dim1;
#line 280 "dggev.f"
    a -= a_offset;
#line 280 "dggev.f"
    b_dim1 = *ldb;
#line 280 "dggev.f"
    b_offset = 1 + b_dim1;
#line 280 "dggev.f"
    b -= b_offset;
#line 280 "dggev.f"
    --alphar;
#line 280 "dggev.f"
    --alphai;
#line 280 "dggev.f"
    --beta;
#line 280 "dggev.f"
    vl_dim1 = *ldvl;
#line 280 "dggev.f"
    vl_offset = 1 + vl_dim1;
#line 280 "dggev.f"
    vl -= vl_offset;
#line 280 "dggev.f"
    vr_dim1 = *ldvr;
#line 280 "dggev.f"
    vr_offset = 1 + vr_dim1;
#line 280 "dggev.f"
    vr -= vr_offset;
#line 280 "dggev.f"
    --work;
#line 280 "dggev.f"

#line 280 "dggev.f"
    /* Function Body */
#line 280 "dggev.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 281 "dggev.f"
	ijobvl = 1;
#line 282 "dggev.f"
	ilvl = FALSE_;
#line 283 "dggev.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 284 "dggev.f"
	ijobvl = 2;
#line 285 "dggev.f"
	ilvl = TRUE_;
#line 286 "dggev.f"
    } else {
#line 287 "dggev.f"
	ijobvl = -1;
#line 288 "dggev.f"
	ilvl = FALSE_;
#line 289 "dggev.f"
    }

#line 291 "dggev.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 292 "dggev.f"
	ijobvr = 1;
#line 293 "dggev.f"
	ilvr = FALSE_;
#line 294 "dggev.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 295 "dggev.f"
	ijobvr = 2;
#line 296 "dggev.f"
	ilvr = TRUE_;
#line 297 "dggev.f"
    } else {
#line 298 "dggev.f"
	ijobvr = -1;
#line 299 "dggev.f"
	ilvr = FALSE_;
#line 300 "dggev.f"
    }
#line 301 "dggev.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

#line 305 "dggev.f"
    *info = 0;
#line 306 "dggev.f"
    lquery = *lwork == -1;
#line 307 "dggev.f"
    if (ijobvl <= 0) {
#line 308 "dggev.f"
	*info = -1;
#line 309 "dggev.f"
    } else if (ijobvr <= 0) {
#line 310 "dggev.f"
	*info = -2;
#line 311 "dggev.f"
    } else if (*n < 0) {
#line 312 "dggev.f"
	*info = -3;
#line 313 "dggev.f"
    } else if (*lda < max(1,*n)) {
#line 314 "dggev.f"
	*info = -5;
#line 315 "dggev.f"
    } else if (*ldb < max(1,*n)) {
#line 316 "dggev.f"
	*info = -7;
#line 317 "dggev.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 318 "dggev.f"
	*info = -12;
#line 319 "dggev.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 320 "dggev.f"
	*info = -14;
#line 321 "dggev.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. The workspace is */
/*       computed assuming ILO = 1 and IHI = N, the worst case.) */

#line 331 "dggev.f"
    if (*info == 0) {
/* Computing MAX */
#line 332 "dggev.f"
	i__1 = 1, i__2 = *n << 3;
#line 332 "dggev.f"
	minwrk = max(i__1,i__2);
/* Computing MAX */
#line 333 "dggev.f"
	i__1 = 1, i__2 = *n * (ilaenv_(&c__1, "DGEQRF", " ", n, &c__1, n, &
		c__0, (ftnlen)6, (ftnlen)1) + 7);
#line 333 "dggev.f"
	maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 335 "dggev.f"
	i__1 = maxwrk, i__2 = *n * (ilaenv_(&c__1, "DORMQR", " ", n, &c__1, n,
		 &c__0, (ftnlen)6, (ftnlen)1) + 7);
#line 335 "dggev.f"
	maxwrk = max(i__1,i__2);
#line 337 "dggev.f"
	if (ilvl) {
/* Computing MAX */
#line 338 "dggev.f"
	    i__1 = maxwrk, i__2 = *n * (ilaenv_(&c__1, "DORGQR", " ", n, &
		    c__1, n, &c_n1, (ftnlen)6, (ftnlen)1) + 7);
#line 338 "dggev.f"
	    maxwrk = max(i__1,i__2);
#line 340 "dggev.f"
	}
#line 341 "dggev.f"
	work[1] = (doublereal) maxwrk;

#line 343 "dggev.f"
	if (*lwork < minwrk && ! lquery) {
#line 343 "dggev.f"
	    *info = -16;
#line 343 "dggev.f"
	}
#line 345 "dggev.f"
    }

#line 347 "dggev.f"
    if (*info != 0) {
#line 348 "dggev.f"
	i__1 = -(*info);
#line 348 "dggev.f"
	xerbla_("DGGEV ", &i__1, (ftnlen)6);
#line 349 "dggev.f"
	return 0;
#line 350 "dggev.f"
    } else if (lquery) {
#line 351 "dggev.f"
	return 0;
#line 352 "dggev.f"
    }

/*     Quick return if possible */

#line 356 "dggev.f"
    if (*n == 0) {
#line 356 "dggev.f"
	return 0;
#line 356 "dggev.f"
    }

/*     Get machine constants */

#line 361 "dggev.f"
    eps = dlamch_("P", (ftnlen)1);
#line 362 "dggev.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 363 "dggev.f"
    bignum = 1. / smlnum;
#line 364 "dggev.f"
    dlabad_(&smlnum, &bignum);
#line 365 "dggev.f"
    smlnum = sqrt(smlnum) / eps;
#line 366 "dggev.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 370 "dggev.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 371 "dggev.f"
    ilascl = FALSE_;
#line 372 "dggev.f"
    if (anrm > 0. && anrm < smlnum) {
#line 373 "dggev.f"
	anrmto = smlnum;
#line 374 "dggev.f"
	ilascl = TRUE_;
#line 375 "dggev.f"
    } else if (anrm > bignum) {
#line 376 "dggev.f"
	anrmto = bignum;
#line 377 "dggev.f"
	ilascl = TRUE_;
#line 378 "dggev.f"
    }
#line 379 "dggev.f"
    if (ilascl) {
#line 379 "dggev.f"
	dlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 379 "dggev.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 384 "dggev.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 385 "dggev.f"
    ilbscl = FALSE_;
#line 386 "dggev.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 387 "dggev.f"
	bnrmto = smlnum;
#line 388 "dggev.f"
	ilbscl = TRUE_;
#line 389 "dggev.f"
    } else if (bnrm > bignum) {
#line 390 "dggev.f"
	bnrmto = bignum;
#line 391 "dggev.f"
	ilbscl = TRUE_;
#line 392 "dggev.f"
    }
#line 393 "dggev.f"
    if (ilbscl) {
#line 393 "dggev.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 393 "dggev.f"
    }

/*     Permute the matrices A, B to isolate eigenvalues if possible */
/*     (Workspace: need 6*N) */

#line 399 "dggev.f"
    ileft = 1;
#line 400 "dggev.f"
    iright = *n + 1;
#line 401 "dggev.f"
    iwrk = iright + *n;
#line 402 "dggev.f"
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Workspace: need N, prefer N*NB) */

#line 408 "dggev.f"
    irows = ihi + 1 - ilo;
#line 409 "dggev.f"
    if (ilv) {
#line 410 "dggev.f"
	icols = *n + 1 - ilo;
#line 411 "dggev.f"
    } else {
#line 412 "dggev.f"
	icols = irows;
#line 413 "dggev.f"
    }
#line 414 "dggev.f"
    itau = iwrk;
#line 415 "dggev.f"
    iwrk = itau + irows;
#line 416 "dggev.f"
    i__1 = *lwork + 1 - iwrk;
#line 416 "dggev.f"
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Workspace: need N, prefer N*NB) */

#line 422 "dggev.f"
    i__1 = *lwork + 1 - iwrk;
#line 422 "dggev.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL */
/*     (Workspace: need N, prefer N*NB) */

#line 429 "dggev.f"
    if (ilvl) {
#line 430 "dggev.f"
	dlaset_("Full", n, n, &c_b36, &c_b37, &vl[vl_offset], ldvl, (ftnlen)4)
		;
#line 431 "dggev.f"
	if (irows > 1) {
#line 432 "dggev.f"
	    i__1 = irows - 1;
#line 432 "dggev.f"
	    i__2 = irows - 1;
#line 432 "dggev.f"
	    dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[
		    ilo + 1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 434 "dggev.f"
	}
#line 435 "dggev.f"
	i__1 = *lwork + 1 - iwrk;
#line 435 "dggev.f"
	dorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwrk], &i__1, &ierr);
#line 437 "dggev.f"
    }

/*     Initialize VR */

#line 441 "dggev.f"
    if (ilvr) {
#line 441 "dggev.f"
	dlaset_("Full", n, n, &c_b36, &c_b37, &vr[vr_offset], ldvr, (ftnlen)4)
		;
#line 441 "dggev.f"
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

#line 447 "dggev.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 451 "dggev.f"
	dgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr, (
		ftnlen)1, (ftnlen)1);
#line 453 "dggev.f"
    } else {
#line 454 "dggev.f"
	dgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 456 "dggev.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur forms and Schur vectors) */
/*     (Workspace: need N) */

#line 462 "dggev.f"
    iwrk = itau;
#line 463 "dggev.f"
    if (ilv) {
#line 464 "dggev.f"
	*(unsigned char *)chtemp = 'S';
#line 465 "dggev.f"
    } else {
#line 466 "dggev.f"
	*(unsigned char *)chtemp = 'E';
#line 467 "dggev.f"
    }
#line 468 "dggev.f"
    i__1 = *lwork + 1 - iwrk;
#line 468 "dggev.f"
    dhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], 
	    ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1, &ierr, (ftnlen)1, 
	    (ftnlen)1, (ftnlen)1);
#line 471 "dggev.f"
    if (ierr != 0) {
#line 472 "dggev.f"
	if (ierr > 0 && ierr <= *n) {
#line 473 "dggev.f"
	    *info = ierr;
#line 474 "dggev.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 475 "dggev.f"
	    *info = ierr - *n;
#line 476 "dggev.f"
	} else {
#line 477 "dggev.f"
	    *info = *n + 1;
#line 478 "dggev.f"
	}
#line 479 "dggev.f"
	goto L110;
#line 480 "dggev.f"
    }

/*     Compute Eigenvectors */
/*     (Workspace: need 6*N) */

#line 485 "dggev.f"
    if (ilv) {
#line 486 "dggev.f"
	if (ilvl) {
#line 487 "dggev.f"
	    if (ilvr) {
#line 488 "dggev.f"
		*(unsigned char *)chtemp = 'B';
#line 489 "dggev.f"
	    } else {
#line 490 "dggev.f"
		*(unsigned char *)chtemp = 'L';
#line 491 "dggev.f"
	    }
#line 492 "dggev.f"
	} else {
#line 493 "dggev.f"
	    *(unsigned char *)chtemp = 'R';
#line 494 "dggev.f"
	}
#line 495 "dggev.f"
	dtgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwrk], &ierr, (ftnlen)1, (ftnlen)1);
#line 497 "dggev.f"
	if (ierr != 0) {
#line 498 "dggev.f"
	    *info = *n + 2;
#line 499 "dggev.f"
	    goto L110;
#line 500 "dggev.f"
	}

/*        Undo balancing on VL and VR and normalization */
/*        (Workspace: none needed) */

#line 505 "dggev.f"
	if (ilvl) {
#line 506 "dggev.f"
	    dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vl[vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);
#line 508 "dggev.f"
	    i__1 = *n;
#line 508 "dggev.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 509 "dggev.f"
		if (alphai[jc] < 0.) {
#line 509 "dggev.f"
		    goto L50;
#line 509 "dggev.f"
		}
#line 511 "dggev.f"
		temp = 0.;
#line 512 "dggev.f"
		if (alphai[jc] == 0.) {
#line 513 "dggev.f"
		    i__2 = *n;
#line 513 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 514 "dggev.f"
			d__2 = temp, d__3 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1));
#line 514 "dggev.f"
			temp = max(d__2,d__3);
#line 515 "dggev.f"
/* L10: */
#line 515 "dggev.f"
		    }
#line 516 "dggev.f"
		} else {
#line 517 "dggev.f"
		    i__2 = *n;
#line 517 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 518 "dggev.f"
			d__3 = temp, d__4 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1)) + (d__2 = vl[jr + (jc + 1) * 
				vl_dim1], abs(d__2));
#line 518 "dggev.f"
			temp = max(d__3,d__4);
#line 520 "dggev.f"
/* L20: */
#line 520 "dggev.f"
		    }
#line 521 "dggev.f"
		}
#line 522 "dggev.f"
		if (temp < smlnum) {
#line 522 "dggev.f"
		    goto L50;
#line 522 "dggev.f"
		}
#line 524 "dggev.f"
		temp = 1. / temp;
#line 525 "dggev.f"
		if (alphai[jc] == 0.) {
#line 526 "dggev.f"
		    i__2 = *n;
#line 526 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 527 "dggev.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 528 "dggev.f"
/* L30: */
#line 528 "dggev.f"
		    }
#line 529 "dggev.f"
		} else {
#line 530 "dggev.f"
		    i__2 = *n;
#line 530 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 531 "dggev.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 532 "dggev.f"
			vl[jr + (jc + 1) * vl_dim1] *= temp;
#line 533 "dggev.f"
/* L40: */
#line 533 "dggev.f"
		    }
#line 534 "dggev.f"
		}
#line 535 "dggev.f"
L50:
#line 535 "dggev.f"
		;
#line 535 "dggev.f"
	    }
#line 536 "dggev.f"
	}
#line 537 "dggev.f"
	if (ilvr) {
#line 538 "dggev.f"
	    dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vr[vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 540 "dggev.f"
	    i__1 = *n;
#line 540 "dggev.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 541 "dggev.f"
		if (alphai[jc] < 0.) {
#line 541 "dggev.f"
		    goto L100;
#line 541 "dggev.f"
		}
#line 543 "dggev.f"
		temp = 0.;
#line 544 "dggev.f"
		if (alphai[jc] == 0.) {
#line 545 "dggev.f"
		    i__2 = *n;
#line 545 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 546 "dggev.f"
			d__2 = temp, d__3 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1));
#line 546 "dggev.f"
			temp = max(d__2,d__3);
#line 547 "dggev.f"
/* L60: */
#line 547 "dggev.f"
		    }
#line 548 "dggev.f"
		} else {
#line 549 "dggev.f"
		    i__2 = *n;
#line 549 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 550 "dggev.f"
			d__3 = temp, d__4 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1)) + (d__2 = vr[jr + (jc + 1) * 
				vr_dim1], abs(d__2));
#line 550 "dggev.f"
			temp = max(d__3,d__4);
#line 552 "dggev.f"
/* L70: */
#line 552 "dggev.f"
		    }
#line 553 "dggev.f"
		}
#line 554 "dggev.f"
		if (temp < smlnum) {
#line 554 "dggev.f"
		    goto L100;
#line 554 "dggev.f"
		}
#line 556 "dggev.f"
		temp = 1. / temp;
#line 557 "dggev.f"
		if (alphai[jc] == 0.) {
#line 558 "dggev.f"
		    i__2 = *n;
#line 558 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 559 "dggev.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 560 "dggev.f"
/* L80: */
#line 560 "dggev.f"
		    }
#line 561 "dggev.f"
		} else {
#line 562 "dggev.f"
		    i__2 = *n;
#line 562 "dggev.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 563 "dggev.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 564 "dggev.f"
			vr[jr + (jc + 1) * vr_dim1] *= temp;
#line 565 "dggev.f"
/* L90: */
#line 565 "dggev.f"
		    }
#line 566 "dggev.f"
		}
#line 567 "dggev.f"
L100:
#line 567 "dggev.f"
		;
#line 567 "dggev.f"
	    }
#line 568 "dggev.f"
	}

/*        End of eigenvector calculation */

#line 572 "dggev.f"
    }

/*     Undo scaling if necessary */

#line 576 "dggev.f"
L110:

#line 578 "dggev.f"
    if (ilascl) {
#line 579 "dggev.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 580 "dggev.f"
	dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 581 "dggev.f"
    }

#line 583 "dggev.f"
    if (ilbscl) {
#line 584 "dggev.f"
	dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 585 "dggev.f"
    }

#line 587 "dggev.f"
    work[1] = (doublereal) maxwrk;
#line 588 "dggev.f"
    return 0;

/*     End of DGGEV */

} /* dggev_ */

