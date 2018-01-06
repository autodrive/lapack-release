#line 1 "sggev3.f"
/* sggev3.f -- translated by f2c (version 20100827).
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

#line 1 "sggev3.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b34 = 0.;
static doublereal c_b35 = 1.;

/* > \brief <b> SGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices (blocked algorithm)</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGEV3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggev3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggev3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggev3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, */
/*      $                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, */
/*      $                   INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGEV3 computes for a pair of N-by-N real nonsymmetric matrices (A,B) */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          B is REAL array, dimension (LDB, N) */
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
/* >          alpha/beta.  However, ALPHAR and ALPHAI will be always less */
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
/* >          VR is REAL array, dimension (LDVR,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* >          > N:  =N+1: other than QZ iteration failed in SHGEQZ. */
/* >                =N+2: error return from STGEVC. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2015 */

/* > \ingroup realGEeigen */

/*  ===================================================================== */
/* Subroutine */ int sggev3_(char *jobvl, char *jobvr, integer *n, doublereal 
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
    static integer ileft, icols, irows;
    extern /* Subroutine */ int sgghd3_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), slabad_(doublereal *, 
	    doublereal *), sggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), sggbal_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ilascl, ilbscl;
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical ldumma[1];
    static char chtemp[1];
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer ijobvl, iright;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer ijobvr;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), stgevc_(char *, char *, logical 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal anrmto, bnrmto;
    extern /* Subroutine */ int shgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int sorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 278 "sggev3.f"
    /* Parameter adjustments */
#line 278 "sggev3.f"
    a_dim1 = *lda;
#line 278 "sggev3.f"
    a_offset = 1 + a_dim1;
#line 278 "sggev3.f"
    a -= a_offset;
#line 278 "sggev3.f"
    b_dim1 = *ldb;
#line 278 "sggev3.f"
    b_offset = 1 + b_dim1;
#line 278 "sggev3.f"
    b -= b_offset;
#line 278 "sggev3.f"
    --alphar;
#line 278 "sggev3.f"
    --alphai;
#line 278 "sggev3.f"
    --beta;
#line 278 "sggev3.f"
    vl_dim1 = *ldvl;
#line 278 "sggev3.f"
    vl_offset = 1 + vl_dim1;
#line 278 "sggev3.f"
    vl -= vl_offset;
#line 278 "sggev3.f"
    vr_dim1 = *ldvr;
#line 278 "sggev3.f"
    vr_offset = 1 + vr_dim1;
#line 278 "sggev3.f"
    vr -= vr_offset;
#line 278 "sggev3.f"
    --work;
#line 278 "sggev3.f"

#line 278 "sggev3.f"
    /* Function Body */
#line 278 "sggev3.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 279 "sggev3.f"
	ijobvl = 1;
#line 280 "sggev3.f"
	ilvl = FALSE_;
#line 281 "sggev3.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 282 "sggev3.f"
	ijobvl = 2;
#line 283 "sggev3.f"
	ilvl = TRUE_;
#line 284 "sggev3.f"
    } else {
#line 285 "sggev3.f"
	ijobvl = -1;
#line 286 "sggev3.f"
	ilvl = FALSE_;
#line 287 "sggev3.f"
    }

#line 289 "sggev3.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 290 "sggev3.f"
	ijobvr = 1;
#line 291 "sggev3.f"
	ilvr = FALSE_;
#line 292 "sggev3.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 293 "sggev3.f"
	ijobvr = 2;
#line 294 "sggev3.f"
	ilvr = TRUE_;
#line 295 "sggev3.f"
    } else {
#line 296 "sggev3.f"
	ijobvr = -1;
#line 297 "sggev3.f"
	ilvr = FALSE_;
#line 298 "sggev3.f"
    }
#line 299 "sggev3.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

#line 303 "sggev3.f"
    *info = 0;
#line 304 "sggev3.f"
    lquery = *lwork == -1;
#line 305 "sggev3.f"
    if (ijobvl <= 0) {
#line 306 "sggev3.f"
	*info = -1;
#line 307 "sggev3.f"
    } else if (ijobvr <= 0) {
#line 308 "sggev3.f"
	*info = -2;
#line 309 "sggev3.f"
    } else if (*n < 0) {
#line 310 "sggev3.f"
	*info = -3;
#line 311 "sggev3.f"
    } else if (*lda < max(1,*n)) {
#line 312 "sggev3.f"
	*info = -5;
#line 313 "sggev3.f"
    } else if (*ldb < max(1,*n)) {
#line 314 "sggev3.f"
	*info = -7;
#line 315 "sggev3.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 316 "sggev3.f"
	*info = -12;
#line 317 "sggev3.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 318 "sggev3.f"
	*info = -14;
#line 319 "sggev3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 319 "sggev3.f"
	i__1 = 1, i__2 = *n << 3;
#line 319 "sggev3.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 320 "sggev3.f"
	    *info = -16;
#line 321 "sggev3.f"
	}
#line 321 "sggev3.f"
    }

/*     Compute workspace */

#line 325 "sggev3.f"
    if (*info == 0) {
#line 326 "sggev3.f"
	sgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
/* Computing MAX */
#line 327 "sggev3.f"
	i__1 = 1, i__2 = *n << 3, i__1 = max(i__1,i__2), i__2 = *n * 3 + (
		integer) work[1];
#line 327 "sggev3.f"
	lwkopt = max(i__1,i__2);
#line 328 "sggev3.f"
	sormqr_("L", "T", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], 
		lda, &work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 330 "sggev3.f"
	i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 330 "sggev3.f"
	lwkopt = max(i__1,i__2);
#line 331 "sggev3.f"
	sgghd3_(jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[1], &
		c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 333 "sggev3.f"
	i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 333 "sggev3.f"
	lwkopt = max(i__1,i__2);
#line 334 "sggev3.f"
	if (ilvl) {
#line 335 "sggev3.f"
	    sorgqr_(n, n, n, &vl[vl_offset], ldvl, &work[1], &work[1], &c_n1, 
		    &ierr);
/* Computing MAX */
#line 336 "sggev3.f"
	    i__1 = lwkopt, i__2 = *n * 3 + (integer) work[1];
#line 336 "sggev3.f"
	    lwkopt = max(i__1,i__2);
#line 337 "sggev3.f"
	    shgeqz_("S", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[
		    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[
		    vl_offset], ldvl, &vr[vr_offset], ldvr, &work[1], &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 340 "sggev3.f"
	    i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 340 "sggev3.f"
	    lwkopt = max(i__1,i__2);
#line 341 "sggev3.f"
	} else {
#line 342 "sggev3.f"
	    shgeqz_("E", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[
		    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[
		    vl_offset], ldvl, &vr[vr_offset], ldvr, &work[1], &c_n1, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 345 "sggev3.f"
	    i__1 = lwkopt, i__2 = (*n << 1) + (integer) work[1];
#line 345 "sggev3.f"
	    lwkopt = max(i__1,i__2);
#line 346 "sggev3.f"
	}
#line 347 "sggev3.f"
	work[1] = (doublereal) lwkopt;

#line 349 "sggev3.f"
    }

#line 351 "sggev3.f"
    if (*info != 0) {
#line 352 "sggev3.f"
	i__1 = -(*info);
#line 352 "sggev3.f"
	xerbla_("SGGEV3 ", &i__1, (ftnlen)7);
#line 353 "sggev3.f"
	return 0;
#line 354 "sggev3.f"
    } else if (lquery) {
#line 355 "sggev3.f"
	return 0;
#line 356 "sggev3.f"
    }

/*     Quick return if possible */

#line 360 "sggev3.f"
    if (*n == 0) {
#line 360 "sggev3.f"
	return 0;
#line 360 "sggev3.f"
    }

/*     Get machine constants */

#line 365 "sggev3.f"
    eps = slamch_("P", (ftnlen)1);
#line 366 "sggev3.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 367 "sggev3.f"
    bignum = 1. / smlnum;
#line 368 "sggev3.f"
    slabad_(&smlnum, &bignum);
#line 369 "sggev3.f"
    smlnum = sqrt(smlnum) / eps;
#line 370 "sggev3.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 374 "sggev3.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 375 "sggev3.f"
    ilascl = FALSE_;
#line 376 "sggev3.f"
    if (anrm > 0. && anrm < smlnum) {
#line 377 "sggev3.f"
	anrmto = smlnum;
#line 378 "sggev3.f"
	ilascl = TRUE_;
#line 379 "sggev3.f"
    } else if (anrm > bignum) {
#line 380 "sggev3.f"
	anrmto = bignum;
#line 381 "sggev3.f"
	ilascl = TRUE_;
#line 382 "sggev3.f"
    }
#line 383 "sggev3.f"
    if (ilascl) {
#line 383 "sggev3.f"
	slascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 383 "sggev3.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 388 "sggev3.f"
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 389 "sggev3.f"
    ilbscl = FALSE_;
#line 390 "sggev3.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 391 "sggev3.f"
	bnrmto = smlnum;
#line 392 "sggev3.f"
	ilbscl = TRUE_;
#line 393 "sggev3.f"
    } else if (bnrm > bignum) {
#line 394 "sggev3.f"
	bnrmto = bignum;
#line 395 "sggev3.f"
	ilbscl = TRUE_;
#line 396 "sggev3.f"
    }
#line 397 "sggev3.f"
    if (ilbscl) {
#line 397 "sggev3.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 397 "sggev3.f"
    }

/*     Permute the matrices A, B to isolate eigenvalues if possible */

#line 402 "sggev3.f"
    ileft = 1;
#line 403 "sggev3.f"
    iright = *n + 1;
#line 404 "sggev3.f"
    iwrk = iright + *n;
#line 405 "sggev3.f"
    sggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */

#line 410 "sggev3.f"
    irows = ihi + 1 - ilo;
#line 411 "sggev3.f"
    if (ilv) {
#line 412 "sggev3.f"
	icols = *n + 1 - ilo;
#line 413 "sggev3.f"
    } else {
#line 414 "sggev3.f"
	icols = irows;
#line 415 "sggev3.f"
    }
#line 416 "sggev3.f"
    itau = iwrk;
#line 417 "sggev3.f"
    iwrk = itau + irows;
#line 418 "sggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 418 "sggev3.f"
    sgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */

#line 423 "sggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 423 "sggev3.f"
    sormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL */

#line 429 "sggev3.f"
    if (ilvl) {
#line 430 "sggev3.f"
	slaset_("Full", n, n, &c_b34, &c_b35, &vl[vl_offset], ldvl, (ftnlen)4)
		;
#line 431 "sggev3.f"
	if (irows > 1) {
#line 432 "sggev3.f"
	    i__1 = irows - 1;
#line 432 "sggev3.f"
	    i__2 = irows - 1;
#line 432 "sggev3.f"
	    slacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[
		    ilo + 1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 434 "sggev3.f"
	}
#line 435 "sggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 435 "sggev3.f"
	sorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwrk], &i__1, &ierr);
#line 437 "sggev3.f"
    }

/*     Initialize VR */

#line 441 "sggev3.f"
    if (ilvr) {
#line 441 "sggev3.f"
	slaset_("Full", n, n, &c_b34, &c_b35, &vr[vr_offset], ldvr, (ftnlen)4)
		;
#line 441 "sggev3.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 446 "sggev3.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 450 "sggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 450 "sggev3.f"
	sgghd3_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[iwrk], 
		&i__1, &ierr, (ftnlen)1, (ftnlen)1);
#line 452 "sggev3.f"
    } else {
#line 453 "sggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 453 "sggev3.f"
	sgghd3_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, &ierr, (ftnlen)1, (
		ftnlen)1);
#line 456 "sggev3.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur forms and Schur vectors) */

#line 461 "sggev3.f"
    iwrk = itau;
#line 462 "sggev3.f"
    if (ilv) {
#line 463 "sggev3.f"
	*(unsigned char *)chtemp = 'S';
#line 464 "sggev3.f"
    } else {
#line 465 "sggev3.f"
	*(unsigned char *)chtemp = 'E';
#line 466 "sggev3.f"
    }
#line 467 "sggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 467 "sggev3.f"
    shgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], 
	    ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1, &ierr, (ftnlen)1, 
	    (ftnlen)1, (ftnlen)1);
#line 470 "sggev3.f"
    if (ierr != 0) {
#line 471 "sggev3.f"
	if (ierr > 0 && ierr <= *n) {
#line 472 "sggev3.f"
	    *info = ierr;
#line 473 "sggev3.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 474 "sggev3.f"
	    *info = ierr - *n;
#line 475 "sggev3.f"
	} else {
#line 476 "sggev3.f"
	    *info = *n + 1;
#line 477 "sggev3.f"
	}
#line 478 "sggev3.f"
	goto L110;
#line 479 "sggev3.f"
    }

/*     Compute Eigenvectors */

#line 483 "sggev3.f"
    if (ilv) {
#line 484 "sggev3.f"
	if (ilvl) {
#line 485 "sggev3.f"
	    if (ilvr) {
#line 486 "sggev3.f"
		*(unsigned char *)chtemp = 'B';
#line 487 "sggev3.f"
	    } else {
#line 488 "sggev3.f"
		*(unsigned char *)chtemp = 'L';
#line 489 "sggev3.f"
	    }
#line 490 "sggev3.f"
	} else {
#line 491 "sggev3.f"
	    *(unsigned char *)chtemp = 'R';
#line 492 "sggev3.f"
	}
#line 493 "sggev3.f"
	stgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwrk], &ierr, (ftnlen)1, (ftnlen)1);
#line 495 "sggev3.f"
	if (ierr != 0) {
#line 496 "sggev3.f"
	    *info = *n + 2;
#line 497 "sggev3.f"
	    goto L110;
#line 498 "sggev3.f"
	}

/*        Undo balancing on VL and VR and normalization */

#line 502 "sggev3.f"
	if (ilvl) {
#line 503 "sggev3.f"
	    sggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vl[vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);
#line 505 "sggev3.f"
	    i__1 = *n;
#line 505 "sggev3.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 506 "sggev3.f"
		if (alphai[jc] < 0.) {
#line 506 "sggev3.f"
		    goto L50;
#line 506 "sggev3.f"
		}
#line 508 "sggev3.f"
		temp = 0.;
#line 509 "sggev3.f"
		if (alphai[jc] == 0.) {
#line 510 "sggev3.f"
		    i__2 = *n;
#line 510 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 511 "sggev3.f"
			d__2 = temp, d__3 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1));
#line 511 "sggev3.f"
			temp = max(d__2,d__3);
#line 512 "sggev3.f"
/* L10: */
#line 512 "sggev3.f"
		    }
#line 513 "sggev3.f"
		} else {
#line 514 "sggev3.f"
		    i__2 = *n;
#line 514 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 515 "sggev3.f"
			d__3 = temp, d__4 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1)) + (d__2 = vl[jr + (jc + 1) * 
				vl_dim1], abs(d__2));
#line 515 "sggev3.f"
			temp = max(d__3,d__4);
#line 517 "sggev3.f"
/* L20: */
#line 517 "sggev3.f"
		    }
#line 518 "sggev3.f"
		}
#line 519 "sggev3.f"
		if (temp < smlnum) {
#line 519 "sggev3.f"
		    goto L50;
#line 519 "sggev3.f"
		}
#line 521 "sggev3.f"
		temp = 1. / temp;
#line 522 "sggev3.f"
		if (alphai[jc] == 0.) {
#line 523 "sggev3.f"
		    i__2 = *n;
#line 523 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 524 "sggev3.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 525 "sggev3.f"
/* L30: */
#line 525 "sggev3.f"
		    }
#line 526 "sggev3.f"
		} else {
#line 527 "sggev3.f"
		    i__2 = *n;
#line 527 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 528 "sggev3.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 529 "sggev3.f"
			vl[jr + (jc + 1) * vl_dim1] *= temp;
#line 530 "sggev3.f"
/* L40: */
#line 530 "sggev3.f"
		    }
#line 531 "sggev3.f"
		}
#line 532 "sggev3.f"
L50:
#line 532 "sggev3.f"
		;
#line 532 "sggev3.f"
	    }
#line 533 "sggev3.f"
	}
#line 534 "sggev3.f"
	if (ilvr) {
#line 535 "sggev3.f"
	    sggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vr[vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 537 "sggev3.f"
	    i__1 = *n;
#line 537 "sggev3.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 538 "sggev3.f"
		if (alphai[jc] < 0.) {
#line 538 "sggev3.f"
		    goto L100;
#line 538 "sggev3.f"
		}
#line 540 "sggev3.f"
		temp = 0.;
#line 541 "sggev3.f"
		if (alphai[jc] == 0.) {
#line 542 "sggev3.f"
		    i__2 = *n;
#line 542 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 543 "sggev3.f"
			d__2 = temp, d__3 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1));
#line 543 "sggev3.f"
			temp = max(d__2,d__3);
#line 544 "sggev3.f"
/* L60: */
#line 544 "sggev3.f"
		    }
#line 545 "sggev3.f"
		} else {
#line 546 "sggev3.f"
		    i__2 = *n;
#line 546 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 547 "sggev3.f"
			d__3 = temp, d__4 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1)) + (d__2 = vr[jr + (jc + 1) * 
				vr_dim1], abs(d__2));
#line 547 "sggev3.f"
			temp = max(d__3,d__4);
#line 549 "sggev3.f"
/* L70: */
#line 549 "sggev3.f"
		    }
#line 550 "sggev3.f"
		}
#line 551 "sggev3.f"
		if (temp < smlnum) {
#line 551 "sggev3.f"
		    goto L100;
#line 551 "sggev3.f"
		}
#line 553 "sggev3.f"
		temp = 1. / temp;
#line 554 "sggev3.f"
		if (alphai[jc] == 0.) {
#line 555 "sggev3.f"
		    i__2 = *n;
#line 555 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 556 "sggev3.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 557 "sggev3.f"
/* L80: */
#line 557 "sggev3.f"
		    }
#line 558 "sggev3.f"
		} else {
#line 559 "sggev3.f"
		    i__2 = *n;
#line 559 "sggev3.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 560 "sggev3.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 561 "sggev3.f"
			vr[jr + (jc + 1) * vr_dim1] *= temp;
#line 562 "sggev3.f"
/* L90: */
#line 562 "sggev3.f"
		    }
#line 563 "sggev3.f"
		}
#line 564 "sggev3.f"
L100:
#line 564 "sggev3.f"
		;
#line 564 "sggev3.f"
	    }
#line 565 "sggev3.f"
	}

/*        End of eigenvector calculation */

#line 569 "sggev3.f"
    }

/*     Undo scaling if necessary */

#line 573 "sggev3.f"
L110:

#line 575 "sggev3.f"
    if (ilascl) {
#line 576 "sggev3.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr, (ftnlen)1);
#line 577 "sggev3.f"
	slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr, (ftnlen)1);
#line 578 "sggev3.f"
    }

#line 580 "sggev3.f"
    if (ilbscl) {
#line 581 "sggev3.f"
	slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 582 "sggev3.f"
    }

#line 584 "sggev3.f"
    work[1] = (doublereal) lwkopt;
#line 585 "sggev3.f"
    return 0;

/*     End of SGGEV3 */

} /* sggev3_ */

