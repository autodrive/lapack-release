#line 1 "zggev3.f"
/* zggev3.f -- translated by f2c (version 20100827).
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

#line 1 "zggev3.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;

/* > \brief <b> ZGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices (blocked algorithm)</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGEV3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggev3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggev3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggev3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, */
/*                          VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGEV3 computes for a pair of N-by-N complex nonsymmetric matrices */
/* > (A,B), the generalized eigenvalues, and optionally, the left and/or */
/* > right generalized eigenvectors. */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar */
/* > lambda or a ratio alpha/beta = lambda, such that A - lambda*B is */
/* > singular. It is usually represented as the pair (alpha,beta), as */
/* > there is a reasonable interpretation for beta=0, and even for both */
/* > being zero. */
/* > */
/* > The right generalized eigenvector v(j) corresponding to the */
/* > generalized eigenvalue lambda(j) of (A,B) satisfies */
/* > */
/* >              A * v(j) = lambda(j) * B * v(j). */
/* > */
/* > The left generalized eigenvector u(j) corresponding to the */
/* > generalized eigenvalues lambda(j) of (A,B) satisfies */
/* > */
/* >              u(j)**H * A = lambda(j) * u(j)**H * B */
/* > */
/* > where u(j)**H is the conjugate-transpose of u(j). */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB, N) */
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
/* > \param[out] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 array, dimension (N) */
/* >          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the */
/* >          generalized eigenvalues. */
/* > */
/* >          Note: the quotients ALPHA(j)/BETA(j) may easily over- or */
/* >          underflow, and BETA(j) may even be zero.  Thus, the user */
/* >          should avoid naively computing the ratio alpha/beta. */
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
/* >          Each eigenvector is scaled so the largest component has */
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
/* >          VR is COMPLEX*16 array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right generalized eigenvectors v(j) are */
/* >          stored one after another in the columns of VR, in the same */
/* >          order as their eigenvalues. */
/* >          Each eigenvector is scaled so the largest component has */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (8*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          =1,...,N: */
/* >                The QZ iteration failed.  No eigenvectors have been */
/* >                calculated, but ALPHA(j) and BETA(j) should be */
/* >                correct for j=INFO+1,...,N. */
/* >          > N:  =N+1: other then QZ iteration failed in DHGEQZ, */
/* >                =N+2: error return from DTGEVC. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2015 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zggev3_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info, ftnlen jobvl_len, ftnlen 
	jobvr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);

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
    static integer ileft, icols, irwrk, irows;
    extern /* Subroutine */ int zgghd3_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), zggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical ldumma[1];
    static char chtemp[1];
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static integer ijobvl, iright;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static doublereal anrmto, bnrmto;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), ztgevc_(
	    char *, char *, logical *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     doublereal *, integer *, ftnlen, ftnlen), zhgeqz_(char *, char *,
	     char *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal smlnum;
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.6.1) -- */
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the input arguments */

#line 280 "zggev3.f"
    /* Parameter adjustments */
#line 280 "zggev3.f"
    a_dim1 = *lda;
#line 280 "zggev3.f"
    a_offset = 1 + a_dim1;
#line 280 "zggev3.f"
    a -= a_offset;
#line 280 "zggev3.f"
    b_dim1 = *ldb;
#line 280 "zggev3.f"
    b_offset = 1 + b_dim1;
#line 280 "zggev3.f"
    b -= b_offset;
#line 280 "zggev3.f"
    --alpha;
#line 280 "zggev3.f"
    --beta;
#line 280 "zggev3.f"
    vl_dim1 = *ldvl;
#line 280 "zggev3.f"
    vl_offset = 1 + vl_dim1;
#line 280 "zggev3.f"
    vl -= vl_offset;
#line 280 "zggev3.f"
    vr_dim1 = *ldvr;
#line 280 "zggev3.f"
    vr_offset = 1 + vr_dim1;
#line 280 "zggev3.f"
    vr -= vr_offset;
#line 280 "zggev3.f"
    --work;
#line 280 "zggev3.f"
    --rwork;
#line 280 "zggev3.f"

#line 280 "zggev3.f"
    /* Function Body */
#line 280 "zggev3.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 281 "zggev3.f"
	ijobvl = 1;
#line 282 "zggev3.f"
	ilvl = FALSE_;
#line 283 "zggev3.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 284 "zggev3.f"
	ijobvl = 2;
#line 285 "zggev3.f"
	ilvl = TRUE_;
#line 286 "zggev3.f"
    } else {
#line 287 "zggev3.f"
	ijobvl = -1;
#line 288 "zggev3.f"
	ilvl = FALSE_;
#line 289 "zggev3.f"
    }

#line 291 "zggev3.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 292 "zggev3.f"
	ijobvr = 1;
#line 293 "zggev3.f"
	ilvr = FALSE_;
#line 294 "zggev3.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 295 "zggev3.f"
	ijobvr = 2;
#line 296 "zggev3.f"
	ilvr = TRUE_;
#line 297 "zggev3.f"
    } else {
#line 298 "zggev3.f"
	ijobvr = -1;
#line 299 "zggev3.f"
	ilvr = FALSE_;
#line 300 "zggev3.f"
    }
#line 301 "zggev3.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

#line 305 "zggev3.f"
    *info = 0;
#line 306 "zggev3.f"
    lquery = *lwork == -1;
#line 307 "zggev3.f"
    if (ijobvl <= 0) {
#line 308 "zggev3.f"
	*info = -1;
#line 309 "zggev3.f"
    } else if (ijobvr <= 0) {
#line 310 "zggev3.f"
	*info = -2;
#line 311 "zggev3.f"
    } else if (*n < 0) {
#line 312 "zggev3.f"
	*info = -3;
#line 313 "zggev3.f"
    } else if (*lda < max(1,*n)) {
#line 314 "zggev3.f"
	*info = -5;
#line 315 "zggev3.f"
    } else if (*ldb < max(1,*n)) {
#line 316 "zggev3.f"
	*info = -7;
#line 317 "zggev3.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 318 "zggev3.f"
	*info = -11;
#line 319 "zggev3.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 320 "zggev3.f"
	*info = -13;
#line 321 "zggev3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 321 "zggev3.f"
	i__1 = 1, i__2 = *n << 1;
#line 321 "zggev3.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 322 "zggev3.f"
	    *info = -15;
#line 323 "zggev3.f"
	}
#line 323 "zggev3.f"
    }

/*     Compute workspace */

#line 327 "zggev3.f"
    if (*info == 0) {
#line 328 "zggev3.f"
	zgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
/* Computing MAX */
#line 329 "zggev3.f"
	i__1 = 1, i__2 = *n + (integer) work[1].r;
#line 329 "zggev3.f"
	lwkopt = max(i__1,i__2);
#line 330 "zggev3.f"
	zunmqr_("L", "C", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], 
		lda, &work[1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 332 "zggev3.f"
	i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 332 "zggev3.f"
	lwkopt = max(i__1,i__2);
#line 333 "zggev3.f"
	if (ilvl) {
#line 334 "zggev3.f"
	    zungqr_(n, n, n, &vl[vl_offset], ldvl, &work[1], &work[1], &c_n1, 
		    &ierr);
/* Computing MAX */
#line 335 "zggev3.f"
	    i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 335 "zggev3.f"
	    lwkopt = max(i__1,i__2);
#line 336 "zggev3.f"
	}
#line 337 "zggev3.f"
	if (ilv) {
#line 338 "zggev3.f"
	    zgghd3_(jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[b_offset]
		    , ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[
		    1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 340 "zggev3.f"
	    i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 340 "zggev3.f"
	    lwkopt = max(i__1,i__2);
#line 341 "zggev3.f"
	    zhgeqz_("S", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[
		    b_offset], ldb, &alpha[1], &beta[1], &vl[vl_offset], ldvl,
		     &vr[vr_offset], ldvr, &work[1], &c_n1, &rwork[1], &ierr, 
		    (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 344 "zggev3.f"
	    i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 344 "zggev3.f"
	    lwkopt = max(i__1,i__2);
#line 345 "zggev3.f"
	} else {
#line 346 "zggev3.f"
	    zgghd3_(jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[b_offset]
		    , ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[
		    1], &c_n1, &ierr, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 348 "zggev3.f"
	    i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 348 "zggev3.f"
	    lwkopt = max(i__1,i__2);
#line 349 "zggev3.f"
	    zhgeqz_("E", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[
		    b_offset], ldb, &alpha[1], &beta[1], &vl[vl_offset], ldvl,
		     &vr[vr_offset], ldvr, &work[1], &c_n1, &rwork[1], &ierr, 
		    (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 352 "zggev3.f"
	    i__1 = lwkopt, i__2 = *n + (integer) work[1].r;
#line 352 "zggev3.f"
	    lwkopt = max(i__1,i__2);
#line 353 "zggev3.f"
	}
#line 354 "zggev3.f"
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 354 "zggev3.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 355 "zggev3.f"
    }

#line 357 "zggev3.f"
    if (*info != 0) {
#line 358 "zggev3.f"
	i__1 = -(*info);
#line 358 "zggev3.f"
	xerbla_("ZGGEV3 ", &i__1, (ftnlen)7);
#line 359 "zggev3.f"
	return 0;
#line 360 "zggev3.f"
    } else if (lquery) {
#line 361 "zggev3.f"
	return 0;
#line 362 "zggev3.f"
    }

/*     Quick return if possible */

#line 366 "zggev3.f"
    if (*n == 0) {
#line 366 "zggev3.f"
	return 0;
#line 366 "zggev3.f"
    }

/*     Get machine constants */

#line 371 "zggev3.f"
    eps = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 372 "zggev3.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 373 "zggev3.f"
    bignum = 1. / smlnum;
#line 374 "zggev3.f"
    dlabad_(&smlnum, &bignum);
#line 375 "zggev3.f"
    smlnum = sqrt(smlnum) / eps;
#line 376 "zggev3.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 380 "zggev3.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 381 "zggev3.f"
    ilascl = FALSE_;
#line 382 "zggev3.f"
    if (anrm > 0. && anrm < smlnum) {
#line 383 "zggev3.f"
	anrmto = smlnum;
#line 384 "zggev3.f"
	ilascl = TRUE_;
#line 385 "zggev3.f"
    } else if (anrm > bignum) {
#line 386 "zggev3.f"
	anrmto = bignum;
#line 387 "zggev3.f"
	ilascl = TRUE_;
#line 388 "zggev3.f"
    }
#line 389 "zggev3.f"
    if (ilascl) {
#line 389 "zggev3.f"
	zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 389 "zggev3.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 394 "zggev3.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 395 "zggev3.f"
    ilbscl = FALSE_;
#line 396 "zggev3.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 397 "zggev3.f"
	bnrmto = smlnum;
#line 398 "zggev3.f"
	ilbscl = TRUE_;
#line 399 "zggev3.f"
    } else if (bnrm > bignum) {
#line 400 "zggev3.f"
	bnrmto = bignum;
#line 401 "zggev3.f"
	ilbscl = TRUE_;
#line 402 "zggev3.f"
    }
#line 403 "zggev3.f"
    if (ilbscl) {
#line 403 "zggev3.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 403 "zggev3.f"
    }

/*     Permute the matrices A, B to isolate eigenvalues if possible */

#line 408 "zggev3.f"
    ileft = 1;
#line 409 "zggev3.f"
    iright = *n + 1;
#line 410 "zggev3.f"
    irwrk = iright + *n;
#line 411 "zggev3.f"
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */

#line 416 "zggev3.f"
    irows = ihi + 1 - ilo;
#line 417 "zggev3.f"
    if (ilv) {
#line 418 "zggev3.f"
	icols = *n + 1 - ilo;
#line 419 "zggev3.f"
    } else {
#line 420 "zggev3.f"
	icols = irows;
#line 421 "zggev3.f"
    }
#line 422 "zggev3.f"
    itau = 1;
#line 423 "zggev3.f"
    iwrk = itau + irows;
#line 424 "zggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 424 "zggev3.f"
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */

#line 429 "zggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 429 "zggev3.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL */

#line 435 "zggev3.f"
    if (ilvl) {
#line 436 "zggev3.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vl[vl_offset], ldvl, (ftnlen)4);
#line 437 "zggev3.f"
	if (irows > 1) {
#line 438 "zggev3.f"
	    i__1 = irows - 1;
#line 438 "zggev3.f"
	    i__2 = irows - 1;
#line 438 "zggev3.f"
	    zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[
		    ilo + 1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 440 "zggev3.f"
	}
#line 441 "zggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 441 "zggev3.f"
	zungqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwrk], &i__1, &ierr);
#line 443 "zggev3.f"
    }

/*     Initialize VR */

#line 447 "zggev3.f"
    if (ilvr) {
#line 447 "zggev3.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vr[vr_offset], ldvr, (ftnlen)4);
#line 447 "zggev3.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 452 "zggev3.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 456 "zggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 456 "zggev3.f"
	zgghd3_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[iwrk], 
		&i__1, &ierr, (ftnlen)1, (ftnlen)1);
#line 458 "zggev3.f"
    } else {
#line 459 "zggev3.f"
	i__1 = *lwork + 1 - iwrk;
#line 459 "zggev3.f"
	zgghd3_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, &ierr, (ftnlen)1, (
		ftnlen)1);
#line 462 "zggev3.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur form and Schur vectors) */

#line 467 "zggev3.f"
    iwrk = itau;
#line 468 "zggev3.f"
    if (ilv) {
#line 469 "zggev3.f"
	*(unsigned char *)chtemp = 'S';
#line 470 "zggev3.f"
    } else {
#line 471 "zggev3.f"
	*(unsigned char *)chtemp = 'E';
#line 472 "zggev3.f"
    }
#line 473 "zggev3.f"
    i__1 = *lwork + 1 - iwrk;
#line 473 "zggev3.f"
    zhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vl[vl_offset], ldvl, &vr[
	    vr_offset], ldvr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 476 "zggev3.f"
    if (ierr != 0) {
#line 477 "zggev3.f"
	if (ierr > 0 && ierr <= *n) {
#line 478 "zggev3.f"
	    *info = ierr;
#line 479 "zggev3.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 480 "zggev3.f"
	    *info = ierr - *n;
#line 481 "zggev3.f"
	} else {
#line 482 "zggev3.f"
	    *info = *n + 1;
#line 483 "zggev3.f"
	}
#line 484 "zggev3.f"
	goto L70;
#line 485 "zggev3.f"
    }

/*     Compute Eigenvectors */

#line 489 "zggev3.f"
    if (ilv) {
#line 490 "zggev3.f"
	if (ilvl) {
#line 491 "zggev3.f"
	    if (ilvr) {
#line 492 "zggev3.f"
		*(unsigned char *)chtemp = 'B';
#line 493 "zggev3.f"
	    } else {
#line 494 "zggev3.f"
		*(unsigned char *)chtemp = 'L';
#line 495 "zggev3.f"
	    }
#line 496 "zggev3.f"
	} else {
#line 497 "zggev3.f"
	    *(unsigned char *)chtemp = 'R';
#line 498 "zggev3.f"
	}

#line 500 "zggev3.f"
	ztgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwrk], &rwork[irwrk], &ierr, (ftnlen)1, (ftnlen)1);
#line 503 "zggev3.f"
	if (ierr != 0) {
#line 504 "zggev3.f"
	    *info = *n + 2;
#line 505 "zggev3.f"
	    goto L70;
#line 506 "zggev3.f"
	}

/*        Undo balancing on VL and VR and normalization */

#line 510 "zggev3.f"
	if (ilvl) {
#line 511 "zggev3.f"
	    zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
		     &vl[vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);
#line 513 "zggev3.f"
	    i__1 = *n;
#line 513 "zggev3.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 514 "zggev3.f"
		temp = 0.;
#line 515 "zggev3.f"
		i__2 = *n;
#line 515 "zggev3.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 516 "zggev3.f"
		    i__3 = jr + jc * vl_dim1;
#line 516 "zggev3.f"
		    d__3 = temp, d__4 = (d__1 = vl[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&vl[jr + jc * vl_dim1]), abs(d__2));
#line 516 "zggev3.f"
		    temp = max(d__3,d__4);
#line 517 "zggev3.f"
/* L10: */
#line 517 "zggev3.f"
		}
#line 518 "zggev3.f"
		if (temp < smlnum) {
#line 518 "zggev3.f"
		    goto L30;
#line 518 "zggev3.f"
		}
#line 520 "zggev3.f"
		temp = 1. / temp;
#line 521 "zggev3.f"
		i__2 = *n;
#line 521 "zggev3.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 522 "zggev3.f"
		    i__3 = jr + jc * vl_dim1;
#line 522 "zggev3.f"
		    i__4 = jr + jc * vl_dim1;
#line 522 "zggev3.f"
		    z__1.r = temp * vl[i__4].r, z__1.i = temp * vl[i__4].i;
#line 522 "zggev3.f"
		    vl[i__3].r = z__1.r, vl[i__3].i = z__1.i;
#line 523 "zggev3.f"
/* L20: */
#line 523 "zggev3.f"
		}
#line 524 "zggev3.f"
L30:
#line 524 "zggev3.f"
		;
#line 524 "zggev3.f"
	    }
#line 525 "zggev3.f"
	}
#line 526 "zggev3.f"
	if (ilvr) {
#line 527 "zggev3.f"
	    zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
		     &vr[vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 529 "zggev3.f"
	    i__1 = *n;
#line 529 "zggev3.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 530 "zggev3.f"
		temp = 0.;
#line 531 "zggev3.f"
		i__2 = *n;
#line 531 "zggev3.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 532 "zggev3.f"
		    i__3 = jr + jc * vr_dim1;
#line 532 "zggev3.f"
		    d__3 = temp, d__4 = (d__1 = vr[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&vr[jr + jc * vr_dim1]), abs(d__2));
#line 532 "zggev3.f"
		    temp = max(d__3,d__4);
#line 533 "zggev3.f"
/* L40: */
#line 533 "zggev3.f"
		}
#line 534 "zggev3.f"
		if (temp < smlnum) {
#line 534 "zggev3.f"
		    goto L60;
#line 534 "zggev3.f"
		}
#line 536 "zggev3.f"
		temp = 1. / temp;
#line 537 "zggev3.f"
		i__2 = *n;
#line 537 "zggev3.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 538 "zggev3.f"
		    i__3 = jr + jc * vr_dim1;
#line 538 "zggev3.f"
		    i__4 = jr + jc * vr_dim1;
#line 538 "zggev3.f"
		    z__1.r = temp * vr[i__4].r, z__1.i = temp * vr[i__4].i;
#line 538 "zggev3.f"
		    vr[i__3].r = z__1.r, vr[i__3].i = z__1.i;
#line 539 "zggev3.f"
/* L50: */
#line 539 "zggev3.f"
		}
#line 540 "zggev3.f"
L60:
#line 540 "zggev3.f"
		;
#line 540 "zggev3.f"
	    }
#line 541 "zggev3.f"
	}
#line 542 "zggev3.f"
    }

/*     Undo scaling if necessary */

#line 546 "zggev3.f"
L70:

#line 548 "zggev3.f"
    if (ilascl) {
#line 548 "zggev3.f"
	zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 548 "zggev3.f"
    }

#line 551 "zggev3.f"
    if (ilbscl) {
#line 551 "zggev3.f"
	zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 551 "zggev3.f"
    }

#line 554 "zggev3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 554 "zggev3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 555 "zggev3.f"
    return 0;

/*     End of ZGGEV3 */

} /* zggev3_ */

