#line 1 "zggev.f"
/* zggev.f -- translated by f2c (version 20100827).
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

#line 1 "zggev.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> ZGGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, */
/*                         VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO ) */

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
/* > ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,2*N). */
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

/* > \date April 2012 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zggev_(char *jobvl, char *jobvr, integer *n, 
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
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), zggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical ldumma[1];
    static char chtemp[1];
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static integer ijobvl, iright;
    extern /* Subroutine */ int zgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static doublereal anrmto;
    static integer lwkmin;
    static doublereal bnrmto;
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the input arguments */

#line 282 "zggev.f"
    /* Parameter adjustments */
#line 282 "zggev.f"
    a_dim1 = *lda;
#line 282 "zggev.f"
    a_offset = 1 + a_dim1;
#line 282 "zggev.f"
    a -= a_offset;
#line 282 "zggev.f"
    b_dim1 = *ldb;
#line 282 "zggev.f"
    b_offset = 1 + b_dim1;
#line 282 "zggev.f"
    b -= b_offset;
#line 282 "zggev.f"
    --alpha;
#line 282 "zggev.f"
    --beta;
#line 282 "zggev.f"
    vl_dim1 = *ldvl;
#line 282 "zggev.f"
    vl_offset = 1 + vl_dim1;
#line 282 "zggev.f"
    vl -= vl_offset;
#line 282 "zggev.f"
    vr_dim1 = *ldvr;
#line 282 "zggev.f"
    vr_offset = 1 + vr_dim1;
#line 282 "zggev.f"
    vr -= vr_offset;
#line 282 "zggev.f"
    --work;
#line 282 "zggev.f"
    --rwork;
#line 282 "zggev.f"

#line 282 "zggev.f"
    /* Function Body */
#line 282 "zggev.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 283 "zggev.f"
	ijobvl = 1;
#line 284 "zggev.f"
	ilvl = FALSE_;
#line 285 "zggev.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 286 "zggev.f"
	ijobvl = 2;
#line 287 "zggev.f"
	ilvl = TRUE_;
#line 288 "zggev.f"
    } else {
#line 289 "zggev.f"
	ijobvl = -1;
#line 290 "zggev.f"
	ilvl = FALSE_;
#line 291 "zggev.f"
    }

#line 293 "zggev.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 294 "zggev.f"
	ijobvr = 1;
#line 295 "zggev.f"
	ilvr = FALSE_;
#line 296 "zggev.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 297 "zggev.f"
	ijobvr = 2;
#line 298 "zggev.f"
	ilvr = TRUE_;
#line 299 "zggev.f"
    } else {
#line 300 "zggev.f"
	ijobvr = -1;
#line 301 "zggev.f"
	ilvr = FALSE_;
#line 302 "zggev.f"
    }
#line 303 "zggev.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

#line 307 "zggev.f"
    *info = 0;
#line 308 "zggev.f"
    lquery = *lwork == -1;
#line 309 "zggev.f"
    if (ijobvl <= 0) {
#line 310 "zggev.f"
	*info = -1;
#line 311 "zggev.f"
    } else if (ijobvr <= 0) {
#line 312 "zggev.f"
	*info = -2;
#line 313 "zggev.f"
    } else if (*n < 0) {
#line 314 "zggev.f"
	*info = -3;
#line 315 "zggev.f"
    } else if (*lda < max(1,*n)) {
#line 316 "zggev.f"
	*info = -5;
#line 317 "zggev.f"
    } else if (*ldb < max(1,*n)) {
#line 318 "zggev.f"
	*info = -7;
#line 319 "zggev.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 320 "zggev.f"
	*info = -11;
#line 321 "zggev.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 322 "zggev.f"
	*info = -13;
#line 323 "zggev.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. The workspace is */
/*       computed assuming ILO = 1 and IHI = N, the worst case.) */

#line 333 "zggev.f"
    if (*info == 0) {
/* Computing MAX */
#line 334 "zggev.f"
	i__1 = 1, i__2 = *n << 1;
#line 334 "zggev.f"
	lwkmin = max(i__1,i__2);
/* Computing MAX */
#line 335 "zggev.f"
	i__1 = 1, i__2 = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", n, &c__1, n, 
		&c__0, (ftnlen)6, (ftnlen)1);
#line 335 "zggev.f"
	lwkopt = max(i__1,i__2);
/* Computing MAX */
#line 336 "zggev.f"
	i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "ZUNMQR", " ", n, &
		c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
#line 336 "zggev.f"
	lwkopt = max(i__1,i__2);
#line 338 "zggev.f"
	if (ilvl) {
/* Computing MAX */
#line 339 "zggev.f"
	    i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", " ", n, &
		    c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 339 "zggev.f"
	    lwkopt = max(i__1,i__2);
#line 341 "zggev.f"
	}
#line 342 "zggev.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 344 "zggev.f"
	if (*lwork < lwkmin && ! lquery) {
#line 344 "zggev.f"
	    *info = -15;
#line 344 "zggev.f"
	}
#line 346 "zggev.f"
    }

#line 348 "zggev.f"
    if (*info != 0) {
#line 349 "zggev.f"
	i__1 = -(*info);
#line 349 "zggev.f"
	xerbla_("ZGGEV ", &i__1, (ftnlen)6);
#line 350 "zggev.f"
	return 0;
#line 351 "zggev.f"
    } else if (lquery) {
#line 352 "zggev.f"
	return 0;
#line 353 "zggev.f"
    }

/*     Quick return if possible */

#line 357 "zggev.f"
    if (*n == 0) {
#line 357 "zggev.f"
	return 0;
#line 357 "zggev.f"
    }

/*     Get machine constants */

#line 362 "zggev.f"
    eps = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 363 "zggev.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 364 "zggev.f"
    bignum = 1. / smlnum;
#line 365 "zggev.f"
    dlabad_(&smlnum, &bignum);
#line 366 "zggev.f"
    smlnum = sqrt(smlnum) / eps;
#line 367 "zggev.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 371 "zggev.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 372 "zggev.f"
    ilascl = FALSE_;
#line 373 "zggev.f"
    if (anrm > 0. && anrm < smlnum) {
#line 374 "zggev.f"
	anrmto = smlnum;
#line 375 "zggev.f"
	ilascl = TRUE_;
#line 376 "zggev.f"
    } else if (anrm > bignum) {
#line 377 "zggev.f"
	anrmto = bignum;
#line 378 "zggev.f"
	ilascl = TRUE_;
#line 379 "zggev.f"
    }
#line 380 "zggev.f"
    if (ilascl) {
#line 380 "zggev.f"
	zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 380 "zggev.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 385 "zggev.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 386 "zggev.f"
    ilbscl = FALSE_;
#line 387 "zggev.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 388 "zggev.f"
	bnrmto = smlnum;
#line 389 "zggev.f"
	ilbscl = TRUE_;
#line 390 "zggev.f"
    } else if (bnrm > bignum) {
#line 391 "zggev.f"
	bnrmto = bignum;
#line 392 "zggev.f"
	ilbscl = TRUE_;
#line 393 "zggev.f"
    }
#line 394 "zggev.f"
    if (ilbscl) {
#line 394 "zggev.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr, (ftnlen)1);
#line 394 "zggev.f"
    }

/*     Permute the matrices A, B to isolate eigenvalues if possible */
/*     (Real Workspace: need 6*N) */

#line 400 "zggev.f"
    ileft = 1;
#line 401 "zggev.f"
    iright = *n + 1;
#line 402 "zggev.f"
    irwrk = iright + *n;
#line 403 "zggev.f"
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr, (ftnlen)1);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 409 "zggev.f"
    irows = ihi + 1 - ilo;
#line 410 "zggev.f"
    if (ilv) {
#line 411 "zggev.f"
	icols = *n + 1 - ilo;
#line 412 "zggev.f"
    } else {
#line 413 "zggev.f"
	icols = irows;
#line 414 "zggev.f"
    }
#line 415 "zggev.f"
    itau = 1;
#line 416 "zggev.f"
    iwrk = itau + irows;
#line 417 "zggev.f"
    i__1 = *lwork + 1 - iwrk;
#line 417 "zggev.f"
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 423 "zggev.f"
    i__1 = *lwork + 1 - iwrk;
#line 423 "zggev.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr, (ftnlen)1, (ftnlen)1);

/*     Initialize VL */
/*     (Complex Workspace: need N, prefer N*NB) */

#line 430 "zggev.f"
    if (ilvl) {
#line 431 "zggev.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vl[vl_offset], ldvl, (ftnlen)4);
#line 432 "zggev.f"
	if (irows > 1) {
#line 433 "zggev.f"
	    i__1 = irows - 1;
#line 433 "zggev.f"
	    i__2 = irows - 1;
#line 433 "zggev.f"
	    zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[
		    ilo + 1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 435 "zggev.f"
	}
#line 436 "zggev.f"
	i__1 = *lwork + 1 - iwrk;
#line 436 "zggev.f"
	zungqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwrk], &i__1, &ierr);
#line 438 "zggev.f"
    }

/*     Initialize VR */

#line 442 "zggev.f"
    if (ilvr) {
#line 442 "zggev.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vr[vr_offset], ldvr, (ftnlen)4);
#line 442 "zggev.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 447 "zggev.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 451 "zggev.f"
	zgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr, (
		ftnlen)1, (ftnlen)1);
#line 453 "zggev.f"
    } else {
#line 454 "zggev.f"
	zgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 456 "zggev.f"
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the */
/*     Schur form and Schur vectors) */
/*     (Complex Workspace: need N) */
/*     (Real Workspace: need N) */

#line 463 "zggev.f"
    iwrk = itau;
#line 464 "zggev.f"
    if (ilv) {
#line 465 "zggev.f"
	*(unsigned char *)chtemp = 'S';
#line 466 "zggev.f"
    } else {
#line 467 "zggev.f"
	*(unsigned char *)chtemp = 'E';
#line 468 "zggev.f"
    }
#line 469 "zggev.f"
    i__1 = *lwork + 1 - iwrk;
#line 469 "zggev.f"
    zhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vl[vl_offset], ldvl, &vr[
	    vr_offset], ldvr, &work[iwrk], &i__1, &rwork[irwrk], &ierr, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 472 "zggev.f"
    if (ierr != 0) {
#line 473 "zggev.f"
	if (ierr > 0 && ierr <= *n) {
#line 474 "zggev.f"
	    *info = ierr;
#line 475 "zggev.f"
	} else if (ierr > *n && ierr <= *n << 1) {
#line 476 "zggev.f"
	    *info = ierr - *n;
#line 477 "zggev.f"
	} else {
#line 478 "zggev.f"
	    *info = *n + 1;
#line 479 "zggev.f"
	}
#line 480 "zggev.f"
	goto L70;
#line 481 "zggev.f"
    }

/*     Compute Eigenvectors */
/*     (Real Workspace: need 2*N) */
/*     (Complex Workspace: need 2*N) */

#line 487 "zggev.f"
    if (ilv) {
#line 488 "zggev.f"
	if (ilvl) {
#line 489 "zggev.f"
	    if (ilvr) {
#line 490 "zggev.f"
		*(unsigned char *)chtemp = 'B';
#line 491 "zggev.f"
	    } else {
#line 492 "zggev.f"
		*(unsigned char *)chtemp = 'L';
#line 493 "zggev.f"
	    }
#line 494 "zggev.f"
	} else {
#line 495 "zggev.f"
	    *(unsigned char *)chtemp = 'R';
#line 496 "zggev.f"
	}

#line 498 "zggev.f"
	ztgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwrk], &rwork[irwrk], &ierr, (ftnlen)1, (ftnlen)1);
#line 501 "zggev.f"
	if (ierr != 0) {
#line 502 "zggev.f"
	    *info = *n + 2;
#line 503 "zggev.f"
	    goto L70;
#line 504 "zggev.f"
	}

/*        Undo balancing on VL and VR and normalization */
/*        (Workspace: none needed) */

#line 509 "zggev.f"
	if (ilvl) {
#line 510 "zggev.f"
	    zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
		     &vl[vl_offset], ldvl, &ierr, (ftnlen)1, (ftnlen)1);
#line 512 "zggev.f"
	    i__1 = *n;
#line 512 "zggev.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 513 "zggev.f"
		temp = 0.;
#line 514 "zggev.f"
		i__2 = *n;
#line 514 "zggev.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 515 "zggev.f"
		    i__3 = jr + jc * vl_dim1;
#line 515 "zggev.f"
		    d__3 = temp, d__4 = (d__1 = vl[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&vl[jr + jc * vl_dim1]), abs(d__2));
#line 515 "zggev.f"
		    temp = max(d__3,d__4);
#line 516 "zggev.f"
/* L10: */
#line 516 "zggev.f"
		}
#line 517 "zggev.f"
		if (temp < smlnum) {
#line 517 "zggev.f"
		    goto L30;
#line 517 "zggev.f"
		}
#line 519 "zggev.f"
		temp = 1. / temp;
#line 520 "zggev.f"
		i__2 = *n;
#line 520 "zggev.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 521 "zggev.f"
		    i__3 = jr + jc * vl_dim1;
#line 521 "zggev.f"
		    i__4 = jr + jc * vl_dim1;
#line 521 "zggev.f"
		    z__1.r = temp * vl[i__4].r, z__1.i = temp * vl[i__4].i;
#line 521 "zggev.f"
		    vl[i__3].r = z__1.r, vl[i__3].i = z__1.i;
#line 522 "zggev.f"
/* L20: */
#line 522 "zggev.f"
		}
#line 523 "zggev.f"
L30:
#line 523 "zggev.f"
		;
#line 523 "zggev.f"
	    }
#line 524 "zggev.f"
	}
#line 525 "zggev.f"
	if (ilvr) {
#line 526 "zggev.f"
	    zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
		     &vr[vr_offset], ldvr, &ierr, (ftnlen)1, (ftnlen)1);
#line 528 "zggev.f"
	    i__1 = *n;
#line 528 "zggev.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 529 "zggev.f"
		temp = 0.;
#line 530 "zggev.f"
		i__2 = *n;
#line 530 "zggev.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 531 "zggev.f"
		    i__3 = jr + jc * vr_dim1;
#line 531 "zggev.f"
		    d__3 = temp, d__4 = (d__1 = vr[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&vr[jr + jc * vr_dim1]), abs(d__2));
#line 531 "zggev.f"
		    temp = max(d__3,d__4);
#line 532 "zggev.f"
/* L40: */
#line 532 "zggev.f"
		}
#line 533 "zggev.f"
		if (temp < smlnum) {
#line 533 "zggev.f"
		    goto L60;
#line 533 "zggev.f"
		}
#line 535 "zggev.f"
		temp = 1. / temp;
#line 536 "zggev.f"
		i__2 = *n;
#line 536 "zggev.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 537 "zggev.f"
		    i__3 = jr + jc * vr_dim1;
#line 537 "zggev.f"
		    i__4 = jr + jc * vr_dim1;
#line 537 "zggev.f"
		    z__1.r = temp * vr[i__4].r, z__1.i = temp * vr[i__4].i;
#line 537 "zggev.f"
		    vr[i__3].r = z__1.r, vr[i__3].i = z__1.i;
#line 538 "zggev.f"
/* L50: */
#line 538 "zggev.f"
		}
#line 539 "zggev.f"
L60:
#line 539 "zggev.f"
		;
#line 539 "zggev.f"
	    }
#line 540 "zggev.f"
	}
#line 541 "zggev.f"
    }

/*     Undo scaling if necessary */

#line 545 "zggev.f"
L70:

#line 547 "zggev.f"
    if (ilascl) {
#line 547 "zggev.f"
	zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr, (ftnlen)1);
#line 547 "zggev.f"
    }

#line 550 "zggev.f"
    if (ilbscl) {
#line 550 "zggev.f"
	zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr, (ftnlen)1);
#line 550 "zggev.f"
    }

#line 553 "zggev.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 554 "zggev.f"
    return 0;

/*     End of ZGGEV */

} /* zggev_ */

