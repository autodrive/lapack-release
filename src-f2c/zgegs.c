#line 1 "zgegs.f"
/* zgegs.f -- translated by f2c (version 20100827).
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

#line 1 "zgegs.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> ZGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgegs.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgegs.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgegs.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA, */
/*                         VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine ZGGES. */
/* > */
/* > ZGEGS computes the eigenvalues, Schur form, and, optionally, the */
/* > left and or/right Schur vectors of a complex matrix pair (A,B). */
/* > Given two square matrices A and B, the generalized Schur */
/* > factorization has the form */
/* > */
/* >    A = Q*S*Z**H,  B = Q*T*Z**H */
/* > */
/* > where Q and Z are unitary matrices and S and T are upper triangular. */
/* > The columns of Q are the left Schur vectors */
/* > and the columns of Z are the right Schur vectors. */
/* > */
/* > If only the eigenvalues of (A,B) are needed, the driver routine */
/* > ZGEGV should be used instead.  See ZGEGV for a description of the */
/* > eigenvalues of the generalized nonsymmetric eigenvalue problem */
/* > (GNEP). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBVSL */
/* > \verbatim */
/* >          JOBVSL is CHARACTER*1 */
/* >          = 'N':  do not compute the left Schur vectors; */
/* >          = 'V':  compute the left Schur vectors (returned in VSL). */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVSR */
/* > \verbatim */
/* >          JOBVSR is CHARACTER*1 */
/* >          = 'N':  do not compute the right Schur vectors; */
/* >          = 'V':  compute the right Schur vectors (returned in VSR). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A, B, VSL, and VSR.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the matrix A. */
/* >          On exit, the upper triangular matrix S from the generalized */
/* >          Schur factorization. */
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
/* >          On entry, the matrix B. */
/* >          On exit, the upper triangular matrix T from the generalized */
/* >          Schur factorization. */
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
/* >          The complex scalars alpha that define the eigenvalues of */
/* >          GNEP.  ALPHA(j) = S(j,j), the diagonal element of the Schur */
/* >          form of A. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 array, dimension (N) */
/* >          The non-negative real scalars beta that define the */
/* >          eigenvalues of GNEP.  BETA(j) = T(j,j), the diagonal element */
/* >          of the triangular factor T. */
/* > */
/* >          Together, the quantities alpha = ALPHA(j) and beta = BETA(j) */
/* >          represent the j-th eigenvalue of the matrix pair (A,B), in */
/* >          one of the forms lambda = alpha/beta or mu = beta/alpha. */
/* >          Since either lambda or mu may overflow, they should not, */
/* >          in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* >          VSL is COMPLEX*16 array, dimension (LDVSL,N) */
/* >          If JOBVSL = 'V', the matrix of left Schur vectors Q. */
/* >          Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* >          LDVSL is INTEGER */
/* >          The leading dimension of the matrix VSL. LDVSL >= 1, and */
/* >          if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* >          VSR is COMPLEX*16 array, dimension (LDVSR,N) */
/* >          If JOBVSR = 'V', the matrix of right Schur vectors Z. */
/* >          Not referenced if JOBVSR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSR */
/* > \verbatim */
/* >          LDVSR is INTEGER */
/* >          The leading dimension of the matrix VSR. LDVSR >= 1, and */
/* >          if JOBVSR = 'V', LDVSR >= N. */
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
/* >          To compute the optimal value of LWORK, call ILAENV to get */
/* >          blocksizes (for ZGEQRF, ZUNMQR, and CUNGQR.)  Then compute: */
/* >          NB  -- MAX of the blocksizes for ZGEQRF, ZUNMQR, and CUNGQR; */
/* >          the optimal LWORK is N*(NB+1). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          =1,...,N: */
/* >                The QZ iteration failed.  (A,B) are not in Schur */
/* >                form, but ALPHA(j) and BETA(j) should be correct for */
/* >                j=INFO+1,...,N. */
/* >          > N:  errors that usually indicate LAPACK problems: */
/* >                =N+1: error return from ZGGBAL */
/* >                =N+2: error return from ZGEQRF */
/* >                =N+3: error return from ZUNMQR */
/* >                =N+4: error return from ZUNGQR */
/* >                =N+5: error return from ZGGHRD */
/* >                =N+6: error return from ZHGEQZ (other than failed */
/* >                                               iteration) */
/* >                =N+7: error return from ZGGBAK (computing VSL) */
/* >                =N+8: error return from ZGGBAK (computing VSR) */
/* >                =N+9: error return from ZLASCL (various places) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zgegs_(char *jobvsl, char *jobvsr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, 
	integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *info, ftnlen 
	jobvsl_len, ftnlen jobvsr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer nb, nb1, nb2, nb3, ihi, ilo;
    static doublereal eps, anrm, bnrm;
    static integer itau, lopt;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ileft, iinfo, icols;
    static logical ilvsl;
    static integer iwork;
    static logical ilvsr;
    static integer irows;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), zggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ilascl, ilbscl;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
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
	    doublecomplex *, doublecomplex *, integer *, ftnlen), zhgeqz_(
	    char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal smlnum;
    static integer irwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the input arguments */

#line 278 "zgegs.f"
    /* Parameter adjustments */
#line 278 "zgegs.f"
    a_dim1 = *lda;
#line 278 "zgegs.f"
    a_offset = 1 + a_dim1;
#line 278 "zgegs.f"
    a -= a_offset;
#line 278 "zgegs.f"
    b_dim1 = *ldb;
#line 278 "zgegs.f"
    b_offset = 1 + b_dim1;
#line 278 "zgegs.f"
    b -= b_offset;
#line 278 "zgegs.f"
    --alpha;
#line 278 "zgegs.f"
    --beta;
#line 278 "zgegs.f"
    vsl_dim1 = *ldvsl;
#line 278 "zgegs.f"
    vsl_offset = 1 + vsl_dim1;
#line 278 "zgegs.f"
    vsl -= vsl_offset;
#line 278 "zgegs.f"
    vsr_dim1 = *ldvsr;
#line 278 "zgegs.f"
    vsr_offset = 1 + vsr_dim1;
#line 278 "zgegs.f"
    vsr -= vsr_offset;
#line 278 "zgegs.f"
    --work;
#line 278 "zgegs.f"
    --rwork;
#line 278 "zgegs.f"

#line 278 "zgegs.f"
    /* Function Body */
#line 278 "zgegs.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 279 "zgegs.f"
	ijobvl = 1;
#line 280 "zgegs.f"
	ilvsl = FALSE_;
#line 281 "zgegs.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 282 "zgegs.f"
	ijobvl = 2;
#line 283 "zgegs.f"
	ilvsl = TRUE_;
#line 284 "zgegs.f"
    } else {
#line 285 "zgegs.f"
	ijobvl = -1;
#line 286 "zgegs.f"
	ilvsl = FALSE_;
#line 287 "zgegs.f"
    }

#line 289 "zgegs.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 290 "zgegs.f"
	ijobvr = 1;
#line 291 "zgegs.f"
	ilvsr = FALSE_;
#line 292 "zgegs.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 293 "zgegs.f"
	ijobvr = 2;
#line 294 "zgegs.f"
	ilvsr = TRUE_;
#line 295 "zgegs.f"
    } else {
#line 296 "zgegs.f"
	ijobvr = -1;
#line 297 "zgegs.f"
	ilvsr = FALSE_;
#line 298 "zgegs.f"
    }

/*     Test the input arguments */

/* Computing MAX */
#line 302 "zgegs.f"
    i__1 = *n << 1;
#line 302 "zgegs.f"
    lwkmin = max(i__1,1);
#line 303 "zgegs.f"
    lwkopt = lwkmin;
#line 304 "zgegs.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 305 "zgegs.f"
    lquery = *lwork == -1;
#line 306 "zgegs.f"
    *info = 0;
#line 307 "zgegs.f"
    if (ijobvl <= 0) {
#line 308 "zgegs.f"
	*info = -1;
#line 309 "zgegs.f"
    } else if (ijobvr <= 0) {
#line 310 "zgegs.f"
	*info = -2;
#line 311 "zgegs.f"
    } else if (*n < 0) {
#line 312 "zgegs.f"
	*info = -3;
#line 313 "zgegs.f"
    } else if (*lda < max(1,*n)) {
#line 314 "zgegs.f"
	*info = -5;
#line 315 "zgegs.f"
    } else if (*ldb < max(1,*n)) {
#line 316 "zgegs.f"
	*info = -7;
#line 317 "zgegs.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 318 "zgegs.f"
	*info = -11;
#line 319 "zgegs.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 320 "zgegs.f"
	*info = -13;
#line 321 "zgegs.f"
    } else if (*lwork < lwkmin && ! lquery) {
#line 322 "zgegs.f"
	*info = -15;
#line 323 "zgegs.f"
    }

#line 325 "zgegs.f"
    if (*info == 0) {
#line 326 "zgegs.f"
	nb1 = ilaenv_(&c__1, "ZGEQRF", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 327 "zgegs.f"
	nb2 = ilaenv_(&c__1, "ZUNMQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 328 "zgegs.f"
	nb3 = ilaenv_(&c__1, "ZUNGQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
/* Computing MAX */
#line 329 "zgegs.f"
	i__1 = max(nb1,nb2);
#line 329 "zgegs.f"
	nb = max(i__1,nb3);
#line 330 "zgegs.f"
	lopt = *n * (nb + 1);
#line 331 "zgegs.f"
	work[1].r = (doublereal) lopt, work[1].i = 0.;
#line 332 "zgegs.f"
    }

#line 334 "zgegs.f"
    if (*info != 0) {
#line 335 "zgegs.f"
	i__1 = -(*info);
#line 335 "zgegs.f"
	xerbla_("ZGEGS ", &i__1, (ftnlen)6);
#line 336 "zgegs.f"
	return 0;
#line 337 "zgegs.f"
    } else if (lquery) {
#line 338 "zgegs.f"
	return 0;
#line 339 "zgegs.f"
    }

/*     Quick return if possible */

#line 343 "zgegs.f"
    if (*n == 0) {
#line 343 "zgegs.f"
	return 0;
#line 343 "zgegs.f"
    }

/*     Get machine constants */

#line 348 "zgegs.f"
    eps = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 349 "zgegs.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 350 "zgegs.f"
    smlnum = *n * safmin / eps;
#line 351 "zgegs.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 355 "zgegs.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 356 "zgegs.f"
    ilascl = FALSE_;
#line 357 "zgegs.f"
    if (anrm > 0. && anrm < smlnum) {
#line 358 "zgegs.f"
	anrmto = smlnum;
#line 359 "zgegs.f"
	ilascl = TRUE_;
#line 360 "zgegs.f"
    } else if (anrm > bignum) {
#line 361 "zgegs.f"
	anrmto = bignum;
#line 362 "zgegs.f"
	ilascl = TRUE_;
#line 363 "zgegs.f"
    }

#line 365 "zgegs.f"
    if (ilascl) {
#line 366 "zgegs.f"
	zlascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 367 "zgegs.f"
	if (iinfo != 0) {
#line 368 "zgegs.f"
	    *info = *n + 9;
#line 369 "zgegs.f"
	    return 0;
#line 370 "zgegs.f"
	}
#line 371 "zgegs.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 375 "zgegs.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 376 "zgegs.f"
    ilbscl = FALSE_;
#line 377 "zgegs.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 378 "zgegs.f"
	bnrmto = smlnum;
#line 379 "zgegs.f"
	ilbscl = TRUE_;
#line 380 "zgegs.f"
    } else if (bnrm > bignum) {
#line 381 "zgegs.f"
	bnrmto = bignum;
#line 382 "zgegs.f"
	ilbscl = TRUE_;
#line 383 "zgegs.f"
    }

#line 385 "zgegs.f"
    if (ilbscl) {
#line 386 "zgegs.f"
	zlascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 387 "zgegs.f"
	if (iinfo != 0) {
#line 388 "zgegs.f"
	    *info = *n + 9;
#line 389 "zgegs.f"
	    return 0;
#line 390 "zgegs.f"
	}
#line 391 "zgegs.f"
    }

/*     Permute the matrix to make it more nearly triangular */

#line 395 "zgegs.f"
    ileft = 1;
#line 396 "zgegs.f"
    iright = *n + 1;
#line 397 "zgegs.f"
    irwork = iright + *n;
#line 398 "zgegs.f"
    iwork = 1;
#line 399 "zgegs.f"
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwork], &iinfo, (ftnlen)1);
#line 401 "zgegs.f"
    if (iinfo != 0) {
#line 402 "zgegs.f"
	*info = *n + 1;
#line 403 "zgegs.f"
	goto L10;
#line 404 "zgegs.f"
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR */

#line 408 "zgegs.f"
    irows = ihi + 1 - ilo;
#line 409 "zgegs.f"
    icols = *n + 1 - ilo;
#line 410 "zgegs.f"
    itau = iwork;
#line 411 "zgegs.f"
    iwork = itau + irows;
#line 412 "zgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 412 "zgegs.f"
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
#line 414 "zgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 414 "zgegs.f"
	i__3 = iwork;
#line 414 "zgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 414 "zgegs.f"
	lwkopt = max(i__1,i__2);
#line 414 "zgegs.f"
    }
#line 416 "zgegs.f"
    if (iinfo != 0) {
#line 417 "zgegs.f"
	*info = *n + 2;
#line 418 "zgegs.f"
	goto L10;
#line 419 "zgegs.f"
    }

#line 421 "zgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 421 "zgegs.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo, (ftnlen)1, (ftnlen)1);
#line 424 "zgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 424 "zgegs.f"
	i__3 = iwork;
#line 424 "zgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 424 "zgegs.f"
	lwkopt = max(i__1,i__2);
#line 424 "zgegs.f"
    }
#line 426 "zgegs.f"
    if (iinfo != 0) {
#line 427 "zgegs.f"
	*info = *n + 3;
#line 428 "zgegs.f"
	goto L10;
#line 429 "zgegs.f"
    }

#line 431 "zgegs.f"
    if (ilvsl) {
#line 432 "zgegs.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 433 "zgegs.f"
	i__1 = irows - 1;
#line 433 "zgegs.f"
	i__2 = irows - 1;
#line 433 "zgegs.f"
	zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[ilo 
		+ 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 435 "zgegs.f"
	i__1 = *lwork + 1 - iwork;
#line 435 "zgegs.f"
	zungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwork], &i__1, &iinfo);
#line 438 "zgegs.f"
	if (iinfo >= 0) {
/* Computing MAX */
#line 438 "zgegs.f"
	    i__3 = iwork;
#line 438 "zgegs.f"
	    i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 438 "zgegs.f"
	    lwkopt = max(i__1,i__2);
#line 438 "zgegs.f"
	}
#line 440 "zgegs.f"
	if (iinfo != 0) {
#line 441 "zgegs.f"
	    *info = *n + 4;
#line 442 "zgegs.f"
	    goto L10;
#line 443 "zgegs.f"
	}
#line 444 "zgegs.f"
    }

#line 446 "zgegs.f"
    if (ilvsr) {
#line 446 "zgegs.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 446 "zgegs.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 451 "zgegs.f"
    zgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &iinfo, (
	    ftnlen)1, (ftnlen)1);
#line 453 "zgegs.f"
    if (iinfo != 0) {
#line 454 "zgegs.f"
	*info = *n + 5;
#line 455 "zgegs.f"
	goto L10;
#line 456 "zgegs.f"
    }

/*     Perform QZ algorithm, computing Schur vectors if desired */

#line 460 "zgegs.f"
    iwork = itau;
#line 461 "zgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 461 "zgegs.f"
    zhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwork], &i__1, &rwork[irwork], &
	    iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 464 "zgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 464 "zgegs.f"
	i__3 = iwork;
#line 464 "zgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 464 "zgegs.f"
	lwkopt = max(i__1,i__2);
#line 464 "zgegs.f"
    }
#line 466 "zgegs.f"
    if (iinfo != 0) {
#line 467 "zgegs.f"
	if (iinfo > 0 && iinfo <= *n) {
#line 468 "zgegs.f"
	    *info = iinfo;
#line 469 "zgegs.f"
	} else if (iinfo > *n && iinfo <= *n << 1) {
#line 470 "zgegs.f"
	    *info = iinfo - *n;
#line 471 "zgegs.f"
	} else {
#line 472 "zgegs.f"
	    *info = *n + 6;
#line 473 "zgegs.f"
	}
#line 474 "zgegs.f"
	goto L10;
#line 475 "zgegs.f"
    }

/*     Apply permutation to VSL and VSR */

#line 479 "zgegs.f"
    if (ilvsl) {
#line 480 "zgegs.f"
	zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &iinfo, (ftnlen)1, (ftnlen)1);
#line 482 "zgegs.f"
	if (iinfo != 0) {
#line 483 "zgegs.f"
	    *info = *n + 7;
#line 484 "zgegs.f"
	    goto L10;
#line 485 "zgegs.f"
	}
#line 486 "zgegs.f"
    }
#line 487 "zgegs.f"
    if (ilvsr) {
#line 488 "zgegs.f"
	zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 490 "zgegs.f"
	if (iinfo != 0) {
#line 491 "zgegs.f"
	    *info = *n + 8;
#line 492 "zgegs.f"
	    goto L10;
#line 493 "zgegs.f"
	}
#line 494 "zgegs.f"
    }

/*     Undo scaling */

#line 498 "zgegs.f"
    if (ilascl) {
#line 499 "zgegs.f"
	zlascl_("U", &c_n1, &c_n1, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 500 "zgegs.f"
	if (iinfo != 0) {
#line 501 "zgegs.f"
	    *info = *n + 9;
#line 502 "zgegs.f"
	    return 0;
#line 503 "zgegs.f"
	}
#line 504 "zgegs.f"
	zlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		iinfo, (ftnlen)1);
#line 505 "zgegs.f"
	if (iinfo != 0) {
#line 506 "zgegs.f"
	    *info = *n + 9;
#line 507 "zgegs.f"
	    return 0;
#line 508 "zgegs.f"
	}
#line 509 "zgegs.f"
    }

#line 511 "zgegs.f"
    if (ilbscl) {
#line 512 "zgegs.f"
	zlascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 513 "zgegs.f"
	if (iinfo != 0) {
#line 514 "zgegs.f"
	    *info = *n + 9;
#line 515 "zgegs.f"
	    return 0;
#line 516 "zgegs.f"
	}
#line 517 "zgegs.f"
	zlascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		iinfo, (ftnlen)1);
#line 518 "zgegs.f"
	if (iinfo != 0) {
#line 519 "zgegs.f"
	    *info = *n + 9;
#line 520 "zgegs.f"
	    return 0;
#line 521 "zgegs.f"
	}
#line 522 "zgegs.f"
    }

#line 524 "zgegs.f"
L10:
#line 525 "zgegs.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 527 "zgegs.f"
    return 0;

/*     End of ZGEGS */

} /* zgegs_ */

