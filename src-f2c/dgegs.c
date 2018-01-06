#line 1 "dgegs.f"
/* dgegs.f -- translated by f2c (version 20100827).
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

#line 1 "dgegs.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b36 = 0.;
static doublereal c_b37 = 1.;

/* > \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgegs.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgegs.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgegs.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR, */
/*                         ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, */
/*                         LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ), */
/*      $                   VSR( LDVSR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine DGGES. */
/* > */
/* > DGEGS computes the eigenvalues, real Schur form, and, optionally, */
/* > left and or/right Schur vectors of a real matrix pair (A,B). */
/* > Given two square matrices A and B, the generalized real Schur */
/* > factorization has the form */
/* > */
/* >   A = Q*S*Z**T,  B = Q*T*Z**T */
/* > */
/* > where Q and Z are orthogonal matrices, T is upper triangular, and S */
/* > is an upper quasi-triangular matrix with 1-by-1 and 2-by-2 diagonal */
/* > blocks, the 2-by-2 blocks corresponding to complex conjugate pairs */
/* > of eigenvalues of (A,B).  The columns of Q are the left Schur vectors */
/* > and the columns of Z are the right Schur vectors. */
/* > */
/* > If only the eigenvalues of (A,B) are needed, the driver routine */
/* > DGEGV should be used instead.  See DGEGV for a description of the */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the matrix A. */
/* >          On exit, the upper quasi-triangular matrix S from the */
/* >          generalized real Schur factorization. */
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
/* >          On entry, the matrix B. */
/* >          On exit, the upper triangular matrix T from the generalized */
/* >          real Schur factorization. */
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
/* >          The real parts of each scalar alpha defining an eigenvalue */
/* >          of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is DOUBLE PRECISION array, dimension (N) */
/* >          The imaginary parts of each scalar alpha defining an */
/* >          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th */
/* >          eigenvalue is real; if positive, then the j-th and (j+1)-st */
/* >          eigenvalues are a complex conjugate pair, with */
/* >          ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION array, dimension (N) */
/* >          The scalars beta that define the eigenvalues of GNEP. */
/* >          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/* >          beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* >          pair (A,B), in one of the forms lambda = alpha/beta or */
/* >          mu = beta/alpha.  Since either lambda or mu may overflow, */
/* >          they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* >          VSL is DOUBLE PRECISION array, dimension (LDVSL,N) */
/* >          If JOBVSL = 'V', the matrix of left Schur vectors Q. */
/* >          Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* >          LDVSL is INTEGER */
/* >          The leading dimension of the matrix VSL. LDVSL >=1, and */
/* >          if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* >          VSR is DOUBLE PRECISION array, dimension (LDVSR,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,4*N). */
/* >          For good performance, LWORK must generally be larger. */
/* >          To compute the optimal value of LWORK, call ILAENV to get */
/* >          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute: */
/* >          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR */
/* >          The optimal LWORK is  2*N + N*(NB+1). */
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
/* >                The QZ iteration failed.  (A,B) are not in Schur */
/* >                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should */
/* >                be correct for j=INFO+1,...,N. */
/* >          > N:  errors that usually indicate LAPACK problems: */
/* >                =N+1: error return from DGGBAL */
/* >                =N+2: error return from DGEQRF */
/* >                =N+3: error return from DORMQR */
/* >                =N+4: error return from DORGQR */
/* >                =N+5: error return from DGGHRD */
/* >                =N+6: error return from DHGEQZ (other than failed */
/* >                                                iteration) */
/* >                =N+7: error return from DGGBAK (computing VSL) */
/* >                =N+8: error return from DGGBAK (computing VSR) */
/* >                =N+9: error return from DLASCL (various places) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dgegs_(char *jobvsl, char *jobvsr, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, 
	integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work, 
	integer *lwork, integer *info, ftnlen jobvsl_len, ftnlen jobvsr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;

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
    extern /* Subroutine */ int dggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), dggbal_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
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
	    doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer ijobvl, iright, ijobvr;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal anrmto;
    static integer lwkmin;
    static doublereal bnrmto;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal smlnum;
    static integer lwkopt;
    static logical lquery;


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

#line 276 "dgegs.f"
    /* Parameter adjustments */
#line 276 "dgegs.f"
    a_dim1 = *lda;
#line 276 "dgegs.f"
    a_offset = 1 + a_dim1;
#line 276 "dgegs.f"
    a -= a_offset;
#line 276 "dgegs.f"
    b_dim1 = *ldb;
#line 276 "dgegs.f"
    b_offset = 1 + b_dim1;
#line 276 "dgegs.f"
    b -= b_offset;
#line 276 "dgegs.f"
    --alphar;
#line 276 "dgegs.f"
    --alphai;
#line 276 "dgegs.f"
    --beta;
#line 276 "dgegs.f"
    vsl_dim1 = *ldvsl;
#line 276 "dgegs.f"
    vsl_offset = 1 + vsl_dim1;
#line 276 "dgegs.f"
    vsl -= vsl_offset;
#line 276 "dgegs.f"
    vsr_dim1 = *ldvsr;
#line 276 "dgegs.f"
    vsr_offset = 1 + vsr_dim1;
#line 276 "dgegs.f"
    vsr -= vsr_offset;
#line 276 "dgegs.f"
    --work;
#line 276 "dgegs.f"

#line 276 "dgegs.f"
    /* Function Body */
#line 276 "dgegs.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 277 "dgegs.f"
	ijobvl = 1;
#line 278 "dgegs.f"
	ilvsl = FALSE_;
#line 279 "dgegs.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 280 "dgegs.f"
	ijobvl = 2;
#line 281 "dgegs.f"
	ilvsl = TRUE_;
#line 282 "dgegs.f"
    } else {
#line 283 "dgegs.f"
	ijobvl = -1;
#line 284 "dgegs.f"
	ilvsl = FALSE_;
#line 285 "dgegs.f"
    }

#line 287 "dgegs.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 288 "dgegs.f"
	ijobvr = 1;
#line 289 "dgegs.f"
	ilvsr = FALSE_;
#line 290 "dgegs.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 291 "dgegs.f"
	ijobvr = 2;
#line 292 "dgegs.f"
	ilvsr = TRUE_;
#line 293 "dgegs.f"
    } else {
#line 294 "dgegs.f"
	ijobvr = -1;
#line 295 "dgegs.f"
	ilvsr = FALSE_;
#line 296 "dgegs.f"
    }

/*     Test the input arguments */

/* Computing MAX */
#line 300 "dgegs.f"
    i__1 = *n << 2;
#line 300 "dgegs.f"
    lwkmin = max(i__1,1);
#line 301 "dgegs.f"
    lwkopt = lwkmin;
#line 302 "dgegs.f"
    work[1] = (doublereal) lwkopt;
#line 303 "dgegs.f"
    lquery = *lwork == -1;
#line 304 "dgegs.f"
    *info = 0;
#line 305 "dgegs.f"
    if (ijobvl <= 0) {
#line 306 "dgegs.f"
	*info = -1;
#line 307 "dgegs.f"
    } else if (ijobvr <= 0) {
#line 308 "dgegs.f"
	*info = -2;
#line 309 "dgegs.f"
    } else if (*n < 0) {
#line 310 "dgegs.f"
	*info = -3;
#line 311 "dgegs.f"
    } else if (*lda < max(1,*n)) {
#line 312 "dgegs.f"
	*info = -5;
#line 313 "dgegs.f"
    } else if (*ldb < max(1,*n)) {
#line 314 "dgegs.f"
	*info = -7;
#line 315 "dgegs.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 316 "dgegs.f"
	*info = -12;
#line 317 "dgegs.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 318 "dgegs.f"
	*info = -14;
#line 319 "dgegs.f"
    } else if (*lwork < lwkmin && ! lquery) {
#line 320 "dgegs.f"
	*info = -16;
#line 321 "dgegs.f"
    }

#line 323 "dgegs.f"
    if (*info == 0) {
#line 324 "dgegs.f"
	nb1 = ilaenv_(&c__1, "DGEQRF", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 325 "dgegs.f"
	nb2 = ilaenv_(&c__1, "DORMQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 326 "dgegs.f"
	nb3 = ilaenv_(&c__1, "DORGQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
/* Computing MAX */
#line 327 "dgegs.f"
	i__1 = max(nb1,nb2);
#line 327 "dgegs.f"
	nb = max(i__1,nb3);
#line 328 "dgegs.f"
	lopt = (*n << 1) + *n * (nb + 1);
#line 329 "dgegs.f"
	work[1] = (doublereal) lopt;
#line 330 "dgegs.f"
    }

#line 332 "dgegs.f"
    if (*info != 0) {
#line 333 "dgegs.f"
	i__1 = -(*info);
#line 333 "dgegs.f"
	xerbla_("DGEGS ", &i__1, (ftnlen)6);
#line 334 "dgegs.f"
	return 0;
#line 335 "dgegs.f"
    } else if (lquery) {
#line 336 "dgegs.f"
	return 0;
#line 337 "dgegs.f"
    }

/*     Quick return if possible */

#line 341 "dgegs.f"
    if (*n == 0) {
#line 341 "dgegs.f"
	return 0;
#line 341 "dgegs.f"
    }

/*     Get machine constants */

#line 346 "dgegs.f"
    eps = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 347 "dgegs.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 348 "dgegs.f"
    smlnum = *n * safmin / eps;
#line 349 "dgegs.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 353 "dgegs.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 354 "dgegs.f"
    ilascl = FALSE_;
#line 355 "dgegs.f"
    if (anrm > 0. && anrm < smlnum) {
#line 356 "dgegs.f"
	anrmto = smlnum;
#line 357 "dgegs.f"
	ilascl = TRUE_;
#line 358 "dgegs.f"
    } else if (anrm > bignum) {
#line 359 "dgegs.f"
	anrmto = bignum;
#line 360 "dgegs.f"
	ilascl = TRUE_;
#line 361 "dgegs.f"
    }

#line 363 "dgegs.f"
    if (ilascl) {
#line 364 "dgegs.f"
	dlascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 365 "dgegs.f"
	if (iinfo != 0) {
#line 366 "dgegs.f"
	    *info = *n + 9;
#line 367 "dgegs.f"
	    return 0;
#line 368 "dgegs.f"
	}
#line 369 "dgegs.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 373 "dgegs.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 374 "dgegs.f"
    ilbscl = FALSE_;
#line 375 "dgegs.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 376 "dgegs.f"
	bnrmto = smlnum;
#line 377 "dgegs.f"
	ilbscl = TRUE_;
#line 378 "dgegs.f"
    } else if (bnrm > bignum) {
#line 379 "dgegs.f"
	bnrmto = bignum;
#line 380 "dgegs.f"
	ilbscl = TRUE_;
#line 381 "dgegs.f"
    }

#line 383 "dgegs.f"
    if (ilbscl) {
#line 384 "dgegs.f"
	dlascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 385 "dgegs.f"
	if (iinfo != 0) {
#line 386 "dgegs.f"
	    *info = *n + 9;
#line 387 "dgegs.f"
	    return 0;
#line 388 "dgegs.f"
	}
#line 389 "dgegs.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     Workspace layout:  (2*N words -- "work..." not actually used) */
/*        left_permutation, right_permutation, work... */

#line 395 "dgegs.f"
    ileft = 1;
#line 396 "dgegs.f"
    iright = *n + 1;
#line 397 "dgegs.f"
    iwork = iright + *n;
#line 398 "dgegs.f"
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwork], &iinfo, (ftnlen)1);
#line 400 "dgegs.f"
    if (iinfo != 0) {
#line 401 "dgegs.f"
	*info = *n + 1;
#line 402 "dgegs.f"
	goto L10;
#line 403 "dgegs.f"
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR */
/*     Workspace layout:  ("work..." must have at least N words) */
/*        left_permutation, right_permutation, tau, work... */

#line 409 "dgegs.f"
    irows = ihi + 1 - ilo;
#line 410 "dgegs.f"
    icols = *n + 1 - ilo;
#line 411 "dgegs.f"
    itau = iwork;
#line 412 "dgegs.f"
    iwork = itau + irows;
#line 413 "dgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 413 "dgegs.f"
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
#line 415 "dgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 415 "dgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 415 "dgegs.f"
	lwkopt = max(i__1,i__2);
#line 415 "dgegs.f"
    }
#line 417 "dgegs.f"
    if (iinfo != 0) {
#line 418 "dgegs.f"
	*info = *n + 2;
#line 419 "dgegs.f"
	goto L10;
#line 420 "dgegs.f"
    }

#line 422 "dgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 422 "dgegs.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo, (ftnlen)1, (ftnlen)1);
#line 425 "dgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 425 "dgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 425 "dgegs.f"
	lwkopt = max(i__1,i__2);
#line 425 "dgegs.f"
    }
#line 427 "dgegs.f"
    if (iinfo != 0) {
#line 428 "dgegs.f"
	*info = *n + 3;
#line 429 "dgegs.f"
	goto L10;
#line 430 "dgegs.f"
    }

#line 432 "dgegs.f"
    if (ilvsl) {
#line 433 "dgegs.f"
	dlaset_("Full", n, n, &c_b36, &c_b37, &vsl[vsl_offset], ldvsl, (
		ftnlen)4);
#line 434 "dgegs.f"
	i__1 = irows - 1;
#line 434 "dgegs.f"
	i__2 = irows - 1;
#line 434 "dgegs.f"
	dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[ilo 
		+ 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 436 "dgegs.f"
	i__1 = *lwork + 1 - iwork;
#line 436 "dgegs.f"
	dorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwork], &i__1, &iinfo);
#line 439 "dgegs.f"
	if (iinfo >= 0) {
/* Computing MAX */
#line 439 "dgegs.f"
	    i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 439 "dgegs.f"
	    lwkopt = max(i__1,i__2);
#line 439 "dgegs.f"
	}
#line 441 "dgegs.f"
	if (iinfo != 0) {
#line 442 "dgegs.f"
	    *info = *n + 4;
#line 443 "dgegs.f"
	    goto L10;
#line 444 "dgegs.f"
	}
#line 445 "dgegs.f"
    }

#line 447 "dgegs.f"
    if (ilvsr) {
#line 447 "dgegs.f"
	dlaset_("Full", n, n, &c_b36, &c_b37, &vsr[vsr_offset], ldvsr, (
		ftnlen)4);
#line 447 "dgegs.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 452 "dgegs.f"
    dgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &iinfo, (
	    ftnlen)1, (ftnlen)1);
#line 454 "dgegs.f"
    if (iinfo != 0) {
#line 455 "dgegs.f"
	*info = *n + 5;
#line 456 "dgegs.f"
	goto L10;
#line 457 "dgegs.f"
    }

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     Workspace layout:  ("work..." must have at least 1 word) */
/*        left_permutation, right_permutation, work... */

#line 463 "dgegs.f"
    iwork = itau;
#line 464 "dgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 464 "dgegs.f"
    dhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
	    , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwork], &i__1, &iinfo, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 467 "dgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 467 "dgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 467 "dgegs.f"
	lwkopt = max(i__1,i__2);
#line 467 "dgegs.f"
    }
#line 469 "dgegs.f"
    if (iinfo != 0) {
#line 470 "dgegs.f"
	if (iinfo > 0 && iinfo <= *n) {
#line 471 "dgegs.f"
	    *info = iinfo;
#line 472 "dgegs.f"
	} else if (iinfo > *n && iinfo <= *n << 1) {
#line 473 "dgegs.f"
	    *info = iinfo - *n;
#line 474 "dgegs.f"
	} else {
#line 475 "dgegs.f"
	    *info = *n + 6;
#line 476 "dgegs.f"
	}
#line 477 "dgegs.f"
	goto L10;
#line 478 "dgegs.f"
    }

/*     Apply permutation to VSL and VSR */

#line 482 "dgegs.f"
    if (ilvsl) {
#line 483 "dgegs.f"
	dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &iinfo, (ftnlen)1, (ftnlen)1);
#line 485 "dgegs.f"
	if (iinfo != 0) {
#line 486 "dgegs.f"
	    *info = *n + 7;
#line 487 "dgegs.f"
	    goto L10;
#line 488 "dgegs.f"
	}
#line 489 "dgegs.f"
    }
#line 490 "dgegs.f"
    if (ilvsr) {
#line 491 "dgegs.f"
	dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 493 "dgegs.f"
	if (iinfo != 0) {
#line 494 "dgegs.f"
	    *info = *n + 8;
#line 495 "dgegs.f"
	    goto L10;
#line 496 "dgegs.f"
	}
#line 497 "dgegs.f"
    }

/*     Undo scaling */

#line 501 "dgegs.f"
    if (ilascl) {
#line 502 "dgegs.f"
	dlascl_("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 503 "dgegs.f"
	if (iinfo != 0) {
#line 504 "dgegs.f"
	    *info = *n + 9;
#line 505 "dgegs.f"
	    return 0;
#line 506 "dgegs.f"
	}
#line 507 "dgegs.f"
	dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		iinfo, (ftnlen)1);
#line 509 "dgegs.f"
	if (iinfo != 0) {
#line 510 "dgegs.f"
	    *info = *n + 9;
#line 511 "dgegs.f"
	    return 0;
#line 512 "dgegs.f"
	}
#line 513 "dgegs.f"
	dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		iinfo, (ftnlen)1);
#line 515 "dgegs.f"
	if (iinfo != 0) {
#line 516 "dgegs.f"
	    *info = *n + 9;
#line 517 "dgegs.f"
	    return 0;
#line 518 "dgegs.f"
	}
#line 519 "dgegs.f"
    }

#line 521 "dgegs.f"
    if (ilbscl) {
#line 522 "dgegs.f"
	dlascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 523 "dgegs.f"
	if (iinfo != 0) {
#line 524 "dgegs.f"
	    *info = *n + 9;
#line 525 "dgegs.f"
	    return 0;
#line 526 "dgegs.f"
	}
#line 527 "dgegs.f"
	dlascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		iinfo, (ftnlen)1);
#line 528 "dgegs.f"
	if (iinfo != 0) {
#line 529 "dgegs.f"
	    *info = *n + 9;
#line 530 "dgegs.f"
	    return 0;
#line 531 "dgegs.f"
	}
#line 532 "dgegs.f"
    }

#line 534 "dgegs.f"
L10:
#line 535 "dgegs.f"
    work[1] = (doublereal) lwkopt;

#line 537 "dgegs.f"
    return 0;

/*     End of DGEGS */

} /* dgegs_ */

