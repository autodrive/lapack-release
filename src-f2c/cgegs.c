#line 1 "cgegs.f"
/* cgegs.f -- translated by f2c (version 20100827).
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

#line 1 "cgegs.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> CGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgegs.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgegs.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgegs.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA, */
/*                         VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVSL, JOBVSR */
/*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ) */
/*       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine CGGES. */
/* > */
/* > CGEGS computes the eigenvalues, Schur form, and, optionally, the */
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
/* > CGEGV should be used instead.  See CGEGV for a description of the */
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
/* >          A is COMPLEX array, dimension (LDA, N) */
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
/* >          B is COMPLEX array, dimension (LDB, N) */
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
/* >          ALPHA is COMPLEX array, dimension (N) */
/* >          The complex scalars alpha that define the eigenvalues of */
/* >          GNEP.  ALPHA(j) = S(j,j), the diagonal element of the Schur */
/* >          form of A. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX array, dimension (N) */
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
/* >          VSL is COMPLEX array, dimension (LDVSL,N) */
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
/* >          VSR is COMPLEX array, dimension (LDVSR,N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,2*N). */
/* >          For good performance, LWORK must generally be larger. */
/* >          To compute the optimal value of LWORK, call ILAENV to get */
/* >          blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute: */
/* >          NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR; */
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
/* >          RWORK is REAL array, dimension (3*N) */
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
/* >                =N+1: error return from CGGBAL */
/* >                =N+2: error return from CGEQRF */
/* >                =N+3: error return from CUNMQR */
/* >                =N+4: error return from CUNGQR */
/* >                =N+5: error return from CGGHRD */
/* >                =N+6: error return from CHGEQZ (other than failed */
/* >                                               iteration) */
/* >                =N+7: error return from CGGBAK (computing VSL) */
/* >                =N+8: error return from CGGBAK (computing VSR) */
/* >                =N+9: error return from CLASCL (various places) */
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
/* Subroutine */ int cgegs_(char *jobvsl, char *jobvsr, integer *n, 
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
    extern /* Subroutine */ int cggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), cggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static logical ilascl, ilbscl;
    extern /* Subroutine */ int cgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int chgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static integer ijobvl, iright, ijobvr;
    static doublereal anrmto;
    static integer lwkmin;
    static doublereal bnrmto;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), cunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static doublereal smlnum;
    static integer irwork, lwkopt;
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

#line 278 "cgegs.f"
    /* Parameter adjustments */
#line 278 "cgegs.f"
    a_dim1 = *lda;
#line 278 "cgegs.f"
    a_offset = 1 + a_dim1;
#line 278 "cgegs.f"
    a -= a_offset;
#line 278 "cgegs.f"
    b_dim1 = *ldb;
#line 278 "cgegs.f"
    b_offset = 1 + b_dim1;
#line 278 "cgegs.f"
    b -= b_offset;
#line 278 "cgegs.f"
    --alpha;
#line 278 "cgegs.f"
    --beta;
#line 278 "cgegs.f"
    vsl_dim1 = *ldvsl;
#line 278 "cgegs.f"
    vsl_offset = 1 + vsl_dim1;
#line 278 "cgegs.f"
    vsl -= vsl_offset;
#line 278 "cgegs.f"
    vsr_dim1 = *ldvsr;
#line 278 "cgegs.f"
    vsr_offset = 1 + vsr_dim1;
#line 278 "cgegs.f"
    vsr -= vsr_offset;
#line 278 "cgegs.f"
    --work;
#line 278 "cgegs.f"
    --rwork;
#line 278 "cgegs.f"

#line 278 "cgegs.f"
    /* Function Body */
#line 278 "cgegs.f"
    if (lsame_(jobvsl, "N", (ftnlen)1, (ftnlen)1)) {
#line 279 "cgegs.f"
	ijobvl = 1;
#line 280 "cgegs.f"
	ilvsl = FALSE_;
#line 281 "cgegs.f"
    } else if (lsame_(jobvsl, "V", (ftnlen)1, (ftnlen)1)) {
#line 282 "cgegs.f"
	ijobvl = 2;
#line 283 "cgegs.f"
	ilvsl = TRUE_;
#line 284 "cgegs.f"
    } else {
#line 285 "cgegs.f"
	ijobvl = -1;
#line 286 "cgegs.f"
	ilvsl = FALSE_;
#line 287 "cgegs.f"
    }

#line 289 "cgegs.f"
    if (lsame_(jobvsr, "N", (ftnlen)1, (ftnlen)1)) {
#line 290 "cgegs.f"
	ijobvr = 1;
#line 291 "cgegs.f"
	ilvsr = FALSE_;
#line 292 "cgegs.f"
    } else if (lsame_(jobvsr, "V", (ftnlen)1, (ftnlen)1)) {
#line 293 "cgegs.f"
	ijobvr = 2;
#line 294 "cgegs.f"
	ilvsr = TRUE_;
#line 295 "cgegs.f"
    } else {
#line 296 "cgegs.f"
	ijobvr = -1;
#line 297 "cgegs.f"
	ilvsr = FALSE_;
#line 298 "cgegs.f"
    }

/*     Test the input arguments */

/* Computing MAX */
#line 302 "cgegs.f"
    i__1 = *n << 1;
#line 302 "cgegs.f"
    lwkmin = max(i__1,1);
#line 303 "cgegs.f"
    lwkopt = lwkmin;
#line 304 "cgegs.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 305 "cgegs.f"
    lquery = *lwork == -1;
#line 306 "cgegs.f"
    *info = 0;
#line 307 "cgegs.f"
    if (ijobvl <= 0) {
#line 308 "cgegs.f"
	*info = -1;
#line 309 "cgegs.f"
    } else if (ijobvr <= 0) {
#line 310 "cgegs.f"
	*info = -2;
#line 311 "cgegs.f"
    } else if (*n < 0) {
#line 312 "cgegs.f"
	*info = -3;
#line 313 "cgegs.f"
    } else if (*lda < max(1,*n)) {
#line 314 "cgegs.f"
	*info = -5;
#line 315 "cgegs.f"
    } else if (*ldb < max(1,*n)) {
#line 316 "cgegs.f"
	*info = -7;
#line 317 "cgegs.f"
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
#line 318 "cgegs.f"
	*info = -11;
#line 319 "cgegs.f"
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
#line 320 "cgegs.f"
	*info = -13;
#line 321 "cgegs.f"
    } else if (*lwork < lwkmin && ! lquery) {
#line 322 "cgegs.f"
	*info = -15;
#line 323 "cgegs.f"
    }

#line 325 "cgegs.f"
    if (*info == 0) {
#line 326 "cgegs.f"
	nb1 = ilaenv_(&c__1, "CGEQRF", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 327 "cgegs.f"
	nb2 = ilaenv_(&c__1, "CUNMQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 328 "cgegs.f"
	nb3 = ilaenv_(&c__1, "CUNGQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
/* Computing MAX */
#line 329 "cgegs.f"
	i__1 = max(nb1,nb2);
#line 329 "cgegs.f"
	nb = max(i__1,nb3);
#line 330 "cgegs.f"
	lopt = *n * (nb + 1);
#line 331 "cgegs.f"
	work[1].r = (doublereal) lopt, work[1].i = 0.;
#line 332 "cgegs.f"
    }

#line 334 "cgegs.f"
    if (*info != 0) {
#line 335 "cgegs.f"
	i__1 = -(*info);
#line 335 "cgegs.f"
	xerbla_("CGEGS ", &i__1, (ftnlen)6);
#line 336 "cgegs.f"
	return 0;
#line 337 "cgegs.f"
    } else if (lquery) {
#line 338 "cgegs.f"
	return 0;
#line 339 "cgegs.f"
    }

/*     Quick return if possible */

#line 343 "cgegs.f"
    if (*n == 0) {
#line 343 "cgegs.f"
	return 0;
#line 343 "cgegs.f"
    }

/*     Get machine constants */

#line 348 "cgegs.f"
    eps = slamch_("E", (ftnlen)1) * slamch_("B", (ftnlen)1);
#line 349 "cgegs.f"
    safmin = slamch_("S", (ftnlen)1);
#line 350 "cgegs.f"
    smlnum = *n * safmin / eps;
#line 351 "cgegs.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 355 "cgegs.f"
    anrm = clange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 356 "cgegs.f"
    ilascl = FALSE_;
#line 357 "cgegs.f"
    if (anrm > 0. && anrm < smlnum) {
#line 358 "cgegs.f"
	anrmto = smlnum;
#line 359 "cgegs.f"
	ilascl = TRUE_;
#line 360 "cgegs.f"
    } else if (anrm > bignum) {
#line 361 "cgegs.f"
	anrmto = bignum;
#line 362 "cgegs.f"
	ilascl = TRUE_;
#line 363 "cgegs.f"
    }

#line 365 "cgegs.f"
    if (ilascl) {
#line 366 "cgegs.f"
	clascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 367 "cgegs.f"
	if (iinfo != 0) {
#line 368 "cgegs.f"
	    *info = *n + 9;
#line 369 "cgegs.f"
	    return 0;
#line 370 "cgegs.f"
	}
#line 371 "cgegs.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 375 "cgegs.f"
    bnrm = clange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 376 "cgegs.f"
    ilbscl = FALSE_;
#line 377 "cgegs.f"
    if (bnrm > 0. && bnrm < smlnum) {
#line 378 "cgegs.f"
	bnrmto = smlnum;
#line 379 "cgegs.f"
	ilbscl = TRUE_;
#line 380 "cgegs.f"
    } else if (bnrm > bignum) {
#line 381 "cgegs.f"
	bnrmto = bignum;
#line 382 "cgegs.f"
	ilbscl = TRUE_;
#line 383 "cgegs.f"
    }

#line 385 "cgegs.f"
    if (ilbscl) {
#line 386 "cgegs.f"
	clascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 387 "cgegs.f"
	if (iinfo != 0) {
#line 388 "cgegs.f"
	    *info = *n + 9;
#line 389 "cgegs.f"
	    return 0;
#line 390 "cgegs.f"
	}
#line 391 "cgegs.f"
    }

/*     Permute the matrix to make it more nearly triangular */

#line 395 "cgegs.f"
    ileft = 1;
#line 396 "cgegs.f"
    iright = *n + 1;
#line 397 "cgegs.f"
    irwork = iright + *n;
#line 398 "cgegs.f"
    iwork = 1;
#line 399 "cgegs.f"
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwork], &iinfo, (ftnlen)1);
#line 401 "cgegs.f"
    if (iinfo != 0) {
#line 402 "cgegs.f"
	*info = *n + 1;
#line 403 "cgegs.f"
	goto L10;
#line 404 "cgegs.f"
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR */

#line 408 "cgegs.f"
    irows = ihi + 1 - ilo;
#line 409 "cgegs.f"
    icols = *n + 1 - ilo;
#line 410 "cgegs.f"
    itau = iwork;
#line 411 "cgegs.f"
    iwork = itau + irows;
#line 412 "cgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 412 "cgegs.f"
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
#line 414 "cgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 414 "cgegs.f"
	i__3 = iwork;
#line 414 "cgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 414 "cgegs.f"
	lwkopt = max(i__1,i__2);
#line 414 "cgegs.f"
    }
#line 416 "cgegs.f"
    if (iinfo != 0) {
#line 417 "cgegs.f"
	*info = *n + 2;
#line 418 "cgegs.f"
	goto L10;
#line 419 "cgegs.f"
    }

#line 421 "cgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 421 "cgegs.f"
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo, (ftnlen)1, (ftnlen)1);
#line 424 "cgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 424 "cgegs.f"
	i__3 = iwork;
#line 424 "cgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 424 "cgegs.f"
	lwkopt = max(i__1,i__2);
#line 424 "cgegs.f"
    }
#line 426 "cgegs.f"
    if (iinfo != 0) {
#line 427 "cgegs.f"
	*info = *n + 3;
#line 428 "cgegs.f"
	goto L10;
#line 429 "cgegs.f"
    }

#line 431 "cgegs.f"
    if (ilvsl) {
#line 432 "cgegs.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl, (ftnlen)
		4);
#line 433 "cgegs.f"
	i__1 = irows - 1;
#line 433 "cgegs.f"
	i__2 = irows - 1;
#line 433 "cgegs.f"
	clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[ilo 
		+ 1 + ilo * vsl_dim1], ldvsl, (ftnlen)1);
#line 435 "cgegs.f"
	i__1 = *lwork + 1 - iwork;
#line 435 "cgegs.f"
	cungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwork], &i__1, &iinfo);
#line 438 "cgegs.f"
	if (iinfo >= 0) {
/* Computing MAX */
#line 438 "cgegs.f"
	    i__3 = iwork;
#line 438 "cgegs.f"
	    i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 438 "cgegs.f"
	    lwkopt = max(i__1,i__2);
#line 438 "cgegs.f"
	}
#line 440 "cgegs.f"
	if (iinfo != 0) {
#line 441 "cgegs.f"
	    *info = *n + 4;
#line 442 "cgegs.f"
	    goto L10;
#line 443 "cgegs.f"
	}
#line 444 "cgegs.f"
    }

#line 446 "cgegs.f"
    if (ilvsr) {
#line 446 "cgegs.f"
	claset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr, (ftnlen)
		4);
#line 446 "cgegs.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 451 "cgegs.f"
    cgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &iinfo, (
	    ftnlen)1, (ftnlen)1);
#line 453 "cgegs.f"
    if (iinfo != 0) {
#line 454 "cgegs.f"
	*info = *n + 5;
#line 455 "cgegs.f"
	goto L10;
#line 456 "cgegs.f"
    }

/*     Perform QZ algorithm, computing Schur vectors if desired */

#line 460 "cgegs.f"
    iwork = itau;
#line 461 "cgegs.f"
    i__1 = *lwork + 1 - iwork;
#line 461 "cgegs.f"
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwork], &i__1, &rwork[irwork], &
	    iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 464 "cgegs.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 464 "cgegs.f"
	i__3 = iwork;
#line 464 "cgegs.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 464 "cgegs.f"
	lwkopt = max(i__1,i__2);
#line 464 "cgegs.f"
    }
#line 466 "cgegs.f"
    if (iinfo != 0) {
#line 467 "cgegs.f"
	if (iinfo > 0 && iinfo <= *n) {
#line 468 "cgegs.f"
	    *info = iinfo;
#line 469 "cgegs.f"
	} else if (iinfo > *n && iinfo <= *n << 1) {
#line 470 "cgegs.f"
	    *info = iinfo - *n;
#line 471 "cgegs.f"
	} else {
#line 472 "cgegs.f"
	    *info = *n + 6;
#line 473 "cgegs.f"
	}
#line 474 "cgegs.f"
	goto L10;
#line 475 "cgegs.f"
    }

/*     Apply permutation to VSL and VSR */

#line 479 "cgegs.f"
    if (ilvsl) {
#line 480 "cgegs.f"
	cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &iinfo, (ftnlen)1, (ftnlen)1);
#line 482 "cgegs.f"
	if (iinfo != 0) {
#line 483 "cgegs.f"
	    *info = *n + 7;
#line 484 "cgegs.f"
	    goto L10;
#line 485 "cgegs.f"
	}
#line 486 "cgegs.f"
    }
#line 487 "cgegs.f"
    if (ilvsr) {
#line 488 "cgegs.f"
	cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 490 "cgegs.f"
	if (iinfo != 0) {
#line 491 "cgegs.f"
	    *info = *n + 8;
#line 492 "cgegs.f"
	    goto L10;
#line 493 "cgegs.f"
	}
#line 494 "cgegs.f"
    }

/*     Undo scaling */

#line 498 "cgegs.f"
    if (ilascl) {
#line 499 "cgegs.f"
	clascl_("U", &c_n1, &c_n1, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 500 "cgegs.f"
	if (iinfo != 0) {
#line 501 "cgegs.f"
	    *info = *n + 9;
#line 502 "cgegs.f"
	    return 0;
#line 503 "cgegs.f"
	}
#line 504 "cgegs.f"
	clascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		iinfo, (ftnlen)1);
#line 505 "cgegs.f"
	if (iinfo != 0) {
#line 506 "cgegs.f"
	    *info = *n + 9;
#line 507 "cgegs.f"
	    return 0;
#line 508 "cgegs.f"
	}
#line 509 "cgegs.f"
    }

#line 511 "cgegs.f"
    if (ilbscl) {
#line 512 "cgegs.f"
	clascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 513 "cgegs.f"
	if (iinfo != 0) {
#line 514 "cgegs.f"
	    *info = *n + 9;
#line 515 "cgegs.f"
	    return 0;
#line 516 "cgegs.f"
	}
#line 517 "cgegs.f"
	clascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		iinfo, (ftnlen)1);
#line 518 "cgegs.f"
	if (iinfo != 0) {
#line 519 "cgegs.f"
	    *info = *n + 9;
#line 520 "cgegs.f"
	    return 0;
#line 521 "cgegs.f"
	}
#line 522 "cgegs.f"
    }

#line 524 "cgegs.f"
L10:
#line 525 "cgegs.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 527 "cgegs.f"
    return 0;

/*     End of CGEGS */

} /* cgegs_ */

