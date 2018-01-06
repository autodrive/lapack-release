#line 1 "dsygvx.f"
/* dsygvx.f -- translated by f2c (version 20100827).
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

#line 1 "dsygvx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b19 = 1.;

/* > \brief \b DSYGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, */
/*                          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
/*                          LWORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYGVX computes selected eigenvalues, and optionally, eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A */
/* > and B are assumed to be symmetric and B is also positive definite. */
/* > Eigenvalues and eigenvectors can be selected by specifying either a */
/* > range of values or a range of indices for the desired eigenvalues. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          Specifies the problem type to be solved: */
/* >          = 1:  A*x = (lambda)*B*x */
/* >          = 2:  A*B*x = (lambda)*x */
/* >          = 3:  B*A*x = (lambda)*x */
/* > \endverbatim */
/* > */
/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only; */
/* >          = 'V':  Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RANGE */
/* > \verbatim */
/* >          RANGE is CHARACTER*1 */
/* >          = 'A': all eigenvalues will be found. */
/* >          = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* >                 will be found. */
/* >          = 'I': the IL-th through IU-th eigenvalues will be found. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A and B are stored; */
/* >          = 'L':  Lower triangle of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix pencil (A,B).  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* > */
/* >          On exit, the lower triangle (if UPLO='L') or the upper */
/* >          triangle (if UPLO='U') of A, including the diagonal, is */
/* >          destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
/* >          On entry, the symmetric matrix B.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of B contains the */
/* >          upper triangular part of the matrix B.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of B contains */
/* >          the lower triangular part of the matrix B. */
/* > */
/* >          On exit, if INFO <= N, the part of B containing the matrix is */
/* >          overwritten by the triangular factor U or L from the Cholesky */
/* >          factorization B = U**T*U or B = L*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          If RANGE='V', the lower and upper bounds of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the indices (in ascending order) of the */
/* >          smallest and largest eigenvalues to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* >          ABSTOL is DOUBLE PRECISION */
/* >          The absolute error tolerance for the eigenvalues. */
/* >          An approximate eigenvalue is accepted as converged */
/* >          when it is determined to lie in an interval [a,b] */
/* >          of width less than or equal to */
/* > */
/* >                  ABSTOL + EPS *   max( |a|,|b| ) , */
/* > */
/* >          where EPS is the machine precision.  If ABSTOL is less than */
/* >          or equal to zero, then  EPS*|T|  will be used in its place, */
/* >          where |T| is the 1-norm of the tridiagonal matrix obtained */
/* >          by reducing C to tridiagonal form, where C is the symmetric */
/* >          matrix of the standard symmetric problem to which the */
/* >          generalized problem is transformed. */
/* > */
/* >          Eigenvalues will be computed most accurately when ABSTOL is */
/* >          set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
/* >          If this routine returns with INFO>0, indicating that some */
/* >          eigenvectors did not converge, try setting ABSTOL to */
/* >          2*DLAMCH('S'). */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The total number of eigenvalues found.  0 <= M <= N. */
/* >          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M)) */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          The eigenvectors are normalized as follows: */
/* >          if ITYPE = 1 or 2, Z**T*B*Z = I; */
/* >          if ITYPE = 3, Z**T*inv(B)*Z = I. */
/* > */
/* >          If an eigenvector fails to converge, then that column of Z */
/* >          contains the latest approximation to the eigenvector, and the */
/* >          index of the eigenvector is returned in IFAIL. */
/* >          Note: the user must ensure that at least max(1,M) columns are */
/* >          supplied in the array Z; if RANGE = 'V', the exact value of M */
/* >          is not known in advance and an upper bound must be used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= max(1,N). */
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
/* >          The length of the array WORK.  LWORK >= max(1,8*N). */
/* >          For optimal efficiency, LWORK >= (NB+3)*N, */
/* >          where NB is the blocksize for DSYTRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* >          IFAIL is INTEGER array, dimension (N) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M elements of */
/* >          IFAIL are zero.  If INFO > 0, then IFAIL contains the */
/* >          indices of the eigenvectors that failed to converge. */
/* >          If JOBZ = 'N', then IFAIL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  DPOTRF or DSYEVX returned an error code: */
/* >             <= N:  if INFO = i, DSYEVX failed to converge; */
/* >                    i eigenvectors failed to converge.  Their indices */
/* >                    are stored in array IFAIL. */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then the leading */
/* >                    minor of order i of B is not positive definite. */
/* >                    The factorization of B could not be completed and */
/* >                    no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleSYeigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int dsygvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *a, integer *lda, doublereal *b, integer 
	*ldb, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *ifail, integer *info, ftnlen jobz_len, ftnlen range_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static char trans[1];
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper, wantz, alleig, indeig, valeig;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer lwkmin;
    extern /* Subroutine */ int dsygst_(integer *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int dsyevx_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 335 "dsygvx.f"
    /* Parameter adjustments */
#line 335 "dsygvx.f"
    a_dim1 = *lda;
#line 335 "dsygvx.f"
    a_offset = 1 + a_dim1;
#line 335 "dsygvx.f"
    a -= a_offset;
#line 335 "dsygvx.f"
    b_dim1 = *ldb;
#line 335 "dsygvx.f"
    b_offset = 1 + b_dim1;
#line 335 "dsygvx.f"
    b -= b_offset;
#line 335 "dsygvx.f"
    --w;
#line 335 "dsygvx.f"
    z_dim1 = *ldz;
#line 335 "dsygvx.f"
    z_offset = 1 + z_dim1;
#line 335 "dsygvx.f"
    z__ -= z_offset;
#line 335 "dsygvx.f"
    --work;
#line 335 "dsygvx.f"
    --iwork;
#line 335 "dsygvx.f"
    --ifail;
#line 335 "dsygvx.f"

#line 335 "dsygvx.f"
    /* Function Body */
#line 335 "dsygvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 336 "dsygvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 337 "dsygvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 338 "dsygvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 339 "dsygvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 340 "dsygvx.f"
    lquery = *lwork == -1;

#line 342 "dsygvx.f"
    *info = 0;
#line 343 "dsygvx.f"
    if (*itype < 1 || *itype > 3) {
#line 344 "dsygvx.f"
	*info = -1;
#line 345 "dsygvx.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 346 "dsygvx.f"
	*info = -2;
#line 347 "dsygvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 348 "dsygvx.f"
	*info = -3;
#line 349 "dsygvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 350 "dsygvx.f"
	*info = -4;
#line 351 "dsygvx.f"
    } else if (*n < 0) {
#line 352 "dsygvx.f"
	*info = -5;
#line 353 "dsygvx.f"
    } else if (*lda < max(1,*n)) {
#line 354 "dsygvx.f"
	*info = -7;
#line 355 "dsygvx.f"
    } else if (*ldb < max(1,*n)) {
#line 356 "dsygvx.f"
	*info = -9;
#line 357 "dsygvx.f"
    } else {
#line 358 "dsygvx.f"
	if (valeig) {
#line 359 "dsygvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 359 "dsygvx.f"
		*info = -11;
#line 359 "dsygvx.f"
	    }
#line 361 "dsygvx.f"
	} else if (indeig) {
#line 362 "dsygvx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 363 "dsygvx.f"
		*info = -12;
#line 364 "dsygvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 365 "dsygvx.f"
		*info = -13;
#line 366 "dsygvx.f"
	    }
#line 367 "dsygvx.f"
	}
#line 368 "dsygvx.f"
    }
#line 369 "dsygvx.f"
    if (*info == 0) {
#line 370 "dsygvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 371 "dsygvx.f"
	    *info = -18;
#line 372 "dsygvx.f"
	}
#line 373 "dsygvx.f"
    }

#line 375 "dsygvx.f"
    if (*info == 0) {
/* Computing MAX */
#line 376 "dsygvx.f"
	i__1 = 1, i__2 = *n << 3;
#line 376 "dsygvx.f"
	lwkmin = max(i__1,i__2);
#line 377 "dsygvx.f"
	nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 378 "dsygvx.f"
	i__1 = lwkmin, i__2 = (nb + 3) * *n;
#line 378 "dsygvx.f"
	lwkopt = max(i__1,i__2);
#line 379 "dsygvx.f"
	work[1] = (doublereal) lwkopt;

#line 381 "dsygvx.f"
	if (*lwork < lwkmin && ! lquery) {
#line 382 "dsygvx.f"
	    *info = -20;
#line 383 "dsygvx.f"
	}
#line 384 "dsygvx.f"
    }

#line 386 "dsygvx.f"
    if (*info != 0) {
#line 387 "dsygvx.f"
	i__1 = -(*info);
#line 387 "dsygvx.f"
	xerbla_("DSYGVX", &i__1, (ftnlen)6);
#line 388 "dsygvx.f"
	return 0;
#line 389 "dsygvx.f"
    } else if (lquery) {
#line 390 "dsygvx.f"
	return 0;
#line 391 "dsygvx.f"
    }

/*     Quick return if possible */

#line 395 "dsygvx.f"
    *m = 0;
#line 396 "dsygvx.f"
    if (*n == 0) {
#line 397 "dsygvx.f"
	return 0;
#line 398 "dsygvx.f"
    }

/*     Form a Cholesky factorization of B. */

#line 402 "dsygvx.f"
    dpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 403 "dsygvx.f"
    if (*info != 0) {
#line 404 "dsygvx.f"
	*info = *n + *info;
#line 405 "dsygvx.f"
	return 0;
#line 406 "dsygvx.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 410 "dsygvx.f"
    dsygst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 411 "dsygvx.f"
    dsyevx_(jobz, range, uplo, n, &a[a_offset], lda, vl, vu, il, iu, abstol, 
	    m, &w[1], &z__[z_offset], ldz, &work[1], lwork, &iwork[1], &ifail[
	    1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 414 "dsygvx.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 418 "dsygvx.f"
	if (*info > 0) {
#line 418 "dsygvx.f"
	    *m = *info - 1;
#line 418 "dsygvx.f"
	}
#line 420 "dsygvx.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */

#line 425 "dsygvx.f"
	    if (upper) {
#line 426 "dsygvx.f"
		*(unsigned char *)trans = 'N';
#line 427 "dsygvx.f"
	    } else {
#line 428 "dsygvx.f"
		*(unsigned char *)trans = 'T';
#line 429 "dsygvx.f"
	    }

#line 431 "dsygvx.f"
	    dtrsm_("Left", uplo, trans, "Non-unit", n, m, &c_b19, &b[b_offset]
		    , ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);

#line 434 "dsygvx.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T*y */

#line 439 "dsygvx.f"
	    if (upper) {
#line 440 "dsygvx.f"
		*(unsigned char *)trans = 'T';
#line 441 "dsygvx.f"
	    } else {
#line 442 "dsygvx.f"
		*(unsigned char *)trans = 'N';
#line 443 "dsygvx.f"
	    }

#line 445 "dsygvx.f"
	    dtrmm_("Left", uplo, trans, "Non-unit", n, m, &c_b19, &b[b_offset]
		    , ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);
#line 447 "dsygvx.f"
	}
#line 448 "dsygvx.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 452 "dsygvx.f"
    work[1] = (doublereal) lwkopt;

#line 454 "dsygvx.f"
    return 0;

/*     End of DSYGVX */

} /* dsygvx_ */

