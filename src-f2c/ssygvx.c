#line 1 "ssygvx.f"
/* ssygvx.f -- translated by f2c (version 20100827).
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

#line 1 "ssygvx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b19 = 1.;

/* > \brief \b SSYGVX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, */
/*                          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
/*                          LWORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYGVX computes selected eigenvalues, and optionally, eigenvectors */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          B is REAL array, dimension (LDB, N) */
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
/* >          VL is REAL */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* >          ABSTOL is REAL */
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
/* >          2*SLAMCH('S'). */
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
/* >          W is REAL array, dimension (N) */
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= max(1,8*N). */
/* >          For optimal efficiency, LWORK >= (NB+3)*N, */
/* >          where NB is the blocksize for SSYTRD returned by ILAENV. */
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
/* >          > 0:  SPOTRF or SSYEVX returned an error code: */
/* >             <= N:  if INFO = i, SSYEVX failed to converge; */
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

/* > \date June 2016 */

/* > \ingroup realSYeigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int ssygvx_(integer *itype, char *jobz, char *range, char *
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
    static char trans[1];
    static logical upper;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical wantz;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical alleig, indeig, valeig;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkmin;
    extern /* Subroutine */ int spotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int ssygst_(integer *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), ssyevx_(char *, char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 342 "ssygvx.f"
    /* Parameter adjustments */
#line 342 "ssygvx.f"
    a_dim1 = *lda;
#line 342 "ssygvx.f"
    a_offset = 1 + a_dim1;
#line 342 "ssygvx.f"
    a -= a_offset;
#line 342 "ssygvx.f"
    b_dim1 = *ldb;
#line 342 "ssygvx.f"
    b_offset = 1 + b_dim1;
#line 342 "ssygvx.f"
    b -= b_offset;
#line 342 "ssygvx.f"
    --w;
#line 342 "ssygvx.f"
    z_dim1 = *ldz;
#line 342 "ssygvx.f"
    z_offset = 1 + z_dim1;
#line 342 "ssygvx.f"
    z__ -= z_offset;
#line 342 "ssygvx.f"
    --work;
#line 342 "ssygvx.f"
    --iwork;
#line 342 "ssygvx.f"
    --ifail;
#line 342 "ssygvx.f"

#line 342 "ssygvx.f"
    /* Function Body */
#line 342 "ssygvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 343 "ssygvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 344 "ssygvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 345 "ssygvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 346 "ssygvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 347 "ssygvx.f"
    lquery = *lwork == -1;

#line 349 "ssygvx.f"
    *info = 0;
#line 350 "ssygvx.f"
    if (*itype < 1 || *itype > 3) {
#line 351 "ssygvx.f"
	*info = -1;
#line 352 "ssygvx.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 353 "ssygvx.f"
	*info = -2;
#line 354 "ssygvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 355 "ssygvx.f"
	*info = -3;
#line 356 "ssygvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 357 "ssygvx.f"
	*info = -4;
#line 358 "ssygvx.f"
    } else if (*n < 0) {
#line 359 "ssygvx.f"
	*info = -5;
#line 360 "ssygvx.f"
    } else if (*lda < max(1,*n)) {
#line 361 "ssygvx.f"
	*info = -7;
#line 362 "ssygvx.f"
    } else if (*ldb < max(1,*n)) {
#line 363 "ssygvx.f"
	*info = -9;
#line 364 "ssygvx.f"
    } else {
#line 365 "ssygvx.f"
	if (valeig) {
#line 366 "ssygvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 366 "ssygvx.f"
		*info = -11;
#line 366 "ssygvx.f"
	    }
#line 368 "ssygvx.f"
	} else if (indeig) {
#line 369 "ssygvx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 370 "ssygvx.f"
		*info = -12;
#line 371 "ssygvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 372 "ssygvx.f"
		*info = -13;
#line 373 "ssygvx.f"
	    }
#line 374 "ssygvx.f"
	}
#line 375 "ssygvx.f"
    }
#line 376 "ssygvx.f"
    if (*info == 0) {
#line 377 "ssygvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 378 "ssygvx.f"
	    *info = -18;
#line 379 "ssygvx.f"
	}
#line 380 "ssygvx.f"
    }

#line 382 "ssygvx.f"
    if (*info == 0) {
/* Computing MAX */
#line 383 "ssygvx.f"
	i__1 = 1, i__2 = *n << 3;
#line 383 "ssygvx.f"
	lwkmin = max(i__1,i__2);
#line 384 "ssygvx.f"
	nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 385 "ssygvx.f"
	i__1 = lwkmin, i__2 = (nb + 3) * *n;
#line 385 "ssygvx.f"
	lwkopt = max(i__1,i__2);
#line 386 "ssygvx.f"
	work[1] = (doublereal) lwkopt;

#line 388 "ssygvx.f"
	if (*lwork < lwkmin && ! lquery) {
#line 389 "ssygvx.f"
	    *info = -20;
#line 390 "ssygvx.f"
	}
#line 391 "ssygvx.f"
    }

#line 393 "ssygvx.f"
    if (*info != 0) {
#line 394 "ssygvx.f"
	i__1 = -(*info);
#line 394 "ssygvx.f"
	xerbla_("SSYGVX", &i__1, (ftnlen)6);
#line 395 "ssygvx.f"
	return 0;
#line 396 "ssygvx.f"
    } else if (lquery) {
#line 397 "ssygvx.f"
	return 0;
#line 398 "ssygvx.f"
    }

/*     Quick return if possible */

#line 402 "ssygvx.f"
    *m = 0;
#line 403 "ssygvx.f"
    if (*n == 0) {
#line 404 "ssygvx.f"
	return 0;
#line 405 "ssygvx.f"
    }

/*     Form a Cholesky factorization of B. */

#line 409 "ssygvx.f"
    spotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 410 "ssygvx.f"
    if (*info != 0) {
#line 411 "ssygvx.f"
	*info = *n + *info;
#line 412 "ssygvx.f"
	return 0;
#line 413 "ssygvx.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 417 "ssygvx.f"
    ssygst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 418 "ssygvx.f"
    ssyevx_(jobz, range, uplo, n, &a[a_offset], lda, vl, vu, il, iu, abstol, 
	    m, &w[1], &z__[z_offset], ldz, &work[1], lwork, &iwork[1], &ifail[
	    1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 421 "ssygvx.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 425 "ssygvx.f"
	if (*info > 0) {
#line 425 "ssygvx.f"
	    *m = *info - 1;
#line 425 "ssygvx.f"
	}
#line 427 "ssygvx.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */

#line 432 "ssygvx.f"
	    if (upper) {
#line 433 "ssygvx.f"
		*(unsigned char *)trans = 'N';
#line 434 "ssygvx.f"
	    } else {
#line 435 "ssygvx.f"
		*(unsigned char *)trans = 'T';
#line 436 "ssygvx.f"
	    }

#line 438 "ssygvx.f"
	    strsm_("Left", uplo, trans, "Non-unit", n, m, &c_b19, &b[b_offset]
		    , ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);

#line 441 "ssygvx.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T*y */

#line 446 "ssygvx.f"
	    if (upper) {
#line 447 "ssygvx.f"
		*(unsigned char *)trans = 'T';
#line 448 "ssygvx.f"
	    } else {
#line 449 "ssygvx.f"
		*(unsigned char *)trans = 'N';
#line 450 "ssygvx.f"
	    }

#line 452 "ssygvx.f"
	    strmm_("Left", uplo, trans, "Non-unit", n, m, &c_b19, &b[b_offset]
		    , ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);
#line 454 "ssygvx.f"
	}
#line 455 "ssygvx.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 459 "ssygvx.f"
    work[1] = (doublereal) lwkopt;

#line 461 "ssygvx.f"
    return 0;

/*     End of SSYGVX */

} /* ssygvx_ */

