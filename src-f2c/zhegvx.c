#line 1 "zhegvx.f"
/* zhegvx.f -- translated by f2c (version 20100827).
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

#line 1 "zhegvx.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZHEGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, */
/*                          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
/*                          LWORK, RWORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEGVX computes selected eigenvalues, and optionally, eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and */
/* > B are assumed to be Hermitian and B is also positive definite. */
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
/* >          = 'U':  Upper triangles of A and B are stored; */
/* >          = 'L':  Lower triangles of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* > */
/* >          On exit,  the lower triangle (if UPLO='L') or the upper */
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
/* >          B is COMPLEX*16 array, dimension (LDB, N) */
/* >          On entry, the Hermitian matrix B.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of B contains the */
/* >          upper triangular part of the matrix B.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of B contains */
/* >          the lower triangular part of the matrix B. */
/* > */
/* >          On exit, if INFO <= N, the part of B containing the matrix is */
/* >          overwritten by the triangular factor U or L from the Cholesky */
/* >          factorization B = U**H*U or B = L*L**H. */
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
/* > */
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
/* > */
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
/* >          The first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= max(1,2*N). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the blocksize for ZHETRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (7*N) */
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
/* >          > 0:  ZPOTRF or ZHEEVX returned an error code: */
/* >             <= N:  if INFO = i, ZHEEVX failed to converge; */
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

/* > \ingroup complex16HEeigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int zhegvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
	 integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork,
	 integer *iwork, integer *ifail, integer *info, ftnlen jobz_len, 
	ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char trans[1];
    static logical upper, wantz;
    extern /* Subroutine */ int ztrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    ztrsm_(char *, char *, char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical alleig, indeig, valeig;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zhegst_(integer *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen), zheevx_(char *, char *, char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, doublereal *, doublecomplex *
	    , integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zpotrf_(char *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 344 "zhegvx.f"
    /* Parameter adjustments */
#line 344 "zhegvx.f"
    a_dim1 = *lda;
#line 344 "zhegvx.f"
    a_offset = 1 + a_dim1;
#line 344 "zhegvx.f"
    a -= a_offset;
#line 344 "zhegvx.f"
    b_dim1 = *ldb;
#line 344 "zhegvx.f"
    b_offset = 1 + b_dim1;
#line 344 "zhegvx.f"
    b -= b_offset;
#line 344 "zhegvx.f"
    --w;
#line 344 "zhegvx.f"
    z_dim1 = *ldz;
#line 344 "zhegvx.f"
    z_offset = 1 + z_dim1;
#line 344 "zhegvx.f"
    z__ -= z_offset;
#line 344 "zhegvx.f"
    --work;
#line 344 "zhegvx.f"
    --rwork;
#line 344 "zhegvx.f"
    --iwork;
#line 344 "zhegvx.f"
    --ifail;
#line 344 "zhegvx.f"

#line 344 "zhegvx.f"
    /* Function Body */
#line 344 "zhegvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 345 "zhegvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 346 "zhegvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 347 "zhegvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 348 "zhegvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 349 "zhegvx.f"
    lquery = *lwork == -1;

#line 351 "zhegvx.f"
    *info = 0;
#line 352 "zhegvx.f"
    if (*itype < 1 || *itype > 3) {
#line 353 "zhegvx.f"
	*info = -1;
#line 354 "zhegvx.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 355 "zhegvx.f"
	*info = -2;
#line 356 "zhegvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 357 "zhegvx.f"
	*info = -3;
#line 358 "zhegvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 359 "zhegvx.f"
	*info = -4;
#line 360 "zhegvx.f"
    } else if (*n < 0) {
#line 361 "zhegvx.f"
	*info = -5;
#line 362 "zhegvx.f"
    } else if (*lda < max(1,*n)) {
#line 363 "zhegvx.f"
	*info = -7;
#line 364 "zhegvx.f"
    } else if (*ldb < max(1,*n)) {
#line 365 "zhegvx.f"
	*info = -9;
#line 366 "zhegvx.f"
    } else {
#line 367 "zhegvx.f"
	if (valeig) {
#line 368 "zhegvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 368 "zhegvx.f"
		*info = -11;
#line 368 "zhegvx.f"
	    }
#line 370 "zhegvx.f"
	} else if (indeig) {
#line 371 "zhegvx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 372 "zhegvx.f"
		*info = -12;
#line 373 "zhegvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 374 "zhegvx.f"
		*info = -13;
#line 375 "zhegvx.f"
	    }
#line 376 "zhegvx.f"
	}
#line 377 "zhegvx.f"
    }
#line 378 "zhegvx.f"
    if (*info == 0) {
#line 379 "zhegvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 380 "zhegvx.f"
	    *info = -18;
#line 381 "zhegvx.f"
	}
#line 382 "zhegvx.f"
    }

#line 384 "zhegvx.f"
    if (*info == 0) {
#line 385 "zhegvx.f"
	nb = ilaenv_(&c__1, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 386 "zhegvx.f"
	i__1 = 1, i__2 = (nb + 1) * *n;
#line 386 "zhegvx.f"
	lwkopt = max(i__1,i__2);
#line 387 "zhegvx.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

/* Computing MAX */
#line 389 "zhegvx.f"
	i__1 = 1, i__2 = *n << 1;
#line 389 "zhegvx.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 390 "zhegvx.f"
	    *info = -20;
#line 391 "zhegvx.f"
	}
#line 392 "zhegvx.f"
    }

#line 394 "zhegvx.f"
    if (*info != 0) {
#line 395 "zhegvx.f"
	i__1 = -(*info);
#line 395 "zhegvx.f"
	xerbla_("ZHEGVX", &i__1, (ftnlen)6);
#line 396 "zhegvx.f"
	return 0;
#line 397 "zhegvx.f"
    } else if (lquery) {
#line 398 "zhegvx.f"
	return 0;
#line 399 "zhegvx.f"
    }

/*     Quick return if possible */

#line 403 "zhegvx.f"
    *m = 0;
#line 404 "zhegvx.f"
    if (*n == 0) {
#line 405 "zhegvx.f"
	return 0;
#line 406 "zhegvx.f"
    }

/*     Form a Cholesky factorization of B. */

#line 410 "zhegvx.f"
    zpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 411 "zhegvx.f"
    if (*info != 0) {
#line 412 "zhegvx.f"
	*info = *n + *info;
#line 413 "zhegvx.f"
	return 0;
#line 414 "zhegvx.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 418 "zhegvx.f"
    zhegst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 419 "zhegvx.f"
    zheevx_(jobz, range, uplo, n, &a[a_offset], lda, vl, vu, il, iu, abstol, 
	    m, &w[1], &z__[z_offset], ldz, &work[1], lwork, &rwork[1], &iwork[
	    1], &ifail[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 423 "zhegvx.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 427 "zhegvx.f"
	if (*info > 0) {
#line 427 "zhegvx.f"
	    *m = *info - 1;
#line 427 "zhegvx.f"
	}
#line 429 "zhegvx.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y */

#line 434 "zhegvx.f"
	    if (upper) {
#line 435 "zhegvx.f"
		*(unsigned char *)trans = 'N';
#line 436 "zhegvx.f"
	    } else {
#line 437 "zhegvx.f"
		*(unsigned char *)trans = 'C';
#line 438 "zhegvx.f"
	    }

#line 440 "zhegvx.f"
	    ztrsm_("Left", uplo, trans, "Non-unit", n, m, &c_b1, &b[b_offset],
		     ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);

#line 443 "zhegvx.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H *y */

#line 448 "zhegvx.f"
	    if (upper) {
#line 449 "zhegvx.f"
		*(unsigned char *)trans = 'C';
#line 450 "zhegvx.f"
	    } else {
#line 451 "zhegvx.f"
		*(unsigned char *)trans = 'N';
#line 452 "zhegvx.f"
	    }

#line 454 "zhegvx.f"
	    ztrmm_("Left", uplo, trans, "Non-unit", n, m, &c_b1, &b[b_offset],
		     ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);
#line 456 "zhegvx.f"
	}
#line 457 "zhegvx.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 461 "zhegvx.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 463 "zhegvx.f"
    return 0;

/*     End of ZHEGVX */

} /* zhegvx_ */

