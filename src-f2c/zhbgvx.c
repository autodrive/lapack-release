#line 1 "zhbgvx.f"
/* zhbgvx.f -- translated by f2c (version 20100827).
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

#line 1 "zhbgvx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHBGVX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHBGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbgvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbgvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbgvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, */
/*                          LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, */
/*                          LDZ, WORK, RWORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M, */
/*      $                   N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBGVX computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite banded eigenproblem, of */
/* > the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian */
/* > and banded, and B is also positive definite.  Eigenvalues and */
/* > eigenvectors can be selected by specifying either all eigenvalues, */
/* > a range of values or a range of indices for the desired eigenvalues. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* >          = 'A': all eigenvalues will be found; */
/* >          = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* >                 will be found; */
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
/* > \param[in] KA */
/* > \verbatim */
/* >          KA is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'. KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of superdiagonals of the matrix B if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'. KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first ka+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). */
/* > */
/* >          On exit, the contents of AB are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KA+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BB */
/* > \verbatim */
/* >          BB is COMPLEX*16 array, dimension (LDBB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix B, stored in the first kb+1 rows of the array.  The */
/* >          j-th column of B is stored in the j-th column of the array BB */
/* >          as follows: */
/* >          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; */
/* >          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). */
/* > */
/* >          On exit, the factor S from the split Cholesky factorization */
/* >          B = S**H*S, as returned by ZPBSTF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* >          LDBB is INTEGER */
/* >          The leading dimension of the array BB.  LDBB >= KB+1. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ, N) */
/* >          If JOBZ = 'V', the n-by-n matrix used in the reduction of */
/* >          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x, */
/* >          and consequently C to tridiagonal form. */
/* >          If JOBZ = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  If JOBZ = 'N', */
/* >          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION */
/* > */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* > */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* > */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* > */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
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
/* >          by reducing AP to tridiagonal form. */
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
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors, with the i-th column of Z holding the */
/* >          eigenvector associated with W(i). The eigenvectors are */
/* >          normalized so that Z**H*B*Z = I. */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
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
/* >          > 0:  if INFO = i, and i is: */
/* >             <= N:  then i eigenvectors failed to converge.  Their */
/* >                    indices are stored in array IFAIL. */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF */
/* >                    returned INFO = i: B is not positive definite. */
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

/* > \ingroup complex16OTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int zhbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, doublecomplex *ab, integer *ldab, 
	doublecomplex *bb, integer *ldbb, doublecomplex *q, integer *ldq, 
	doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *
	abstol, integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info, ftnlen jobz_len, ftnlen range_len, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, jj;
    static doublereal tmp1;
    static integer indd, inde;
    static char vect[1];
    static logical test;
    static integer itmp1, indee;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical upper, wantz;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static logical alleig, indeig;
    static integer indibl;
    static logical valeig;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer indiwk, indisp;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dstebz_(char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    zhbtrd_(char *, char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int zhbgst_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, ftnlen, ftnlen), zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer nsplit;
    extern /* Subroutine */ int zpbstf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen), zstein_(integer *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    integer *, integer *), zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 353 "zhbgvx.f"
    /* Parameter adjustments */
#line 353 "zhbgvx.f"
    ab_dim1 = *ldab;
#line 353 "zhbgvx.f"
    ab_offset = 1 + ab_dim1;
#line 353 "zhbgvx.f"
    ab -= ab_offset;
#line 353 "zhbgvx.f"
    bb_dim1 = *ldbb;
#line 353 "zhbgvx.f"
    bb_offset = 1 + bb_dim1;
#line 353 "zhbgvx.f"
    bb -= bb_offset;
#line 353 "zhbgvx.f"
    q_dim1 = *ldq;
#line 353 "zhbgvx.f"
    q_offset = 1 + q_dim1;
#line 353 "zhbgvx.f"
    q -= q_offset;
#line 353 "zhbgvx.f"
    --w;
#line 353 "zhbgvx.f"
    z_dim1 = *ldz;
#line 353 "zhbgvx.f"
    z_offset = 1 + z_dim1;
#line 353 "zhbgvx.f"
    z__ -= z_offset;
#line 353 "zhbgvx.f"
    --work;
#line 353 "zhbgvx.f"
    --rwork;
#line 353 "zhbgvx.f"
    --iwork;
#line 353 "zhbgvx.f"
    --ifail;
#line 353 "zhbgvx.f"

#line 353 "zhbgvx.f"
    /* Function Body */
#line 353 "zhbgvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 354 "zhbgvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 355 "zhbgvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 356 "zhbgvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 357 "zhbgvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 359 "zhbgvx.f"
    *info = 0;
#line 360 "zhbgvx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 361 "zhbgvx.f"
	*info = -1;
#line 362 "zhbgvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 363 "zhbgvx.f"
	*info = -2;
#line 364 "zhbgvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 365 "zhbgvx.f"
	*info = -3;
#line 366 "zhbgvx.f"
    } else if (*n < 0) {
#line 367 "zhbgvx.f"
	*info = -4;
#line 368 "zhbgvx.f"
    } else if (*ka < 0) {
#line 369 "zhbgvx.f"
	*info = -5;
#line 370 "zhbgvx.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 371 "zhbgvx.f"
	*info = -6;
#line 372 "zhbgvx.f"
    } else if (*ldab < *ka + 1) {
#line 373 "zhbgvx.f"
	*info = -8;
#line 374 "zhbgvx.f"
    } else if (*ldbb < *kb + 1) {
#line 375 "zhbgvx.f"
	*info = -10;
#line 376 "zhbgvx.f"
    } else if (*ldq < 1 || wantz && *ldq < *n) {
#line 377 "zhbgvx.f"
	*info = -12;
#line 378 "zhbgvx.f"
    } else {
#line 379 "zhbgvx.f"
	if (valeig) {
#line 380 "zhbgvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 380 "zhbgvx.f"
		*info = -14;
#line 380 "zhbgvx.f"
	    }
#line 382 "zhbgvx.f"
	} else if (indeig) {
#line 383 "zhbgvx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 384 "zhbgvx.f"
		*info = -15;
#line 385 "zhbgvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 386 "zhbgvx.f"
		*info = -16;
#line 387 "zhbgvx.f"
	    }
#line 388 "zhbgvx.f"
	}
#line 389 "zhbgvx.f"
    }
#line 390 "zhbgvx.f"
    if (*info == 0) {
#line 391 "zhbgvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 392 "zhbgvx.f"
	    *info = -21;
#line 393 "zhbgvx.f"
	}
#line 394 "zhbgvx.f"
    }

#line 396 "zhbgvx.f"
    if (*info != 0) {
#line 397 "zhbgvx.f"
	i__1 = -(*info);
#line 397 "zhbgvx.f"
	xerbla_("ZHBGVX", &i__1, (ftnlen)6);
#line 398 "zhbgvx.f"
	return 0;
#line 399 "zhbgvx.f"
    }

/*     Quick return if possible */

#line 403 "zhbgvx.f"
    *m = 0;
#line 404 "zhbgvx.f"
    if (*n == 0) {
#line 404 "zhbgvx.f"
	return 0;
#line 404 "zhbgvx.f"
    }

/*     Form a split Cholesky factorization of B. */

#line 409 "zhbgvx.f"
    zpbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info, (ftnlen)1);
#line 410 "zhbgvx.f"
    if (*info != 0) {
#line 411 "zhbgvx.f"
	*info = *n + *info;
#line 412 "zhbgvx.f"
	return 0;
#line 413 "zhbgvx.f"
    }

/*     Transform problem to standard eigenvalue problem. */

#line 417 "zhbgvx.f"
    zhbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
	     &q[q_offset], ldq, &work[1], &rwork[1], &iinfo, (ftnlen)1, (
	    ftnlen)1);

/*     Solve the standard eigenvalue problem. */
/*     Reduce Hermitian band matrix to tridiagonal form. */

#line 423 "zhbgvx.f"
    indd = 1;
#line 424 "zhbgvx.f"
    inde = indd + *n;
#line 425 "zhbgvx.f"
    indrwk = inde + *n;
#line 426 "zhbgvx.f"
    indwrk = 1;
#line 427 "zhbgvx.f"
    if (wantz) {
#line 428 "zhbgvx.f"
	*(unsigned char *)vect = 'U';
#line 429 "zhbgvx.f"
    } else {
#line 430 "zhbgvx.f"
	*(unsigned char *)vect = 'N';
#line 431 "zhbgvx.f"
    }
#line 432 "zhbgvx.f"
    zhbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &rwork[indd], &rwork[
	    inde], &q[q_offset], ldq, &work[indwrk], &iinfo, (ftnlen)1, (
	    ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call DSTERF or ZSTEQR.  If this fails for some */
/*     eigenvalue, then try DSTEBZ. */

#line 439 "zhbgvx.f"
    test = FALSE_;
#line 440 "zhbgvx.f"
    if (indeig) {
#line 441 "zhbgvx.f"
	if (*il == 1 && *iu == *n) {
#line 442 "zhbgvx.f"
	    test = TRUE_;
#line 443 "zhbgvx.f"
	}
#line 444 "zhbgvx.f"
    }
#line 445 "zhbgvx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 446 "zhbgvx.f"
	dcopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 447 "zhbgvx.f"
	indee = indrwk + (*n << 1);
#line 448 "zhbgvx.f"
	i__1 = *n - 1;
#line 448 "zhbgvx.f"
	dcopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 449 "zhbgvx.f"
	if (! wantz) {
#line 450 "zhbgvx.f"
	    dsterf_(n, &w[1], &rwork[indee], info);
#line 451 "zhbgvx.f"
	} else {
#line 452 "zhbgvx.f"
	    zlacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 453 "zhbgvx.f"
	    zsteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 455 "zhbgvx.f"
	    if (*info == 0) {
#line 456 "zhbgvx.f"
		i__1 = *n;
#line 456 "zhbgvx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 457 "zhbgvx.f"
		    ifail[i__] = 0;
#line 458 "zhbgvx.f"
/* L10: */
#line 458 "zhbgvx.f"
		}
#line 459 "zhbgvx.f"
	    }
#line 460 "zhbgvx.f"
	}
#line 461 "zhbgvx.f"
	if (*info == 0) {
#line 462 "zhbgvx.f"
	    *m = *n;
#line 463 "zhbgvx.f"
	    goto L30;
#line 464 "zhbgvx.f"
	}
#line 465 "zhbgvx.f"
	*info = 0;
#line 466 "zhbgvx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, */
/*     call ZSTEIN. */

#line 471 "zhbgvx.f"
    if (wantz) {
#line 472 "zhbgvx.f"
	*(unsigned char *)order = 'B';
#line 473 "zhbgvx.f"
    } else {
#line 474 "zhbgvx.f"
	*(unsigned char *)order = 'E';
#line 475 "zhbgvx.f"
    }
#line 476 "zhbgvx.f"
    indibl = 1;
#line 477 "zhbgvx.f"
    indisp = indibl + *n;
#line 478 "zhbgvx.f"
    indiwk = indisp + *n;
#line 479 "zhbgvx.f"
    dstebz_(range, order, n, vl, vu, il, iu, abstol, &rwork[indd], &rwork[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &rwork[
	    indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 484 "zhbgvx.f"
    if (wantz) {
#line 485 "zhbgvx.f"
	zstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by ZSTEIN. */

#line 492 "zhbgvx.f"
	i__1 = *m;
#line 492 "zhbgvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 493 "zhbgvx.f"
	    zcopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 494 "zhbgvx.f"
	    zgemv_("N", n, n, &c_b2, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b1, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 496 "zhbgvx.f"
/* L20: */
#line 496 "zhbgvx.f"
	}
#line 497 "zhbgvx.f"
    }

#line 499 "zhbgvx.f"
L30:

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 504 "zhbgvx.f"
    if (wantz) {
#line 505 "zhbgvx.f"
	i__1 = *m - 1;
#line 505 "zhbgvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 506 "zhbgvx.f"
	    i__ = 0;
#line 507 "zhbgvx.f"
	    tmp1 = w[j];
#line 508 "zhbgvx.f"
	    i__2 = *m;
#line 508 "zhbgvx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 509 "zhbgvx.f"
		if (w[jj] < tmp1) {
#line 510 "zhbgvx.f"
		    i__ = jj;
#line 511 "zhbgvx.f"
		    tmp1 = w[jj];
#line 512 "zhbgvx.f"
		}
#line 513 "zhbgvx.f"
/* L40: */
#line 513 "zhbgvx.f"
	    }

#line 515 "zhbgvx.f"
	    if (i__ != 0) {
#line 516 "zhbgvx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 517 "zhbgvx.f"
		w[i__] = w[j];
#line 518 "zhbgvx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 519 "zhbgvx.f"
		w[j] = tmp1;
#line 520 "zhbgvx.f"
		iwork[indibl + j - 1] = itmp1;
#line 521 "zhbgvx.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 522 "zhbgvx.f"
		if (*info != 0) {
#line 523 "zhbgvx.f"
		    itmp1 = ifail[i__];
#line 524 "zhbgvx.f"
		    ifail[i__] = ifail[j];
#line 525 "zhbgvx.f"
		    ifail[j] = itmp1;
#line 526 "zhbgvx.f"
		}
#line 527 "zhbgvx.f"
	    }
#line 528 "zhbgvx.f"
/* L50: */
#line 528 "zhbgvx.f"
	}
#line 529 "zhbgvx.f"
    }

#line 531 "zhbgvx.f"
    return 0;

/*     End of ZHBGVX */

} /* zhbgvx_ */

