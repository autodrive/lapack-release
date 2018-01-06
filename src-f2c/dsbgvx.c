#line 1 "dsbgvx.f"
/* dsbgvx.f -- translated by f2c (version 20100827).
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

#line 1 "dsbgvx.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b25 = 1.;
static doublereal c_b27 = 0.;

/* > \brief \b DSBGVX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, */
/*                          LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, */
/*                          LDZ, WORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M, */
/*      $                   N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ), */
/*      $                   W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBGVX computes selected eigenvalues, and optionally, eigenvectors */
/* > of a real generalized symmetric-definite banded eigenproblem, of */
/* > the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric */
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
/* > \param[in] KA */
/* > \verbatim */
/* >          KA is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of superdiagonals of the matrix B if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
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
/* >          BB is DOUBLE PRECISION array, dimension (LDBB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix B, stored in the first kb+1 rows of the array.  The */
/* >          j-th column of B is stored in the j-th column of the array BB */
/* >          as follows: */
/* >          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; */
/* >          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). */
/* > */
/* >          On exit, the factor S from the split Cholesky factorization */
/* >          B = S**T*S, as returned by DPBSTF. */
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
/* >          Q is DOUBLE PRECISION array, dimension (LDQ, N) */
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
/* >          by reducing A to tridiagonal form. */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors, with the i-th column of Z holding the */
/* >          eigenvector associated with W(i).  The eigenvectors are */
/* >          normalized so Z**T*B*Z = I. */
/* >          If JOBZ = 'N', then Z is not referenced. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (7*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* >          IFAIL is INTEGER array, dimension (M) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M elements of */
/* >          IFAIL are zero.  If INFO > 0, then IFAIL contains the */
/* >          indices of the eigenvalues that failed to converge. */
/* >          If JOBZ = 'N', then IFAIL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0 : successful exit */
/* >          < 0 : if INFO = -i, the i-th argument had an illegal value */
/* >          <= N: if INFO = i, then i eigenvectors failed to converge. */
/* >                  Their indices are stored in IFAIL. */
/* >          > N : DPBSTF returned an error code; i.e., */
/* >                if INFO = N + i, for 1 <= i <= N, then the leading */
/* >                minor of order i of B is not positive definite. */
/* >                The factorization of B could not be completed and */
/* >                no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup doubleOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int dsbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *
	bb, integer *ldbb, doublereal *q, integer *ldq, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *iwork, integer *ifail, integer *info, ftnlen jobz_len, 
	ftnlen range_len, ftnlen uplo_len)
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
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical upper, wantz, alleig, indeig;
    static integer indibl;
    static logical valeig;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpbstf_(char *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen), dsbtrd_(
	    char *, char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static integer indisp;
    extern /* Subroutine */ int dsbgst_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen),
	     dstein_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *);
    static integer indiwo;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dstebz_(char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer nsplit;


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 333 "dsbgvx.f"
    /* Parameter adjustments */
#line 333 "dsbgvx.f"
    ab_dim1 = *ldab;
#line 333 "dsbgvx.f"
    ab_offset = 1 + ab_dim1;
#line 333 "dsbgvx.f"
    ab -= ab_offset;
#line 333 "dsbgvx.f"
    bb_dim1 = *ldbb;
#line 333 "dsbgvx.f"
    bb_offset = 1 + bb_dim1;
#line 333 "dsbgvx.f"
    bb -= bb_offset;
#line 333 "dsbgvx.f"
    q_dim1 = *ldq;
#line 333 "dsbgvx.f"
    q_offset = 1 + q_dim1;
#line 333 "dsbgvx.f"
    q -= q_offset;
#line 333 "dsbgvx.f"
    --w;
#line 333 "dsbgvx.f"
    z_dim1 = *ldz;
#line 333 "dsbgvx.f"
    z_offset = 1 + z_dim1;
#line 333 "dsbgvx.f"
    z__ -= z_offset;
#line 333 "dsbgvx.f"
    --work;
#line 333 "dsbgvx.f"
    --iwork;
#line 333 "dsbgvx.f"
    --ifail;
#line 333 "dsbgvx.f"

#line 333 "dsbgvx.f"
    /* Function Body */
#line 333 "dsbgvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 334 "dsbgvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 335 "dsbgvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 336 "dsbgvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 337 "dsbgvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 339 "dsbgvx.f"
    *info = 0;
#line 340 "dsbgvx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 341 "dsbgvx.f"
	*info = -1;
#line 342 "dsbgvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 343 "dsbgvx.f"
	*info = -2;
#line 344 "dsbgvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 345 "dsbgvx.f"
	*info = -3;
#line 346 "dsbgvx.f"
    } else if (*n < 0) {
#line 347 "dsbgvx.f"
	*info = -4;
#line 348 "dsbgvx.f"
    } else if (*ka < 0) {
#line 349 "dsbgvx.f"
	*info = -5;
#line 350 "dsbgvx.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 351 "dsbgvx.f"
	*info = -6;
#line 352 "dsbgvx.f"
    } else if (*ldab < *ka + 1) {
#line 353 "dsbgvx.f"
	*info = -8;
#line 354 "dsbgvx.f"
    } else if (*ldbb < *kb + 1) {
#line 355 "dsbgvx.f"
	*info = -10;
#line 356 "dsbgvx.f"
    } else if (*ldq < 1 || wantz && *ldq < *n) {
#line 357 "dsbgvx.f"
	*info = -12;
#line 358 "dsbgvx.f"
    } else {
#line 359 "dsbgvx.f"
	if (valeig) {
#line 360 "dsbgvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 360 "dsbgvx.f"
		*info = -14;
#line 360 "dsbgvx.f"
	    }
#line 362 "dsbgvx.f"
	} else if (indeig) {
#line 363 "dsbgvx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 364 "dsbgvx.f"
		*info = -15;
#line 365 "dsbgvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 366 "dsbgvx.f"
		*info = -16;
#line 367 "dsbgvx.f"
	    }
#line 368 "dsbgvx.f"
	}
#line 369 "dsbgvx.f"
    }
#line 370 "dsbgvx.f"
    if (*info == 0) {
#line 371 "dsbgvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 372 "dsbgvx.f"
	    *info = -21;
#line 373 "dsbgvx.f"
	}
#line 374 "dsbgvx.f"
    }

#line 376 "dsbgvx.f"
    if (*info != 0) {
#line 377 "dsbgvx.f"
	i__1 = -(*info);
#line 377 "dsbgvx.f"
	xerbla_("DSBGVX", &i__1, (ftnlen)6);
#line 378 "dsbgvx.f"
	return 0;
#line 379 "dsbgvx.f"
    }

/*     Quick return if possible */

#line 383 "dsbgvx.f"
    *m = 0;
#line 384 "dsbgvx.f"
    if (*n == 0) {
#line 384 "dsbgvx.f"
	return 0;
#line 384 "dsbgvx.f"
    }

/*     Form a split Cholesky factorization of B. */

#line 389 "dsbgvx.f"
    dpbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info, (ftnlen)1);
#line 390 "dsbgvx.f"
    if (*info != 0) {
#line 391 "dsbgvx.f"
	*info = *n + *info;
#line 392 "dsbgvx.f"
	return 0;
#line 393 "dsbgvx.f"
    }

/*     Transform problem to standard eigenvalue problem. */

#line 397 "dsbgvx.f"
    dsbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
	     &q[q_offset], ldq, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);

/*     Reduce symmetric band matrix to tridiagonal form. */

#line 402 "dsbgvx.f"
    indd = 1;
#line 403 "dsbgvx.f"
    inde = indd + *n;
#line 404 "dsbgvx.f"
    indwrk = inde + *n;
#line 405 "dsbgvx.f"
    if (wantz) {
#line 406 "dsbgvx.f"
	*(unsigned char *)vect = 'U';
#line 407 "dsbgvx.f"
    } else {
#line 408 "dsbgvx.f"
	*(unsigned char *)vect = 'N';
#line 409 "dsbgvx.f"
    }
#line 410 "dsbgvx.f"
    dsbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &work[indd], &work[inde],
	     &q[q_offset], ldq, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call DSTERF or SSTEQR.  If this fails for some */
/*     eigenvalue, then try DSTEBZ. */

#line 417 "dsbgvx.f"
    test = FALSE_;
#line 418 "dsbgvx.f"
    if (indeig) {
#line 419 "dsbgvx.f"
	if (*il == 1 && *iu == *n) {
#line 420 "dsbgvx.f"
	    test = TRUE_;
#line 421 "dsbgvx.f"
	}
#line 422 "dsbgvx.f"
    }
#line 423 "dsbgvx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 424 "dsbgvx.f"
	dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 425 "dsbgvx.f"
	indee = indwrk + (*n << 1);
#line 426 "dsbgvx.f"
	i__1 = *n - 1;
#line 426 "dsbgvx.f"
	dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 427 "dsbgvx.f"
	if (! wantz) {
#line 428 "dsbgvx.f"
	    dsterf_(n, &w[1], &work[indee], info);
#line 429 "dsbgvx.f"
	} else {
#line 430 "dsbgvx.f"
	    dlacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 431 "dsbgvx.f"
	    dsteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 433 "dsbgvx.f"
	    if (*info == 0) {
#line 434 "dsbgvx.f"
		i__1 = *n;
#line 434 "dsbgvx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 435 "dsbgvx.f"
		    ifail[i__] = 0;
#line 436 "dsbgvx.f"
/* L10: */
#line 436 "dsbgvx.f"
		}
#line 437 "dsbgvx.f"
	    }
#line 438 "dsbgvx.f"
	}
#line 439 "dsbgvx.f"
	if (*info == 0) {
#line 440 "dsbgvx.f"
	    *m = *n;
#line 441 "dsbgvx.f"
	    goto L30;
#line 442 "dsbgvx.f"
	}
#line 443 "dsbgvx.f"
	*info = 0;
#line 444 "dsbgvx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, */
/*     call DSTEIN. */

#line 449 "dsbgvx.f"
    if (wantz) {
#line 450 "dsbgvx.f"
	*(unsigned char *)order = 'B';
#line 451 "dsbgvx.f"
    } else {
#line 452 "dsbgvx.f"
	*(unsigned char *)order = 'E';
#line 453 "dsbgvx.f"
    }
#line 454 "dsbgvx.f"
    indibl = 1;
#line 455 "dsbgvx.f"
    indisp = indibl + *n;
#line 456 "dsbgvx.f"
    indiwo = indisp + *n;
#line 457 "dsbgvx.f"
    dstebz_(range, order, n, vl, vu, il, iu, abstol, &work[indd], &work[inde],
	     m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[indwrk],
	     &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 462 "dsbgvx.f"
    if (wantz) {
#line 463 "dsbgvx.f"
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply transformation matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

#line 470 "dsbgvx.f"
	i__1 = *m;
#line 470 "dsbgvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 471 "dsbgvx.f"
	    dcopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 472 "dsbgvx.f"
	    dgemv_("N", n, n, &c_b25, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b27, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 474 "dsbgvx.f"
/* L20: */
#line 474 "dsbgvx.f"
	}
#line 475 "dsbgvx.f"
    }

#line 477 "dsbgvx.f"
L30:

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 482 "dsbgvx.f"
    if (wantz) {
#line 483 "dsbgvx.f"
	i__1 = *m - 1;
#line 483 "dsbgvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 484 "dsbgvx.f"
	    i__ = 0;
#line 485 "dsbgvx.f"
	    tmp1 = w[j];
#line 486 "dsbgvx.f"
	    i__2 = *m;
#line 486 "dsbgvx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 487 "dsbgvx.f"
		if (w[jj] < tmp1) {
#line 488 "dsbgvx.f"
		    i__ = jj;
#line 489 "dsbgvx.f"
		    tmp1 = w[jj];
#line 490 "dsbgvx.f"
		}
#line 491 "dsbgvx.f"
/* L40: */
#line 491 "dsbgvx.f"
	    }

#line 493 "dsbgvx.f"
	    if (i__ != 0) {
#line 494 "dsbgvx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 495 "dsbgvx.f"
		w[i__] = w[j];
#line 496 "dsbgvx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 497 "dsbgvx.f"
		w[j] = tmp1;
#line 498 "dsbgvx.f"
		iwork[indibl + j - 1] = itmp1;
#line 499 "dsbgvx.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 500 "dsbgvx.f"
		if (*info != 0) {
#line 501 "dsbgvx.f"
		    itmp1 = ifail[i__];
#line 502 "dsbgvx.f"
		    ifail[i__] = ifail[j];
#line 503 "dsbgvx.f"
		    ifail[j] = itmp1;
#line 504 "dsbgvx.f"
		}
#line 505 "dsbgvx.f"
	    }
#line 506 "dsbgvx.f"
/* L50: */
#line 506 "dsbgvx.f"
	}
#line 507 "dsbgvx.f"
    }

#line 509 "dsbgvx.f"
    return 0;

/*     End of DSBGVX */

} /* dsbgvx_ */

