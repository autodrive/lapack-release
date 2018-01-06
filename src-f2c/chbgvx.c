#line 1 "chbgvx.f"
/* chbgvx.f -- translated by f2c (version 20100827).
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

#line 1 "chbgvx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CHBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbgvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbgvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbgvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, */
/*                          LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, */
/*                          LDZ, WORK, RWORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M, */
/*      $                   N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBGVX computes all the eigenvalues, and optionally, the eigenvectors */
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
/* >          AB is COMPLEX array, dimension (LDAB, N) */
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
/* >          BB is COMPLEX array, dimension (LDBB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix B, stored in the first kb+1 rows of the array.  The */
/* >          j-th column of B is stored in the j-th column of the array BB */
/* >          as follows: */
/* >          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; */
/* >          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). */
/* > */
/* >          On exit, the factor S from the split Cholesky factorization */
/* >          B = S**H*S, as returned by CPBSTF. */
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
/* >          Q is COMPLEX array, dimension (LDQ, N) */
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
/* >          VL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
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
/* >          by reducing AP to tridiagonal form. */
/* > */
/* >          Eigenvalues will be computed most accurately when ABSTOL is */
/* >          set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
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
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
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
/* >          WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (7*N) */
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
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then CPBSTF */
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

/* > \date November 2011 */

/* > \ingroup complexOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int chbgvx_(char *jobz, char *range, char *uplo, integer *n, 
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
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantz, alleig, indeig;
    static integer indibl;
    extern /* Subroutine */ int chbtrd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    static logical valeig;
    extern /* Subroutine */ int chbgst_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, ftnlen, ftnlen), clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), cpbstf_(char *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, ftnlen);
    static integer indiwk, indisp;
    extern /* Subroutine */ int cstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), ssterf_(integer *, doublereal *, doublereal *, integer *
	    );
    static integer nsplit;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 344 "chbgvx.f"
    /* Parameter adjustments */
#line 344 "chbgvx.f"
    ab_dim1 = *ldab;
#line 344 "chbgvx.f"
    ab_offset = 1 + ab_dim1;
#line 344 "chbgvx.f"
    ab -= ab_offset;
#line 344 "chbgvx.f"
    bb_dim1 = *ldbb;
#line 344 "chbgvx.f"
    bb_offset = 1 + bb_dim1;
#line 344 "chbgvx.f"
    bb -= bb_offset;
#line 344 "chbgvx.f"
    q_dim1 = *ldq;
#line 344 "chbgvx.f"
    q_offset = 1 + q_dim1;
#line 344 "chbgvx.f"
    q -= q_offset;
#line 344 "chbgvx.f"
    --w;
#line 344 "chbgvx.f"
    z_dim1 = *ldz;
#line 344 "chbgvx.f"
    z_offset = 1 + z_dim1;
#line 344 "chbgvx.f"
    z__ -= z_offset;
#line 344 "chbgvx.f"
    --work;
#line 344 "chbgvx.f"
    --rwork;
#line 344 "chbgvx.f"
    --iwork;
#line 344 "chbgvx.f"
    --ifail;
#line 344 "chbgvx.f"

#line 344 "chbgvx.f"
    /* Function Body */
#line 344 "chbgvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 345 "chbgvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 346 "chbgvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 347 "chbgvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 348 "chbgvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 350 "chbgvx.f"
    *info = 0;
#line 351 "chbgvx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 352 "chbgvx.f"
	*info = -1;
#line 353 "chbgvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 354 "chbgvx.f"
	*info = -2;
#line 355 "chbgvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 356 "chbgvx.f"
	*info = -3;
#line 357 "chbgvx.f"
    } else if (*n < 0) {
#line 358 "chbgvx.f"
	*info = -4;
#line 359 "chbgvx.f"
    } else if (*ka < 0) {
#line 360 "chbgvx.f"
	*info = -5;
#line 361 "chbgvx.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 362 "chbgvx.f"
	*info = -6;
#line 363 "chbgvx.f"
    } else if (*ldab < *ka + 1) {
#line 364 "chbgvx.f"
	*info = -8;
#line 365 "chbgvx.f"
    } else if (*ldbb < *kb + 1) {
#line 366 "chbgvx.f"
	*info = -10;
#line 367 "chbgvx.f"
    } else if (*ldq < 1 || wantz && *ldq < *n) {
#line 368 "chbgvx.f"
	*info = -12;
#line 369 "chbgvx.f"
    } else {
#line 370 "chbgvx.f"
	if (valeig) {
#line 371 "chbgvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 371 "chbgvx.f"
		*info = -14;
#line 371 "chbgvx.f"
	    }
#line 373 "chbgvx.f"
	} else if (indeig) {
#line 374 "chbgvx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 375 "chbgvx.f"
		*info = -15;
#line 376 "chbgvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 377 "chbgvx.f"
		*info = -16;
#line 378 "chbgvx.f"
	    }
#line 379 "chbgvx.f"
	}
#line 380 "chbgvx.f"
    }
#line 381 "chbgvx.f"
    if (*info == 0) {
#line 382 "chbgvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 383 "chbgvx.f"
	    *info = -21;
#line 384 "chbgvx.f"
	}
#line 385 "chbgvx.f"
    }

#line 387 "chbgvx.f"
    if (*info != 0) {
#line 388 "chbgvx.f"
	i__1 = -(*info);
#line 388 "chbgvx.f"
	xerbla_("CHBGVX", &i__1, (ftnlen)6);
#line 389 "chbgvx.f"
	return 0;
#line 390 "chbgvx.f"
    }

/*     Quick return if possible */

#line 394 "chbgvx.f"
    *m = 0;
#line 395 "chbgvx.f"
    if (*n == 0) {
#line 395 "chbgvx.f"
	return 0;
#line 395 "chbgvx.f"
    }

/*     Form a split Cholesky factorization of B. */

#line 400 "chbgvx.f"
    cpbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info, (ftnlen)1);
#line 401 "chbgvx.f"
    if (*info != 0) {
#line 402 "chbgvx.f"
	*info = *n + *info;
#line 403 "chbgvx.f"
	return 0;
#line 404 "chbgvx.f"
    }

/*     Transform problem to standard eigenvalue problem. */

#line 408 "chbgvx.f"
    chbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
	     &q[q_offset], ldq, &work[1], &rwork[1], &iinfo, (ftnlen)1, (
	    ftnlen)1);

/*     Solve the standard eigenvalue problem. */
/*     Reduce Hermitian band matrix to tridiagonal form. */

#line 414 "chbgvx.f"
    indd = 1;
#line 415 "chbgvx.f"
    inde = indd + *n;
#line 416 "chbgvx.f"
    indrwk = inde + *n;
#line 417 "chbgvx.f"
    indwrk = 1;
#line 418 "chbgvx.f"
    if (wantz) {
#line 419 "chbgvx.f"
	*(unsigned char *)vect = 'U';
#line 420 "chbgvx.f"
    } else {
#line 421 "chbgvx.f"
	*(unsigned char *)vect = 'N';
#line 422 "chbgvx.f"
    }
#line 423 "chbgvx.f"
    chbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &rwork[indd], &rwork[
	    inde], &q[q_offset], ldq, &work[indwrk], &iinfo, (ftnlen)1, (
	    ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call SSTERF or CSTEQR.  If this fails for some */
/*     eigenvalue, then try SSTEBZ. */

#line 430 "chbgvx.f"
    test = FALSE_;
#line 431 "chbgvx.f"
    if (indeig) {
#line 432 "chbgvx.f"
	if (*il == 1 && *iu == *n) {
#line 433 "chbgvx.f"
	    test = TRUE_;
#line 434 "chbgvx.f"
	}
#line 435 "chbgvx.f"
    }
#line 436 "chbgvx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 437 "chbgvx.f"
	scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 438 "chbgvx.f"
	indee = indrwk + (*n << 1);
#line 439 "chbgvx.f"
	i__1 = *n - 1;
#line 439 "chbgvx.f"
	scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 440 "chbgvx.f"
	if (! wantz) {
#line 441 "chbgvx.f"
	    ssterf_(n, &w[1], &rwork[indee], info);
#line 442 "chbgvx.f"
	} else {
#line 443 "chbgvx.f"
	    clacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 444 "chbgvx.f"
	    csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 446 "chbgvx.f"
	    if (*info == 0) {
#line 447 "chbgvx.f"
		i__1 = *n;
#line 447 "chbgvx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 448 "chbgvx.f"
		    ifail[i__] = 0;
#line 449 "chbgvx.f"
/* L10: */
#line 449 "chbgvx.f"
		}
#line 450 "chbgvx.f"
	    }
#line 451 "chbgvx.f"
	}
#line 452 "chbgvx.f"
	if (*info == 0) {
#line 453 "chbgvx.f"
	    *m = *n;
#line 454 "chbgvx.f"
	    goto L30;
#line 455 "chbgvx.f"
	}
#line 456 "chbgvx.f"
	*info = 0;
#line 457 "chbgvx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, */
/*     call CSTEIN. */

#line 462 "chbgvx.f"
    if (wantz) {
#line 463 "chbgvx.f"
	*(unsigned char *)order = 'B';
#line 464 "chbgvx.f"
    } else {
#line 465 "chbgvx.f"
	*(unsigned char *)order = 'E';
#line 466 "chbgvx.f"
    }
#line 467 "chbgvx.f"
    indibl = 1;
#line 468 "chbgvx.f"
    indisp = indibl + *n;
#line 469 "chbgvx.f"
    indiwk = indisp + *n;
#line 470 "chbgvx.f"
    sstebz_(range, order, n, vl, vu, il, iu, abstol, &rwork[indd], &rwork[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &rwork[
	    indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 475 "chbgvx.f"
    if (wantz) {
#line 476 "chbgvx.f"
	cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by CSTEIN. */

#line 483 "chbgvx.f"
	i__1 = *m;
#line 483 "chbgvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 484 "chbgvx.f"
	    ccopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 485 "chbgvx.f"
	    cgemv_("N", n, n, &c_b2, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b1, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 487 "chbgvx.f"
/* L20: */
#line 487 "chbgvx.f"
	}
#line 488 "chbgvx.f"
    }

#line 490 "chbgvx.f"
L30:

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 495 "chbgvx.f"
    if (wantz) {
#line 496 "chbgvx.f"
	i__1 = *m - 1;
#line 496 "chbgvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 497 "chbgvx.f"
	    i__ = 0;
#line 498 "chbgvx.f"
	    tmp1 = w[j];
#line 499 "chbgvx.f"
	    i__2 = *m;
#line 499 "chbgvx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 500 "chbgvx.f"
		if (w[jj] < tmp1) {
#line 501 "chbgvx.f"
		    i__ = jj;
#line 502 "chbgvx.f"
		    tmp1 = w[jj];
#line 503 "chbgvx.f"
		}
#line 504 "chbgvx.f"
/* L40: */
#line 504 "chbgvx.f"
	    }

#line 506 "chbgvx.f"
	    if (i__ != 0) {
#line 507 "chbgvx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 508 "chbgvx.f"
		w[i__] = w[j];
#line 509 "chbgvx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 510 "chbgvx.f"
		w[j] = tmp1;
#line 511 "chbgvx.f"
		iwork[indibl + j - 1] = itmp1;
#line 512 "chbgvx.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 513 "chbgvx.f"
		if (*info != 0) {
#line 514 "chbgvx.f"
		    itmp1 = ifail[i__];
#line 515 "chbgvx.f"
		    ifail[i__] = ifail[j];
#line 516 "chbgvx.f"
		    ifail[j] = itmp1;
#line 517 "chbgvx.f"
		}
#line 518 "chbgvx.f"
	    }
#line 519 "chbgvx.f"
/* L50: */
#line 519 "chbgvx.f"
	}
#line 520 "chbgvx.f"
    }

#line 522 "chbgvx.f"
    return 0;

/*     End of CHBGVX */

} /* chbgvx_ */

