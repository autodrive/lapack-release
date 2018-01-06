#line 1 "chegvx.f"
/* chegvx.f -- translated by f2c (version 20100827).
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

#line 1 "chegvx.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CHEGVX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, */
/*                          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
/*                          LWORK, RWORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEGVX computes selected eigenvalues, and optionally, eigenvectors */
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
/* >          A is COMPLEX array, dimension (LDA, N) */
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
/* >          B is COMPLEX array, dimension (LDB, N) */
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
/* >          VL is REAL */
/* > */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
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
/* >          The first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= max(1,2*N). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the blocksize for CHETRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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
/* >          > 0:  CPOTRF or CHEEVX returned an error code: */
/* >             <= N:  if INFO = i, CHEEVX failed to converge; */
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

/* > \ingroup complexHEeigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int chegvx_(integer *itype, char *jobz, char *range, char *
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
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static char trans[1];
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper, wantz, alleig, indeig, valeig;
    extern /* Subroutine */ int chegst_(integer *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cheevx_(
	    char *, char *, char *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), cpotrf_(char *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;


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

#line 353 "chegvx.f"
    /* Parameter adjustments */
#line 353 "chegvx.f"
    a_dim1 = *lda;
#line 353 "chegvx.f"
    a_offset = 1 + a_dim1;
#line 353 "chegvx.f"
    a -= a_offset;
#line 353 "chegvx.f"
    b_dim1 = *ldb;
#line 353 "chegvx.f"
    b_offset = 1 + b_dim1;
#line 353 "chegvx.f"
    b -= b_offset;
#line 353 "chegvx.f"
    --w;
#line 353 "chegvx.f"
    z_dim1 = *ldz;
#line 353 "chegvx.f"
    z_offset = 1 + z_dim1;
#line 353 "chegvx.f"
    z__ -= z_offset;
#line 353 "chegvx.f"
    --work;
#line 353 "chegvx.f"
    --rwork;
#line 353 "chegvx.f"
    --iwork;
#line 353 "chegvx.f"
    --ifail;
#line 353 "chegvx.f"

#line 353 "chegvx.f"
    /* Function Body */
#line 353 "chegvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 354 "chegvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 355 "chegvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 356 "chegvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 357 "chegvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 358 "chegvx.f"
    lquery = *lwork == -1;

#line 360 "chegvx.f"
    *info = 0;
#line 361 "chegvx.f"
    if (*itype < 1 || *itype > 3) {
#line 362 "chegvx.f"
	*info = -1;
#line 363 "chegvx.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 364 "chegvx.f"
	*info = -2;
#line 365 "chegvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 366 "chegvx.f"
	*info = -3;
#line 367 "chegvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 368 "chegvx.f"
	*info = -4;
#line 369 "chegvx.f"
    } else if (*n < 0) {
#line 370 "chegvx.f"
	*info = -5;
#line 371 "chegvx.f"
    } else if (*lda < max(1,*n)) {
#line 372 "chegvx.f"
	*info = -7;
#line 373 "chegvx.f"
    } else if (*ldb < max(1,*n)) {
#line 374 "chegvx.f"
	*info = -9;
#line 375 "chegvx.f"
    } else {
#line 376 "chegvx.f"
	if (valeig) {
#line 377 "chegvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 377 "chegvx.f"
		*info = -11;
#line 377 "chegvx.f"
	    }
#line 379 "chegvx.f"
	} else if (indeig) {
#line 380 "chegvx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 381 "chegvx.f"
		*info = -12;
#line 382 "chegvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 383 "chegvx.f"
		*info = -13;
#line 384 "chegvx.f"
	    }
#line 385 "chegvx.f"
	}
#line 386 "chegvx.f"
    }
#line 387 "chegvx.f"
    if (*info == 0) {
#line 388 "chegvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 389 "chegvx.f"
	    *info = -18;
#line 390 "chegvx.f"
	}
#line 391 "chegvx.f"
    }

#line 393 "chegvx.f"
    if (*info == 0) {
#line 394 "chegvx.f"
	nb = ilaenv_(&c__1, "CHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 395 "chegvx.f"
	i__1 = 1, i__2 = (nb + 1) * *n;
#line 395 "chegvx.f"
	lwkopt = max(i__1,i__2);
#line 396 "chegvx.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

/* Computing MAX */
#line 398 "chegvx.f"
	i__1 = 1, i__2 = *n << 1;
#line 398 "chegvx.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 399 "chegvx.f"
	    *info = -20;
#line 400 "chegvx.f"
	}
#line 401 "chegvx.f"
    }

#line 403 "chegvx.f"
    if (*info != 0) {
#line 404 "chegvx.f"
	i__1 = -(*info);
#line 404 "chegvx.f"
	xerbla_("CHEGVX", &i__1, (ftnlen)6);
#line 405 "chegvx.f"
	return 0;
#line 406 "chegvx.f"
    } else if (lquery) {
#line 407 "chegvx.f"
	return 0;
#line 408 "chegvx.f"
    }

/*     Quick return if possible */

#line 412 "chegvx.f"
    *m = 0;
#line 413 "chegvx.f"
    if (*n == 0) {
#line 414 "chegvx.f"
	return 0;
#line 415 "chegvx.f"
    }

/*     Form a Cholesky factorization of B. */

#line 419 "chegvx.f"
    cpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 420 "chegvx.f"
    if (*info != 0) {
#line 421 "chegvx.f"
	*info = *n + *info;
#line 422 "chegvx.f"
	return 0;
#line 423 "chegvx.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 427 "chegvx.f"
    chegst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 428 "chegvx.f"
    cheevx_(jobz, range, uplo, n, &a[a_offset], lda, vl, vu, il, iu, abstol, 
	    m, &w[1], &z__[z_offset], ldz, &work[1], lwork, &rwork[1], &iwork[
	    1], &ifail[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 432 "chegvx.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 436 "chegvx.f"
	if (*info > 0) {
#line 436 "chegvx.f"
	    *m = *info - 1;
#line 436 "chegvx.f"
	}
#line 438 "chegvx.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y */

#line 443 "chegvx.f"
	    if (upper) {
#line 444 "chegvx.f"
		*(unsigned char *)trans = 'N';
#line 445 "chegvx.f"
	    } else {
#line 446 "chegvx.f"
		*(unsigned char *)trans = 'C';
#line 447 "chegvx.f"
	    }

#line 449 "chegvx.f"
	    ctrsm_("Left", uplo, trans, "Non-unit", n, m, &c_b1, &b[b_offset],
		     ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);

#line 452 "chegvx.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H*y */

#line 457 "chegvx.f"
	    if (upper) {
#line 458 "chegvx.f"
		*(unsigned char *)trans = 'C';
#line 459 "chegvx.f"
	    } else {
#line 460 "chegvx.f"
		*(unsigned char *)trans = 'N';
#line 461 "chegvx.f"
	    }

#line 463 "chegvx.f"
	    ctrmm_("Left", uplo, trans, "Non-unit", n, m, &c_b1, &b[b_offset],
		     ldb, &z__[z_offset], ldz, (ftnlen)4, (ftnlen)1, (ftnlen)
		    1, (ftnlen)8);
#line 465 "chegvx.f"
	}
#line 466 "chegvx.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 470 "chegvx.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 472 "chegvx.f"
    return 0;

/*     End of CHEGVX */

} /* chegvx_ */

