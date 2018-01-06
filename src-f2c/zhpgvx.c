#line 1 "zhpgvx.f"
/* zhpgvx.f -- translated by f2c (version 20100827).
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

#line 1 "zhpgvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZHPGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpgvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpgvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpgvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, */
/*                          IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPGVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and */
/* > B are assumed to be Hermitian, stored in packed format, and B is also */
/* > positive definite.  Eigenvalues and eigenvectors can be selected by */
/* > specifying either a range of values or a range of indices for the */
/* > desired eigenvalues. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the contents of AP are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BP */
/* > \verbatim */
/* >          BP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          B, packed columnwise in a linear array.  The j-th column of B */
/* >          is stored in the array BP as follows: */
/* >          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the triangular factor U or L from the Cholesky */
/* >          factorization B = U**H*U or B = L*L**H, in the same storage */
/* >          format as B. */
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
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          The eigenvectors are normalized as follows: */
/* >          if ITYPE = 1 or 2, Z**H*B*Z = I; */
/* >          if ITYPE = 3, Z**H*inv(B)*Z = I. */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
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
/* >          > 0:  ZPPTRF or ZHPEVX returned an error code: */
/* >             <= N:  if INFO = i, ZHPEVX failed to converge; */
/* >                    i eigenvectors failed to converge.  Their indices */
/* >                    are stored in array IFAIL. */
/* >             > N:   if INFO = N + i, for 1 <= i <= n, then the leading */
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

/* > \ingroup complex16OTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int zhpgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *
	vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info, ftnlen jobz_len, ftnlen range_len, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;

    /* Local variables */
    static integer j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char trans[1];
    static logical upper, wantz;
    extern /* Subroutine */ int ztpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen), ztpsv_(char *, char *, char *, integer *, doublecomplex *
	    , doublecomplex *, integer *, ftnlen, ftnlen, ftnlen);
    static logical alleig, indeig, valeig;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zhpgst_(
	    integer *, char *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen), zhpevx_(char *, char *, char *, integer *, 
	    doublecomplex *, doublereal *, doublereal *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublecomplex *, integer *
	    , doublecomplex *, doublereal *, integer *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), zpptrf_(char *, integer *, doublecomplex 
	    *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 308 "zhpgvx.f"
    /* Parameter adjustments */
#line 308 "zhpgvx.f"
    --ap;
#line 308 "zhpgvx.f"
    --bp;
#line 308 "zhpgvx.f"
    --w;
#line 308 "zhpgvx.f"
    z_dim1 = *ldz;
#line 308 "zhpgvx.f"
    z_offset = 1 + z_dim1;
#line 308 "zhpgvx.f"
    z__ -= z_offset;
#line 308 "zhpgvx.f"
    --work;
#line 308 "zhpgvx.f"
    --rwork;
#line 308 "zhpgvx.f"
    --iwork;
#line 308 "zhpgvx.f"
    --ifail;
#line 308 "zhpgvx.f"

#line 308 "zhpgvx.f"
    /* Function Body */
#line 308 "zhpgvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 309 "zhpgvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 310 "zhpgvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 311 "zhpgvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 312 "zhpgvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 314 "zhpgvx.f"
    *info = 0;
#line 315 "zhpgvx.f"
    if (*itype < 1 || *itype > 3) {
#line 316 "zhpgvx.f"
	*info = -1;
#line 317 "zhpgvx.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 318 "zhpgvx.f"
	*info = -2;
#line 319 "zhpgvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 320 "zhpgvx.f"
	*info = -3;
#line 321 "zhpgvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 322 "zhpgvx.f"
	*info = -4;
#line 323 "zhpgvx.f"
    } else if (*n < 0) {
#line 324 "zhpgvx.f"
	*info = -5;
#line 325 "zhpgvx.f"
    } else {
#line 326 "zhpgvx.f"
	if (valeig) {
#line 327 "zhpgvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 328 "zhpgvx.f"
		*info = -9;
#line 329 "zhpgvx.f"
	    }
#line 330 "zhpgvx.f"
	} else if (indeig) {
#line 331 "zhpgvx.f"
	    if (*il < 1) {
#line 332 "zhpgvx.f"
		*info = -10;
#line 333 "zhpgvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 334 "zhpgvx.f"
		*info = -11;
#line 335 "zhpgvx.f"
	    }
#line 336 "zhpgvx.f"
	}
#line 337 "zhpgvx.f"
    }
#line 338 "zhpgvx.f"
    if (*info == 0) {
#line 339 "zhpgvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 340 "zhpgvx.f"
	    *info = -16;
#line 341 "zhpgvx.f"
	}
#line 342 "zhpgvx.f"
    }

#line 344 "zhpgvx.f"
    if (*info != 0) {
#line 345 "zhpgvx.f"
	i__1 = -(*info);
#line 345 "zhpgvx.f"
	xerbla_("ZHPGVX", &i__1, (ftnlen)6);
#line 346 "zhpgvx.f"
	return 0;
#line 347 "zhpgvx.f"
    }

/*     Quick return if possible */

#line 351 "zhpgvx.f"
    if (*n == 0) {
#line 351 "zhpgvx.f"
	return 0;
#line 351 "zhpgvx.f"
    }

/*     Form a Cholesky factorization of B. */

#line 356 "zhpgvx.f"
    zpptrf_(uplo, n, &bp[1], info, (ftnlen)1);
#line 357 "zhpgvx.f"
    if (*info != 0) {
#line 358 "zhpgvx.f"
	*info = *n + *info;
#line 359 "zhpgvx.f"
	return 0;
#line 360 "zhpgvx.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 364 "zhpgvx.f"
    zhpgst_(itype, uplo, n, &ap[1], &bp[1], info, (ftnlen)1);
#line 365 "zhpgvx.f"
    zhpevx_(jobz, range, uplo, n, &ap[1], vl, vu, il, iu, abstol, m, &w[1], &
	    z__[z_offset], ldz, &work[1], &rwork[1], &iwork[1], &ifail[1], 
	    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 368 "zhpgvx.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 372 "zhpgvx.f"
	if (*info > 0) {
#line 372 "zhpgvx.f"
	    *m = *info - 1;
#line 372 "zhpgvx.f"
	}
#line 374 "zhpgvx.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y */

#line 379 "zhpgvx.f"
	    if (upper) {
#line 380 "zhpgvx.f"
		*(unsigned char *)trans = 'N';
#line 381 "zhpgvx.f"
	    } else {
#line 382 "zhpgvx.f"
		*(unsigned char *)trans = 'C';
#line 383 "zhpgvx.f"
	    }

#line 385 "zhpgvx.f"
	    i__1 = *m;
#line 385 "zhpgvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 386 "zhpgvx.f"
		ztpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 388 "zhpgvx.f"
/* L10: */
#line 388 "zhpgvx.f"
	    }

#line 390 "zhpgvx.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H *y */

#line 395 "zhpgvx.f"
	    if (upper) {
#line 396 "zhpgvx.f"
		*(unsigned char *)trans = 'C';
#line 397 "zhpgvx.f"
	    } else {
#line 398 "zhpgvx.f"
		*(unsigned char *)trans = 'N';
#line 399 "zhpgvx.f"
	    }

#line 401 "zhpgvx.f"
	    i__1 = *m;
#line 401 "zhpgvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 402 "zhpgvx.f"
		ztpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 404 "zhpgvx.f"
/* L20: */
#line 404 "zhpgvx.f"
	    }
#line 405 "zhpgvx.f"
	}
#line 406 "zhpgvx.f"
    }

#line 408 "zhpgvx.f"
    return 0;

/*     End of ZHPGVX */

} /* zhpgvx_ */

