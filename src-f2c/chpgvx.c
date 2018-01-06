#line 1 "chpgvx.f"
/* chpgvx.f -- translated by f2c (version 20100827).
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

#line 1 "chpgvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHPGVX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpgvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpgvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpgvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, */
/*                          IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AP( * ), BP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPGVX computes selected eigenvalues and, optionally, eigenvectors */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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
/* >          BP is COMPLEX array, dimension (N*(N+1)/2) */
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
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
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
/* >          WORK is COMPLEX array, dimension (2*N) */
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
/* >          > 0:  CPPTRF or CHPEVX returned an error code: */
/* >             <= N:  if INFO = i, CHPEVX failed to converge; */
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

/* > \date June 2016 */

/* > \ingroup complexOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int chpgvx_(integer *itype, char *jobz, char *range, char *
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
    extern /* Subroutine */ int ctpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical upper;
    extern /* Subroutine */ int ctpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical wantz, alleig, indeig, valeig;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), chpgst_(
	    integer *, char *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen), chpevx_(char *, char *, char *, integer *, 
	    doublecomplex *, doublereal *, doublereal *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublecomplex *, integer *
	    , doublecomplex *, doublereal *, integer *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), cpptrf_(char *, integer *, doublecomplex 
	    *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 317 "chpgvx.f"
    /* Parameter adjustments */
#line 317 "chpgvx.f"
    --ap;
#line 317 "chpgvx.f"
    --bp;
#line 317 "chpgvx.f"
    --w;
#line 317 "chpgvx.f"
    z_dim1 = *ldz;
#line 317 "chpgvx.f"
    z_offset = 1 + z_dim1;
#line 317 "chpgvx.f"
    z__ -= z_offset;
#line 317 "chpgvx.f"
    --work;
#line 317 "chpgvx.f"
    --rwork;
#line 317 "chpgvx.f"
    --iwork;
#line 317 "chpgvx.f"
    --ifail;
#line 317 "chpgvx.f"

#line 317 "chpgvx.f"
    /* Function Body */
#line 317 "chpgvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 318 "chpgvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 319 "chpgvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 320 "chpgvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 321 "chpgvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 323 "chpgvx.f"
    *info = 0;
#line 324 "chpgvx.f"
    if (*itype < 1 || *itype > 3) {
#line 325 "chpgvx.f"
	*info = -1;
#line 326 "chpgvx.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 327 "chpgvx.f"
	*info = -2;
#line 328 "chpgvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 329 "chpgvx.f"
	*info = -3;
#line 330 "chpgvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 331 "chpgvx.f"
	*info = -4;
#line 332 "chpgvx.f"
    } else if (*n < 0) {
#line 333 "chpgvx.f"
	*info = -5;
#line 334 "chpgvx.f"
    } else {
#line 335 "chpgvx.f"
	if (valeig) {
#line 336 "chpgvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 337 "chpgvx.f"
		*info = -9;
#line 338 "chpgvx.f"
	    }
#line 339 "chpgvx.f"
	} else if (indeig) {
#line 340 "chpgvx.f"
	    if (*il < 1) {
#line 341 "chpgvx.f"
		*info = -10;
#line 342 "chpgvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 343 "chpgvx.f"
		*info = -11;
#line 344 "chpgvx.f"
	    }
#line 345 "chpgvx.f"
	}
#line 346 "chpgvx.f"
    }
#line 347 "chpgvx.f"
    if (*info == 0) {
#line 348 "chpgvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 349 "chpgvx.f"
	    *info = -16;
#line 350 "chpgvx.f"
	}
#line 351 "chpgvx.f"
    }

#line 353 "chpgvx.f"
    if (*info != 0) {
#line 354 "chpgvx.f"
	i__1 = -(*info);
#line 354 "chpgvx.f"
	xerbla_("CHPGVX", &i__1, (ftnlen)6);
#line 355 "chpgvx.f"
	return 0;
#line 356 "chpgvx.f"
    }

/*     Quick return if possible */

#line 360 "chpgvx.f"
    if (*n == 0) {
#line 360 "chpgvx.f"
	return 0;
#line 360 "chpgvx.f"
    }

/*     Form a Cholesky factorization of B. */

#line 365 "chpgvx.f"
    cpptrf_(uplo, n, &bp[1], info, (ftnlen)1);
#line 366 "chpgvx.f"
    if (*info != 0) {
#line 367 "chpgvx.f"
	*info = *n + *info;
#line 368 "chpgvx.f"
	return 0;
#line 369 "chpgvx.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 373 "chpgvx.f"
    chpgst_(itype, uplo, n, &ap[1], &bp[1], info, (ftnlen)1);
#line 374 "chpgvx.f"
    chpevx_(jobz, range, uplo, n, &ap[1], vl, vu, il, iu, abstol, m, &w[1], &
	    z__[z_offset], ldz, &work[1], &rwork[1], &iwork[1], &ifail[1], 
	    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 377 "chpgvx.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 381 "chpgvx.f"
	if (*info > 0) {
#line 381 "chpgvx.f"
	    *m = *info - 1;
#line 381 "chpgvx.f"
	}
#line 383 "chpgvx.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y */

#line 388 "chpgvx.f"
	    if (upper) {
#line 389 "chpgvx.f"
		*(unsigned char *)trans = 'N';
#line 390 "chpgvx.f"
	    } else {
#line 391 "chpgvx.f"
		*(unsigned char *)trans = 'C';
#line 392 "chpgvx.f"
	    }

#line 394 "chpgvx.f"
	    i__1 = *m;
#line 394 "chpgvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 395 "chpgvx.f"
		ctpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 397 "chpgvx.f"
/* L10: */
#line 397 "chpgvx.f"
	    }

#line 399 "chpgvx.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H*y */

#line 404 "chpgvx.f"
	    if (upper) {
#line 405 "chpgvx.f"
		*(unsigned char *)trans = 'C';
#line 406 "chpgvx.f"
	    } else {
#line 407 "chpgvx.f"
		*(unsigned char *)trans = 'N';
#line 408 "chpgvx.f"
	    }

#line 410 "chpgvx.f"
	    i__1 = *m;
#line 410 "chpgvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 411 "chpgvx.f"
		ctpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 413 "chpgvx.f"
/* L20: */
#line 413 "chpgvx.f"
	    }
#line 414 "chpgvx.f"
	}
#line 415 "chpgvx.f"
    }

#line 417 "chpgvx.f"
    return 0;

/*     End of CHPGVX */

} /* chpgvx_ */

