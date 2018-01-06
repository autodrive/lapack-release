#line 1 "sspgvx.f"
/* sspgvx.f -- translated by f2c (version 20100827).
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

#line 1 "sspgvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SSPGVX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPGVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, */
/*                          IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, */
/*                          IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               AP( * ), BP( * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPGVX computes selected eigenvalues, and optionally, eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A */
/* > and B are assumed to be symmetric, stored in packed storage, and B */
/* > is also positive definite.  Eigenvalues and eigenvectors can be */
/* > selected by specifying either a range of values or a range of indices */
/* > for the desired eigenvalues. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
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
/* >          BP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          B, packed columnwise in a linear array.  The j-th column of B */
/* >          is stored in the array BP as follows: */
/* >          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the triangular factor U or L from the Cholesky */
/* >          factorization B = U**T*U or B = L*L**T, in the same storage */
/* >          format as B. */
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
/* >          by reducing A to tridiagonal form. */
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
/* >          WORK is REAL array, dimension (8*N) */
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
/* >          > 0:  SPPTRF or SSPEVX returned an error code: */
/* >             <= N:  if INFO = i, SSPEVX failed to converge; */
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

/* > \date November 2015 */

/* > \ingroup realOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int sspgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *ap, doublereal *bp, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *iwork, integer *ifail, integer *info, ftnlen jobz_len, 
	ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;

    /* Local variables */
    static integer j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char trans[1];
    static logical upper, wantz;
    extern /* Subroutine */ int stpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    stpsv_(char *, char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static logical alleig, indeig, valeig;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), spptrf_(
	    char *, integer *, doublereal *, integer *, ftnlen), sspgst_(
	    integer *, char *, integer *, doublereal *, doublereal *, integer 
	    *, ftnlen), sspevx_(char *, char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

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

#line 303 "sspgvx.f"
    /* Parameter adjustments */
#line 303 "sspgvx.f"
    --ap;
#line 303 "sspgvx.f"
    --bp;
#line 303 "sspgvx.f"
    --w;
#line 303 "sspgvx.f"
    z_dim1 = *ldz;
#line 303 "sspgvx.f"
    z_offset = 1 + z_dim1;
#line 303 "sspgvx.f"
    z__ -= z_offset;
#line 303 "sspgvx.f"
    --work;
#line 303 "sspgvx.f"
    --iwork;
#line 303 "sspgvx.f"
    --ifail;
#line 303 "sspgvx.f"

#line 303 "sspgvx.f"
    /* Function Body */
#line 303 "sspgvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 304 "sspgvx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 305 "sspgvx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 306 "sspgvx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 307 "sspgvx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 309 "sspgvx.f"
    *info = 0;
#line 310 "sspgvx.f"
    if (*itype < 1 || *itype > 3) {
#line 311 "sspgvx.f"
	*info = -1;
#line 312 "sspgvx.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 313 "sspgvx.f"
	*info = -2;
#line 314 "sspgvx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 315 "sspgvx.f"
	*info = -3;
#line 316 "sspgvx.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 317 "sspgvx.f"
	*info = -4;
#line 318 "sspgvx.f"
    } else if (*n < 0) {
#line 319 "sspgvx.f"
	*info = -5;
#line 320 "sspgvx.f"
    } else {
#line 321 "sspgvx.f"
	if (valeig) {
#line 322 "sspgvx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 323 "sspgvx.f"
		*info = -9;
#line 324 "sspgvx.f"
	    }
#line 325 "sspgvx.f"
	} else if (indeig) {
#line 326 "sspgvx.f"
	    if (*il < 1) {
#line 327 "sspgvx.f"
		*info = -10;
#line 328 "sspgvx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 329 "sspgvx.f"
		*info = -11;
#line 330 "sspgvx.f"
	    }
#line 331 "sspgvx.f"
	}
#line 332 "sspgvx.f"
    }
#line 333 "sspgvx.f"
    if (*info == 0) {
#line 334 "sspgvx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 335 "sspgvx.f"
	    *info = -16;
#line 336 "sspgvx.f"
	}
#line 337 "sspgvx.f"
    }

#line 339 "sspgvx.f"
    if (*info != 0) {
#line 340 "sspgvx.f"
	i__1 = -(*info);
#line 340 "sspgvx.f"
	xerbla_("SSPGVX", &i__1, (ftnlen)6);
#line 341 "sspgvx.f"
	return 0;
#line 342 "sspgvx.f"
    }

/*     Quick return if possible */

#line 346 "sspgvx.f"
    *m = 0;
#line 347 "sspgvx.f"
    if (*n == 0) {
#line 347 "sspgvx.f"
	return 0;
#line 347 "sspgvx.f"
    }

/*     Form a Cholesky factorization of B. */

#line 352 "sspgvx.f"
    spptrf_(uplo, n, &bp[1], info, (ftnlen)1);
#line 353 "sspgvx.f"
    if (*info != 0) {
#line 354 "sspgvx.f"
	*info = *n + *info;
#line 355 "sspgvx.f"
	return 0;
#line 356 "sspgvx.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 360 "sspgvx.f"
    sspgst_(itype, uplo, n, &ap[1], &bp[1], info, (ftnlen)1);
#line 361 "sspgvx.f"
    sspevx_(jobz, range, uplo, n, &ap[1], vl, vu, il, iu, abstol, m, &w[1], &
	    z__[z_offset], ldz, &work[1], &iwork[1], &ifail[1], info, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1);

#line 364 "sspgvx.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 368 "sspgvx.f"
	if (*info > 0) {
#line 368 "sspgvx.f"
	    *m = *info - 1;
#line 368 "sspgvx.f"
	}
#line 370 "sspgvx.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */

#line 375 "sspgvx.f"
	    if (upper) {
#line 376 "sspgvx.f"
		*(unsigned char *)trans = 'N';
#line 377 "sspgvx.f"
	    } else {
#line 378 "sspgvx.f"
		*(unsigned char *)trans = 'T';
#line 379 "sspgvx.f"
	    }

#line 381 "sspgvx.f"
	    i__1 = *m;
#line 381 "sspgvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 382 "sspgvx.f"
		stpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 384 "sspgvx.f"
/* L10: */
#line 384 "sspgvx.f"
	    }

#line 386 "sspgvx.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T*y */

#line 391 "sspgvx.f"
	    if (upper) {
#line 392 "sspgvx.f"
		*(unsigned char *)trans = 'T';
#line 393 "sspgvx.f"
	    } else {
#line 394 "sspgvx.f"
		*(unsigned char *)trans = 'N';
#line 395 "sspgvx.f"
	    }

#line 397 "sspgvx.f"
	    i__1 = *m;
#line 397 "sspgvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 398 "sspgvx.f"
		stpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 400 "sspgvx.f"
/* L20: */
#line 400 "sspgvx.f"
	    }
#line 401 "sspgvx.f"
	}
#line 402 "sspgvx.f"
    }

#line 404 "sspgvx.f"
    return 0;

/*     End of SSPGVX */

} /* sspgvx_ */

