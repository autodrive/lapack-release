#line 1 "dsbevx.f"
/* dsbevx.f -- translated by f2c (version 20100827).
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

#line 1 "dsbevx.f"
/* Table of constant values */

static doublereal c_b14 = 1.;
static integer c__1 = 1;
static doublereal c_b34 = 0.;

/* > \brief <b> DSBEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, */
/*                          VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, */
/*                          IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric band matrix A.  Eigenvalues and eigenvectors can */
/* > be selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
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
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, AB is overwritten by values generated during the */
/* >          reduction to tridiagonal form.  If UPLO = 'U', the first */
/* >          superdiagonal and the diagonal of the tridiagonal matrix T */
/* >          are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
/* >          the diagonal and first subdiagonal of T are returned in the */
/* >          first two rows of AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD + 1. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* >          If JOBZ = 'V', the N-by-N orthogonal matrix used in the */
/* >                         reduction to tridiagonal form. */
/* >          If JOBZ = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  If JOBZ = 'V', then */
/* >          LDQ >= max(1,N). */
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
/* >          by reducing AB to tridiagonal form. */
/* > */
/* >          Eigenvalues will be computed most accurately when ABSTOL is */
/* >          set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
/* >          If this routine returns with INFO>0, indicating that some */
/* >          eigenvectors did not converge, try setting ABSTOL to */
/* >          2*DLAMCH('S'). */
/* > */
/* >          See "Computing Small Singular Values of Bidiagonal Matrices */
/* >          with Guaranteed High Relative Accuracy," by Demmel and */
/* >          Kahan, LAPACK Working Note #3. */
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
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M)) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          If an eigenvector fails to converge, then that column of Z */
/* >          contains the latest approximation to the eigenvector, and the */
/* >          index of the eigenvector is returned in IFAIL. */
/* >          If JOBZ = 'N', then Z is not referenced. */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, then i eigenvectors failed to converge. */
/* >                Their indices are stored in array IFAIL. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dsbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *
	ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	integer *info, ftnlen jobz_len, ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, 
	    i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, jj;
    static doublereal eps, vll, vuu, tmp1;
    static integer indd, inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static integer itmp1, indee;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lower, wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern doublereal dlansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    static logical valeig;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    extern /* Subroutine */ int dsbtrd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen, ftnlen);
    static integer indisp;
    extern /* Subroutine */ int dstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dsterf_(integer *, doublereal *, doublereal *, integer *);
    static integer indiwo;
    extern /* Subroutine */ int dstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer nsplit;
    static doublereal smlnum;


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

#line 308 "dsbevx.f"
    /* Parameter adjustments */
#line 308 "dsbevx.f"
    ab_dim1 = *ldab;
#line 308 "dsbevx.f"
    ab_offset = 1 + ab_dim1;
#line 308 "dsbevx.f"
    ab -= ab_offset;
#line 308 "dsbevx.f"
    q_dim1 = *ldq;
#line 308 "dsbevx.f"
    q_offset = 1 + q_dim1;
#line 308 "dsbevx.f"
    q -= q_offset;
#line 308 "dsbevx.f"
    --w;
#line 308 "dsbevx.f"
    z_dim1 = *ldz;
#line 308 "dsbevx.f"
    z_offset = 1 + z_dim1;
#line 308 "dsbevx.f"
    z__ -= z_offset;
#line 308 "dsbevx.f"
    --work;
#line 308 "dsbevx.f"
    --iwork;
#line 308 "dsbevx.f"
    --ifail;
#line 308 "dsbevx.f"

#line 308 "dsbevx.f"
    /* Function Body */
#line 308 "dsbevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 309 "dsbevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 310 "dsbevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 311 "dsbevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 312 "dsbevx.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

#line 314 "dsbevx.f"
    *info = 0;
#line 315 "dsbevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 316 "dsbevx.f"
	*info = -1;
#line 317 "dsbevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 318 "dsbevx.f"
	*info = -2;
#line 319 "dsbevx.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 320 "dsbevx.f"
	*info = -3;
#line 321 "dsbevx.f"
    } else if (*n < 0) {
#line 322 "dsbevx.f"
	*info = -4;
#line 323 "dsbevx.f"
    } else if (*kd < 0) {
#line 324 "dsbevx.f"
	*info = -5;
#line 325 "dsbevx.f"
    } else if (*ldab < *kd + 1) {
#line 326 "dsbevx.f"
	*info = -7;
#line 327 "dsbevx.f"
    } else if (wantz && *ldq < max(1,*n)) {
#line 328 "dsbevx.f"
	*info = -9;
#line 329 "dsbevx.f"
    } else {
#line 330 "dsbevx.f"
	if (valeig) {
#line 331 "dsbevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 331 "dsbevx.f"
		*info = -11;
#line 331 "dsbevx.f"
	    }
#line 333 "dsbevx.f"
	} else if (indeig) {
#line 334 "dsbevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 335 "dsbevx.f"
		*info = -12;
#line 336 "dsbevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 337 "dsbevx.f"
		*info = -13;
#line 338 "dsbevx.f"
	    }
#line 339 "dsbevx.f"
	}
#line 340 "dsbevx.f"
    }
#line 341 "dsbevx.f"
    if (*info == 0) {
#line 342 "dsbevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 342 "dsbevx.f"
	    *info = -18;
#line 342 "dsbevx.f"
	}
#line 344 "dsbevx.f"
    }

#line 346 "dsbevx.f"
    if (*info != 0) {
#line 347 "dsbevx.f"
	i__1 = -(*info);
#line 347 "dsbevx.f"
	xerbla_("DSBEVX", &i__1, (ftnlen)6);
#line 348 "dsbevx.f"
	return 0;
#line 349 "dsbevx.f"
    }

/*     Quick return if possible */

#line 353 "dsbevx.f"
    *m = 0;
#line 354 "dsbevx.f"
    if (*n == 0) {
#line 354 "dsbevx.f"
	return 0;
#line 354 "dsbevx.f"
    }

#line 357 "dsbevx.f"
    if (*n == 1) {
#line 358 "dsbevx.f"
	*m = 1;
#line 359 "dsbevx.f"
	if (lower) {
#line 360 "dsbevx.f"
	    tmp1 = ab[ab_dim1 + 1];
#line 361 "dsbevx.f"
	} else {
#line 362 "dsbevx.f"
	    tmp1 = ab[*kd + 1 + ab_dim1];
#line 363 "dsbevx.f"
	}
#line 364 "dsbevx.f"
	if (valeig) {
#line 365 "dsbevx.f"
	    if (! (*vl < tmp1 && *vu >= tmp1)) {
#line 365 "dsbevx.f"
		*m = 0;
#line 365 "dsbevx.f"
	    }
#line 367 "dsbevx.f"
	}
#line 368 "dsbevx.f"
	if (*m == 1) {
#line 369 "dsbevx.f"
	    w[1] = tmp1;
#line 370 "dsbevx.f"
	    if (wantz) {
#line 370 "dsbevx.f"
		z__[z_dim1 + 1] = 1.;
#line 370 "dsbevx.f"
	    }
#line 372 "dsbevx.f"
	}
#line 373 "dsbevx.f"
	return 0;
#line 374 "dsbevx.f"
    }

/*     Get machine constants. */

#line 378 "dsbevx.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 379 "dsbevx.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 380 "dsbevx.f"
    smlnum = safmin / eps;
#line 381 "dsbevx.f"
    bignum = 1. / smlnum;
#line 382 "dsbevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 383 "dsbevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 383 "dsbevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 387 "dsbevx.f"
    iscale = 0;
#line 388 "dsbevx.f"
    abstll = *abstol;
#line 389 "dsbevx.f"
    if (valeig) {
#line 390 "dsbevx.f"
	vll = *vl;
#line 391 "dsbevx.f"
	vuu = *vu;
#line 392 "dsbevx.f"
    } else {
#line 393 "dsbevx.f"
	vll = 0.;
#line 394 "dsbevx.f"
	vuu = 0.;
#line 395 "dsbevx.f"
    }
#line 396 "dsbevx.f"
    anrm = dlansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 397 "dsbevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 398 "dsbevx.f"
	iscale = 1;
#line 399 "dsbevx.f"
	sigma = rmin / anrm;
#line 400 "dsbevx.f"
    } else if (anrm > rmax) {
#line 401 "dsbevx.f"
	iscale = 1;
#line 402 "dsbevx.f"
	sigma = rmax / anrm;
#line 403 "dsbevx.f"
    }
#line 404 "dsbevx.f"
    if (iscale == 1) {
#line 405 "dsbevx.f"
	if (lower) {
#line 406 "dsbevx.f"
	    dlascl_("B", kd, kd, &c_b14, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 407 "dsbevx.f"
	} else {
#line 408 "dsbevx.f"
	    dlascl_("Q", kd, kd, &c_b14, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 409 "dsbevx.f"
	}
#line 410 "dsbevx.f"
	if (*abstol > 0.) {
#line 410 "dsbevx.f"
	    abstll = *abstol * sigma;
#line 410 "dsbevx.f"
	}
#line 412 "dsbevx.f"
	if (valeig) {
#line 413 "dsbevx.f"
	    vll = *vl * sigma;
#line 414 "dsbevx.f"
	    vuu = *vu * sigma;
#line 415 "dsbevx.f"
	}
#line 416 "dsbevx.f"
    }

/*     Call DSBTRD to reduce symmetric band matrix to tridiagonal form. */

#line 420 "dsbevx.f"
    indd = 1;
#line 421 "dsbevx.f"
    inde = indd + *n;
#line 422 "dsbevx.f"
    indwrk = inde + *n;
#line 423 "dsbevx.f"
    dsbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &work[indd], &work[inde],
	     &q[q_offset], ldq, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call DSTERF or SSTEQR.  If this fails for some */
/*     eigenvalue, then try DSTEBZ. */

#line 430 "dsbevx.f"
    test = FALSE_;
#line 431 "dsbevx.f"
    if (indeig) {
#line 432 "dsbevx.f"
	if (*il == 1 && *iu == *n) {
#line 433 "dsbevx.f"
	    test = TRUE_;
#line 434 "dsbevx.f"
	}
#line 435 "dsbevx.f"
    }
#line 436 "dsbevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 437 "dsbevx.f"
	dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 438 "dsbevx.f"
	indee = indwrk + (*n << 1);
#line 439 "dsbevx.f"
	if (! wantz) {
#line 440 "dsbevx.f"
	    i__1 = *n - 1;
#line 440 "dsbevx.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 441 "dsbevx.f"
	    dsterf_(n, &w[1], &work[indee], info);
#line 442 "dsbevx.f"
	} else {
#line 443 "dsbevx.f"
	    dlacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 444 "dsbevx.f"
	    i__1 = *n - 1;
#line 444 "dsbevx.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 445 "dsbevx.f"
	    dsteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 447 "dsbevx.f"
	    if (*info == 0) {
#line 448 "dsbevx.f"
		i__1 = *n;
#line 448 "dsbevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 449 "dsbevx.f"
		    ifail[i__] = 0;
#line 450 "dsbevx.f"
/* L10: */
#line 450 "dsbevx.f"
		}
#line 451 "dsbevx.f"
	    }
#line 452 "dsbevx.f"
	}
#line 453 "dsbevx.f"
	if (*info == 0) {
#line 454 "dsbevx.f"
	    *m = *n;
#line 455 "dsbevx.f"
	    goto L30;
#line 456 "dsbevx.f"
	}
#line 457 "dsbevx.f"
	*info = 0;
#line 458 "dsbevx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 462 "dsbevx.f"
    if (wantz) {
#line 463 "dsbevx.f"
	*(unsigned char *)order = 'B';
#line 464 "dsbevx.f"
    } else {
#line 465 "dsbevx.f"
	*(unsigned char *)order = 'E';
#line 466 "dsbevx.f"
    }
#line 467 "dsbevx.f"
    indibl = 1;
#line 468 "dsbevx.f"
    indisp = indibl + *n;
#line 469 "dsbevx.f"
    indiwo = indisp + *n;
#line 470 "dsbevx.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 475 "dsbevx.f"
    if (wantz) {
#line 476 "dsbevx.f"
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

#line 483 "dsbevx.f"
	i__1 = *m;
#line 483 "dsbevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 484 "dsbevx.f"
	    dcopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 485 "dsbevx.f"
	    dgemv_("N", n, n, &c_b14, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b34, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 487 "dsbevx.f"
/* L20: */
#line 487 "dsbevx.f"
	}
#line 488 "dsbevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 492 "dsbevx.f"
L30:
#line 493 "dsbevx.f"
    if (iscale == 1) {
#line 494 "dsbevx.f"
	if (*info == 0) {
#line 495 "dsbevx.f"
	    imax = *m;
#line 496 "dsbevx.f"
	} else {
#line 497 "dsbevx.f"
	    imax = *info - 1;
#line 498 "dsbevx.f"
	}
#line 499 "dsbevx.f"
	d__1 = 1. / sigma;
#line 499 "dsbevx.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 500 "dsbevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 505 "dsbevx.f"
    if (wantz) {
#line 506 "dsbevx.f"
	i__1 = *m - 1;
#line 506 "dsbevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 507 "dsbevx.f"
	    i__ = 0;
#line 508 "dsbevx.f"
	    tmp1 = w[j];
#line 509 "dsbevx.f"
	    i__2 = *m;
#line 509 "dsbevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 510 "dsbevx.f"
		if (w[jj] < tmp1) {
#line 511 "dsbevx.f"
		    i__ = jj;
#line 512 "dsbevx.f"
		    tmp1 = w[jj];
#line 513 "dsbevx.f"
		}
#line 514 "dsbevx.f"
/* L40: */
#line 514 "dsbevx.f"
	    }

#line 516 "dsbevx.f"
	    if (i__ != 0) {
#line 517 "dsbevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 518 "dsbevx.f"
		w[i__] = w[j];
#line 519 "dsbevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 520 "dsbevx.f"
		w[j] = tmp1;
#line 521 "dsbevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 522 "dsbevx.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 523 "dsbevx.f"
		if (*info != 0) {
#line 524 "dsbevx.f"
		    itmp1 = ifail[i__];
#line 525 "dsbevx.f"
		    ifail[i__] = ifail[j];
#line 526 "dsbevx.f"
		    ifail[j] = itmp1;
#line 527 "dsbevx.f"
		}
#line 528 "dsbevx.f"
	    }
#line 529 "dsbevx.f"
/* L50: */
#line 529 "dsbevx.f"
	}
#line 530 "dsbevx.f"
    }

#line 532 "dsbevx.f"
    return 0;

/*     End of DSBEVX */

} /* dsbevx_ */

