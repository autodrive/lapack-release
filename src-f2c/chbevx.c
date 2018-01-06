#line 1 "chbevx.f"
/* chbevx.f -- translated by f2c (version 20100827).
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

#line 1 "chbevx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static doublereal c_b16 = 1.;
static integer c__1 = 1;

/* > \brief <b> CHBEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, */
/*                          VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian band matrix A.  Eigenvalues and eigenvectors */
/* > can be selected by specifying either a range of values or a range of */
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
/* >          AB is COMPLEX array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, AB is overwritten by values generated during the */
/* >          reduction to tridiagonal form. */
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
/* >          Q is COMPLEX array, dimension (LDQ, N) */
/* >          If JOBZ = 'V', the N-by-N unitary matrix used in the */
/* >                          reduction to tridiagonal form. */
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
/* >          by reducing AB to tridiagonal form. */
/* > */
/* >          Eigenvalues will be computed most accurately when ABSTOL is */
/* >          set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
/* >          If this routine returns with INFO>0, indicating that some */
/* >          eigenvectors did not converge, try setting ABSTOL to */
/* >          2*SLAMCH('S'). */
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
/* >          W is REAL array, dimension (N) */
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, max(1,M)) */
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
/* >          > 0:  if INFO = i, then i eigenvectors failed to converge. */
/* >                Their indices are stored in array IFAIL. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int chbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublecomplex *q, 
	integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
	 integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork,
	 integer *ifail, integer *info, ftnlen jobz_len, ftnlen range_len, 
	ftnlen uplo_len)
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
    static doublecomplex ctmp1;
    static integer itmp1, indee;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static logical lower;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantz;
    extern doublereal clanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), chbtrd_(char *, char *, integer *,
	     integer *, doublecomplex *, integer *, doublereal *, doublereal *
	    , doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
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
    static doublereal smlnum;


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

#line 323 "chbevx.f"
    /* Parameter adjustments */
#line 323 "chbevx.f"
    ab_dim1 = *ldab;
#line 323 "chbevx.f"
    ab_offset = 1 + ab_dim1;
#line 323 "chbevx.f"
    ab -= ab_offset;
#line 323 "chbevx.f"
    q_dim1 = *ldq;
#line 323 "chbevx.f"
    q_offset = 1 + q_dim1;
#line 323 "chbevx.f"
    q -= q_offset;
#line 323 "chbevx.f"
    --w;
#line 323 "chbevx.f"
    z_dim1 = *ldz;
#line 323 "chbevx.f"
    z_offset = 1 + z_dim1;
#line 323 "chbevx.f"
    z__ -= z_offset;
#line 323 "chbevx.f"
    --work;
#line 323 "chbevx.f"
    --rwork;
#line 323 "chbevx.f"
    --iwork;
#line 323 "chbevx.f"
    --ifail;
#line 323 "chbevx.f"

#line 323 "chbevx.f"
    /* Function Body */
#line 323 "chbevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 324 "chbevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 325 "chbevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 326 "chbevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 327 "chbevx.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

#line 329 "chbevx.f"
    *info = 0;
#line 330 "chbevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 331 "chbevx.f"
	*info = -1;
#line 332 "chbevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 333 "chbevx.f"
	*info = -2;
#line 334 "chbevx.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 335 "chbevx.f"
	*info = -3;
#line 336 "chbevx.f"
    } else if (*n < 0) {
#line 337 "chbevx.f"
	*info = -4;
#line 338 "chbevx.f"
    } else if (*kd < 0) {
#line 339 "chbevx.f"
	*info = -5;
#line 340 "chbevx.f"
    } else if (*ldab < *kd + 1) {
#line 341 "chbevx.f"
	*info = -7;
#line 342 "chbevx.f"
    } else if (wantz && *ldq < max(1,*n)) {
#line 343 "chbevx.f"
	*info = -9;
#line 344 "chbevx.f"
    } else {
#line 345 "chbevx.f"
	if (valeig) {
#line 346 "chbevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 346 "chbevx.f"
		*info = -11;
#line 346 "chbevx.f"
	    }
#line 348 "chbevx.f"
	} else if (indeig) {
#line 349 "chbevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 350 "chbevx.f"
		*info = -12;
#line 351 "chbevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 352 "chbevx.f"
		*info = -13;
#line 353 "chbevx.f"
	    }
#line 354 "chbevx.f"
	}
#line 355 "chbevx.f"
    }
#line 356 "chbevx.f"
    if (*info == 0) {
#line 357 "chbevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 357 "chbevx.f"
	    *info = -18;
#line 357 "chbevx.f"
	}
#line 359 "chbevx.f"
    }

#line 361 "chbevx.f"
    if (*info != 0) {
#line 362 "chbevx.f"
	i__1 = -(*info);
#line 362 "chbevx.f"
	xerbla_("CHBEVX", &i__1, (ftnlen)6);
#line 363 "chbevx.f"
	return 0;
#line 364 "chbevx.f"
    }

/*     Quick return if possible */

#line 368 "chbevx.f"
    *m = 0;
#line 369 "chbevx.f"
    if (*n == 0) {
#line 369 "chbevx.f"
	return 0;
#line 369 "chbevx.f"
    }

#line 372 "chbevx.f"
    if (*n == 1) {
#line 373 "chbevx.f"
	*m = 1;
#line 374 "chbevx.f"
	if (lower) {
#line 375 "chbevx.f"
	    i__1 = ab_dim1 + 1;
#line 375 "chbevx.f"
	    ctmp1.r = ab[i__1].r, ctmp1.i = ab[i__1].i;
#line 376 "chbevx.f"
	} else {
#line 377 "chbevx.f"
	    i__1 = *kd + 1 + ab_dim1;
#line 377 "chbevx.f"
	    ctmp1.r = ab[i__1].r, ctmp1.i = ab[i__1].i;
#line 378 "chbevx.f"
	}
#line 379 "chbevx.f"
	tmp1 = ctmp1.r;
#line 380 "chbevx.f"
	if (valeig) {
#line 381 "chbevx.f"
	    if (! (*vl < tmp1 && *vu >= tmp1)) {
#line 381 "chbevx.f"
		*m = 0;
#line 381 "chbevx.f"
	    }
#line 383 "chbevx.f"
	}
#line 384 "chbevx.f"
	if (*m == 1) {
#line 385 "chbevx.f"
	    w[1] = ctmp1.r;
#line 386 "chbevx.f"
	    if (wantz) {
#line 386 "chbevx.f"
		i__1 = z_dim1 + 1;
#line 386 "chbevx.f"
		z__[i__1].r = 1., z__[i__1].i = 0.;
#line 386 "chbevx.f"
	    }
#line 388 "chbevx.f"
	}
#line 389 "chbevx.f"
	return 0;
#line 390 "chbevx.f"
    }

/*     Get machine constants. */

#line 394 "chbevx.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 395 "chbevx.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 396 "chbevx.f"
    smlnum = safmin / eps;
#line 397 "chbevx.f"
    bignum = 1. / smlnum;
#line 398 "chbevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 399 "chbevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 399 "chbevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 403 "chbevx.f"
    iscale = 0;
#line 404 "chbevx.f"
    abstll = *abstol;
#line 405 "chbevx.f"
    if (valeig) {
#line 406 "chbevx.f"
	vll = *vl;
#line 407 "chbevx.f"
	vuu = *vu;
#line 408 "chbevx.f"
    } else {
#line 409 "chbevx.f"
	vll = 0.;
#line 410 "chbevx.f"
	vuu = 0.;
#line 411 "chbevx.f"
    }
#line 412 "chbevx.f"
    anrm = clanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1, (ftnlen)1);
#line 413 "chbevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 414 "chbevx.f"
	iscale = 1;
#line 415 "chbevx.f"
	sigma = rmin / anrm;
#line 416 "chbevx.f"
    } else if (anrm > rmax) {
#line 417 "chbevx.f"
	iscale = 1;
#line 418 "chbevx.f"
	sigma = rmax / anrm;
#line 419 "chbevx.f"
    }
#line 420 "chbevx.f"
    if (iscale == 1) {
#line 421 "chbevx.f"
	if (lower) {
#line 422 "chbevx.f"
	    clascl_("B", kd, kd, &c_b16, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 423 "chbevx.f"
	} else {
#line 424 "chbevx.f"
	    clascl_("Q", kd, kd, &c_b16, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 425 "chbevx.f"
	}
#line 426 "chbevx.f"
	if (*abstol > 0.) {
#line 426 "chbevx.f"
	    abstll = *abstol * sigma;
#line 426 "chbevx.f"
	}
#line 428 "chbevx.f"
	if (valeig) {
#line 429 "chbevx.f"
	    vll = *vl * sigma;
#line 430 "chbevx.f"
	    vuu = *vu * sigma;
#line 431 "chbevx.f"
	}
#line 432 "chbevx.f"
    }

/*     Call CHBTRD to reduce Hermitian band matrix to tridiagonal form. */

#line 436 "chbevx.f"
    indd = 1;
#line 437 "chbevx.f"
    inde = indd + *n;
#line 438 "chbevx.f"
    indrwk = inde + *n;
#line 439 "chbevx.f"
    indwrk = 1;
#line 440 "chbevx.f"
    chbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &rwork[indd], &rwork[
	    inde], &q[q_offset], ldq, &work[indwrk], &iinfo, (ftnlen)1, (
	    ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call SSTERF or CSTEQR.  If this fails for some */
/*     eigenvalue, then try SSTEBZ. */

#line 447 "chbevx.f"
    test = FALSE_;
#line 448 "chbevx.f"
    if (indeig) {
#line 449 "chbevx.f"
	if (*il == 1 && *iu == *n) {
#line 450 "chbevx.f"
	    test = TRUE_;
#line 451 "chbevx.f"
	}
#line 452 "chbevx.f"
    }
#line 453 "chbevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 454 "chbevx.f"
	scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 455 "chbevx.f"
	indee = indrwk + (*n << 1);
#line 456 "chbevx.f"
	if (! wantz) {
#line 457 "chbevx.f"
	    i__1 = *n - 1;
#line 457 "chbevx.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 458 "chbevx.f"
	    ssterf_(n, &w[1], &rwork[indee], info);
#line 459 "chbevx.f"
	} else {
#line 460 "chbevx.f"
	    clacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 461 "chbevx.f"
	    i__1 = *n - 1;
#line 461 "chbevx.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 462 "chbevx.f"
	    csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 464 "chbevx.f"
	    if (*info == 0) {
#line 465 "chbevx.f"
		i__1 = *n;
#line 465 "chbevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 466 "chbevx.f"
		    ifail[i__] = 0;
#line 467 "chbevx.f"
/* L10: */
#line 467 "chbevx.f"
		}
#line 468 "chbevx.f"
	    }
#line 469 "chbevx.f"
	}
#line 470 "chbevx.f"
	if (*info == 0) {
#line 471 "chbevx.f"
	    *m = *n;
#line 472 "chbevx.f"
	    goto L30;
#line 473 "chbevx.f"
	}
#line 474 "chbevx.f"
	*info = 0;
#line 475 "chbevx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */

#line 479 "chbevx.f"
    if (wantz) {
#line 480 "chbevx.f"
	*(unsigned char *)order = 'B';
#line 481 "chbevx.f"
    } else {
#line 482 "chbevx.f"
	*(unsigned char *)order = 'E';
#line 483 "chbevx.f"
    }
#line 484 "chbevx.f"
    indibl = 1;
#line 485 "chbevx.f"
    indisp = indibl + *n;
#line 486 "chbevx.f"
    indiwk = indisp + *n;
#line 487 "chbevx.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 492 "chbevx.f"
    if (wantz) {
#line 493 "chbevx.f"
	cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by CSTEIN. */

#line 500 "chbevx.f"
	i__1 = *m;
#line 500 "chbevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 501 "chbevx.f"
	    ccopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 502 "chbevx.f"
	    cgemv_("N", n, n, &c_b2, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b1, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 504 "chbevx.f"
/* L20: */
#line 504 "chbevx.f"
	}
#line 505 "chbevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 509 "chbevx.f"
L30:
#line 510 "chbevx.f"
    if (iscale == 1) {
#line 511 "chbevx.f"
	if (*info == 0) {
#line 512 "chbevx.f"
	    imax = *m;
#line 513 "chbevx.f"
	} else {
#line 514 "chbevx.f"
	    imax = *info - 1;
#line 515 "chbevx.f"
	}
#line 516 "chbevx.f"
	d__1 = 1. / sigma;
#line 516 "chbevx.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 517 "chbevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 522 "chbevx.f"
    if (wantz) {
#line 523 "chbevx.f"
	i__1 = *m - 1;
#line 523 "chbevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 524 "chbevx.f"
	    i__ = 0;
#line 525 "chbevx.f"
	    tmp1 = w[j];
#line 526 "chbevx.f"
	    i__2 = *m;
#line 526 "chbevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 527 "chbevx.f"
		if (w[jj] < tmp1) {
#line 528 "chbevx.f"
		    i__ = jj;
#line 529 "chbevx.f"
		    tmp1 = w[jj];
#line 530 "chbevx.f"
		}
#line 531 "chbevx.f"
/* L40: */
#line 531 "chbevx.f"
	    }

#line 533 "chbevx.f"
	    if (i__ != 0) {
#line 534 "chbevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 535 "chbevx.f"
		w[i__] = w[j];
#line 536 "chbevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 537 "chbevx.f"
		w[j] = tmp1;
#line 538 "chbevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 539 "chbevx.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 540 "chbevx.f"
		if (*info != 0) {
#line 541 "chbevx.f"
		    itmp1 = ifail[i__];
#line 542 "chbevx.f"
		    ifail[i__] = ifail[j];
#line 543 "chbevx.f"
		    ifail[j] = itmp1;
#line 544 "chbevx.f"
		}
#line 545 "chbevx.f"
	    }
#line 546 "chbevx.f"
/* L50: */
#line 546 "chbevx.f"
	}
#line 547 "chbevx.f"
    }

#line 549 "chbevx.f"
    return 0;

/*     End of CHBEVX */

} /* chbevx_ */

