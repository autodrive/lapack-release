#line 1 "ssyevx.f"
/* ssyevx.f -- translated by f2c (version 20100827).
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

#line 1 "ssyevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> SSYEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, */
/*                          IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A.  Eigenvalues and eigenvectors can be */
/* > selected by specifying either a range of values or a range of indices */
/* > for the desired eigenvalues. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA, N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* >          On exit, the lower triangle (if UPLO='L') or the upper */
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
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
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
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= 1, when N <= 1; */
/* >          otherwise 8*N. */
/* >          For optimal efficiency, LWORK >= (NB+3)*N, */
/* >          where NB is the max of the blocksize for SSYTRD and SORMTR */
/* >          returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \date November 2011 */

/* > \ingroup realSYeigen */

/*  ===================================================================== */
/* Subroutine */ int ssyevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
	integer *iwork, integer *ifail, integer *info, ftnlen jobz_len, 
	ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nb, jj;
    static doublereal eps, vll, vuu, tmp1;
    static integer indd, inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static integer itmp1, indee;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    static logical lower;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical wantz, alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp, indiwo, indwkn;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer indwrk, lwkmin;
    extern /* Subroutine */ int sstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    ssterf_(integer *, doublereal *, doublereal *, integer *);
    static integer llwrkn, llwork, nsplit;
    static doublereal smlnum;
    extern doublereal slansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sorgtr_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), ssteqr_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    sormtr_(char *, char *, char *, integer *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), ssytrd_(char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

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

#line 298 "ssyevx.f"
    /* Parameter adjustments */
#line 298 "ssyevx.f"
    a_dim1 = *lda;
#line 298 "ssyevx.f"
    a_offset = 1 + a_dim1;
#line 298 "ssyevx.f"
    a -= a_offset;
#line 298 "ssyevx.f"
    --w;
#line 298 "ssyevx.f"
    z_dim1 = *ldz;
#line 298 "ssyevx.f"
    z_offset = 1 + z_dim1;
#line 298 "ssyevx.f"
    z__ -= z_offset;
#line 298 "ssyevx.f"
    --work;
#line 298 "ssyevx.f"
    --iwork;
#line 298 "ssyevx.f"
    --ifail;
#line 298 "ssyevx.f"

#line 298 "ssyevx.f"
    /* Function Body */
#line 298 "ssyevx.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 299 "ssyevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 300 "ssyevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 301 "ssyevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 302 "ssyevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 303 "ssyevx.f"
    lquery = *lwork == -1;

#line 305 "ssyevx.f"
    *info = 0;
#line 306 "ssyevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 307 "ssyevx.f"
	*info = -1;
#line 308 "ssyevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 309 "ssyevx.f"
	*info = -2;
#line 310 "ssyevx.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 311 "ssyevx.f"
	*info = -3;
#line 312 "ssyevx.f"
    } else if (*n < 0) {
#line 313 "ssyevx.f"
	*info = -4;
#line 314 "ssyevx.f"
    } else if (*lda < max(1,*n)) {
#line 315 "ssyevx.f"
	*info = -6;
#line 316 "ssyevx.f"
    } else {
#line 317 "ssyevx.f"
	if (valeig) {
#line 318 "ssyevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 318 "ssyevx.f"
		*info = -8;
#line 318 "ssyevx.f"
	    }
#line 320 "ssyevx.f"
	} else if (indeig) {
#line 321 "ssyevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 322 "ssyevx.f"
		*info = -9;
#line 323 "ssyevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 324 "ssyevx.f"
		*info = -10;
#line 325 "ssyevx.f"
	    }
#line 326 "ssyevx.f"
	}
#line 327 "ssyevx.f"
    }
#line 328 "ssyevx.f"
    if (*info == 0) {
#line 329 "ssyevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 330 "ssyevx.f"
	    *info = -15;
#line 331 "ssyevx.f"
	}
#line 332 "ssyevx.f"
    }

#line 334 "ssyevx.f"
    if (*info == 0) {
#line 335 "ssyevx.f"
	if (*n <= 1) {
#line 336 "ssyevx.f"
	    lwkmin = 1;
#line 337 "ssyevx.f"
	    work[1] = (doublereal) lwkmin;
#line 338 "ssyevx.f"
	} else {
#line 339 "ssyevx.f"
	    lwkmin = *n << 3;
#line 340 "ssyevx.f"
	    nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 341 "ssyevx.f"
	    i__1 = nb, i__2 = ilaenv_(&c__1, "SORMTR", uplo, n, &c_n1, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1);
#line 341 "ssyevx.f"
	    nb = max(i__1,i__2);
/* Computing MAX */
#line 342 "ssyevx.f"
	    i__1 = lwkmin, i__2 = (nb + 3) * *n;
#line 342 "ssyevx.f"
	    lwkopt = max(i__1,i__2);
#line 343 "ssyevx.f"
	    work[1] = (doublereal) lwkopt;
#line 344 "ssyevx.f"
	}

#line 346 "ssyevx.f"
	if (*lwork < lwkmin && ! lquery) {
#line 346 "ssyevx.f"
	    *info = -17;
#line 346 "ssyevx.f"
	}
#line 348 "ssyevx.f"
    }

#line 350 "ssyevx.f"
    if (*info != 0) {
#line 351 "ssyevx.f"
	i__1 = -(*info);
#line 351 "ssyevx.f"
	xerbla_("SSYEVX", &i__1, (ftnlen)6);
#line 352 "ssyevx.f"
	return 0;
#line 353 "ssyevx.f"
    } else if (lquery) {
#line 354 "ssyevx.f"
	return 0;
#line 355 "ssyevx.f"
    }

/*     Quick return if possible */

#line 359 "ssyevx.f"
    *m = 0;
#line 360 "ssyevx.f"
    if (*n == 0) {
#line 361 "ssyevx.f"
	return 0;
#line 362 "ssyevx.f"
    }

#line 364 "ssyevx.f"
    if (*n == 1) {
#line 365 "ssyevx.f"
	if (alleig || indeig) {
#line 366 "ssyevx.f"
	    *m = 1;
#line 367 "ssyevx.f"
	    w[1] = a[a_dim1 + 1];
#line 368 "ssyevx.f"
	} else {
#line 369 "ssyevx.f"
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
#line 370 "ssyevx.f"
		*m = 1;
#line 371 "ssyevx.f"
		w[1] = a[a_dim1 + 1];
#line 372 "ssyevx.f"
	    }
#line 373 "ssyevx.f"
	}
#line 374 "ssyevx.f"
	if (wantz) {
#line 374 "ssyevx.f"
	    z__[z_dim1 + 1] = 1.;
#line 374 "ssyevx.f"
	}
#line 376 "ssyevx.f"
	return 0;
#line 377 "ssyevx.f"
    }

/*     Get machine constants. */

#line 381 "ssyevx.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 382 "ssyevx.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 383 "ssyevx.f"
    smlnum = safmin / eps;
#line 384 "ssyevx.f"
    bignum = 1. / smlnum;
#line 385 "ssyevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 386 "ssyevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 386 "ssyevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 390 "ssyevx.f"
    iscale = 0;
#line 391 "ssyevx.f"
    abstll = *abstol;
#line 392 "ssyevx.f"
    if (valeig) {
#line 393 "ssyevx.f"
	vll = *vl;
#line 394 "ssyevx.f"
	vuu = *vu;
#line 395 "ssyevx.f"
    }
#line 396 "ssyevx.f"
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 397 "ssyevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 398 "ssyevx.f"
	iscale = 1;
#line 399 "ssyevx.f"
	sigma = rmin / anrm;
#line 400 "ssyevx.f"
    } else if (anrm > rmax) {
#line 401 "ssyevx.f"
	iscale = 1;
#line 402 "ssyevx.f"
	sigma = rmax / anrm;
#line 403 "ssyevx.f"
    }
#line 404 "ssyevx.f"
    if (iscale == 1) {
#line 405 "ssyevx.f"
	if (lower) {
#line 406 "ssyevx.f"
	    i__1 = *n;
#line 406 "ssyevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 407 "ssyevx.f"
		i__2 = *n - j + 1;
#line 407 "ssyevx.f"
		sscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 408 "ssyevx.f"
/* L10: */
#line 408 "ssyevx.f"
	    }
#line 409 "ssyevx.f"
	} else {
#line 410 "ssyevx.f"
	    i__1 = *n;
#line 410 "ssyevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 411 "ssyevx.f"
		sscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 412 "ssyevx.f"
/* L20: */
#line 412 "ssyevx.f"
	    }
#line 413 "ssyevx.f"
	}
#line 414 "ssyevx.f"
	if (*abstol > 0.) {
#line 414 "ssyevx.f"
	    abstll = *abstol * sigma;
#line 414 "ssyevx.f"
	}
#line 416 "ssyevx.f"
	if (valeig) {
#line 417 "ssyevx.f"
	    vll = *vl * sigma;
#line 418 "ssyevx.f"
	    vuu = *vu * sigma;
#line 419 "ssyevx.f"
	}
#line 420 "ssyevx.f"
    }

/*     Call SSYTRD to reduce symmetric matrix to tridiagonal form. */

#line 424 "ssyevx.f"
    indtau = 1;
#line 425 "ssyevx.f"
    inde = indtau + *n;
#line 426 "ssyevx.f"
    indd = inde + *n;
#line 427 "ssyevx.f"
    indwrk = indd + *n;
#line 428 "ssyevx.f"
    llwork = *lwork - indwrk + 1;
#line 429 "ssyevx.f"
    ssytrd_(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call SSTERF or SORGTR and SSTEQR.  If this fails for */
/*     some eigenvalue, then try SSTEBZ. */

#line 436 "ssyevx.f"
    test = FALSE_;
#line 437 "ssyevx.f"
    if (indeig) {
#line 438 "ssyevx.f"
	if (*il == 1 && *iu == *n) {
#line 439 "ssyevx.f"
	    test = TRUE_;
#line 440 "ssyevx.f"
	}
#line 441 "ssyevx.f"
    }
#line 442 "ssyevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 443 "ssyevx.f"
	scopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 444 "ssyevx.f"
	indee = indwrk + (*n << 1);
#line 445 "ssyevx.f"
	if (! wantz) {
#line 446 "ssyevx.f"
	    i__1 = *n - 1;
#line 446 "ssyevx.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 447 "ssyevx.f"
	    ssterf_(n, &w[1], &work[indee], info);
#line 448 "ssyevx.f"
	} else {
#line 449 "ssyevx.f"
	    slacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 450 "ssyevx.f"
	    sorgtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
#line 452 "ssyevx.f"
	    i__1 = *n - 1;
#line 452 "ssyevx.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 453 "ssyevx.f"
	    ssteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 455 "ssyevx.f"
	    if (*info == 0) {
#line 456 "ssyevx.f"
		i__1 = *n;
#line 456 "ssyevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 457 "ssyevx.f"
		    ifail[i__] = 0;
#line 458 "ssyevx.f"
/* L30: */
#line 458 "ssyevx.f"
		}
#line 459 "ssyevx.f"
	    }
#line 460 "ssyevx.f"
	}
#line 461 "ssyevx.f"
	if (*info == 0) {
#line 462 "ssyevx.f"
	    *m = *n;
#line 463 "ssyevx.f"
	    goto L40;
#line 464 "ssyevx.f"
	}
#line 465 "ssyevx.f"
	*info = 0;
#line 466 "ssyevx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 470 "ssyevx.f"
    if (wantz) {
#line 471 "ssyevx.f"
	*(unsigned char *)order = 'B';
#line 472 "ssyevx.f"
    } else {
#line 473 "ssyevx.f"
	*(unsigned char *)order = 'E';
#line 474 "ssyevx.f"
    }
#line 475 "ssyevx.f"
    indibl = 1;
#line 476 "ssyevx.f"
    indisp = indibl + *n;
#line 477 "ssyevx.f"
    indiwo = indisp + *n;
#line 478 "ssyevx.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 483 "ssyevx.f"
    if (wantz) {
#line 484 "ssyevx.f"
	sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEIN. */

#line 491 "ssyevx.f"
	indwkn = inde;
#line 492 "ssyevx.f"
	llwrkn = *lwork - indwkn + 1;
#line 493 "ssyevx.f"
	sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 495 "ssyevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 499 "ssyevx.f"
L40:
#line 500 "ssyevx.f"
    if (iscale == 1) {
#line 501 "ssyevx.f"
	if (*info == 0) {
#line 502 "ssyevx.f"
	    imax = *m;
#line 503 "ssyevx.f"
	} else {
#line 504 "ssyevx.f"
	    imax = *info - 1;
#line 505 "ssyevx.f"
	}
#line 506 "ssyevx.f"
	d__1 = 1. / sigma;
#line 506 "ssyevx.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 507 "ssyevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 512 "ssyevx.f"
    if (wantz) {
#line 513 "ssyevx.f"
	i__1 = *m - 1;
#line 513 "ssyevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 514 "ssyevx.f"
	    i__ = 0;
#line 515 "ssyevx.f"
	    tmp1 = w[j];
#line 516 "ssyevx.f"
	    i__2 = *m;
#line 516 "ssyevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 517 "ssyevx.f"
		if (w[jj] < tmp1) {
#line 518 "ssyevx.f"
		    i__ = jj;
#line 519 "ssyevx.f"
		    tmp1 = w[jj];
#line 520 "ssyevx.f"
		}
#line 521 "ssyevx.f"
/* L50: */
#line 521 "ssyevx.f"
	    }

#line 523 "ssyevx.f"
	    if (i__ != 0) {
#line 524 "ssyevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 525 "ssyevx.f"
		w[i__] = w[j];
#line 526 "ssyevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 527 "ssyevx.f"
		w[j] = tmp1;
#line 528 "ssyevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 529 "ssyevx.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 530 "ssyevx.f"
		if (*info != 0) {
#line 531 "ssyevx.f"
		    itmp1 = ifail[i__];
#line 532 "ssyevx.f"
		    ifail[i__] = ifail[j];
#line 533 "ssyevx.f"
		    ifail[j] = itmp1;
#line 534 "ssyevx.f"
		}
#line 535 "ssyevx.f"
	    }
#line 536 "ssyevx.f"
/* L60: */
#line 536 "ssyevx.f"
	}
#line 537 "ssyevx.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 541 "ssyevx.f"
    work[1] = (doublereal) lwkopt;

#line 543 "ssyevx.f"
    return 0;

/*     End of SSYEVX */

} /* ssyevx_ */

