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

/* > \date June 2016 */

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


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 305 "ssyevx.f"
    /* Parameter adjustments */
#line 305 "ssyevx.f"
    a_dim1 = *lda;
#line 305 "ssyevx.f"
    a_offset = 1 + a_dim1;
#line 305 "ssyevx.f"
    a -= a_offset;
#line 305 "ssyevx.f"
    --w;
#line 305 "ssyevx.f"
    z_dim1 = *ldz;
#line 305 "ssyevx.f"
    z_offset = 1 + z_dim1;
#line 305 "ssyevx.f"
    z__ -= z_offset;
#line 305 "ssyevx.f"
    --work;
#line 305 "ssyevx.f"
    --iwork;
#line 305 "ssyevx.f"
    --ifail;
#line 305 "ssyevx.f"

#line 305 "ssyevx.f"
    /* Function Body */
#line 305 "ssyevx.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 306 "ssyevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 307 "ssyevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 308 "ssyevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 309 "ssyevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 310 "ssyevx.f"
    lquery = *lwork == -1;

#line 312 "ssyevx.f"
    *info = 0;
#line 313 "ssyevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 314 "ssyevx.f"
	*info = -1;
#line 315 "ssyevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 316 "ssyevx.f"
	*info = -2;
#line 317 "ssyevx.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 318 "ssyevx.f"
	*info = -3;
#line 319 "ssyevx.f"
    } else if (*n < 0) {
#line 320 "ssyevx.f"
	*info = -4;
#line 321 "ssyevx.f"
    } else if (*lda < max(1,*n)) {
#line 322 "ssyevx.f"
	*info = -6;
#line 323 "ssyevx.f"
    } else {
#line 324 "ssyevx.f"
	if (valeig) {
#line 325 "ssyevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 325 "ssyevx.f"
		*info = -8;
#line 325 "ssyevx.f"
	    }
#line 327 "ssyevx.f"
	} else if (indeig) {
#line 328 "ssyevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 329 "ssyevx.f"
		*info = -9;
#line 330 "ssyevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 331 "ssyevx.f"
		*info = -10;
#line 332 "ssyevx.f"
	    }
#line 333 "ssyevx.f"
	}
#line 334 "ssyevx.f"
    }
#line 335 "ssyevx.f"
    if (*info == 0) {
#line 336 "ssyevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 337 "ssyevx.f"
	    *info = -15;
#line 338 "ssyevx.f"
	}
#line 339 "ssyevx.f"
    }

#line 341 "ssyevx.f"
    if (*info == 0) {
#line 342 "ssyevx.f"
	if (*n <= 1) {
#line 343 "ssyevx.f"
	    lwkmin = 1;
#line 344 "ssyevx.f"
	    work[1] = (doublereal) lwkmin;
#line 345 "ssyevx.f"
	} else {
#line 346 "ssyevx.f"
	    lwkmin = *n << 3;
#line 347 "ssyevx.f"
	    nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 348 "ssyevx.f"
	    i__1 = nb, i__2 = ilaenv_(&c__1, "SORMTR", uplo, n, &c_n1, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1);
#line 348 "ssyevx.f"
	    nb = max(i__1,i__2);
/* Computing MAX */
#line 349 "ssyevx.f"
	    i__1 = lwkmin, i__2 = (nb + 3) * *n;
#line 349 "ssyevx.f"
	    lwkopt = max(i__1,i__2);
#line 350 "ssyevx.f"
	    work[1] = (doublereal) lwkopt;
#line 351 "ssyevx.f"
	}

#line 353 "ssyevx.f"
	if (*lwork < lwkmin && ! lquery) {
#line 353 "ssyevx.f"
	    *info = -17;
#line 353 "ssyevx.f"
	}
#line 355 "ssyevx.f"
    }

#line 357 "ssyevx.f"
    if (*info != 0) {
#line 358 "ssyevx.f"
	i__1 = -(*info);
#line 358 "ssyevx.f"
	xerbla_("SSYEVX", &i__1, (ftnlen)6);
#line 359 "ssyevx.f"
	return 0;
#line 360 "ssyevx.f"
    } else if (lquery) {
#line 361 "ssyevx.f"
	return 0;
#line 362 "ssyevx.f"
    }

/*     Quick return if possible */

#line 366 "ssyevx.f"
    *m = 0;
#line 367 "ssyevx.f"
    if (*n == 0) {
#line 368 "ssyevx.f"
	return 0;
#line 369 "ssyevx.f"
    }

#line 371 "ssyevx.f"
    if (*n == 1) {
#line 372 "ssyevx.f"
	if (alleig || indeig) {
#line 373 "ssyevx.f"
	    *m = 1;
#line 374 "ssyevx.f"
	    w[1] = a[a_dim1 + 1];
#line 375 "ssyevx.f"
	} else {
#line 376 "ssyevx.f"
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
#line 377 "ssyevx.f"
		*m = 1;
#line 378 "ssyevx.f"
		w[1] = a[a_dim1 + 1];
#line 379 "ssyevx.f"
	    }
#line 380 "ssyevx.f"
	}
#line 381 "ssyevx.f"
	if (wantz) {
#line 381 "ssyevx.f"
	    z__[z_dim1 + 1] = 1.;
#line 381 "ssyevx.f"
	}
#line 383 "ssyevx.f"
	return 0;
#line 384 "ssyevx.f"
    }

/*     Get machine constants. */

#line 388 "ssyevx.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 389 "ssyevx.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 390 "ssyevx.f"
    smlnum = safmin / eps;
#line 391 "ssyevx.f"
    bignum = 1. / smlnum;
#line 392 "ssyevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 393 "ssyevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 393 "ssyevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 397 "ssyevx.f"
    iscale = 0;
#line 398 "ssyevx.f"
    abstll = *abstol;
#line 399 "ssyevx.f"
    if (valeig) {
#line 400 "ssyevx.f"
	vll = *vl;
#line 401 "ssyevx.f"
	vuu = *vu;
#line 402 "ssyevx.f"
    }
#line 403 "ssyevx.f"
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 404 "ssyevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 405 "ssyevx.f"
	iscale = 1;
#line 406 "ssyevx.f"
	sigma = rmin / anrm;
#line 407 "ssyevx.f"
    } else if (anrm > rmax) {
#line 408 "ssyevx.f"
	iscale = 1;
#line 409 "ssyevx.f"
	sigma = rmax / anrm;
#line 410 "ssyevx.f"
    }
#line 411 "ssyevx.f"
    if (iscale == 1) {
#line 412 "ssyevx.f"
	if (lower) {
#line 413 "ssyevx.f"
	    i__1 = *n;
#line 413 "ssyevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 414 "ssyevx.f"
		i__2 = *n - j + 1;
#line 414 "ssyevx.f"
		sscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 415 "ssyevx.f"
/* L10: */
#line 415 "ssyevx.f"
	    }
#line 416 "ssyevx.f"
	} else {
#line 417 "ssyevx.f"
	    i__1 = *n;
#line 417 "ssyevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 418 "ssyevx.f"
		sscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 419 "ssyevx.f"
/* L20: */
#line 419 "ssyevx.f"
	    }
#line 420 "ssyevx.f"
	}
#line 421 "ssyevx.f"
	if (*abstol > 0.) {
#line 421 "ssyevx.f"
	    abstll = *abstol * sigma;
#line 421 "ssyevx.f"
	}
#line 423 "ssyevx.f"
	if (valeig) {
#line 424 "ssyevx.f"
	    vll = *vl * sigma;
#line 425 "ssyevx.f"
	    vuu = *vu * sigma;
#line 426 "ssyevx.f"
	}
#line 427 "ssyevx.f"
    }

/*     Call SSYTRD to reduce symmetric matrix to tridiagonal form. */

#line 431 "ssyevx.f"
    indtau = 1;
#line 432 "ssyevx.f"
    inde = indtau + *n;
#line 433 "ssyevx.f"
    indd = inde + *n;
#line 434 "ssyevx.f"
    indwrk = indd + *n;
#line 435 "ssyevx.f"
    llwork = *lwork - indwrk + 1;
#line 436 "ssyevx.f"
    ssytrd_(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call SSTERF or SORGTR and SSTEQR.  If this fails for */
/*     some eigenvalue, then try SSTEBZ. */

#line 443 "ssyevx.f"
    test = FALSE_;
#line 444 "ssyevx.f"
    if (indeig) {
#line 445 "ssyevx.f"
	if (*il == 1 && *iu == *n) {
#line 446 "ssyevx.f"
	    test = TRUE_;
#line 447 "ssyevx.f"
	}
#line 448 "ssyevx.f"
    }
#line 449 "ssyevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 450 "ssyevx.f"
	scopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 451 "ssyevx.f"
	indee = indwrk + (*n << 1);
#line 452 "ssyevx.f"
	if (! wantz) {
#line 453 "ssyevx.f"
	    i__1 = *n - 1;
#line 453 "ssyevx.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 454 "ssyevx.f"
	    ssterf_(n, &w[1], &work[indee], info);
#line 455 "ssyevx.f"
	} else {
#line 456 "ssyevx.f"
	    slacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 457 "ssyevx.f"
	    sorgtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
#line 459 "ssyevx.f"
	    i__1 = *n - 1;
#line 459 "ssyevx.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 460 "ssyevx.f"
	    ssteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 462 "ssyevx.f"
	    if (*info == 0) {
#line 463 "ssyevx.f"
		i__1 = *n;
#line 463 "ssyevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 464 "ssyevx.f"
		    ifail[i__] = 0;
#line 465 "ssyevx.f"
/* L30: */
#line 465 "ssyevx.f"
		}
#line 466 "ssyevx.f"
	    }
#line 467 "ssyevx.f"
	}
#line 468 "ssyevx.f"
	if (*info == 0) {
#line 469 "ssyevx.f"
	    *m = *n;
#line 470 "ssyevx.f"
	    goto L40;
#line 471 "ssyevx.f"
	}
#line 472 "ssyevx.f"
	*info = 0;
#line 473 "ssyevx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 477 "ssyevx.f"
    if (wantz) {
#line 478 "ssyevx.f"
	*(unsigned char *)order = 'B';
#line 479 "ssyevx.f"
    } else {
#line 480 "ssyevx.f"
	*(unsigned char *)order = 'E';
#line 481 "ssyevx.f"
    }
#line 482 "ssyevx.f"
    indibl = 1;
#line 483 "ssyevx.f"
    indisp = indibl + *n;
#line 484 "ssyevx.f"
    indiwo = indisp + *n;
#line 485 "ssyevx.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 490 "ssyevx.f"
    if (wantz) {
#line 491 "ssyevx.f"
	sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEIN. */

#line 498 "ssyevx.f"
	indwkn = inde;
#line 499 "ssyevx.f"
	llwrkn = *lwork - indwkn + 1;
#line 500 "ssyevx.f"
	sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 502 "ssyevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 506 "ssyevx.f"
L40:
#line 507 "ssyevx.f"
    if (iscale == 1) {
#line 508 "ssyevx.f"
	if (*info == 0) {
#line 509 "ssyevx.f"
	    imax = *m;
#line 510 "ssyevx.f"
	} else {
#line 511 "ssyevx.f"
	    imax = *info - 1;
#line 512 "ssyevx.f"
	}
#line 513 "ssyevx.f"
	d__1 = 1. / sigma;
#line 513 "ssyevx.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 514 "ssyevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 519 "ssyevx.f"
    if (wantz) {
#line 520 "ssyevx.f"
	i__1 = *m - 1;
#line 520 "ssyevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 521 "ssyevx.f"
	    i__ = 0;
#line 522 "ssyevx.f"
	    tmp1 = w[j];
#line 523 "ssyevx.f"
	    i__2 = *m;
#line 523 "ssyevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 524 "ssyevx.f"
		if (w[jj] < tmp1) {
#line 525 "ssyevx.f"
		    i__ = jj;
#line 526 "ssyevx.f"
		    tmp1 = w[jj];
#line 527 "ssyevx.f"
		}
#line 528 "ssyevx.f"
/* L50: */
#line 528 "ssyevx.f"
	    }

#line 530 "ssyevx.f"
	    if (i__ != 0) {
#line 531 "ssyevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 532 "ssyevx.f"
		w[i__] = w[j];
#line 533 "ssyevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 534 "ssyevx.f"
		w[j] = tmp1;
#line 535 "ssyevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 536 "ssyevx.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 537 "ssyevx.f"
		if (*info != 0) {
#line 538 "ssyevx.f"
		    itmp1 = ifail[i__];
#line 539 "ssyevx.f"
		    ifail[i__] = ifail[j];
#line 540 "ssyevx.f"
		    ifail[j] = itmp1;
#line 541 "ssyevx.f"
		}
#line 542 "ssyevx.f"
	    }
#line 543 "ssyevx.f"
/* L60: */
#line 543 "ssyevx.f"
	}
#line 544 "ssyevx.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 548 "ssyevx.f"
    work[1] = (doublereal) lwkopt;

#line 550 "ssyevx.f"
    return 0;

/*     End of SSYEVX */

} /* ssyevx_ */

