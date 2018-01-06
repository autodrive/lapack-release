#line 1 "dsyevx.f"
/* dsyevx.f -- translated by f2c (version 20100827).
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

#line 1 "dsyevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> DSYEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, */
/*                          IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYEVX computes selected eigenvalues and, optionally, eigenvectors */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
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
/* >          by reducing A to tridiagonal form. */
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
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= 1, when N <= 1; */
/* >          otherwise 8*N. */
/* >          For optimal efficiency, LWORK >= (NB+3)*N, */
/* >          where NB is the max of the blocksize for DSYTRD and DORMTR */
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

/* > \ingroup doubleSYeigen */

/*  ===================================================================== */
/* Subroutine */ int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lower, wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp;
    extern /* Subroutine */ int dstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dsterf_(integer *, doublereal *, doublereal *, integer *);
    static integer indiwo, indwkn;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indwrk, lwkmin;
    extern /* Subroutine */ int dorgtr_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), dsteqr_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dormtr_(char *, char *, char *, integer *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer llwrkn, llwork, nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int dsytrd_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;


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

#line 298 "dsyevx.f"
    /* Parameter adjustments */
#line 298 "dsyevx.f"
    a_dim1 = *lda;
#line 298 "dsyevx.f"
    a_offset = 1 + a_dim1;
#line 298 "dsyevx.f"
    a -= a_offset;
#line 298 "dsyevx.f"
    --w;
#line 298 "dsyevx.f"
    z_dim1 = *ldz;
#line 298 "dsyevx.f"
    z_offset = 1 + z_dim1;
#line 298 "dsyevx.f"
    z__ -= z_offset;
#line 298 "dsyevx.f"
    --work;
#line 298 "dsyevx.f"
    --iwork;
#line 298 "dsyevx.f"
    --ifail;
#line 298 "dsyevx.f"

#line 298 "dsyevx.f"
    /* Function Body */
#line 298 "dsyevx.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 299 "dsyevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 300 "dsyevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 301 "dsyevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 302 "dsyevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 303 "dsyevx.f"
    lquery = *lwork == -1;

#line 305 "dsyevx.f"
    *info = 0;
#line 306 "dsyevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 307 "dsyevx.f"
	*info = -1;
#line 308 "dsyevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 309 "dsyevx.f"
	*info = -2;
#line 310 "dsyevx.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 311 "dsyevx.f"
	*info = -3;
#line 312 "dsyevx.f"
    } else if (*n < 0) {
#line 313 "dsyevx.f"
	*info = -4;
#line 314 "dsyevx.f"
    } else if (*lda < max(1,*n)) {
#line 315 "dsyevx.f"
	*info = -6;
#line 316 "dsyevx.f"
    } else {
#line 317 "dsyevx.f"
	if (valeig) {
#line 318 "dsyevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 318 "dsyevx.f"
		*info = -8;
#line 318 "dsyevx.f"
	    }
#line 320 "dsyevx.f"
	} else if (indeig) {
#line 321 "dsyevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 322 "dsyevx.f"
		*info = -9;
#line 323 "dsyevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 324 "dsyevx.f"
		*info = -10;
#line 325 "dsyevx.f"
	    }
#line 326 "dsyevx.f"
	}
#line 327 "dsyevx.f"
    }
#line 328 "dsyevx.f"
    if (*info == 0) {
#line 329 "dsyevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 330 "dsyevx.f"
	    *info = -15;
#line 331 "dsyevx.f"
	}
#line 332 "dsyevx.f"
    }

#line 334 "dsyevx.f"
    if (*info == 0) {
#line 335 "dsyevx.f"
	if (*n <= 1) {
#line 336 "dsyevx.f"
	    lwkmin = 1;
#line 337 "dsyevx.f"
	    work[1] = (doublereal) lwkmin;
#line 338 "dsyevx.f"
	} else {
#line 339 "dsyevx.f"
	    lwkmin = *n << 3;
#line 340 "dsyevx.f"
	    nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 341 "dsyevx.f"
	    i__1 = nb, i__2 = ilaenv_(&c__1, "DORMTR", uplo, n, &c_n1, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1);
#line 341 "dsyevx.f"
	    nb = max(i__1,i__2);
/* Computing MAX */
#line 342 "dsyevx.f"
	    i__1 = lwkmin, i__2 = (nb + 3) * *n;
#line 342 "dsyevx.f"
	    lwkopt = max(i__1,i__2);
#line 343 "dsyevx.f"
	    work[1] = (doublereal) lwkopt;
#line 344 "dsyevx.f"
	}

#line 346 "dsyevx.f"
	if (*lwork < lwkmin && ! lquery) {
#line 346 "dsyevx.f"
	    *info = -17;
#line 346 "dsyevx.f"
	}
#line 348 "dsyevx.f"
    }

#line 350 "dsyevx.f"
    if (*info != 0) {
#line 351 "dsyevx.f"
	i__1 = -(*info);
#line 351 "dsyevx.f"
	xerbla_("DSYEVX", &i__1, (ftnlen)6);
#line 352 "dsyevx.f"
	return 0;
#line 353 "dsyevx.f"
    } else if (lquery) {
#line 354 "dsyevx.f"
	return 0;
#line 355 "dsyevx.f"
    }

/*     Quick return if possible */

#line 359 "dsyevx.f"
    *m = 0;
#line 360 "dsyevx.f"
    if (*n == 0) {
#line 361 "dsyevx.f"
	return 0;
#line 362 "dsyevx.f"
    }

#line 364 "dsyevx.f"
    if (*n == 1) {
#line 365 "dsyevx.f"
	if (alleig || indeig) {
#line 366 "dsyevx.f"
	    *m = 1;
#line 367 "dsyevx.f"
	    w[1] = a[a_dim1 + 1];
#line 368 "dsyevx.f"
	} else {
#line 369 "dsyevx.f"
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
#line 370 "dsyevx.f"
		*m = 1;
#line 371 "dsyevx.f"
		w[1] = a[a_dim1 + 1];
#line 372 "dsyevx.f"
	    }
#line 373 "dsyevx.f"
	}
#line 374 "dsyevx.f"
	if (wantz) {
#line 374 "dsyevx.f"
	    z__[z_dim1 + 1] = 1.;
#line 374 "dsyevx.f"
	}
#line 376 "dsyevx.f"
	return 0;
#line 377 "dsyevx.f"
    }

/*     Get machine constants. */

#line 381 "dsyevx.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 382 "dsyevx.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 383 "dsyevx.f"
    smlnum = safmin / eps;
#line 384 "dsyevx.f"
    bignum = 1. / smlnum;
#line 385 "dsyevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 386 "dsyevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 386 "dsyevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 390 "dsyevx.f"
    iscale = 0;
#line 391 "dsyevx.f"
    abstll = *abstol;
#line 392 "dsyevx.f"
    if (valeig) {
#line 393 "dsyevx.f"
	vll = *vl;
#line 394 "dsyevx.f"
	vuu = *vu;
#line 395 "dsyevx.f"
    }
#line 396 "dsyevx.f"
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 397 "dsyevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 398 "dsyevx.f"
	iscale = 1;
#line 399 "dsyevx.f"
	sigma = rmin / anrm;
#line 400 "dsyevx.f"
    } else if (anrm > rmax) {
#line 401 "dsyevx.f"
	iscale = 1;
#line 402 "dsyevx.f"
	sigma = rmax / anrm;
#line 403 "dsyevx.f"
    }
#line 404 "dsyevx.f"
    if (iscale == 1) {
#line 405 "dsyevx.f"
	if (lower) {
#line 406 "dsyevx.f"
	    i__1 = *n;
#line 406 "dsyevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 407 "dsyevx.f"
		i__2 = *n - j + 1;
#line 407 "dsyevx.f"
		dscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 408 "dsyevx.f"
/* L10: */
#line 408 "dsyevx.f"
	    }
#line 409 "dsyevx.f"
	} else {
#line 410 "dsyevx.f"
	    i__1 = *n;
#line 410 "dsyevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 411 "dsyevx.f"
		dscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 412 "dsyevx.f"
/* L20: */
#line 412 "dsyevx.f"
	    }
#line 413 "dsyevx.f"
	}
#line 414 "dsyevx.f"
	if (*abstol > 0.) {
#line 414 "dsyevx.f"
	    abstll = *abstol * sigma;
#line 414 "dsyevx.f"
	}
#line 416 "dsyevx.f"
	if (valeig) {
#line 417 "dsyevx.f"
	    vll = *vl * sigma;
#line 418 "dsyevx.f"
	    vuu = *vu * sigma;
#line 419 "dsyevx.f"
	}
#line 420 "dsyevx.f"
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

#line 424 "dsyevx.f"
    indtau = 1;
#line 425 "dsyevx.f"
    inde = indtau + *n;
#line 426 "dsyevx.f"
    indd = inde + *n;
#line 427 "dsyevx.f"
    indwrk = indd + *n;
#line 428 "dsyevx.f"
    llwork = *lwork - indwrk + 1;
#line 429 "dsyevx.f"
    dsytrd_(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call DSTERF or DORGTR and SSTEQR.  If this fails for */
/*     some eigenvalue, then try DSTEBZ. */

#line 436 "dsyevx.f"
    test = FALSE_;
#line 437 "dsyevx.f"
    if (indeig) {
#line 438 "dsyevx.f"
	if (*il == 1 && *iu == *n) {
#line 439 "dsyevx.f"
	    test = TRUE_;
#line 440 "dsyevx.f"
	}
#line 441 "dsyevx.f"
    }
#line 442 "dsyevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 443 "dsyevx.f"
	dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 444 "dsyevx.f"
	indee = indwrk + (*n << 1);
#line 445 "dsyevx.f"
	if (! wantz) {
#line 446 "dsyevx.f"
	    i__1 = *n - 1;
#line 446 "dsyevx.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 447 "dsyevx.f"
	    dsterf_(n, &w[1], &work[indee], info);
#line 448 "dsyevx.f"
	} else {
#line 449 "dsyevx.f"
	    dlacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 450 "dsyevx.f"
	    dorgtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
#line 452 "dsyevx.f"
	    i__1 = *n - 1;
#line 452 "dsyevx.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 453 "dsyevx.f"
	    dsteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 455 "dsyevx.f"
	    if (*info == 0) {
#line 456 "dsyevx.f"
		i__1 = *n;
#line 456 "dsyevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 457 "dsyevx.f"
		    ifail[i__] = 0;
#line 458 "dsyevx.f"
/* L30: */
#line 458 "dsyevx.f"
		}
#line 459 "dsyevx.f"
	    }
#line 460 "dsyevx.f"
	}
#line 461 "dsyevx.f"
	if (*info == 0) {
#line 462 "dsyevx.f"
	    *m = *n;
#line 463 "dsyevx.f"
	    goto L40;
#line 464 "dsyevx.f"
	}
#line 465 "dsyevx.f"
	*info = 0;
#line 466 "dsyevx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 470 "dsyevx.f"
    if (wantz) {
#line 471 "dsyevx.f"
	*(unsigned char *)order = 'B';
#line 472 "dsyevx.f"
    } else {
#line 473 "dsyevx.f"
	*(unsigned char *)order = 'E';
#line 474 "dsyevx.f"
    }
#line 475 "dsyevx.f"
    indibl = 1;
#line 476 "dsyevx.f"
    indisp = indibl + *n;
#line 477 "dsyevx.f"
    indiwo = indisp + *n;
#line 478 "dsyevx.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 483 "dsyevx.f"
    if (wantz) {
#line 484 "dsyevx.f"
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

#line 491 "dsyevx.f"
	indwkn = inde;
#line 492 "dsyevx.f"
	llwrkn = *lwork - indwkn + 1;
#line 493 "dsyevx.f"
	dormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 495 "dsyevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 499 "dsyevx.f"
L40:
#line 500 "dsyevx.f"
    if (iscale == 1) {
#line 501 "dsyevx.f"
	if (*info == 0) {
#line 502 "dsyevx.f"
	    imax = *m;
#line 503 "dsyevx.f"
	} else {
#line 504 "dsyevx.f"
	    imax = *info - 1;
#line 505 "dsyevx.f"
	}
#line 506 "dsyevx.f"
	d__1 = 1. / sigma;
#line 506 "dsyevx.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 507 "dsyevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 512 "dsyevx.f"
    if (wantz) {
#line 513 "dsyevx.f"
	i__1 = *m - 1;
#line 513 "dsyevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 514 "dsyevx.f"
	    i__ = 0;
#line 515 "dsyevx.f"
	    tmp1 = w[j];
#line 516 "dsyevx.f"
	    i__2 = *m;
#line 516 "dsyevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 517 "dsyevx.f"
		if (w[jj] < tmp1) {
#line 518 "dsyevx.f"
		    i__ = jj;
#line 519 "dsyevx.f"
		    tmp1 = w[jj];
#line 520 "dsyevx.f"
		}
#line 521 "dsyevx.f"
/* L50: */
#line 521 "dsyevx.f"
	    }

#line 523 "dsyevx.f"
	    if (i__ != 0) {
#line 524 "dsyevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 525 "dsyevx.f"
		w[i__] = w[j];
#line 526 "dsyevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 527 "dsyevx.f"
		w[j] = tmp1;
#line 528 "dsyevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 529 "dsyevx.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 530 "dsyevx.f"
		if (*info != 0) {
#line 531 "dsyevx.f"
		    itmp1 = ifail[i__];
#line 532 "dsyevx.f"
		    ifail[i__] = ifail[j];
#line 533 "dsyevx.f"
		    ifail[j] = itmp1;
#line 534 "dsyevx.f"
		}
#line 535 "dsyevx.f"
	    }
#line 536 "dsyevx.f"
/* L60: */
#line 536 "dsyevx.f"
	}
#line 537 "dsyevx.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 541 "dsyevx.f"
    work[1] = (doublereal) lwkopt;

#line 543 "dsyevx.f"
    return 0;

/*     End of DSYEVX */

} /* dsyevx_ */

