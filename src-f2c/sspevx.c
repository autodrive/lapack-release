#line 1 "sspevx.f"
/* sspevx.f -- translated by f2c (version 20100827).
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

#line 1 "sspevx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SSPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDZ, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A in packed storage.  Eigenvalues/vectors */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, AP is overwritten by values generated during the */
/* >          reduction to tridiagonal form.  If UPLO = 'U', the diagonal */
/* >          and first superdiagonal of the tridiagonal matrix T overwrite */
/* >          the corresponding elements of A, and if UPLO = 'L', the */
/* >          diagonal and first subdiagonal of T overwrite the */
/* >          corresponding elements of A. */
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
/* >          by reducing AP to tridiagonal form. */
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
/* >          If INFO = 0, the selected eigenvalues in ascending order. */
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

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int sspevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *ap, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	integer *info, ftnlen jobz_len, ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
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
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical wantz, alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp, indiwo, indwrk;
    extern doublereal slansp_(char *, char *, integer *, doublereal *, 
	    doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int sstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    ssterf_(integer *, doublereal *, doublereal *, integer *);
    static integer nsplit;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int sopgtr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), ssptrd_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), ssteqr_(char *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), sopmtr_(char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);


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

#line 283 "sspevx.f"
    /* Parameter adjustments */
#line 283 "sspevx.f"
    --ap;
#line 283 "sspevx.f"
    --w;
#line 283 "sspevx.f"
    z_dim1 = *ldz;
#line 283 "sspevx.f"
    z_offset = 1 + z_dim1;
#line 283 "sspevx.f"
    z__ -= z_offset;
#line 283 "sspevx.f"
    --work;
#line 283 "sspevx.f"
    --iwork;
#line 283 "sspevx.f"
    --ifail;
#line 283 "sspevx.f"

#line 283 "sspevx.f"
    /* Function Body */
#line 283 "sspevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 284 "sspevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 285 "sspevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 286 "sspevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 288 "sspevx.f"
    *info = 0;
#line 289 "sspevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 290 "sspevx.f"
	*info = -1;
#line 291 "sspevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 292 "sspevx.f"
	*info = -2;
#line 293 "sspevx.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 295 "sspevx.f"
	*info = -3;
#line 296 "sspevx.f"
    } else if (*n < 0) {
#line 297 "sspevx.f"
	*info = -4;
#line 298 "sspevx.f"
    } else {
#line 299 "sspevx.f"
	if (valeig) {
#line 300 "sspevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 300 "sspevx.f"
		*info = -7;
#line 300 "sspevx.f"
	    }
#line 302 "sspevx.f"
	} else if (indeig) {
#line 303 "sspevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 304 "sspevx.f"
		*info = -8;
#line 305 "sspevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 306 "sspevx.f"
		*info = -9;
#line 307 "sspevx.f"
	    }
#line 308 "sspevx.f"
	}
#line 309 "sspevx.f"
    }
#line 310 "sspevx.f"
    if (*info == 0) {
#line 311 "sspevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 311 "sspevx.f"
	    *info = -14;
#line 311 "sspevx.f"
	}
#line 313 "sspevx.f"
    }

#line 315 "sspevx.f"
    if (*info != 0) {
#line 316 "sspevx.f"
	i__1 = -(*info);
#line 316 "sspevx.f"
	xerbla_("SSPEVX", &i__1, (ftnlen)6);
#line 317 "sspevx.f"
	return 0;
#line 318 "sspevx.f"
    }

/*     Quick return if possible */

#line 322 "sspevx.f"
    *m = 0;
#line 323 "sspevx.f"
    if (*n == 0) {
#line 323 "sspevx.f"
	return 0;
#line 323 "sspevx.f"
    }

#line 326 "sspevx.f"
    if (*n == 1) {
#line 327 "sspevx.f"
	if (alleig || indeig) {
#line 328 "sspevx.f"
	    *m = 1;
#line 329 "sspevx.f"
	    w[1] = ap[1];
#line 330 "sspevx.f"
	} else {
#line 331 "sspevx.f"
	    if (*vl < ap[1] && *vu >= ap[1]) {
#line 332 "sspevx.f"
		*m = 1;
#line 333 "sspevx.f"
		w[1] = ap[1];
#line 334 "sspevx.f"
	    }
#line 335 "sspevx.f"
	}
#line 336 "sspevx.f"
	if (wantz) {
#line 336 "sspevx.f"
	    z__[z_dim1 + 1] = 1.;
#line 336 "sspevx.f"
	}
#line 338 "sspevx.f"
	return 0;
#line 339 "sspevx.f"
    }

/*     Get machine constants. */

#line 343 "sspevx.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 344 "sspevx.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 345 "sspevx.f"
    smlnum = safmin / eps;
#line 346 "sspevx.f"
    bignum = 1. / smlnum;
#line 347 "sspevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 348 "sspevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 348 "sspevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 352 "sspevx.f"
    iscale = 0;
#line 353 "sspevx.f"
    abstll = *abstol;
#line 354 "sspevx.f"
    if (valeig) {
#line 355 "sspevx.f"
	vll = *vl;
#line 356 "sspevx.f"
	vuu = *vu;
#line 357 "sspevx.f"
    } else {
#line 358 "sspevx.f"
	vll = 0.;
#line 359 "sspevx.f"
	vuu = 0.;
#line 360 "sspevx.f"
    }
#line 361 "sspevx.f"
    anrm = slansp_("M", uplo, n, &ap[1], &work[1], (ftnlen)1, (ftnlen)1);
#line 362 "sspevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 363 "sspevx.f"
	iscale = 1;
#line 364 "sspevx.f"
	sigma = rmin / anrm;
#line 365 "sspevx.f"
    } else if (anrm > rmax) {
#line 366 "sspevx.f"
	iscale = 1;
#line 367 "sspevx.f"
	sigma = rmax / anrm;
#line 368 "sspevx.f"
    }
#line 369 "sspevx.f"
    if (iscale == 1) {
#line 370 "sspevx.f"
	i__1 = *n * (*n + 1) / 2;
#line 370 "sspevx.f"
	sscal_(&i__1, &sigma, &ap[1], &c__1);
#line 371 "sspevx.f"
	if (*abstol > 0.) {
#line 371 "sspevx.f"
	    abstll = *abstol * sigma;
#line 371 "sspevx.f"
	}
#line 373 "sspevx.f"
	if (valeig) {
#line 374 "sspevx.f"
	    vll = *vl * sigma;
#line 375 "sspevx.f"
	    vuu = *vu * sigma;
#line 376 "sspevx.f"
	}
#line 377 "sspevx.f"
    }

/*     Call SSPTRD to reduce symmetric packed matrix to tridiagonal form. */

#line 381 "sspevx.f"
    indtau = 1;
#line 382 "sspevx.f"
    inde = indtau + *n;
#line 383 "sspevx.f"
    indd = inde + *n;
#line 384 "sspevx.f"
    indwrk = indd + *n;
#line 385 "sspevx.f"
    ssptrd_(uplo, n, &ap[1], &work[indd], &work[inde], &work[indtau], &iinfo, 
	    (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call SSTERF or SOPGTR and SSTEQR.  If this fails */
/*     for some eigenvalue, then try SSTEBZ. */

#line 392 "sspevx.f"
    test = FALSE_;
#line 393 "sspevx.f"
    if (indeig) {
#line 394 "sspevx.f"
	if (*il == 1 && *iu == *n) {
#line 395 "sspevx.f"
	    test = TRUE_;
#line 396 "sspevx.f"
	}
#line 397 "sspevx.f"
    }
#line 398 "sspevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 399 "sspevx.f"
	scopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 400 "sspevx.f"
	indee = indwrk + (*n << 1);
#line 401 "sspevx.f"
	if (! wantz) {
#line 402 "sspevx.f"
	    i__1 = *n - 1;
#line 402 "sspevx.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 403 "sspevx.f"
	    ssterf_(n, &w[1], &work[indee], info);
#line 404 "sspevx.f"
	} else {
#line 405 "sspevx.f"
	    sopgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &
		    work[indwrk], &iinfo, (ftnlen)1);
#line 407 "sspevx.f"
	    i__1 = *n - 1;
#line 407 "sspevx.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 408 "sspevx.f"
	    ssteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 410 "sspevx.f"
	    if (*info == 0) {
#line 411 "sspevx.f"
		i__1 = *n;
#line 411 "sspevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 412 "sspevx.f"
		    ifail[i__] = 0;
#line 413 "sspevx.f"
/* L10: */
#line 413 "sspevx.f"
		}
#line 414 "sspevx.f"
	    }
#line 415 "sspevx.f"
	}
#line 416 "sspevx.f"
	if (*info == 0) {
#line 417 "sspevx.f"
	    *m = *n;
#line 418 "sspevx.f"
	    goto L20;
#line 419 "sspevx.f"
	}
#line 420 "sspevx.f"
	*info = 0;
#line 421 "sspevx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 425 "sspevx.f"
    if (wantz) {
#line 426 "sspevx.f"
	*(unsigned char *)order = 'B';
#line 427 "sspevx.f"
    } else {
#line 428 "sspevx.f"
	*(unsigned char *)order = 'E';
#line 429 "sspevx.f"
    }
#line 430 "sspevx.f"
    indibl = 1;
#line 431 "sspevx.f"
    indisp = indibl + *n;
#line 432 "sspevx.f"
    indiwo = indisp + *n;
#line 433 "sspevx.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 438 "sspevx.f"
    if (wantz) {
#line 439 "sspevx.f"
	sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEIN. */

#line 446 "sspevx.f"
	sopmtr_("L", uplo, "N", n, m, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 448 "sspevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 452 "sspevx.f"
L20:
#line 453 "sspevx.f"
    if (iscale == 1) {
#line 454 "sspevx.f"
	if (*info == 0) {
#line 455 "sspevx.f"
	    imax = *m;
#line 456 "sspevx.f"
	} else {
#line 457 "sspevx.f"
	    imax = *info - 1;
#line 458 "sspevx.f"
	}
#line 459 "sspevx.f"
	d__1 = 1. / sigma;
#line 459 "sspevx.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 460 "sspevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 465 "sspevx.f"
    if (wantz) {
#line 466 "sspevx.f"
	i__1 = *m - 1;
#line 466 "sspevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 467 "sspevx.f"
	    i__ = 0;
#line 468 "sspevx.f"
	    tmp1 = w[j];
#line 469 "sspevx.f"
	    i__2 = *m;
#line 469 "sspevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 470 "sspevx.f"
		if (w[jj] < tmp1) {
#line 471 "sspevx.f"
		    i__ = jj;
#line 472 "sspevx.f"
		    tmp1 = w[jj];
#line 473 "sspevx.f"
		}
#line 474 "sspevx.f"
/* L30: */
#line 474 "sspevx.f"
	    }

#line 476 "sspevx.f"
	    if (i__ != 0) {
#line 477 "sspevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 478 "sspevx.f"
		w[i__] = w[j];
#line 479 "sspevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 480 "sspevx.f"
		w[j] = tmp1;
#line 481 "sspevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 482 "sspevx.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 483 "sspevx.f"
		if (*info != 0) {
#line 484 "sspevx.f"
		    itmp1 = ifail[i__];
#line 485 "sspevx.f"
		    ifail[i__] = ifail[j];
#line 486 "sspevx.f"
		    ifail[j] = itmp1;
#line 487 "sspevx.f"
		}
#line 488 "sspevx.f"
	    }
#line 489 "sspevx.f"
/* L40: */
#line 489 "sspevx.f"
	}
#line 490 "sspevx.f"
    }

#line 492 "sspevx.f"
    return 0;

/*     End of SSPEVX */

} /* sspevx_ */

