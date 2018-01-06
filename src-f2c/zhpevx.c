#line 1 "zhpevx.f"
/* zhpevx.f -- translated by f2c (version 20100827).
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

#line 1 "zhpevx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> ZHPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, */
/*                          IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDZ, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian matrix A in packed storage. */
/* > Eigenvalues/vectors can be selected by specifying either a range of */
/* > values or a range of indices for the desired eigenvalues. */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
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
/* >          VL is DOUBLE PRECISION */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
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
/* >          If INFO = 0, the selected eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, max(1,M)) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          If an eigenvector fails to converge, then that column of Z */
/* >          contains the latest approximation to the eigenvector, and */
/* >          the index of the eigenvector is returned in IFAIL. */
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

/* > \ingroup complex16OTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int zhpevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, doublereal *
	rwork, integer *iwork, integer *ifail, integer *info, ftnlen jobz_len,
	 ftnlen range_len, ftnlen uplo_len)
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantz;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal abstll, bignum;
    static integer indiwk, indisp, indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dstebz_(char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen);
    static integer indrwk, indwrk, nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int zhptrd_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, ftnlen), 
	    zstein_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, integer *, integer *), zsteqr_(char *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen), zupgtr_(char *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen), zupmtr_(char *, char *, char 
	    *, integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen);


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

#line 292 "zhpevx.f"
    /* Parameter adjustments */
#line 292 "zhpevx.f"
    --ap;
#line 292 "zhpevx.f"
    --w;
#line 292 "zhpevx.f"
    z_dim1 = *ldz;
#line 292 "zhpevx.f"
    z_offset = 1 + z_dim1;
#line 292 "zhpevx.f"
    z__ -= z_offset;
#line 292 "zhpevx.f"
    --work;
#line 292 "zhpevx.f"
    --rwork;
#line 292 "zhpevx.f"
    --iwork;
#line 292 "zhpevx.f"
    --ifail;
#line 292 "zhpevx.f"

#line 292 "zhpevx.f"
    /* Function Body */
#line 292 "zhpevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 293 "zhpevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 294 "zhpevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 295 "zhpevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 297 "zhpevx.f"
    *info = 0;
#line 298 "zhpevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 299 "zhpevx.f"
	*info = -1;
#line 300 "zhpevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 301 "zhpevx.f"
	*info = -2;
#line 302 "zhpevx.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 304 "zhpevx.f"
	*info = -3;
#line 305 "zhpevx.f"
    } else if (*n < 0) {
#line 306 "zhpevx.f"
	*info = -4;
#line 307 "zhpevx.f"
    } else {
#line 308 "zhpevx.f"
	if (valeig) {
#line 309 "zhpevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 309 "zhpevx.f"
		*info = -7;
#line 309 "zhpevx.f"
	    }
#line 311 "zhpevx.f"
	} else if (indeig) {
#line 312 "zhpevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 313 "zhpevx.f"
		*info = -8;
#line 314 "zhpevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 315 "zhpevx.f"
		*info = -9;
#line 316 "zhpevx.f"
	    }
#line 317 "zhpevx.f"
	}
#line 318 "zhpevx.f"
    }
#line 319 "zhpevx.f"
    if (*info == 0) {
#line 320 "zhpevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 320 "zhpevx.f"
	    *info = -14;
#line 320 "zhpevx.f"
	}
#line 322 "zhpevx.f"
    }

#line 324 "zhpevx.f"
    if (*info != 0) {
#line 325 "zhpevx.f"
	i__1 = -(*info);
#line 325 "zhpevx.f"
	xerbla_("ZHPEVX", &i__1, (ftnlen)6);
#line 326 "zhpevx.f"
	return 0;
#line 327 "zhpevx.f"
    }

/*     Quick return if possible */

#line 331 "zhpevx.f"
    *m = 0;
#line 332 "zhpevx.f"
    if (*n == 0) {
#line 332 "zhpevx.f"
	return 0;
#line 332 "zhpevx.f"
    }

#line 335 "zhpevx.f"
    if (*n == 1) {
#line 336 "zhpevx.f"
	if (alleig || indeig) {
#line 337 "zhpevx.f"
	    *m = 1;
#line 338 "zhpevx.f"
	    w[1] = ap[1].r;
#line 339 "zhpevx.f"
	} else {
#line 340 "zhpevx.f"
	    if (*vl < ap[1].r && *vu >= ap[1].r) {
#line 341 "zhpevx.f"
		*m = 1;
#line 342 "zhpevx.f"
		w[1] = ap[1].r;
#line 343 "zhpevx.f"
	    }
#line 344 "zhpevx.f"
	}
#line 345 "zhpevx.f"
	if (wantz) {
#line 345 "zhpevx.f"
	    i__1 = z_dim1 + 1;
#line 345 "zhpevx.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 345 "zhpevx.f"
	}
#line 347 "zhpevx.f"
	return 0;
#line 348 "zhpevx.f"
    }

/*     Get machine constants. */

#line 352 "zhpevx.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 353 "zhpevx.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 354 "zhpevx.f"
    smlnum = safmin / eps;
#line 355 "zhpevx.f"
    bignum = 1. / smlnum;
#line 356 "zhpevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 357 "zhpevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 357 "zhpevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 361 "zhpevx.f"
    iscale = 0;
#line 362 "zhpevx.f"
    abstll = *abstol;
#line 363 "zhpevx.f"
    if (valeig) {
#line 364 "zhpevx.f"
	vll = *vl;
#line 365 "zhpevx.f"
	vuu = *vu;
#line 366 "zhpevx.f"
    } else {
#line 367 "zhpevx.f"
	vll = 0.;
#line 368 "zhpevx.f"
	vuu = 0.;
#line 369 "zhpevx.f"
    }
#line 370 "zhpevx.f"
    anrm = zlanhp_("M", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);
#line 371 "zhpevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 372 "zhpevx.f"
	iscale = 1;
#line 373 "zhpevx.f"
	sigma = rmin / anrm;
#line 374 "zhpevx.f"
    } else if (anrm > rmax) {
#line 375 "zhpevx.f"
	iscale = 1;
#line 376 "zhpevx.f"
	sigma = rmax / anrm;
#line 377 "zhpevx.f"
    }
#line 378 "zhpevx.f"
    if (iscale == 1) {
#line 379 "zhpevx.f"
	i__1 = *n * (*n + 1) / 2;
#line 379 "zhpevx.f"
	zdscal_(&i__1, &sigma, &ap[1], &c__1);
#line 380 "zhpevx.f"
	if (*abstol > 0.) {
#line 380 "zhpevx.f"
	    abstll = *abstol * sigma;
#line 380 "zhpevx.f"
	}
#line 382 "zhpevx.f"
	if (valeig) {
#line 383 "zhpevx.f"
	    vll = *vl * sigma;
#line 384 "zhpevx.f"
	    vuu = *vu * sigma;
#line 385 "zhpevx.f"
	}
#line 386 "zhpevx.f"
    }

/*     Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form. */

#line 390 "zhpevx.f"
    indd = 1;
#line 391 "zhpevx.f"
    inde = indd + *n;
#line 392 "zhpevx.f"
    indrwk = inde + *n;
#line 393 "zhpevx.f"
    indtau = 1;
#line 394 "zhpevx.f"
    indwrk = indtau + *n;
#line 395 "zhpevx.f"
    zhptrd_(uplo, n, &ap[1], &rwork[indd], &rwork[inde], &work[indtau], &
	    iinfo, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call DSTERF or ZUPGTR and ZSTEQR.  If this fails */
/*     for some eigenvalue, then try DSTEBZ. */

#line 402 "zhpevx.f"
    test = FALSE_;
#line 403 "zhpevx.f"
    if (indeig) {
#line 404 "zhpevx.f"
	if (*il == 1 && *iu == *n) {
#line 405 "zhpevx.f"
	    test = TRUE_;
#line 406 "zhpevx.f"
	}
#line 407 "zhpevx.f"
    }
#line 408 "zhpevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 409 "zhpevx.f"
	dcopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 410 "zhpevx.f"
	indee = indrwk + (*n << 1);
#line 411 "zhpevx.f"
	if (! wantz) {
#line 412 "zhpevx.f"
	    i__1 = *n - 1;
#line 412 "zhpevx.f"
	    dcopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 413 "zhpevx.f"
	    dsterf_(n, &w[1], &rwork[indee], info);
#line 414 "zhpevx.f"
	} else {
#line 415 "zhpevx.f"
	    zupgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &
		    work[indwrk], &iinfo, (ftnlen)1);
#line 417 "zhpevx.f"
	    i__1 = *n - 1;
#line 417 "zhpevx.f"
	    dcopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 418 "zhpevx.f"
	    zsteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 420 "zhpevx.f"
	    if (*info == 0) {
#line 421 "zhpevx.f"
		i__1 = *n;
#line 421 "zhpevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 422 "zhpevx.f"
		    ifail[i__] = 0;
#line 423 "zhpevx.f"
/* L10: */
#line 423 "zhpevx.f"
		}
#line 424 "zhpevx.f"
	    }
#line 425 "zhpevx.f"
	}
#line 426 "zhpevx.f"
	if (*info == 0) {
#line 427 "zhpevx.f"
	    *m = *n;
#line 428 "zhpevx.f"
	    goto L20;
#line 429 "zhpevx.f"
	}
#line 430 "zhpevx.f"
	*info = 0;
#line 431 "zhpevx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN. */

#line 435 "zhpevx.f"
    if (wantz) {
#line 436 "zhpevx.f"
	*(unsigned char *)order = 'B';
#line 437 "zhpevx.f"
    } else {
#line 438 "zhpevx.f"
	*(unsigned char *)order = 'E';
#line 439 "zhpevx.f"
    }
#line 440 "zhpevx.f"
    indibl = 1;
#line 441 "zhpevx.f"
    indisp = indibl + *n;
#line 442 "zhpevx.f"
    indiwk = indisp + *n;
#line 443 "zhpevx.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 448 "zhpevx.f"
    if (wantz) {
#line 449 "zhpevx.f"
	zstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by ZSTEIN. */

#line 456 "zhpevx.f"
	indwrk = indtau + *n;
#line 457 "zhpevx.f"
	zupmtr_("L", uplo, "N", n, m, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 459 "zhpevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 463 "zhpevx.f"
L20:
#line 464 "zhpevx.f"
    if (iscale == 1) {
#line 465 "zhpevx.f"
	if (*info == 0) {
#line 466 "zhpevx.f"
	    imax = *m;
#line 467 "zhpevx.f"
	} else {
#line 468 "zhpevx.f"
	    imax = *info - 1;
#line 469 "zhpevx.f"
	}
#line 470 "zhpevx.f"
	d__1 = 1. / sigma;
#line 470 "zhpevx.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 471 "zhpevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 476 "zhpevx.f"
    if (wantz) {
#line 477 "zhpevx.f"
	i__1 = *m - 1;
#line 477 "zhpevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 478 "zhpevx.f"
	    i__ = 0;
#line 479 "zhpevx.f"
	    tmp1 = w[j];
#line 480 "zhpevx.f"
	    i__2 = *m;
#line 480 "zhpevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 481 "zhpevx.f"
		if (w[jj] < tmp1) {
#line 482 "zhpevx.f"
		    i__ = jj;
#line 483 "zhpevx.f"
		    tmp1 = w[jj];
#line 484 "zhpevx.f"
		}
#line 485 "zhpevx.f"
/* L30: */
#line 485 "zhpevx.f"
	    }

#line 487 "zhpevx.f"
	    if (i__ != 0) {
#line 488 "zhpevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 489 "zhpevx.f"
		w[i__] = w[j];
#line 490 "zhpevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 491 "zhpevx.f"
		w[j] = tmp1;
#line 492 "zhpevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 493 "zhpevx.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 494 "zhpevx.f"
		if (*info != 0) {
#line 495 "zhpevx.f"
		    itmp1 = ifail[i__];
#line 496 "zhpevx.f"
		    ifail[i__] = ifail[j];
#line 497 "zhpevx.f"
		    ifail[j] = itmp1;
#line 498 "zhpevx.f"
		}
#line 499 "zhpevx.f"
	    }
#line 500 "zhpevx.f"
/* L40: */
#line 500 "zhpevx.f"
	}
#line 501 "zhpevx.f"
    }

#line 503 "zhpevx.f"
    return 0;

/*     End of ZHPEVX */

} /* zhpevx_ */

