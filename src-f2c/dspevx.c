#line 1 "dspevx.f"
/* dspevx.f -- translated by f2c (version 20100827).
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

#line 1 "dspevx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DSPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDZ, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPEVX computes selected eigenvalues and, optionally, eigenvectors */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (8*N) */
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

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dspevx_(char *jobz, char *range, char *uplo, integer *n, 
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    extern doublereal dlansp_(char *, char *, integer *, doublereal *, 
	    doublereal *, ftnlen, ftnlen);
    static integer indtau, indisp;
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
    extern /* Subroutine */ int dopgtr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dsptrd_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), dsteqr_(char *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dopmtr_(char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static integer nsplit;
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

#line 283 "dspevx.f"
    /* Parameter adjustments */
#line 283 "dspevx.f"
    --ap;
#line 283 "dspevx.f"
    --w;
#line 283 "dspevx.f"
    z_dim1 = *ldz;
#line 283 "dspevx.f"
    z_offset = 1 + z_dim1;
#line 283 "dspevx.f"
    z__ -= z_offset;
#line 283 "dspevx.f"
    --work;
#line 283 "dspevx.f"
    --iwork;
#line 283 "dspevx.f"
    --ifail;
#line 283 "dspevx.f"

#line 283 "dspevx.f"
    /* Function Body */
#line 283 "dspevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 284 "dspevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 285 "dspevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 286 "dspevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 288 "dspevx.f"
    *info = 0;
#line 289 "dspevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 290 "dspevx.f"
	*info = -1;
#line 291 "dspevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 292 "dspevx.f"
	*info = -2;
#line 293 "dspevx.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 295 "dspevx.f"
	*info = -3;
#line 296 "dspevx.f"
    } else if (*n < 0) {
#line 297 "dspevx.f"
	*info = -4;
#line 298 "dspevx.f"
    } else {
#line 299 "dspevx.f"
	if (valeig) {
#line 300 "dspevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 300 "dspevx.f"
		*info = -7;
#line 300 "dspevx.f"
	    }
#line 302 "dspevx.f"
	} else if (indeig) {
#line 303 "dspevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 304 "dspevx.f"
		*info = -8;
#line 305 "dspevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 306 "dspevx.f"
		*info = -9;
#line 307 "dspevx.f"
	    }
#line 308 "dspevx.f"
	}
#line 309 "dspevx.f"
    }
#line 310 "dspevx.f"
    if (*info == 0) {
#line 311 "dspevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 311 "dspevx.f"
	    *info = -14;
#line 311 "dspevx.f"
	}
#line 313 "dspevx.f"
    }

#line 315 "dspevx.f"
    if (*info != 0) {
#line 316 "dspevx.f"
	i__1 = -(*info);
#line 316 "dspevx.f"
	xerbla_("DSPEVX", &i__1, (ftnlen)6);
#line 317 "dspevx.f"
	return 0;
#line 318 "dspevx.f"
    }

/*     Quick return if possible */

#line 322 "dspevx.f"
    *m = 0;
#line 323 "dspevx.f"
    if (*n == 0) {
#line 323 "dspevx.f"
	return 0;
#line 323 "dspevx.f"
    }

#line 326 "dspevx.f"
    if (*n == 1) {
#line 327 "dspevx.f"
	if (alleig || indeig) {
#line 328 "dspevx.f"
	    *m = 1;
#line 329 "dspevx.f"
	    w[1] = ap[1];
#line 330 "dspevx.f"
	} else {
#line 331 "dspevx.f"
	    if (*vl < ap[1] && *vu >= ap[1]) {
#line 332 "dspevx.f"
		*m = 1;
#line 333 "dspevx.f"
		w[1] = ap[1];
#line 334 "dspevx.f"
	    }
#line 335 "dspevx.f"
	}
#line 336 "dspevx.f"
	if (wantz) {
#line 336 "dspevx.f"
	    z__[z_dim1 + 1] = 1.;
#line 336 "dspevx.f"
	}
#line 338 "dspevx.f"
	return 0;
#line 339 "dspevx.f"
    }

/*     Get machine constants. */

#line 343 "dspevx.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 344 "dspevx.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 345 "dspevx.f"
    smlnum = safmin / eps;
#line 346 "dspevx.f"
    bignum = 1. / smlnum;
#line 347 "dspevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 348 "dspevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 348 "dspevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 352 "dspevx.f"
    iscale = 0;
#line 353 "dspevx.f"
    abstll = *abstol;
#line 354 "dspevx.f"
    if (valeig) {
#line 355 "dspevx.f"
	vll = *vl;
#line 356 "dspevx.f"
	vuu = *vu;
#line 357 "dspevx.f"
    } else {
#line 358 "dspevx.f"
	vll = 0.;
#line 359 "dspevx.f"
	vuu = 0.;
#line 360 "dspevx.f"
    }
#line 361 "dspevx.f"
    anrm = dlansp_("M", uplo, n, &ap[1], &work[1], (ftnlen)1, (ftnlen)1);
#line 362 "dspevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 363 "dspevx.f"
	iscale = 1;
#line 364 "dspevx.f"
	sigma = rmin / anrm;
#line 365 "dspevx.f"
    } else if (anrm > rmax) {
#line 366 "dspevx.f"
	iscale = 1;
#line 367 "dspevx.f"
	sigma = rmax / anrm;
#line 368 "dspevx.f"
    }
#line 369 "dspevx.f"
    if (iscale == 1) {
#line 370 "dspevx.f"
	i__1 = *n * (*n + 1) / 2;
#line 370 "dspevx.f"
	dscal_(&i__1, &sigma, &ap[1], &c__1);
#line 371 "dspevx.f"
	if (*abstol > 0.) {
#line 371 "dspevx.f"
	    abstll = *abstol * sigma;
#line 371 "dspevx.f"
	}
#line 373 "dspevx.f"
	if (valeig) {
#line 374 "dspevx.f"
	    vll = *vl * sigma;
#line 375 "dspevx.f"
	    vuu = *vu * sigma;
#line 376 "dspevx.f"
	}
#line 377 "dspevx.f"
    }

/*     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form. */

#line 381 "dspevx.f"
    indtau = 1;
#line 382 "dspevx.f"
    inde = indtau + *n;
#line 383 "dspevx.f"
    indd = inde + *n;
#line 384 "dspevx.f"
    indwrk = indd + *n;
#line 385 "dspevx.f"
    dsptrd_(uplo, n, &ap[1], &work[indd], &work[inde], &work[indtau], &iinfo, 
	    (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call DSTERF or DOPGTR and SSTEQR.  If this fails */
/*     for some eigenvalue, then try DSTEBZ. */

#line 392 "dspevx.f"
    test = FALSE_;
#line 393 "dspevx.f"
    if (indeig) {
#line 394 "dspevx.f"
	if (*il == 1 && *iu == *n) {
#line 395 "dspevx.f"
	    test = TRUE_;
#line 396 "dspevx.f"
	}
#line 397 "dspevx.f"
    }
#line 398 "dspevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 399 "dspevx.f"
	dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 400 "dspevx.f"
	indee = indwrk + (*n << 1);
#line 401 "dspevx.f"
	if (! wantz) {
#line 402 "dspevx.f"
	    i__1 = *n - 1;
#line 402 "dspevx.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 403 "dspevx.f"
	    dsterf_(n, &w[1], &work[indee], info);
#line 404 "dspevx.f"
	} else {
#line 405 "dspevx.f"
	    dopgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &
		    work[indwrk], &iinfo, (ftnlen)1);
#line 407 "dspevx.f"
	    i__1 = *n - 1;
#line 407 "dspevx.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 408 "dspevx.f"
	    dsteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 410 "dspevx.f"
	    if (*info == 0) {
#line 411 "dspevx.f"
		i__1 = *n;
#line 411 "dspevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 412 "dspevx.f"
		    ifail[i__] = 0;
#line 413 "dspevx.f"
/* L10: */
#line 413 "dspevx.f"
		}
#line 414 "dspevx.f"
	    }
#line 415 "dspevx.f"
	}
#line 416 "dspevx.f"
	if (*info == 0) {
#line 417 "dspevx.f"
	    *m = *n;
#line 418 "dspevx.f"
	    goto L20;
#line 419 "dspevx.f"
	}
#line 420 "dspevx.f"
	*info = 0;
#line 421 "dspevx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 425 "dspevx.f"
    if (wantz) {
#line 426 "dspevx.f"
	*(unsigned char *)order = 'B';
#line 427 "dspevx.f"
    } else {
#line 428 "dspevx.f"
	*(unsigned char *)order = 'E';
#line 429 "dspevx.f"
    }
#line 430 "dspevx.f"
    indibl = 1;
#line 431 "dspevx.f"
    indisp = indibl + *n;
#line 432 "dspevx.f"
    indiwo = indisp + *n;
#line 433 "dspevx.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 438 "dspevx.f"
    if (wantz) {
#line 439 "dspevx.f"
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

#line 446 "dspevx.f"
	dopmtr_("L", uplo, "N", n, m, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 448 "dspevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 452 "dspevx.f"
L20:
#line 453 "dspevx.f"
    if (iscale == 1) {
#line 454 "dspevx.f"
	if (*info == 0) {
#line 455 "dspevx.f"
	    imax = *m;
#line 456 "dspevx.f"
	} else {
#line 457 "dspevx.f"
	    imax = *info - 1;
#line 458 "dspevx.f"
	}
#line 459 "dspevx.f"
	d__1 = 1. / sigma;
#line 459 "dspevx.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 460 "dspevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 465 "dspevx.f"
    if (wantz) {
#line 466 "dspevx.f"
	i__1 = *m - 1;
#line 466 "dspevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 467 "dspevx.f"
	    i__ = 0;
#line 468 "dspevx.f"
	    tmp1 = w[j];
#line 469 "dspevx.f"
	    i__2 = *m;
#line 469 "dspevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 470 "dspevx.f"
		if (w[jj] < tmp1) {
#line 471 "dspevx.f"
		    i__ = jj;
#line 472 "dspevx.f"
		    tmp1 = w[jj];
#line 473 "dspevx.f"
		}
#line 474 "dspevx.f"
/* L30: */
#line 474 "dspevx.f"
	    }

#line 476 "dspevx.f"
	    if (i__ != 0) {
#line 477 "dspevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 478 "dspevx.f"
		w[i__] = w[j];
#line 479 "dspevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 480 "dspevx.f"
		w[j] = tmp1;
#line 481 "dspevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 482 "dspevx.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 483 "dspevx.f"
		if (*info != 0) {
#line 484 "dspevx.f"
		    itmp1 = ifail[i__];
#line 485 "dspevx.f"
		    ifail[i__] = ifail[j];
#line 486 "dspevx.f"
		    ifail[j] = itmp1;
#line 487 "dspevx.f"
		}
#line 488 "dspevx.f"
	    }
#line 489 "dspevx.f"
/* L40: */
#line 489 "dspevx.f"
	}
#line 490 "dspevx.f"
    }

#line 492 "dspevx.f"
    return 0;

/*     End of DSPEVX */

} /* dspevx_ */

