#line 1 "chpevx.f"
/* chpevx.f -- translated by f2c (version 20100827).
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

#line 1 "chpevx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CHPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, */
/*                          IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDZ, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPEVX computes selected eigenvalues and, optionally, eigenvectors */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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
/* >          Z is COMPLEX array, dimension (LDZ, max(1,M)) */
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

/* > \ingroup complexOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int chpevx_(char *jobz, char *range, char *uplo, integer *n, 
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
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), scopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical wantz, alleig, indeig;
    static integer iscale, indibl;
    extern doublereal clanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen);
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indiwk, indisp, indtau;
    extern /* Subroutine */ int chptrd_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, ftnlen), 
	    cstein_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, integer *, integer *);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), cupgtr_(char *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen), ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer nsplit;
    extern /* Subroutine */ int cupmtr_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 285 "chpevx.f"
    /* Parameter adjustments */
#line 285 "chpevx.f"
    --ap;
#line 285 "chpevx.f"
    --w;
#line 285 "chpevx.f"
    z_dim1 = *ldz;
#line 285 "chpevx.f"
    z_offset = 1 + z_dim1;
#line 285 "chpevx.f"
    z__ -= z_offset;
#line 285 "chpevx.f"
    --work;
#line 285 "chpevx.f"
    --rwork;
#line 285 "chpevx.f"
    --iwork;
#line 285 "chpevx.f"
    --ifail;
#line 285 "chpevx.f"

#line 285 "chpevx.f"
    /* Function Body */
#line 285 "chpevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 286 "chpevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 287 "chpevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 288 "chpevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 290 "chpevx.f"
    *info = 0;
#line 291 "chpevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 292 "chpevx.f"
	*info = -1;
#line 293 "chpevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 294 "chpevx.f"
	*info = -2;
#line 295 "chpevx.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 297 "chpevx.f"
	*info = -3;
#line 298 "chpevx.f"
    } else if (*n < 0) {
#line 299 "chpevx.f"
	*info = -4;
#line 300 "chpevx.f"
    } else {
#line 301 "chpevx.f"
	if (valeig) {
#line 302 "chpevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 302 "chpevx.f"
		*info = -7;
#line 302 "chpevx.f"
	    }
#line 304 "chpevx.f"
	} else if (indeig) {
#line 305 "chpevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 306 "chpevx.f"
		*info = -8;
#line 307 "chpevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 308 "chpevx.f"
		*info = -9;
#line 309 "chpevx.f"
	    }
#line 310 "chpevx.f"
	}
#line 311 "chpevx.f"
    }
#line 312 "chpevx.f"
    if (*info == 0) {
#line 313 "chpevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 313 "chpevx.f"
	    *info = -14;
#line 313 "chpevx.f"
	}
#line 315 "chpevx.f"
    }

#line 317 "chpevx.f"
    if (*info != 0) {
#line 318 "chpevx.f"
	i__1 = -(*info);
#line 318 "chpevx.f"
	xerbla_("CHPEVX", &i__1, (ftnlen)6);
#line 319 "chpevx.f"
	return 0;
#line 320 "chpevx.f"
    }

/*     Quick return if possible */

#line 324 "chpevx.f"
    *m = 0;
#line 325 "chpevx.f"
    if (*n == 0) {
#line 325 "chpevx.f"
	return 0;
#line 325 "chpevx.f"
    }

#line 328 "chpevx.f"
    if (*n == 1) {
#line 329 "chpevx.f"
	if (alleig || indeig) {
#line 330 "chpevx.f"
	    *m = 1;
#line 331 "chpevx.f"
	    w[1] = ap[1].r;
#line 332 "chpevx.f"
	} else {
#line 333 "chpevx.f"
	    if (*vl < ap[1].r && *vu >= ap[1].r) {
#line 334 "chpevx.f"
		*m = 1;
#line 335 "chpevx.f"
		w[1] = ap[1].r;
#line 336 "chpevx.f"
	    }
#line 337 "chpevx.f"
	}
#line 338 "chpevx.f"
	if (wantz) {
#line 338 "chpevx.f"
	    i__1 = z_dim1 + 1;
#line 338 "chpevx.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 338 "chpevx.f"
	}
#line 340 "chpevx.f"
	return 0;
#line 341 "chpevx.f"
    }

/*     Get machine constants. */

#line 345 "chpevx.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 346 "chpevx.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 347 "chpevx.f"
    smlnum = safmin / eps;
#line 348 "chpevx.f"
    bignum = 1. / smlnum;
#line 349 "chpevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 350 "chpevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 350 "chpevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 354 "chpevx.f"
    iscale = 0;
#line 355 "chpevx.f"
    abstll = *abstol;
#line 356 "chpevx.f"
    if (valeig) {
#line 357 "chpevx.f"
	vll = *vl;
#line 358 "chpevx.f"
	vuu = *vu;
#line 359 "chpevx.f"
    } else {
#line 360 "chpevx.f"
	vll = 0.;
#line 361 "chpevx.f"
	vuu = 0.;
#line 362 "chpevx.f"
    }
#line 363 "chpevx.f"
    anrm = clanhp_("M", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);
#line 364 "chpevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 365 "chpevx.f"
	iscale = 1;
#line 366 "chpevx.f"
	sigma = rmin / anrm;
#line 367 "chpevx.f"
    } else if (anrm > rmax) {
#line 368 "chpevx.f"
	iscale = 1;
#line 369 "chpevx.f"
	sigma = rmax / anrm;
#line 370 "chpevx.f"
    }
#line 371 "chpevx.f"
    if (iscale == 1) {
#line 372 "chpevx.f"
	i__1 = *n * (*n + 1) / 2;
#line 372 "chpevx.f"
	csscal_(&i__1, &sigma, &ap[1], &c__1);
#line 373 "chpevx.f"
	if (*abstol > 0.) {
#line 373 "chpevx.f"
	    abstll = *abstol * sigma;
#line 373 "chpevx.f"
	}
#line 375 "chpevx.f"
	if (valeig) {
#line 376 "chpevx.f"
	    vll = *vl * sigma;
#line 377 "chpevx.f"
	    vuu = *vu * sigma;
#line 378 "chpevx.f"
	}
#line 379 "chpevx.f"
    }

/*     Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form. */

#line 383 "chpevx.f"
    indd = 1;
#line 384 "chpevx.f"
    inde = indd + *n;
#line 385 "chpevx.f"
    indrwk = inde + *n;
#line 386 "chpevx.f"
    indtau = 1;
#line 387 "chpevx.f"
    indwrk = indtau + *n;
#line 388 "chpevx.f"
    chptrd_(uplo, n, &ap[1], &rwork[indd], &rwork[inde], &work[indtau], &
	    iinfo, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call SSTERF or CUPGTR and CSTEQR.  If this fails */
/*     for some eigenvalue, then try SSTEBZ. */

#line 395 "chpevx.f"
    test = FALSE_;
#line 396 "chpevx.f"
    if (indeig) {
#line 397 "chpevx.f"
	if (*il == 1 && *iu == *n) {
#line 398 "chpevx.f"
	    test = TRUE_;
#line 399 "chpevx.f"
	}
#line 400 "chpevx.f"
    }
#line 401 "chpevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 402 "chpevx.f"
	scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 403 "chpevx.f"
	indee = indrwk + (*n << 1);
#line 404 "chpevx.f"
	if (! wantz) {
#line 405 "chpevx.f"
	    i__1 = *n - 1;
#line 405 "chpevx.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 406 "chpevx.f"
	    ssterf_(n, &w[1], &rwork[indee], info);
#line 407 "chpevx.f"
	} else {
#line 408 "chpevx.f"
	    cupgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &
		    work[indwrk], &iinfo, (ftnlen)1);
#line 410 "chpevx.f"
	    i__1 = *n - 1;
#line 410 "chpevx.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 411 "chpevx.f"
	    csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 413 "chpevx.f"
	    if (*info == 0) {
#line 414 "chpevx.f"
		i__1 = *n;
#line 414 "chpevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 415 "chpevx.f"
		    ifail[i__] = 0;
#line 416 "chpevx.f"
/* L10: */
#line 416 "chpevx.f"
		}
#line 417 "chpevx.f"
	    }
#line 418 "chpevx.f"
	}
#line 419 "chpevx.f"
	if (*info == 0) {
#line 420 "chpevx.f"
	    *m = *n;
#line 421 "chpevx.f"
	    goto L20;
#line 422 "chpevx.f"
	}
#line 423 "chpevx.f"
	*info = 0;
#line 424 "chpevx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */

#line 428 "chpevx.f"
    if (wantz) {
#line 429 "chpevx.f"
	*(unsigned char *)order = 'B';
#line 430 "chpevx.f"
    } else {
#line 431 "chpevx.f"
	*(unsigned char *)order = 'E';
#line 432 "chpevx.f"
    }
#line 433 "chpevx.f"
    indibl = 1;
#line 434 "chpevx.f"
    indisp = indibl + *n;
#line 435 "chpevx.f"
    indiwk = indisp + *n;
#line 436 "chpevx.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 441 "chpevx.f"
    if (wantz) {
#line 442 "chpevx.f"
	cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by CSTEIN. */

#line 449 "chpevx.f"
	indwrk = indtau + *n;
#line 450 "chpevx.f"
	cupmtr_("L", uplo, "N", n, m, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 452 "chpevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 456 "chpevx.f"
L20:
#line 457 "chpevx.f"
    if (iscale == 1) {
#line 458 "chpevx.f"
	if (*info == 0) {
#line 459 "chpevx.f"
	    imax = *m;
#line 460 "chpevx.f"
	} else {
#line 461 "chpevx.f"
	    imax = *info - 1;
#line 462 "chpevx.f"
	}
#line 463 "chpevx.f"
	d__1 = 1. / sigma;
#line 463 "chpevx.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 464 "chpevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 469 "chpevx.f"
    if (wantz) {
#line 470 "chpevx.f"
	i__1 = *m - 1;
#line 470 "chpevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 471 "chpevx.f"
	    i__ = 0;
#line 472 "chpevx.f"
	    tmp1 = w[j];
#line 473 "chpevx.f"
	    i__2 = *m;
#line 473 "chpevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 474 "chpevx.f"
		if (w[jj] < tmp1) {
#line 475 "chpevx.f"
		    i__ = jj;
#line 476 "chpevx.f"
		    tmp1 = w[jj];
#line 477 "chpevx.f"
		}
#line 478 "chpevx.f"
/* L30: */
#line 478 "chpevx.f"
	    }

#line 480 "chpevx.f"
	    if (i__ != 0) {
#line 481 "chpevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 482 "chpevx.f"
		w[i__] = w[j];
#line 483 "chpevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 484 "chpevx.f"
		w[j] = tmp1;
#line 485 "chpevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 486 "chpevx.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 487 "chpevx.f"
		if (*info != 0) {
#line 488 "chpevx.f"
		    itmp1 = ifail[i__];
#line 489 "chpevx.f"
		    ifail[i__] = ifail[j];
#line 490 "chpevx.f"
		    ifail[j] = itmp1;
#line 491 "chpevx.f"
		}
#line 492 "chpevx.f"
	    }
#line 493 "chpevx.f"
/* L40: */
#line 493 "chpevx.f"
	}
#line 494 "chpevx.f"
    }

#line 496 "chpevx.f"
    return 0;

/*     End of CHPEVX */

} /* chpevx_ */

