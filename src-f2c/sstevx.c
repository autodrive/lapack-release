#line 1 "sstevx.f"
/* sstevx.f -- translated by f2c (version 20100827).
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

#line 1 "sstevx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SSTEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, */
/*                          M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE */
/*       INTEGER            IL, INFO, IU, LDZ, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric tridiagonal matrix A.  Eigenvalues and */
/* > eigenvectors can be selected by specifying either a range of values */
/* > or a range of indices for the desired eigenvalues. */
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
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A. */
/* >          On exit, D may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (max(1,N-1)) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A in elements 1 to N-1 of E. */
/* >          On exit, E may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
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
/* >          where EPS is the machine precision.  If ABSTOL is less */
/* >          than or equal to zero, then  EPS*|T|  will be used in */
/* >          its place, where |T| is the 1-norm of the tridiagonal */
/* >          matrix. */
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
/* >          Z is REAL array, dimension (LDZ, max(1,M) ) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          If an eigenvector fails to converge (INFO > 0), then that */
/* >          column of Z contains the latest approximation to the */
/* >          eigenvector, and the index of the eigenvector is returned */
/* >          in IFAIL.  If JOBZ = 'N', then Z is not referenced. */
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
/* >          WORK is REAL array, dimension (5*N) */
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

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int sstevx_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info, ftnlen jobz_len, ftnlen range_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, jj;
    static doublereal eps, vll, vuu, tmp1;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static doublereal tnrm;
    static integer itmp1;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
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
    static doublereal bignum;
    static integer indisp, indiwo, indwrk;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
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
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);


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

#line 268 "sstevx.f"
    /* Parameter adjustments */
#line 268 "sstevx.f"
    --d__;
#line 268 "sstevx.f"
    --e;
#line 268 "sstevx.f"
    --w;
#line 268 "sstevx.f"
    z_dim1 = *ldz;
#line 268 "sstevx.f"
    z_offset = 1 + z_dim1;
#line 268 "sstevx.f"
    z__ -= z_offset;
#line 268 "sstevx.f"
    --work;
#line 268 "sstevx.f"
    --iwork;
#line 268 "sstevx.f"
    --ifail;
#line 268 "sstevx.f"

#line 268 "sstevx.f"
    /* Function Body */
#line 268 "sstevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 269 "sstevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 270 "sstevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 271 "sstevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 273 "sstevx.f"
    *info = 0;
#line 274 "sstevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 275 "sstevx.f"
	*info = -1;
#line 276 "sstevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 277 "sstevx.f"
	*info = -2;
#line 278 "sstevx.f"
    } else if (*n < 0) {
#line 279 "sstevx.f"
	*info = -3;
#line 280 "sstevx.f"
    } else {
#line 281 "sstevx.f"
	if (valeig) {
#line 282 "sstevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 282 "sstevx.f"
		*info = -7;
#line 282 "sstevx.f"
	    }
#line 284 "sstevx.f"
	} else if (indeig) {
#line 285 "sstevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 286 "sstevx.f"
		*info = -8;
#line 287 "sstevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 288 "sstevx.f"
		*info = -9;
#line 289 "sstevx.f"
	    }
#line 290 "sstevx.f"
	}
#line 291 "sstevx.f"
    }
#line 292 "sstevx.f"
    if (*info == 0) {
#line 293 "sstevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 293 "sstevx.f"
	    *info = -14;
#line 293 "sstevx.f"
	}
#line 295 "sstevx.f"
    }

#line 297 "sstevx.f"
    if (*info != 0) {
#line 298 "sstevx.f"
	i__1 = -(*info);
#line 298 "sstevx.f"
	xerbla_("SSTEVX", &i__1, (ftnlen)6);
#line 299 "sstevx.f"
	return 0;
#line 300 "sstevx.f"
    }

/*     Quick return if possible */

#line 304 "sstevx.f"
    *m = 0;
#line 305 "sstevx.f"
    if (*n == 0) {
#line 305 "sstevx.f"
	return 0;
#line 305 "sstevx.f"
    }

#line 308 "sstevx.f"
    if (*n == 1) {
#line 309 "sstevx.f"
	if (alleig || indeig) {
#line 310 "sstevx.f"
	    *m = 1;
#line 311 "sstevx.f"
	    w[1] = d__[1];
#line 312 "sstevx.f"
	} else {
#line 313 "sstevx.f"
	    if (*vl < d__[1] && *vu >= d__[1]) {
#line 314 "sstevx.f"
		*m = 1;
#line 315 "sstevx.f"
		w[1] = d__[1];
#line 316 "sstevx.f"
	    }
#line 317 "sstevx.f"
	}
#line 318 "sstevx.f"
	if (wantz) {
#line 318 "sstevx.f"
	    z__[z_dim1 + 1] = 1.;
#line 318 "sstevx.f"
	}
#line 320 "sstevx.f"
	return 0;
#line 321 "sstevx.f"
    }

/*     Get machine constants. */

#line 325 "sstevx.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 326 "sstevx.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 327 "sstevx.f"
    smlnum = safmin / eps;
#line 328 "sstevx.f"
    bignum = 1. / smlnum;
#line 329 "sstevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 330 "sstevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 330 "sstevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 334 "sstevx.f"
    iscale = 0;
#line 335 "sstevx.f"
    if (valeig) {
#line 336 "sstevx.f"
	vll = *vl;
#line 337 "sstevx.f"
	vuu = *vu;
#line 338 "sstevx.f"
    } else {
#line 339 "sstevx.f"
	vll = 0.;
#line 340 "sstevx.f"
	vuu = 0.;
#line 341 "sstevx.f"
    }
#line 342 "sstevx.f"
    tnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 343 "sstevx.f"
    if (tnrm > 0. && tnrm < rmin) {
#line 344 "sstevx.f"
	iscale = 1;
#line 345 "sstevx.f"
	sigma = rmin / tnrm;
#line 346 "sstevx.f"
    } else if (tnrm > rmax) {
#line 347 "sstevx.f"
	iscale = 1;
#line 348 "sstevx.f"
	sigma = rmax / tnrm;
#line 349 "sstevx.f"
    }
#line 350 "sstevx.f"
    if (iscale == 1) {
#line 351 "sstevx.f"
	sscal_(n, &sigma, &d__[1], &c__1);
#line 352 "sstevx.f"
	i__1 = *n - 1;
#line 352 "sstevx.f"
	sscal_(&i__1, &sigma, &e[1], &c__1);
#line 353 "sstevx.f"
	if (valeig) {
#line 354 "sstevx.f"
	    vll = *vl * sigma;
#line 355 "sstevx.f"
	    vuu = *vu * sigma;
#line 356 "sstevx.f"
	}
#line 357 "sstevx.f"
    }

/*     If all eigenvalues are desired and ABSTOL is less than zero, then */
/*     call SSTERF or SSTEQR.  If this fails for some eigenvalue, then */
/*     try SSTEBZ. */

#line 363 "sstevx.f"
    test = FALSE_;
#line 364 "sstevx.f"
    if (indeig) {
#line 365 "sstevx.f"
	if (*il == 1 && *iu == *n) {
#line 366 "sstevx.f"
	    test = TRUE_;
#line 367 "sstevx.f"
	}
#line 368 "sstevx.f"
    }
#line 369 "sstevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 370 "sstevx.f"
	scopy_(n, &d__[1], &c__1, &w[1], &c__1);
#line 371 "sstevx.f"
	i__1 = *n - 1;
#line 371 "sstevx.f"
	scopy_(&i__1, &e[1], &c__1, &work[1], &c__1);
#line 372 "sstevx.f"
	indwrk = *n + 1;
#line 373 "sstevx.f"
	if (! wantz) {
#line 374 "sstevx.f"
	    ssterf_(n, &w[1], &work[1], info);
#line 375 "sstevx.f"
	} else {
#line 376 "sstevx.f"
	    ssteqr_("I", n, &w[1], &work[1], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 377 "sstevx.f"
	    if (*info == 0) {
#line 378 "sstevx.f"
		i__1 = *n;
#line 378 "sstevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 379 "sstevx.f"
		    ifail[i__] = 0;
#line 380 "sstevx.f"
/* L10: */
#line 380 "sstevx.f"
		}
#line 381 "sstevx.f"
	    }
#line 382 "sstevx.f"
	}
#line 383 "sstevx.f"
	if (*info == 0) {
#line 384 "sstevx.f"
	    *m = *n;
#line 385 "sstevx.f"
	    goto L20;
#line 386 "sstevx.f"
	}
#line 387 "sstevx.f"
	*info = 0;
#line 388 "sstevx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 392 "sstevx.f"
    if (wantz) {
#line 393 "sstevx.f"
	*(unsigned char *)order = 'B';
#line 394 "sstevx.f"
    } else {
#line 395 "sstevx.f"
	*(unsigned char *)order = 'E';
#line 396 "sstevx.f"
    }
#line 397 "sstevx.f"
    indwrk = 1;
#line 398 "sstevx.f"
    indibl = 1;
#line 399 "sstevx.f"
    indisp = indibl + *n;
#line 400 "sstevx.f"
    indiwo = indisp + *n;
#line 401 "sstevx.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &d__[1], &e[1], m, &
	    nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[indwrk], &
	    iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 405 "sstevx.f"
    if (wantz) {
#line 406 "sstevx.f"
	sstein_(n, &d__[1], &e[1], m, &w[1], &iwork[indibl], &iwork[indisp], &
		z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &ifail[1], 
		info);
#line 409 "sstevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 413 "sstevx.f"
L20:
#line 414 "sstevx.f"
    if (iscale == 1) {
#line 415 "sstevx.f"
	if (*info == 0) {
#line 416 "sstevx.f"
	    imax = *m;
#line 417 "sstevx.f"
	} else {
#line 418 "sstevx.f"
	    imax = *info - 1;
#line 419 "sstevx.f"
	}
#line 420 "sstevx.f"
	d__1 = 1. / sigma;
#line 420 "sstevx.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 421 "sstevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 426 "sstevx.f"
    if (wantz) {
#line 427 "sstevx.f"
	i__1 = *m - 1;
#line 427 "sstevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 428 "sstevx.f"
	    i__ = 0;
#line 429 "sstevx.f"
	    tmp1 = w[j];
#line 430 "sstevx.f"
	    i__2 = *m;
#line 430 "sstevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 431 "sstevx.f"
		if (w[jj] < tmp1) {
#line 432 "sstevx.f"
		    i__ = jj;
#line 433 "sstevx.f"
		    tmp1 = w[jj];
#line 434 "sstevx.f"
		}
#line 435 "sstevx.f"
/* L30: */
#line 435 "sstevx.f"
	    }

#line 437 "sstevx.f"
	    if (i__ != 0) {
#line 438 "sstevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 439 "sstevx.f"
		w[i__] = w[j];
#line 440 "sstevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 441 "sstevx.f"
		w[j] = tmp1;
#line 442 "sstevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 443 "sstevx.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 444 "sstevx.f"
		if (*info != 0) {
#line 445 "sstevx.f"
		    itmp1 = ifail[i__];
#line 446 "sstevx.f"
		    ifail[i__] = ifail[j];
#line 447 "sstevx.f"
		    ifail[j] = itmp1;
#line 448 "sstevx.f"
		}
#line 449 "sstevx.f"
	    }
#line 450 "sstevx.f"
/* L40: */
#line 450 "sstevx.f"
	}
#line 451 "sstevx.f"
    }

#line 453 "sstevx.f"
    return 0;

/*     End of SSTEVX */

} /* sstevx_ */

