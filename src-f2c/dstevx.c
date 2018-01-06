#line 1 "dstevx.f"
/* dstevx.f -- translated by f2c (version 20100827).
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

#line 1 "dstevx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DSTEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, */
/*                          M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE */
/*       INTEGER            IL, INFO, IU, LDZ, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEVX computes selected eigenvalues and, optionally, eigenvectors */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A. */
/* >          On exit, D may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (max(1,N-1)) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A in elements 1 to N-1 of E. */
/* >          On exit, E may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
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
/* >          where EPS is the machine precision.  If ABSTOL is less */
/* >          than or equal to zero, then  EPS*|T|  will be used in */
/* >          its place, where |T| is the 1-norm of the tridiagonal */
/* >          matrix. */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) ) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (5*N) */
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
/* Subroutine */ int dstevx_(char *jobz, char *range, integer *n, doublereal *
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
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
    static doublereal bignum;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
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

#line 275 "dstevx.f"
    /* Parameter adjustments */
#line 275 "dstevx.f"
    --d__;
#line 275 "dstevx.f"
    --e;
#line 275 "dstevx.f"
    --w;
#line 275 "dstevx.f"
    z_dim1 = *ldz;
#line 275 "dstevx.f"
    z_offset = 1 + z_dim1;
#line 275 "dstevx.f"
    z__ -= z_offset;
#line 275 "dstevx.f"
    --work;
#line 275 "dstevx.f"
    --iwork;
#line 275 "dstevx.f"
    --ifail;
#line 275 "dstevx.f"

#line 275 "dstevx.f"
    /* Function Body */
#line 275 "dstevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 276 "dstevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 277 "dstevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 278 "dstevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 280 "dstevx.f"
    *info = 0;
#line 281 "dstevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 282 "dstevx.f"
	*info = -1;
#line 283 "dstevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 284 "dstevx.f"
	*info = -2;
#line 285 "dstevx.f"
    } else if (*n < 0) {
#line 286 "dstevx.f"
	*info = -3;
#line 287 "dstevx.f"
    } else {
#line 288 "dstevx.f"
	if (valeig) {
#line 289 "dstevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 289 "dstevx.f"
		*info = -7;
#line 289 "dstevx.f"
	    }
#line 291 "dstevx.f"
	} else if (indeig) {
#line 292 "dstevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 293 "dstevx.f"
		*info = -8;
#line 294 "dstevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 295 "dstevx.f"
		*info = -9;
#line 296 "dstevx.f"
	    }
#line 297 "dstevx.f"
	}
#line 298 "dstevx.f"
    }
#line 299 "dstevx.f"
    if (*info == 0) {
#line 300 "dstevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 300 "dstevx.f"
	    *info = -14;
#line 300 "dstevx.f"
	}
#line 302 "dstevx.f"
    }

#line 304 "dstevx.f"
    if (*info != 0) {
#line 305 "dstevx.f"
	i__1 = -(*info);
#line 305 "dstevx.f"
	xerbla_("DSTEVX", &i__1, (ftnlen)6);
#line 306 "dstevx.f"
	return 0;
#line 307 "dstevx.f"
    }

/*     Quick return if possible */

#line 311 "dstevx.f"
    *m = 0;
#line 312 "dstevx.f"
    if (*n == 0) {
#line 312 "dstevx.f"
	return 0;
#line 312 "dstevx.f"
    }

#line 315 "dstevx.f"
    if (*n == 1) {
#line 316 "dstevx.f"
	if (alleig || indeig) {
#line 317 "dstevx.f"
	    *m = 1;
#line 318 "dstevx.f"
	    w[1] = d__[1];
#line 319 "dstevx.f"
	} else {
#line 320 "dstevx.f"
	    if (*vl < d__[1] && *vu >= d__[1]) {
#line 321 "dstevx.f"
		*m = 1;
#line 322 "dstevx.f"
		w[1] = d__[1];
#line 323 "dstevx.f"
	    }
#line 324 "dstevx.f"
	}
#line 325 "dstevx.f"
	if (wantz) {
#line 325 "dstevx.f"
	    z__[z_dim1 + 1] = 1.;
#line 325 "dstevx.f"
	}
#line 327 "dstevx.f"
	return 0;
#line 328 "dstevx.f"
    }

/*     Get machine constants. */

#line 332 "dstevx.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 333 "dstevx.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 334 "dstevx.f"
    smlnum = safmin / eps;
#line 335 "dstevx.f"
    bignum = 1. / smlnum;
#line 336 "dstevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 337 "dstevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 337 "dstevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 341 "dstevx.f"
    iscale = 0;
#line 342 "dstevx.f"
    if (valeig) {
#line 343 "dstevx.f"
	vll = *vl;
#line 344 "dstevx.f"
	vuu = *vu;
#line 345 "dstevx.f"
    } else {
#line 346 "dstevx.f"
	vll = 0.;
#line 347 "dstevx.f"
	vuu = 0.;
#line 348 "dstevx.f"
    }
#line 349 "dstevx.f"
    tnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 350 "dstevx.f"
    if (tnrm > 0. && tnrm < rmin) {
#line 351 "dstevx.f"
	iscale = 1;
#line 352 "dstevx.f"
	sigma = rmin / tnrm;
#line 353 "dstevx.f"
    } else if (tnrm > rmax) {
#line 354 "dstevx.f"
	iscale = 1;
#line 355 "dstevx.f"
	sigma = rmax / tnrm;
#line 356 "dstevx.f"
    }
#line 357 "dstevx.f"
    if (iscale == 1) {
#line 358 "dstevx.f"
	dscal_(n, &sigma, &d__[1], &c__1);
#line 359 "dstevx.f"
	i__1 = *n - 1;
#line 359 "dstevx.f"
	dscal_(&i__1, &sigma, &e[1], &c__1);
#line 360 "dstevx.f"
	if (valeig) {
#line 361 "dstevx.f"
	    vll = *vl * sigma;
#line 362 "dstevx.f"
	    vuu = *vu * sigma;
#line 363 "dstevx.f"
	}
#line 364 "dstevx.f"
    }

/*     If all eigenvalues are desired and ABSTOL is less than zero, then */
/*     call DSTERF or SSTEQR.  If this fails for some eigenvalue, then */
/*     try DSTEBZ. */

#line 370 "dstevx.f"
    test = FALSE_;
#line 371 "dstevx.f"
    if (indeig) {
#line 372 "dstevx.f"
	if (*il == 1 && *iu == *n) {
#line 373 "dstevx.f"
	    test = TRUE_;
#line 374 "dstevx.f"
	}
#line 375 "dstevx.f"
    }
#line 376 "dstevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 377 "dstevx.f"
	dcopy_(n, &d__[1], &c__1, &w[1], &c__1);
#line 378 "dstevx.f"
	i__1 = *n - 1;
#line 378 "dstevx.f"
	dcopy_(&i__1, &e[1], &c__1, &work[1], &c__1);
#line 379 "dstevx.f"
	indwrk = *n + 1;
#line 380 "dstevx.f"
	if (! wantz) {
#line 381 "dstevx.f"
	    dsterf_(n, &w[1], &work[1], info);
#line 382 "dstevx.f"
	} else {
#line 383 "dstevx.f"
	    dsteqr_("I", n, &w[1], &work[1], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 384 "dstevx.f"
	    if (*info == 0) {
#line 385 "dstevx.f"
		i__1 = *n;
#line 385 "dstevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "dstevx.f"
		    ifail[i__] = 0;
#line 387 "dstevx.f"
/* L10: */
#line 387 "dstevx.f"
		}
#line 388 "dstevx.f"
	    }
#line 389 "dstevx.f"
	}
#line 390 "dstevx.f"
	if (*info == 0) {
#line 391 "dstevx.f"
	    *m = *n;
#line 392 "dstevx.f"
	    goto L20;
#line 393 "dstevx.f"
	}
#line 394 "dstevx.f"
	*info = 0;
#line 395 "dstevx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 399 "dstevx.f"
    if (wantz) {
#line 400 "dstevx.f"
	*(unsigned char *)order = 'B';
#line 401 "dstevx.f"
    } else {
#line 402 "dstevx.f"
	*(unsigned char *)order = 'E';
#line 403 "dstevx.f"
    }
#line 404 "dstevx.f"
    indwrk = 1;
#line 405 "dstevx.f"
    indibl = 1;
#line 406 "dstevx.f"
    indisp = indibl + *n;
#line 407 "dstevx.f"
    indiwo = indisp + *n;
#line 408 "dstevx.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &d__[1], &e[1], m, &
	    nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[indwrk], &
	    iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 412 "dstevx.f"
    if (wantz) {
#line 413 "dstevx.f"
	dstein_(n, &d__[1], &e[1], m, &w[1], &iwork[indibl], &iwork[indisp], &
		z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &ifail[1], 
		info);
#line 416 "dstevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 420 "dstevx.f"
L20:
#line 421 "dstevx.f"
    if (iscale == 1) {
#line 422 "dstevx.f"
	if (*info == 0) {
#line 423 "dstevx.f"
	    imax = *m;
#line 424 "dstevx.f"
	} else {
#line 425 "dstevx.f"
	    imax = *info - 1;
#line 426 "dstevx.f"
	}
#line 427 "dstevx.f"
	d__1 = 1. / sigma;
#line 427 "dstevx.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 428 "dstevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 433 "dstevx.f"
    if (wantz) {
#line 434 "dstevx.f"
	i__1 = *m - 1;
#line 434 "dstevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 435 "dstevx.f"
	    i__ = 0;
#line 436 "dstevx.f"
	    tmp1 = w[j];
#line 437 "dstevx.f"
	    i__2 = *m;
#line 437 "dstevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 438 "dstevx.f"
		if (w[jj] < tmp1) {
#line 439 "dstevx.f"
		    i__ = jj;
#line 440 "dstevx.f"
		    tmp1 = w[jj];
#line 441 "dstevx.f"
		}
#line 442 "dstevx.f"
/* L30: */
#line 442 "dstevx.f"
	    }

#line 444 "dstevx.f"
	    if (i__ != 0) {
#line 445 "dstevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 446 "dstevx.f"
		w[i__] = w[j];
#line 447 "dstevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 448 "dstevx.f"
		w[j] = tmp1;
#line 449 "dstevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 450 "dstevx.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 451 "dstevx.f"
		if (*info != 0) {
#line 452 "dstevx.f"
		    itmp1 = ifail[i__];
#line 453 "dstevx.f"
		    ifail[i__] = ifail[j];
#line 454 "dstevx.f"
		    ifail[j] = itmp1;
#line 455 "dstevx.f"
		}
#line 456 "dstevx.f"
	    }
#line 457 "dstevx.f"
/* L40: */
#line 457 "dstevx.f"
	}
#line 458 "dstevx.f"
    }

#line 460 "dstevx.f"
    return 0;

/*     End of DSTEVX */

} /* dstevx_ */

