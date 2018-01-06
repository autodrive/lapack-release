#line 1 "zheevx.f"
/* zheevx.f -- translated by f2c (version 20100827).
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

#line 1 "zheevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> ZHEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can */
/* > be selected by specifying either a range of values or a range of */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= 1, when N <= 1; */
/* >          otherwise 2*N. */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the max of the blocksize for ZHETRD and for */
/* >          ZUNMTR as returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \date November 2011 */

/* > \ingroup complex16HEeigen */

/*  ===================================================================== */
/* Subroutine */ int zheevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *iwork, integer *ifail, integer *
	info, ftnlen jobz_len, ftnlen range_len, ftnlen uplo_len)
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
	    doublereal *, integer *);
    static logical lower, wantz;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal abstll, bignum;
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer indiwk, indisp, indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dstebz_(char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int zhetrd_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer lwkmin;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer llwork, nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int zstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), zungtr_(char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen), 
	    zunmtr_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, ftnlen, ftnlen, 
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

#line 308 "zheevx.f"
    /* Parameter adjustments */
#line 308 "zheevx.f"
    a_dim1 = *lda;
#line 308 "zheevx.f"
    a_offset = 1 + a_dim1;
#line 308 "zheevx.f"
    a -= a_offset;
#line 308 "zheevx.f"
    --w;
#line 308 "zheevx.f"
    z_dim1 = *ldz;
#line 308 "zheevx.f"
    z_offset = 1 + z_dim1;
#line 308 "zheevx.f"
    z__ -= z_offset;
#line 308 "zheevx.f"
    --work;
#line 308 "zheevx.f"
    --rwork;
#line 308 "zheevx.f"
    --iwork;
#line 308 "zheevx.f"
    --ifail;
#line 308 "zheevx.f"

#line 308 "zheevx.f"
    /* Function Body */
#line 308 "zheevx.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 309 "zheevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 310 "zheevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 311 "zheevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 312 "zheevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 313 "zheevx.f"
    lquery = *lwork == -1;

#line 315 "zheevx.f"
    *info = 0;
#line 316 "zheevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 317 "zheevx.f"
	*info = -1;
#line 318 "zheevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 319 "zheevx.f"
	*info = -2;
#line 320 "zheevx.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 321 "zheevx.f"
	*info = -3;
#line 322 "zheevx.f"
    } else if (*n < 0) {
#line 323 "zheevx.f"
	*info = -4;
#line 324 "zheevx.f"
    } else if (*lda < max(1,*n)) {
#line 325 "zheevx.f"
	*info = -6;
#line 326 "zheevx.f"
    } else {
#line 327 "zheevx.f"
	if (valeig) {
#line 328 "zheevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 328 "zheevx.f"
		*info = -8;
#line 328 "zheevx.f"
	    }
#line 330 "zheevx.f"
	} else if (indeig) {
#line 331 "zheevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 332 "zheevx.f"
		*info = -9;
#line 333 "zheevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 334 "zheevx.f"
		*info = -10;
#line 335 "zheevx.f"
	    }
#line 336 "zheevx.f"
	}
#line 337 "zheevx.f"
    }
#line 338 "zheevx.f"
    if (*info == 0) {
#line 339 "zheevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 340 "zheevx.f"
	    *info = -15;
#line 341 "zheevx.f"
	}
#line 342 "zheevx.f"
    }

#line 344 "zheevx.f"
    if (*info == 0) {
#line 345 "zheevx.f"
	if (*n <= 1) {
#line 346 "zheevx.f"
	    lwkmin = 1;
#line 347 "zheevx.f"
	    work[1].r = (doublereal) lwkmin, work[1].i = 0.;
#line 348 "zheevx.f"
	} else {
#line 349 "zheevx.f"
	    lwkmin = *n << 1;
#line 350 "zheevx.f"
	    nb = ilaenv_(&c__1, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 351 "zheevx.f"
	    i__1 = nb, i__2 = ilaenv_(&c__1, "ZUNMTR", uplo, n, &c_n1, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1);
#line 351 "zheevx.f"
	    nb = max(i__1,i__2);
/* Computing MAX */
#line 352 "zheevx.f"
	    i__1 = 1, i__2 = (nb + 1) * *n;
#line 352 "zheevx.f"
	    lwkopt = max(i__1,i__2);
#line 353 "zheevx.f"
	    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 354 "zheevx.f"
	}

#line 356 "zheevx.f"
	if (*lwork < lwkmin && ! lquery) {
#line 356 "zheevx.f"
	    *info = -17;
#line 356 "zheevx.f"
	}
#line 358 "zheevx.f"
    }

#line 360 "zheevx.f"
    if (*info != 0) {
#line 361 "zheevx.f"
	i__1 = -(*info);
#line 361 "zheevx.f"
	xerbla_("ZHEEVX", &i__1, (ftnlen)6);
#line 362 "zheevx.f"
	return 0;
#line 363 "zheevx.f"
    } else if (lquery) {
#line 364 "zheevx.f"
	return 0;
#line 365 "zheevx.f"
    }

/*     Quick return if possible */

#line 369 "zheevx.f"
    *m = 0;
#line 370 "zheevx.f"
    if (*n == 0) {
#line 371 "zheevx.f"
	return 0;
#line 372 "zheevx.f"
    }

#line 374 "zheevx.f"
    if (*n == 1) {
#line 375 "zheevx.f"
	if (alleig || indeig) {
#line 376 "zheevx.f"
	    *m = 1;
#line 377 "zheevx.f"
	    i__1 = a_dim1 + 1;
#line 377 "zheevx.f"
	    w[1] = a[i__1].r;
#line 378 "zheevx.f"
	} else if (valeig) {
#line 379 "zheevx.f"
	    i__1 = a_dim1 + 1;
#line 379 "zheevx.f"
	    i__2 = a_dim1 + 1;
#line 379 "zheevx.f"
	    if (*vl < a[i__1].r && *vu >= a[i__2].r) {
#line 381 "zheevx.f"
		*m = 1;
#line 382 "zheevx.f"
		i__1 = a_dim1 + 1;
#line 382 "zheevx.f"
		w[1] = a[i__1].r;
#line 383 "zheevx.f"
	    }
#line 384 "zheevx.f"
	}
#line 385 "zheevx.f"
	if (wantz) {
#line 385 "zheevx.f"
	    i__1 = z_dim1 + 1;
#line 385 "zheevx.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 385 "zheevx.f"
	}
#line 387 "zheevx.f"
	return 0;
#line 388 "zheevx.f"
    }

/*     Get machine constants. */

#line 392 "zheevx.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 393 "zheevx.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 394 "zheevx.f"
    smlnum = safmin / eps;
#line 395 "zheevx.f"
    bignum = 1. / smlnum;
#line 396 "zheevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 397 "zheevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 397 "zheevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 401 "zheevx.f"
    iscale = 0;
#line 402 "zheevx.f"
    abstll = *abstol;
#line 403 "zheevx.f"
    if (valeig) {
#line 404 "zheevx.f"
	vll = *vl;
#line 405 "zheevx.f"
	vuu = *vu;
#line 406 "zheevx.f"
    }
#line 407 "zheevx.f"
    anrm = zlanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 408 "zheevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 409 "zheevx.f"
	iscale = 1;
#line 410 "zheevx.f"
	sigma = rmin / anrm;
#line 411 "zheevx.f"
    } else if (anrm > rmax) {
#line 412 "zheevx.f"
	iscale = 1;
#line 413 "zheevx.f"
	sigma = rmax / anrm;
#line 414 "zheevx.f"
    }
#line 415 "zheevx.f"
    if (iscale == 1) {
#line 416 "zheevx.f"
	if (lower) {
#line 417 "zheevx.f"
	    i__1 = *n;
#line 417 "zheevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 418 "zheevx.f"
		i__2 = *n - j + 1;
#line 418 "zheevx.f"
		zdscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 419 "zheevx.f"
/* L10: */
#line 419 "zheevx.f"
	    }
#line 420 "zheevx.f"
	} else {
#line 421 "zheevx.f"
	    i__1 = *n;
#line 421 "zheevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 422 "zheevx.f"
		zdscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 423 "zheevx.f"
/* L20: */
#line 423 "zheevx.f"
	    }
#line 424 "zheevx.f"
	}
#line 425 "zheevx.f"
	if (*abstol > 0.) {
#line 425 "zheevx.f"
	    abstll = *abstol * sigma;
#line 425 "zheevx.f"
	}
#line 427 "zheevx.f"
	if (valeig) {
#line 428 "zheevx.f"
	    vll = *vl * sigma;
#line 429 "zheevx.f"
	    vuu = *vu * sigma;
#line 430 "zheevx.f"
	}
#line 431 "zheevx.f"
    }

/*     Call ZHETRD to reduce Hermitian matrix to tridiagonal form. */

#line 435 "zheevx.f"
    indd = 1;
#line 436 "zheevx.f"
    inde = indd + *n;
#line 437 "zheevx.f"
    indrwk = inde + *n;
#line 438 "zheevx.f"
    indtau = 1;
#line 439 "zheevx.f"
    indwrk = indtau + *n;
#line 440 "zheevx.f"
    llwork = *lwork - indwrk + 1;
#line 441 "zheevx.f"
    zhetrd_(uplo, n, &a[a_offset], lda, &rwork[indd], &rwork[inde], &work[
	    indtau], &work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call DSTERF or ZUNGTR and ZSTEQR.  If this fails for */
/*     some eigenvalue, then try DSTEBZ. */

#line 448 "zheevx.f"
    test = FALSE_;
#line 449 "zheevx.f"
    if (indeig) {
#line 450 "zheevx.f"
	if (*il == 1 && *iu == *n) {
#line 451 "zheevx.f"
	    test = TRUE_;
#line 452 "zheevx.f"
	}
#line 453 "zheevx.f"
    }
#line 454 "zheevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 455 "zheevx.f"
	dcopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 456 "zheevx.f"
	indee = indrwk + (*n << 1);
#line 457 "zheevx.f"
	if (! wantz) {
#line 458 "zheevx.f"
	    i__1 = *n - 1;
#line 458 "zheevx.f"
	    dcopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 459 "zheevx.f"
	    dsterf_(n, &w[1], &rwork[indee], info);
#line 460 "zheevx.f"
	} else {
#line 461 "zheevx.f"
	    zlacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 462 "zheevx.f"
	    zungtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
#line 464 "zheevx.f"
	    i__1 = *n - 1;
#line 464 "zheevx.f"
	    dcopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 465 "zheevx.f"
	    zsteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 467 "zheevx.f"
	    if (*info == 0) {
#line 468 "zheevx.f"
		i__1 = *n;
#line 468 "zheevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 469 "zheevx.f"
		    ifail[i__] = 0;
#line 470 "zheevx.f"
/* L30: */
#line 470 "zheevx.f"
		}
#line 471 "zheevx.f"
	    }
#line 472 "zheevx.f"
	}
#line 473 "zheevx.f"
	if (*info == 0) {
#line 474 "zheevx.f"
	    *m = *n;
#line 475 "zheevx.f"
	    goto L40;
#line 476 "zheevx.f"
	}
#line 477 "zheevx.f"
	*info = 0;
#line 478 "zheevx.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN. */

#line 482 "zheevx.f"
    if (wantz) {
#line 483 "zheevx.f"
	*(unsigned char *)order = 'B';
#line 484 "zheevx.f"
    } else {
#line 485 "zheevx.f"
	*(unsigned char *)order = 'E';
#line 486 "zheevx.f"
    }
#line 487 "zheevx.f"
    indibl = 1;
#line 488 "zheevx.f"
    indisp = indibl + *n;
#line 489 "zheevx.f"
    indiwk = indisp + *n;
#line 490 "zheevx.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 495 "zheevx.f"
    if (wantz) {
#line 496 "zheevx.f"
	zstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by ZSTEIN. */

#line 503 "zheevx.f"
	zunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwrk], &llwork, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 505 "zheevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 509 "zheevx.f"
L40:
#line 510 "zheevx.f"
    if (iscale == 1) {
#line 511 "zheevx.f"
	if (*info == 0) {
#line 512 "zheevx.f"
	    imax = *m;
#line 513 "zheevx.f"
	} else {
#line 514 "zheevx.f"
	    imax = *info - 1;
#line 515 "zheevx.f"
	}
#line 516 "zheevx.f"
	d__1 = 1. / sigma;
#line 516 "zheevx.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 517 "zheevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 522 "zheevx.f"
    if (wantz) {
#line 523 "zheevx.f"
	i__1 = *m - 1;
#line 523 "zheevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 524 "zheevx.f"
	    i__ = 0;
#line 525 "zheevx.f"
	    tmp1 = w[j];
#line 526 "zheevx.f"
	    i__2 = *m;
#line 526 "zheevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 527 "zheevx.f"
		if (w[jj] < tmp1) {
#line 528 "zheevx.f"
		    i__ = jj;
#line 529 "zheevx.f"
		    tmp1 = w[jj];
#line 530 "zheevx.f"
		}
#line 531 "zheevx.f"
/* L50: */
#line 531 "zheevx.f"
	    }

#line 533 "zheevx.f"
	    if (i__ != 0) {
#line 534 "zheevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 535 "zheevx.f"
		w[i__] = w[j];
#line 536 "zheevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 537 "zheevx.f"
		w[j] = tmp1;
#line 538 "zheevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 539 "zheevx.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 540 "zheevx.f"
		if (*info != 0) {
#line 541 "zheevx.f"
		    itmp1 = ifail[i__];
#line 542 "zheevx.f"
		    ifail[i__] = ifail[j];
#line 543 "zheevx.f"
		    ifail[j] = itmp1;
#line 544 "zheevx.f"
		}
#line 545 "zheevx.f"
	    }
#line 546 "zheevx.f"
/* L60: */
#line 546 "zheevx.f"
	}
#line 547 "zheevx.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 551 "zheevx.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 553 "zheevx.f"
    return 0;

/*     End of ZHEEVX */

} /* zheevx_ */

