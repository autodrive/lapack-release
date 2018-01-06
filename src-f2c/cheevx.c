#line 1 "cheevx.f"
/* cheevx.f -- translated by f2c (version 20100827).
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

#line 1 "cheevx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> CHEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheevx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheevx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheevx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEEVX computes selected eigenvalues and, optionally, eigenvectors */
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
/* >          A is COMPLEX array, dimension (LDA, N) */
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
/* >          Z is COMPLEX array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= 1, when N <= 1; */
/* >          otherwise 2*N. */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the max of the blocksize for CHETRD and for */
/* >          CUNMTR as returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \date June 2016 */

/* > \ingroup complexHEeigen */

/*  ===================================================================== */
/* Subroutine */ int cheevx_(char *jobz, char *range, char *uplo, integer *n, 
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
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical lower;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantz;
    extern doublereal clanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), chetrd_(char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen), 
	    clacpy_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indiwk, indisp, indtau;
    extern /* Subroutine */ int cstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static integer indrwk, indwrk, lwkmin;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), cungtr_(char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen), 
	    ssterf_(integer *, doublereal *, doublereal *, integer *), 
	    cunmtr_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer nsplit, llwork;
    static doublereal smlnum;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;


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

#line 315 "cheevx.f"
    /* Parameter adjustments */
#line 315 "cheevx.f"
    a_dim1 = *lda;
#line 315 "cheevx.f"
    a_offset = 1 + a_dim1;
#line 315 "cheevx.f"
    a -= a_offset;
#line 315 "cheevx.f"
    --w;
#line 315 "cheevx.f"
    z_dim1 = *ldz;
#line 315 "cheevx.f"
    z_offset = 1 + z_dim1;
#line 315 "cheevx.f"
    z__ -= z_offset;
#line 315 "cheevx.f"
    --work;
#line 315 "cheevx.f"
    --rwork;
#line 315 "cheevx.f"
    --iwork;
#line 315 "cheevx.f"
    --ifail;
#line 315 "cheevx.f"

#line 315 "cheevx.f"
    /* Function Body */
#line 315 "cheevx.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 316 "cheevx.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 317 "cheevx.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 318 "cheevx.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 319 "cheevx.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 320 "cheevx.f"
    lquery = *lwork == -1;

#line 322 "cheevx.f"
    *info = 0;
#line 323 "cheevx.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 324 "cheevx.f"
	*info = -1;
#line 325 "cheevx.f"
    } else if (! (alleig || valeig || indeig)) {
#line 326 "cheevx.f"
	*info = -2;
#line 327 "cheevx.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 328 "cheevx.f"
	*info = -3;
#line 329 "cheevx.f"
    } else if (*n < 0) {
#line 330 "cheevx.f"
	*info = -4;
#line 331 "cheevx.f"
    } else if (*lda < max(1,*n)) {
#line 332 "cheevx.f"
	*info = -6;
#line 333 "cheevx.f"
    } else {
#line 334 "cheevx.f"
	if (valeig) {
#line 335 "cheevx.f"
	    if (*n > 0 && *vu <= *vl) {
#line 335 "cheevx.f"
		*info = -8;
#line 335 "cheevx.f"
	    }
#line 337 "cheevx.f"
	} else if (indeig) {
#line 338 "cheevx.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 339 "cheevx.f"
		*info = -9;
#line 340 "cheevx.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 341 "cheevx.f"
		*info = -10;
#line 342 "cheevx.f"
	    }
#line 343 "cheevx.f"
	}
#line 344 "cheevx.f"
    }
#line 345 "cheevx.f"
    if (*info == 0) {
#line 346 "cheevx.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 347 "cheevx.f"
	    *info = -15;
#line 348 "cheevx.f"
	}
#line 349 "cheevx.f"
    }

#line 351 "cheevx.f"
    if (*info == 0) {
#line 352 "cheevx.f"
	if (*n <= 1) {
#line 353 "cheevx.f"
	    lwkmin = 1;
#line 354 "cheevx.f"
	    work[1].r = (doublereal) lwkmin, work[1].i = 0.;
#line 355 "cheevx.f"
	} else {
#line 356 "cheevx.f"
	    lwkmin = *n << 1;
#line 357 "cheevx.f"
	    nb = ilaenv_(&c__1, "CHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 358 "cheevx.f"
	    i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMTR", uplo, n, &c_n1, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1);
#line 358 "cheevx.f"
	    nb = max(i__1,i__2);
/* Computing MAX */
#line 359 "cheevx.f"
	    i__1 = 1, i__2 = (nb + 1) * *n;
#line 359 "cheevx.f"
	    lwkopt = max(i__1,i__2);
#line 360 "cheevx.f"
	    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 361 "cheevx.f"
	}

#line 363 "cheevx.f"
	if (*lwork < lwkmin && ! lquery) {
#line 363 "cheevx.f"
	    *info = -17;
#line 363 "cheevx.f"
	}
#line 365 "cheevx.f"
    }

#line 367 "cheevx.f"
    if (*info != 0) {
#line 368 "cheevx.f"
	i__1 = -(*info);
#line 368 "cheevx.f"
	xerbla_("CHEEVX", &i__1, (ftnlen)6);
#line 369 "cheevx.f"
	return 0;
#line 370 "cheevx.f"
    } else if (lquery) {
#line 371 "cheevx.f"
	return 0;
#line 372 "cheevx.f"
    }

/*     Quick return if possible */

#line 376 "cheevx.f"
    *m = 0;
#line 377 "cheevx.f"
    if (*n == 0) {
#line 378 "cheevx.f"
	return 0;
#line 379 "cheevx.f"
    }

#line 381 "cheevx.f"
    if (*n == 1) {
#line 382 "cheevx.f"
	if (alleig || indeig) {
#line 383 "cheevx.f"
	    *m = 1;
#line 384 "cheevx.f"
	    i__1 = a_dim1 + 1;
#line 384 "cheevx.f"
	    w[1] = a[i__1].r;
#line 385 "cheevx.f"
	} else if (valeig) {
#line 386 "cheevx.f"
	    i__1 = a_dim1 + 1;
#line 386 "cheevx.f"
	    i__2 = a_dim1 + 1;
#line 386 "cheevx.f"
	    if (*vl < a[i__1].r && *vu >= a[i__2].r) {
#line 388 "cheevx.f"
		*m = 1;
#line 389 "cheevx.f"
		i__1 = a_dim1 + 1;
#line 389 "cheevx.f"
		w[1] = a[i__1].r;
#line 390 "cheevx.f"
	    }
#line 391 "cheevx.f"
	}
#line 392 "cheevx.f"
	if (wantz) {
#line 392 "cheevx.f"
	    i__1 = z_dim1 + 1;
#line 392 "cheevx.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 392 "cheevx.f"
	}
#line 394 "cheevx.f"
	return 0;
#line 395 "cheevx.f"
    }

/*     Get machine constants. */

#line 399 "cheevx.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 400 "cheevx.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 401 "cheevx.f"
    smlnum = safmin / eps;
#line 402 "cheevx.f"
    bignum = 1. / smlnum;
#line 403 "cheevx.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 404 "cheevx.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 404 "cheevx.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 408 "cheevx.f"
    iscale = 0;
#line 409 "cheevx.f"
    abstll = *abstol;
#line 410 "cheevx.f"
    if (valeig) {
#line 411 "cheevx.f"
	vll = *vl;
#line 412 "cheevx.f"
	vuu = *vu;
#line 413 "cheevx.f"
    }
#line 414 "cheevx.f"
    anrm = clanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 415 "cheevx.f"
    if (anrm > 0. && anrm < rmin) {
#line 416 "cheevx.f"
	iscale = 1;
#line 417 "cheevx.f"
	sigma = rmin / anrm;
#line 418 "cheevx.f"
    } else if (anrm > rmax) {
#line 419 "cheevx.f"
	iscale = 1;
#line 420 "cheevx.f"
	sigma = rmax / anrm;
#line 421 "cheevx.f"
    }
#line 422 "cheevx.f"
    if (iscale == 1) {
#line 423 "cheevx.f"
	if (lower) {
#line 424 "cheevx.f"
	    i__1 = *n;
#line 424 "cheevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 425 "cheevx.f"
		i__2 = *n - j + 1;
#line 425 "cheevx.f"
		csscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 426 "cheevx.f"
/* L10: */
#line 426 "cheevx.f"
	    }
#line 427 "cheevx.f"
	} else {
#line 428 "cheevx.f"
	    i__1 = *n;
#line 428 "cheevx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 429 "cheevx.f"
		csscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 430 "cheevx.f"
/* L20: */
#line 430 "cheevx.f"
	    }
#line 431 "cheevx.f"
	}
#line 432 "cheevx.f"
	if (*abstol > 0.) {
#line 432 "cheevx.f"
	    abstll = *abstol * sigma;
#line 432 "cheevx.f"
	}
#line 434 "cheevx.f"
	if (valeig) {
#line 435 "cheevx.f"
	    vll = *vl * sigma;
#line 436 "cheevx.f"
	    vuu = *vu * sigma;
#line 437 "cheevx.f"
	}
#line 438 "cheevx.f"
    }

/*     Call CHETRD to reduce Hermitian matrix to tridiagonal form. */

#line 442 "cheevx.f"
    indd = 1;
#line 443 "cheevx.f"
    inde = indd + *n;
#line 444 "cheevx.f"
    indrwk = inde + *n;
#line 445 "cheevx.f"
    indtau = 1;
#line 446 "cheevx.f"
    indwrk = indtau + *n;
#line 447 "cheevx.f"
    llwork = *lwork - indwrk + 1;
#line 448 "cheevx.f"
    chetrd_(uplo, n, &a[a_offset], lda, &rwork[indd], &rwork[inde], &work[
	    indtau], &work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call SSTERF or CUNGTR and CSTEQR.  If this fails for */
/*     some eigenvalue, then try SSTEBZ. */

#line 455 "cheevx.f"
    test = FALSE_;
#line 456 "cheevx.f"
    if (indeig) {
#line 457 "cheevx.f"
	if (*il == 1 && *iu == *n) {
#line 458 "cheevx.f"
	    test = TRUE_;
#line 459 "cheevx.f"
	}
#line 460 "cheevx.f"
    }
#line 461 "cheevx.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 462 "cheevx.f"
	scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 463 "cheevx.f"
	indee = indrwk + (*n << 1);
#line 464 "cheevx.f"
	if (! wantz) {
#line 465 "cheevx.f"
	    i__1 = *n - 1;
#line 465 "cheevx.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 466 "cheevx.f"
	    ssterf_(n, &w[1], &rwork[indee], info);
#line 467 "cheevx.f"
	} else {
#line 468 "cheevx.f"
	    clacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 469 "cheevx.f"
	    cungtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
#line 471 "cheevx.f"
	    i__1 = *n - 1;
#line 471 "cheevx.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 472 "cheevx.f"
	    csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 474 "cheevx.f"
	    if (*info == 0) {
#line 475 "cheevx.f"
		i__1 = *n;
#line 475 "cheevx.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 476 "cheevx.f"
		    ifail[i__] = 0;
#line 477 "cheevx.f"
/* L30: */
#line 477 "cheevx.f"
		}
#line 478 "cheevx.f"
	    }
#line 479 "cheevx.f"
	}
#line 480 "cheevx.f"
	if (*info == 0) {
#line 481 "cheevx.f"
	    *m = *n;
#line 482 "cheevx.f"
	    goto L40;
#line 483 "cheevx.f"
	}
#line 484 "cheevx.f"
	*info = 0;
#line 485 "cheevx.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */

#line 489 "cheevx.f"
    if (wantz) {
#line 490 "cheevx.f"
	*(unsigned char *)order = 'B';
#line 491 "cheevx.f"
    } else {
#line 492 "cheevx.f"
	*(unsigned char *)order = 'E';
#line 493 "cheevx.f"
    }
#line 494 "cheevx.f"
    indibl = 1;
#line 495 "cheevx.f"
    indisp = indibl + *n;
#line 496 "cheevx.f"
    indiwk = indisp + *n;
#line 497 "cheevx.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 502 "cheevx.f"
    if (wantz) {
#line 503 "cheevx.f"
	cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by CSTEIN. */

#line 510 "cheevx.f"
	cunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwrk], &llwork, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 512 "cheevx.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 516 "cheevx.f"
L40:
#line 517 "cheevx.f"
    if (iscale == 1) {
#line 518 "cheevx.f"
	if (*info == 0) {
#line 519 "cheevx.f"
	    imax = *m;
#line 520 "cheevx.f"
	} else {
#line 521 "cheevx.f"
	    imax = *info - 1;
#line 522 "cheevx.f"
	}
#line 523 "cheevx.f"
	d__1 = 1. / sigma;
#line 523 "cheevx.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 524 "cheevx.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 529 "cheevx.f"
    if (wantz) {
#line 530 "cheevx.f"
	i__1 = *m - 1;
#line 530 "cheevx.f"
	for (j = 1; j <= i__1; ++j) {
#line 531 "cheevx.f"
	    i__ = 0;
#line 532 "cheevx.f"
	    tmp1 = w[j];
#line 533 "cheevx.f"
	    i__2 = *m;
#line 533 "cheevx.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 534 "cheevx.f"
		if (w[jj] < tmp1) {
#line 535 "cheevx.f"
		    i__ = jj;
#line 536 "cheevx.f"
		    tmp1 = w[jj];
#line 537 "cheevx.f"
		}
#line 538 "cheevx.f"
/* L50: */
#line 538 "cheevx.f"
	    }

#line 540 "cheevx.f"
	    if (i__ != 0) {
#line 541 "cheevx.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 542 "cheevx.f"
		w[i__] = w[j];
#line 543 "cheevx.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 544 "cheevx.f"
		w[j] = tmp1;
#line 545 "cheevx.f"
		iwork[indibl + j - 1] = itmp1;
#line 546 "cheevx.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 547 "cheevx.f"
		if (*info != 0) {
#line 548 "cheevx.f"
		    itmp1 = ifail[i__];
#line 549 "cheevx.f"
		    ifail[i__] = ifail[j];
#line 550 "cheevx.f"
		    ifail[j] = itmp1;
#line 551 "cheevx.f"
		}
#line 552 "cheevx.f"
	    }
#line 553 "cheevx.f"
/* L60: */
#line 553 "cheevx.f"
	}
#line 554 "cheevx.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 558 "cheevx.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 560 "cheevx.f"
    return 0;

/*     End of CHEEVX */

} /* cheevx_ */

