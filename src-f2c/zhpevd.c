#line 1 "zhpevd.f"
/* zhpevd.f -- translated by f2c (version 20100827).
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

#line 1 "zhpevd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> ZHPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, */
/*                          RWORK, LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPEVD computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian matrix A in packed storage.  If eigenvectors are */
/* > desired, it uses a divide and conquer algorithm. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
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
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* >          eigenvectors of the matrix A, with the i-th column of Z */
/* >          holding the eigenvector associated with W(i). */
/* >          If JOBZ = 'N', then Z is not referenced. */
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
/* >          On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of array WORK. */
/* >          If N <= 1,               LWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK must be at least N. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the required sizes of the WORK, RWORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, */
/* >                                         dimension (LRWORK) */
/* >          On exit, if INFO = 0, RWORK(1) returns the required LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of array RWORK. */
/* >          If N <= 1,               LRWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LRWORK must be at least N. */
/* >          If JOBZ = 'V' and N > 1, LRWORK must be at least */
/* >                    1 + 5*N + 2*N**2. */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the required LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of array IWORK. */
/* >          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the algorithm failed to converge; i */
/* >                off-diagonal elements of an intermediate tridiagonal */
/* >                form did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16OTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int zhpevd_(char *jobz, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info, ftnlen 
	jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo, lwmin, llrwk, llwrk;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal zlanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int zstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen);
    static integer indrwk, indwrk, liwmin, lrwmin;
    static doublereal smlnum;
    extern /* Subroutine */ int zhptrd_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int zupmtr_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen);


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

#line 250 "zhpevd.f"
    /* Parameter adjustments */
#line 250 "zhpevd.f"
    --ap;
#line 250 "zhpevd.f"
    --w;
#line 250 "zhpevd.f"
    z_dim1 = *ldz;
#line 250 "zhpevd.f"
    z_offset = 1 + z_dim1;
#line 250 "zhpevd.f"
    z__ -= z_offset;
#line 250 "zhpevd.f"
    --work;
#line 250 "zhpevd.f"
    --rwork;
#line 250 "zhpevd.f"
    --iwork;
#line 250 "zhpevd.f"

#line 250 "zhpevd.f"
    /* Function Body */
#line 250 "zhpevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 251 "zhpevd.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 253 "zhpevd.f"
    *info = 0;
#line 254 "zhpevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 255 "zhpevd.f"
	*info = -1;
#line 256 "zhpevd.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 258 "zhpevd.f"
	*info = -2;
#line 259 "zhpevd.f"
    } else if (*n < 0) {
#line 260 "zhpevd.f"
	*info = -3;
#line 261 "zhpevd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 262 "zhpevd.f"
	*info = -7;
#line 263 "zhpevd.f"
    }

#line 265 "zhpevd.f"
    if (*info == 0) {
#line 266 "zhpevd.f"
	if (*n <= 1) {
#line 267 "zhpevd.f"
	    lwmin = 1;
#line 268 "zhpevd.f"
	    liwmin = 1;
#line 269 "zhpevd.f"
	    lrwmin = 1;
#line 270 "zhpevd.f"
	} else {
#line 271 "zhpevd.f"
	    if (wantz) {
#line 272 "zhpevd.f"
		lwmin = *n << 1;
/* Computing 2nd power */
#line 273 "zhpevd.f"
		i__1 = *n;
#line 273 "zhpevd.f"
		lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 274 "zhpevd.f"
		liwmin = *n * 5 + 3;
#line 275 "zhpevd.f"
	    } else {
#line 276 "zhpevd.f"
		lwmin = *n;
#line 277 "zhpevd.f"
		lrwmin = *n;
#line 278 "zhpevd.f"
		liwmin = 1;
#line 279 "zhpevd.f"
	    }
#line 280 "zhpevd.f"
	}
#line 281 "zhpevd.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 282 "zhpevd.f"
	rwork[1] = (doublereal) lrwmin;
#line 283 "zhpevd.f"
	iwork[1] = liwmin;

#line 285 "zhpevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 286 "zhpevd.f"
	    *info = -9;
#line 287 "zhpevd.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 288 "zhpevd.f"
	    *info = -11;
#line 289 "zhpevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 290 "zhpevd.f"
	    *info = -13;
#line 291 "zhpevd.f"
	}
#line 292 "zhpevd.f"
    }

#line 294 "zhpevd.f"
    if (*info != 0) {
#line 295 "zhpevd.f"
	i__1 = -(*info);
#line 295 "zhpevd.f"
	xerbla_("ZHPEVD", &i__1, (ftnlen)6);
#line 296 "zhpevd.f"
	return 0;
#line 297 "zhpevd.f"
    } else if (lquery) {
#line 298 "zhpevd.f"
	return 0;
#line 299 "zhpevd.f"
    }

/*     Quick return if possible */

#line 303 "zhpevd.f"
    if (*n == 0) {
#line 303 "zhpevd.f"
	return 0;
#line 303 "zhpevd.f"
    }

#line 306 "zhpevd.f"
    if (*n == 1) {
#line 307 "zhpevd.f"
	w[1] = ap[1].r;
#line 308 "zhpevd.f"
	if (wantz) {
#line 308 "zhpevd.f"
	    i__1 = z_dim1 + 1;
#line 308 "zhpevd.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 308 "zhpevd.f"
	}
#line 310 "zhpevd.f"
	return 0;
#line 311 "zhpevd.f"
    }

/*     Get machine constants. */

#line 315 "zhpevd.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 316 "zhpevd.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 317 "zhpevd.f"
    smlnum = safmin / eps;
#line 318 "zhpevd.f"
    bignum = 1. / smlnum;
#line 319 "zhpevd.f"
    rmin = sqrt(smlnum);
#line 320 "zhpevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 324 "zhpevd.f"
    anrm = zlanhp_("M", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);
#line 325 "zhpevd.f"
    iscale = 0;
#line 326 "zhpevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 327 "zhpevd.f"
	iscale = 1;
#line 328 "zhpevd.f"
	sigma = rmin / anrm;
#line 329 "zhpevd.f"
    } else if (anrm > rmax) {
#line 330 "zhpevd.f"
	iscale = 1;
#line 331 "zhpevd.f"
	sigma = rmax / anrm;
#line 332 "zhpevd.f"
    }
#line 333 "zhpevd.f"
    if (iscale == 1) {
#line 334 "zhpevd.f"
	i__1 = *n * (*n + 1) / 2;
#line 334 "zhpevd.f"
	zdscal_(&i__1, &sigma, &ap[1], &c__1);
#line 335 "zhpevd.f"
    }

/*     Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form. */

#line 339 "zhpevd.f"
    inde = 1;
#line 340 "zhpevd.f"
    indtau = 1;
#line 341 "zhpevd.f"
    indrwk = inde + *n;
#line 342 "zhpevd.f"
    indwrk = indtau + *n;
#line 343 "zhpevd.f"
    llwrk = *lwork - indwrk + 1;
#line 344 "zhpevd.f"
    llrwk = *lrwork - indrwk + 1;
#line 345 "zhpevd.f"
    zhptrd_(uplo, n, &ap[1], &w[1], &rwork[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     ZUPGTR to generate the orthogonal matrix, then call ZSTEDC. */

#line 351 "zhpevd.f"
    if (! wantz) {
#line 352 "zhpevd.f"
	dsterf_(n, &w[1], &rwork[inde], info);
#line 353 "zhpevd.f"
    } else {
#line 354 "zhpevd.f"
	zstedc_("I", n, &w[1], &rwork[inde], &z__[z_offset], ldz, &work[
		indwrk], &llwrk, &rwork[indrwk], &llrwk, &iwork[1], liwork, 
		info, (ftnlen)1);
#line 357 "zhpevd.f"
	zupmtr_("L", uplo, "N", n, n, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 359 "zhpevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 363 "zhpevd.f"
    if (iscale == 1) {
#line 364 "zhpevd.f"
	if (*info == 0) {
#line 365 "zhpevd.f"
	    imax = *n;
#line 366 "zhpevd.f"
	} else {
#line 367 "zhpevd.f"
	    imax = *info - 1;
#line 368 "zhpevd.f"
	}
#line 369 "zhpevd.f"
	d__1 = 1. / sigma;
#line 369 "zhpevd.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 370 "zhpevd.f"
    }

#line 372 "zhpevd.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 373 "zhpevd.f"
    rwork[1] = (doublereal) lrwmin;
#line 374 "zhpevd.f"
    iwork[1] = liwmin;
#line 375 "zhpevd.f"
    return 0;

/*     End of ZHPEVD */

} /* zhpevd_ */

