#line 1 "chbevd.f"
/* chbevd.f -- translated by f2c (version 20100827).
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

#line 1 "chbevd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static doublereal c_b13 = 1.;
static integer c__1 = 1;

/* > \brief <b> CHBEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/*                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEVD computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian band matrix A.  If eigenvectors are desired, it */
/* > uses a divide and conquer algorithm. */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, AB is overwritten by values generated during the */
/* >          reduction to tridiagonal form.  If UPLO = 'U', the first */
/* >          superdiagonal and the diagonal of the tridiagonal matrix T */
/* >          are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
/* >          the diagonal and first subdiagonal of T are returned in the */
/* >          first two rows of AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD + 1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= 1,               LWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK must be at least N. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N**2. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal sizes of the WORK, RWORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, */
/* >                                         dimension (LRWORK) */
/* >          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of array RWORK. */
/* >          If N <= 1,               LRWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LRWORK must be at least N. */
/* >          If JOBZ = 'V' and N > 1, LRWORK must be at least */
/* >                        1 + 5*N + 2*N**2. */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of array IWORK. */
/* >          If JOBZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* >          If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N . */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
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

/* > \ingroup complexOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int chbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static integer llwk2;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lwmin;
    static logical lower;
    static integer llrwk;
    static logical wantz;
    static integer indwk2;
    extern doublereal clanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), cstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, ftnlen), chbtrd_(char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indwrk, liwmin;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer lrwmin;
    static doublereal smlnum;
    static logical lquery;


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

#line 265 "chbevd.f"
    /* Parameter adjustments */
#line 265 "chbevd.f"
    ab_dim1 = *ldab;
#line 265 "chbevd.f"
    ab_offset = 1 + ab_dim1;
#line 265 "chbevd.f"
    ab -= ab_offset;
#line 265 "chbevd.f"
    --w;
#line 265 "chbevd.f"
    z_dim1 = *ldz;
#line 265 "chbevd.f"
    z_offset = 1 + z_dim1;
#line 265 "chbevd.f"
    z__ -= z_offset;
#line 265 "chbevd.f"
    --work;
#line 265 "chbevd.f"
    --rwork;
#line 265 "chbevd.f"
    --iwork;
#line 265 "chbevd.f"

#line 265 "chbevd.f"
    /* Function Body */
#line 265 "chbevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 266 "chbevd.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 267 "chbevd.f"
    lquery = *lwork == -1 || *liwork == -1 || *lrwork == -1;

#line 269 "chbevd.f"
    *info = 0;
#line 270 "chbevd.f"
    if (*n <= 1) {
#line 271 "chbevd.f"
	lwmin = 1;
#line 272 "chbevd.f"
	lrwmin = 1;
#line 273 "chbevd.f"
	liwmin = 1;
#line 274 "chbevd.f"
    } else {
#line 275 "chbevd.f"
	if (wantz) {
/* Computing 2nd power */
#line 276 "chbevd.f"
	    i__1 = *n;
#line 276 "chbevd.f"
	    lwmin = i__1 * i__1 << 1;
/* Computing 2nd power */
#line 277 "chbevd.f"
	    i__1 = *n;
#line 277 "chbevd.f"
	    lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 278 "chbevd.f"
	    liwmin = *n * 5 + 3;
#line 279 "chbevd.f"
	} else {
#line 280 "chbevd.f"
	    lwmin = *n;
#line 281 "chbevd.f"
	    lrwmin = *n;
#line 282 "chbevd.f"
	    liwmin = 1;
#line 283 "chbevd.f"
	}
#line 284 "chbevd.f"
    }
#line 285 "chbevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 286 "chbevd.f"
	*info = -1;
#line 287 "chbevd.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 288 "chbevd.f"
	*info = -2;
#line 289 "chbevd.f"
    } else if (*n < 0) {
#line 290 "chbevd.f"
	*info = -3;
#line 291 "chbevd.f"
    } else if (*kd < 0) {
#line 292 "chbevd.f"
	*info = -4;
#line 293 "chbevd.f"
    } else if (*ldab < *kd + 1) {
#line 294 "chbevd.f"
	*info = -6;
#line 295 "chbevd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 296 "chbevd.f"
	*info = -9;
#line 297 "chbevd.f"
    }

#line 299 "chbevd.f"
    if (*info == 0) {
#line 300 "chbevd.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 301 "chbevd.f"
	rwork[1] = (doublereal) lrwmin;
#line 302 "chbevd.f"
	iwork[1] = liwmin;

#line 304 "chbevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 305 "chbevd.f"
	    *info = -11;
#line 306 "chbevd.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 307 "chbevd.f"
	    *info = -13;
#line 308 "chbevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 309 "chbevd.f"
	    *info = -15;
#line 310 "chbevd.f"
	}
#line 311 "chbevd.f"
    }

#line 313 "chbevd.f"
    if (*info != 0) {
#line 314 "chbevd.f"
	i__1 = -(*info);
#line 314 "chbevd.f"
	xerbla_("CHBEVD", &i__1, (ftnlen)6);
#line 315 "chbevd.f"
	return 0;
#line 316 "chbevd.f"
    } else if (lquery) {
#line 317 "chbevd.f"
	return 0;
#line 318 "chbevd.f"
    }

/*     Quick return if possible */

#line 322 "chbevd.f"
    if (*n == 0) {
#line 322 "chbevd.f"
	return 0;
#line 322 "chbevd.f"
    }

#line 325 "chbevd.f"
    if (*n == 1) {
#line 326 "chbevd.f"
	i__1 = ab_dim1 + 1;
#line 326 "chbevd.f"
	w[1] = ab[i__1].r;
#line 327 "chbevd.f"
	if (wantz) {
#line 327 "chbevd.f"
	    i__1 = z_dim1 + 1;
#line 327 "chbevd.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 327 "chbevd.f"
	}
#line 329 "chbevd.f"
	return 0;
#line 330 "chbevd.f"
    }

/*     Get machine constants. */

#line 334 "chbevd.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 335 "chbevd.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 336 "chbevd.f"
    smlnum = safmin / eps;
#line 337 "chbevd.f"
    bignum = 1. / smlnum;
#line 338 "chbevd.f"
    rmin = sqrt(smlnum);
#line 339 "chbevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 343 "chbevd.f"
    anrm = clanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1, (ftnlen)1);
#line 344 "chbevd.f"
    iscale = 0;
#line 345 "chbevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 346 "chbevd.f"
	iscale = 1;
#line 347 "chbevd.f"
	sigma = rmin / anrm;
#line 348 "chbevd.f"
    } else if (anrm > rmax) {
#line 349 "chbevd.f"
	iscale = 1;
#line 350 "chbevd.f"
	sigma = rmax / anrm;
#line 351 "chbevd.f"
    }
#line 352 "chbevd.f"
    if (iscale == 1) {
#line 353 "chbevd.f"
	if (lower) {
#line 354 "chbevd.f"
	    clascl_("B", kd, kd, &c_b13, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 355 "chbevd.f"
	} else {
#line 356 "chbevd.f"
	    clascl_("Q", kd, kd, &c_b13, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 357 "chbevd.f"
	}
#line 358 "chbevd.f"
    }

/*     Call CHBTRD to reduce Hermitian band matrix to tridiagonal form. */

#line 362 "chbevd.f"
    inde = 1;
#line 363 "chbevd.f"
    indwrk = inde + *n;
#line 364 "chbevd.f"
    indwk2 = *n * *n + 1;
#line 365 "chbevd.f"
    llwk2 = *lwork - indwk2 + 1;
#line 366 "chbevd.f"
    llrwk = *lrwork - indwrk + 1;
#line 367 "chbevd.f"
    chbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &rwork[inde], &
	    z__[z_offset], ldz, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEDC. */

#line 372 "chbevd.f"
    if (! wantz) {
#line 373 "chbevd.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 374 "chbevd.f"
    } else {
#line 375 "chbevd.f"
	cstedc_("I", n, &w[1], &rwork[inde], &work[1], n, &work[indwk2], &
		llwk2, &rwork[indwrk], &llrwk, &iwork[1], liwork, info, (
		ftnlen)1);
#line 378 "chbevd.f"
	cgemm_("N", "N", n, n, n, &c_b2, &z__[z_offset], ldz, &work[1], n, &
		c_b1, &work[indwk2], n, (ftnlen)1, (ftnlen)1);
#line 380 "chbevd.f"
	clacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz, (ftnlen)1);
#line 381 "chbevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 385 "chbevd.f"
    if (iscale == 1) {
#line 386 "chbevd.f"
	if (*info == 0) {
#line 387 "chbevd.f"
	    imax = *n;
#line 388 "chbevd.f"
	} else {
#line 389 "chbevd.f"
	    imax = *info - 1;
#line 390 "chbevd.f"
	}
#line 391 "chbevd.f"
	d__1 = 1. / sigma;
#line 391 "chbevd.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 392 "chbevd.f"
    }

#line 394 "chbevd.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 395 "chbevd.f"
    rwork[1] = (doublereal) lrwmin;
#line 396 "chbevd.f"
    iwork[1] = liwmin;
#line 397 "chbevd.f"
    return 0;

/*     End of CHBEVD */

} /* chbevd_ */

