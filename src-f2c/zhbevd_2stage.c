#line 1 "zhbevd_2stage.f"
/* zhbevd_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "zhbevd_2stage.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__4 = 4;
static doublereal c_b23 = 1.;
static integer c__1 = 1;

/* > \brief <b> ZHBEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 OTHER matrices</b> */

/*  @precisions fortran z -> s d c */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHBEVD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbevd_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbevd_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbevd_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/*                                 WORK, LWORK, RWORK, LRWORK, IWORK, */
/*                                 LIWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBEVD_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian band matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  If eigenvectors are desired, it */
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
/* >                  Not available in this release. */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB, N) */
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
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK. LWORK >= 1, when N <= 1; */
/* >          otherwise */
/* >          If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* >                                   LWORK = MAX(1, dimension) where */
/* >                                   dimension = (2KD+1)*N + KD*NTHREADS */
/* >                                   where KD is the size of the band. */
/* >                                   NTHREADS is the number of threads used when */
/* >                                   openMP compilation is enabled, otherwise =1. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available. */
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
/* >          RWORK is DOUBLE PRECISION array, */
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

/* > \date November 2017 */

/* > \ingroup complex16OTHEReigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  All details about the 2stage techniques are available in: */
/* > */
/* >  Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* >  Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* >  using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* >  of 2011 International Conference for High Performance Computing, */
/* >  Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* >  Article 8 , 11 pages. */
/* >  http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* >  A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* >  An improved parallel singular value algorithm and its implementation */
/* >  for multicore hardware, In Proceedings of 2013 International Conference */
/* >  for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* >  Denver, Colorado, USA, 2013. */
/* >  Article 90, 12 pages. */
/* >  http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* >  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* >  A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* >  calculations based on fine-grained memory aware tasks. */
/* >  International Journal of High Performance Computing Applications. */
/* >  Volume 28 Issue 2, Pages 196-209, May 2014. */
/* >  http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zhbevd_2stage__(char *jobz, char *uplo, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublereal *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork,
	 doublereal *rwork, integer *lrwork, integer *iwork, integer *liwork, 
	integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ib;
    static doublereal eps;
    static integer inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static doublereal anrm;
    static integer imax;
    extern /* Subroutine */ int zhetrd_hb2st__(char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal rmin, rmax;
    static integer llwk2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo, indwk, lhtrd;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer lwmin;
    static logical lower;
    static integer lwtrd, llrwk;
    static logical wantz;
    static integer indwk2;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    static doublereal safmin;
    extern doublereal zlanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), zlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen), zstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen);
    static integer indrwk, liwmin;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer lrwmin, llwork;
    static doublereal smlnum;
    static logical lquery;
    static integer indhous;



/*  -- LAPACK driver routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 314 "zhbevd_2stage.f"
    /* Parameter adjustments */
#line 314 "zhbevd_2stage.f"
    ab_dim1 = *ldab;
#line 314 "zhbevd_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 314 "zhbevd_2stage.f"
    ab -= ab_offset;
#line 314 "zhbevd_2stage.f"
    --w;
#line 314 "zhbevd_2stage.f"
    z_dim1 = *ldz;
#line 314 "zhbevd_2stage.f"
    z_offset = 1 + z_dim1;
#line 314 "zhbevd_2stage.f"
    z__ -= z_offset;
#line 314 "zhbevd_2stage.f"
    --work;
#line 314 "zhbevd_2stage.f"
    --rwork;
#line 314 "zhbevd_2stage.f"
    --iwork;
#line 314 "zhbevd_2stage.f"

#line 314 "zhbevd_2stage.f"
    /* Function Body */
#line 314 "zhbevd_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 315 "zhbevd_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 316 "zhbevd_2stage.f"
    lquery = *lwork == -1 || *liwork == -1 || *lrwork == -1;

#line 318 "zhbevd_2stage.f"
    *info = 0;
#line 319 "zhbevd_2stage.f"
    if (*n <= 1) {
#line 320 "zhbevd_2stage.f"
	lwmin = 1;
#line 321 "zhbevd_2stage.f"
	lrwmin = 1;
#line 322 "zhbevd_2stage.f"
	liwmin = 1;
#line 323 "zhbevd_2stage.f"
    } else {
#line 324 "zhbevd_2stage.f"
	ib = ilaenv2stage_(&c__2, "ZHETRD_HB2ST", jobz, n, kd, &c_n1, &c_n1, (
		ftnlen)12, (ftnlen)1);
#line 325 "zhbevd_2stage.f"
	lhtrd = ilaenv2stage_(&c__3, "ZHETRD_HB2ST", jobz, n, kd, &ib, &c_n1, 
		(ftnlen)12, (ftnlen)1);
#line 326 "zhbevd_2stage.f"
	lwtrd = ilaenv2stage_(&c__4, "ZHETRD_HB2ST", jobz, n, kd, &ib, &c_n1, 
		(ftnlen)12, (ftnlen)1);
#line 327 "zhbevd_2stage.f"
	if (wantz) {
/* Computing 2nd power */
#line 328 "zhbevd_2stage.f"
	    i__1 = *n;
#line 328 "zhbevd_2stage.f"
	    lwmin = i__1 * i__1 << 1;
/* Computing 2nd power */
#line 329 "zhbevd_2stage.f"
	    i__1 = *n;
#line 329 "zhbevd_2stage.f"
	    lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 330 "zhbevd_2stage.f"
	    liwmin = *n * 5 + 3;
#line 331 "zhbevd_2stage.f"
	} else {
/* Computing MAX */
#line 332 "zhbevd_2stage.f"
	    i__1 = *n, i__2 = lhtrd + lwtrd;
#line 332 "zhbevd_2stage.f"
	    lwmin = max(i__1,i__2);
#line 333 "zhbevd_2stage.f"
	    lrwmin = *n;
#line 334 "zhbevd_2stage.f"
	    liwmin = 1;
#line 335 "zhbevd_2stage.f"
	}
#line 336 "zhbevd_2stage.f"
    }
#line 337 "zhbevd_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 338 "zhbevd_2stage.f"
	*info = -1;
#line 339 "zhbevd_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 340 "zhbevd_2stage.f"
	*info = -2;
#line 341 "zhbevd_2stage.f"
    } else if (*n < 0) {
#line 342 "zhbevd_2stage.f"
	*info = -3;
#line 343 "zhbevd_2stage.f"
    } else if (*kd < 0) {
#line 344 "zhbevd_2stage.f"
	*info = -4;
#line 345 "zhbevd_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 346 "zhbevd_2stage.f"
	*info = -6;
#line 347 "zhbevd_2stage.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 348 "zhbevd_2stage.f"
	*info = -9;
#line 349 "zhbevd_2stage.f"
    }

#line 351 "zhbevd_2stage.f"
    if (*info == 0) {
#line 352 "zhbevd_2stage.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 353 "zhbevd_2stage.f"
	rwork[1] = (doublereal) lrwmin;
#line 354 "zhbevd_2stage.f"
	iwork[1] = liwmin;

#line 356 "zhbevd_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 357 "zhbevd_2stage.f"
	    *info = -11;
#line 358 "zhbevd_2stage.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 359 "zhbevd_2stage.f"
	    *info = -13;
#line 360 "zhbevd_2stage.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 361 "zhbevd_2stage.f"
	    *info = -15;
#line 362 "zhbevd_2stage.f"
	}
#line 363 "zhbevd_2stage.f"
    }

#line 365 "zhbevd_2stage.f"
    if (*info != 0) {
#line 366 "zhbevd_2stage.f"
	i__1 = -(*info);
#line 366 "zhbevd_2stage.f"
	xerbla_("ZHBEVD_2STAGE", &i__1, (ftnlen)13);
#line 367 "zhbevd_2stage.f"
	return 0;
#line 368 "zhbevd_2stage.f"
    } else if (lquery) {
#line 369 "zhbevd_2stage.f"
	return 0;
#line 370 "zhbevd_2stage.f"
    }

/*     Quick return if possible */

#line 374 "zhbevd_2stage.f"
    if (*n == 0) {
#line 374 "zhbevd_2stage.f"
	return 0;
#line 374 "zhbevd_2stage.f"
    }

#line 377 "zhbevd_2stage.f"
    if (*n == 1) {
#line 378 "zhbevd_2stage.f"
	i__1 = ab_dim1 + 1;
#line 378 "zhbevd_2stage.f"
	w[1] = ab[i__1].r;
#line 379 "zhbevd_2stage.f"
	if (wantz) {
#line 379 "zhbevd_2stage.f"
	    i__1 = z_dim1 + 1;
#line 379 "zhbevd_2stage.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 379 "zhbevd_2stage.f"
	}
#line 381 "zhbevd_2stage.f"
	return 0;
#line 382 "zhbevd_2stage.f"
    }

/*     Get machine constants. */

#line 386 "zhbevd_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 387 "zhbevd_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 388 "zhbevd_2stage.f"
    smlnum = safmin / eps;
#line 389 "zhbevd_2stage.f"
    bignum = 1. / smlnum;
#line 390 "zhbevd_2stage.f"
    rmin = sqrt(smlnum);
#line 391 "zhbevd_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 395 "zhbevd_2stage.f"
    anrm = zlanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1, (ftnlen)1);
#line 396 "zhbevd_2stage.f"
    iscale = 0;
#line 397 "zhbevd_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 398 "zhbevd_2stage.f"
	iscale = 1;
#line 399 "zhbevd_2stage.f"
	sigma = rmin / anrm;
#line 400 "zhbevd_2stage.f"
    } else if (anrm > rmax) {
#line 401 "zhbevd_2stage.f"
	iscale = 1;
#line 402 "zhbevd_2stage.f"
	sigma = rmax / anrm;
#line 403 "zhbevd_2stage.f"
    }
#line 404 "zhbevd_2stage.f"
    if (iscale == 1) {
#line 405 "zhbevd_2stage.f"
	if (lower) {
#line 406 "zhbevd_2stage.f"
	    zlascl_("B", kd, kd, &c_b23, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 407 "zhbevd_2stage.f"
	} else {
#line 408 "zhbevd_2stage.f"
	    zlascl_("Q", kd, kd, &c_b23, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 409 "zhbevd_2stage.f"
	}
#line 410 "zhbevd_2stage.f"
    }

/*     Call ZHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form. */

#line 414 "zhbevd_2stage.f"
    inde = 1;
#line 415 "zhbevd_2stage.f"
    indrwk = inde + *n;
#line 416 "zhbevd_2stage.f"
    llrwk = *lrwork - indrwk + 1;
#line 417 "zhbevd_2stage.f"
    indhous = 1;
#line 418 "zhbevd_2stage.f"
    indwk = indhous + lhtrd;
#line 419 "zhbevd_2stage.f"
    llwork = *lwork - indwk + 1;
#line 420 "zhbevd_2stage.f"
    indwk2 = indwk + *n * *n;
#line 421 "zhbevd_2stage.f"
    llwk2 = *lwork - indwk2 + 1;

#line 423 "zhbevd_2stage.f"
    zhetrd_hb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &
	    rwork[inde], &work[indhous], &lhtrd, &work[indwk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC. */

#line 429 "zhbevd_2stage.f"
    if (! wantz) {
#line 430 "zhbevd_2stage.f"
	dsterf_(n, &w[1], &rwork[inde], info);
#line 431 "zhbevd_2stage.f"
    } else {
#line 432 "zhbevd_2stage.f"
	zstedc_("I", n, &w[1], &rwork[inde], &work[1], n, &work[indwk2], &
		llwk2, &rwork[indrwk], &llrwk, &iwork[1], liwork, info, (
		ftnlen)1);
#line 435 "zhbevd_2stage.f"
	zgemm_("N", "N", n, n, n, &c_b2, &z__[z_offset], ldz, &work[1], n, &
		c_b1, &work[indwk2], n, (ftnlen)1, (ftnlen)1);
#line 437 "zhbevd_2stage.f"
	zlacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz, (ftnlen)1);
#line 438 "zhbevd_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 442 "zhbevd_2stage.f"
    if (iscale == 1) {
#line 443 "zhbevd_2stage.f"
	if (*info == 0) {
#line 444 "zhbevd_2stage.f"
	    imax = *n;
#line 445 "zhbevd_2stage.f"
	} else {
#line 446 "zhbevd_2stage.f"
	    imax = *info - 1;
#line 447 "zhbevd_2stage.f"
	}
#line 448 "zhbevd_2stage.f"
	d__1 = 1. / sigma;
#line 448 "zhbevd_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 449 "zhbevd_2stage.f"
    }

#line 451 "zhbevd_2stage.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 452 "zhbevd_2stage.f"
    rwork[1] = (doublereal) lrwmin;
#line 453 "zhbevd_2stage.f"
    iwork[1] = liwmin;
#line 454 "zhbevd_2stage.f"
    return 0;

/*     End of ZHBEVD_2STAGE */

} /* zhbevd_2stage__ */

