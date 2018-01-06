#line 1 "cheevd_2stage.f"
/* cheevd_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "cheevd_2stage.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__0 = 0;
static doublereal c_b28 = 1.;

/* > \brief <b> CHEEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 HE matrices</b> */

/*  @generated from zheevd_2stage.f, fortran z -> c, Sat Nov  5 23:18:14 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEEVD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheevd_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheevd_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheevd_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEEVD_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, */
/*                          RWORK, LRWORK, IWORK, LIWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEEVD_2STAGE computes all eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  If eigenvectors are desired, it uses a */
/* > divide and conquer algorithm. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* >          On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* >          orthonormal eigenvectors of the matrix A. */
/* >          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') */
/* >          or the upper triangle (if UPLO='U') of A, including the */
/* >          diagonal, is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
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
/* >          If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* >                                   LWORK = MAX(1, dimension) where */
/* >                                   dimension = max(stage1,stage2) + (KD+1)*N + N+1 */
/* >                                             = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                                               + max(2*KD*KD, KD*NTHREADS) */
/* >                                               + (KD+1)*N + N+1 */
/* >                                   where KD is the blocking size of the reduction, */
/* >                                   FACTOPTNB is the blocking used by the QR or LQ */
/* >                                   algorithm, usually FACTOPTNB=128 is a good choice */
/* >                                   NTHREADS is the number of threads used when */
/* >                                   openMP compilation is enabled, otherwise =1. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N + N**2 */
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
/* >          The dimension of the array RWORK. */
/* >          If N <= 1,                LRWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LRWORK must be at least N. */
/* >          If JOBZ  = 'V' and N > 1, LRWORK must be at least */
/* >                         1 + 5*N + 2*N**2. */
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
/* >          The dimension of the array IWORK. */
/* >          If N <= 1,                LIWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
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
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed */
/* >                to converge; i off-diagonal elements of an intermediate */
/* >                tridiagonal form did not converge to zero; */
/* >                if INFO = i and JOBZ = 'V', then the algorithm failed */
/* >                to compute an eigenvalue while working on the submatrix */
/* >                lying in rows and columns INFO/(N+1) through */
/* >                mod(INFO,N+1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complexHEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  Modified description of INFO. Sven, 16 Feb 05. */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* > */
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
/* Subroutine */ int cheevd_2stage__(char *jobz, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ib, kd;
    static doublereal eps;
    static integer inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    extern /* Subroutine */ int chetrd_2stage__(char *, char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lhtrd, lwmin;
    static logical lower;
    static integer llrwk, lwtrd;
    static logical wantz;
    static integer indwk2, llwrk2;
    extern doublereal clanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), cstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau, indrwk, indwrk, liwmin;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer lrwmin;
    extern /* Subroutine */ int cunmtr_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen, ftnlen);
    static integer llwork;
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

#line 309 "cheevd_2stage.f"
    /* Parameter adjustments */
#line 309 "cheevd_2stage.f"
    a_dim1 = *lda;
#line 309 "cheevd_2stage.f"
    a_offset = 1 + a_dim1;
#line 309 "cheevd_2stage.f"
    a -= a_offset;
#line 309 "cheevd_2stage.f"
    --w;
#line 309 "cheevd_2stage.f"
    --work;
#line 309 "cheevd_2stage.f"
    --rwork;
#line 309 "cheevd_2stage.f"
    --iwork;
#line 309 "cheevd_2stage.f"

#line 309 "cheevd_2stage.f"
    /* Function Body */
#line 309 "cheevd_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 310 "cheevd_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 311 "cheevd_2stage.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 313 "cheevd_2stage.f"
    *info = 0;
#line 314 "cheevd_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 315 "cheevd_2stage.f"
	*info = -1;
#line 316 "cheevd_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 317 "cheevd_2stage.f"
	*info = -2;
#line 318 "cheevd_2stage.f"
    } else if (*n < 0) {
#line 319 "cheevd_2stage.f"
	*info = -3;
#line 320 "cheevd_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 321 "cheevd_2stage.f"
	*info = -5;
#line 322 "cheevd_2stage.f"
    }

#line 324 "cheevd_2stage.f"
    if (*info == 0) {
#line 325 "cheevd_2stage.f"
	if (*n <= 1) {
#line 326 "cheevd_2stage.f"
	    lwmin = 1;
#line 327 "cheevd_2stage.f"
	    lrwmin = 1;
#line 328 "cheevd_2stage.f"
	    liwmin = 1;
#line 329 "cheevd_2stage.f"
	} else {
#line 330 "cheevd_2stage.f"
	    kd = ilaenv2stage_(&c__1, "CHETRD_2STAGE", jobz, n, &c_n1, &c_n1, 
		    &c_n1, (ftnlen)13, (ftnlen)1);
#line 332 "cheevd_2stage.f"
	    ib = ilaenv2stage_(&c__2, "CHETRD_2STAGE", jobz, n, &kd, &c_n1, &
		    c_n1, (ftnlen)13, (ftnlen)1);
#line 334 "cheevd_2stage.f"
	    lhtrd = ilaenv2stage_(&c__3, "CHETRD_2STAGE", jobz, n, &kd, &ib, &
		    c_n1, (ftnlen)13, (ftnlen)1);
#line 336 "cheevd_2stage.f"
	    lwtrd = ilaenv2stage_(&c__4, "CHETRD_2STAGE", jobz, n, &kd, &ib, &
		    c_n1, (ftnlen)13, (ftnlen)1);
#line 338 "cheevd_2stage.f"
	    if (wantz) {
#line 339 "cheevd_2stage.f"
		lwmin = (*n << 1) + *n * *n;
/* Computing 2nd power */
#line 340 "cheevd_2stage.f"
		i__1 = *n;
#line 340 "cheevd_2stage.f"
		lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 341 "cheevd_2stage.f"
		liwmin = *n * 5 + 3;
#line 342 "cheevd_2stage.f"
	    } else {
#line 343 "cheevd_2stage.f"
		lwmin = *n + 1 + lhtrd + lwtrd;
#line 344 "cheevd_2stage.f"
		lrwmin = *n;
#line 345 "cheevd_2stage.f"
		liwmin = 1;
#line 346 "cheevd_2stage.f"
	    }
#line 347 "cheevd_2stage.f"
	}
#line 348 "cheevd_2stage.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 349 "cheevd_2stage.f"
	rwork[1] = (doublereal) lrwmin;
#line 350 "cheevd_2stage.f"
	iwork[1] = liwmin;

#line 352 "cheevd_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 353 "cheevd_2stage.f"
	    *info = -8;
#line 354 "cheevd_2stage.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 355 "cheevd_2stage.f"
	    *info = -10;
#line 356 "cheevd_2stage.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 357 "cheevd_2stage.f"
	    *info = -12;
#line 358 "cheevd_2stage.f"
	}
#line 359 "cheevd_2stage.f"
    }

#line 361 "cheevd_2stage.f"
    if (*info != 0) {
#line 362 "cheevd_2stage.f"
	i__1 = -(*info);
#line 362 "cheevd_2stage.f"
	xerbla_("CHEEVD_2STAGE", &i__1, (ftnlen)13);
#line 363 "cheevd_2stage.f"
	return 0;
#line 364 "cheevd_2stage.f"
    } else if (lquery) {
#line 365 "cheevd_2stage.f"
	return 0;
#line 366 "cheevd_2stage.f"
    }

/*     Quick return if possible */

#line 370 "cheevd_2stage.f"
    if (*n == 0) {
#line 370 "cheevd_2stage.f"
	return 0;
#line 370 "cheevd_2stage.f"
    }

#line 373 "cheevd_2stage.f"
    if (*n == 1) {
#line 374 "cheevd_2stage.f"
	i__1 = a_dim1 + 1;
#line 374 "cheevd_2stage.f"
	w[1] = a[i__1].r;
#line 375 "cheevd_2stage.f"
	if (wantz) {
#line 375 "cheevd_2stage.f"
	    i__1 = a_dim1 + 1;
#line 375 "cheevd_2stage.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 375 "cheevd_2stage.f"
	}
#line 377 "cheevd_2stage.f"
	return 0;
#line 378 "cheevd_2stage.f"
    }

/*     Get machine constants. */

#line 382 "cheevd_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 383 "cheevd_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 384 "cheevd_2stage.f"
    smlnum = safmin / eps;
#line 385 "cheevd_2stage.f"
    bignum = 1. / smlnum;
#line 386 "cheevd_2stage.f"
    rmin = sqrt(smlnum);
#line 387 "cheevd_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 391 "cheevd_2stage.f"
    anrm = clanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 392 "cheevd_2stage.f"
    iscale = 0;
#line 393 "cheevd_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 394 "cheevd_2stage.f"
	iscale = 1;
#line 395 "cheevd_2stage.f"
	sigma = rmin / anrm;
#line 396 "cheevd_2stage.f"
    } else if (anrm > rmax) {
#line 397 "cheevd_2stage.f"
	iscale = 1;
#line 398 "cheevd_2stage.f"
	sigma = rmax / anrm;
#line 399 "cheevd_2stage.f"
    }
#line 400 "cheevd_2stage.f"
    if (iscale == 1) {
#line 400 "cheevd_2stage.f"
	clascl_(uplo, &c__0, &c__0, &c_b28, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 400 "cheevd_2stage.f"
    }

/*     Call CHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form. */

#line 405 "cheevd_2stage.f"
    inde = 1;
#line 406 "cheevd_2stage.f"
    indrwk = inde + *n;
#line 407 "cheevd_2stage.f"
    llrwk = *lrwork - indrwk + 1;
#line 408 "cheevd_2stage.f"
    indtau = 1;
#line 409 "cheevd_2stage.f"
    indhous = indtau + *n;
#line 410 "cheevd_2stage.f"
    indwrk = indhous + lhtrd;
#line 411 "cheevd_2stage.f"
    llwork = *lwork - indwrk + 1;
#line 412 "cheevd_2stage.f"
    indwk2 = indwrk + *n * *n;
#line 413 "cheevd_2stage.f"
    llwrk2 = *lwork - indwk2 + 1;

#line 415 "cheevd_2stage.f"
    chetrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &
	    work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     CSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
/*     tridiagonal matrix, then call CUNMTR to multiply it to the */
/*     Householder transformations represented as Householder vectors in */
/*     A. */

#line 425 "cheevd_2stage.f"
    if (! wantz) {
#line 426 "cheevd_2stage.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 427 "cheevd_2stage.f"
    } else {
#line 428 "cheevd_2stage.f"
	cstedc_("I", n, &w[1], &rwork[inde], &work[indwrk], n, &work[indwk2], 
		&llwrk2, &rwork[indrwk], &llrwk, &iwork[1], liwork, info, (
		ftnlen)1);
#line 431 "cheevd_2stage.f"
	cunmtr_("L", uplo, "N", n, n, &a[a_offset], lda, &work[indtau], &work[
		indwrk], n, &work[indwk2], &llwrk2, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 433 "cheevd_2stage.f"
	clacpy_("A", n, n, &work[indwrk], n, &a[a_offset], lda, (ftnlen)1);
#line 434 "cheevd_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 438 "cheevd_2stage.f"
    if (iscale == 1) {
#line 439 "cheevd_2stage.f"
	if (*info == 0) {
#line 440 "cheevd_2stage.f"
	    imax = *n;
#line 441 "cheevd_2stage.f"
	} else {
#line 442 "cheevd_2stage.f"
	    imax = *info - 1;
#line 443 "cheevd_2stage.f"
	}
#line 444 "cheevd_2stage.f"
	d__1 = 1. / sigma;
#line 444 "cheevd_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 445 "cheevd_2stage.f"
    }

#line 447 "cheevd_2stage.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 448 "cheevd_2stage.f"
    rwork[1] = (doublereal) lrwmin;
#line 449 "cheevd_2stage.f"
    iwork[1] = liwmin;

#line 451 "cheevd_2stage.f"
    return 0;

/*     End of CHEEVD_2STAGE */

} /* cheevd_2stage__ */

