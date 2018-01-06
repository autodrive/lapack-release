#line 1 "dsyevd_2stage.f"
/* dsyevd_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsyevd_2stage.f"
/* Table of constant values */

static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__0 = 0;
static doublereal c_b27 = 1.;
static integer c__1 = 1;

/* > \brief <b> DSYEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 SY matrices</b> */

/*  @precisions fortran d -> s */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYEVD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyevd_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyevd_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyevd_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYEVD_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, */
/*                                IWORK, LIWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYEVD_2STAGE computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. If eigenvectors are desired, it uses a */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, */
/* >                                         dimension (LWORK) */
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
/* >                                   dimension = max(stage1,stage2) + (KD+1)*N + 2*N+1 */
/* >                                             = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                                               + max(2*KD*KD, KD*NTHREADS) */
/* >                                               + (KD+1)*N + 2*N+1 */
/* >                                   where KD is the blocking size of the reduction, */
/* >                                   FACTOPTNB is the blocking used by the QR or LQ */
/* >                                   algorithm, usually FACTOPTNB=128 is a good choice */
/* >                                   NTHREADS is the number of threads used when */
/* >                                   openMP compilation is enabled, otherwise =1. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be at least */
/* >                                                1 + 6*N + 2*N**2. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal sizes of the WORK and IWORK */
/* >          arrays, returns these values as the first entries of the WORK */
/* >          and IWORK arrays, and no error message related to LWORK or */
/* >          LIWORK is issued by XERBLA. */
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
/* >          routine only calculates the optimal sizes of the WORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK and IWORK arrays, and no error message related to */
/* >          LWORK or LIWORK is issued by XERBLA. */
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

/* > \date December 2016 */

/* > \ingroup doubleSYeigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* >  Modified by Francoise Tisseur, University of Tennessee \n */
/* >  Modified description of INFO. Sven, 16 Feb 05. \n */
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
/* Subroutine */ int dsyevd_2stage__(char *jobz, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *w, doublereal *work, integer 
	*lwork, integer *iwork, integer *liwork, integer *info, ftnlen 
	jobz_len, ftnlen uplo_len)
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
    static doublereal anrm, rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int dsytrd_2stage__(char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    static integer indwk2, llwrk2;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen), dlacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer indwrk, liwmin;
    extern /* Subroutine */ int dormtr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer llwork;
    static doublereal smlnum;
    static logical lquery;
    static integer indhous;



/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 278 "dsyevd_2stage.f"
    /* Parameter adjustments */
#line 278 "dsyevd_2stage.f"
    a_dim1 = *lda;
#line 278 "dsyevd_2stage.f"
    a_offset = 1 + a_dim1;
#line 278 "dsyevd_2stage.f"
    a -= a_offset;
#line 278 "dsyevd_2stage.f"
    --w;
#line 278 "dsyevd_2stage.f"
    --work;
#line 278 "dsyevd_2stage.f"
    --iwork;
#line 278 "dsyevd_2stage.f"

#line 278 "dsyevd_2stage.f"
    /* Function Body */
#line 278 "dsyevd_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 279 "dsyevd_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 280 "dsyevd_2stage.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 282 "dsyevd_2stage.f"
    *info = 0;
#line 283 "dsyevd_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 284 "dsyevd_2stage.f"
	*info = -1;
#line 285 "dsyevd_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 286 "dsyevd_2stage.f"
	*info = -2;
#line 287 "dsyevd_2stage.f"
    } else if (*n < 0) {
#line 288 "dsyevd_2stage.f"
	*info = -3;
#line 289 "dsyevd_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 290 "dsyevd_2stage.f"
	*info = -5;
#line 291 "dsyevd_2stage.f"
    }

#line 293 "dsyevd_2stage.f"
    if (*info == 0) {
#line 294 "dsyevd_2stage.f"
	if (*n <= 1) {
#line 295 "dsyevd_2stage.f"
	    liwmin = 1;
#line 296 "dsyevd_2stage.f"
	    lwmin = 1;
#line 297 "dsyevd_2stage.f"
	} else {
#line 298 "dsyevd_2stage.f"
	    kd = ilaenv_(&c__17, "DSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &
		    c_n1, (ftnlen)13, (ftnlen)1);
#line 299 "dsyevd_2stage.f"
	    ib = ilaenv_(&c__18, "DSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, 
		    (ftnlen)13, (ftnlen)1);
#line 300 "dsyevd_2stage.f"
	    lhtrd = ilaenv_(&c__19, "DSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1,
		     (ftnlen)13, (ftnlen)1);
#line 301 "dsyevd_2stage.f"
	    lwtrd = ilaenv_(&c__20, "DSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1,
		     (ftnlen)13, (ftnlen)1);
#line 302 "dsyevd_2stage.f"
	    if (wantz) {
#line 303 "dsyevd_2stage.f"
		liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 304 "dsyevd_2stage.f"
		i__1 = *n;
#line 304 "dsyevd_2stage.f"
		lwmin = *n * 6 + 1 + (i__1 * i__1 << 1);
#line 305 "dsyevd_2stage.f"
	    } else {
#line 306 "dsyevd_2stage.f"
		liwmin = 1;
#line 307 "dsyevd_2stage.f"
		lwmin = (*n << 1) + 1 + lhtrd + lwtrd;
#line 308 "dsyevd_2stage.f"
	    }
#line 309 "dsyevd_2stage.f"
	}
#line 310 "dsyevd_2stage.f"
	work[1] = (doublereal) lwmin;
#line 311 "dsyevd_2stage.f"
	iwork[1] = liwmin;

#line 313 "dsyevd_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 314 "dsyevd_2stage.f"
	    *info = -8;
#line 315 "dsyevd_2stage.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 316 "dsyevd_2stage.f"
	    *info = -10;
#line 317 "dsyevd_2stage.f"
	}
#line 318 "dsyevd_2stage.f"
    }

#line 320 "dsyevd_2stage.f"
    if (*info != 0) {
#line 321 "dsyevd_2stage.f"
	i__1 = -(*info);
#line 321 "dsyevd_2stage.f"
	xerbla_("DSYEVD_2STAGE", &i__1, (ftnlen)13);
#line 322 "dsyevd_2stage.f"
	return 0;
#line 323 "dsyevd_2stage.f"
    } else if (lquery) {
#line 324 "dsyevd_2stage.f"
	return 0;
#line 325 "dsyevd_2stage.f"
    }

/*     Quick return if possible */

#line 329 "dsyevd_2stage.f"
    if (*n == 0) {
#line 329 "dsyevd_2stage.f"
	return 0;
#line 329 "dsyevd_2stage.f"
    }

#line 332 "dsyevd_2stage.f"
    if (*n == 1) {
#line 333 "dsyevd_2stage.f"
	w[1] = a[a_dim1 + 1];
#line 334 "dsyevd_2stage.f"
	if (wantz) {
#line 334 "dsyevd_2stage.f"
	    a[a_dim1 + 1] = 1.;
#line 334 "dsyevd_2stage.f"
	}
#line 336 "dsyevd_2stage.f"
	return 0;
#line 337 "dsyevd_2stage.f"
    }

/*     Get machine constants. */

#line 341 "dsyevd_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 342 "dsyevd_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 343 "dsyevd_2stage.f"
    smlnum = safmin / eps;
#line 344 "dsyevd_2stage.f"
    bignum = 1. / smlnum;
#line 345 "dsyevd_2stage.f"
    rmin = sqrt(smlnum);
#line 346 "dsyevd_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 350 "dsyevd_2stage.f"
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 351 "dsyevd_2stage.f"
    iscale = 0;
#line 352 "dsyevd_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 353 "dsyevd_2stage.f"
	iscale = 1;
#line 354 "dsyevd_2stage.f"
	sigma = rmin / anrm;
#line 355 "dsyevd_2stage.f"
    } else if (anrm > rmax) {
#line 356 "dsyevd_2stage.f"
	iscale = 1;
#line 357 "dsyevd_2stage.f"
	sigma = rmax / anrm;
#line 358 "dsyevd_2stage.f"
    }
#line 359 "dsyevd_2stage.f"
    if (iscale == 1) {
#line 359 "dsyevd_2stage.f"
	dlascl_(uplo, &c__0, &c__0, &c_b27, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 359 "dsyevd_2stage.f"
    }

/*     Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */

#line 364 "dsyevd_2stage.f"
    inde = 1;
#line 365 "dsyevd_2stage.f"
    indtau = inde + *n;
#line 366 "dsyevd_2stage.f"
    indhous = indtau + *n;
#line 367 "dsyevd_2stage.f"
    indwrk = indhous + lhtrd;
#line 368 "dsyevd_2stage.f"
    llwork = *lwork - indwrk + 1;
#line 369 "dsyevd_2stage.f"
    indwk2 = indwrk + *n * *n;
#line 370 "dsyevd_2stage.f"
    llwrk2 = *lwork - indwk2 + 1;

#line 372 "dsyevd_2stage.f"
    dsytrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[inde], &
	    work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
/*     tridiagonal matrix, then call DORMTR to multiply it by the */
/*     Householder transformations stored in A. */

#line 381 "dsyevd_2stage.f"
    if (! wantz) {
#line 382 "dsyevd_2stage.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 383 "dsyevd_2stage.f"
    } else {
/*        Not available in this release, and agrument checking should not */
/*        let it getting here */
#line 386 "dsyevd_2stage.f"
	return 0;
#line 387 "dsyevd_2stage.f"
	dstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &
		llwrk2, &iwork[1], liwork, info, (ftnlen)1);
#line 389 "dsyevd_2stage.f"
	dormtr_("L", uplo, "N", n, n, &a[a_offset], lda, &work[indtau], &work[
		indwrk], n, &work[indwk2], &llwrk2, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 391 "dsyevd_2stage.f"
	dlacpy_("A", n, n, &work[indwrk], n, &a[a_offset], lda, (ftnlen)1);
#line 392 "dsyevd_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 396 "dsyevd_2stage.f"
    if (iscale == 1) {
#line 396 "dsyevd_2stage.f"
	d__1 = 1. / sigma;
#line 396 "dsyevd_2stage.f"
	dscal_(n, &d__1, &w[1], &c__1);
#line 396 "dsyevd_2stage.f"
    }

#line 399 "dsyevd_2stage.f"
    work[1] = (doublereal) lwmin;
#line 400 "dsyevd_2stage.f"
    iwork[1] = liwmin;

#line 402 "dsyevd_2stage.f"
    return 0;

/*     End of DSYEVD_2STAGE */

} /* dsyevd_2stage__ */

