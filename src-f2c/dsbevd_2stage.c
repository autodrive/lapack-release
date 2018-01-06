#line 1 "dsbevd_2stage.f"
/* dsbevd_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsbevd_2stage.f"
/* Table of constant values */

static integer c__18 = 18;
static integer c_n1 = -1;
static integer c__19 = 19;
static integer c__20 = 20;
static doublereal c_b21 = 1.;
static doublereal c_b29 = 0.;
static integer c__1 = 1;

/* > \brief <b> DSBEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 OTHER matrices</b> */

/*  @precisions fortran d -> s */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBEVD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbevd_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbevd_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbevd_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/*                                 WORK, LWORK, IWORK, LIWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBEVD_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
/* > a real symmetric band matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. If eigenvectors are desired, it uses */
/* > a divide and conquer algorithm. */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension LWORK */
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
/* >                                   dimension = (2KD+1)*N + KD*NTHREADS + N */
/* >                                   where KD is the size of the band. */
/* >                                   NTHREADS is the number of threads used when */
/* >                                   openMP compilation is enabled, otherwise =1. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available. */
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
/* >          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N. */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHEReigen */

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
/* Subroutine */ int dsbevd_2stage__(char *jobz, char *uplo, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *w, doublereal 
	*z__, integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
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
    static doublereal anrm, rmin, rmax;
    extern /* Subroutine */ int dsytrd_sb2st__(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dscal_(integer *, doublereal *
	    , doublereal *, integer *), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo, lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    static integer indwk2, llwrk2;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern doublereal dlansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, ftnlen), dlacpy_(char *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer indwrk, liwmin, llwork;
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

#line 284 "dsbevd_2stage.f"
    /* Parameter adjustments */
#line 284 "dsbevd_2stage.f"
    ab_dim1 = *ldab;
#line 284 "dsbevd_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 284 "dsbevd_2stage.f"
    ab -= ab_offset;
#line 284 "dsbevd_2stage.f"
    --w;
#line 284 "dsbevd_2stage.f"
    z_dim1 = *ldz;
#line 284 "dsbevd_2stage.f"
    z_offset = 1 + z_dim1;
#line 284 "dsbevd_2stage.f"
    z__ -= z_offset;
#line 284 "dsbevd_2stage.f"
    --work;
#line 284 "dsbevd_2stage.f"
    --iwork;
#line 284 "dsbevd_2stage.f"

#line 284 "dsbevd_2stage.f"
    /* Function Body */
#line 284 "dsbevd_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 285 "dsbevd_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 286 "dsbevd_2stage.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 288 "dsbevd_2stage.f"
    *info = 0;
#line 289 "dsbevd_2stage.f"
    if (*n <= 1) {
#line 290 "dsbevd_2stage.f"
	liwmin = 1;
#line 291 "dsbevd_2stage.f"
	lwmin = 1;
#line 292 "dsbevd_2stage.f"
    } else {
#line 293 "dsbevd_2stage.f"
	ib = ilaenv_(&c__18, "DSYTRD_SB2ST", jobz, n, kd, &c_n1, &c_n1, (
		ftnlen)12, (ftnlen)1);
#line 294 "dsbevd_2stage.f"
	lhtrd = ilaenv_(&c__19, "DSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1, (
		ftnlen)12, (ftnlen)1);
#line 295 "dsbevd_2stage.f"
	lwtrd = ilaenv_(&c__20, "DSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1, (
		ftnlen)12, (ftnlen)1);
#line 296 "dsbevd_2stage.f"
	if (wantz) {
#line 297 "dsbevd_2stage.f"
	    liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 298 "dsbevd_2stage.f"
	    i__1 = *n;
#line 298 "dsbevd_2stage.f"
	    lwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 299 "dsbevd_2stage.f"
	} else {
#line 300 "dsbevd_2stage.f"
	    liwmin = 1;
/* Computing MAX */
#line 301 "dsbevd_2stage.f"
	    i__1 = *n << 1, i__2 = *n + lhtrd + lwtrd;
#line 301 "dsbevd_2stage.f"
	    lwmin = max(i__1,i__2);
#line 302 "dsbevd_2stage.f"
	}
#line 303 "dsbevd_2stage.f"
    }
#line 304 "dsbevd_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 305 "dsbevd_2stage.f"
	*info = -1;
#line 306 "dsbevd_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 307 "dsbevd_2stage.f"
	*info = -2;
#line 308 "dsbevd_2stage.f"
    } else if (*n < 0) {
#line 309 "dsbevd_2stage.f"
	*info = -3;
#line 310 "dsbevd_2stage.f"
    } else if (*kd < 0) {
#line 311 "dsbevd_2stage.f"
	*info = -4;
#line 312 "dsbevd_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 313 "dsbevd_2stage.f"
	*info = -6;
#line 314 "dsbevd_2stage.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 315 "dsbevd_2stage.f"
	*info = -9;
#line 316 "dsbevd_2stage.f"
    }

#line 318 "dsbevd_2stage.f"
    if (*info == 0) {
#line 319 "dsbevd_2stage.f"
	work[1] = (doublereal) lwmin;
#line 320 "dsbevd_2stage.f"
	iwork[1] = liwmin;

#line 322 "dsbevd_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 323 "dsbevd_2stage.f"
	    *info = -11;
#line 324 "dsbevd_2stage.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 325 "dsbevd_2stage.f"
	    *info = -13;
#line 326 "dsbevd_2stage.f"
	}
#line 327 "dsbevd_2stage.f"
    }

#line 329 "dsbevd_2stage.f"
    if (*info != 0) {
#line 330 "dsbevd_2stage.f"
	i__1 = -(*info);
#line 330 "dsbevd_2stage.f"
	xerbla_("DSBEVD_2STAGE", &i__1, (ftnlen)13);
#line 331 "dsbevd_2stage.f"
	return 0;
#line 332 "dsbevd_2stage.f"
    } else if (lquery) {
#line 333 "dsbevd_2stage.f"
	return 0;
#line 334 "dsbevd_2stage.f"
    }

/*     Quick return if possible */

#line 338 "dsbevd_2stage.f"
    if (*n == 0) {
#line 338 "dsbevd_2stage.f"
	return 0;
#line 338 "dsbevd_2stage.f"
    }

#line 341 "dsbevd_2stage.f"
    if (*n == 1) {
#line 342 "dsbevd_2stage.f"
	w[1] = ab[ab_dim1 + 1];
#line 343 "dsbevd_2stage.f"
	if (wantz) {
#line 343 "dsbevd_2stage.f"
	    z__[z_dim1 + 1] = 1.;
#line 343 "dsbevd_2stage.f"
	}
#line 345 "dsbevd_2stage.f"
	return 0;
#line 346 "dsbevd_2stage.f"
    }

/*     Get machine constants. */

#line 350 "dsbevd_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 351 "dsbevd_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 352 "dsbevd_2stage.f"
    smlnum = safmin / eps;
#line 353 "dsbevd_2stage.f"
    bignum = 1. / smlnum;
#line 354 "dsbevd_2stage.f"
    rmin = sqrt(smlnum);
#line 355 "dsbevd_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 359 "dsbevd_2stage.f"
    anrm = dlansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 360 "dsbevd_2stage.f"
    iscale = 0;
#line 361 "dsbevd_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 362 "dsbevd_2stage.f"
	iscale = 1;
#line 363 "dsbevd_2stage.f"
	sigma = rmin / anrm;
#line 364 "dsbevd_2stage.f"
    } else if (anrm > rmax) {
#line 365 "dsbevd_2stage.f"
	iscale = 1;
#line 366 "dsbevd_2stage.f"
	sigma = rmax / anrm;
#line 367 "dsbevd_2stage.f"
    }
#line 368 "dsbevd_2stage.f"
    if (iscale == 1) {
#line 369 "dsbevd_2stage.f"
	if (lower) {
#line 370 "dsbevd_2stage.f"
	    dlascl_("B", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 371 "dsbevd_2stage.f"
	} else {
#line 372 "dsbevd_2stage.f"
	    dlascl_("Q", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 373 "dsbevd_2stage.f"
	}
#line 374 "dsbevd_2stage.f"
    }

/*     Call DSYTRD_SB2ST to reduce band symmetric matrix to tridiagonal form. */

#line 378 "dsbevd_2stage.f"
    inde = 1;
#line 379 "dsbevd_2stage.f"
    indhous = inde + *n;
#line 380 "dsbevd_2stage.f"
    indwrk = indhous + lhtrd;
#line 381 "dsbevd_2stage.f"
    llwork = *lwork - indwrk + 1;
#line 382 "dsbevd_2stage.f"
    indwk2 = indwrk + *n * *n;
#line 383 "dsbevd_2stage.f"
    llwrk2 = *lwork - indwk2 + 1;

#line 385 "dsbevd_2stage.f"
    dsytrd_sb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[
	    inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &iinfo, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC. */

#line 391 "dsbevd_2stage.f"
    if (! wantz) {
#line 392 "dsbevd_2stage.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 393 "dsbevd_2stage.f"
    } else {
#line 394 "dsbevd_2stage.f"
	dstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &
		llwrk2, &iwork[1], liwork, info, (ftnlen)1);
#line 396 "dsbevd_2stage.f"
	dgemm_("N", "N", n, n, n, &c_b21, &z__[z_offset], ldz, &work[indwrk], 
		n, &c_b29, &work[indwk2], n, (ftnlen)1, (ftnlen)1);
#line 398 "dsbevd_2stage.f"
	dlacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz, (ftnlen)1);
#line 399 "dsbevd_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 403 "dsbevd_2stage.f"
    if (iscale == 1) {
#line 403 "dsbevd_2stage.f"
	d__1 = 1. / sigma;
#line 403 "dsbevd_2stage.f"
	dscal_(n, &d__1, &w[1], &c__1);
#line 403 "dsbevd_2stage.f"
    }

#line 406 "dsbevd_2stage.f"
    work[1] = (doublereal) lwmin;
#line 407 "dsbevd_2stage.f"
    iwork[1] = liwmin;
#line 408 "dsbevd_2stage.f"
    return 0;

/*     End of DSBEVD_2STAGE */

} /* dsbevd_2stage__ */

