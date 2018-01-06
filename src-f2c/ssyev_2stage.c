#line 1 "ssyev_2stage.f"
/* ssyev_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "ssyev_2stage.f"
/* Table of constant values */

static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__0 = 0;
static doublereal c_b27 = 1.;
static integer c__1 = 1;

/* > \brief <b> SSYEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for 
SY matrices</b> */

/*  @generated from dsyev_2stage.f, fortran d -> s, Sat Nov  5 23:55:51 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYEV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevd_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevd_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevd_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, */
/*                                INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension LWORK */
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
/* >                                   dimension = max(stage1,stage2) + (KD+1)*N + 2*N */
/* >                                             = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                                               + max(2*KD*KD, KD*NTHREADS) */
/* >                                               + (KD+1)*N + 2*N */
/* >                                   where KD is the blocking size of the reduction, */
/* >                                   FACTOPTNB is the blocking used by the QR or LQ */
/* >                                   algorithm, usually FACTOPTNB=128 is a good choice */
/* >                                   NTHREADS is the number of threads used when */
/* >                                   openMP compilation is enabled, otherwise =1. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \ingroup realSYeigen */

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
/* Subroutine */ int ssyev_2stage__(char *jobz, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *w, doublereal *work, integer 
	*lwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
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
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lhtrd, lwmin;
    static logical lower;
    extern /* Subroutine */ int ssytrd_2stage__(char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer lwtrd;
    static logical wantz;
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer indtau, indwrk;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal slansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer llwork;
    static doublereal smlnum;
    extern /* Subroutine */ int sorgtr_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static logical lquery;
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
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

#line 231 "ssyev_2stage.f"
    /* Parameter adjustments */
#line 231 "ssyev_2stage.f"
    a_dim1 = *lda;
#line 231 "ssyev_2stage.f"
    a_offset = 1 + a_dim1;
#line 231 "ssyev_2stage.f"
    a -= a_offset;
#line 231 "ssyev_2stage.f"
    --w;
#line 231 "ssyev_2stage.f"
    --work;
#line 231 "ssyev_2stage.f"

#line 231 "ssyev_2stage.f"
    /* Function Body */
#line 231 "ssyev_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 232 "ssyev_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 233 "ssyev_2stage.f"
    lquery = *lwork == -1;

#line 235 "ssyev_2stage.f"
    *info = 0;
#line 236 "ssyev_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 237 "ssyev_2stage.f"
	*info = -1;
#line 238 "ssyev_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 239 "ssyev_2stage.f"
	*info = -2;
#line 240 "ssyev_2stage.f"
    } else if (*n < 0) {
#line 241 "ssyev_2stage.f"
	*info = -3;
#line 242 "ssyev_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 243 "ssyev_2stage.f"
	*info = -5;
#line 244 "ssyev_2stage.f"
    }

#line 246 "ssyev_2stage.f"
    if (*info == 0) {
#line 247 "ssyev_2stage.f"
	kd = ilaenv_(&c__17, "SSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 248 "ssyev_2stage.f"
	ib = ilaenv_(&c__18, "SSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 249 "ssyev_2stage.f"
	lhtrd = ilaenv_(&c__19, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 250 "ssyev_2stage.f"
	lwtrd = ilaenv_(&c__20, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 251 "ssyev_2stage.f"
	lwmin = (*n << 1) + lhtrd + lwtrd;
#line 252 "ssyev_2stage.f"
	work[1] = (doublereal) lwmin;

#line 254 "ssyev_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 254 "ssyev_2stage.f"
	    *info = -8;
#line 254 "ssyev_2stage.f"
	}
#line 256 "ssyev_2stage.f"
    }

#line 258 "ssyev_2stage.f"
    if (*info != 0) {
#line 259 "ssyev_2stage.f"
	i__1 = -(*info);
#line 259 "ssyev_2stage.f"
	xerbla_("SSYEV_2STAGE ", &i__1, (ftnlen)13);
#line 260 "ssyev_2stage.f"
	return 0;
#line 261 "ssyev_2stage.f"
    } else if (lquery) {
#line 262 "ssyev_2stage.f"
	return 0;
#line 263 "ssyev_2stage.f"
    }

/*     Quick return if possible */

#line 267 "ssyev_2stage.f"
    if (*n == 0) {
#line 268 "ssyev_2stage.f"
	return 0;
#line 269 "ssyev_2stage.f"
    }

#line 271 "ssyev_2stage.f"
    if (*n == 1) {
#line 272 "ssyev_2stage.f"
	w[1] = a[a_dim1 + 1];
#line 273 "ssyev_2stage.f"
	work[1] = 2.;
#line 274 "ssyev_2stage.f"
	if (wantz) {
#line 274 "ssyev_2stage.f"
	    a[a_dim1 + 1] = 1.;
#line 274 "ssyev_2stage.f"
	}
#line 276 "ssyev_2stage.f"
	return 0;
#line 277 "ssyev_2stage.f"
    }

/*     Get machine constants. */

#line 281 "ssyev_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 282 "ssyev_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 283 "ssyev_2stage.f"
    smlnum = safmin / eps;
#line 284 "ssyev_2stage.f"
    bignum = 1. / smlnum;
#line 285 "ssyev_2stage.f"
    rmin = sqrt(smlnum);
#line 286 "ssyev_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 290 "ssyev_2stage.f"
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 291 "ssyev_2stage.f"
    iscale = 0;
#line 292 "ssyev_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 293 "ssyev_2stage.f"
	iscale = 1;
#line 294 "ssyev_2stage.f"
	sigma = rmin / anrm;
#line 295 "ssyev_2stage.f"
    } else if (anrm > rmax) {
#line 296 "ssyev_2stage.f"
	iscale = 1;
#line 297 "ssyev_2stage.f"
	sigma = rmax / anrm;
#line 298 "ssyev_2stage.f"
    }
#line 299 "ssyev_2stage.f"
    if (iscale == 1) {
#line 299 "ssyev_2stage.f"
	slascl_(uplo, &c__0, &c__0, &c_b27, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 299 "ssyev_2stage.f"
    }

/*     Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */

#line 304 "ssyev_2stage.f"
    inde = 1;
#line 305 "ssyev_2stage.f"
    indtau = inde + *n;
#line 306 "ssyev_2stage.f"
    indhous = indtau + *n;
#line 307 "ssyev_2stage.f"
    indwrk = indhous + lhtrd;
#line 308 "ssyev_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 310 "ssyev_2stage.f"
    ssytrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[inde], &
	    work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     SORGTR to generate the orthogonal matrix, then call SSTEQR. */

#line 317 "ssyev_2stage.f"
    if (! wantz) {
#line 318 "ssyev_2stage.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 319 "ssyev_2stage.f"
    } else {
/*        Not available in this release, and agrument checking should not */
/*        let it getting here */
#line 322 "ssyev_2stage.f"
	return 0;
#line 323 "ssyev_2stage.f"
	sorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo, (ftnlen)1);
#line 325 "ssyev_2stage.f"
	ssteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
		 info, (ftnlen)1);
#line 327 "ssyev_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 331 "ssyev_2stage.f"
    if (iscale == 1) {
#line 332 "ssyev_2stage.f"
	if (*info == 0) {
#line 333 "ssyev_2stage.f"
	    imax = *n;
#line 334 "ssyev_2stage.f"
	} else {
#line 335 "ssyev_2stage.f"
	    imax = *info - 1;
#line 336 "ssyev_2stage.f"
	}
#line 337 "ssyev_2stage.f"
	d__1 = 1. / sigma;
#line 337 "ssyev_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 338 "ssyev_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 342 "ssyev_2stage.f"
    work[1] = (doublereal) lwmin;

#line 344 "ssyev_2stage.f"
    return 0;

/*     End of SSYEV_2STAGE */

} /* ssyev_2stage__ */

