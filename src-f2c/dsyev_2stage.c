#line 1 "dsyev_2stage.f"
/* dsyev_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsyev_2stage.f"
/* Table of constant values */

static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__0 = 0;
static doublereal c_b27 = 1.;
static integer c__1 = 1;

/* > \brief <b> DSYEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for 
SY matrices</b> */

/*  @precisions fortran d -> s */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYEV_2STAGE + dependencies */
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

/*       SUBROUTINE DSYEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, */
/*                                INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of a */
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

/* > \ingroup doubleSYeigen */

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
/* Subroutine */ int dsyev_2stage__(char *jobz, char *uplo, integer *n, 
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
    static integer inde, imax;
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
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
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
    static integer indwrk;
    extern /* Subroutine */ int dorgtr_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), dsteqr_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
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

#line 231 "dsyev_2stage.f"
    /* Parameter adjustments */
#line 231 "dsyev_2stage.f"
    a_dim1 = *lda;
#line 231 "dsyev_2stage.f"
    a_offset = 1 + a_dim1;
#line 231 "dsyev_2stage.f"
    a -= a_offset;
#line 231 "dsyev_2stage.f"
    --w;
#line 231 "dsyev_2stage.f"
    --work;
#line 231 "dsyev_2stage.f"

#line 231 "dsyev_2stage.f"
    /* Function Body */
#line 231 "dsyev_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 232 "dsyev_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 233 "dsyev_2stage.f"
    lquery = *lwork == -1;

#line 235 "dsyev_2stage.f"
    *info = 0;
#line 236 "dsyev_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 237 "dsyev_2stage.f"
	*info = -1;
#line 238 "dsyev_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 239 "dsyev_2stage.f"
	*info = -2;
#line 240 "dsyev_2stage.f"
    } else if (*n < 0) {
#line 241 "dsyev_2stage.f"
	*info = -3;
#line 242 "dsyev_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 243 "dsyev_2stage.f"
	*info = -5;
#line 244 "dsyev_2stage.f"
    }

#line 246 "dsyev_2stage.f"
    if (*info == 0) {
#line 247 "dsyev_2stage.f"
	kd = ilaenv_(&c__17, "DSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 248 "dsyev_2stage.f"
	ib = ilaenv_(&c__18, "DSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 249 "dsyev_2stage.f"
	lhtrd = ilaenv_(&c__19, "DSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 250 "dsyev_2stage.f"
	lwtrd = ilaenv_(&c__20, "DSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 251 "dsyev_2stage.f"
	lwmin = (*n << 1) + lhtrd + lwtrd;
#line 252 "dsyev_2stage.f"
	work[1] = (doublereal) lwmin;

#line 254 "dsyev_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 254 "dsyev_2stage.f"
	    *info = -8;
#line 254 "dsyev_2stage.f"
	}
#line 256 "dsyev_2stage.f"
    }

#line 258 "dsyev_2stage.f"
    if (*info != 0) {
#line 259 "dsyev_2stage.f"
	i__1 = -(*info);
#line 259 "dsyev_2stage.f"
	xerbla_("DSYEV_2STAGE ", &i__1, (ftnlen)13);
#line 260 "dsyev_2stage.f"
	return 0;
#line 261 "dsyev_2stage.f"
    } else if (lquery) {
#line 262 "dsyev_2stage.f"
	return 0;
#line 263 "dsyev_2stage.f"
    }

/*     Quick return if possible */

#line 267 "dsyev_2stage.f"
    if (*n == 0) {
#line 268 "dsyev_2stage.f"
	return 0;
#line 269 "dsyev_2stage.f"
    }

#line 271 "dsyev_2stage.f"
    if (*n == 1) {
#line 272 "dsyev_2stage.f"
	w[1] = a[a_dim1 + 1];
#line 273 "dsyev_2stage.f"
	work[1] = 2.;
#line 274 "dsyev_2stage.f"
	if (wantz) {
#line 274 "dsyev_2stage.f"
	    a[a_dim1 + 1] = 1.;
#line 274 "dsyev_2stage.f"
	}
#line 276 "dsyev_2stage.f"
	return 0;
#line 277 "dsyev_2stage.f"
    }

/*     Get machine constants. */

#line 281 "dsyev_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 282 "dsyev_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 283 "dsyev_2stage.f"
    smlnum = safmin / eps;
#line 284 "dsyev_2stage.f"
    bignum = 1. / smlnum;
#line 285 "dsyev_2stage.f"
    rmin = sqrt(smlnum);
#line 286 "dsyev_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 290 "dsyev_2stage.f"
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 291 "dsyev_2stage.f"
    iscale = 0;
#line 292 "dsyev_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 293 "dsyev_2stage.f"
	iscale = 1;
#line 294 "dsyev_2stage.f"
	sigma = rmin / anrm;
#line 295 "dsyev_2stage.f"
    } else if (anrm > rmax) {
#line 296 "dsyev_2stage.f"
	iscale = 1;
#line 297 "dsyev_2stage.f"
	sigma = rmax / anrm;
#line 298 "dsyev_2stage.f"
    }
#line 299 "dsyev_2stage.f"
    if (iscale == 1) {
#line 299 "dsyev_2stage.f"
	dlascl_(uplo, &c__0, &c__0, &c_b27, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 299 "dsyev_2stage.f"
    }

/*     Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */

#line 304 "dsyev_2stage.f"
    inde = 1;
#line 305 "dsyev_2stage.f"
    indtau = inde + *n;
#line 306 "dsyev_2stage.f"
    indhous = indtau + *n;
#line 307 "dsyev_2stage.f"
    indwrk = indhous + lhtrd;
#line 308 "dsyev_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 310 "dsyev_2stage.f"
    dsytrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[inde], &
	    work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     DORGTR to generate the orthogonal matrix, then call DSTEQR. */

#line 317 "dsyev_2stage.f"
    if (! wantz) {
#line 318 "dsyev_2stage.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 319 "dsyev_2stage.f"
    } else {
/*        Not available in this release, and agrument checking should not */
/*        let it getting here */
#line 322 "dsyev_2stage.f"
	return 0;
#line 323 "dsyev_2stage.f"
	dorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo, (ftnlen)1);
#line 325 "dsyev_2stage.f"
	dsteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
		 info, (ftnlen)1);
#line 327 "dsyev_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 331 "dsyev_2stage.f"
    if (iscale == 1) {
#line 332 "dsyev_2stage.f"
	if (*info == 0) {
#line 333 "dsyev_2stage.f"
	    imax = *n;
#line 334 "dsyev_2stage.f"
	} else {
#line 335 "dsyev_2stage.f"
	    imax = *info - 1;
#line 336 "dsyev_2stage.f"
	}
#line 337 "dsyev_2stage.f"
	d__1 = 1. / sigma;
#line 337 "dsyev_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 338 "dsyev_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 342 "dsyev_2stage.f"
    work[1] = (doublereal) lwmin;

#line 344 "dsyev_2stage.f"
    return 0;

/*     End of DSYEV_2STAGE */

} /* dsyev_2stage__ */

