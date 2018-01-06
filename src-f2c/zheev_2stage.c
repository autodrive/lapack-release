#line 1 "zheev_2stage.f"
/* zheev_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "zheev_2stage.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__0 = 0;
static doublereal c_b28 = 1.;

/* > \brief <b> ZHEEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for 
HE matrices</b> */

/*  @precisions fortran z -> s d c */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheev_2
stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheev_2
stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheev_2
stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, */
/*                                RWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix A using the 2stage technique for */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
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
/* >                                   dimension = max(stage1,stage2) + (KD+1)*N + N */
/* >                                             = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                                               + max(2*KD*KD, KD*NTHREADS) */
/* >                                               + (KD+1)*N + N */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2)) */
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

/* > \date November 2017 */

/* > \ingroup complex16HEeigen */

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
/* Subroutine */ int zheev_2stage__(char *jobz, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info, ftnlen jobz_len, 
	ftnlen uplo_len)
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zhetrd_2stage__(char *, char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer iinfo, lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), zlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen);
    static integer indwrk, llwork;
    static doublereal smlnum;
    static logical lquery;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), zungtr_(char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen);
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

#line 240 "zheev_2stage.f"
    /* Parameter adjustments */
#line 240 "zheev_2stage.f"
    a_dim1 = *lda;
#line 240 "zheev_2stage.f"
    a_offset = 1 + a_dim1;
#line 240 "zheev_2stage.f"
    a -= a_offset;
#line 240 "zheev_2stage.f"
    --w;
#line 240 "zheev_2stage.f"
    --work;
#line 240 "zheev_2stage.f"
    --rwork;
#line 240 "zheev_2stage.f"

#line 240 "zheev_2stage.f"
    /* Function Body */
#line 240 "zheev_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 241 "zheev_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 242 "zheev_2stage.f"
    lquery = *lwork == -1;

#line 244 "zheev_2stage.f"
    *info = 0;
#line 245 "zheev_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 246 "zheev_2stage.f"
	*info = -1;
#line 247 "zheev_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 248 "zheev_2stage.f"
	*info = -2;
#line 249 "zheev_2stage.f"
    } else if (*n < 0) {
#line 250 "zheev_2stage.f"
	*info = -3;
#line 251 "zheev_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 252 "zheev_2stage.f"
	*info = -5;
#line 253 "zheev_2stage.f"
    }

#line 255 "zheev_2stage.f"
    if (*info == 0) {
#line 256 "zheev_2stage.f"
	kd = ilaenv2stage_(&c__1, "ZHETRD_2STAGE", jobz, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)13, (ftnlen)1);
#line 257 "zheev_2stage.f"
	ib = ilaenv2stage_(&c__2, "ZHETRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1,
		 (ftnlen)13, (ftnlen)1);
#line 258 "zheev_2stage.f"
	lhtrd = ilaenv2stage_(&c__3, "ZHETRD_2STAGE", jobz, n, &kd, &ib, &
		c_n1, (ftnlen)13, (ftnlen)1);
#line 259 "zheev_2stage.f"
	lwtrd = ilaenv2stage_(&c__4, "ZHETRD_2STAGE", jobz, n, &kd, &ib, &
		c_n1, (ftnlen)13, (ftnlen)1);
#line 260 "zheev_2stage.f"
	lwmin = *n + lhtrd + lwtrd;
#line 261 "zheev_2stage.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 263 "zheev_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 263 "zheev_2stage.f"
	    *info = -8;
#line 263 "zheev_2stage.f"
	}
#line 265 "zheev_2stage.f"
    }

#line 267 "zheev_2stage.f"
    if (*info != 0) {
#line 268 "zheev_2stage.f"
	i__1 = -(*info);
#line 268 "zheev_2stage.f"
	xerbla_("ZHEEV_2STAGE ", &i__1, (ftnlen)13);
#line 269 "zheev_2stage.f"
	return 0;
#line 270 "zheev_2stage.f"
    } else if (lquery) {
#line 271 "zheev_2stage.f"
	return 0;
#line 272 "zheev_2stage.f"
    }

/*     Quick return if possible */

#line 276 "zheev_2stage.f"
    if (*n == 0) {
#line 277 "zheev_2stage.f"
	return 0;
#line 278 "zheev_2stage.f"
    }

#line 280 "zheev_2stage.f"
    if (*n == 1) {
#line 281 "zheev_2stage.f"
	i__1 = a_dim1 + 1;
#line 281 "zheev_2stage.f"
	w[1] = a[i__1].r;
#line 282 "zheev_2stage.f"
	work[1].r = 1., work[1].i = 0.;
#line 283 "zheev_2stage.f"
	if (wantz) {
#line 283 "zheev_2stage.f"
	    i__1 = a_dim1 + 1;
#line 283 "zheev_2stage.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 283 "zheev_2stage.f"
	}
#line 285 "zheev_2stage.f"
	return 0;
#line 286 "zheev_2stage.f"
    }

/*     Get machine constants. */

#line 290 "zheev_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 291 "zheev_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 292 "zheev_2stage.f"
    smlnum = safmin / eps;
#line 293 "zheev_2stage.f"
    bignum = 1. / smlnum;
#line 294 "zheev_2stage.f"
    rmin = sqrt(smlnum);
#line 295 "zheev_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 299 "zheev_2stage.f"
    anrm = zlanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 300 "zheev_2stage.f"
    iscale = 0;
#line 301 "zheev_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 302 "zheev_2stage.f"
	iscale = 1;
#line 303 "zheev_2stage.f"
	sigma = rmin / anrm;
#line 304 "zheev_2stage.f"
    } else if (anrm > rmax) {
#line 305 "zheev_2stage.f"
	iscale = 1;
#line 306 "zheev_2stage.f"
	sigma = rmax / anrm;
#line 307 "zheev_2stage.f"
    }
#line 308 "zheev_2stage.f"
    if (iscale == 1) {
#line 308 "zheev_2stage.f"
	zlascl_(uplo, &c__0, &c__0, &c_b28, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 308 "zheev_2stage.f"
    }

/*     Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form. */

#line 313 "zheev_2stage.f"
    inde = 1;
#line 314 "zheev_2stage.f"
    indtau = 1;
#line 315 "zheev_2stage.f"
    indhous = indtau + *n;
#line 316 "zheev_2stage.f"
    indwrk = indhous + lhtrd;
#line 317 "zheev_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 319 "zheev_2stage.f"
    zhetrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &
	    work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     ZUNGTR to generate the unitary matrix, then call ZSTEQR. */

#line 326 "zheev_2stage.f"
    if (! wantz) {
#line 327 "zheev_2stage.f"
	dsterf_(n, &w[1], &rwork[inde], info);
#line 328 "zheev_2stage.f"
    } else {
#line 329 "zheev_2stage.f"
	zungtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo, (ftnlen)1);
#line 331 "zheev_2stage.f"
	indwrk = inde + *n;
#line 332 "zheev_2stage.f"
	zsteqr_(jobz, n, &w[1], &rwork[inde], &a[a_offset], lda, &rwork[
		indwrk], info, (ftnlen)1);
#line 334 "zheev_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 338 "zheev_2stage.f"
    if (iscale == 1) {
#line 339 "zheev_2stage.f"
	if (*info == 0) {
#line 340 "zheev_2stage.f"
	    imax = *n;
#line 341 "zheev_2stage.f"
	} else {
#line 342 "zheev_2stage.f"
	    imax = *info - 1;
#line 343 "zheev_2stage.f"
	}
#line 344 "zheev_2stage.f"
	d__1 = 1. / sigma;
#line 344 "zheev_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 345 "zheev_2stage.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 349 "zheev_2stage.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 351 "zheev_2stage.f"
    return 0;

/*     End of ZHEEV_2STAGE */

} /* zheev_2stage__ */

