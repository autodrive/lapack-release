#line 1 "chbev_2stage.f"
/* chbev_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "chbev_2stage.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__4 = 4;
static doublereal c_b21 = 1.;
static integer c__1 = 1;

/* > \brief <b> CHBEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for 
OTHER matrices</b> */

/*  @generated from zhbev_2stage.f, fortran z -> c, Sat Nov  5 23:18:20 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBEV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbev_2
stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbev_2
stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbev_2
stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/*                                WORK, LWORK, RWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, N, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEV_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian band matrix A using the 2stage technique for */
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
/* >          WORK is COMPLEX array, dimension LWORK */
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
/* >          RWORK is REAL array, dimension (max(1,3*N-2)) */
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

/* > \ingroup complexOTHEReigen */

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
/* Subroutine */ int chbev_2stage__(char *jobz, char *uplo, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublereal *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork,
	 doublereal *rwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ib;
    static doublereal eps;
    extern /* Subroutine */ int chetrd_hb2st__(char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    extern doublereal clanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indwrk, indrwk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), ssterf_(integer *, doublereal *, doublereal *, integer *
	    );
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

#line 260 "chbev_2stage.f"
    /* Parameter adjustments */
#line 260 "chbev_2stage.f"
    ab_dim1 = *ldab;
#line 260 "chbev_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 260 "chbev_2stage.f"
    ab -= ab_offset;
#line 260 "chbev_2stage.f"
    --w;
#line 260 "chbev_2stage.f"
    z_dim1 = *ldz;
#line 260 "chbev_2stage.f"
    z_offset = 1 + z_dim1;
#line 260 "chbev_2stage.f"
    z__ -= z_offset;
#line 260 "chbev_2stage.f"
    --work;
#line 260 "chbev_2stage.f"
    --rwork;
#line 260 "chbev_2stage.f"

#line 260 "chbev_2stage.f"
    /* Function Body */
#line 260 "chbev_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 261 "chbev_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 262 "chbev_2stage.f"
    lquery = *lwork == -1;

#line 264 "chbev_2stage.f"
    *info = 0;
#line 265 "chbev_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 266 "chbev_2stage.f"
	*info = -1;
#line 267 "chbev_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 268 "chbev_2stage.f"
	*info = -2;
#line 269 "chbev_2stage.f"
    } else if (*n < 0) {
#line 270 "chbev_2stage.f"
	*info = -3;
#line 271 "chbev_2stage.f"
    } else if (*kd < 0) {
#line 272 "chbev_2stage.f"
	*info = -4;
#line 273 "chbev_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 274 "chbev_2stage.f"
	*info = -6;
#line 275 "chbev_2stage.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 276 "chbev_2stage.f"
	*info = -9;
#line 277 "chbev_2stage.f"
    }

#line 279 "chbev_2stage.f"
    if (*info == 0) {
#line 280 "chbev_2stage.f"
	if (*n <= 1) {
#line 281 "chbev_2stage.f"
	    lwmin = 1;
#line 282 "chbev_2stage.f"
	    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 283 "chbev_2stage.f"
	} else {
#line 284 "chbev_2stage.f"
	    ib = ilaenv2stage_(&c__2, "CHETRD_HB2ST", jobz, n, kd, &c_n1, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 286 "chbev_2stage.f"
	    lhtrd = ilaenv2stage_(&c__3, "CHETRD_HB2ST", jobz, n, kd, &ib, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 288 "chbev_2stage.f"
	    lwtrd = ilaenv2stage_(&c__4, "CHETRD_HB2ST", jobz, n, kd, &ib, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 290 "chbev_2stage.f"
	    lwmin = lhtrd + lwtrd;
#line 291 "chbev_2stage.f"
	    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 292 "chbev_2stage.f"
	}

#line 294 "chbev_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 294 "chbev_2stage.f"
	    *info = -11;
#line 294 "chbev_2stage.f"
	}
#line 296 "chbev_2stage.f"
    }

#line 298 "chbev_2stage.f"
    if (*info != 0) {
#line 299 "chbev_2stage.f"
	i__1 = -(*info);
#line 299 "chbev_2stage.f"
	xerbla_("CHBEV_2STAGE ", &i__1, (ftnlen)13);
#line 300 "chbev_2stage.f"
	return 0;
#line 301 "chbev_2stage.f"
    } else if (lquery) {
#line 302 "chbev_2stage.f"
	return 0;
#line 303 "chbev_2stage.f"
    }

/*     Quick return if possible */

#line 307 "chbev_2stage.f"
    if (*n == 0) {
#line 307 "chbev_2stage.f"
	return 0;
#line 307 "chbev_2stage.f"
    }

#line 310 "chbev_2stage.f"
    if (*n == 1) {
#line 311 "chbev_2stage.f"
	if (lower) {
#line 312 "chbev_2stage.f"
	    i__1 = ab_dim1 + 1;
#line 312 "chbev_2stage.f"
	    w[1] = ab[i__1].r;
#line 313 "chbev_2stage.f"
	} else {
#line 314 "chbev_2stage.f"
	    i__1 = *kd + 1 + ab_dim1;
#line 314 "chbev_2stage.f"
	    w[1] = ab[i__1].r;
#line 315 "chbev_2stage.f"
	}
#line 316 "chbev_2stage.f"
	if (wantz) {
#line 316 "chbev_2stage.f"
	    i__1 = z_dim1 + 1;
#line 316 "chbev_2stage.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 316 "chbev_2stage.f"
	}
#line 318 "chbev_2stage.f"
	return 0;
#line 319 "chbev_2stage.f"
    }

/*     Get machine constants. */

#line 323 "chbev_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 324 "chbev_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 325 "chbev_2stage.f"
    smlnum = safmin / eps;
#line 326 "chbev_2stage.f"
    bignum = 1. / smlnum;
#line 327 "chbev_2stage.f"
    rmin = sqrt(smlnum);
#line 328 "chbev_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 332 "chbev_2stage.f"
    anrm = clanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1, (ftnlen)1);
#line 333 "chbev_2stage.f"
    iscale = 0;
#line 334 "chbev_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 335 "chbev_2stage.f"
	iscale = 1;
#line 336 "chbev_2stage.f"
	sigma = rmin / anrm;
#line 337 "chbev_2stage.f"
    } else if (anrm > rmax) {
#line 338 "chbev_2stage.f"
	iscale = 1;
#line 339 "chbev_2stage.f"
	sigma = rmax / anrm;
#line 340 "chbev_2stage.f"
    }
#line 341 "chbev_2stage.f"
    if (iscale == 1) {
#line 342 "chbev_2stage.f"
	if (lower) {
#line 343 "chbev_2stage.f"
	    clascl_("B", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 344 "chbev_2stage.f"
	} else {
#line 345 "chbev_2stage.f"
	    clascl_("Q", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 346 "chbev_2stage.f"
	}
#line 347 "chbev_2stage.f"
    }

/*     Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form. */

#line 351 "chbev_2stage.f"
    inde = 1;
#line 352 "chbev_2stage.f"
    indhous = 1;
#line 353 "chbev_2stage.f"
    indwrk = indhous + lhtrd;
#line 354 "chbev_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 356 "chbev_2stage.f"
    chetrd_hb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &
	    rwork[inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR. */

#line 362 "chbev_2stage.f"
    if (! wantz) {
#line 363 "chbev_2stage.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 364 "chbev_2stage.f"
    } else {
#line 365 "chbev_2stage.f"
	indrwk = inde + *n;
#line 366 "chbev_2stage.f"
	csteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[
		indrwk], info, (ftnlen)1);
#line 368 "chbev_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 372 "chbev_2stage.f"
    if (iscale == 1) {
#line 373 "chbev_2stage.f"
	if (*info == 0) {
#line 374 "chbev_2stage.f"
	    imax = *n;
#line 375 "chbev_2stage.f"
	} else {
#line 376 "chbev_2stage.f"
	    imax = *info - 1;
#line 377 "chbev_2stage.f"
	}
#line 378 "chbev_2stage.f"
	d__1 = 1. / sigma;
#line 378 "chbev_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 379 "chbev_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 383 "chbev_2stage.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 385 "chbev_2stage.f"
    return 0;

/*     End of CHBEV_2STAGE */

} /* chbev_2stage__ */

