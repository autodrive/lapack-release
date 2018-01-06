#line 1 "dsbev_2stage.f"
/* dsbev_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsbev_2stage.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__4 = 4;
static doublereal c_b21 = 1.;
static integer c__1 = 1;

/* > \brief <b> DSBEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for 
OTHER matrices</b> */

/*  @precisions fortran d -> s */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBEV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbev_2
stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbev_2
stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbev_2
stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/*                                WORK, LWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, N, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBEV_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
/* > a real symmetric band matrix A using the 2stage technique for */
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

/* > \date November 2017 */

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
/* Subroutine */ int dsbev_2stage__(char *jobz, char *uplo, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *w, doublereal 
	*z__, integer *ldz, doublereal *work, integer *lwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
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
    static doublereal rmin, rmax;
    extern /* Subroutine */ int dsytrd_sb2st__(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dscal_(integer *, doublereal *
	    , doublereal *, integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo, lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern doublereal dlansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer indwrk;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
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

#line 252 "dsbev_2stage.f"
    /* Parameter adjustments */
#line 252 "dsbev_2stage.f"
    ab_dim1 = *ldab;
#line 252 "dsbev_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 252 "dsbev_2stage.f"
    ab -= ab_offset;
#line 252 "dsbev_2stage.f"
    --w;
#line 252 "dsbev_2stage.f"
    z_dim1 = *ldz;
#line 252 "dsbev_2stage.f"
    z_offset = 1 + z_dim1;
#line 252 "dsbev_2stage.f"
    z__ -= z_offset;
#line 252 "dsbev_2stage.f"
    --work;
#line 252 "dsbev_2stage.f"

#line 252 "dsbev_2stage.f"
    /* Function Body */
#line 252 "dsbev_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 253 "dsbev_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 254 "dsbev_2stage.f"
    lquery = *lwork == -1;

#line 256 "dsbev_2stage.f"
    *info = 0;
#line 257 "dsbev_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 258 "dsbev_2stage.f"
	*info = -1;
#line 259 "dsbev_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 260 "dsbev_2stage.f"
	*info = -2;
#line 261 "dsbev_2stage.f"
    } else if (*n < 0) {
#line 262 "dsbev_2stage.f"
	*info = -3;
#line 263 "dsbev_2stage.f"
    } else if (*kd < 0) {
#line 264 "dsbev_2stage.f"
	*info = -4;
#line 265 "dsbev_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 266 "dsbev_2stage.f"
	*info = -6;
#line 267 "dsbev_2stage.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 268 "dsbev_2stage.f"
	*info = -9;
#line 269 "dsbev_2stage.f"
    }

#line 271 "dsbev_2stage.f"
    if (*info == 0) {
#line 272 "dsbev_2stage.f"
	if (*n <= 1) {
#line 273 "dsbev_2stage.f"
	    lwmin = 1;
#line 274 "dsbev_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 275 "dsbev_2stage.f"
	} else {
#line 276 "dsbev_2stage.f"
	    ib = ilaenv2stage_(&c__2, "DSYTRD_SB2ST", jobz, n, kd, &c_n1, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 278 "dsbev_2stage.f"
	    lhtrd = ilaenv2stage_(&c__3, "DSYTRD_SB2ST", jobz, n, kd, &ib, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 280 "dsbev_2stage.f"
	    lwtrd = ilaenv2stage_(&c__4, "DSYTRD_SB2ST", jobz, n, kd, &ib, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 282 "dsbev_2stage.f"
	    lwmin = *n + lhtrd + lwtrd;
#line 283 "dsbev_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 284 "dsbev_2stage.f"
	}

#line 286 "dsbev_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 286 "dsbev_2stage.f"
	    *info = -11;
#line 286 "dsbev_2stage.f"
	}
#line 288 "dsbev_2stage.f"
    }

#line 290 "dsbev_2stage.f"
    if (*info != 0) {
#line 291 "dsbev_2stage.f"
	i__1 = -(*info);
#line 291 "dsbev_2stage.f"
	xerbla_("DSBEV_2STAGE ", &i__1, (ftnlen)13);
#line 292 "dsbev_2stage.f"
	return 0;
#line 293 "dsbev_2stage.f"
    } else if (lquery) {
#line 294 "dsbev_2stage.f"
	return 0;
#line 295 "dsbev_2stage.f"
    }

/*     Quick return if possible */

#line 299 "dsbev_2stage.f"
    if (*n == 0) {
#line 299 "dsbev_2stage.f"
	return 0;
#line 299 "dsbev_2stage.f"
    }

#line 302 "dsbev_2stage.f"
    if (*n == 1) {
#line 303 "dsbev_2stage.f"
	if (lower) {
#line 304 "dsbev_2stage.f"
	    w[1] = ab[ab_dim1 + 1];
#line 305 "dsbev_2stage.f"
	} else {
#line 306 "dsbev_2stage.f"
	    w[1] = ab[*kd + 1 + ab_dim1];
#line 307 "dsbev_2stage.f"
	}
#line 308 "dsbev_2stage.f"
	if (wantz) {
#line 308 "dsbev_2stage.f"
	    z__[z_dim1 + 1] = 1.;
#line 308 "dsbev_2stage.f"
	}
#line 310 "dsbev_2stage.f"
	return 0;
#line 311 "dsbev_2stage.f"
    }

/*     Get machine constants. */

#line 315 "dsbev_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 316 "dsbev_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 317 "dsbev_2stage.f"
    smlnum = safmin / eps;
#line 318 "dsbev_2stage.f"
    bignum = 1. / smlnum;
#line 319 "dsbev_2stage.f"
    rmin = sqrt(smlnum);
#line 320 "dsbev_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 324 "dsbev_2stage.f"
    anrm = dlansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 325 "dsbev_2stage.f"
    iscale = 0;
#line 326 "dsbev_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 327 "dsbev_2stage.f"
	iscale = 1;
#line 328 "dsbev_2stage.f"
	sigma = rmin / anrm;
#line 329 "dsbev_2stage.f"
    } else if (anrm > rmax) {
#line 330 "dsbev_2stage.f"
	iscale = 1;
#line 331 "dsbev_2stage.f"
	sigma = rmax / anrm;
#line 332 "dsbev_2stage.f"
    }
#line 333 "dsbev_2stage.f"
    if (iscale == 1) {
#line 334 "dsbev_2stage.f"
	if (lower) {
#line 335 "dsbev_2stage.f"
	    dlascl_("B", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 336 "dsbev_2stage.f"
	} else {
#line 337 "dsbev_2stage.f"
	    dlascl_("Q", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 338 "dsbev_2stage.f"
	}
#line 339 "dsbev_2stage.f"
    }

/*     Call DSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form. */

#line 343 "dsbev_2stage.f"
    inde = 1;
#line 344 "dsbev_2stage.f"
    indhous = inde + *n;
#line 345 "dsbev_2stage.f"
    indwrk = indhous + lhtrd;
#line 346 "dsbev_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 348 "dsbev_2stage.f"
    dsytrd_sb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[
	    inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &iinfo, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR. */

#line 354 "dsbev_2stage.f"
    if (! wantz) {
#line 355 "dsbev_2stage.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 356 "dsbev_2stage.f"
    } else {
#line 357 "dsbev_2stage.f"
	dsteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[
		indwrk], info, (ftnlen)1);
#line 359 "dsbev_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 363 "dsbev_2stage.f"
    if (iscale == 1) {
#line 364 "dsbev_2stage.f"
	if (*info == 0) {
#line 365 "dsbev_2stage.f"
	    imax = *n;
#line 366 "dsbev_2stage.f"
	} else {
#line 367 "dsbev_2stage.f"
	    imax = *info - 1;
#line 368 "dsbev_2stage.f"
	}
#line 369 "dsbev_2stage.f"
	d__1 = 1. / sigma;
#line 369 "dsbev_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 370 "dsbev_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 374 "dsbev_2stage.f"
    work[1] = (doublereal) lwmin;

#line 376 "dsbev_2stage.f"
    return 0;

/*     End of DSBEV_2STAGE */

} /* dsbev_2stage__ */

