#line 1 "ssbev_2stage.f"
/* ssbev_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "ssbev_2stage.f"
/* Table of constant values */

static integer c__18 = 18;
static integer c_n1 = -1;
static integer c__19 = 19;
static integer c__20 = 20;
static doublereal c_b21 = 1.;
static integer c__1 = 1;

/* > \brief <b> SSBEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for 
OTHER matrices</b> */

/*  @generated from dsbev_2stage.f, fortran d -> s, Sat Nov  5 23:58:09 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSBEV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbev_2
stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbev_2
stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbev_2
stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/*                                WORK, LWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, N, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBEV_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
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
/* >          AB is REAL array, dimension (LDAB, N) */
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
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, N) */
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

/* > \date December 2016 */

/* > \ingroup realOTHEReigen */

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
/* Subroutine */ int ssbev_2stage__(char *jobz, char *uplo, integer *n, 
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
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    extern /* Subroutine */ int ssytrd_sb2st__(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal slansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer llwork;
    static doublereal smlnum;
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

#line 252 "ssbev_2stage.f"
    /* Parameter adjustments */
#line 252 "ssbev_2stage.f"
    ab_dim1 = *ldab;
#line 252 "ssbev_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 252 "ssbev_2stage.f"
    ab -= ab_offset;
#line 252 "ssbev_2stage.f"
    --w;
#line 252 "ssbev_2stage.f"
    z_dim1 = *ldz;
#line 252 "ssbev_2stage.f"
    z_offset = 1 + z_dim1;
#line 252 "ssbev_2stage.f"
    z__ -= z_offset;
#line 252 "ssbev_2stage.f"
    --work;
#line 252 "ssbev_2stage.f"

#line 252 "ssbev_2stage.f"
    /* Function Body */
#line 252 "ssbev_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 253 "ssbev_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 254 "ssbev_2stage.f"
    lquery = *lwork == -1;

#line 256 "ssbev_2stage.f"
    *info = 0;
#line 257 "ssbev_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 258 "ssbev_2stage.f"
	*info = -1;
#line 259 "ssbev_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 260 "ssbev_2stage.f"
	*info = -2;
#line 261 "ssbev_2stage.f"
    } else if (*n < 0) {
#line 262 "ssbev_2stage.f"
	*info = -3;
#line 263 "ssbev_2stage.f"
    } else if (*kd < 0) {
#line 264 "ssbev_2stage.f"
	*info = -4;
#line 265 "ssbev_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 266 "ssbev_2stage.f"
	*info = -6;
#line 267 "ssbev_2stage.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 268 "ssbev_2stage.f"
	*info = -9;
#line 269 "ssbev_2stage.f"
    }

#line 271 "ssbev_2stage.f"
    if (*info == 0) {
#line 272 "ssbev_2stage.f"
	if (*n <= 1) {
#line 273 "ssbev_2stage.f"
	    lwmin = 1;
#line 274 "ssbev_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 275 "ssbev_2stage.f"
	} else {
#line 276 "ssbev_2stage.f"
	    ib = ilaenv_(&c__18, "SSYTRD_SB2ST", jobz, n, kd, &c_n1, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 277 "ssbev_2stage.f"
	    lhtrd = ilaenv_(&c__19, "SSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 278 "ssbev_2stage.f"
	    lwtrd = ilaenv_(&c__20, "SSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 279 "ssbev_2stage.f"
	    lwmin = *n + lhtrd + lwtrd;
#line 280 "ssbev_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 281 "ssbev_2stage.f"
	}

#line 283 "ssbev_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 283 "ssbev_2stage.f"
	    *info = -11;
#line 283 "ssbev_2stage.f"
	}
#line 285 "ssbev_2stage.f"
    }

#line 287 "ssbev_2stage.f"
    if (*info != 0) {
#line 288 "ssbev_2stage.f"
	i__1 = -(*info);
#line 288 "ssbev_2stage.f"
	xerbla_("SSBEV_2STAGE ", &i__1, (ftnlen)13);
#line 289 "ssbev_2stage.f"
	return 0;
#line 290 "ssbev_2stage.f"
    } else if (lquery) {
#line 291 "ssbev_2stage.f"
	return 0;
#line 292 "ssbev_2stage.f"
    }

/*     Quick return if possible */

#line 296 "ssbev_2stage.f"
    if (*n == 0) {
#line 296 "ssbev_2stage.f"
	return 0;
#line 296 "ssbev_2stage.f"
    }

#line 299 "ssbev_2stage.f"
    if (*n == 1) {
#line 300 "ssbev_2stage.f"
	if (lower) {
#line 301 "ssbev_2stage.f"
	    w[1] = ab[ab_dim1 + 1];
#line 302 "ssbev_2stage.f"
	} else {
#line 303 "ssbev_2stage.f"
	    w[1] = ab[*kd + 1 + ab_dim1];
#line 304 "ssbev_2stage.f"
	}
#line 305 "ssbev_2stage.f"
	if (wantz) {
#line 305 "ssbev_2stage.f"
	    z__[z_dim1 + 1] = 1.;
#line 305 "ssbev_2stage.f"
	}
#line 307 "ssbev_2stage.f"
	return 0;
#line 308 "ssbev_2stage.f"
    }

/*     Get machine constants. */

#line 312 "ssbev_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 313 "ssbev_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 314 "ssbev_2stage.f"
    smlnum = safmin / eps;
#line 315 "ssbev_2stage.f"
    bignum = 1. / smlnum;
#line 316 "ssbev_2stage.f"
    rmin = sqrt(smlnum);
#line 317 "ssbev_2stage.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 321 "ssbev_2stage.f"
    anrm = slansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 322 "ssbev_2stage.f"
    iscale = 0;
#line 323 "ssbev_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 324 "ssbev_2stage.f"
	iscale = 1;
#line 325 "ssbev_2stage.f"
	sigma = rmin / anrm;
#line 326 "ssbev_2stage.f"
    } else if (anrm > rmax) {
#line 327 "ssbev_2stage.f"
	iscale = 1;
#line 328 "ssbev_2stage.f"
	sigma = rmax / anrm;
#line 329 "ssbev_2stage.f"
    }
#line 330 "ssbev_2stage.f"
    if (iscale == 1) {
#line 331 "ssbev_2stage.f"
	if (lower) {
#line 332 "ssbev_2stage.f"
	    slascl_("B", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 333 "ssbev_2stage.f"
	} else {
#line 334 "ssbev_2stage.f"
	    slascl_("Q", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 335 "ssbev_2stage.f"
	}
#line 336 "ssbev_2stage.f"
    }

/*     Call SSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form. */

#line 340 "ssbev_2stage.f"
    inde = 1;
#line 341 "ssbev_2stage.f"
    indhous = inde + *n;
#line 342 "ssbev_2stage.f"
    indwrk = indhous + lhtrd;
#line 343 "ssbev_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 345 "ssbev_2stage.f"
    ssytrd_sb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[
	    inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &iinfo, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, call SSTEQR. */

#line 351 "ssbev_2stage.f"
    if (! wantz) {
#line 352 "ssbev_2stage.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 353 "ssbev_2stage.f"
    } else {
#line 354 "ssbev_2stage.f"
	ssteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[
		indwrk], info, (ftnlen)1);
#line 356 "ssbev_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 360 "ssbev_2stage.f"
    if (iscale == 1) {
#line 361 "ssbev_2stage.f"
	if (*info == 0) {
#line 362 "ssbev_2stage.f"
	    imax = *n;
#line 363 "ssbev_2stage.f"
	} else {
#line 364 "ssbev_2stage.f"
	    imax = *info - 1;
#line 365 "ssbev_2stage.f"
	}
#line 366 "ssbev_2stage.f"
	d__1 = 1. / sigma;
#line 366 "ssbev_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 367 "ssbev_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 371 "ssbev_2stage.f"
    work[1] = (doublereal) lwmin;

#line 373 "ssbev_2stage.f"
    return 0;

/*     End of SSBEV_2STAGE */

} /* ssbev_2stage__ */

