#line 1 "ssyevx_2stage.f"
/* ssyevx_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "ssyevx_2stage.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;

/* > \brief <b> SSYEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 SY matrices</b> */

/*  @generated from dsyevx_2stage.f, fortran d -> s, Sat Nov  5 23:55:46 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYEVX_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevx_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevx_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevx_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYEVX_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, */
/*                                 IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
/*                                 LWORK, IWORK, IFAIL, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  Eigenvalues and eigenvectors can be */
/* > selected by specifying either a range of values or a range of indices */
/* > for the desired eigenvalues. */
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
/* > \param[in] RANGE */
/* > \verbatim */
/* >          RANGE is CHARACTER*1 */
/* >          = 'A': all eigenvalues will be found. */
/* >          = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* >                 will be found. */
/* >          = 'I': the IL-th through IU-th eigenvalues will be found. */
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
/* >          On exit, the lower triangle (if UPLO='L') or the upper */
/* >          triangle (if UPLO='U') of A, including the diagonal, is */
/* >          destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is REAL */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* >          ABSTOL is REAL */
/* >          The absolute error tolerance for the eigenvalues. */
/* >          An approximate eigenvalue is accepted as converged */
/* >          when it is determined to lie in an interval [a,b] */
/* >          of width less than or equal to */
/* > */
/* >                  ABSTOL + EPS *   max( |a|,|b| ) , */
/* > */
/* >          where EPS is the machine precision.  If ABSTOL is less than */
/* >          or equal to zero, then  EPS*|T|  will be used in its place, */
/* >          where |T| is the 1-norm of the tridiagonal matrix obtained */
/* >          by reducing A to tridiagonal form. */
/* > */
/* >          Eigenvalues will be computed most accurately when ABSTOL is */
/* >          set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
/* >          If this routine returns with INFO>0, indicating that some */
/* >          eigenvectors did not converge, try setting ABSTOL to */
/* >          2*SLAMCH('S'). */
/* > */
/* >          See "Computing Small Singular Values of Bidiagonal Matrices */
/* >          with Guaranteed High Relative Accuracy," by Demmel and */
/* >          Kahan, LAPACK Working Note #3. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The total number of eigenvalues found.  0 <= M <= N. */
/* >          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, max(1,M)) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          If an eigenvector fails to converge, then that column of Z */
/* >          contains the latest approximation to the eigenvector, and the */
/* >          index of the eigenvector is returned in IFAIL. */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* >          Note: the user must ensure that at least max(1,M) columns are */
/* >          supplied in the array Z; if RANGE = 'V', the exact value of M */
/* >          is not known in advance and an upper bound must be used. */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK. LWORK >= 1, when N <= 1; */
/* >          otherwise */
/* >          If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* >                                   LWORK = MAX(1, 8*N, dimension) where */
/* >                                   dimension = max(stage1,stage2) + (KD+1)*N + 3*N */
/* >                                             = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                                               + max(2*KD*KD, KD*NTHREADS) */
/* >                                               + (KD+1)*N + 3*N */
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
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* >          IFAIL is INTEGER array, dimension (N) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M elements of */
/* >          IFAIL are zero.  If INFO > 0, then IFAIL contains the */
/* >          indices of the eigenvectors that failed to converge. */
/* >          If JOBZ = 'N', then IFAIL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, then i eigenvectors failed to converge. */
/* >                Their indices are stored in array IFAIL. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

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
/* Subroutine */ int ssyevx_2stage__(char *jobz, char *range, char *uplo, 
	integer *n, doublereal *a, integer *lda, doublereal *vl, doublereal *
	vu, integer *il, integer *iu, doublereal *abstol, integer *m, 
	doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *ifail, integer *info, ftnlen 
	jobz_len, ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, ib, kd, jj;
    static doublereal eps, vll, vuu, tmp1;
    static integer indd, inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static integer itmp1, indee;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    static integer lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), ssytrd_2stage__(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static logical wantz, alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp, indiwo, indwkn;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int sstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    ssterf_(integer *, doublereal *, doublereal *, integer *);
    static integer llwrkn, llwork, nsplit;
    static doublereal smlnum;
    extern doublereal slansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    sorgtr_(char *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), sormtr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer indhous;



/*  -- LAPACK driver routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

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

#line 355 "ssyevx_2stage.f"
    /* Parameter adjustments */
#line 355 "ssyevx_2stage.f"
    a_dim1 = *lda;
#line 355 "ssyevx_2stage.f"
    a_offset = 1 + a_dim1;
#line 355 "ssyevx_2stage.f"
    a -= a_offset;
#line 355 "ssyevx_2stage.f"
    --w;
#line 355 "ssyevx_2stage.f"
    z_dim1 = *ldz;
#line 355 "ssyevx_2stage.f"
    z_offset = 1 + z_dim1;
#line 355 "ssyevx_2stage.f"
    z__ -= z_offset;
#line 355 "ssyevx_2stage.f"
    --work;
#line 355 "ssyevx_2stage.f"
    --iwork;
#line 355 "ssyevx_2stage.f"
    --ifail;
#line 355 "ssyevx_2stage.f"

#line 355 "ssyevx_2stage.f"
    /* Function Body */
#line 355 "ssyevx_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 356 "ssyevx_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 357 "ssyevx_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 358 "ssyevx_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 359 "ssyevx_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 360 "ssyevx_2stage.f"
    lquery = *lwork == -1;

#line 362 "ssyevx_2stage.f"
    *info = 0;
#line 363 "ssyevx_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 364 "ssyevx_2stage.f"
	*info = -1;
#line 365 "ssyevx_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 366 "ssyevx_2stage.f"
	*info = -2;
#line 367 "ssyevx_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 368 "ssyevx_2stage.f"
	*info = -3;
#line 369 "ssyevx_2stage.f"
    } else if (*n < 0) {
#line 370 "ssyevx_2stage.f"
	*info = -4;
#line 371 "ssyevx_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 372 "ssyevx_2stage.f"
	*info = -6;
#line 373 "ssyevx_2stage.f"
    } else {
#line 374 "ssyevx_2stage.f"
	if (valeig) {
#line 375 "ssyevx_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 375 "ssyevx_2stage.f"
		*info = -8;
#line 375 "ssyevx_2stage.f"
	    }
#line 377 "ssyevx_2stage.f"
	} else if (indeig) {
#line 378 "ssyevx_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 379 "ssyevx_2stage.f"
		*info = -9;
#line 380 "ssyevx_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 381 "ssyevx_2stage.f"
		*info = -10;
#line 382 "ssyevx_2stage.f"
	    }
#line 383 "ssyevx_2stage.f"
	}
#line 384 "ssyevx_2stage.f"
    }
#line 385 "ssyevx_2stage.f"
    if (*info == 0) {
#line 386 "ssyevx_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 387 "ssyevx_2stage.f"
	    *info = -15;
#line 388 "ssyevx_2stage.f"
	}
#line 389 "ssyevx_2stage.f"
    }

#line 391 "ssyevx_2stage.f"
    if (*info == 0) {
#line 392 "ssyevx_2stage.f"
	if (*n <= 1) {
#line 393 "ssyevx_2stage.f"
	    lwmin = 1;
#line 394 "ssyevx_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 395 "ssyevx_2stage.f"
	} else {
#line 396 "ssyevx_2stage.f"
	    kd = ilaenv2stage_(&c__1, "SSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, 
		    &c_n1, (ftnlen)13, (ftnlen)1);
#line 398 "ssyevx_2stage.f"
	    ib = ilaenv2stage_(&c__2, "SSYTRD_2STAGE", jobz, n, &kd, &c_n1, &
		    c_n1, (ftnlen)13, (ftnlen)1);
#line 400 "ssyevx_2stage.f"
	    lhtrd = ilaenv2stage_(&c__3, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &
		    c_n1, (ftnlen)13, (ftnlen)1);
#line 402 "ssyevx_2stage.f"
	    lwtrd = ilaenv2stage_(&c__4, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &
		    c_n1, (ftnlen)13, (ftnlen)1);
/* Computing MAX */
#line 404 "ssyevx_2stage.f"
	    i__1 = *n << 3, i__2 = *n * 3 + lhtrd + lwtrd;
#line 404 "ssyevx_2stage.f"
	    lwmin = max(i__1,i__2);
#line 405 "ssyevx_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 406 "ssyevx_2stage.f"
	}

#line 408 "ssyevx_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 408 "ssyevx_2stage.f"
	    *info = -17;
#line 408 "ssyevx_2stage.f"
	}
#line 410 "ssyevx_2stage.f"
    }

#line 412 "ssyevx_2stage.f"
    if (*info != 0) {
#line 413 "ssyevx_2stage.f"
	i__1 = -(*info);
#line 413 "ssyevx_2stage.f"
	xerbla_("SSYEVX_2STAGE", &i__1, (ftnlen)13);
#line 414 "ssyevx_2stage.f"
	return 0;
#line 415 "ssyevx_2stage.f"
    } else if (lquery) {
#line 416 "ssyevx_2stage.f"
	return 0;
#line 417 "ssyevx_2stage.f"
    }

/*     Quick return if possible */

#line 421 "ssyevx_2stage.f"
    *m = 0;
#line 422 "ssyevx_2stage.f"
    if (*n == 0) {
#line 423 "ssyevx_2stage.f"
	return 0;
#line 424 "ssyevx_2stage.f"
    }

#line 426 "ssyevx_2stage.f"
    if (*n == 1) {
#line 427 "ssyevx_2stage.f"
	if (alleig || indeig) {
#line 428 "ssyevx_2stage.f"
	    *m = 1;
#line 429 "ssyevx_2stage.f"
	    w[1] = a[a_dim1 + 1];
#line 430 "ssyevx_2stage.f"
	} else {
#line 431 "ssyevx_2stage.f"
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
#line 432 "ssyevx_2stage.f"
		*m = 1;
#line 433 "ssyevx_2stage.f"
		w[1] = a[a_dim1 + 1];
#line 434 "ssyevx_2stage.f"
	    }
#line 435 "ssyevx_2stage.f"
	}
#line 436 "ssyevx_2stage.f"
	if (wantz) {
#line 436 "ssyevx_2stage.f"
	    z__[z_dim1 + 1] = 1.;
#line 436 "ssyevx_2stage.f"
	}
#line 438 "ssyevx_2stage.f"
	return 0;
#line 439 "ssyevx_2stage.f"
    }

/*     Get machine constants. */

#line 443 "ssyevx_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 444 "ssyevx_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 445 "ssyevx_2stage.f"
    smlnum = safmin / eps;
#line 446 "ssyevx_2stage.f"
    bignum = 1. / smlnum;
#line 447 "ssyevx_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 448 "ssyevx_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 448 "ssyevx_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 452 "ssyevx_2stage.f"
    iscale = 0;
#line 453 "ssyevx_2stage.f"
    abstll = *abstol;
#line 454 "ssyevx_2stage.f"
    if (valeig) {
#line 455 "ssyevx_2stage.f"
	vll = *vl;
#line 456 "ssyevx_2stage.f"
	vuu = *vu;
#line 457 "ssyevx_2stage.f"
    }
#line 458 "ssyevx_2stage.f"
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 459 "ssyevx_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 460 "ssyevx_2stage.f"
	iscale = 1;
#line 461 "ssyevx_2stage.f"
	sigma = rmin / anrm;
#line 462 "ssyevx_2stage.f"
    } else if (anrm > rmax) {
#line 463 "ssyevx_2stage.f"
	iscale = 1;
#line 464 "ssyevx_2stage.f"
	sigma = rmax / anrm;
#line 465 "ssyevx_2stage.f"
    }
#line 466 "ssyevx_2stage.f"
    if (iscale == 1) {
#line 467 "ssyevx_2stage.f"
	if (lower) {
#line 468 "ssyevx_2stage.f"
	    i__1 = *n;
#line 468 "ssyevx_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 469 "ssyevx_2stage.f"
		i__2 = *n - j + 1;
#line 469 "ssyevx_2stage.f"
		sscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 470 "ssyevx_2stage.f"
/* L10: */
#line 470 "ssyevx_2stage.f"
	    }
#line 471 "ssyevx_2stage.f"
	} else {
#line 472 "ssyevx_2stage.f"
	    i__1 = *n;
#line 472 "ssyevx_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 473 "ssyevx_2stage.f"
		sscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 474 "ssyevx_2stage.f"
/* L20: */
#line 474 "ssyevx_2stage.f"
	    }
#line 475 "ssyevx_2stage.f"
	}
#line 476 "ssyevx_2stage.f"
	if (*abstol > 0.) {
#line 476 "ssyevx_2stage.f"
	    abstll = *abstol * sigma;
#line 476 "ssyevx_2stage.f"
	}
#line 478 "ssyevx_2stage.f"
	if (valeig) {
#line 479 "ssyevx_2stage.f"
	    vll = *vl * sigma;
#line 480 "ssyevx_2stage.f"
	    vuu = *vu * sigma;
#line 481 "ssyevx_2stage.f"
	}
#line 482 "ssyevx_2stage.f"
    }

/*     Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */

#line 486 "ssyevx_2stage.f"
    indtau = 1;
#line 487 "ssyevx_2stage.f"
    inde = indtau + *n;
#line 488 "ssyevx_2stage.f"
    indd = inde + *n;
#line 489 "ssyevx_2stage.f"
    indhous = indd + *n;
#line 490 "ssyevx_2stage.f"
    indwrk = indhous + lhtrd;
#line 491 "ssyevx_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 493 "ssyevx_2stage.f"
    ssytrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &work[indd], &work[inde]
	    , &work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call SSTERF or SORGTR and SSTEQR.  If this fails for */
/*     some eigenvalue, then try SSTEBZ. */

#line 501 "ssyevx_2stage.f"
    test = FALSE_;
#line 502 "ssyevx_2stage.f"
    if (indeig) {
#line 503 "ssyevx_2stage.f"
	if (*il == 1 && *iu == *n) {
#line 504 "ssyevx_2stage.f"
	    test = TRUE_;
#line 505 "ssyevx_2stage.f"
	}
#line 506 "ssyevx_2stage.f"
    }
#line 507 "ssyevx_2stage.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 508 "ssyevx_2stage.f"
	scopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 509 "ssyevx_2stage.f"
	indee = indwrk + (*n << 1);
#line 510 "ssyevx_2stage.f"
	if (! wantz) {
#line 511 "ssyevx_2stage.f"
	    i__1 = *n - 1;
#line 511 "ssyevx_2stage.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 512 "ssyevx_2stage.f"
	    ssterf_(n, &w[1], &work[indee], info);
#line 513 "ssyevx_2stage.f"
	} else {
#line 514 "ssyevx_2stage.f"
	    slacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 515 "ssyevx_2stage.f"
	    sorgtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
#line 517 "ssyevx_2stage.f"
	    i__1 = *n - 1;
#line 517 "ssyevx_2stage.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 518 "ssyevx_2stage.f"
	    ssteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 520 "ssyevx_2stage.f"
	    if (*info == 0) {
#line 521 "ssyevx_2stage.f"
		i__1 = *n;
#line 521 "ssyevx_2stage.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 522 "ssyevx_2stage.f"
		    ifail[i__] = 0;
#line 523 "ssyevx_2stage.f"
/* L30: */
#line 523 "ssyevx_2stage.f"
		}
#line 524 "ssyevx_2stage.f"
	    }
#line 525 "ssyevx_2stage.f"
	}
#line 526 "ssyevx_2stage.f"
	if (*info == 0) {
#line 527 "ssyevx_2stage.f"
	    *m = *n;
#line 528 "ssyevx_2stage.f"
	    goto L40;
#line 529 "ssyevx_2stage.f"
	}
#line 530 "ssyevx_2stage.f"
	*info = 0;
#line 531 "ssyevx_2stage.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 535 "ssyevx_2stage.f"
    if (wantz) {
#line 536 "ssyevx_2stage.f"
	*(unsigned char *)order = 'B';
#line 537 "ssyevx_2stage.f"
    } else {
#line 538 "ssyevx_2stage.f"
	*(unsigned char *)order = 'E';
#line 539 "ssyevx_2stage.f"
    }
#line 540 "ssyevx_2stage.f"
    indibl = 1;
#line 541 "ssyevx_2stage.f"
    indisp = indibl + *n;
#line 542 "ssyevx_2stage.f"
    indiwo = indisp + *n;
#line 543 "ssyevx_2stage.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 548 "ssyevx_2stage.f"
    if (wantz) {
#line 549 "ssyevx_2stage.f"
	sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEIN. */

#line 556 "ssyevx_2stage.f"
	indwkn = inde;
#line 557 "ssyevx_2stage.f"
	llwrkn = *lwork - indwkn + 1;
#line 558 "ssyevx_2stage.f"
	sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 560 "ssyevx_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 564 "ssyevx_2stage.f"
L40:
#line 565 "ssyevx_2stage.f"
    if (iscale == 1) {
#line 566 "ssyevx_2stage.f"
	if (*info == 0) {
#line 567 "ssyevx_2stage.f"
	    imax = *m;
#line 568 "ssyevx_2stage.f"
	} else {
#line 569 "ssyevx_2stage.f"
	    imax = *info - 1;
#line 570 "ssyevx_2stage.f"
	}
#line 571 "ssyevx_2stage.f"
	d__1 = 1. / sigma;
#line 571 "ssyevx_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 572 "ssyevx_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 577 "ssyevx_2stage.f"
    if (wantz) {
#line 578 "ssyevx_2stage.f"
	i__1 = *m - 1;
#line 578 "ssyevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 579 "ssyevx_2stage.f"
	    i__ = 0;
#line 580 "ssyevx_2stage.f"
	    tmp1 = w[j];
#line 581 "ssyevx_2stage.f"
	    i__2 = *m;
#line 581 "ssyevx_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 582 "ssyevx_2stage.f"
		if (w[jj] < tmp1) {
#line 583 "ssyevx_2stage.f"
		    i__ = jj;
#line 584 "ssyevx_2stage.f"
		    tmp1 = w[jj];
#line 585 "ssyevx_2stage.f"
		}
#line 586 "ssyevx_2stage.f"
/* L50: */
#line 586 "ssyevx_2stage.f"
	    }

#line 588 "ssyevx_2stage.f"
	    if (i__ != 0) {
#line 589 "ssyevx_2stage.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 590 "ssyevx_2stage.f"
		w[i__] = w[j];
#line 591 "ssyevx_2stage.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 592 "ssyevx_2stage.f"
		w[j] = tmp1;
#line 593 "ssyevx_2stage.f"
		iwork[indibl + j - 1] = itmp1;
#line 594 "ssyevx_2stage.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 595 "ssyevx_2stage.f"
		if (*info != 0) {
#line 596 "ssyevx_2stage.f"
		    itmp1 = ifail[i__];
#line 597 "ssyevx_2stage.f"
		    ifail[i__] = ifail[j];
#line 598 "ssyevx_2stage.f"
		    ifail[j] = itmp1;
#line 599 "ssyevx_2stage.f"
		}
#line 600 "ssyevx_2stage.f"
	    }
#line 601 "ssyevx_2stage.f"
/* L60: */
#line 601 "ssyevx_2stage.f"
	}
#line 602 "ssyevx_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 606 "ssyevx_2stage.f"
    work[1] = (doublereal) lwmin;

#line 608 "ssyevx_2stage.f"
    return 0;

/*     End of SSYEVX_2STAGE */

} /* ssyevx_2stage__ */

