#line 1 "zheevx_2stage.f"
/* zheevx_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "zheevx_2stage.f"
/* Table of constant values */

static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__1 = 1;

/* > \brief <b> ZHEEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 HE matrices</b> */

/*  @precisions fortran z -> s d c */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEVX_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevx_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevx_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevx_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEVX_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, */
/*                                 IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
/*                                 LWORK, RWORK, IWORK, IFAIL, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  Eigenvalues and eigenvectors can */
/* > be selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
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
/* >          VL is DOUBLE PRECISION */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
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
/* >          ABSTOL is DOUBLE PRECISION */
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
/* >          set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
/* >          If this routine returns with INFO>0, indicating that some */
/* >          eigenvectors did not converge, try setting ABSTOL to */
/* >          2*DLAMCH('S'). */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          On normal exit, the first M elements contain the selected */
/* >          eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, max(1,M)) */
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
/* >                                   LWORK = MAX(1, 8*N, dimension) where */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (7*N) */
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
/* Subroutine */ int zheevx_2stage__(char *jobz, char *range, char *uplo, 
	integer *n, doublecomplex *a, integer *lda, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info, ftnlen jobz_len, ftnlen range_len, ftnlen 
	uplo_len)
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
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static integer itmp1, indee;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zhetrd_2stage__(char *, char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer iinfo;
    static char order[1];
    static integer lhtrd;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal abstll, bignum;
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer indiwk, indisp, indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dstebz_(char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer nsplit, llwork;
    static doublereal smlnum;
    extern /* Subroutine */ int zstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static logical lquery;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), zungtr_(char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen), 
	    zunmtr_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer indhous;



/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 364 "zheevx_2stage.f"
    /* Parameter adjustments */
#line 364 "zheevx_2stage.f"
    a_dim1 = *lda;
#line 364 "zheevx_2stage.f"
    a_offset = 1 + a_dim1;
#line 364 "zheevx_2stage.f"
    a -= a_offset;
#line 364 "zheevx_2stage.f"
    --w;
#line 364 "zheevx_2stage.f"
    z_dim1 = *ldz;
#line 364 "zheevx_2stage.f"
    z_offset = 1 + z_dim1;
#line 364 "zheevx_2stage.f"
    z__ -= z_offset;
#line 364 "zheevx_2stage.f"
    --work;
#line 364 "zheevx_2stage.f"
    --rwork;
#line 364 "zheevx_2stage.f"
    --iwork;
#line 364 "zheevx_2stage.f"
    --ifail;
#line 364 "zheevx_2stage.f"

#line 364 "zheevx_2stage.f"
    /* Function Body */
#line 364 "zheevx_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 365 "zheevx_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 366 "zheevx_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 367 "zheevx_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 368 "zheevx_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 369 "zheevx_2stage.f"
    lquery = *lwork == -1;

#line 371 "zheevx_2stage.f"
    *info = 0;
#line 372 "zheevx_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 373 "zheevx_2stage.f"
	*info = -1;
#line 374 "zheevx_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 375 "zheevx_2stage.f"
	*info = -2;
#line 376 "zheevx_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 377 "zheevx_2stage.f"
	*info = -3;
#line 378 "zheevx_2stage.f"
    } else if (*n < 0) {
#line 379 "zheevx_2stage.f"
	*info = -4;
#line 380 "zheevx_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 381 "zheevx_2stage.f"
	*info = -6;
#line 382 "zheevx_2stage.f"
    } else {
#line 383 "zheevx_2stage.f"
	if (valeig) {
#line 384 "zheevx_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 384 "zheevx_2stage.f"
		*info = -8;
#line 384 "zheevx_2stage.f"
	    }
#line 386 "zheevx_2stage.f"
	} else if (indeig) {
#line 387 "zheevx_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 388 "zheevx_2stage.f"
		*info = -9;
#line 389 "zheevx_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 390 "zheevx_2stage.f"
		*info = -10;
#line 391 "zheevx_2stage.f"
	    }
#line 392 "zheevx_2stage.f"
	}
#line 393 "zheevx_2stage.f"
    }
#line 394 "zheevx_2stage.f"
    if (*info == 0) {
#line 395 "zheevx_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 396 "zheevx_2stage.f"
	    *info = -15;
#line 397 "zheevx_2stage.f"
	}
#line 398 "zheevx_2stage.f"
    }

#line 400 "zheevx_2stage.f"
    if (*info == 0) {
#line 401 "zheevx_2stage.f"
	if (*n <= 1) {
#line 402 "zheevx_2stage.f"
	    lwmin = 1;
#line 403 "zheevx_2stage.f"
	    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 404 "zheevx_2stage.f"
	} else {
#line 405 "zheevx_2stage.f"
	    kd = ilaenv_(&c__17, "ZHETRD_2STAGE", jobz, n, &c_n1, &c_n1, &
		    c_n1, (ftnlen)13, (ftnlen)1);
#line 406 "zheevx_2stage.f"
	    ib = ilaenv_(&c__18, "ZHETRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, 
		    (ftnlen)13, (ftnlen)1);
#line 407 "zheevx_2stage.f"
	    lhtrd = ilaenv_(&c__19, "ZHETRD_2STAGE", jobz, n, &kd, &ib, &c_n1,
		     (ftnlen)13, (ftnlen)1);
#line 408 "zheevx_2stage.f"
	    lwtrd = ilaenv_(&c__20, "ZHETRD_2STAGE", jobz, n, &kd, &ib, &c_n1,
		     (ftnlen)13, (ftnlen)1);
#line 409 "zheevx_2stage.f"
	    lwmin = *n + lhtrd + lwtrd;
#line 410 "zheevx_2stage.f"
	    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 411 "zheevx_2stage.f"
	}

#line 413 "zheevx_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 413 "zheevx_2stage.f"
	    *info = -17;
#line 413 "zheevx_2stage.f"
	}
#line 415 "zheevx_2stage.f"
    }

#line 417 "zheevx_2stage.f"
    if (*info != 0) {
#line 418 "zheevx_2stage.f"
	i__1 = -(*info);
#line 418 "zheevx_2stage.f"
	xerbla_("ZHEEVX_2STAGE", &i__1, (ftnlen)13);
#line 419 "zheevx_2stage.f"
	return 0;
#line 420 "zheevx_2stage.f"
    } else if (lquery) {
#line 421 "zheevx_2stage.f"
	return 0;
#line 422 "zheevx_2stage.f"
    }

/*     Quick return if possible */

#line 426 "zheevx_2stage.f"
    *m = 0;
#line 427 "zheevx_2stage.f"
    if (*n == 0) {
#line 428 "zheevx_2stage.f"
	return 0;
#line 429 "zheevx_2stage.f"
    }

#line 431 "zheevx_2stage.f"
    if (*n == 1) {
#line 432 "zheevx_2stage.f"
	if (alleig || indeig) {
#line 433 "zheevx_2stage.f"
	    *m = 1;
#line 434 "zheevx_2stage.f"
	    i__1 = a_dim1 + 1;
#line 434 "zheevx_2stage.f"
	    w[1] = a[i__1].r;
#line 435 "zheevx_2stage.f"
	} else if (valeig) {
#line 436 "zheevx_2stage.f"
	    i__1 = a_dim1 + 1;
#line 436 "zheevx_2stage.f"
	    i__2 = a_dim1 + 1;
#line 436 "zheevx_2stage.f"
	    if (*vl < a[i__1].r && *vu >= a[i__2].r) {
#line 438 "zheevx_2stage.f"
		*m = 1;
#line 439 "zheevx_2stage.f"
		i__1 = a_dim1 + 1;
#line 439 "zheevx_2stage.f"
		w[1] = a[i__1].r;
#line 440 "zheevx_2stage.f"
	    }
#line 441 "zheevx_2stage.f"
	}
#line 442 "zheevx_2stage.f"
	if (wantz) {
#line 442 "zheevx_2stage.f"
	    i__1 = z_dim1 + 1;
#line 442 "zheevx_2stage.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 442 "zheevx_2stage.f"
	}
#line 444 "zheevx_2stage.f"
	return 0;
#line 445 "zheevx_2stage.f"
    }

/*     Get machine constants. */

#line 449 "zheevx_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 450 "zheevx_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 451 "zheevx_2stage.f"
    smlnum = safmin / eps;
#line 452 "zheevx_2stage.f"
    bignum = 1. / smlnum;
#line 453 "zheevx_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 454 "zheevx_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 454 "zheevx_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 458 "zheevx_2stage.f"
    iscale = 0;
#line 459 "zheevx_2stage.f"
    abstll = *abstol;
#line 460 "zheevx_2stage.f"
    if (valeig) {
#line 461 "zheevx_2stage.f"
	vll = *vl;
#line 462 "zheevx_2stage.f"
	vuu = *vu;
#line 463 "zheevx_2stage.f"
    }
#line 464 "zheevx_2stage.f"
    anrm = zlanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 465 "zheevx_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 466 "zheevx_2stage.f"
	iscale = 1;
#line 467 "zheevx_2stage.f"
	sigma = rmin / anrm;
#line 468 "zheevx_2stage.f"
    } else if (anrm > rmax) {
#line 469 "zheevx_2stage.f"
	iscale = 1;
#line 470 "zheevx_2stage.f"
	sigma = rmax / anrm;
#line 471 "zheevx_2stage.f"
    }
#line 472 "zheevx_2stage.f"
    if (iscale == 1) {
#line 473 "zheevx_2stage.f"
	if (lower) {
#line 474 "zheevx_2stage.f"
	    i__1 = *n;
#line 474 "zheevx_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 475 "zheevx_2stage.f"
		i__2 = *n - j + 1;
#line 475 "zheevx_2stage.f"
		zdscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 476 "zheevx_2stage.f"
/* L10: */
#line 476 "zheevx_2stage.f"
	    }
#line 477 "zheevx_2stage.f"
	} else {
#line 478 "zheevx_2stage.f"
	    i__1 = *n;
#line 478 "zheevx_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 479 "zheevx_2stage.f"
		zdscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 480 "zheevx_2stage.f"
/* L20: */
#line 480 "zheevx_2stage.f"
	    }
#line 481 "zheevx_2stage.f"
	}
#line 482 "zheevx_2stage.f"
	if (*abstol > 0.) {
#line 482 "zheevx_2stage.f"
	    abstll = *abstol * sigma;
#line 482 "zheevx_2stage.f"
	}
#line 484 "zheevx_2stage.f"
	if (valeig) {
#line 485 "zheevx_2stage.f"
	    vll = *vl * sigma;
#line 486 "zheevx_2stage.f"
	    vuu = *vu * sigma;
#line 487 "zheevx_2stage.f"
	}
#line 488 "zheevx_2stage.f"
    }

/*     Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form. */

#line 492 "zheevx_2stage.f"
    indd = 1;
#line 493 "zheevx_2stage.f"
    inde = indd + *n;
#line 494 "zheevx_2stage.f"
    indrwk = inde + *n;
#line 495 "zheevx_2stage.f"
    indtau = 1;
#line 496 "zheevx_2stage.f"
    indhous = indtau + *n;
#line 497 "zheevx_2stage.f"
    indwrk = indhous + lhtrd;
#line 498 "zheevx_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 500 "zheevx_2stage.f"
    zhetrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &rwork[indd], &rwork[
	    inde], &work[indtau], &work[indhous], &lhtrd, &work[indwrk], &
	    llwork, &iinfo, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call DSTERF or ZUNGTR and ZSTEQR.  If this fails for */
/*     some eigenvalue, then try DSTEBZ. */

#line 509 "zheevx_2stage.f"
    test = FALSE_;
#line 510 "zheevx_2stage.f"
    if (indeig) {
#line 511 "zheevx_2stage.f"
	if (*il == 1 && *iu == *n) {
#line 512 "zheevx_2stage.f"
	    test = TRUE_;
#line 513 "zheevx_2stage.f"
	}
#line 514 "zheevx_2stage.f"
    }
#line 515 "zheevx_2stage.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 516 "zheevx_2stage.f"
	dcopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 517 "zheevx_2stage.f"
	indee = indrwk + (*n << 1);
#line 518 "zheevx_2stage.f"
	if (! wantz) {
#line 519 "zheevx_2stage.f"
	    i__1 = *n - 1;
#line 519 "zheevx_2stage.f"
	    dcopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 520 "zheevx_2stage.f"
	    dsterf_(n, &w[1], &rwork[indee], info);
#line 521 "zheevx_2stage.f"
	} else {
#line 522 "zheevx_2stage.f"
	    zlacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 523 "zheevx_2stage.f"
	    zungtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
#line 525 "zheevx_2stage.f"
	    i__1 = *n - 1;
#line 525 "zheevx_2stage.f"
	    dcopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 526 "zheevx_2stage.f"
	    zsteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 528 "zheevx_2stage.f"
	    if (*info == 0) {
#line 529 "zheevx_2stage.f"
		i__1 = *n;
#line 529 "zheevx_2stage.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 530 "zheevx_2stage.f"
		    ifail[i__] = 0;
#line 531 "zheevx_2stage.f"
/* L30: */
#line 531 "zheevx_2stage.f"
		}
#line 532 "zheevx_2stage.f"
	    }
#line 533 "zheevx_2stage.f"
	}
#line 534 "zheevx_2stage.f"
	if (*info == 0) {
#line 535 "zheevx_2stage.f"
	    *m = *n;
#line 536 "zheevx_2stage.f"
	    goto L40;
#line 537 "zheevx_2stage.f"
	}
#line 538 "zheevx_2stage.f"
	*info = 0;
#line 539 "zheevx_2stage.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN. */

#line 543 "zheevx_2stage.f"
    if (wantz) {
#line 544 "zheevx_2stage.f"
	*(unsigned char *)order = 'B';
#line 545 "zheevx_2stage.f"
    } else {
#line 546 "zheevx_2stage.f"
	*(unsigned char *)order = 'E';
#line 547 "zheevx_2stage.f"
    }
#line 548 "zheevx_2stage.f"
    indibl = 1;
#line 549 "zheevx_2stage.f"
    indisp = indibl + *n;
#line 550 "zheevx_2stage.f"
    indiwk = indisp + *n;
#line 551 "zheevx_2stage.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 556 "zheevx_2stage.f"
    if (wantz) {
#line 557 "zheevx_2stage.f"
	zstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by ZSTEIN. */

#line 564 "zheevx_2stage.f"
	zunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwrk], &llwork, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 566 "zheevx_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 570 "zheevx_2stage.f"
L40:
#line 571 "zheevx_2stage.f"
    if (iscale == 1) {
#line 572 "zheevx_2stage.f"
	if (*info == 0) {
#line 573 "zheevx_2stage.f"
	    imax = *m;
#line 574 "zheevx_2stage.f"
	} else {
#line 575 "zheevx_2stage.f"
	    imax = *info - 1;
#line 576 "zheevx_2stage.f"
	}
#line 577 "zheevx_2stage.f"
	d__1 = 1. / sigma;
#line 577 "zheevx_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 578 "zheevx_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 583 "zheevx_2stage.f"
    if (wantz) {
#line 584 "zheevx_2stage.f"
	i__1 = *m - 1;
#line 584 "zheevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 585 "zheevx_2stage.f"
	    i__ = 0;
#line 586 "zheevx_2stage.f"
	    tmp1 = w[j];
#line 587 "zheevx_2stage.f"
	    i__2 = *m;
#line 587 "zheevx_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 588 "zheevx_2stage.f"
		if (w[jj] < tmp1) {
#line 589 "zheevx_2stage.f"
		    i__ = jj;
#line 590 "zheevx_2stage.f"
		    tmp1 = w[jj];
#line 591 "zheevx_2stage.f"
		}
#line 592 "zheevx_2stage.f"
/* L50: */
#line 592 "zheevx_2stage.f"
	    }

#line 594 "zheevx_2stage.f"
	    if (i__ != 0) {
#line 595 "zheevx_2stage.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 596 "zheevx_2stage.f"
		w[i__] = w[j];
#line 597 "zheevx_2stage.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 598 "zheevx_2stage.f"
		w[j] = tmp1;
#line 599 "zheevx_2stage.f"
		iwork[indibl + j - 1] = itmp1;
#line 600 "zheevx_2stage.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 601 "zheevx_2stage.f"
		if (*info != 0) {
#line 602 "zheevx_2stage.f"
		    itmp1 = ifail[i__];
#line 603 "zheevx_2stage.f"
		    ifail[i__] = ifail[j];
#line 604 "zheevx_2stage.f"
		    ifail[j] = itmp1;
#line 605 "zheevx_2stage.f"
		}
#line 606 "zheevx_2stage.f"
	    }
#line 607 "zheevx_2stage.f"
/* L60: */
#line 607 "zheevx_2stage.f"
	}
#line 608 "zheevx_2stage.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 612 "zheevx_2stage.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 614 "zheevx_2stage.f"
    return 0;

/*     End of ZHEEVX_2STAGE */

} /* zheevx_2stage__ */

