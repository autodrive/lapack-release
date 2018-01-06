#line 1 "dsbevx_2stage.f"
/* dsbevx_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsbevx_2stage.f"
/* Table of constant values */

static integer c__18 = 18;
static integer c_n1 = -1;
static integer c__19 = 19;
static integer c__20 = 20;
static doublereal c_b24 = 1.;
static integer c__1 = 1;
static doublereal c_b45 = 0.;

/* > \brief <b> DSBEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 OTHER matrices</b> */

/*  @precisions fortran d -> s */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBEVX_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbevx_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbevx_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbevx_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, */
/*                                 LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, */
/*                                 LDZ, WORK, LWORK, IWORK, IFAIL, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric band matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. Eigenvalues and eigenvectors can */
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
/* >          = 'A': all eigenvalues will be found; */
/* >          = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* >                 will be found; */
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
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* >          If JOBZ = 'V', the N-by-N orthogonal matrix used in the */
/* >                         reduction to tridiagonal form. */
/* >          If JOBZ = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  If JOBZ = 'V', then */
/* >          LDQ >= max(1,N). */
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
/* >          by reducing AB to tridiagonal form. */
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
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK. LWORK >= 1, when N <= 1; */
/* >          otherwise */
/* >          If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* >                                   LWORK = MAX(1, 7*N, dimension) where */
/* >                                   dimension = (2KD+1)*N + KD*NTHREADS + 2*N */
/* >                                   where KD is the size of the band. */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
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
/* Subroutine */ int dsbevx_2stage__(char *jobz, char *range, char *uplo, 
	integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *q,
	 integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *ifail, integer *info, ftnlen jobz_len, ftnlen range_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, 
	    i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, ib, jj;
    static doublereal eps, vll, vuu, tmp1;
    static integer indd, inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    extern /* Subroutine */ int dsytrd_sb2st__(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer itmp1, indee;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer iinfo;
    static char order[1];
    static integer lhtrd;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern doublereal dlansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    static logical valeig;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indisp;
    extern /* Subroutine */ int dstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dsterf_(integer *, doublereal *, doublereal *, integer *);
    static integer indiwo;
    extern /* Subroutine */ int dstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer nsplit, llwork;
    static doublereal smlnum;
    static logical lquery;
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

#line 378 "dsbevx_2stage.f"
    /* Parameter adjustments */
#line 378 "dsbevx_2stage.f"
    ab_dim1 = *ldab;
#line 378 "dsbevx_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 378 "dsbevx_2stage.f"
    ab -= ab_offset;
#line 378 "dsbevx_2stage.f"
    q_dim1 = *ldq;
#line 378 "dsbevx_2stage.f"
    q_offset = 1 + q_dim1;
#line 378 "dsbevx_2stage.f"
    q -= q_offset;
#line 378 "dsbevx_2stage.f"
    --w;
#line 378 "dsbevx_2stage.f"
    z_dim1 = *ldz;
#line 378 "dsbevx_2stage.f"
    z_offset = 1 + z_dim1;
#line 378 "dsbevx_2stage.f"
    z__ -= z_offset;
#line 378 "dsbevx_2stage.f"
    --work;
#line 378 "dsbevx_2stage.f"
    --iwork;
#line 378 "dsbevx_2stage.f"
    --ifail;
#line 378 "dsbevx_2stage.f"

#line 378 "dsbevx_2stage.f"
    /* Function Body */
#line 378 "dsbevx_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 379 "dsbevx_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 380 "dsbevx_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 381 "dsbevx_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 382 "dsbevx_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 383 "dsbevx_2stage.f"
    lquery = *lwork == -1;

#line 385 "dsbevx_2stage.f"
    *info = 0;
#line 386 "dsbevx_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 387 "dsbevx_2stage.f"
	*info = -1;
#line 388 "dsbevx_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 389 "dsbevx_2stage.f"
	*info = -2;
#line 390 "dsbevx_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 391 "dsbevx_2stage.f"
	*info = -3;
#line 392 "dsbevx_2stage.f"
    } else if (*n < 0) {
#line 393 "dsbevx_2stage.f"
	*info = -4;
#line 394 "dsbevx_2stage.f"
    } else if (*kd < 0) {
#line 395 "dsbevx_2stage.f"
	*info = -5;
#line 396 "dsbevx_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 397 "dsbevx_2stage.f"
	*info = -7;
#line 398 "dsbevx_2stage.f"
    } else if (wantz && *ldq < max(1,*n)) {
#line 399 "dsbevx_2stage.f"
	*info = -9;
#line 400 "dsbevx_2stage.f"
    } else {
#line 401 "dsbevx_2stage.f"
	if (valeig) {
#line 402 "dsbevx_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 402 "dsbevx_2stage.f"
		*info = -11;
#line 402 "dsbevx_2stage.f"
	    }
#line 404 "dsbevx_2stage.f"
	} else if (indeig) {
#line 405 "dsbevx_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 406 "dsbevx_2stage.f"
		*info = -12;
#line 407 "dsbevx_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 408 "dsbevx_2stage.f"
		*info = -13;
#line 409 "dsbevx_2stage.f"
	    }
#line 410 "dsbevx_2stage.f"
	}
#line 411 "dsbevx_2stage.f"
    }
#line 412 "dsbevx_2stage.f"
    if (*info == 0) {
#line 413 "dsbevx_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 413 "dsbevx_2stage.f"
	    *info = -18;
#line 413 "dsbevx_2stage.f"
	}
#line 415 "dsbevx_2stage.f"
    }

#line 417 "dsbevx_2stage.f"
    if (*info == 0) {
#line 418 "dsbevx_2stage.f"
	if (*n <= 1) {
#line 419 "dsbevx_2stage.f"
	    lwmin = 1;
#line 420 "dsbevx_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 421 "dsbevx_2stage.f"
	} else {
#line 422 "dsbevx_2stage.f"
	    ib = ilaenv_(&c__18, "DSYTRD_SB2ST", jobz, n, kd, &c_n1, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 423 "dsbevx_2stage.f"
	    lhtrd = ilaenv_(&c__19, "DSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 424 "dsbevx_2stage.f"
	    lwtrd = ilaenv_(&c__20, "DSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 425 "dsbevx_2stage.f"
	    lwmin = (*n << 1) + lhtrd + lwtrd;
#line 426 "dsbevx_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 427 "dsbevx_2stage.f"
	}

#line 429 "dsbevx_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 429 "dsbevx_2stage.f"
	    *info = -20;
#line 429 "dsbevx_2stage.f"
	}
#line 431 "dsbevx_2stage.f"
    }

#line 433 "dsbevx_2stage.f"
    if (*info != 0) {
#line 434 "dsbevx_2stage.f"
	i__1 = -(*info);
#line 434 "dsbevx_2stage.f"
	xerbla_("DSBEVX_2STAGE ", &i__1, (ftnlen)14);
#line 435 "dsbevx_2stage.f"
	return 0;
#line 436 "dsbevx_2stage.f"
    } else if (lquery) {
#line 437 "dsbevx_2stage.f"
	return 0;
#line 438 "dsbevx_2stage.f"
    }

/*     Quick return if possible */

#line 442 "dsbevx_2stage.f"
    *m = 0;
#line 443 "dsbevx_2stage.f"
    if (*n == 0) {
#line 443 "dsbevx_2stage.f"
	return 0;
#line 443 "dsbevx_2stage.f"
    }

#line 446 "dsbevx_2stage.f"
    if (*n == 1) {
#line 447 "dsbevx_2stage.f"
	*m = 1;
#line 448 "dsbevx_2stage.f"
	if (lower) {
#line 449 "dsbevx_2stage.f"
	    tmp1 = ab[ab_dim1 + 1];
#line 450 "dsbevx_2stage.f"
	} else {
#line 451 "dsbevx_2stage.f"
	    tmp1 = ab[*kd + 1 + ab_dim1];
#line 452 "dsbevx_2stage.f"
	}
#line 453 "dsbevx_2stage.f"
	if (valeig) {
#line 454 "dsbevx_2stage.f"
	    if (! (*vl < tmp1 && *vu >= tmp1)) {
#line 454 "dsbevx_2stage.f"
		*m = 0;
#line 454 "dsbevx_2stage.f"
	    }
#line 456 "dsbevx_2stage.f"
	}
#line 457 "dsbevx_2stage.f"
	if (*m == 1) {
#line 458 "dsbevx_2stage.f"
	    w[1] = tmp1;
#line 459 "dsbevx_2stage.f"
	    if (wantz) {
#line 459 "dsbevx_2stage.f"
		z__[z_dim1 + 1] = 1.;
#line 459 "dsbevx_2stage.f"
	    }
#line 461 "dsbevx_2stage.f"
	}
#line 462 "dsbevx_2stage.f"
	return 0;
#line 463 "dsbevx_2stage.f"
    }

/*     Get machine constants. */

#line 467 "dsbevx_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 468 "dsbevx_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 469 "dsbevx_2stage.f"
    smlnum = safmin / eps;
#line 470 "dsbevx_2stage.f"
    bignum = 1. / smlnum;
#line 471 "dsbevx_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 472 "dsbevx_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 472 "dsbevx_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 476 "dsbevx_2stage.f"
    iscale = 0;
#line 477 "dsbevx_2stage.f"
    abstll = *abstol;
#line 478 "dsbevx_2stage.f"
    if (valeig) {
#line 479 "dsbevx_2stage.f"
	vll = *vl;
#line 480 "dsbevx_2stage.f"
	vuu = *vu;
#line 481 "dsbevx_2stage.f"
    } else {
#line 482 "dsbevx_2stage.f"
	vll = 0.;
#line 483 "dsbevx_2stage.f"
	vuu = 0.;
#line 484 "dsbevx_2stage.f"
    }
#line 485 "dsbevx_2stage.f"
    anrm = dlansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 486 "dsbevx_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 487 "dsbevx_2stage.f"
	iscale = 1;
#line 488 "dsbevx_2stage.f"
	sigma = rmin / anrm;
#line 489 "dsbevx_2stage.f"
    } else if (anrm > rmax) {
#line 490 "dsbevx_2stage.f"
	iscale = 1;
#line 491 "dsbevx_2stage.f"
	sigma = rmax / anrm;
#line 492 "dsbevx_2stage.f"
    }
#line 493 "dsbevx_2stage.f"
    if (iscale == 1) {
#line 494 "dsbevx_2stage.f"
	if (lower) {
#line 495 "dsbevx_2stage.f"
	    dlascl_("B", kd, kd, &c_b24, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 496 "dsbevx_2stage.f"
	} else {
#line 497 "dsbevx_2stage.f"
	    dlascl_("Q", kd, kd, &c_b24, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 498 "dsbevx_2stage.f"
	}
#line 499 "dsbevx_2stage.f"
	if (*abstol > 0.) {
#line 499 "dsbevx_2stage.f"
	    abstll = *abstol * sigma;
#line 499 "dsbevx_2stage.f"
	}
#line 501 "dsbevx_2stage.f"
	if (valeig) {
#line 502 "dsbevx_2stage.f"
	    vll = *vl * sigma;
#line 503 "dsbevx_2stage.f"
	    vuu = *vu * sigma;
#line 504 "dsbevx_2stage.f"
	}
#line 505 "dsbevx_2stage.f"
    }

/*     Call DSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form. */

#line 509 "dsbevx_2stage.f"
    indd = 1;
#line 510 "dsbevx_2stage.f"
    inde = indd + *n;
#line 511 "dsbevx_2stage.f"
    indhous = inde + *n;
#line 512 "dsbevx_2stage.f"
    indwrk = indhous + lhtrd;
#line 513 "dsbevx_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 515 "dsbevx_2stage.f"
    dsytrd_sb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &work[indd], 
	    &work[inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call DSTERF or SSTEQR.  If this fails for some */
/*     eigenvalue, then try DSTEBZ. */

#line 523 "dsbevx_2stage.f"
    test = FALSE_;
#line 524 "dsbevx_2stage.f"
    if (indeig) {
#line 525 "dsbevx_2stage.f"
	if (*il == 1 && *iu == *n) {
#line 526 "dsbevx_2stage.f"
	    test = TRUE_;
#line 527 "dsbevx_2stage.f"
	}
#line 528 "dsbevx_2stage.f"
    }
#line 529 "dsbevx_2stage.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 530 "dsbevx_2stage.f"
	dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 531 "dsbevx_2stage.f"
	indee = indwrk + (*n << 1);
#line 532 "dsbevx_2stage.f"
	if (! wantz) {
#line 533 "dsbevx_2stage.f"
	    i__1 = *n - 1;
#line 533 "dsbevx_2stage.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 534 "dsbevx_2stage.f"
	    dsterf_(n, &w[1], &work[indee], info);
#line 535 "dsbevx_2stage.f"
	} else {
#line 536 "dsbevx_2stage.f"
	    dlacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 537 "dsbevx_2stage.f"
	    i__1 = *n - 1;
#line 537 "dsbevx_2stage.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 538 "dsbevx_2stage.f"
	    dsteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 540 "dsbevx_2stage.f"
	    if (*info == 0) {
#line 541 "dsbevx_2stage.f"
		i__1 = *n;
#line 541 "dsbevx_2stage.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 542 "dsbevx_2stage.f"
		    ifail[i__] = 0;
#line 543 "dsbevx_2stage.f"
/* L10: */
#line 543 "dsbevx_2stage.f"
		}
#line 544 "dsbevx_2stage.f"
	    }
#line 545 "dsbevx_2stage.f"
	}
#line 546 "dsbevx_2stage.f"
	if (*info == 0) {
#line 547 "dsbevx_2stage.f"
	    *m = *n;
#line 548 "dsbevx_2stage.f"
	    goto L30;
#line 549 "dsbevx_2stage.f"
	}
#line 550 "dsbevx_2stage.f"
	*info = 0;
#line 551 "dsbevx_2stage.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 555 "dsbevx_2stage.f"
    if (wantz) {
#line 556 "dsbevx_2stage.f"
	*(unsigned char *)order = 'B';
#line 557 "dsbevx_2stage.f"
    } else {
#line 558 "dsbevx_2stage.f"
	*(unsigned char *)order = 'E';
#line 559 "dsbevx_2stage.f"
    }
#line 560 "dsbevx_2stage.f"
    indibl = 1;
#line 561 "dsbevx_2stage.f"
    indisp = indibl + *n;
#line 562 "dsbevx_2stage.f"
    indiwo = indisp + *n;
#line 563 "dsbevx_2stage.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 568 "dsbevx_2stage.f"
    if (wantz) {
#line 569 "dsbevx_2stage.f"
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

#line 576 "dsbevx_2stage.f"
	i__1 = *m;
#line 576 "dsbevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 577 "dsbevx_2stage.f"
	    dcopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 578 "dsbevx_2stage.f"
	    dgemv_("N", n, n, &c_b24, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b45, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 580 "dsbevx_2stage.f"
/* L20: */
#line 580 "dsbevx_2stage.f"
	}
#line 581 "dsbevx_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 585 "dsbevx_2stage.f"
L30:
#line 586 "dsbevx_2stage.f"
    if (iscale == 1) {
#line 587 "dsbevx_2stage.f"
	if (*info == 0) {
#line 588 "dsbevx_2stage.f"
	    imax = *m;
#line 589 "dsbevx_2stage.f"
	} else {
#line 590 "dsbevx_2stage.f"
	    imax = *info - 1;
#line 591 "dsbevx_2stage.f"
	}
#line 592 "dsbevx_2stage.f"
	d__1 = 1. / sigma;
#line 592 "dsbevx_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 593 "dsbevx_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 598 "dsbevx_2stage.f"
    if (wantz) {
#line 599 "dsbevx_2stage.f"
	i__1 = *m - 1;
#line 599 "dsbevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 600 "dsbevx_2stage.f"
	    i__ = 0;
#line 601 "dsbevx_2stage.f"
	    tmp1 = w[j];
#line 602 "dsbevx_2stage.f"
	    i__2 = *m;
#line 602 "dsbevx_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 603 "dsbevx_2stage.f"
		if (w[jj] < tmp1) {
#line 604 "dsbevx_2stage.f"
		    i__ = jj;
#line 605 "dsbevx_2stage.f"
		    tmp1 = w[jj];
#line 606 "dsbevx_2stage.f"
		}
#line 607 "dsbevx_2stage.f"
/* L40: */
#line 607 "dsbevx_2stage.f"
	    }

#line 609 "dsbevx_2stage.f"
	    if (i__ != 0) {
#line 610 "dsbevx_2stage.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 611 "dsbevx_2stage.f"
		w[i__] = w[j];
#line 612 "dsbevx_2stage.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 613 "dsbevx_2stage.f"
		w[j] = tmp1;
#line 614 "dsbevx_2stage.f"
		iwork[indibl + j - 1] = itmp1;
#line 615 "dsbevx_2stage.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 616 "dsbevx_2stage.f"
		if (*info != 0) {
#line 617 "dsbevx_2stage.f"
		    itmp1 = ifail[i__];
#line 618 "dsbevx_2stage.f"
		    ifail[i__] = ifail[j];
#line 619 "dsbevx_2stage.f"
		    ifail[j] = itmp1;
#line 620 "dsbevx_2stage.f"
		}
#line 621 "dsbevx_2stage.f"
	    }
#line 622 "dsbevx_2stage.f"
/* L50: */
#line 622 "dsbevx_2stage.f"
	}
#line 623 "dsbevx_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 627 "dsbevx_2stage.f"
    work[1] = (doublereal) lwmin;

#line 629 "dsbevx_2stage.f"
    return 0;

/*     End of DSBEVX_2STAGE */

} /* dsbevx_2stage__ */

