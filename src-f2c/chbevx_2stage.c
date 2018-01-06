#line 1 "chbevx_2stage.f"
/* chbevx_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "chbevx_2stage.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__18 = 18;
static integer c_n1 = -1;
static integer c__19 = 19;
static integer c__20 = 20;
static doublereal c_b26 = 1.;
static integer c__1 = 1;

/* > \brief <b> CHBEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 OTHER matrices</b> */

/*  @generated from zhbevx_2stage.f, fortran z -> c, Sat Nov  5 23:18:22 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBEVX_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevx_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevx_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevx_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, */
/*                                 Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, */
/*                                 Z, LDZ, WORK, LWORK, RWORK, IWORK, */
/*                                 IFAIL, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian band matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  Eigenvalues and eigenvectors */
/* > can be selected by specifying either a range of values or a range of */
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
/* >          AB is COMPLEX array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, AB is overwritten by values generated during the */
/* >          reduction to tridiagonal form. */
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
/* >          Q is COMPLEX array, dimension (LDQ, N) */
/* >          If JOBZ = 'V', the N-by-N unitary matrix used in the */
/* >                          reduction to tridiagonal form. */
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
/* >          by reducing AB to tridiagonal form. */
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
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, max(1,M)) */
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
/* >          WORK is COMPLEX array, dimension (LWORK) */
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
/* >          RWORK is REAL array, dimension (7*N) */
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
/* Subroutine */ int chbevx_2stage__(char *jobz, char *range, char *uplo, 
	integer *n, integer *kd, doublecomplex *ab, integer *ldab, 
	doublecomplex *q, integer *ldq, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *iwork, integer *ifail, integer *
	info, ftnlen jobz_len, ftnlen range_len, ftnlen uplo_len)
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
    extern /* Subroutine */ int chetrd_hb2st__(char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer indd, inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static doublecomplex ctmp1;
    static integer itmp1, indee;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    static integer lhtrd;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static integer lwmin;
    static logical lower;
    static integer lwtrd;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantz;
    extern doublereal clanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indiwk, indisp;
    extern /* Subroutine */ int cstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static integer indrwk, indwrk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), ssterf_(integer *, doublereal *, doublereal *, integer *
	    );
    static integer nsplit, llwork;
    static doublereal smlnum;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
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

#line 388 "chbevx_2stage.f"
    /* Parameter adjustments */
#line 388 "chbevx_2stage.f"
    ab_dim1 = *ldab;
#line 388 "chbevx_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 388 "chbevx_2stage.f"
    ab -= ab_offset;
#line 388 "chbevx_2stage.f"
    q_dim1 = *ldq;
#line 388 "chbevx_2stage.f"
    q_offset = 1 + q_dim1;
#line 388 "chbevx_2stage.f"
    q -= q_offset;
#line 388 "chbevx_2stage.f"
    --w;
#line 388 "chbevx_2stage.f"
    z_dim1 = *ldz;
#line 388 "chbevx_2stage.f"
    z_offset = 1 + z_dim1;
#line 388 "chbevx_2stage.f"
    z__ -= z_offset;
#line 388 "chbevx_2stage.f"
    --work;
#line 388 "chbevx_2stage.f"
    --rwork;
#line 388 "chbevx_2stage.f"
    --iwork;
#line 388 "chbevx_2stage.f"
    --ifail;
#line 388 "chbevx_2stage.f"

#line 388 "chbevx_2stage.f"
    /* Function Body */
#line 388 "chbevx_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 389 "chbevx_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 390 "chbevx_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 391 "chbevx_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 392 "chbevx_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 393 "chbevx_2stage.f"
    lquery = *lwork == -1;

#line 395 "chbevx_2stage.f"
    *info = 0;
#line 396 "chbevx_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 397 "chbevx_2stage.f"
	*info = -1;
#line 398 "chbevx_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 399 "chbevx_2stage.f"
	*info = -2;
#line 400 "chbevx_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 401 "chbevx_2stage.f"
	*info = -3;
#line 402 "chbevx_2stage.f"
    } else if (*n < 0) {
#line 403 "chbevx_2stage.f"
	*info = -4;
#line 404 "chbevx_2stage.f"
    } else if (*kd < 0) {
#line 405 "chbevx_2stage.f"
	*info = -5;
#line 406 "chbevx_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 407 "chbevx_2stage.f"
	*info = -7;
#line 408 "chbevx_2stage.f"
    } else if (wantz && *ldq < max(1,*n)) {
#line 409 "chbevx_2stage.f"
	*info = -9;
#line 410 "chbevx_2stage.f"
    } else {
#line 411 "chbevx_2stage.f"
	if (valeig) {
#line 412 "chbevx_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 412 "chbevx_2stage.f"
		*info = -11;
#line 412 "chbevx_2stage.f"
	    }
#line 414 "chbevx_2stage.f"
	} else if (indeig) {
#line 415 "chbevx_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 416 "chbevx_2stage.f"
		*info = -12;
#line 417 "chbevx_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 418 "chbevx_2stage.f"
		*info = -13;
#line 419 "chbevx_2stage.f"
	    }
#line 420 "chbevx_2stage.f"
	}
#line 421 "chbevx_2stage.f"
    }
#line 422 "chbevx_2stage.f"
    if (*info == 0) {
#line 423 "chbevx_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 423 "chbevx_2stage.f"
	    *info = -18;
#line 423 "chbevx_2stage.f"
	}
#line 425 "chbevx_2stage.f"
    }

#line 427 "chbevx_2stage.f"
    if (*info == 0) {
#line 428 "chbevx_2stage.f"
	if (*n <= 1) {
#line 429 "chbevx_2stage.f"
	    lwmin = 1;
#line 430 "chbevx_2stage.f"
	    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 431 "chbevx_2stage.f"
	} else {
#line 432 "chbevx_2stage.f"
	    ib = ilaenv_(&c__18, "CHETRD_HB2ST", jobz, n, kd, &c_n1, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 433 "chbevx_2stage.f"
	    lhtrd = ilaenv_(&c__19, "CHETRD_HB2ST", jobz, n, kd, &ib, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 434 "chbevx_2stage.f"
	    lwtrd = ilaenv_(&c__20, "CHETRD_HB2ST", jobz, n, kd, &ib, &c_n1, (
		    ftnlen)12, (ftnlen)1);
#line 435 "chbevx_2stage.f"
	    lwmin = lhtrd + lwtrd;
#line 436 "chbevx_2stage.f"
	    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 437 "chbevx_2stage.f"
	}

#line 439 "chbevx_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 439 "chbevx_2stage.f"
	    *info = -20;
#line 439 "chbevx_2stage.f"
	}
#line 441 "chbevx_2stage.f"
    }

#line 443 "chbevx_2stage.f"
    if (*info != 0) {
#line 444 "chbevx_2stage.f"
	i__1 = -(*info);
#line 444 "chbevx_2stage.f"
	xerbla_("CHBEVX_2STAGE", &i__1, (ftnlen)13);
#line 445 "chbevx_2stage.f"
	return 0;
#line 446 "chbevx_2stage.f"
    } else if (lquery) {
#line 447 "chbevx_2stage.f"
	return 0;
#line 448 "chbevx_2stage.f"
    }

/*     Quick return if possible */

#line 452 "chbevx_2stage.f"
    *m = 0;
#line 453 "chbevx_2stage.f"
    if (*n == 0) {
#line 453 "chbevx_2stage.f"
	return 0;
#line 453 "chbevx_2stage.f"
    }

#line 456 "chbevx_2stage.f"
    if (*n == 1) {
#line 457 "chbevx_2stage.f"
	*m = 1;
#line 458 "chbevx_2stage.f"
	if (lower) {
#line 459 "chbevx_2stage.f"
	    i__1 = ab_dim1 + 1;
#line 459 "chbevx_2stage.f"
	    ctmp1.r = ab[i__1].r, ctmp1.i = ab[i__1].i;
#line 460 "chbevx_2stage.f"
	} else {
#line 461 "chbevx_2stage.f"
	    i__1 = *kd + 1 + ab_dim1;
#line 461 "chbevx_2stage.f"
	    ctmp1.r = ab[i__1].r, ctmp1.i = ab[i__1].i;
#line 462 "chbevx_2stage.f"
	}
#line 463 "chbevx_2stage.f"
	tmp1 = ctmp1.r;
#line 464 "chbevx_2stage.f"
	if (valeig) {
#line 465 "chbevx_2stage.f"
	    if (! (*vl < tmp1 && *vu >= tmp1)) {
#line 465 "chbevx_2stage.f"
		*m = 0;
#line 465 "chbevx_2stage.f"
	    }
#line 467 "chbevx_2stage.f"
	}
#line 468 "chbevx_2stage.f"
	if (*m == 1) {
#line 469 "chbevx_2stage.f"
	    w[1] = ctmp1.r;
#line 470 "chbevx_2stage.f"
	    if (wantz) {
#line 470 "chbevx_2stage.f"
		i__1 = z_dim1 + 1;
#line 470 "chbevx_2stage.f"
		z__[i__1].r = 1., z__[i__1].i = 0.;
#line 470 "chbevx_2stage.f"
	    }
#line 472 "chbevx_2stage.f"
	}
#line 473 "chbevx_2stage.f"
	return 0;
#line 474 "chbevx_2stage.f"
    }

/*     Get machine constants. */

#line 478 "chbevx_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 479 "chbevx_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 480 "chbevx_2stage.f"
    smlnum = safmin / eps;
#line 481 "chbevx_2stage.f"
    bignum = 1. / smlnum;
#line 482 "chbevx_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 483 "chbevx_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 483 "chbevx_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 487 "chbevx_2stage.f"
    iscale = 0;
#line 488 "chbevx_2stage.f"
    abstll = *abstol;
#line 489 "chbevx_2stage.f"
    if (valeig) {
#line 490 "chbevx_2stage.f"
	vll = *vl;
#line 491 "chbevx_2stage.f"
	vuu = *vu;
#line 492 "chbevx_2stage.f"
    } else {
#line 493 "chbevx_2stage.f"
	vll = 0.;
#line 494 "chbevx_2stage.f"
	vuu = 0.;
#line 495 "chbevx_2stage.f"
    }
#line 496 "chbevx_2stage.f"
    anrm = clanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1, (ftnlen)1);
#line 497 "chbevx_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 498 "chbevx_2stage.f"
	iscale = 1;
#line 499 "chbevx_2stage.f"
	sigma = rmin / anrm;
#line 500 "chbevx_2stage.f"
    } else if (anrm > rmax) {
#line 501 "chbevx_2stage.f"
	iscale = 1;
#line 502 "chbevx_2stage.f"
	sigma = rmax / anrm;
#line 503 "chbevx_2stage.f"
    }
#line 504 "chbevx_2stage.f"
    if (iscale == 1) {
#line 505 "chbevx_2stage.f"
	if (lower) {
#line 506 "chbevx_2stage.f"
	    clascl_("B", kd, kd, &c_b26, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 507 "chbevx_2stage.f"
	} else {
#line 508 "chbevx_2stage.f"
	    clascl_("Q", kd, kd, &c_b26, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 509 "chbevx_2stage.f"
	}
#line 510 "chbevx_2stage.f"
	if (*abstol > 0.) {
#line 510 "chbevx_2stage.f"
	    abstll = *abstol * sigma;
#line 510 "chbevx_2stage.f"
	}
#line 512 "chbevx_2stage.f"
	if (valeig) {
#line 513 "chbevx_2stage.f"
	    vll = *vl * sigma;
#line 514 "chbevx_2stage.f"
	    vuu = *vu * sigma;
#line 515 "chbevx_2stage.f"
	}
#line 516 "chbevx_2stage.f"
    }

/*     Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form. */

#line 520 "chbevx_2stage.f"
    indd = 1;
#line 521 "chbevx_2stage.f"
    inde = indd + *n;
#line 522 "chbevx_2stage.f"
    indrwk = inde + *n;

#line 524 "chbevx_2stage.f"
    indhous = 1;
#line 525 "chbevx_2stage.f"
    indwrk = indhous + lhtrd;
#line 526 "chbevx_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 528 "chbevx_2stage.f"
    chetrd_hb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &rwork[indd],
	     &rwork[inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call SSTERF or CSTEQR.  If this fails for some */
/*     eigenvalue, then try SSTEBZ. */

#line 536 "chbevx_2stage.f"
    test = FALSE_;
#line 537 "chbevx_2stage.f"
    if (indeig) {
#line 538 "chbevx_2stage.f"
	if (*il == 1 && *iu == *n) {
#line 539 "chbevx_2stage.f"
	    test = TRUE_;
#line 540 "chbevx_2stage.f"
	}
#line 541 "chbevx_2stage.f"
    }
#line 542 "chbevx_2stage.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 543 "chbevx_2stage.f"
	scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
#line 544 "chbevx_2stage.f"
	indee = indrwk + (*n << 1);
#line 545 "chbevx_2stage.f"
	if (! wantz) {
#line 546 "chbevx_2stage.f"
	    i__1 = *n - 1;
#line 546 "chbevx_2stage.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 547 "chbevx_2stage.f"
	    ssterf_(n, &w[1], &rwork[indee], info);
#line 548 "chbevx_2stage.f"
	} else {
#line 549 "chbevx_2stage.f"
	    clacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 550 "chbevx_2stage.f"
	    i__1 = *n - 1;
#line 550 "chbevx_2stage.f"
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
#line 551 "chbevx_2stage.f"
	    csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info, (ftnlen)1);
#line 553 "chbevx_2stage.f"
	    if (*info == 0) {
#line 554 "chbevx_2stage.f"
		i__1 = *n;
#line 554 "chbevx_2stage.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 555 "chbevx_2stage.f"
		    ifail[i__] = 0;
#line 556 "chbevx_2stage.f"
/* L10: */
#line 556 "chbevx_2stage.f"
		}
#line 557 "chbevx_2stage.f"
	    }
#line 558 "chbevx_2stage.f"
	}
#line 559 "chbevx_2stage.f"
	if (*info == 0) {
#line 560 "chbevx_2stage.f"
	    *m = *n;
#line 561 "chbevx_2stage.f"
	    goto L30;
#line 562 "chbevx_2stage.f"
	}
#line 563 "chbevx_2stage.f"
	*info = 0;
#line 564 "chbevx_2stage.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */

#line 568 "chbevx_2stage.f"
    if (wantz) {
#line 569 "chbevx_2stage.f"
	*(unsigned char *)order = 'B';
#line 570 "chbevx_2stage.f"
    } else {
#line 571 "chbevx_2stage.f"
	*(unsigned char *)order = 'E';
#line 572 "chbevx_2stage.f"
    }
#line 573 "chbevx_2stage.f"
    indibl = 1;
#line 574 "chbevx_2stage.f"
    indisp = indibl + *n;
#line 575 "chbevx_2stage.f"
    indiwk = indisp + *n;
#line 576 "chbevx_2stage.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], info, (ftnlen)1, (ftnlen)1);

#line 581 "chbevx_2stage.f"
    if (wantz) {
#line 582 "chbevx_2stage.f"
	cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by CSTEIN. */

#line 589 "chbevx_2stage.f"
	i__1 = *m;
#line 589 "chbevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 590 "chbevx_2stage.f"
	    ccopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 591 "chbevx_2stage.f"
	    cgemv_("N", n, n, &c_b2, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b1, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 593 "chbevx_2stage.f"
/* L20: */
#line 593 "chbevx_2stage.f"
	}
#line 594 "chbevx_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 598 "chbevx_2stage.f"
L30:
#line 599 "chbevx_2stage.f"
    if (iscale == 1) {
#line 600 "chbevx_2stage.f"
	if (*info == 0) {
#line 601 "chbevx_2stage.f"
	    imax = *m;
#line 602 "chbevx_2stage.f"
	} else {
#line 603 "chbevx_2stage.f"
	    imax = *info - 1;
#line 604 "chbevx_2stage.f"
	}
#line 605 "chbevx_2stage.f"
	d__1 = 1. / sigma;
#line 605 "chbevx_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 606 "chbevx_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 611 "chbevx_2stage.f"
    if (wantz) {
#line 612 "chbevx_2stage.f"
	i__1 = *m - 1;
#line 612 "chbevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 613 "chbevx_2stage.f"
	    i__ = 0;
#line 614 "chbevx_2stage.f"
	    tmp1 = w[j];
#line 615 "chbevx_2stage.f"
	    i__2 = *m;
#line 615 "chbevx_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 616 "chbevx_2stage.f"
		if (w[jj] < tmp1) {
#line 617 "chbevx_2stage.f"
		    i__ = jj;
#line 618 "chbevx_2stage.f"
		    tmp1 = w[jj];
#line 619 "chbevx_2stage.f"
		}
#line 620 "chbevx_2stage.f"
/* L40: */
#line 620 "chbevx_2stage.f"
	    }

#line 622 "chbevx_2stage.f"
	    if (i__ != 0) {
#line 623 "chbevx_2stage.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 624 "chbevx_2stage.f"
		w[i__] = w[j];
#line 625 "chbevx_2stage.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 626 "chbevx_2stage.f"
		w[j] = tmp1;
#line 627 "chbevx_2stage.f"
		iwork[indibl + j - 1] = itmp1;
#line 628 "chbevx_2stage.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 629 "chbevx_2stage.f"
		if (*info != 0) {
#line 630 "chbevx_2stage.f"
		    itmp1 = ifail[i__];
#line 631 "chbevx_2stage.f"
		    ifail[i__] = ifail[j];
#line 632 "chbevx_2stage.f"
		    ifail[j] = itmp1;
#line 633 "chbevx_2stage.f"
		}
#line 634 "chbevx_2stage.f"
	    }
#line 635 "chbevx_2stage.f"
/* L50: */
#line 635 "chbevx_2stage.f"
	}
#line 636 "chbevx_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 640 "chbevx_2stage.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 642 "chbevx_2stage.f"
    return 0;

/*     End of CHBEVX_2STAGE */

} /* chbevx_2stage__ */

