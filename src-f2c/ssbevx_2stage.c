#line 1 "ssbevx_2stage.f"
/* ssbevx_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "ssbevx_2stage.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__4 = 4;
static doublereal c_b24 = 1.;
static integer c__1 = 1;
static doublereal c_b45 = 0.;

/* > \brief <b> SSBEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 OTHER matrices</b> */

/*  @generated from dsbevx_2stage.f, fortran d -> s, Sat Nov  5 23:58:06 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSBEVX_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbevx_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbevx_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbevx_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, */
/*                                 LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, */
/*                                 LDZ, WORK, LWORK, IWORK, IFAIL, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IFAIL( * ), IWORK( * ) */
/*       REAL               AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
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
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is REAL array, dimension (LDQ, N) */
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
/* >          WORK is REAL array, dimension (LWORK) */
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
/* Subroutine */ int ssbevx_2stage__(char *jobz, char *range, char *uplo, 
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
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    extern /* Subroutine */ int ssytrd_sb2st__(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer itmp1, indee;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    static integer lhtrd;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer lwmin;
    static logical lower;
    static integer lwtrd;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical wantz, alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    extern doublereal slansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer indisp, indiwo;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int sstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    ssterf_(integer *, doublereal *, doublereal *, integer *);
    static integer nsplit, llwork;
    static doublereal smlnum;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer indhous;



/*  -- LAPACK driver routine (version 3.8.0) -- */
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

#line 378 "ssbevx_2stage.f"
    /* Parameter adjustments */
#line 378 "ssbevx_2stage.f"
    ab_dim1 = *ldab;
#line 378 "ssbevx_2stage.f"
    ab_offset = 1 + ab_dim1;
#line 378 "ssbevx_2stage.f"
    ab -= ab_offset;
#line 378 "ssbevx_2stage.f"
    q_dim1 = *ldq;
#line 378 "ssbevx_2stage.f"
    q_offset = 1 + q_dim1;
#line 378 "ssbevx_2stage.f"
    q -= q_offset;
#line 378 "ssbevx_2stage.f"
    --w;
#line 378 "ssbevx_2stage.f"
    z_dim1 = *ldz;
#line 378 "ssbevx_2stage.f"
    z_offset = 1 + z_dim1;
#line 378 "ssbevx_2stage.f"
    z__ -= z_offset;
#line 378 "ssbevx_2stage.f"
    --work;
#line 378 "ssbevx_2stage.f"
    --iwork;
#line 378 "ssbevx_2stage.f"
    --ifail;
#line 378 "ssbevx_2stage.f"

#line 378 "ssbevx_2stage.f"
    /* Function Body */
#line 378 "ssbevx_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 379 "ssbevx_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 380 "ssbevx_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 381 "ssbevx_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
#line 382 "ssbevx_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 383 "ssbevx_2stage.f"
    lquery = *lwork == -1;

#line 385 "ssbevx_2stage.f"
    *info = 0;
#line 386 "ssbevx_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 387 "ssbevx_2stage.f"
	*info = -1;
#line 388 "ssbevx_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 389 "ssbevx_2stage.f"
	*info = -2;
#line 390 "ssbevx_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 391 "ssbevx_2stage.f"
	*info = -3;
#line 392 "ssbevx_2stage.f"
    } else if (*n < 0) {
#line 393 "ssbevx_2stage.f"
	*info = -4;
#line 394 "ssbevx_2stage.f"
    } else if (*kd < 0) {
#line 395 "ssbevx_2stage.f"
	*info = -5;
#line 396 "ssbevx_2stage.f"
    } else if (*ldab < *kd + 1) {
#line 397 "ssbevx_2stage.f"
	*info = -7;
#line 398 "ssbevx_2stage.f"
    } else if (wantz && *ldq < max(1,*n)) {
#line 399 "ssbevx_2stage.f"
	*info = -9;
#line 400 "ssbevx_2stage.f"
    } else {
#line 401 "ssbevx_2stage.f"
	if (valeig) {
#line 402 "ssbevx_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 402 "ssbevx_2stage.f"
		*info = -11;
#line 402 "ssbevx_2stage.f"
	    }
#line 404 "ssbevx_2stage.f"
	} else if (indeig) {
#line 405 "ssbevx_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 406 "ssbevx_2stage.f"
		*info = -12;
#line 407 "ssbevx_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 408 "ssbevx_2stage.f"
		*info = -13;
#line 409 "ssbevx_2stage.f"
	    }
#line 410 "ssbevx_2stage.f"
	}
#line 411 "ssbevx_2stage.f"
    }
#line 412 "ssbevx_2stage.f"
    if (*info == 0) {
#line 413 "ssbevx_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 413 "ssbevx_2stage.f"
	    *info = -18;
#line 413 "ssbevx_2stage.f"
	}
#line 415 "ssbevx_2stage.f"
    }

#line 417 "ssbevx_2stage.f"
    if (*info == 0) {
#line 418 "ssbevx_2stage.f"
	if (*n <= 1) {
#line 419 "ssbevx_2stage.f"
	    lwmin = 1;
#line 420 "ssbevx_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 421 "ssbevx_2stage.f"
	} else {
#line 422 "ssbevx_2stage.f"
	    ib = ilaenv2stage_(&c__2, "SSYTRD_SB2ST", jobz, n, kd, &c_n1, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 424 "ssbevx_2stage.f"
	    lhtrd = ilaenv2stage_(&c__3, "SSYTRD_SB2ST", jobz, n, kd, &ib, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 426 "ssbevx_2stage.f"
	    lwtrd = ilaenv2stage_(&c__4, "SSYTRD_SB2ST", jobz, n, kd, &ib, &
		    c_n1, (ftnlen)12, (ftnlen)1);
#line 428 "ssbevx_2stage.f"
	    lwmin = (*n << 1) + lhtrd + lwtrd;
#line 429 "ssbevx_2stage.f"
	    work[1] = (doublereal) lwmin;
#line 430 "ssbevx_2stage.f"
	}

#line 432 "ssbevx_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 432 "ssbevx_2stage.f"
	    *info = -20;
#line 432 "ssbevx_2stage.f"
	}
#line 434 "ssbevx_2stage.f"
    }

#line 436 "ssbevx_2stage.f"
    if (*info != 0) {
#line 437 "ssbevx_2stage.f"
	i__1 = -(*info);
#line 437 "ssbevx_2stage.f"
	xerbla_("SSBEVX_2STAGE ", &i__1, (ftnlen)14);
#line 438 "ssbevx_2stage.f"
	return 0;
#line 439 "ssbevx_2stage.f"
    } else if (lquery) {
#line 440 "ssbevx_2stage.f"
	return 0;
#line 441 "ssbevx_2stage.f"
    }

/*     Quick return if possible */

#line 445 "ssbevx_2stage.f"
    *m = 0;
#line 446 "ssbevx_2stage.f"
    if (*n == 0) {
#line 446 "ssbevx_2stage.f"
	return 0;
#line 446 "ssbevx_2stage.f"
    }

#line 449 "ssbevx_2stage.f"
    if (*n == 1) {
#line 450 "ssbevx_2stage.f"
	*m = 1;
#line 451 "ssbevx_2stage.f"
	if (lower) {
#line 452 "ssbevx_2stage.f"
	    tmp1 = ab[ab_dim1 + 1];
#line 453 "ssbevx_2stage.f"
	} else {
#line 454 "ssbevx_2stage.f"
	    tmp1 = ab[*kd + 1 + ab_dim1];
#line 455 "ssbevx_2stage.f"
	}
#line 456 "ssbevx_2stage.f"
	if (valeig) {
#line 457 "ssbevx_2stage.f"
	    if (! (*vl < tmp1 && *vu >= tmp1)) {
#line 457 "ssbevx_2stage.f"
		*m = 0;
#line 457 "ssbevx_2stage.f"
	    }
#line 459 "ssbevx_2stage.f"
	}
#line 460 "ssbevx_2stage.f"
	if (*m == 1) {
#line 461 "ssbevx_2stage.f"
	    w[1] = tmp1;
#line 462 "ssbevx_2stage.f"
	    if (wantz) {
#line 462 "ssbevx_2stage.f"
		z__[z_dim1 + 1] = 1.;
#line 462 "ssbevx_2stage.f"
	    }
#line 464 "ssbevx_2stage.f"
	}
#line 465 "ssbevx_2stage.f"
	return 0;
#line 466 "ssbevx_2stage.f"
    }

/*     Get machine constants. */

#line 470 "ssbevx_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 471 "ssbevx_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 472 "ssbevx_2stage.f"
    smlnum = safmin / eps;
#line 473 "ssbevx_2stage.f"
    bignum = 1. / smlnum;
#line 474 "ssbevx_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 475 "ssbevx_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 475 "ssbevx_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 479 "ssbevx_2stage.f"
    iscale = 0;
#line 480 "ssbevx_2stage.f"
    abstll = *abstol;
#line 481 "ssbevx_2stage.f"
    if (valeig) {
#line 482 "ssbevx_2stage.f"
	vll = *vl;
#line 483 "ssbevx_2stage.f"
	vuu = *vu;
#line 484 "ssbevx_2stage.f"
    } else {
#line 485 "ssbevx_2stage.f"
	vll = 0.;
#line 486 "ssbevx_2stage.f"
	vuu = 0.;
#line 487 "ssbevx_2stage.f"
    }
#line 488 "ssbevx_2stage.f"
    anrm = slansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 489 "ssbevx_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 490 "ssbevx_2stage.f"
	iscale = 1;
#line 491 "ssbevx_2stage.f"
	sigma = rmin / anrm;
#line 492 "ssbevx_2stage.f"
    } else if (anrm > rmax) {
#line 493 "ssbevx_2stage.f"
	iscale = 1;
#line 494 "ssbevx_2stage.f"
	sigma = rmax / anrm;
#line 495 "ssbevx_2stage.f"
    }
#line 496 "ssbevx_2stage.f"
    if (iscale == 1) {
#line 497 "ssbevx_2stage.f"
	if (lower) {
#line 498 "ssbevx_2stage.f"
	    slascl_("B", kd, kd, &c_b24, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 499 "ssbevx_2stage.f"
	} else {
#line 500 "ssbevx_2stage.f"
	    slascl_("Q", kd, kd, &c_b24, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 501 "ssbevx_2stage.f"
	}
#line 502 "ssbevx_2stage.f"
	if (*abstol > 0.) {
#line 502 "ssbevx_2stage.f"
	    abstll = *abstol * sigma;
#line 502 "ssbevx_2stage.f"
	}
#line 504 "ssbevx_2stage.f"
	if (valeig) {
#line 505 "ssbevx_2stage.f"
	    vll = *vl * sigma;
#line 506 "ssbevx_2stage.f"
	    vuu = *vu * sigma;
#line 507 "ssbevx_2stage.f"
	}
#line 508 "ssbevx_2stage.f"
    }

/*     Call SSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form. */

#line 512 "ssbevx_2stage.f"
    indd = 1;
#line 513 "ssbevx_2stage.f"
    inde = indd + *n;
#line 514 "ssbevx_2stage.f"
    indhous = inde + *n;
#line 515 "ssbevx_2stage.f"
    indwrk = indhous + lhtrd;
#line 516 "ssbevx_2stage.f"
    llwork = *lwork - indwrk + 1;

#line 518 "ssbevx_2stage.f"
    ssytrd_sb2st__("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &work[indd], 
	    &work[inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired and ABSTOL is less than or equal */
/*     to zero, then call SSTERF or SSTEQR.  If this fails for some */
/*     eigenvalue, then try SSTEBZ. */

#line 526 "ssbevx_2stage.f"
    test = FALSE_;
#line 527 "ssbevx_2stage.f"
    if (indeig) {
#line 528 "ssbevx_2stage.f"
	if (*il == 1 && *iu == *n) {
#line 529 "ssbevx_2stage.f"
	    test = TRUE_;
#line 530 "ssbevx_2stage.f"
	}
#line 531 "ssbevx_2stage.f"
    }
#line 532 "ssbevx_2stage.f"
    if ((alleig || test) && *abstol <= 0.) {
#line 533 "ssbevx_2stage.f"
	scopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 534 "ssbevx_2stage.f"
	indee = indwrk + (*n << 1);
#line 535 "ssbevx_2stage.f"
	if (! wantz) {
#line 536 "ssbevx_2stage.f"
	    i__1 = *n - 1;
#line 536 "ssbevx_2stage.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 537 "ssbevx_2stage.f"
	    ssterf_(n, &w[1], &work[indee], info);
#line 538 "ssbevx_2stage.f"
	} else {
#line 539 "ssbevx_2stage.f"
	    slacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz, (
		    ftnlen)1);
#line 540 "ssbevx_2stage.f"
	    i__1 = *n - 1;
#line 540 "ssbevx_2stage.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 541 "ssbevx_2stage.f"
	    ssteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
#line 543 "ssbevx_2stage.f"
	    if (*info == 0) {
#line 544 "ssbevx_2stage.f"
		i__1 = *n;
#line 544 "ssbevx_2stage.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 545 "ssbevx_2stage.f"
		    ifail[i__] = 0;
#line 546 "ssbevx_2stage.f"
/* L10: */
#line 546 "ssbevx_2stage.f"
		}
#line 547 "ssbevx_2stage.f"
	    }
#line 548 "ssbevx_2stage.f"
	}
#line 549 "ssbevx_2stage.f"
	if (*info == 0) {
#line 550 "ssbevx_2stage.f"
	    *m = *n;
#line 551 "ssbevx_2stage.f"
	    goto L30;
#line 552 "ssbevx_2stage.f"
	}
#line 553 "ssbevx_2stage.f"
	*info = 0;
#line 554 "ssbevx_2stage.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 558 "ssbevx_2stage.f"
    if (wantz) {
#line 559 "ssbevx_2stage.f"
	*(unsigned char *)order = 'B';
#line 560 "ssbevx_2stage.f"
    } else {
#line 561 "ssbevx_2stage.f"
	*(unsigned char *)order = 'E';
#line 562 "ssbevx_2stage.f"
    }
#line 563 "ssbevx_2stage.f"
    indibl = 1;
#line 564 "ssbevx_2stage.f"
    indisp = indibl + *n;
#line 565 "ssbevx_2stage.f"
    indiwo = indisp + *n;
#line 566 "ssbevx_2stage.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 571 "ssbevx_2stage.f"
    if (wantz) {
#line 572 "ssbevx_2stage.f"
	sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEIN. */

#line 579 "ssbevx_2stage.f"
	i__1 = *m;
#line 579 "ssbevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 580 "ssbevx_2stage.f"
	    scopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
#line 581 "ssbevx_2stage.f"
	    sgemv_("N", n, n, &c_b24, &q[q_offset], ldq, &work[1], &c__1, &
		    c_b45, &z__[j * z_dim1 + 1], &c__1, (ftnlen)1);
#line 583 "ssbevx_2stage.f"
/* L20: */
#line 583 "ssbevx_2stage.f"
	}
#line 584 "ssbevx_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 588 "ssbevx_2stage.f"
L30:
#line 589 "ssbevx_2stage.f"
    if (iscale == 1) {
#line 590 "ssbevx_2stage.f"
	if (*info == 0) {
#line 591 "ssbevx_2stage.f"
	    imax = *m;
#line 592 "ssbevx_2stage.f"
	} else {
#line 593 "ssbevx_2stage.f"
	    imax = *info - 1;
#line 594 "ssbevx_2stage.f"
	}
#line 595 "ssbevx_2stage.f"
	d__1 = 1. / sigma;
#line 595 "ssbevx_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 596 "ssbevx_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 601 "ssbevx_2stage.f"
    if (wantz) {
#line 602 "ssbevx_2stage.f"
	i__1 = *m - 1;
#line 602 "ssbevx_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 603 "ssbevx_2stage.f"
	    i__ = 0;
#line 604 "ssbevx_2stage.f"
	    tmp1 = w[j];
#line 605 "ssbevx_2stage.f"
	    i__2 = *m;
#line 605 "ssbevx_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 606 "ssbevx_2stage.f"
		if (w[jj] < tmp1) {
#line 607 "ssbevx_2stage.f"
		    i__ = jj;
#line 608 "ssbevx_2stage.f"
		    tmp1 = w[jj];
#line 609 "ssbevx_2stage.f"
		}
#line 610 "ssbevx_2stage.f"
/* L40: */
#line 610 "ssbevx_2stage.f"
	    }

#line 612 "ssbevx_2stage.f"
	    if (i__ != 0) {
#line 613 "ssbevx_2stage.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 614 "ssbevx_2stage.f"
		w[i__] = w[j];
#line 615 "ssbevx_2stage.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 616 "ssbevx_2stage.f"
		w[j] = tmp1;
#line 617 "ssbevx_2stage.f"
		iwork[indibl + j - 1] = itmp1;
#line 618 "ssbevx_2stage.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 619 "ssbevx_2stage.f"
		if (*info != 0) {
#line 620 "ssbevx_2stage.f"
		    itmp1 = ifail[i__];
#line 621 "ssbevx_2stage.f"
		    ifail[i__] = ifail[j];
#line 622 "ssbevx_2stage.f"
		    ifail[j] = itmp1;
#line 623 "ssbevx_2stage.f"
		}
#line 624 "ssbevx_2stage.f"
	    }
#line 625 "ssbevx_2stage.f"
/* L50: */
#line 625 "ssbevx_2stage.f"
	}
#line 626 "ssbevx_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 630 "ssbevx_2stage.f"
    work[1] = (doublereal) lwmin;

#line 632 "ssbevx_2stage.f"
    return 0;

/*     End of SSBEVX_2STAGE */

} /* ssbevx_2stage__ */

