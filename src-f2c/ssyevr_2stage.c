#line 1 "ssyevr_2stage.f"
/* ssyevr_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "ssyevr_2stage.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;

/* > \brief <b> SSYEVR_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 SY matrices</b> */

/*  @generated from dsyevr_2stage.f, fortran d -> s, Sat Nov  5 23:50:10 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYEVR_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevr_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevr_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevr_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, */
/*                          IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, */
/*                          LWORK, IWORK, LIWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEVR_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  Eigenvalues and eigenvectors can be */
/* > selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
/* > */
/* > SSYEVR_2STAGE first reduces the matrix A to tridiagonal form T with a call */
/* > to SSYTRD.  Then, whenever possible, SSYEVR_2STAGE calls SSTEMR to compute */
/* > the eigenspectrum using Relatively Robust Representations.  SSTEMR */
/* > computes eigenvalues by the dqds algorithm, while orthogonal */
/* > eigenvectors are computed from various "good" L D L^T representations */
/* > (also known as Relatively Robust Representations). Gram-Schmidt */
/* > orthogonalization is avoided as far as possible. More specifically, */
/* > the various steps of the algorithm are as follows. */
/* > */
/* > For each unreduced block (submatrix) of T, */
/* >    (a) Compute T - sigma I  = L D L^T, so that L and D */
/* >        define all the wanted eigenvalues to high relative accuracy. */
/* >        This means that small relative changes in the entries of D and L */
/* >        cause only small relative changes in the eigenvalues and */
/* >        eigenvectors. The standard (unfactored) representation of the */
/* >        tridiagonal matrix T does not have this property in general. */
/* >    (b) Compute the eigenvalues to suitable accuracy. */
/* >        If the eigenvectors are desired, the algorithm attains full */
/* >        accuracy of the computed eigenvalues only right before */
/* >        the corresponding vectors have to be computed, see steps c) and d). */
/* >    (c) For each cluster of close eigenvalues, select a new */
/* >        shift close to the cluster, find a new factorization, and refine */
/* >        the shifted eigenvalues to suitable accuracy. */
/* >    (d) For each eigenvalue with a large enough relative separation compute */
/* >        the corresponding eigenvector by forming a rank revealing twisted */
/* >        factorization. Go back to (c) for any clusters that remain. */
/* > */
/* > The desired accuracy of the output can be specified by the input */
/* > parameter ABSTOL. */
/* > */
/* > For more details, see SSTEMR's documentation and: */
/* > - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations */
/* >   to compute orthogonal eigenvectors of symmetric tridiagonal matrices," */
/* >   Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004. */
/* > - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and */
/* >   Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25, */
/* >   2004.  Also LAPACK Working Note 154. */
/* > - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric */
/* >   tridiagonal eigenvalue/eigenvector problem", */
/* >   Computer Science Division Technical Report No. UCB/CSD-97-971, */
/* >   UC Berkeley, May 1997. */
/* > */
/* > */
/* > Note 1 : SSYEVR_2STAGE calls SSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > SSYEVR_2STAGE calls SSTEBZ and SSTEIN on non-ieee machines and */
/* > when partial spectrum requests are made. */
/* > */
/* > Normal execution of SSTEMR may create NaNs and infinities and */
/* > hence may abort due to a floating point exception in environments */
/* > which do not handle NaNs and infinities in the ieee standard default */
/* > manner. */
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
/* >          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and */
/* >          SSTEIN are called */
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
/* >          See "Computing Small Singular Values of Bidiagonal Matrices */
/* >          with Guaranteed High Relative Accuracy," by Demmel and */
/* >          Kahan, LAPACK Working Note #3. */
/* > */
/* >          If high relative accuracy is important, set ABSTOL to */
/* >          SLAMCH( 'Safe minimum' ).  Doing so will guarantee that */
/* >          eigenvalues are computed to high relative accuracy when */
/* >          possible in future releases.  The current code does not */
/* >          make any guarantees about high relative accuracy, but */
/* >          future releases will. See J. Barlow and J. Demmel, */
/* >          "Computing Accurate Eigensystems of Scaled Diagonally */
/* >          Dominant Matrices", LAPACK Working Note #7, for a discussion */
/* >          of which matrices define their eigenvalues to high relative */
/* >          accuracy. */
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
/* >          If JOBZ = 'N', then Z is not referenced. */
/* >          Note: the user must ensure that at least max(1,M) columns are */
/* >          supplied in the array Z; if RANGE = 'V', the exact value of M */
/* >          is not known in advance and an upper bound must be used. */
/* >          Supplying N columns is always safe. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* >          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) ) */
/* >          The support of the eigenvectors in Z, i.e., the indices */
/* >          indicating the nonzero elements in Z. The i-th eigenvector */
/* >          is nonzero only in elements ISUPPZ( 2*i-1 ) through */
/* >          ISUPPZ( 2*i ). This is an output of SSTEMR (tridiagonal */
/* >          matrix). The support of the eigenvectors of A is typically */
/* >          1:N because of the orthogonal transformations applied by SORMTR. */
/* >          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 */
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
/* >          The dimension of the array WORK. */
/* >          If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* >                                   LWORK = MAX(1, 26*N, dimension) where */
/* >                                   dimension = max(stage1,stage2) + (KD+1)*N + 5*N */
/* >                                             = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                                               + max(2*KD*KD, KD*NTHREADS) */
/* >                                               + (KD+1)*N + 5*N */
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
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK.  LIWORK >= max(1,10*N). */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of the IWORK array, */
/* >          returns this value as the first entry of the IWORK array, and */
/* >          no error message related to LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  Internal error */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup realSYeigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Inderjit Dhillon, IBM Almaden, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */
/* >     Ken Stanley, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Jason Riedy, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* > */
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
/* Subroutine */ int ssyevr_2stage__(char *jobz, char *range, char *uplo, 
	integer *n, doublereal *a, integer *lda, doublereal *vl, doublereal *
	vu, integer *il, integer *iu, doublereal *abstol, integer *m, 
	doublereal *w, doublereal *z__, integer *ldz, integer *isuppz, 
	doublereal *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info, ftnlen jobz_len, ftnlen range_len, ftnlen uplo_len)
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
    static integer inddd, indee;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    static integer indwk, lhtrd, lwmin;
    static logical lower;
    static integer lwtrd;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), ssytrd_2stage__(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static logical wantz, alleig, indeig;
    static integer iscale, ieeeok, indibl, indifl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp, indiwo, indwkn, liwmin;
    static logical tryrac;
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
	    sstemr_(char *, char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int sormtr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
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

#line 436 "ssyevr_2stage.f"
    /* Parameter adjustments */
#line 436 "ssyevr_2stage.f"
    a_dim1 = *lda;
#line 436 "ssyevr_2stage.f"
    a_offset = 1 + a_dim1;
#line 436 "ssyevr_2stage.f"
    a -= a_offset;
#line 436 "ssyevr_2stage.f"
    --w;
#line 436 "ssyevr_2stage.f"
    z_dim1 = *ldz;
#line 436 "ssyevr_2stage.f"
    z_offset = 1 + z_dim1;
#line 436 "ssyevr_2stage.f"
    z__ -= z_offset;
#line 436 "ssyevr_2stage.f"
    --isuppz;
#line 436 "ssyevr_2stage.f"
    --work;
#line 436 "ssyevr_2stage.f"
    --iwork;
#line 436 "ssyevr_2stage.f"

#line 436 "ssyevr_2stage.f"
    /* Function Body */
#line 436 "ssyevr_2stage.f"
    ieeeok = ilaenv_(&c__10, "SSYEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 438 "ssyevr_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 439 "ssyevr_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 440 "ssyevr_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 441 "ssyevr_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 442 "ssyevr_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 444 "ssyevr_2stage.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 446 "ssyevr_2stage.f"
    kd = ilaenv2stage_(&c__1, "SSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 447 "ssyevr_2stage.f"
    ib = ilaenv2stage_(&c__2, "SSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 448 "ssyevr_2stage.f"
    lhtrd = ilaenv2stage_(&c__3, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 449 "ssyevr_2stage.f"
    lwtrd = ilaenv2stage_(&c__4, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
/* Computing MAX */
#line 450 "ssyevr_2stage.f"
    i__1 = *n * 26, i__2 = *n * 5 + lhtrd + lwtrd;
#line 450 "ssyevr_2stage.f"
    lwmin = max(i__1,i__2);
/* Computing MAX */
#line 451 "ssyevr_2stage.f"
    i__1 = 1, i__2 = *n * 10;
#line 451 "ssyevr_2stage.f"
    liwmin = max(i__1,i__2);

#line 453 "ssyevr_2stage.f"
    *info = 0;
#line 454 "ssyevr_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 455 "ssyevr_2stage.f"
	*info = -1;
#line 456 "ssyevr_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 457 "ssyevr_2stage.f"
	*info = -2;
#line 458 "ssyevr_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 459 "ssyevr_2stage.f"
	*info = -3;
#line 460 "ssyevr_2stage.f"
    } else if (*n < 0) {
#line 461 "ssyevr_2stage.f"
	*info = -4;
#line 462 "ssyevr_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 463 "ssyevr_2stage.f"
	*info = -6;
#line 464 "ssyevr_2stage.f"
    } else {
#line 465 "ssyevr_2stage.f"
	if (valeig) {
#line 466 "ssyevr_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 466 "ssyevr_2stage.f"
		*info = -8;
#line 466 "ssyevr_2stage.f"
	    }
#line 468 "ssyevr_2stage.f"
	} else if (indeig) {
#line 469 "ssyevr_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 470 "ssyevr_2stage.f"
		*info = -9;
#line 471 "ssyevr_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 472 "ssyevr_2stage.f"
		*info = -10;
#line 473 "ssyevr_2stage.f"
	    }
#line 474 "ssyevr_2stage.f"
	}
#line 475 "ssyevr_2stage.f"
    }
#line 476 "ssyevr_2stage.f"
    if (*info == 0) {
#line 477 "ssyevr_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 478 "ssyevr_2stage.f"
	    *info = -15;
#line 479 "ssyevr_2stage.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 480 "ssyevr_2stage.f"
	    *info = -18;
#line 481 "ssyevr_2stage.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 482 "ssyevr_2stage.f"
	    *info = -20;
#line 483 "ssyevr_2stage.f"
	}
#line 484 "ssyevr_2stage.f"
    }

#line 486 "ssyevr_2stage.f"
    if (*info == 0) {
/*         NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 ) */
/*         NB = MAX( NB, ILAENV( 1, 'SORMTR', UPLO, N, -1, -1, -1 ) ) */
/*         LWKOPT = MAX( ( NB+1 )*N, LWMIN ) */
#line 490 "ssyevr_2stage.f"
	work[1] = (doublereal) lwmin;
#line 491 "ssyevr_2stage.f"
	iwork[1] = liwmin;
#line 492 "ssyevr_2stage.f"
    }

#line 494 "ssyevr_2stage.f"
    if (*info != 0) {
#line 495 "ssyevr_2stage.f"
	i__1 = -(*info);
#line 495 "ssyevr_2stage.f"
	xerbla_("SSYEVR_2STAGE", &i__1, (ftnlen)13);
#line 496 "ssyevr_2stage.f"
	return 0;
#line 497 "ssyevr_2stage.f"
    } else if (lquery) {
#line 498 "ssyevr_2stage.f"
	return 0;
#line 499 "ssyevr_2stage.f"
    }

/*     Quick return if possible */

#line 503 "ssyevr_2stage.f"
    *m = 0;
#line 504 "ssyevr_2stage.f"
    if (*n == 0) {
#line 505 "ssyevr_2stage.f"
	work[1] = 1.;
#line 506 "ssyevr_2stage.f"
	return 0;
#line 507 "ssyevr_2stage.f"
    }

#line 509 "ssyevr_2stage.f"
    if (*n == 1) {
#line 510 "ssyevr_2stage.f"
	work[1] = 26.;
#line 511 "ssyevr_2stage.f"
	if (alleig || indeig) {
#line 512 "ssyevr_2stage.f"
	    *m = 1;
#line 513 "ssyevr_2stage.f"
	    w[1] = a[a_dim1 + 1];
#line 514 "ssyevr_2stage.f"
	} else {
#line 515 "ssyevr_2stage.f"
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
#line 516 "ssyevr_2stage.f"
		*m = 1;
#line 517 "ssyevr_2stage.f"
		w[1] = a[a_dim1 + 1];
#line 518 "ssyevr_2stage.f"
	    }
#line 519 "ssyevr_2stage.f"
	}
#line 520 "ssyevr_2stage.f"
	if (wantz) {
#line 521 "ssyevr_2stage.f"
	    z__[z_dim1 + 1] = 1.;
#line 522 "ssyevr_2stage.f"
	    isuppz[1] = 1;
#line 523 "ssyevr_2stage.f"
	    isuppz[2] = 1;
#line 524 "ssyevr_2stage.f"
	}
#line 525 "ssyevr_2stage.f"
	return 0;
#line 526 "ssyevr_2stage.f"
    }

/*     Get machine constants. */

#line 530 "ssyevr_2stage.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 531 "ssyevr_2stage.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 532 "ssyevr_2stage.f"
    smlnum = safmin / eps;
#line 533 "ssyevr_2stage.f"
    bignum = 1. / smlnum;
#line 534 "ssyevr_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 535 "ssyevr_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 535 "ssyevr_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 539 "ssyevr_2stage.f"
    iscale = 0;
#line 540 "ssyevr_2stage.f"
    abstll = *abstol;
#line 541 "ssyevr_2stage.f"
    if (valeig) {
#line 542 "ssyevr_2stage.f"
	vll = *vl;
#line 543 "ssyevr_2stage.f"
	vuu = *vu;
#line 544 "ssyevr_2stage.f"
    }
#line 545 "ssyevr_2stage.f"
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 546 "ssyevr_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 547 "ssyevr_2stage.f"
	iscale = 1;
#line 548 "ssyevr_2stage.f"
	sigma = rmin / anrm;
#line 549 "ssyevr_2stage.f"
    } else if (anrm > rmax) {
#line 550 "ssyevr_2stage.f"
	iscale = 1;
#line 551 "ssyevr_2stage.f"
	sigma = rmax / anrm;
#line 552 "ssyevr_2stage.f"
    }
#line 553 "ssyevr_2stage.f"
    if (iscale == 1) {
#line 554 "ssyevr_2stage.f"
	if (lower) {
#line 555 "ssyevr_2stage.f"
	    i__1 = *n;
#line 555 "ssyevr_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 556 "ssyevr_2stage.f"
		i__2 = *n - j + 1;
#line 556 "ssyevr_2stage.f"
		sscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 557 "ssyevr_2stage.f"
/* L10: */
#line 557 "ssyevr_2stage.f"
	    }
#line 558 "ssyevr_2stage.f"
	} else {
#line 559 "ssyevr_2stage.f"
	    i__1 = *n;
#line 559 "ssyevr_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 560 "ssyevr_2stage.f"
		sscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 561 "ssyevr_2stage.f"
/* L20: */
#line 561 "ssyevr_2stage.f"
	    }
#line 562 "ssyevr_2stage.f"
	}
#line 563 "ssyevr_2stage.f"
	if (*abstol > 0.) {
#line 563 "ssyevr_2stage.f"
	    abstll = *abstol * sigma;
#line 563 "ssyevr_2stage.f"
	}
#line 565 "ssyevr_2stage.f"
	if (valeig) {
#line 566 "ssyevr_2stage.f"
	    vll = *vl * sigma;
#line 567 "ssyevr_2stage.f"
	    vuu = *vu * sigma;
#line 568 "ssyevr_2stage.f"
	}
#line 569 "ssyevr_2stage.f"
    }
/*     Initialize indices into workspaces.  Note: The IWORK indices are */
/*     used only if SSTERF or SSTEMR fail. */
/*     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the */
/*     elementary reflectors used in SSYTRD. */
#line 576 "ssyevr_2stage.f"
    indtau = 1;
/*     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries. */
#line 578 "ssyevr_2stage.f"
    indd = indtau + *n;
/*     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the */
/*     tridiagonal matrix from SSYTRD. */
#line 581 "ssyevr_2stage.f"
    inde = indd + *n;
/*     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over */
/*     -written by SSTEMR (the SSTERF path copies the diagonal to W). */
#line 584 "ssyevr_2stage.f"
    inddd = inde + *n;
/*     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over */
/*     -written while computing the eigenvalues in SSTERF and SSTEMR. */
#line 587 "ssyevr_2stage.f"
    indee = inddd + *n;
/*     INDHOUS is the starting offset Householder storage of stage 2 */
#line 589 "ssyevr_2stage.f"
    indhous = indee + *n;
/*     INDWK is the starting offset of the left-over workspace, and */
/*     LLWORK is the remaining workspace size. */
#line 592 "ssyevr_2stage.f"
    indwk = indhous + lhtrd;
#line 593 "ssyevr_2stage.f"
    llwork = *lwork - indwk + 1;
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 598 "ssyevr_2stage.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 601 "ssyevr_2stage.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     SSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 606 "ssyevr_2stage.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 608 "ssyevr_2stage.f"
    indiwo = indifl + *n;

/*     Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */


#line 614 "ssyevr_2stage.f"
    ssytrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &work[indd], &work[inde]
	    , &work[indtau], &work[indhous], &lhtrd, &work[indwk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired */
/*     then call SSTERF or SSTEMR and SORMTR. */

#line 621 "ssyevr_2stage.f"
    test = FALSE_;
#line 622 "ssyevr_2stage.f"
    if (indeig) {
#line 623 "ssyevr_2stage.f"
	if (*il == 1 && *iu == *n) {
#line 624 "ssyevr_2stage.f"
	    test = TRUE_;
#line 625 "ssyevr_2stage.f"
	}
#line 626 "ssyevr_2stage.f"
    }
#line 627 "ssyevr_2stage.f"
    if ((alleig || test) && ieeeok == 1) {
#line 628 "ssyevr_2stage.f"
	if (! wantz) {
#line 629 "ssyevr_2stage.f"
	    scopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 630 "ssyevr_2stage.f"
	    i__1 = *n - 1;
#line 630 "ssyevr_2stage.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 631 "ssyevr_2stage.f"
	    ssterf_(n, &w[1], &work[indee], info);
#line 632 "ssyevr_2stage.f"
	} else {
#line 633 "ssyevr_2stage.f"
	    i__1 = *n - 1;
#line 633 "ssyevr_2stage.f"
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 634 "ssyevr_2stage.f"
	    scopy_(n, &work[indd], &c__1, &work[inddd], &c__1);

#line 636 "ssyevr_2stage.f"
	    if (*abstol <= *n * 2. * eps) {
#line 637 "ssyevr_2stage.f"
		tryrac = TRUE_;
#line 638 "ssyevr_2stage.f"
	    } else {
#line 639 "ssyevr_2stage.f"
		tryrac = FALSE_;
#line 640 "ssyevr_2stage.f"
	    }
#line 641 "ssyevr_2stage.f"
	    sstemr_(jobz, "A", n, &work[inddd], &work[indee], vl, vu, il, iu, 
		    m, &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &
		    work[indwk], lwork, &iwork[1], liwork, info, (ftnlen)1, (
		    ftnlen)1);



/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEMR. */

#line 651 "ssyevr_2stage.f"
	    if (wantz && *info == 0) {
#line 652 "ssyevr_2stage.f"
		indwkn = inde;
#line 653 "ssyevr_2stage.f"
		llwrkn = *lwork - indwkn + 1;
#line 654 "ssyevr_2stage.f"
		sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
			, &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 657 "ssyevr_2stage.f"
	    }
#line 658 "ssyevr_2stage.f"
	}


#line 661 "ssyevr_2stage.f"
	if (*info == 0) {
/*           Everything worked.  Skip SSTEBZ/SSTEIN.  IWORK(:) are */
/*           undefined. */
#line 664 "ssyevr_2stage.f"
	    *m = *n;
#line 665 "ssyevr_2stage.f"
	    goto L30;
#line 666 "ssyevr_2stage.f"
	}
#line 667 "ssyevr_2stage.f"
	*info = 0;
#line 668 "ssyevr_2stage.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */
/*     Also call SSTEBZ and SSTEIN if SSTEMR fails. */

#line 673 "ssyevr_2stage.f"
    if (wantz) {
#line 674 "ssyevr_2stage.f"
	*(unsigned char *)order = 'B';
#line 675 "ssyevr_2stage.f"
    } else {
#line 676 "ssyevr_2stage.f"
	*(unsigned char *)order = 'E';
#line 677 "ssyevr_2stage.f"
    }
#line 679 "ssyevr_2stage.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 684 "ssyevr_2stage.f"
    if (wantz) {
#line 685 "ssyevr_2stage.f"
	sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwk], &iwork[indiwo], &
		iwork[indifl], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEIN. */

#line 693 "ssyevr_2stage.f"
	indwkn = inde;
#line 694 "ssyevr_2stage.f"
	llwrkn = *lwork - indwkn + 1;
#line 695 "ssyevr_2stage.f"
	sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 697 "ssyevr_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

/*  Jump here if SSTEMR/SSTEIN succeeded. */
#line 702 "ssyevr_2stage.f"
L30:
#line 703 "ssyevr_2stage.f"
    if (iscale == 1) {
#line 704 "ssyevr_2stage.f"
	if (*info == 0) {
#line 705 "ssyevr_2stage.f"
	    imax = *m;
#line 706 "ssyevr_2stage.f"
	} else {
#line 707 "ssyevr_2stage.f"
	    imax = *info - 1;
#line 708 "ssyevr_2stage.f"
	}
#line 709 "ssyevr_2stage.f"
	d__1 = 1. / sigma;
#line 709 "ssyevr_2stage.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 710 "ssyevr_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK. */
/*     It may not be initialized (if SSTEMR/SSTEIN succeeded), and we do */
/*     not return this detailed information to the user. */

#line 717 "ssyevr_2stage.f"
    if (wantz) {
#line 718 "ssyevr_2stage.f"
	i__1 = *m - 1;
#line 718 "ssyevr_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 719 "ssyevr_2stage.f"
	    i__ = 0;
#line 720 "ssyevr_2stage.f"
	    tmp1 = w[j];
#line 721 "ssyevr_2stage.f"
	    i__2 = *m;
#line 721 "ssyevr_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 722 "ssyevr_2stage.f"
		if (w[jj] < tmp1) {
#line 723 "ssyevr_2stage.f"
		    i__ = jj;
#line 724 "ssyevr_2stage.f"
		    tmp1 = w[jj];
#line 725 "ssyevr_2stage.f"
		}
#line 726 "ssyevr_2stage.f"
/* L40: */
#line 726 "ssyevr_2stage.f"
	    }

#line 728 "ssyevr_2stage.f"
	    if (i__ != 0) {
#line 729 "ssyevr_2stage.f"
		w[i__] = w[j];
#line 730 "ssyevr_2stage.f"
		w[j] = tmp1;
#line 731 "ssyevr_2stage.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 732 "ssyevr_2stage.f"
	    }
#line 733 "ssyevr_2stage.f"
/* L50: */
#line 733 "ssyevr_2stage.f"
	}
#line 734 "ssyevr_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 738 "ssyevr_2stage.f"
    work[1] = (doublereal) lwmin;
#line 739 "ssyevr_2stage.f"
    iwork[1] = liwmin;

#line 741 "ssyevr_2stage.f"
    return 0;

/*     End of SSYEVR_2STAGE */

} /* ssyevr_2stage__ */

