#line 1 "dsyevr_2stage.f"
/* dsyevr_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsyevr_2stage.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;

/* > \brief <b> DSYEVR_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 SY matrices</b> */

/*  @precisions fortran d -> s */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYEVR_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyevr_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyevr_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyevr_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, */
/*                          IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, */
/*                          LWORK, IWORK, LIWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYEVR_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  Eigenvalues and eigenvectors can be */
/* > selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
/* > */
/* > DSYEVR_2STAGE first reduces the matrix A to tridiagonal form T with a call */
/* > to DSYTRD.  Then, whenever possible, DSYEVR_2STAGE calls DSTEMR to compute */
/* > the eigenspectrum using Relatively Robust Representations.  DSTEMR */
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
/* > For more details, see DSTEMR's documentation and: */
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
/* > Note 1 : DSYEVR_2STAGE calls DSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > DSYEVR_2STAGE calls DSTEBZ and SSTEIN on non-ieee machines and */
/* > when partial spectrum requests are made. */
/* > */
/* > Normal execution of DSTEMR may create NaNs and infinities and */
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
/* >          For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and */
/* >          DSTEIN are called */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
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
/* >          See "Computing Small Singular Values of Bidiagonal Matrices */
/* >          with Guaranteed High Relative Accuracy," by Demmel and */
/* >          Kahan, LAPACK Working Note #3. */
/* > */
/* >          If high relative accuracy is important, set ABSTOL to */
/* >          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that */
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
/* >          ISUPPZ( 2*i ). This is an output of DSTEMR (tridiagonal */
/* >          matrix). The support of the eigenvectors of A is typically */
/* >          1:N because of the orthogonal transformations applied by DORMTR. */
/* >          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup doubleSYeigen */

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
/* Subroutine */ int dsyevr_2stage__(char *jobz, char *range, char *uplo, 
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
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static integer inddd, indee;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static char order[1];
    static integer indwk, lhtrd;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), dsytrd_2stage__(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static integer lwmin;
    static logical lower;
    static integer lwtrd;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, ieeeok, indibl, indifl;
    static logical valeig;
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp;
    extern /* Subroutine */ int dstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dsterf_(integer *, doublereal *, doublereal *, integer *);
    static integer indiwo, indwkn;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dstemr_(char *, char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen, ftnlen);
    static integer liwmin;
    static logical tryrac;
    extern /* Subroutine */ int dormtr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer llwrkn, llwork, nsplit;
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

#line 436 "dsyevr_2stage.f"
    /* Parameter adjustments */
#line 436 "dsyevr_2stage.f"
    a_dim1 = *lda;
#line 436 "dsyevr_2stage.f"
    a_offset = 1 + a_dim1;
#line 436 "dsyevr_2stage.f"
    a -= a_offset;
#line 436 "dsyevr_2stage.f"
    --w;
#line 436 "dsyevr_2stage.f"
    z_dim1 = *ldz;
#line 436 "dsyevr_2stage.f"
    z_offset = 1 + z_dim1;
#line 436 "dsyevr_2stage.f"
    z__ -= z_offset;
#line 436 "dsyevr_2stage.f"
    --isuppz;
#line 436 "dsyevr_2stage.f"
    --work;
#line 436 "dsyevr_2stage.f"
    --iwork;
#line 436 "dsyevr_2stage.f"

#line 436 "dsyevr_2stage.f"
    /* Function Body */
#line 436 "dsyevr_2stage.f"
    ieeeok = ilaenv_(&c__10, "DSYEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 438 "dsyevr_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 439 "dsyevr_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 440 "dsyevr_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 441 "dsyevr_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 442 "dsyevr_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 444 "dsyevr_2stage.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 446 "dsyevr_2stage.f"
    kd = ilaenv_(&c__17, "DSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 447 "dsyevr_2stage.f"
    ib = ilaenv_(&c__18, "DSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, (ftnlen)
	    13, (ftnlen)1);
#line 448 "dsyevr_2stage.f"
    lhtrd = ilaenv_(&c__19, "DSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 449 "dsyevr_2stage.f"
    lwtrd = ilaenv_(&c__20, "DSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
/* Computing MAX */
#line 450 "dsyevr_2stage.f"
    i__1 = *n * 26, i__2 = *n * 5 + lhtrd + lwtrd;
#line 450 "dsyevr_2stage.f"
    lwmin = max(i__1,i__2);
/* Computing MAX */
#line 451 "dsyevr_2stage.f"
    i__1 = 1, i__2 = *n * 10;
#line 451 "dsyevr_2stage.f"
    liwmin = max(i__1,i__2);

#line 453 "dsyevr_2stage.f"
    *info = 0;
#line 454 "dsyevr_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 455 "dsyevr_2stage.f"
	*info = -1;
#line 456 "dsyevr_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 457 "dsyevr_2stage.f"
	*info = -2;
#line 458 "dsyevr_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 459 "dsyevr_2stage.f"
	*info = -3;
#line 460 "dsyevr_2stage.f"
    } else if (*n < 0) {
#line 461 "dsyevr_2stage.f"
	*info = -4;
#line 462 "dsyevr_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 463 "dsyevr_2stage.f"
	*info = -6;
#line 464 "dsyevr_2stage.f"
    } else {
#line 465 "dsyevr_2stage.f"
	if (valeig) {
#line 466 "dsyevr_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 466 "dsyevr_2stage.f"
		*info = -8;
#line 466 "dsyevr_2stage.f"
	    }
#line 468 "dsyevr_2stage.f"
	} else if (indeig) {
#line 469 "dsyevr_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 470 "dsyevr_2stage.f"
		*info = -9;
#line 471 "dsyevr_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 472 "dsyevr_2stage.f"
		*info = -10;
#line 473 "dsyevr_2stage.f"
	    }
#line 474 "dsyevr_2stage.f"
	}
#line 475 "dsyevr_2stage.f"
    }
#line 476 "dsyevr_2stage.f"
    if (*info == 0) {
#line 477 "dsyevr_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 478 "dsyevr_2stage.f"
	    *info = -15;
#line 479 "dsyevr_2stage.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 480 "dsyevr_2stage.f"
	    *info = -18;
#line 481 "dsyevr_2stage.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 482 "dsyevr_2stage.f"
	    *info = -20;
#line 483 "dsyevr_2stage.f"
	}
#line 484 "dsyevr_2stage.f"
    }

#line 486 "dsyevr_2stage.f"
    if (*info == 0) {
/*         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 ) */
/*         NB = MAX( NB, ILAENV( 1, 'DORMTR', UPLO, N, -1, -1, -1 ) ) */
/*         LWKOPT = MAX( ( NB+1 )*N, LWMIN ) */
#line 490 "dsyevr_2stage.f"
	work[1] = (doublereal) lwmin;
#line 491 "dsyevr_2stage.f"
	iwork[1] = liwmin;
#line 492 "dsyevr_2stage.f"
    }

#line 494 "dsyevr_2stage.f"
    if (*info != 0) {
#line 495 "dsyevr_2stage.f"
	i__1 = -(*info);
#line 495 "dsyevr_2stage.f"
	xerbla_("DSYEVR_2STAGE", &i__1, (ftnlen)13);
#line 496 "dsyevr_2stage.f"
	return 0;
#line 497 "dsyevr_2stage.f"
    } else if (lquery) {
#line 498 "dsyevr_2stage.f"
	return 0;
#line 499 "dsyevr_2stage.f"
    }

/*     Quick return if possible */

#line 503 "dsyevr_2stage.f"
    *m = 0;
#line 504 "dsyevr_2stage.f"
    if (*n == 0) {
#line 505 "dsyevr_2stage.f"
	work[1] = 1.;
#line 506 "dsyevr_2stage.f"
	return 0;
#line 507 "dsyevr_2stage.f"
    }

#line 509 "dsyevr_2stage.f"
    if (*n == 1) {
#line 510 "dsyevr_2stage.f"
	work[1] = 7.;
#line 511 "dsyevr_2stage.f"
	if (alleig || indeig) {
#line 512 "dsyevr_2stage.f"
	    *m = 1;
#line 513 "dsyevr_2stage.f"
	    w[1] = a[a_dim1 + 1];
#line 514 "dsyevr_2stage.f"
	} else {
#line 515 "dsyevr_2stage.f"
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
#line 516 "dsyevr_2stage.f"
		*m = 1;
#line 517 "dsyevr_2stage.f"
		w[1] = a[a_dim1 + 1];
#line 518 "dsyevr_2stage.f"
	    }
#line 519 "dsyevr_2stage.f"
	}
#line 520 "dsyevr_2stage.f"
	if (wantz) {
#line 521 "dsyevr_2stage.f"
	    z__[z_dim1 + 1] = 1.;
#line 522 "dsyevr_2stage.f"
	    isuppz[1] = 1;
#line 523 "dsyevr_2stage.f"
	    isuppz[2] = 1;
#line 524 "dsyevr_2stage.f"
	}
#line 525 "dsyevr_2stage.f"
	return 0;
#line 526 "dsyevr_2stage.f"
    }

/*     Get machine constants. */

#line 530 "dsyevr_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 531 "dsyevr_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 532 "dsyevr_2stage.f"
    smlnum = safmin / eps;
#line 533 "dsyevr_2stage.f"
    bignum = 1. / smlnum;
#line 534 "dsyevr_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 535 "dsyevr_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 535 "dsyevr_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 539 "dsyevr_2stage.f"
    iscale = 0;
#line 540 "dsyevr_2stage.f"
    abstll = *abstol;
#line 541 "dsyevr_2stage.f"
    if (valeig) {
#line 542 "dsyevr_2stage.f"
	vll = *vl;
#line 543 "dsyevr_2stage.f"
	vuu = *vu;
#line 544 "dsyevr_2stage.f"
    }
#line 545 "dsyevr_2stage.f"
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 546 "dsyevr_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 547 "dsyevr_2stage.f"
	iscale = 1;
#line 548 "dsyevr_2stage.f"
	sigma = rmin / anrm;
#line 549 "dsyevr_2stage.f"
    } else if (anrm > rmax) {
#line 550 "dsyevr_2stage.f"
	iscale = 1;
#line 551 "dsyevr_2stage.f"
	sigma = rmax / anrm;
#line 552 "dsyevr_2stage.f"
    }
#line 553 "dsyevr_2stage.f"
    if (iscale == 1) {
#line 554 "dsyevr_2stage.f"
	if (lower) {
#line 555 "dsyevr_2stage.f"
	    i__1 = *n;
#line 555 "dsyevr_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 556 "dsyevr_2stage.f"
		i__2 = *n - j + 1;
#line 556 "dsyevr_2stage.f"
		dscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 557 "dsyevr_2stage.f"
/* L10: */
#line 557 "dsyevr_2stage.f"
	    }
#line 558 "dsyevr_2stage.f"
	} else {
#line 559 "dsyevr_2stage.f"
	    i__1 = *n;
#line 559 "dsyevr_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 560 "dsyevr_2stage.f"
		dscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 561 "dsyevr_2stage.f"
/* L20: */
#line 561 "dsyevr_2stage.f"
	    }
#line 562 "dsyevr_2stage.f"
	}
#line 563 "dsyevr_2stage.f"
	if (*abstol > 0.) {
#line 563 "dsyevr_2stage.f"
	    abstll = *abstol * sigma;
#line 563 "dsyevr_2stage.f"
	}
#line 565 "dsyevr_2stage.f"
	if (valeig) {
#line 566 "dsyevr_2stage.f"
	    vll = *vl * sigma;
#line 567 "dsyevr_2stage.f"
	    vuu = *vu * sigma;
#line 568 "dsyevr_2stage.f"
	}
#line 569 "dsyevr_2stage.f"
    }
/*     Initialize indices into workspaces.  Note: The IWORK indices are */
/*     used only if DSTERF or DSTEMR fail. */
/*     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the */
/*     elementary reflectors used in DSYTRD. */
#line 576 "dsyevr_2stage.f"
    indtau = 1;
/*     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries. */
#line 578 "dsyevr_2stage.f"
    indd = indtau + *n;
/*     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the */
/*     tridiagonal matrix from DSYTRD. */
#line 581 "dsyevr_2stage.f"
    inde = indd + *n;
/*     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over */
/*     -written by DSTEMR (the DSTERF path copies the diagonal to W). */
#line 584 "dsyevr_2stage.f"
    inddd = inde + *n;
/*     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over */
/*     -written while computing the eigenvalues in DSTERF and DSTEMR. */
#line 587 "dsyevr_2stage.f"
    indee = inddd + *n;
/*     INDHOUS is the starting offset Householder storage of stage 2 */
#line 589 "dsyevr_2stage.f"
    indhous = indee + *n;
/*     INDWK is the starting offset of the left-over workspace, and */
/*     LLWORK is the remaining workspace size. */
#line 592 "dsyevr_2stage.f"
    indwk = indhous + lhtrd;
#line 593 "dsyevr_2stage.f"
    llwork = *lwork - indwk + 1;
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 598 "dsyevr_2stage.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 601 "dsyevr_2stage.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     DSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 606 "dsyevr_2stage.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 608 "dsyevr_2stage.f"
    indiwo = indifl + *n;

/*     Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */


#line 614 "dsyevr_2stage.f"
    dsytrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &work[indd], &work[inde]
	    , &work[indtau], &work[indhous], &lhtrd, &work[indwk], &llwork, &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired */
/*     then call DSTERF or DSTEMR and DORMTR. */

#line 621 "dsyevr_2stage.f"
    if ((alleig || indeig && *il == 1 && *iu == *n) && ieeeok == 1) {
#line 623 "dsyevr_2stage.f"
	if (! wantz) {
#line 624 "dsyevr_2stage.f"
	    dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 625 "dsyevr_2stage.f"
	    i__1 = *n - 1;
#line 625 "dsyevr_2stage.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 626 "dsyevr_2stage.f"
	    dsterf_(n, &w[1], &work[indee], info);
#line 627 "dsyevr_2stage.f"
	} else {
#line 628 "dsyevr_2stage.f"
	    i__1 = *n - 1;
#line 628 "dsyevr_2stage.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 629 "dsyevr_2stage.f"
	    dcopy_(n, &work[indd], &c__1, &work[inddd], &c__1);

#line 631 "dsyevr_2stage.f"
	    if (*abstol <= *n * 2. * eps) {
#line 632 "dsyevr_2stage.f"
		tryrac = TRUE_;
#line 633 "dsyevr_2stage.f"
	    } else {
#line 634 "dsyevr_2stage.f"
		tryrac = FALSE_;
#line 635 "dsyevr_2stage.f"
	    }
#line 636 "dsyevr_2stage.f"
	    dstemr_(jobz, "A", n, &work[inddd], &work[indee], vl, vu, il, iu, 
		    m, &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &
		    work[indwk], lwork, &iwork[1], liwork, info, (ftnlen)1, (
		    ftnlen)1);



/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEMR. */

#line 646 "dsyevr_2stage.f"
	    if (wantz && *info == 0) {
#line 647 "dsyevr_2stage.f"
		indwkn = inde;
#line 648 "dsyevr_2stage.f"
		llwrkn = *lwork - indwkn + 1;
#line 649 "dsyevr_2stage.f"
		dormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
			, &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 652 "dsyevr_2stage.f"
	    }
#line 653 "dsyevr_2stage.f"
	}


#line 656 "dsyevr_2stage.f"
	if (*info == 0) {
/*           Everything worked.  Skip DSTEBZ/DSTEIN.  IWORK(:) are */
/*           undefined. */
#line 659 "dsyevr_2stage.f"
	    *m = *n;
#line 660 "dsyevr_2stage.f"
	    goto L30;
#line 661 "dsyevr_2stage.f"
	}
#line 662 "dsyevr_2stage.f"
	*info = 0;
#line 663 "dsyevr_2stage.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN. */
/*     Also call DSTEBZ and DSTEIN if DSTEMR fails. */

#line 668 "dsyevr_2stage.f"
    if (wantz) {
#line 669 "dsyevr_2stage.f"
	*(unsigned char *)order = 'B';
#line 670 "dsyevr_2stage.f"
    } else {
#line 671 "dsyevr_2stage.f"
	*(unsigned char *)order = 'E';
#line 672 "dsyevr_2stage.f"
    }
#line 674 "dsyevr_2stage.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 679 "dsyevr_2stage.f"
    if (wantz) {
#line 680 "dsyevr_2stage.f"
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwk], &iwork[indiwo], &
		iwork[indifl], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

#line 688 "dsyevr_2stage.f"
	indwkn = inde;
#line 689 "dsyevr_2stage.f"
	llwrkn = *lwork - indwkn + 1;
#line 690 "dsyevr_2stage.f"
	dormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 692 "dsyevr_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

/*  Jump here if DSTEMR/DSTEIN succeeded. */
#line 697 "dsyevr_2stage.f"
L30:
#line 698 "dsyevr_2stage.f"
    if (iscale == 1) {
#line 699 "dsyevr_2stage.f"
	if (*info == 0) {
#line 700 "dsyevr_2stage.f"
	    imax = *m;
#line 701 "dsyevr_2stage.f"
	} else {
#line 702 "dsyevr_2stage.f"
	    imax = *info - 1;
#line 703 "dsyevr_2stage.f"
	}
#line 704 "dsyevr_2stage.f"
	d__1 = 1. / sigma;
#line 704 "dsyevr_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 705 "dsyevr_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK. */
/*     It may not be initialized (if DSTEMR/DSTEIN succeeded), and we do */
/*     not return this detailed information to the user. */

#line 712 "dsyevr_2stage.f"
    if (wantz) {
#line 713 "dsyevr_2stage.f"
	i__1 = *m - 1;
#line 713 "dsyevr_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 714 "dsyevr_2stage.f"
	    i__ = 0;
#line 715 "dsyevr_2stage.f"
	    tmp1 = w[j];
#line 716 "dsyevr_2stage.f"
	    i__2 = *m;
#line 716 "dsyevr_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 717 "dsyevr_2stage.f"
		if (w[jj] < tmp1) {
#line 718 "dsyevr_2stage.f"
		    i__ = jj;
#line 719 "dsyevr_2stage.f"
		    tmp1 = w[jj];
#line 720 "dsyevr_2stage.f"
		}
#line 721 "dsyevr_2stage.f"
/* L40: */
#line 721 "dsyevr_2stage.f"
	    }

#line 723 "dsyevr_2stage.f"
	    if (i__ != 0) {
#line 724 "dsyevr_2stage.f"
		w[i__] = w[j];
#line 725 "dsyevr_2stage.f"
		w[j] = tmp1;
#line 726 "dsyevr_2stage.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 727 "dsyevr_2stage.f"
	    }
#line 728 "dsyevr_2stage.f"
/* L50: */
#line 728 "dsyevr_2stage.f"
	}
#line 729 "dsyevr_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 733 "dsyevr_2stage.f"
    work[1] = (doublereal) lwmin;
#line 734 "dsyevr_2stage.f"
    iwork[1] = liwmin;

#line 736 "dsyevr_2stage.f"
    return 0;

/*     End of DSYEVR_2STAGE */

} /* dsyevr_2stage__ */

