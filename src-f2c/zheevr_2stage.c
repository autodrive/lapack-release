#line 1 "zheevr_2stage.f"
/* zheevr_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "zheevr_2stage.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;

/* > \brief <b> ZHEEVR_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for
 HE matrices</b> */

/*  @precisions fortran z -> s d c */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEVR_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevr_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevr_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevr_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, */
/*                                 IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, */
/*                                 WORK, LWORK, RWORK, LRWORK, IWORK, */
/*                                 LIWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, */
/*      $                   M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEVR_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian matrix A using the 2stage technique for */
/* > the reduction to tridiagonal.  Eigenvalues and eigenvectors can */
/* > be selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
/* > */
/* > ZHEEVR_2STAGE first reduces the matrix A to tridiagonal form T with a call */
/* > to ZHETRD.  Then, whenever possible, ZHEEVR_2STAGE calls ZSTEMR to compute */
/* > eigenspectrum using Relatively Robust Representations.  ZSTEMR */
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
/* > Note 1 : ZHEEVR_2STAGE calls ZSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > ZHEEVR_2STAGE calls DSTEBZ and ZSTEIN on non-ieee machines and */
/* > when partial spectrum requests are made. */
/* > */
/* > Normal execution of ZSTEMR may create NaNs and infinities and */
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
/* >          ZSTEIN are called */
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
/* >          See "Computing Small Singular Values of Bidiagonal Matrices */
/* >          with Guaranteed High Relative Accuracy," by Demmel and */
/* >          Kahan, LAPACK Working Note #3. */
/* > */
/* >          If high relative accuracy is important, set ABSTOL to */
/* >          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that */
/* >          eigenvalues are computed to high relative accuracy when */
/* >          possible in future releases.  The current code does not */
/* >          make any guarantees about high relative accuracy, but */
/* >          furutre releases will. See J. Barlow and J. Demmel, */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ, max(1,M)) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
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
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* >          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) ) */
/* >          The support of the eigenvectors in Z, i.e., the indices */
/* >          indicating the nonzero elements in Z. The i-th eigenvector */
/* >          is nonzero only in elements ISUPPZ( 2*i-1 ) through */
/* >          ISUPPZ( 2*i ). This is an output of ZSTEMR (tridiagonal */
/* >          matrix). The support of the eigenvectors of A is typically */
/* >          1:N because of the unitary transformations applied by ZUNMTR. */
/* >          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 */
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
/* >          The dimension of the array WORK. */
/* >          If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* >                                   LWORK = MAX(1, 26*N, dimension) where */
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
/* >          only calculates the optimal sizes of the WORK, RWORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
/* >          On exit, if INFO = 0, RWORK(1) returns the optimal */
/* >          (and minimal) LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The length of the array RWORK.  LRWORK >= max(1,24*N). */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal */
/* >          (and minimal) LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK.  LIWORK >= max(1,10*N). */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
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

/* > \ingroup complex16HEeigen */

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
/* Subroutine */ int zheevr_2stage__(char *jobz, char *range, char *uplo, 
	integer *n, doublecomplex *a, integer *lda, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublecomplex *z__, integer *ldz, integer *isuppz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info, ftnlen 
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
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static integer itmp1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer indrd, indre;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int zhetrd_2stage__(char *, char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char order[1];
    static integer indwk, lhtrd;
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
    static integer iscale, ieeeok, indibl, indrdd, indifl, indree;
    static logical valeig;
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal abstll, bignum;
    static integer indtau, indisp;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer indiwo, indwkn;
    extern /* Subroutine */ int dstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indrwk, liwmin;
    static logical tryrac;
    static integer lrwmin, llwrkn, llwork, nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int zstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static logical lquery;
    extern doublereal zlansy_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int zstemr_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     integer *, doublereal *, doublecomplex *, integer *, integer *, 
	    integer *, logical *, doublereal *, integer *, integer *, integer 
	    *, integer *, ftnlen, ftnlen), zunmtr_(char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen);
    static integer indhous, llrwork;



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

#line 463 "zheevr_2stage.f"
    /* Parameter adjustments */
#line 463 "zheevr_2stage.f"
    a_dim1 = *lda;
#line 463 "zheevr_2stage.f"
    a_offset = 1 + a_dim1;
#line 463 "zheevr_2stage.f"
    a -= a_offset;
#line 463 "zheevr_2stage.f"
    --w;
#line 463 "zheevr_2stage.f"
    z_dim1 = *ldz;
#line 463 "zheevr_2stage.f"
    z_offset = 1 + z_dim1;
#line 463 "zheevr_2stage.f"
    z__ -= z_offset;
#line 463 "zheevr_2stage.f"
    --isuppz;
#line 463 "zheevr_2stage.f"
    --work;
#line 463 "zheevr_2stage.f"
    --rwork;
#line 463 "zheevr_2stage.f"
    --iwork;
#line 463 "zheevr_2stage.f"

#line 463 "zheevr_2stage.f"
    /* Function Body */
#line 463 "zheevr_2stage.f"
    ieeeok = ilaenv_(&c__10, "ZHEEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 465 "zheevr_2stage.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 466 "zheevr_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 467 "zheevr_2stage.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 468 "zheevr_2stage.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 469 "zheevr_2stage.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 471 "zheevr_2stage.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 474 "zheevr_2stage.f"
    kd = ilaenv2stage_(&c__1, "ZHETRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 475 "zheevr_2stage.f"
    ib = ilaenv2stage_(&c__2, "ZHETRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 476 "zheevr_2stage.f"
    lhtrd = ilaenv2stage_(&c__3, "ZHETRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 477 "zheevr_2stage.f"
    lwtrd = ilaenv2stage_(&c__4, "ZHETRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 478 "zheevr_2stage.f"
    lwmin = *n + lhtrd + lwtrd;
/* Computing MAX */
#line 479 "zheevr_2stage.f"
    i__1 = 1, i__2 = *n * 24;
#line 479 "zheevr_2stage.f"
    lrwmin = max(i__1,i__2);
/* Computing MAX */
#line 480 "zheevr_2stage.f"
    i__1 = 1, i__2 = *n * 10;
#line 480 "zheevr_2stage.f"
    liwmin = max(i__1,i__2);

#line 482 "zheevr_2stage.f"
    *info = 0;
#line 483 "zheevr_2stage.f"
    if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 484 "zheevr_2stage.f"
	*info = -1;
#line 485 "zheevr_2stage.f"
    } else if (! (alleig || valeig || indeig)) {
#line 486 "zheevr_2stage.f"
	*info = -2;
#line 487 "zheevr_2stage.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 488 "zheevr_2stage.f"
	*info = -3;
#line 489 "zheevr_2stage.f"
    } else if (*n < 0) {
#line 490 "zheevr_2stage.f"
	*info = -4;
#line 491 "zheevr_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 492 "zheevr_2stage.f"
	*info = -6;
#line 493 "zheevr_2stage.f"
    } else {
#line 494 "zheevr_2stage.f"
	if (valeig) {
#line 495 "zheevr_2stage.f"
	    if (*n > 0 && *vu <= *vl) {
#line 495 "zheevr_2stage.f"
		*info = -8;
#line 495 "zheevr_2stage.f"
	    }
#line 497 "zheevr_2stage.f"
	} else if (indeig) {
#line 498 "zheevr_2stage.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 499 "zheevr_2stage.f"
		*info = -9;
#line 500 "zheevr_2stage.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 501 "zheevr_2stage.f"
		*info = -10;
#line 502 "zheevr_2stage.f"
	    }
#line 503 "zheevr_2stage.f"
	}
#line 504 "zheevr_2stage.f"
    }
#line 505 "zheevr_2stage.f"
    if (*info == 0) {
#line 506 "zheevr_2stage.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 507 "zheevr_2stage.f"
	    *info = -15;
#line 508 "zheevr_2stage.f"
	}
#line 509 "zheevr_2stage.f"
    }

#line 511 "zheevr_2stage.f"
    if (*info == 0) {
#line 512 "zheevr_2stage.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 513 "zheevr_2stage.f"
	rwork[1] = (doublereal) lrwmin;
#line 514 "zheevr_2stage.f"
	iwork[1] = liwmin;

#line 516 "zheevr_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 517 "zheevr_2stage.f"
	    *info = -18;
#line 518 "zheevr_2stage.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 519 "zheevr_2stage.f"
	    *info = -20;
#line 520 "zheevr_2stage.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 521 "zheevr_2stage.f"
	    *info = -22;
#line 522 "zheevr_2stage.f"
	}
#line 523 "zheevr_2stage.f"
    }

#line 525 "zheevr_2stage.f"
    if (*info != 0) {
#line 526 "zheevr_2stage.f"
	i__1 = -(*info);
#line 526 "zheevr_2stage.f"
	xerbla_("ZHEEVR_2STAGE", &i__1, (ftnlen)13);
#line 527 "zheevr_2stage.f"
	return 0;
#line 528 "zheevr_2stage.f"
    } else if (lquery) {
#line 529 "zheevr_2stage.f"
	return 0;
#line 530 "zheevr_2stage.f"
    }

/*     Quick return if possible */

#line 534 "zheevr_2stage.f"
    *m = 0;
#line 535 "zheevr_2stage.f"
    if (*n == 0) {
#line 536 "zheevr_2stage.f"
	work[1].r = 1., work[1].i = 0.;
#line 537 "zheevr_2stage.f"
	return 0;
#line 538 "zheevr_2stage.f"
    }

#line 540 "zheevr_2stage.f"
    if (*n == 1) {
#line 541 "zheevr_2stage.f"
	work[1].r = 2., work[1].i = 0.;
#line 542 "zheevr_2stage.f"
	if (alleig || indeig) {
#line 543 "zheevr_2stage.f"
	    *m = 1;
#line 544 "zheevr_2stage.f"
	    i__1 = a_dim1 + 1;
#line 544 "zheevr_2stage.f"
	    w[1] = a[i__1].r;
#line 545 "zheevr_2stage.f"
	} else {
#line 546 "zheevr_2stage.f"
	    i__1 = a_dim1 + 1;
#line 546 "zheevr_2stage.f"
	    i__2 = a_dim1 + 1;
#line 546 "zheevr_2stage.f"
	    if (*vl < a[i__1].r && *vu >= a[i__2].r) {
#line 548 "zheevr_2stage.f"
		*m = 1;
#line 549 "zheevr_2stage.f"
		i__1 = a_dim1 + 1;
#line 549 "zheevr_2stage.f"
		w[1] = a[i__1].r;
#line 550 "zheevr_2stage.f"
	    }
#line 551 "zheevr_2stage.f"
	}
#line 552 "zheevr_2stage.f"
	if (wantz) {
#line 553 "zheevr_2stage.f"
	    i__1 = z_dim1 + 1;
#line 553 "zheevr_2stage.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 554 "zheevr_2stage.f"
	    isuppz[1] = 1;
#line 555 "zheevr_2stage.f"
	    isuppz[2] = 1;
#line 556 "zheevr_2stage.f"
	}
#line 557 "zheevr_2stage.f"
	return 0;
#line 558 "zheevr_2stage.f"
    }

/*     Get machine constants. */

#line 562 "zheevr_2stage.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 563 "zheevr_2stage.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 564 "zheevr_2stage.f"
    smlnum = safmin / eps;
#line 565 "zheevr_2stage.f"
    bignum = 1. / smlnum;
#line 566 "zheevr_2stage.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 567 "zheevr_2stage.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 567 "zheevr_2stage.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 571 "zheevr_2stage.f"
    iscale = 0;
#line 572 "zheevr_2stage.f"
    abstll = *abstol;
#line 573 "zheevr_2stage.f"
    if (valeig) {
#line 574 "zheevr_2stage.f"
	vll = *vl;
#line 575 "zheevr_2stage.f"
	vuu = *vu;
#line 576 "zheevr_2stage.f"
    }
#line 577 "zheevr_2stage.f"
    anrm = zlansy_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 578 "zheevr_2stage.f"
    if (anrm > 0. && anrm < rmin) {
#line 579 "zheevr_2stage.f"
	iscale = 1;
#line 580 "zheevr_2stage.f"
	sigma = rmin / anrm;
#line 581 "zheevr_2stage.f"
    } else if (anrm > rmax) {
#line 582 "zheevr_2stage.f"
	iscale = 1;
#line 583 "zheevr_2stage.f"
	sigma = rmax / anrm;
#line 584 "zheevr_2stage.f"
    }
#line 585 "zheevr_2stage.f"
    if (iscale == 1) {
#line 586 "zheevr_2stage.f"
	if (lower) {
#line 587 "zheevr_2stage.f"
	    i__1 = *n;
#line 587 "zheevr_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 588 "zheevr_2stage.f"
		i__2 = *n - j + 1;
#line 588 "zheevr_2stage.f"
		zdscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 589 "zheevr_2stage.f"
/* L10: */
#line 589 "zheevr_2stage.f"
	    }
#line 590 "zheevr_2stage.f"
	} else {
#line 591 "zheevr_2stage.f"
	    i__1 = *n;
#line 591 "zheevr_2stage.f"
	    for (j = 1; j <= i__1; ++j) {
#line 592 "zheevr_2stage.f"
		zdscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 593 "zheevr_2stage.f"
/* L20: */
#line 593 "zheevr_2stage.f"
	    }
#line 594 "zheevr_2stage.f"
	}
#line 595 "zheevr_2stage.f"
	if (*abstol > 0.) {
#line 595 "zheevr_2stage.f"
	    abstll = *abstol * sigma;
#line 595 "zheevr_2stage.f"
	}
#line 597 "zheevr_2stage.f"
	if (valeig) {
#line 598 "zheevr_2stage.f"
	    vll = *vl * sigma;
#line 599 "zheevr_2stage.f"
	    vuu = *vu * sigma;
#line 600 "zheevr_2stage.f"
	}
#line 601 "zheevr_2stage.f"
    }
/*     Initialize indices into workspaces.  Note: The IWORK indices are */
/*     used only if DSTERF or ZSTEMR fail. */
/*     WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the */
/*     elementary reflectors used in ZHETRD. */
#line 608 "zheevr_2stage.f"
    indtau = 1;
/*     INDWK is the starting offset of the remaining complex workspace, */
/*     and LLWORK is the remaining complex workspace size. */
#line 611 "zheevr_2stage.f"
    indhous = indtau + *n;
#line 612 "zheevr_2stage.f"
    indwk = indhous + lhtrd;
#line 613 "zheevr_2stage.f"
    llwork = *lwork - indwk + 1;
/*     RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal */
/*     entries. */
#line 617 "zheevr_2stage.f"
    indrd = 1;
/*     RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the */
/*     tridiagonal matrix from ZHETRD. */
#line 620 "zheevr_2stage.f"
    indre = indrd + *n;
/*     RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over */
/*     -written by ZSTEMR (the DSTERF path copies the diagonal to W). */
#line 623 "zheevr_2stage.f"
    indrdd = indre + *n;
/*     RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over */
/*     -written while computing the eigenvalues in DSTERF and ZSTEMR. */
#line 626 "zheevr_2stage.f"
    indree = indrdd + *n;
/*     INDRWK is the starting offset of the left-over real workspace, and */
/*     LLRWORK is the remaining workspace size. */
#line 629 "zheevr_2stage.f"
    indrwk = indree + *n;
#line 630 "zheevr_2stage.f"
    llrwork = *lrwork - indrwk + 1;
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 634 "zheevr_2stage.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 637 "zheevr_2stage.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     ZSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 642 "zheevr_2stage.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 644 "zheevr_2stage.f"
    indiwo = indifl + *n;

/*     Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form. */

#line 649 "zheevr_2stage.f"
    zhetrd_2stage__(jobz, uplo, n, &a[a_offset], lda, &rwork[indrd], &rwork[
	    indre], &work[indtau], &work[indhous], &lhtrd, &work[indwk], &
	    llwork, &iinfo, (ftnlen)1, (ftnlen)1);

/*     If all eigenvalues are desired */
/*     then call DSTERF or ZSTEMR and ZUNMTR. */

#line 657 "zheevr_2stage.f"
    test = FALSE_;
#line 658 "zheevr_2stage.f"
    if (indeig) {
#line 659 "zheevr_2stage.f"
	if (*il == 1 && *iu == *n) {
#line 660 "zheevr_2stage.f"
	    test = TRUE_;
#line 661 "zheevr_2stage.f"
	}
#line 662 "zheevr_2stage.f"
    }
#line 663 "zheevr_2stage.f"
    if ((alleig || test) && ieeeok == 1) {
#line 664 "zheevr_2stage.f"
	if (! wantz) {
#line 665 "zheevr_2stage.f"
	    dcopy_(n, &rwork[indrd], &c__1, &w[1], &c__1);
#line 666 "zheevr_2stage.f"
	    i__1 = *n - 1;
#line 666 "zheevr_2stage.f"
	    dcopy_(&i__1, &rwork[indre], &c__1, &rwork[indree], &c__1);
#line 667 "zheevr_2stage.f"
	    dsterf_(n, &w[1], &rwork[indree], info);
#line 668 "zheevr_2stage.f"
	} else {
#line 669 "zheevr_2stage.f"
	    i__1 = *n - 1;
#line 669 "zheevr_2stage.f"
	    dcopy_(&i__1, &rwork[indre], &c__1, &rwork[indree], &c__1);
#line 670 "zheevr_2stage.f"
	    dcopy_(n, &rwork[indrd], &c__1, &rwork[indrdd], &c__1);

#line 672 "zheevr_2stage.f"
	    if (*abstol <= *n * 2. * eps) {
#line 673 "zheevr_2stage.f"
		tryrac = TRUE_;
#line 674 "zheevr_2stage.f"
	    } else {
#line 675 "zheevr_2stage.f"
		tryrac = FALSE_;
#line 676 "zheevr_2stage.f"
	    }
#line 677 "zheevr_2stage.f"
	    zstemr_(jobz, "A", n, &rwork[indrdd], &rwork[indree], vl, vu, il, 
		    iu, m, &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac,
		     &rwork[indrwk], &llrwork, &iwork[1], liwork, info, (
		    ftnlen)1, (ftnlen)1);

/*           Apply unitary matrix used in reduction to tridiagonal */
/*           form to eigenvectors returned by ZSTEMR. */

#line 686 "zheevr_2stage.f"
	    if (wantz && *info == 0) {
#line 687 "zheevr_2stage.f"
		indwkn = indwk;
#line 688 "zheevr_2stage.f"
		llwrkn = *lwork - indwkn + 1;
#line 689 "zheevr_2stage.f"
		zunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
			, &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 692 "zheevr_2stage.f"
	    }
#line 693 "zheevr_2stage.f"
	}


#line 696 "zheevr_2stage.f"
	if (*info == 0) {
#line 697 "zheevr_2stage.f"
	    *m = *n;
#line 698 "zheevr_2stage.f"
	    goto L30;
#line 699 "zheevr_2stage.f"
	}
#line 700 "zheevr_2stage.f"
	*info = 0;
#line 701 "zheevr_2stage.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN. */
/*     Also call DSTEBZ and ZSTEIN if ZSTEMR fails. */

#line 706 "zheevr_2stage.f"
    if (wantz) {
#line 707 "zheevr_2stage.f"
	*(unsigned char *)order = 'B';
#line 708 "zheevr_2stage.f"
    } else {
#line 709 "zheevr_2stage.f"
	*(unsigned char *)order = 'E';
#line 710 "zheevr_2stage.f"
    }
#line 712 "zheevr_2stage.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indrd], &
	    rwork[indre], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 717 "zheevr_2stage.f"
    if (wantz) {
#line 718 "zheevr_2stage.f"
	zstein_(n, &rwork[indrd], &rwork[indre], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwo], &iwork[indifl], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by ZSTEIN. */

#line 726 "zheevr_2stage.f"
	indwkn = indwk;
#line 727 "zheevr_2stage.f"
	llwrkn = *lwork - indwkn + 1;
#line 728 "zheevr_2stage.f"
	zunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 730 "zheevr_2stage.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 734 "zheevr_2stage.f"
L30:
#line 735 "zheevr_2stage.f"
    if (iscale == 1) {
#line 736 "zheevr_2stage.f"
	if (*info == 0) {
#line 737 "zheevr_2stage.f"
	    imax = *m;
#line 738 "zheevr_2stage.f"
	} else {
#line 739 "zheevr_2stage.f"
	    imax = *info - 1;
#line 740 "zheevr_2stage.f"
	}
#line 741 "zheevr_2stage.f"
	d__1 = 1. / sigma;
#line 741 "zheevr_2stage.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 742 "zheevr_2stage.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 747 "zheevr_2stage.f"
    if (wantz) {
#line 748 "zheevr_2stage.f"
	i__1 = *m - 1;
#line 748 "zheevr_2stage.f"
	for (j = 1; j <= i__1; ++j) {
#line 749 "zheevr_2stage.f"
	    i__ = 0;
#line 750 "zheevr_2stage.f"
	    tmp1 = w[j];
#line 751 "zheevr_2stage.f"
	    i__2 = *m;
#line 751 "zheevr_2stage.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 752 "zheevr_2stage.f"
		if (w[jj] < tmp1) {
#line 753 "zheevr_2stage.f"
		    i__ = jj;
#line 754 "zheevr_2stage.f"
		    tmp1 = w[jj];
#line 755 "zheevr_2stage.f"
		}
#line 756 "zheevr_2stage.f"
/* L40: */
#line 756 "zheevr_2stage.f"
	    }

#line 758 "zheevr_2stage.f"
	    if (i__ != 0) {
#line 759 "zheevr_2stage.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 760 "zheevr_2stage.f"
		w[i__] = w[j];
#line 761 "zheevr_2stage.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 762 "zheevr_2stage.f"
		w[j] = tmp1;
#line 763 "zheevr_2stage.f"
		iwork[indibl + j - 1] = itmp1;
#line 764 "zheevr_2stage.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 765 "zheevr_2stage.f"
	    }
#line 766 "zheevr_2stage.f"
/* L50: */
#line 766 "zheevr_2stage.f"
	}
#line 767 "zheevr_2stage.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 771 "zheevr_2stage.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 772 "zheevr_2stage.f"
    rwork[1] = (doublereal) lrwmin;
#line 773 "zheevr_2stage.f"
    iwork[1] = liwmin;

#line 775 "zheevr_2stage.f"
    return 0;

/*     End of ZHEEVR_2STAGE */

} /* zheevr_2stage__ */

