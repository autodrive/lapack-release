#line 1 "zheevr.f"
/* zheevr.f -- translated by f2c (version 20100827).
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

#line 1 "zheevr.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;

/* > \brief <b> ZHEEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEVR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, */
/*                          RWORK, LRWORK, IWORK, LIWORK, INFO ) */

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
/* > ZHEEVR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can */
/* > be selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
/* > */
/* > ZHEEVR first reduces the matrix A to tridiagonal form T with a call */
/* > to ZHETRD.  Then, whenever possible, ZHEEVR calls ZSTEMR to compute */
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
/* > Note 1 : ZHEEVR calls ZSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > ZHEEVR calls DSTEBZ and ZSTEIN on non-ieee machines and */
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
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          If RANGE='V', the lower and upper bounds of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the indices (in ascending order) of the */
/* >          smallest and largest eigenvalues to be returned. */
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
/* >          ISUPPZ( 2*i ). */
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
/* >          The length of the array WORK.  LWORK >= max(1,2*N). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the max of the blocksize for ZHETRD and for */
/* >          ZUNMTR as returned by ILAENV. */
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

/* > \date September 2012 */

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
/*  ===================================================================== */
/* Subroutine */ int zheevr_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, integer *isuppz, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info, ftnlen jobz_len, ftnlen 
	range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nb, jj;
    static doublereal eps, vll, vuu, tmp1, anrm;
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
    static char order[1];
    static integer indwk;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lwmin;
    static logical lower, wantz;
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
    extern /* Subroutine */ int zhetrd_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static logical tryrac;
    static integer lrwmin, llwrkn, llwork, nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int zstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static logical lquery;
    static integer lwkopt;
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
    static integer llrwork;


/*  -- LAPACK driver routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 403 "zheevr.f"
    /* Parameter adjustments */
#line 403 "zheevr.f"
    a_dim1 = *lda;
#line 403 "zheevr.f"
    a_offset = 1 + a_dim1;
#line 403 "zheevr.f"
    a -= a_offset;
#line 403 "zheevr.f"
    --w;
#line 403 "zheevr.f"
    z_dim1 = *ldz;
#line 403 "zheevr.f"
    z_offset = 1 + z_dim1;
#line 403 "zheevr.f"
    z__ -= z_offset;
#line 403 "zheevr.f"
    --isuppz;
#line 403 "zheevr.f"
    --work;
#line 403 "zheevr.f"
    --rwork;
#line 403 "zheevr.f"
    --iwork;
#line 403 "zheevr.f"

#line 403 "zheevr.f"
    /* Function Body */
#line 403 "zheevr.f"
    ieeeok = ilaenv_(&c__10, "ZHEEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 405 "zheevr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 406 "zheevr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 407 "zheevr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 408 "zheevr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 409 "zheevr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 411 "zheevr.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

/* Computing MAX */
#line 414 "zheevr.f"
    i__1 = 1, i__2 = *n * 24;
#line 414 "zheevr.f"
    lrwmin = max(i__1,i__2);
/* Computing MAX */
#line 415 "zheevr.f"
    i__1 = 1, i__2 = *n * 10;
#line 415 "zheevr.f"
    liwmin = max(i__1,i__2);
/* Computing MAX */
#line 416 "zheevr.f"
    i__1 = 1, i__2 = *n << 1;
#line 416 "zheevr.f"
    lwmin = max(i__1,i__2);

#line 418 "zheevr.f"
    *info = 0;
#line 419 "zheevr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 420 "zheevr.f"
	*info = -1;
#line 421 "zheevr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 422 "zheevr.f"
	*info = -2;
#line 423 "zheevr.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 424 "zheevr.f"
	*info = -3;
#line 425 "zheevr.f"
    } else if (*n < 0) {
#line 426 "zheevr.f"
	*info = -4;
#line 427 "zheevr.f"
    } else if (*lda < max(1,*n)) {
#line 428 "zheevr.f"
	*info = -6;
#line 429 "zheevr.f"
    } else {
#line 430 "zheevr.f"
	if (valeig) {
#line 431 "zheevr.f"
	    if (*n > 0 && *vu <= *vl) {
#line 431 "zheevr.f"
		*info = -8;
#line 431 "zheevr.f"
	    }
#line 433 "zheevr.f"
	} else if (indeig) {
#line 434 "zheevr.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 435 "zheevr.f"
		*info = -9;
#line 436 "zheevr.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 437 "zheevr.f"
		*info = -10;
#line 438 "zheevr.f"
	    }
#line 439 "zheevr.f"
	}
#line 440 "zheevr.f"
    }
#line 441 "zheevr.f"
    if (*info == 0) {
#line 442 "zheevr.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 443 "zheevr.f"
	    *info = -15;
#line 444 "zheevr.f"
	}
#line 445 "zheevr.f"
    }

#line 447 "zheevr.f"
    if (*info == 0) {
#line 448 "zheevr.f"
	nb = ilaenv_(&c__1, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 449 "zheevr.f"
	i__1 = nb, i__2 = ilaenv_(&c__1, "ZUNMTR", uplo, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
#line 449 "zheevr.f"
	nb = max(i__1,i__2);
/* Computing MAX */
#line 450 "zheevr.f"
	i__1 = (nb + 1) * *n;
#line 450 "zheevr.f"
	lwkopt = max(i__1,lwmin);
#line 451 "zheevr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 452 "zheevr.f"
	rwork[1] = (doublereal) lrwmin;
#line 453 "zheevr.f"
	iwork[1] = liwmin;

#line 455 "zheevr.f"
	if (*lwork < lwmin && ! lquery) {
#line 456 "zheevr.f"
	    *info = -18;
#line 457 "zheevr.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 458 "zheevr.f"
	    *info = -20;
#line 459 "zheevr.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 460 "zheevr.f"
	    *info = -22;
#line 461 "zheevr.f"
	}
#line 462 "zheevr.f"
    }

#line 464 "zheevr.f"
    if (*info != 0) {
#line 465 "zheevr.f"
	i__1 = -(*info);
#line 465 "zheevr.f"
	xerbla_("ZHEEVR", &i__1, (ftnlen)6);
#line 466 "zheevr.f"
	return 0;
#line 467 "zheevr.f"
    } else if (lquery) {
#line 468 "zheevr.f"
	return 0;
#line 469 "zheevr.f"
    }

/*     Quick return if possible */

#line 473 "zheevr.f"
    *m = 0;
#line 474 "zheevr.f"
    if (*n == 0) {
#line 475 "zheevr.f"
	work[1].r = 1., work[1].i = 0.;
#line 476 "zheevr.f"
	return 0;
#line 477 "zheevr.f"
    }

#line 479 "zheevr.f"
    if (*n == 1) {
#line 480 "zheevr.f"
	work[1].r = 2., work[1].i = 0.;
#line 481 "zheevr.f"
	if (alleig || indeig) {
#line 482 "zheevr.f"
	    *m = 1;
#line 483 "zheevr.f"
	    i__1 = a_dim1 + 1;
#line 483 "zheevr.f"
	    w[1] = a[i__1].r;
#line 484 "zheevr.f"
	} else {
#line 485 "zheevr.f"
	    i__1 = a_dim1 + 1;
#line 485 "zheevr.f"
	    i__2 = a_dim1 + 1;
#line 485 "zheevr.f"
	    if (*vl < a[i__1].r && *vu >= a[i__2].r) {
#line 487 "zheevr.f"
		*m = 1;
#line 488 "zheevr.f"
		i__1 = a_dim1 + 1;
#line 488 "zheevr.f"
		w[1] = a[i__1].r;
#line 489 "zheevr.f"
	    }
#line 490 "zheevr.f"
	}
#line 491 "zheevr.f"
	if (wantz) {
#line 492 "zheevr.f"
	    i__1 = z_dim1 + 1;
#line 492 "zheevr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 493 "zheevr.f"
	    isuppz[1] = 1;
#line 494 "zheevr.f"
	    isuppz[2] = 1;
#line 495 "zheevr.f"
	}
#line 496 "zheevr.f"
	return 0;
#line 497 "zheevr.f"
    }

/*     Get machine constants. */

#line 501 "zheevr.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 502 "zheevr.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 503 "zheevr.f"
    smlnum = safmin / eps;
#line 504 "zheevr.f"
    bignum = 1. / smlnum;
#line 505 "zheevr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 506 "zheevr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 506 "zheevr.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 510 "zheevr.f"
    iscale = 0;
#line 511 "zheevr.f"
    abstll = *abstol;
#line 512 "zheevr.f"
    if (valeig) {
#line 513 "zheevr.f"
	vll = *vl;
#line 514 "zheevr.f"
	vuu = *vu;
#line 515 "zheevr.f"
    }
#line 516 "zheevr.f"
    anrm = zlansy_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 517 "zheevr.f"
    if (anrm > 0. && anrm < rmin) {
#line 518 "zheevr.f"
	iscale = 1;
#line 519 "zheevr.f"
	sigma = rmin / anrm;
#line 520 "zheevr.f"
    } else if (anrm > rmax) {
#line 521 "zheevr.f"
	iscale = 1;
#line 522 "zheevr.f"
	sigma = rmax / anrm;
#line 523 "zheevr.f"
    }
#line 524 "zheevr.f"
    if (iscale == 1) {
#line 525 "zheevr.f"
	if (lower) {
#line 526 "zheevr.f"
	    i__1 = *n;
#line 526 "zheevr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 527 "zheevr.f"
		i__2 = *n - j + 1;
#line 527 "zheevr.f"
		zdscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 528 "zheevr.f"
/* L10: */
#line 528 "zheevr.f"
	    }
#line 529 "zheevr.f"
	} else {
#line 530 "zheevr.f"
	    i__1 = *n;
#line 530 "zheevr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 531 "zheevr.f"
		zdscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 532 "zheevr.f"
/* L20: */
#line 532 "zheevr.f"
	    }
#line 533 "zheevr.f"
	}
#line 534 "zheevr.f"
	if (*abstol > 0.) {
#line 534 "zheevr.f"
	    abstll = *abstol * sigma;
#line 534 "zheevr.f"
	}
#line 536 "zheevr.f"
	if (valeig) {
#line 537 "zheevr.f"
	    vll = *vl * sigma;
#line 538 "zheevr.f"
	    vuu = *vu * sigma;
#line 539 "zheevr.f"
	}
#line 540 "zheevr.f"
    }
/*     Initialize indices into workspaces.  Note: The IWORK indices are */
/*     used only if DSTERF or ZSTEMR fail. */
/*     WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the */
/*     elementary reflectors used in ZHETRD. */
#line 547 "zheevr.f"
    indtau = 1;
/*     INDWK is the starting offset of the remaining complex workspace, */
/*     and LLWORK is the remaining complex workspace size. */
#line 550 "zheevr.f"
    indwk = indtau + *n;
#line 551 "zheevr.f"
    llwork = *lwork - indwk + 1;
/*     RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal */
/*     entries. */
#line 555 "zheevr.f"
    indrd = 1;
/*     RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the */
/*     tridiagonal matrix from ZHETRD. */
#line 558 "zheevr.f"
    indre = indrd + *n;
/*     RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over */
/*     -written by ZSTEMR (the DSTERF path copies the diagonal to W). */
#line 561 "zheevr.f"
    indrdd = indre + *n;
/*     RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over */
/*     -written while computing the eigenvalues in DSTERF and ZSTEMR. */
#line 564 "zheevr.f"
    indree = indrdd + *n;
/*     INDRWK is the starting offset of the left-over real workspace, and */
/*     LLRWORK is the remaining workspace size. */
#line 567 "zheevr.f"
    indrwk = indree + *n;
#line 568 "zheevr.f"
    llrwork = *lrwork - indrwk + 1;
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 572 "zheevr.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 575 "zheevr.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     DSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 580 "zheevr.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 582 "zheevr.f"
    indiwo = indifl + *n;

/*     Call ZHETRD to reduce Hermitian matrix to tridiagonal form. */

#line 587 "zheevr.f"
    zhetrd_(uplo, n, &a[a_offset], lda, &rwork[indrd], &rwork[indre], &work[
	    indtau], &work[indwk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired */
/*     then call DSTERF or ZSTEMR and ZUNMTR. */

#line 593 "zheevr.f"
    test = FALSE_;
#line 594 "zheevr.f"
    if (indeig) {
#line 595 "zheevr.f"
	if (*il == 1 && *iu == *n) {
#line 596 "zheevr.f"
	    test = TRUE_;
#line 597 "zheevr.f"
	}
#line 598 "zheevr.f"
    }
#line 599 "zheevr.f"
    if ((alleig || test) && ieeeok == 1) {
#line 600 "zheevr.f"
	if (! wantz) {
#line 601 "zheevr.f"
	    dcopy_(n, &rwork[indrd], &c__1, &w[1], &c__1);
#line 602 "zheevr.f"
	    i__1 = *n - 1;
#line 602 "zheevr.f"
	    dcopy_(&i__1, &rwork[indre], &c__1, &rwork[indree], &c__1);
#line 603 "zheevr.f"
	    dsterf_(n, &w[1], &rwork[indree], info);
#line 604 "zheevr.f"
	} else {
#line 605 "zheevr.f"
	    i__1 = *n - 1;
#line 605 "zheevr.f"
	    dcopy_(&i__1, &rwork[indre], &c__1, &rwork[indree], &c__1);
#line 606 "zheevr.f"
	    dcopy_(n, &rwork[indrd], &c__1, &rwork[indrdd], &c__1);

#line 608 "zheevr.f"
	    if (*abstol <= *n * 2. * eps) {
#line 609 "zheevr.f"
		tryrac = TRUE_;
#line 610 "zheevr.f"
	    } else {
#line 611 "zheevr.f"
		tryrac = FALSE_;
#line 612 "zheevr.f"
	    }
#line 613 "zheevr.f"
	    zstemr_(jobz, "A", n, &rwork[indrdd], &rwork[indree], vl, vu, il, 
		    iu, m, &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac,
		     &rwork[indrwk], &llrwork, &iwork[1], liwork, info, (
		    ftnlen)1, (ftnlen)1);

/*           Apply unitary matrix used in reduction to tridiagonal */
/*           form to eigenvectors returned by ZSTEIN. */

#line 622 "zheevr.f"
	    if (wantz && *info == 0) {
#line 623 "zheevr.f"
		indwkn = indwk;
#line 624 "zheevr.f"
		llwrkn = *lwork - indwkn + 1;
#line 625 "zheevr.f"
		zunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
			, &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 628 "zheevr.f"
	    }
#line 629 "zheevr.f"
	}


#line 632 "zheevr.f"
	if (*info == 0) {
#line 633 "zheevr.f"
	    *m = *n;
#line 634 "zheevr.f"
	    goto L30;
#line 635 "zheevr.f"
	}
#line 636 "zheevr.f"
	*info = 0;
#line 637 "zheevr.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN. */
/*     Also call DSTEBZ and ZSTEIN if ZSTEMR fails. */

#line 642 "zheevr.f"
    if (wantz) {
#line 643 "zheevr.f"
	*(unsigned char *)order = 'B';
#line 644 "zheevr.f"
    } else {
#line 645 "zheevr.f"
	*(unsigned char *)order = 'E';
#line 646 "zheevr.f"
    }
#line 648 "zheevr.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indrd], &
	    rwork[indre], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 653 "zheevr.f"
    if (wantz) {
#line 654 "zheevr.f"
	zstein_(n, &rwork[indrd], &rwork[indre], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwo], &iwork[indifl], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by ZSTEIN. */

#line 662 "zheevr.f"
	indwkn = indwk;
#line 663 "zheevr.f"
	llwrkn = *lwork - indwkn + 1;
#line 664 "zheevr.f"
	zunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 666 "zheevr.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 670 "zheevr.f"
L30:
#line 671 "zheevr.f"
    if (iscale == 1) {
#line 672 "zheevr.f"
	if (*info == 0) {
#line 673 "zheevr.f"
	    imax = *m;
#line 674 "zheevr.f"
	} else {
#line 675 "zheevr.f"
	    imax = *info - 1;
#line 676 "zheevr.f"
	}
#line 677 "zheevr.f"
	d__1 = 1. / sigma;
#line 677 "zheevr.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 678 "zheevr.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 683 "zheevr.f"
    if (wantz) {
#line 684 "zheevr.f"
	i__1 = *m - 1;
#line 684 "zheevr.f"
	for (j = 1; j <= i__1; ++j) {
#line 685 "zheevr.f"
	    i__ = 0;
#line 686 "zheevr.f"
	    tmp1 = w[j];
#line 687 "zheevr.f"
	    i__2 = *m;
#line 687 "zheevr.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 688 "zheevr.f"
		if (w[jj] < tmp1) {
#line 689 "zheevr.f"
		    i__ = jj;
#line 690 "zheevr.f"
		    tmp1 = w[jj];
#line 691 "zheevr.f"
		}
#line 692 "zheevr.f"
/* L40: */
#line 692 "zheevr.f"
	    }

#line 694 "zheevr.f"
	    if (i__ != 0) {
#line 695 "zheevr.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 696 "zheevr.f"
		w[i__] = w[j];
#line 697 "zheevr.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 698 "zheevr.f"
		w[j] = tmp1;
#line 699 "zheevr.f"
		iwork[indibl + j - 1] = itmp1;
#line 700 "zheevr.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 701 "zheevr.f"
	    }
#line 702 "zheevr.f"
/* L50: */
#line 702 "zheevr.f"
	}
#line 703 "zheevr.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 707 "zheevr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 708 "zheevr.f"
    rwork[1] = (doublereal) lrwmin;
#line 709 "zheevr.f"
    iwork[1] = liwmin;

#line 711 "zheevr.f"
    return 0;

/*     End of ZHEEVR */

} /* zheevr_ */

