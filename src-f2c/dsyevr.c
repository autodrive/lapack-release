#line 1 "dsyevr.f"
/* dsyevr.f -- translated by f2c (version 20100827).
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

#line 1 "dsyevr.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;

/* > \brief <b> DSYEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYEVR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyevr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyevr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyevr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, */
/*                          IWORK, LIWORK, INFO ) */

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
/* > DSYEVR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A.  Eigenvalues and eigenvectors can be */
/* > selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
/* > */
/* > DSYEVR first reduces the matrix A to tridiagonal form T with a call */
/* > to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute */
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
/* > Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,26*N). */
/* >          For optimal efficiency, LWORK >= (NB+6)*N, */
/* >          where NB is the max of the blocksize for DSYTRD and DORMTR */
/* >          returned by ILAENV. */
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
/*  ===================================================================== */
/* Subroutine */ int dsyevr_(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nb, jj;
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
    static integer indwk;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer lwmin;
    static logical lower, wantz;
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
    extern /* Subroutine */ int dsytrd_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;


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

#line 386 "dsyevr.f"
    /* Parameter adjustments */
#line 386 "dsyevr.f"
    a_dim1 = *lda;
#line 386 "dsyevr.f"
    a_offset = 1 + a_dim1;
#line 386 "dsyevr.f"
    a -= a_offset;
#line 386 "dsyevr.f"
    --w;
#line 386 "dsyevr.f"
    z_dim1 = *ldz;
#line 386 "dsyevr.f"
    z_offset = 1 + z_dim1;
#line 386 "dsyevr.f"
    z__ -= z_offset;
#line 386 "dsyevr.f"
    --isuppz;
#line 386 "dsyevr.f"
    --work;
#line 386 "dsyevr.f"
    --iwork;
#line 386 "dsyevr.f"

#line 386 "dsyevr.f"
    /* Function Body */
#line 386 "dsyevr.f"
    ieeeok = ilaenv_(&c__10, "DSYEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 388 "dsyevr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 389 "dsyevr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 390 "dsyevr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 391 "dsyevr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 392 "dsyevr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 394 "dsyevr.f"
    lquery = *lwork == -1 || *liwork == -1;

/* Computing MAX */
#line 396 "dsyevr.f"
    i__1 = 1, i__2 = *n * 26;
#line 396 "dsyevr.f"
    lwmin = max(i__1,i__2);
/* Computing MAX */
#line 397 "dsyevr.f"
    i__1 = 1, i__2 = *n * 10;
#line 397 "dsyevr.f"
    liwmin = max(i__1,i__2);

#line 399 "dsyevr.f"
    *info = 0;
#line 400 "dsyevr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 401 "dsyevr.f"
	*info = -1;
#line 402 "dsyevr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 403 "dsyevr.f"
	*info = -2;
#line 404 "dsyevr.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 405 "dsyevr.f"
	*info = -3;
#line 406 "dsyevr.f"
    } else if (*n < 0) {
#line 407 "dsyevr.f"
	*info = -4;
#line 408 "dsyevr.f"
    } else if (*lda < max(1,*n)) {
#line 409 "dsyevr.f"
	*info = -6;
#line 410 "dsyevr.f"
    } else {
#line 411 "dsyevr.f"
	if (valeig) {
#line 412 "dsyevr.f"
	    if (*n > 0 && *vu <= *vl) {
#line 412 "dsyevr.f"
		*info = -8;
#line 412 "dsyevr.f"
	    }
#line 414 "dsyevr.f"
	} else if (indeig) {
#line 415 "dsyevr.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 416 "dsyevr.f"
		*info = -9;
#line 417 "dsyevr.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 418 "dsyevr.f"
		*info = -10;
#line 419 "dsyevr.f"
	    }
#line 420 "dsyevr.f"
	}
#line 421 "dsyevr.f"
    }
#line 422 "dsyevr.f"
    if (*info == 0) {
#line 423 "dsyevr.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 424 "dsyevr.f"
	    *info = -15;
#line 425 "dsyevr.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 426 "dsyevr.f"
	    *info = -18;
#line 427 "dsyevr.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 428 "dsyevr.f"
	    *info = -20;
#line 429 "dsyevr.f"
	}
#line 430 "dsyevr.f"
    }

#line 432 "dsyevr.f"
    if (*info == 0) {
#line 433 "dsyevr.f"
	nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 434 "dsyevr.f"
	i__1 = nb, i__2 = ilaenv_(&c__1, "DORMTR", uplo, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
#line 434 "dsyevr.f"
	nb = max(i__1,i__2);
/* Computing MAX */
#line 435 "dsyevr.f"
	i__1 = (nb + 1) * *n;
#line 435 "dsyevr.f"
	lwkopt = max(i__1,lwmin);
#line 436 "dsyevr.f"
	work[1] = (doublereal) lwkopt;
#line 437 "dsyevr.f"
	iwork[1] = liwmin;
#line 438 "dsyevr.f"
    }

#line 440 "dsyevr.f"
    if (*info != 0) {
#line 441 "dsyevr.f"
	i__1 = -(*info);
#line 441 "dsyevr.f"
	xerbla_("DSYEVR", &i__1, (ftnlen)6);
#line 442 "dsyevr.f"
	return 0;
#line 443 "dsyevr.f"
    } else if (lquery) {
#line 444 "dsyevr.f"
	return 0;
#line 445 "dsyevr.f"
    }

/*     Quick return if possible */

#line 449 "dsyevr.f"
    *m = 0;
#line 450 "dsyevr.f"
    if (*n == 0) {
#line 451 "dsyevr.f"
	work[1] = 1.;
#line 452 "dsyevr.f"
	return 0;
#line 453 "dsyevr.f"
    }

#line 455 "dsyevr.f"
    if (*n == 1) {
#line 456 "dsyevr.f"
	work[1] = 7.;
#line 457 "dsyevr.f"
	if (alleig || indeig) {
#line 458 "dsyevr.f"
	    *m = 1;
#line 459 "dsyevr.f"
	    w[1] = a[a_dim1 + 1];
#line 460 "dsyevr.f"
	} else {
#line 461 "dsyevr.f"
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
#line 462 "dsyevr.f"
		*m = 1;
#line 463 "dsyevr.f"
		w[1] = a[a_dim1 + 1];
#line 464 "dsyevr.f"
	    }
#line 465 "dsyevr.f"
	}
#line 466 "dsyevr.f"
	if (wantz) {
#line 467 "dsyevr.f"
	    z__[z_dim1 + 1] = 1.;
#line 468 "dsyevr.f"
	    isuppz[1] = 1;
#line 469 "dsyevr.f"
	    isuppz[2] = 1;
#line 470 "dsyevr.f"
	}
#line 471 "dsyevr.f"
	return 0;
#line 472 "dsyevr.f"
    }

/*     Get machine constants. */

#line 476 "dsyevr.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 477 "dsyevr.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 478 "dsyevr.f"
    smlnum = safmin / eps;
#line 479 "dsyevr.f"
    bignum = 1. / smlnum;
#line 480 "dsyevr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 481 "dsyevr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 481 "dsyevr.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 485 "dsyevr.f"
    iscale = 0;
#line 486 "dsyevr.f"
    abstll = *abstol;
#line 487 "dsyevr.f"
    if (valeig) {
#line 488 "dsyevr.f"
	vll = *vl;
#line 489 "dsyevr.f"
	vuu = *vu;
#line 490 "dsyevr.f"
    }
#line 491 "dsyevr.f"
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 492 "dsyevr.f"
    if (anrm > 0. && anrm < rmin) {
#line 493 "dsyevr.f"
	iscale = 1;
#line 494 "dsyevr.f"
	sigma = rmin / anrm;
#line 495 "dsyevr.f"
    } else if (anrm > rmax) {
#line 496 "dsyevr.f"
	iscale = 1;
#line 497 "dsyevr.f"
	sigma = rmax / anrm;
#line 498 "dsyevr.f"
    }
#line 499 "dsyevr.f"
    if (iscale == 1) {
#line 500 "dsyevr.f"
	if (lower) {
#line 501 "dsyevr.f"
	    i__1 = *n;
#line 501 "dsyevr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 502 "dsyevr.f"
		i__2 = *n - j + 1;
#line 502 "dsyevr.f"
		dscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 503 "dsyevr.f"
/* L10: */
#line 503 "dsyevr.f"
	    }
#line 504 "dsyevr.f"
	} else {
#line 505 "dsyevr.f"
	    i__1 = *n;
#line 505 "dsyevr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 506 "dsyevr.f"
		dscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 507 "dsyevr.f"
/* L20: */
#line 507 "dsyevr.f"
	    }
#line 508 "dsyevr.f"
	}
#line 509 "dsyevr.f"
	if (*abstol > 0.) {
#line 509 "dsyevr.f"
	    abstll = *abstol * sigma;
#line 509 "dsyevr.f"
	}
#line 511 "dsyevr.f"
	if (valeig) {
#line 512 "dsyevr.f"
	    vll = *vl * sigma;
#line 513 "dsyevr.f"
	    vuu = *vu * sigma;
#line 514 "dsyevr.f"
	}
#line 515 "dsyevr.f"
    }
/*     Initialize indices into workspaces.  Note: The IWORK indices are */
/*     used only if DSTERF or DSTEMR fail. */
/*     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the */
/*     elementary reflectors used in DSYTRD. */
#line 522 "dsyevr.f"
    indtau = 1;
/*     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries. */
#line 524 "dsyevr.f"
    indd = indtau + *n;
/*     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the */
/*     tridiagonal matrix from DSYTRD. */
#line 527 "dsyevr.f"
    inde = indd + *n;
/*     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over */
/*     -written by DSTEMR (the DSTERF path copies the diagonal to W). */
#line 530 "dsyevr.f"
    inddd = inde + *n;
/*     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over */
/*     -written while computing the eigenvalues in DSTERF and DSTEMR. */
#line 533 "dsyevr.f"
    indee = inddd + *n;
/*     INDWK is the starting offset of the left-over workspace, and */
/*     LLWORK is the remaining workspace size. */
#line 536 "dsyevr.f"
    indwk = indee + *n;
#line 537 "dsyevr.f"
    llwork = *lwork - indwk + 1;
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 541 "dsyevr.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 544 "dsyevr.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     DSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 549 "dsyevr.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 551 "dsyevr.f"
    indiwo = indifl + *n;

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

#line 556 "dsyevr.f"
    dsytrd_(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired */
/*     then call DSTERF or DSTEMR and DORMTR. */

#line 562 "dsyevr.f"
    if ((alleig || indeig && *il == 1 && *iu == *n) && ieeeok == 1) {
#line 564 "dsyevr.f"
	if (! wantz) {
#line 565 "dsyevr.f"
	    dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
#line 566 "dsyevr.f"
	    i__1 = *n - 1;
#line 566 "dsyevr.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 567 "dsyevr.f"
	    dsterf_(n, &w[1], &work[indee], info);
#line 568 "dsyevr.f"
	} else {
#line 569 "dsyevr.f"
	    i__1 = *n - 1;
#line 569 "dsyevr.f"
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
#line 570 "dsyevr.f"
	    dcopy_(n, &work[indd], &c__1, &work[inddd], &c__1);

#line 572 "dsyevr.f"
	    if (*abstol <= *n * 2. * eps) {
#line 573 "dsyevr.f"
		tryrac = TRUE_;
#line 574 "dsyevr.f"
	    } else {
#line 575 "dsyevr.f"
		tryrac = FALSE_;
#line 576 "dsyevr.f"
	    }
#line 577 "dsyevr.f"
	    dstemr_(jobz, "A", n, &work[inddd], &work[indee], vl, vu, il, iu, 
		    m, &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &
		    work[indwk], lwork, &iwork[1], liwork, info, (ftnlen)1, (
		    ftnlen)1);



/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEMR. */

#line 587 "dsyevr.f"
	    if (wantz && *info == 0) {
#line 588 "dsyevr.f"
		indwkn = inde;
#line 589 "dsyevr.f"
		llwrkn = *lwork - indwkn + 1;
#line 590 "dsyevr.f"
		dormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
			, &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 593 "dsyevr.f"
	    }
#line 594 "dsyevr.f"
	}


#line 597 "dsyevr.f"
	if (*info == 0) {
/*           Everything worked.  Skip DSTEBZ/DSTEIN.  IWORK(:) are */
/*           undefined. */
#line 600 "dsyevr.f"
	    *m = *n;
#line 601 "dsyevr.f"
	    goto L30;
#line 602 "dsyevr.f"
	}
#line 603 "dsyevr.f"
	*info = 0;
#line 604 "dsyevr.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN. */
/*     Also call DSTEBZ and DSTEIN if DSTEMR fails. */

#line 609 "dsyevr.f"
    if (wantz) {
#line 610 "dsyevr.f"
	*(unsigned char *)order = 'B';
#line 611 "dsyevr.f"
    } else {
#line 612 "dsyevr.f"
	*(unsigned char *)order = 'E';
#line 613 "dsyevr.f"
    }
#line 615 "dsyevr.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 620 "dsyevr.f"
    if (wantz) {
#line 621 "dsyevr.f"
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwk], &iwork[indiwo], &
		iwork[indifl], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

#line 629 "dsyevr.f"
	indwkn = inde;
#line 630 "dsyevr.f"
	llwrkn = *lwork - indwkn + 1;
#line 631 "dsyevr.f"
	dormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 633 "dsyevr.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

/*  Jump here if DSTEMR/DSTEIN succeeded. */
#line 638 "dsyevr.f"
L30:
#line 639 "dsyevr.f"
    if (iscale == 1) {
#line 640 "dsyevr.f"
	if (*info == 0) {
#line 641 "dsyevr.f"
	    imax = *m;
#line 642 "dsyevr.f"
	} else {
#line 643 "dsyevr.f"
	    imax = *info - 1;
#line 644 "dsyevr.f"
	}
#line 645 "dsyevr.f"
	d__1 = 1. / sigma;
#line 645 "dsyevr.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 646 "dsyevr.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK. */
/*     It may not be initialized (if DSTEMR/DSTEIN succeeded), and we do */
/*     not return this detailed information to the user. */

#line 653 "dsyevr.f"
    if (wantz) {
#line 654 "dsyevr.f"
	i__1 = *m - 1;
#line 654 "dsyevr.f"
	for (j = 1; j <= i__1; ++j) {
#line 655 "dsyevr.f"
	    i__ = 0;
#line 656 "dsyevr.f"
	    tmp1 = w[j];
#line 657 "dsyevr.f"
	    i__2 = *m;
#line 657 "dsyevr.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 658 "dsyevr.f"
		if (w[jj] < tmp1) {
#line 659 "dsyevr.f"
		    i__ = jj;
#line 660 "dsyevr.f"
		    tmp1 = w[jj];
#line 661 "dsyevr.f"
		}
#line 662 "dsyevr.f"
/* L40: */
#line 662 "dsyevr.f"
	    }

#line 664 "dsyevr.f"
	    if (i__ != 0) {
#line 665 "dsyevr.f"
		w[i__] = w[j];
#line 666 "dsyevr.f"
		w[j] = tmp1;
#line 667 "dsyevr.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 668 "dsyevr.f"
	    }
#line 669 "dsyevr.f"
/* L50: */
#line 669 "dsyevr.f"
	}
#line 670 "dsyevr.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 674 "dsyevr.f"
    work[1] = (doublereal) lwkopt;
#line 675 "dsyevr.f"
    iwork[1] = liwmin;

#line 677 "dsyevr.f"
    return 0;

/*     End of DSYEVR */

} /* dsyevr_ */

