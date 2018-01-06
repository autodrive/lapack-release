#line 1 "cheevr.f"
/* cheevr.f -- translated by f2c (version 20100827).
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

#line 1 "cheevr.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;

/* > \brief <b> CHEEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEEVR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheevr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheevr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheevr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/*                          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, */
/*                          RWORK, LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, */
/*      $                   M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEEVR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can */
/* > be selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
/* > */
/* > CHEEVR first reduces the matrix A to tridiagonal form T with a call */
/* > to CHETRD.  Then, whenever possible, CHEEVR calls CSTEMR to compute */
/* > the eigenspectrum using Relatively Robust Representations.  CSTEMR */
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
/* > Note 1 : CHEEVR calls CSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > CHEEVR calls SSTEBZ and CSTEIN on non-ieee machines and */
/* > when partial spectrum requests are made. */
/* > */
/* > Normal execution of CSTEMR may create NaNs and infinities and */
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
/* >          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and */
/* >          CSTEIN are called */
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
/* >          A is COMPLEX array, dimension (LDA, N) */
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
/* >          ISUPPZ( 2*i ). This is an output of CSTEMR (tridiagonal */
/* >          matrix). The support of the eigenvectors of A is typically */
/* >          1:N because of the unitary transformations applied by CUNMTR. */
/* >          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= max(1,2*N). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the max of the blocksize for CHETRD and for */
/* >          CUNMTR as returned by ILAENV. */
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
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
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

/* > \ingroup complexHEeigen */

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
/* Subroutine */ int cheevr_(char *jobz, char *range, char *uplo, integer *n, 
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
    static integer itmp1, indrd, indre;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    static integer indwk;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer lwmin;
    static logical lower;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantz, alleig, indeig;
    static integer iscale, ieeeok, indibl, indrdd, indifl, indree;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int chetrd_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen), csscal_(integer *,
	     doublereal *, doublecomplex *, integer *);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp;
    extern /* Subroutine */ int cstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static integer indiwo, indwkn;
    extern doublereal clansy_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int cstemr_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     integer *, doublereal *, doublecomplex *, integer *, integer *, 
	    integer *, logical *, doublereal *, integer *, integer *, integer 
	    *, integer *, ftnlen, ftnlen);
    static integer indrwk, liwmin;
    static logical tryrac;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer lrwmin, llwrkn, llwork, nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int cunmtr_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen, ftnlen), sstebz_(char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen, 
	    ftnlen);
    static logical lquery;
    static integer lwkopt, llrwork;


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

#line 412 "cheevr.f"
    /* Parameter adjustments */
#line 412 "cheevr.f"
    a_dim1 = *lda;
#line 412 "cheevr.f"
    a_offset = 1 + a_dim1;
#line 412 "cheevr.f"
    a -= a_offset;
#line 412 "cheevr.f"
    --w;
#line 412 "cheevr.f"
    z_dim1 = *ldz;
#line 412 "cheevr.f"
    z_offset = 1 + z_dim1;
#line 412 "cheevr.f"
    z__ -= z_offset;
#line 412 "cheevr.f"
    --isuppz;
#line 412 "cheevr.f"
    --work;
#line 412 "cheevr.f"
    --rwork;
#line 412 "cheevr.f"
    --iwork;
#line 412 "cheevr.f"

#line 412 "cheevr.f"
    /* Function Body */
#line 412 "cheevr.f"
    ieeeok = ilaenv_(&c__10, "CHEEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 414 "cheevr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 415 "cheevr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 416 "cheevr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 417 "cheevr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 418 "cheevr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 420 "cheevr.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

/* Computing MAX */
#line 423 "cheevr.f"
    i__1 = 1, i__2 = *n * 24;
#line 423 "cheevr.f"
    lrwmin = max(i__1,i__2);
/* Computing MAX */
#line 424 "cheevr.f"
    i__1 = 1, i__2 = *n * 10;
#line 424 "cheevr.f"
    liwmin = max(i__1,i__2);
/* Computing MAX */
#line 425 "cheevr.f"
    i__1 = 1, i__2 = *n << 1;
#line 425 "cheevr.f"
    lwmin = max(i__1,i__2);

#line 427 "cheevr.f"
    *info = 0;
#line 428 "cheevr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 429 "cheevr.f"
	*info = -1;
#line 430 "cheevr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 431 "cheevr.f"
	*info = -2;
#line 432 "cheevr.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 433 "cheevr.f"
	*info = -3;
#line 434 "cheevr.f"
    } else if (*n < 0) {
#line 435 "cheevr.f"
	*info = -4;
#line 436 "cheevr.f"
    } else if (*lda < max(1,*n)) {
#line 437 "cheevr.f"
	*info = -6;
#line 438 "cheevr.f"
    } else {
#line 439 "cheevr.f"
	if (valeig) {
#line 440 "cheevr.f"
	    if (*n > 0 && *vu <= *vl) {
#line 440 "cheevr.f"
		*info = -8;
#line 440 "cheevr.f"
	    }
#line 442 "cheevr.f"
	} else if (indeig) {
#line 443 "cheevr.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 444 "cheevr.f"
		*info = -9;
#line 445 "cheevr.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 446 "cheevr.f"
		*info = -10;
#line 447 "cheevr.f"
	    }
#line 448 "cheevr.f"
	}
#line 449 "cheevr.f"
    }
#line 450 "cheevr.f"
    if (*info == 0) {
#line 451 "cheevr.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 452 "cheevr.f"
	    *info = -15;
#line 453 "cheevr.f"
	}
#line 454 "cheevr.f"
    }

#line 456 "cheevr.f"
    if (*info == 0) {
#line 457 "cheevr.f"
	nb = ilaenv_(&c__1, "CHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 458 "cheevr.f"
	i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMTR", uplo, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
#line 458 "cheevr.f"
	nb = max(i__1,i__2);
/* Computing MAX */
#line 459 "cheevr.f"
	i__1 = (nb + 1) * *n;
#line 459 "cheevr.f"
	lwkopt = max(i__1,lwmin);
#line 460 "cheevr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 461 "cheevr.f"
	rwork[1] = (doublereal) lrwmin;
#line 462 "cheevr.f"
	iwork[1] = liwmin;

#line 464 "cheevr.f"
	if (*lwork < lwmin && ! lquery) {
#line 465 "cheevr.f"
	    *info = -18;
#line 466 "cheevr.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 467 "cheevr.f"
	    *info = -20;
#line 468 "cheevr.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 469 "cheevr.f"
	    *info = -22;
#line 470 "cheevr.f"
	}
#line 471 "cheevr.f"
    }

#line 473 "cheevr.f"
    if (*info != 0) {
#line 474 "cheevr.f"
	i__1 = -(*info);
#line 474 "cheevr.f"
	xerbla_("CHEEVR", &i__1, (ftnlen)6);
#line 475 "cheevr.f"
	return 0;
#line 476 "cheevr.f"
    } else if (lquery) {
#line 477 "cheevr.f"
	return 0;
#line 478 "cheevr.f"
    }

/*     Quick return if possible */

#line 482 "cheevr.f"
    *m = 0;
#line 483 "cheevr.f"
    if (*n == 0) {
#line 484 "cheevr.f"
	work[1].r = 1., work[1].i = 0.;
#line 485 "cheevr.f"
	return 0;
#line 486 "cheevr.f"
    }

#line 488 "cheevr.f"
    if (*n == 1) {
#line 489 "cheevr.f"
	work[1].r = 2., work[1].i = 0.;
#line 490 "cheevr.f"
	if (alleig || indeig) {
#line 491 "cheevr.f"
	    *m = 1;
#line 492 "cheevr.f"
	    i__1 = a_dim1 + 1;
#line 492 "cheevr.f"
	    w[1] = a[i__1].r;
#line 493 "cheevr.f"
	} else {
#line 494 "cheevr.f"
	    i__1 = a_dim1 + 1;
#line 494 "cheevr.f"
	    i__2 = a_dim1 + 1;
#line 494 "cheevr.f"
	    if (*vl < a[i__1].r && *vu >= a[i__2].r) {
#line 496 "cheevr.f"
		*m = 1;
#line 497 "cheevr.f"
		i__1 = a_dim1 + 1;
#line 497 "cheevr.f"
		w[1] = a[i__1].r;
#line 498 "cheevr.f"
	    }
#line 499 "cheevr.f"
	}
#line 500 "cheevr.f"
	if (wantz) {
#line 501 "cheevr.f"
	    i__1 = z_dim1 + 1;
#line 501 "cheevr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 502 "cheevr.f"
	    isuppz[1] = 1;
#line 503 "cheevr.f"
	    isuppz[2] = 1;
#line 504 "cheevr.f"
	}
#line 505 "cheevr.f"
	return 0;
#line 506 "cheevr.f"
    }

/*     Get machine constants. */

#line 510 "cheevr.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 511 "cheevr.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 512 "cheevr.f"
    smlnum = safmin / eps;
#line 513 "cheevr.f"
    bignum = 1. / smlnum;
#line 514 "cheevr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 515 "cheevr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 515 "cheevr.f"
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

#line 519 "cheevr.f"
    iscale = 0;
#line 520 "cheevr.f"
    abstll = *abstol;
#line 521 "cheevr.f"
    if (valeig) {
#line 522 "cheevr.f"
	vll = *vl;
#line 523 "cheevr.f"
	vuu = *vu;
#line 524 "cheevr.f"
    }
#line 525 "cheevr.f"
    anrm = clansy_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 526 "cheevr.f"
    if (anrm > 0. && anrm < rmin) {
#line 527 "cheevr.f"
	iscale = 1;
#line 528 "cheevr.f"
	sigma = rmin / anrm;
#line 529 "cheevr.f"
    } else if (anrm > rmax) {
#line 530 "cheevr.f"
	iscale = 1;
#line 531 "cheevr.f"
	sigma = rmax / anrm;
#line 532 "cheevr.f"
    }
#line 533 "cheevr.f"
    if (iscale == 1) {
#line 534 "cheevr.f"
	if (lower) {
#line 535 "cheevr.f"
	    i__1 = *n;
#line 535 "cheevr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 536 "cheevr.f"
		i__2 = *n - j + 1;
#line 536 "cheevr.f"
		csscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
#line 537 "cheevr.f"
/* L10: */
#line 537 "cheevr.f"
	    }
#line 538 "cheevr.f"
	} else {
#line 539 "cheevr.f"
	    i__1 = *n;
#line 539 "cheevr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 540 "cheevr.f"
		csscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
#line 541 "cheevr.f"
/* L20: */
#line 541 "cheevr.f"
	    }
#line 542 "cheevr.f"
	}
#line 543 "cheevr.f"
	if (*abstol > 0.) {
#line 543 "cheevr.f"
	    abstll = *abstol * sigma;
#line 543 "cheevr.f"
	}
#line 545 "cheevr.f"
	if (valeig) {
#line 546 "cheevr.f"
	    vll = *vl * sigma;
#line 547 "cheevr.f"
	    vuu = *vu * sigma;
#line 548 "cheevr.f"
	}
#line 549 "cheevr.f"
    }
/*     Initialize indices into workspaces.  Note: The IWORK indices are */
/*     used only if SSTERF or CSTEMR fail. */
/*     WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the */
/*     elementary reflectors used in CHETRD. */
#line 556 "cheevr.f"
    indtau = 1;
/*     INDWK is the starting offset of the remaining complex workspace, */
/*     and LLWORK is the remaining complex workspace size. */
#line 559 "cheevr.f"
    indwk = indtau + *n;
#line 560 "cheevr.f"
    llwork = *lwork - indwk + 1;
/*     RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal */
/*     entries. */
#line 564 "cheevr.f"
    indrd = 1;
/*     RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the */
/*     tridiagonal matrix from CHETRD. */
#line 567 "cheevr.f"
    indre = indrd + *n;
/*     RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over */
/*     -written by CSTEMR (the SSTERF path copies the diagonal to W). */
#line 570 "cheevr.f"
    indrdd = indre + *n;
/*     RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over */
/*     -written while computing the eigenvalues in SSTERF and CSTEMR. */
#line 573 "cheevr.f"
    indree = indrdd + *n;
/*     INDRWK is the starting offset of the left-over real workspace, and */
/*     LLRWORK is the remaining workspace size. */
#line 576 "cheevr.f"
    indrwk = indree + *n;
#line 577 "cheevr.f"
    llrwork = *lrwork - indrwk + 1;
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 581 "cheevr.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 584 "cheevr.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     SSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 589 "cheevr.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 591 "cheevr.f"
    indiwo = indifl + *n;

/*     Call CHETRD to reduce Hermitian matrix to tridiagonal form. */

#line 596 "cheevr.f"
    chetrd_(uplo, n, &a[a_offset], lda, &rwork[indrd], &rwork[indre], &work[
	    indtau], &work[indwk], &llwork, &iinfo, (ftnlen)1);

/*     If all eigenvalues are desired */
/*     then call SSTERF or CSTEMR and CUNMTR. */

#line 602 "cheevr.f"
    test = FALSE_;
#line 603 "cheevr.f"
    if (indeig) {
#line 604 "cheevr.f"
	if (*il == 1 && *iu == *n) {
#line 605 "cheevr.f"
	    test = TRUE_;
#line 606 "cheevr.f"
	}
#line 607 "cheevr.f"
    }
#line 608 "cheevr.f"
    if ((alleig || test) && ieeeok == 1) {
#line 609 "cheevr.f"
	if (! wantz) {
#line 610 "cheevr.f"
	    scopy_(n, &rwork[indrd], &c__1, &w[1], &c__1);
#line 611 "cheevr.f"
	    i__1 = *n - 1;
#line 611 "cheevr.f"
	    scopy_(&i__1, &rwork[indre], &c__1, &rwork[indree], &c__1);
#line 612 "cheevr.f"
	    ssterf_(n, &w[1], &rwork[indree], info);
#line 613 "cheevr.f"
	} else {
#line 614 "cheevr.f"
	    i__1 = *n - 1;
#line 614 "cheevr.f"
	    scopy_(&i__1, &rwork[indre], &c__1, &rwork[indree], &c__1);
#line 615 "cheevr.f"
	    scopy_(n, &rwork[indrd], &c__1, &rwork[indrdd], &c__1);

#line 617 "cheevr.f"
	    if (*abstol <= *n * 2. * eps) {
#line 618 "cheevr.f"
		tryrac = TRUE_;
#line 619 "cheevr.f"
	    } else {
#line 620 "cheevr.f"
		tryrac = FALSE_;
#line 621 "cheevr.f"
	    }
#line 622 "cheevr.f"
	    cstemr_(jobz, "A", n, &rwork[indrdd], &rwork[indree], vl, vu, il, 
		    iu, m, &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac,
		     &rwork[indrwk], &llrwork, &iwork[1], liwork, info, (
		    ftnlen)1, (ftnlen)1);

/*           Apply unitary matrix used in reduction to tridiagonal */
/*           form to eigenvectors returned by CSTEMR. */

#line 631 "cheevr.f"
	    if (wantz && *info == 0) {
#line 632 "cheevr.f"
		indwkn = indwk;
#line 633 "cheevr.f"
		llwrkn = *lwork - indwkn + 1;
#line 634 "cheevr.f"
		cunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
			, &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 637 "cheevr.f"
	    }
#line 638 "cheevr.f"
	}


#line 641 "cheevr.f"
	if (*info == 0) {
#line 642 "cheevr.f"
	    *m = *n;
#line 643 "cheevr.f"
	    goto L30;
#line 644 "cheevr.f"
	}
#line 645 "cheevr.f"
	*info = 0;
#line 646 "cheevr.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */
/*     Also call SSTEBZ and CSTEIN if CSTEMR fails. */

#line 651 "cheevr.f"
    if (wantz) {
#line 652 "cheevr.f"
	*(unsigned char *)order = 'B';
#line 653 "cheevr.f"
    } else {
#line 654 "cheevr.f"
	*(unsigned char *)order = 'E';
#line 655 "cheevr.f"
    }
#line 657 "cheevr.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indrd], &
	    rwork[indre], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

#line 662 "cheevr.f"
    if (wantz) {
#line 663 "cheevr.f"
	cstein_(n, &rwork[indrd], &rwork[indre], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwo], &iwork[indifl], info);

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by CSTEIN. */

#line 671 "cheevr.f"
	indwkn = indwk;
#line 672 "cheevr.f"
	llwrkn = *lwork - indwkn + 1;
#line 673 "cheevr.f"
	cunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 675 "cheevr.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 679 "cheevr.f"
L30:
#line 680 "cheevr.f"
    if (iscale == 1) {
#line 681 "cheevr.f"
	if (*info == 0) {
#line 682 "cheevr.f"
	    imax = *m;
#line 683 "cheevr.f"
	} else {
#line 684 "cheevr.f"
	    imax = *info - 1;
#line 685 "cheevr.f"
	}
#line 686 "cheevr.f"
	d__1 = 1. / sigma;
#line 686 "cheevr.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 687 "cheevr.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 692 "cheevr.f"
    if (wantz) {
#line 693 "cheevr.f"
	i__1 = *m - 1;
#line 693 "cheevr.f"
	for (j = 1; j <= i__1; ++j) {
#line 694 "cheevr.f"
	    i__ = 0;
#line 695 "cheevr.f"
	    tmp1 = w[j];
#line 696 "cheevr.f"
	    i__2 = *m;
#line 696 "cheevr.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 697 "cheevr.f"
		if (w[jj] < tmp1) {
#line 698 "cheevr.f"
		    i__ = jj;
#line 699 "cheevr.f"
		    tmp1 = w[jj];
#line 700 "cheevr.f"
		}
#line 701 "cheevr.f"
/* L40: */
#line 701 "cheevr.f"
	    }

#line 703 "cheevr.f"
	    if (i__ != 0) {
#line 704 "cheevr.f"
		itmp1 = iwork[indibl + i__ - 1];
#line 705 "cheevr.f"
		w[i__] = w[j];
#line 706 "cheevr.f"
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
#line 707 "cheevr.f"
		w[j] = tmp1;
#line 708 "cheevr.f"
		iwork[indibl + j - 1] = itmp1;
#line 709 "cheevr.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 710 "cheevr.f"
	    }
#line 711 "cheevr.f"
/* L50: */
#line 711 "cheevr.f"
	}
#line 712 "cheevr.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 716 "cheevr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 717 "cheevr.f"
    rwork[1] = (doublereal) lrwmin;
#line 718 "cheevr.f"
    iwork[1] = liwmin;

#line 720 "cheevr.f"
    return 0;

/*     End of CHEEVR */

} /* cheevr_ */

