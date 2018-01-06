#line 1 "cstemr.f"
/* cstemr.f -- translated by f2c (version 20100827).
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

#line 1 "cstemr.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = .003;

/* > \brief \b CSTEMR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSTEMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cstemr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cstemr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cstemr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, */
/*                          M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, */
/*                          IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE */
/*       LOGICAL            TRYRAC */
/*       INTEGER            IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N */
/*       REAL             VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       REAL               D( * ), E( * ), W( * ), WORK( * ) */
/*       COMPLEX            Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSTEMR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric tridiagonal matrix T. Any such unreduced matrix has */
/* > a well defined set of pairwise different real eigenvalues, the corresponding */
/* > real eigenvectors are pairwise orthogonal. */
/* > */
/* > The spectrum may be computed either completely or partially by specifying */
/* > either an interval (VL,VU] or a range of indices IL:IU for the desired */
/* > eigenvalues. */
/* > */
/* > Depending on the number of desired eigenvalues, these are computed either */
/* > by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are */
/* > computed by the use of various suitable L D L^T factorizations near clusters */
/* > of close eigenvalues (referred to as RRRs, Relatively Robust */
/* > Representations). An informal sketch of the algorithm follows. */
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
/* > For more details, see: */
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
/* > Further Details */
/* > 1.CSTEMR works only on machines which follow IEEE-754 */
/* > floating-point standard in their handling of infinities and NaNs. */
/* > This permits the use of efficient inner loops avoiding a check for */
/* > zero divisors. */
/* > */
/* > 2. LAPACK routines can be used to reduce a complex Hermitean matrix to */
/* > real symmetric tridiagonal form. */
/* > */
/* > (Any complex Hermitean tridiagonal matrix has real values on its diagonal */
/* > and potentially complex numbers on its off-diagonals. By applying a */
/* > similarity transform with an appropriate diagonal matrix */
/* > diag(1,e^{i \phy_1}, ... , e^{i \phy_{n-1}}), the complex Hermitean */
/* > matrix can be transformed into a real symmetric matrix and complex */
/* > arithmetic can be entirely avoided.) */
/* > */
/* > While the eigenvectors of the real symmetric tridiagonal matrix are real, */
/* > the eigenvectors of original complex Hermitean matrix have complex entries */
/* > in general. */
/* > Since LAPACK drivers overwrite the matrix data with the eigenvectors, */
/* > CSTEMR accepts complex workspace to facilitate interoperability */
/* > with CUNMTR or CUPMTR. */
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
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the N diagonal elements of the tridiagonal matrix */
/* >          T. On exit, D is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* >          On entry, the (N-1) subdiagonal elements of the tridiagonal */
/* >          matrix T in elements 1 to N-1 of E. E(N) need not be set on */
/* >          input, but is used internally as workspace. */
/* >          On exit, E is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* > */
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
/* > */
/* >          If RANGE='I', the indices (in ascending order) of the */
/* >          smallest and largest eigenvalues to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
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
/* >          Z is COMPLEX array, dimension (LDZ, max(1,M) ) */
/* >          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix T */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* >          Note: the user must ensure that at least max(1,M) columns are */
/* >          supplied in the array Z; if RANGE = 'V', the exact value of M */
/* >          is not known in advance and can be computed with a workspace */
/* >          query by setting NZC = -1, see below. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', then LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] NZC */
/* > \verbatim */
/* >          NZC is INTEGER */
/* >          The number of eigenvectors to be held in the array Z. */
/* >          If RANGE = 'A', then NZC >= max(1,N). */
/* >          If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU]. */
/* >          If RANGE = 'I', then NZC >= IU-IL+1. */
/* >          If NZC = -1, then a workspace query is assumed; the */
/* >          routine calculates the number of columns of the array Z that */
/* >          are needed to hold the eigenvectors. */
/* >          This value is returned as the first entry of the Z array, and */
/* >          no error message related to NZC is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* >          ISUPPZ is INTEGER ARRAY, dimension ( 2*max(1,M) ) */
/* >          The support of the eigenvectors in Z, i.e., the indices */
/* >          indicating the nonzero elements in Z. The i-th computed eigenvector */
/* >          is nonzero only in elements ISUPPZ( 2*i-1 ) through */
/* >          ISUPPZ( 2*i ). This is relevant in the case when the matrix */
/* >          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] TRYRAC */
/* > \verbatim */
/* >          TRYRAC is LOGICAL */
/* >          If TRYRAC.EQ..TRUE., indicates that the code should check whether */
/* >          the tridiagonal matrix defines its eigenvalues to high relative */
/* >          accuracy.  If so, the code uses relative-accuracy preserving */
/* >          algorithms that might be (a bit) slower depending on the matrix. */
/* >          If the matrix does not define its eigenvalues to high relative */
/* >          accuracy, the code can uses possibly faster algorithms. */
/* >          If TRYRAC.EQ..FALSE., the code is not required to guarantee */
/* >          relatively accurate eigenvalues and can use the fastest possible */
/* >          techniques. */
/* >          On exit, a .TRUE. TRYRAC will be set to .FALSE. if the matrix */
/* >          does not define its eigenvalues to high relative accuracy. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (LWORK) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal */
/* >          (and minimal) LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= max(1,18*N) */
/* >          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'. */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (LIWORK) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK.  LIWORK >= max(1,10*N) */
/* >          if the eigenvectors are desired, and LIWORK >= max(1,8*N) */
/* >          if only the eigenvalues are to be computed. */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of the IWORK array, */
/* >          returns this value as the first entry of the IWORK array, and */
/* >          no error message related to LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          On exit, INFO */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = 1X, internal error in SLARRE, */
/* >                if INFO = 2X, internal error in CLARRV. */
/* >                Here, the digit X = ABS( IINFO ) < 10, where IINFO is */
/* >                the nonzero error code returned by SLARRE or */
/* >                CLARRV, respectively. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int cstemr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, integer *m, doublereal *w, doublecomplex *z__, integer *
	ldz, integer *nzc, integer *isuppz, logical *tryrac, doublereal *work,
	 integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen range_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal r1, r2;
    static integer jj;
    static doublereal cs;
    static integer in;
    static doublereal sn, wl, wu;
    static integer iil, iiu;
    static doublereal eps, tmp;
    static integer indd, iend, jblk, wend;
    static doublereal rmin, rmax;
    static integer itmp;
    static doublereal tnrm;
    static integer inde2;
    extern /* Subroutine */ int slae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static integer itmp2;
    static doublereal rtol1, rtol2, scale;
    static integer indgp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iindw, ilast;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer lwmin;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantz;
    extern /* Subroutine */ int slaev2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static logical alleig;
    static integer ibegin;
    static logical indeig;
    static integer iindbl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    static integer wbegin;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer inderr, iindwk, indgrs, offset;
    extern /* Subroutine */ int slarrc_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, ftnlen), clarrv_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublecomplex *, integer *, integer *, doublereal *, integer *, 
	    integer *), slarre_(char *, integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer iinspl, indwrk, ifirst, liwmin, nzcmin;
    static doublereal pivmin, thresh;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int slarrj_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *);
    static integer nsplit;
    extern /* Subroutine */ int slarrr_(integer *, doublereal *, doublereal *,
	     integer *);
    static doublereal smlnum;
    extern /* Subroutine */ int slasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static logical lquery, zquery;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 387 "cstemr.f"
    /* Parameter adjustments */
#line 387 "cstemr.f"
    --d__;
#line 387 "cstemr.f"
    --e;
#line 387 "cstemr.f"
    --w;
#line 387 "cstemr.f"
    z_dim1 = *ldz;
#line 387 "cstemr.f"
    z_offset = 1 + z_dim1;
#line 387 "cstemr.f"
    z__ -= z_offset;
#line 387 "cstemr.f"
    --isuppz;
#line 387 "cstemr.f"
    --work;
#line 387 "cstemr.f"
    --iwork;
#line 387 "cstemr.f"

#line 387 "cstemr.f"
    /* Function Body */
#line 387 "cstemr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 388 "cstemr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 389 "cstemr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 390 "cstemr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 392 "cstemr.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 393 "cstemr.f"
    zquery = *nzc == -1;
/*     SSTEMR needs WORK of size 6*N, IWORK of size 3*N. */
/*     In addition, SLARRE needs WORK of size 6*N, IWORK of size 5*N. */
/*     Furthermore, CLARRV needs WORK of size 12*N, IWORK of size 7*N. */
#line 398 "cstemr.f"
    if (wantz) {
#line 399 "cstemr.f"
	lwmin = *n * 18;
#line 400 "cstemr.f"
	liwmin = *n * 10;
#line 401 "cstemr.f"
    } else {
/*        need less workspace if only the eigenvalues are wanted */
#line 403 "cstemr.f"
	lwmin = *n * 12;
#line 404 "cstemr.f"
	liwmin = *n << 3;
#line 405 "cstemr.f"
    }
#line 407 "cstemr.f"
    wl = 0.;
#line 408 "cstemr.f"
    wu = 0.;
#line 409 "cstemr.f"
    iil = 0;
#line 410 "cstemr.f"
    iiu = 0;
#line 411 "cstemr.f"
    nsplit = 0;
#line 413 "cstemr.f"
    if (valeig) {
/*        We do not reference VL, VU in the cases RANGE = 'I','A' */
/*        The interval (WL, WU] contains all the wanted eigenvalues. */
/*        It is either given by the user or computed in SLARRE. */
#line 417 "cstemr.f"
	wl = *vl;
#line 418 "cstemr.f"
	wu = *vu;
#line 419 "cstemr.f"
    } else if (indeig) {
/*        We do not reference IL, IU in the cases RANGE = 'V','A' */
#line 421 "cstemr.f"
	iil = *il;
#line 422 "cstemr.f"
	iiu = *iu;
#line 423 "cstemr.f"
    }

#line 425 "cstemr.f"
    *info = 0;
#line 426 "cstemr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 427 "cstemr.f"
	*info = -1;
#line 428 "cstemr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 429 "cstemr.f"
	*info = -2;
#line 430 "cstemr.f"
    } else if (*n < 0) {
#line 431 "cstemr.f"
	*info = -3;
#line 432 "cstemr.f"
    } else if (valeig && *n > 0 && wu <= wl) {
#line 433 "cstemr.f"
	*info = -7;
#line 434 "cstemr.f"
    } else if (indeig && (iil < 1 || iil > *n)) {
#line 435 "cstemr.f"
	*info = -8;
#line 436 "cstemr.f"
    } else if (indeig && (iiu < iil || iiu > *n)) {
#line 437 "cstemr.f"
	*info = -9;
#line 438 "cstemr.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 439 "cstemr.f"
	*info = -13;
#line 440 "cstemr.f"
    } else if (*lwork < lwmin && ! lquery) {
#line 441 "cstemr.f"
	*info = -17;
#line 442 "cstemr.f"
    } else if (*liwork < liwmin && ! lquery) {
#line 443 "cstemr.f"
	*info = -19;
#line 444 "cstemr.f"
    }

/*     Get machine constants. */

#line 448 "cstemr.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 449 "cstemr.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 450 "cstemr.f"
    smlnum = safmin / eps;
#line 451 "cstemr.f"
    bignum = 1. / smlnum;
#line 452 "cstemr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 453 "cstemr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 453 "cstemr.f"
    rmax = min(d__1,d__2);

#line 455 "cstemr.f"
    if (*info == 0) {
#line 456 "cstemr.f"
	work[1] = (doublereal) lwmin;
#line 457 "cstemr.f"
	iwork[1] = liwmin;

#line 459 "cstemr.f"
	if (wantz && alleig) {
#line 460 "cstemr.f"
	    nzcmin = *n;
#line 461 "cstemr.f"
	} else if (wantz && valeig) {
#line 462 "cstemr.f"
	    slarrc_("T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp, &
		    itmp2, info, (ftnlen)1);
#line 464 "cstemr.f"
	} else if (wantz && indeig) {
#line 465 "cstemr.f"
	    nzcmin = iiu - iil + 1;
#line 466 "cstemr.f"
	} else {
/*           WANTZ .EQ. FALSE. */
#line 468 "cstemr.f"
	    nzcmin = 0;
#line 469 "cstemr.f"
	}
#line 470 "cstemr.f"
	if (zquery && *info == 0) {
#line 471 "cstemr.f"
	    i__1 = z_dim1 + 1;
#line 471 "cstemr.f"
	    z__[i__1].r = (doublereal) nzcmin, z__[i__1].i = 0.;
#line 472 "cstemr.f"
	} else if (*nzc < nzcmin && ! zquery) {
#line 473 "cstemr.f"
	    *info = -14;
#line 474 "cstemr.f"
	}
#line 475 "cstemr.f"
    }
#line 477 "cstemr.f"
    if (*info != 0) {

#line 479 "cstemr.f"
	i__1 = -(*info);
#line 479 "cstemr.f"
	xerbla_("CSTEMR", &i__1, (ftnlen)6);

#line 481 "cstemr.f"
	return 0;
#line 482 "cstemr.f"
    } else if (lquery || zquery) {
#line 483 "cstemr.f"
	return 0;
#line 484 "cstemr.f"
    }

/*     Handle N = 0, 1, and 2 cases immediately */

#line 488 "cstemr.f"
    *m = 0;
#line 489 "cstemr.f"
    if (*n == 0) {
#line 489 "cstemr.f"
	return 0;
#line 489 "cstemr.f"
    }

#line 492 "cstemr.f"
    if (*n == 1) {
#line 493 "cstemr.f"
	if (alleig || indeig) {
#line 494 "cstemr.f"
	    *m = 1;
#line 495 "cstemr.f"
	    w[1] = d__[1];
#line 496 "cstemr.f"
	} else {
#line 497 "cstemr.f"
	    if (wl < d__[1] && wu >= d__[1]) {
#line 498 "cstemr.f"
		*m = 1;
#line 499 "cstemr.f"
		w[1] = d__[1];
#line 500 "cstemr.f"
	    }
#line 501 "cstemr.f"
	}
#line 502 "cstemr.f"
	if (wantz && ! zquery) {
#line 503 "cstemr.f"
	    i__1 = z_dim1 + 1;
#line 503 "cstemr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 504 "cstemr.f"
	    isuppz[1] = 1;
#line 505 "cstemr.f"
	    isuppz[2] = 1;
#line 506 "cstemr.f"
	}
#line 507 "cstemr.f"
	return 0;
#line 508 "cstemr.f"
    }

#line 510 "cstemr.f"
    if (*n == 2) {
#line 511 "cstemr.f"
	if (! wantz) {
#line 512 "cstemr.f"
	    slae2_(&d__[1], &e[1], &d__[2], &r1, &r2);
#line 513 "cstemr.f"
	} else if (wantz && ! zquery) {
#line 514 "cstemr.f"
	    slaev2_(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
#line 515 "cstemr.f"
	}
#line 516 "cstemr.f"
	if (alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1) {
#line 520 "cstemr.f"
	    ++(*m);
#line 521 "cstemr.f"
	    w[*m] = r2;
#line 522 "cstemr.f"
	    if (wantz && ! zquery) {
#line 523 "cstemr.f"
		i__1 = *m * z_dim1 + 1;
#line 523 "cstemr.f"
		d__1 = -sn;
#line 523 "cstemr.f"
		z__[i__1].r = d__1, z__[i__1].i = 0.;
#line 524 "cstemr.f"
		i__1 = *m * z_dim1 + 2;
#line 524 "cstemr.f"
		z__[i__1].r = cs, z__[i__1].i = 0.;
/*              Note: At most one of SN and CS can be zero. */
#line 526 "cstemr.f"
		if (sn != 0.) {
#line 527 "cstemr.f"
		    if (cs != 0.) {
#line 528 "cstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 529 "cstemr.f"
			isuppz[*m * 2] = 2;
#line 530 "cstemr.f"
		    } else {
#line 531 "cstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 532 "cstemr.f"
			isuppz[*m * 2] = 1;
#line 533 "cstemr.f"
		    }
#line 534 "cstemr.f"
		} else {
#line 535 "cstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 536 "cstemr.f"
		    isuppz[*m * 2] = 2;
#line 537 "cstemr.f"
		}
#line 538 "cstemr.f"
	    }
#line 539 "cstemr.f"
	}
#line 540 "cstemr.f"
	if (alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2) {
#line 544 "cstemr.f"
	    ++(*m);
#line 545 "cstemr.f"
	    w[*m] = r1;
#line 546 "cstemr.f"
	    if (wantz && ! zquery) {
#line 547 "cstemr.f"
		i__1 = *m * z_dim1 + 1;
#line 547 "cstemr.f"
		z__[i__1].r = cs, z__[i__1].i = 0.;
#line 548 "cstemr.f"
		i__1 = *m * z_dim1 + 2;
#line 548 "cstemr.f"
		z__[i__1].r = sn, z__[i__1].i = 0.;
/*              Note: At most one of SN and CS can be zero. */
#line 550 "cstemr.f"
		if (sn != 0.) {
#line 551 "cstemr.f"
		    if (cs != 0.) {
#line 552 "cstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 553 "cstemr.f"
			isuppz[*m * 2] = 2;
#line 554 "cstemr.f"
		    } else {
#line 555 "cstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 556 "cstemr.f"
			isuppz[*m * 2] = 1;
#line 557 "cstemr.f"
		    }
#line 558 "cstemr.f"
		} else {
#line 559 "cstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 560 "cstemr.f"
		    isuppz[*m * 2] = 2;
#line 561 "cstemr.f"
		}
#line 562 "cstemr.f"
	    }
#line 563 "cstemr.f"
	}
#line 564 "cstemr.f"
    } else {
/*        Continue with general N */
#line 568 "cstemr.f"
	indgrs = 1;
#line 569 "cstemr.f"
	inderr = (*n << 1) + 1;
#line 570 "cstemr.f"
	indgp = *n * 3 + 1;
#line 571 "cstemr.f"
	indd = (*n << 2) + 1;
#line 572 "cstemr.f"
	inde2 = *n * 5 + 1;
#line 573 "cstemr.f"
	indwrk = *n * 6 + 1;

#line 575 "cstemr.f"
	iinspl = 1;
#line 576 "cstemr.f"
	iindbl = *n + 1;
#line 577 "cstemr.f"
	iindw = (*n << 1) + 1;
#line 578 "cstemr.f"
	iindwk = *n * 3 + 1;

/*        Scale matrix to allowable range, if necessary. */
/*        The allowable range is related to the PIVMIN parameter; see the */
/*        comments in SLARRD.  The preference for scaling small values */
/*        up is heuristic; we expect users' matrices not to be close to the */
/*        RMAX threshold. */

#line 586 "cstemr.f"
	scale = 1.;
#line 587 "cstemr.f"
	tnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 588 "cstemr.f"
	if (tnrm > 0. && tnrm < rmin) {
#line 589 "cstemr.f"
	    scale = rmin / tnrm;
#line 590 "cstemr.f"
	} else if (tnrm > rmax) {
#line 591 "cstemr.f"
	    scale = rmax / tnrm;
#line 592 "cstemr.f"
	}
#line 593 "cstemr.f"
	if (scale != 1.) {
#line 594 "cstemr.f"
	    sscal_(n, &scale, &d__[1], &c__1);
#line 595 "cstemr.f"
	    i__1 = *n - 1;
#line 595 "cstemr.f"
	    sscal_(&i__1, &scale, &e[1], &c__1);
#line 596 "cstemr.f"
	    tnrm *= scale;
#line 597 "cstemr.f"
	    if (valeig) {
/*              If eigenvalues in interval have to be found, */
/*              scale (WL, WU] accordingly */
#line 600 "cstemr.f"
		wl *= scale;
#line 601 "cstemr.f"
		wu *= scale;
#line 602 "cstemr.f"
	    }
#line 603 "cstemr.f"
	}

/*        Compute the desired eigenvalues of the tridiagonal after splitting */
/*        into smaller subblocks if the corresponding off-diagonal elements */
/*        are small */
/*        THRESH is the splitting parameter for SLARRE */
/*        A negative THRESH forces the old splitting criterion based on the */
/*        size of the off-diagonal. A positive THRESH switches to splitting */
/*        which preserves relative accuracy. */

#line 613 "cstemr.f"
	if (*tryrac) {
/*           Test whether the matrix warrants the more expensive relative approach. */
#line 615 "cstemr.f"
	    slarrr_(n, &d__[1], &e[1], &iinfo);
#line 616 "cstemr.f"
	} else {
/*           The user does not care about relative accurately eigenvalues */
#line 618 "cstemr.f"
	    iinfo = -1;
#line 619 "cstemr.f"
	}
/*        Set the splitting criterion */
#line 621 "cstemr.f"
	if (iinfo == 0) {
#line 622 "cstemr.f"
	    thresh = eps;
#line 623 "cstemr.f"
	} else {
#line 624 "cstemr.f"
	    thresh = -eps;
/*           relative accuracy is desired but T does not guarantee it */
#line 626 "cstemr.f"
	    *tryrac = FALSE_;
#line 627 "cstemr.f"
	}

#line 629 "cstemr.f"
	if (*tryrac) {
/*           Copy original diagonal, needed to guarantee relative accuracy */
#line 631 "cstemr.f"
	    scopy_(n, &d__[1], &c__1, &work[indd], &c__1);
#line 632 "cstemr.f"
	}
/*        Store the squares of the offdiagonal values of T */
#line 634 "cstemr.f"
	i__1 = *n - 1;
#line 634 "cstemr.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
#line 635 "cstemr.f"
	    d__1 = e[j];
#line 635 "cstemr.f"
	    work[inde2 + j - 1] = d__1 * d__1;
#line 636 "cstemr.f"
/* L5: */
#line 636 "cstemr.f"
	}
/*        Set the tolerance parameters for bisection */
#line 639 "cstemr.f"
	if (! wantz) {
/*           SLARRE computes the eigenvalues to full precision. */
#line 641 "cstemr.f"
	    rtol1 = eps * 4.;
#line 642 "cstemr.f"
	    rtol2 = eps * 4.;
#line 643 "cstemr.f"
	} else {
/*           SLARRE computes the eigenvalues to less than full precision. */
/*           CLARRV will refine the eigenvalue approximations, and we only */
/*           need less accurate initial bisection in SLARRE. */
/*           Note: these settings do only affect the subset case and SLARRE */
/* Computing MAX */
#line 648 "cstemr.f"
	    d__1 = sqrt(eps) * .05, d__2 = eps * 4.;
#line 648 "cstemr.f"
	    rtol1 = max(d__1,d__2);
/* Computing MAX */
#line 649 "cstemr.f"
	    d__1 = sqrt(eps) * .005, d__2 = eps * 4.;
#line 649 "cstemr.f"
	    rtol2 = max(d__1,d__2);
#line 650 "cstemr.f"
	}
#line 651 "cstemr.f"
	slarre_(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2], 
		&rtol1, &rtol2, &thresh, &nsplit, &iwork[iinspl], m, &w[1], &
		work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw], &
		work[indgrs], &pivmin, &work[indwrk], &iwork[iindwk], &iinfo, 
		(ftnlen)1);
#line 657 "cstemr.f"
	if (iinfo != 0) {
#line 658 "cstemr.f"
	    *info = abs(iinfo) + 10;
#line 659 "cstemr.f"
	    return 0;
#line 660 "cstemr.f"
	}
/*        Note that if RANGE .NE. 'V', SLARRE computes bounds on the desired */
/*        part of the spectrum. All desired eigenvalues are contained in */
/*        (WL,WU] */
#line 666 "cstemr.f"
	if (wantz) {

/*           Compute the desired eigenvectors corresponding to the computed */
/*           eigenvalues */

#line 671 "cstemr.f"
	    clarrv_(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &
		    c__1, m, &c_b18, &rtol1, &rtol2, &w[1], &work[inderr], &
		    work[indgp], &iwork[iindbl], &iwork[iindw], &work[indgrs],
		     &z__[z_offset], ldz, &isuppz[1], &work[indwrk], &iwork[
		    iindwk], &iinfo);
#line 677 "cstemr.f"
	    if (iinfo != 0) {
#line 678 "cstemr.f"
		*info = abs(iinfo) + 20;
#line 679 "cstemr.f"
		return 0;
#line 680 "cstemr.f"
	    }
#line 681 "cstemr.f"
	} else {
/*           SLARRE computes eigenvalues of the (shifted) root representation */
/*           CLARRV returns the eigenvalues of the unshifted matrix. */
/*           However, if the eigenvectors are not desired by the user, we need */
/*           to apply the corresponding shifts from SLARRE to obtain the */
/*           eigenvalues of the original matrix. */
#line 687 "cstemr.f"
	    i__1 = *m;
#line 687 "cstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 688 "cstemr.f"
		itmp = iwork[iindbl + j - 1];
#line 689 "cstemr.f"
		w[j] += e[iwork[iinspl + itmp - 1]];
#line 690 "cstemr.f"
/* L20: */
#line 690 "cstemr.f"
	    }
#line 691 "cstemr.f"
	}

#line 694 "cstemr.f"
	if (*tryrac) {
/*           Refine computed eigenvalues so that they are relatively accurate */
/*           with respect to the original matrix T. */
#line 697 "cstemr.f"
	    ibegin = 1;
#line 698 "cstemr.f"
	    wbegin = 1;
#line 699 "cstemr.f"
	    i__1 = iwork[iindbl + *m - 1];
#line 699 "cstemr.f"
	    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 700 "cstemr.f"
		iend = iwork[iinspl + jblk - 1];
#line 701 "cstemr.f"
		in = iend - ibegin + 1;
#line 702 "cstemr.f"
		wend = wbegin - 1;
/*              check if any eigenvalues have to be refined in this block */
#line 704 "cstemr.f"
L36:
#line 705 "cstemr.f"
		if (wend < *m) {
#line 706 "cstemr.f"
		    if (iwork[iindbl + wend] == jblk) {
#line 707 "cstemr.f"
			++wend;
#line 708 "cstemr.f"
			goto L36;
#line 709 "cstemr.f"
		    }
#line 710 "cstemr.f"
		}
#line 711 "cstemr.f"
		if (wend < wbegin) {
#line 712 "cstemr.f"
		    ibegin = iend + 1;
#line 713 "cstemr.f"
		    goto L39;
#line 714 "cstemr.f"
		}
#line 716 "cstemr.f"
		offset = iwork[iindw + wbegin - 1] - 1;
#line 717 "cstemr.f"
		ifirst = iwork[iindw + wbegin - 1];
#line 718 "cstemr.f"
		ilast = iwork[iindw + wend - 1];
#line 719 "cstemr.f"
		rtol2 = eps * 4.;
#line 720 "cstemr.f"
		slarrj_(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 
			1], &ifirst, &ilast, &rtol2, &offset, &w[wbegin], &
			work[inderr + wbegin - 1], &work[indwrk], &iwork[
			iindwk], &pivmin, &tnrm, &iinfo);
#line 726 "cstemr.f"
		ibegin = iend + 1;
#line 727 "cstemr.f"
		wbegin = wend + 1;
#line 728 "cstemr.f"
L39:
#line 728 "cstemr.f"
		;
#line 728 "cstemr.f"
	    }
#line 729 "cstemr.f"
	}

/*        If matrix was scaled, then rescale eigenvalues appropriately. */

#line 733 "cstemr.f"
	if (scale != 1.) {
#line 734 "cstemr.f"
	    d__1 = 1. / scale;
#line 734 "cstemr.f"
	    sscal_(m, &d__1, &w[1], &c__1);
#line 735 "cstemr.f"
	}
#line 736 "cstemr.f"
    }

/*     If eigenvalues are not in increasing order, then sort them, */
/*     possibly along with eigenvectors. */

#line 741 "cstemr.f"
    if (nsplit > 1 || *n == 2) {
#line 742 "cstemr.f"
	if (! wantz) {
#line 743 "cstemr.f"
	    slasrt_("I", m, &w[1], &iinfo, (ftnlen)1);
#line 744 "cstemr.f"
	    if (iinfo != 0) {
#line 745 "cstemr.f"
		*info = 3;
#line 746 "cstemr.f"
		return 0;
#line 747 "cstemr.f"
	    }
#line 748 "cstemr.f"
	} else {
#line 749 "cstemr.f"
	    i__1 = *m - 1;
#line 749 "cstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 750 "cstemr.f"
		i__ = 0;
#line 751 "cstemr.f"
		tmp = w[j];
#line 752 "cstemr.f"
		i__2 = *m;
#line 752 "cstemr.f"
		for (jj = j + 1; jj <= i__2; ++jj) {
#line 753 "cstemr.f"
		    if (w[jj] < tmp) {
#line 754 "cstemr.f"
			i__ = jj;
#line 755 "cstemr.f"
			tmp = w[jj];
#line 756 "cstemr.f"
		    }
#line 757 "cstemr.f"
/* L50: */
#line 757 "cstemr.f"
		}
#line 758 "cstemr.f"
		if (i__ != 0) {
#line 759 "cstemr.f"
		    w[i__] = w[j];
#line 760 "cstemr.f"
		    w[j] = tmp;
#line 761 "cstemr.f"
		    if (wantz) {
#line 762 "cstemr.f"
			cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * 
				z_dim1 + 1], &c__1);
#line 763 "cstemr.f"
			itmp = isuppz[(i__ << 1) - 1];
#line 764 "cstemr.f"
			isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
#line 765 "cstemr.f"
			isuppz[(j << 1) - 1] = itmp;
#line 766 "cstemr.f"
			itmp = isuppz[i__ * 2];
#line 767 "cstemr.f"
			isuppz[i__ * 2] = isuppz[j * 2];
#line 768 "cstemr.f"
			isuppz[j * 2] = itmp;
#line 769 "cstemr.f"
		    }
#line 770 "cstemr.f"
		}
#line 771 "cstemr.f"
/* L60: */
#line 771 "cstemr.f"
	    }
#line 772 "cstemr.f"
	}
#line 773 "cstemr.f"
    }


#line 776 "cstemr.f"
    work[1] = (doublereal) lwmin;
#line 777 "cstemr.f"
    iwork[1] = liwmin;
#line 778 "cstemr.f"
    return 0;

/*     End of CSTEMR */

} /* cstemr_ */

