#line 1 "zstemr.f"
/* zstemr.f -- translated by f2c (version 20100827).
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

#line 1 "zstemr.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = .001;

/* > \brief \b ZSTEMR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSTEMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstemr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstemr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstemr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, */
/*                          M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, */
/*                          IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE */
/*       LOGICAL            TRYRAC */
/*       INTEGER            IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N */
/*       DOUBLE PRECISION VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ) */
/*       COMPLEX*16         Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSTEMR computes selected eigenvalues and, optionally, eigenvectors */
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
/* > 1.ZSTEMR works only on machines which follow IEEE-754 */
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
/* > ZSTEMR accepts complex workspace to facilitate interoperability */
/* > with ZUNMTR or ZUPMTR. */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the N diagonal elements of the tridiagonal matrix */
/* >          T. On exit, D is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the (N-1) subdiagonal elements of the tridiagonal */
/* >          matrix T in elements 1 to N-1 of E. E(N) need not be set on */
/* >          input, but is used internally as workspace. */
/* >          On exit, E is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION */
/* > */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* > */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* > */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* > */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, max(1,M) ) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK) */
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
/* >          > 0:  if INFO = 1X, internal error in DLARRE, */
/* >                if INFO = 2X, internal error in ZLARRV. */
/* >                Here, the digit X = ABS( IINFO ) < 10, where IINFO is */
/* >                the nonzero error code returned by DLARRE or */
/* >                ZLARRV, respectively. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA \n */

/*  ===================================================================== */
/* Subroutine */ int zstemr_(char *jobz, char *range, integer *n, doublereal *
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
    extern /* Subroutine */ int dlae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static integer inde2, itmp2;
    static doublereal rtol1, rtol2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    static integer indgp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo, iindw, ilast;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lwmin;
    static logical wantz;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlaev2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig;
    static integer ibegin;
    static logical indeig;
    static integer iindbl;
    static logical valeig;
    extern /* Subroutine */ int dlarrc_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, ftnlen), dlarre_(char *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static integer wbegin;
    static doublereal safmin;
    extern /* Subroutine */ int dlarrj_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *), xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer inderr, iindwk, indgrs, offset;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dlarrr_(integer *, doublereal *, doublereal *,
	     integer *), dlasrt_(char *, integer *, doublereal *, integer *, 
	    ftnlen);
    static doublereal thresh;
    static integer iinspl, indwrk, ifirst, liwmin, nzcmin;
    static doublereal pivmin;
    static integer nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int zlarrv_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublecomplex *, integer *, integer *, doublereal *,
	     integer *, integer *);
    static logical lquery, zquery;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 396 "zstemr.f"
    /* Parameter adjustments */
#line 396 "zstemr.f"
    --d__;
#line 396 "zstemr.f"
    --e;
#line 396 "zstemr.f"
    --w;
#line 396 "zstemr.f"
    z_dim1 = *ldz;
#line 396 "zstemr.f"
    z_offset = 1 + z_dim1;
#line 396 "zstemr.f"
    z__ -= z_offset;
#line 396 "zstemr.f"
    --isuppz;
#line 396 "zstemr.f"
    --work;
#line 396 "zstemr.f"
    --iwork;
#line 396 "zstemr.f"

#line 396 "zstemr.f"
    /* Function Body */
#line 396 "zstemr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 397 "zstemr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 398 "zstemr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 399 "zstemr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 401 "zstemr.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 402 "zstemr.f"
    zquery = *nzc == -1;
/*     DSTEMR needs WORK of size 6*N, IWORK of size 3*N. */
/*     In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N. */
/*     Furthermore, ZLARRV needs WORK of size 12*N, IWORK of size 7*N. */
#line 407 "zstemr.f"
    if (wantz) {
#line 408 "zstemr.f"
	lwmin = *n * 18;
#line 409 "zstemr.f"
	liwmin = *n * 10;
#line 410 "zstemr.f"
    } else {
/*        need less workspace if only the eigenvalues are wanted */
#line 412 "zstemr.f"
	lwmin = *n * 12;
#line 413 "zstemr.f"
	liwmin = *n << 3;
#line 414 "zstemr.f"
    }
#line 416 "zstemr.f"
    wl = 0.;
#line 417 "zstemr.f"
    wu = 0.;
#line 418 "zstemr.f"
    iil = 0;
#line 419 "zstemr.f"
    iiu = 0;
#line 420 "zstemr.f"
    nsplit = 0;
#line 422 "zstemr.f"
    if (valeig) {
/*        We do not reference VL, VU in the cases RANGE = 'I','A' */
/*        The interval (WL, WU] contains all the wanted eigenvalues. */
/*        It is either given by the user or computed in DLARRE. */
#line 426 "zstemr.f"
	wl = *vl;
#line 427 "zstemr.f"
	wu = *vu;
#line 428 "zstemr.f"
    } else if (indeig) {
/*        We do not reference IL, IU in the cases RANGE = 'V','A' */
#line 430 "zstemr.f"
	iil = *il;
#line 431 "zstemr.f"
	iiu = *iu;
#line 432 "zstemr.f"
    }

#line 434 "zstemr.f"
    *info = 0;
#line 435 "zstemr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 436 "zstemr.f"
	*info = -1;
#line 437 "zstemr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 438 "zstemr.f"
	*info = -2;
#line 439 "zstemr.f"
    } else if (*n < 0) {
#line 440 "zstemr.f"
	*info = -3;
#line 441 "zstemr.f"
    } else if (valeig && *n > 0 && wu <= wl) {
#line 442 "zstemr.f"
	*info = -7;
#line 443 "zstemr.f"
    } else if (indeig && (iil < 1 || iil > *n)) {
#line 444 "zstemr.f"
	*info = -8;
#line 445 "zstemr.f"
    } else if (indeig && (iiu < iil || iiu > *n)) {
#line 446 "zstemr.f"
	*info = -9;
#line 447 "zstemr.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 448 "zstemr.f"
	*info = -13;
#line 449 "zstemr.f"
    } else if (*lwork < lwmin && ! lquery) {
#line 450 "zstemr.f"
	*info = -17;
#line 451 "zstemr.f"
    } else if (*liwork < liwmin && ! lquery) {
#line 452 "zstemr.f"
	*info = -19;
#line 453 "zstemr.f"
    }

/*     Get machine constants. */

#line 457 "zstemr.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 458 "zstemr.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 459 "zstemr.f"
    smlnum = safmin / eps;
#line 460 "zstemr.f"
    bignum = 1. / smlnum;
#line 461 "zstemr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 462 "zstemr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 462 "zstemr.f"
    rmax = min(d__1,d__2);

#line 464 "zstemr.f"
    if (*info == 0) {
#line 465 "zstemr.f"
	work[1] = (doublereal) lwmin;
#line 466 "zstemr.f"
	iwork[1] = liwmin;

#line 468 "zstemr.f"
	if (wantz && alleig) {
#line 469 "zstemr.f"
	    nzcmin = *n;
#line 470 "zstemr.f"
	} else if (wantz && valeig) {
#line 471 "zstemr.f"
	    dlarrc_("T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp, &
		    itmp2, info, (ftnlen)1);
#line 473 "zstemr.f"
	} else if (wantz && indeig) {
#line 474 "zstemr.f"
	    nzcmin = iiu - iil + 1;
#line 475 "zstemr.f"
	} else {
/*           WANTZ .EQ. FALSE. */
#line 477 "zstemr.f"
	    nzcmin = 0;
#line 478 "zstemr.f"
	}
#line 479 "zstemr.f"
	if (zquery && *info == 0) {
#line 480 "zstemr.f"
	    i__1 = z_dim1 + 1;
#line 480 "zstemr.f"
	    z__[i__1].r = (doublereal) nzcmin, z__[i__1].i = 0.;
#line 481 "zstemr.f"
	} else if (*nzc < nzcmin && ! zquery) {
#line 482 "zstemr.f"
	    *info = -14;
#line 483 "zstemr.f"
	}
#line 484 "zstemr.f"
    }
#line 486 "zstemr.f"
    if (*info != 0) {

#line 488 "zstemr.f"
	i__1 = -(*info);
#line 488 "zstemr.f"
	xerbla_("ZSTEMR", &i__1, (ftnlen)6);

#line 490 "zstemr.f"
	return 0;
#line 491 "zstemr.f"
    } else if (lquery || zquery) {
#line 492 "zstemr.f"
	return 0;
#line 493 "zstemr.f"
    }

/*     Handle N = 0, 1, and 2 cases immediately */

#line 497 "zstemr.f"
    *m = 0;
#line 498 "zstemr.f"
    if (*n == 0) {
#line 498 "zstemr.f"
	return 0;
#line 498 "zstemr.f"
    }

#line 501 "zstemr.f"
    if (*n == 1) {
#line 502 "zstemr.f"
	if (alleig || indeig) {
#line 503 "zstemr.f"
	    *m = 1;
#line 504 "zstemr.f"
	    w[1] = d__[1];
#line 505 "zstemr.f"
	} else {
#line 506 "zstemr.f"
	    if (wl < d__[1] && wu >= d__[1]) {
#line 507 "zstemr.f"
		*m = 1;
#line 508 "zstemr.f"
		w[1] = d__[1];
#line 509 "zstemr.f"
	    }
#line 510 "zstemr.f"
	}
#line 511 "zstemr.f"
	if (wantz && ! zquery) {
#line 512 "zstemr.f"
	    i__1 = z_dim1 + 1;
#line 512 "zstemr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 513 "zstemr.f"
	    isuppz[1] = 1;
#line 514 "zstemr.f"
	    isuppz[2] = 1;
#line 515 "zstemr.f"
	}
#line 516 "zstemr.f"
	return 0;
#line 517 "zstemr.f"
    }

#line 519 "zstemr.f"
    if (*n == 2) {
#line 520 "zstemr.f"
	if (! wantz) {
#line 521 "zstemr.f"
	    dlae2_(&d__[1], &e[1], &d__[2], &r1, &r2);
#line 522 "zstemr.f"
	} else if (wantz && ! zquery) {
#line 523 "zstemr.f"
	    dlaev2_(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
#line 524 "zstemr.f"
	}
#line 525 "zstemr.f"
	if (alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1) {
#line 529 "zstemr.f"
	    ++(*m);
#line 530 "zstemr.f"
	    w[*m] = r2;
#line 531 "zstemr.f"
	    if (wantz && ! zquery) {
#line 532 "zstemr.f"
		i__1 = *m * z_dim1 + 1;
#line 532 "zstemr.f"
		d__1 = -sn;
#line 532 "zstemr.f"
		z__[i__1].r = d__1, z__[i__1].i = 0.;
#line 533 "zstemr.f"
		i__1 = *m * z_dim1 + 2;
#line 533 "zstemr.f"
		z__[i__1].r = cs, z__[i__1].i = 0.;
/*              Note: At most one of SN and CS can be zero. */
#line 535 "zstemr.f"
		if (sn != 0.) {
#line 536 "zstemr.f"
		    if (cs != 0.) {
#line 537 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 538 "zstemr.f"
			isuppz[*m * 2] = 2;
#line 539 "zstemr.f"
		    } else {
#line 540 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 541 "zstemr.f"
			isuppz[*m * 2] = 1;
#line 542 "zstemr.f"
		    }
#line 543 "zstemr.f"
		} else {
#line 544 "zstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 545 "zstemr.f"
		    isuppz[*m * 2] = 2;
#line 546 "zstemr.f"
		}
#line 547 "zstemr.f"
	    }
#line 548 "zstemr.f"
	}
#line 549 "zstemr.f"
	if (alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2) {
#line 553 "zstemr.f"
	    ++(*m);
#line 554 "zstemr.f"
	    w[*m] = r1;
#line 555 "zstemr.f"
	    if (wantz && ! zquery) {
#line 556 "zstemr.f"
		i__1 = *m * z_dim1 + 1;
#line 556 "zstemr.f"
		z__[i__1].r = cs, z__[i__1].i = 0.;
#line 557 "zstemr.f"
		i__1 = *m * z_dim1 + 2;
#line 557 "zstemr.f"
		z__[i__1].r = sn, z__[i__1].i = 0.;
/*              Note: At most one of SN and CS can be zero. */
#line 559 "zstemr.f"
		if (sn != 0.) {
#line 560 "zstemr.f"
		    if (cs != 0.) {
#line 561 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 562 "zstemr.f"
			isuppz[*m * 2] = 2;
#line 563 "zstemr.f"
		    } else {
#line 564 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 565 "zstemr.f"
			isuppz[*m * 2] = 1;
#line 566 "zstemr.f"
		    }
#line 567 "zstemr.f"
		} else {
#line 568 "zstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 569 "zstemr.f"
		    isuppz[*m * 2] = 2;
#line 570 "zstemr.f"
		}
#line 571 "zstemr.f"
	    }
#line 572 "zstemr.f"
	}
#line 573 "zstemr.f"
    } else {
/*        Continue with general N */
#line 577 "zstemr.f"
	indgrs = 1;
#line 578 "zstemr.f"
	inderr = (*n << 1) + 1;
#line 579 "zstemr.f"
	indgp = *n * 3 + 1;
#line 580 "zstemr.f"
	indd = (*n << 2) + 1;
#line 581 "zstemr.f"
	inde2 = *n * 5 + 1;
#line 582 "zstemr.f"
	indwrk = *n * 6 + 1;

#line 584 "zstemr.f"
	iinspl = 1;
#line 585 "zstemr.f"
	iindbl = *n + 1;
#line 586 "zstemr.f"
	iindw = (*n << 1) + 1;
#line 587 "zstemr.f"
	iindwk = *n * 3 + 1;

/*        Scale matrix to allowable range, if necessary. */
/*        The allowable range is related to the PIVMIN parameter; see the */
/*        comments in DLARRD.  The preference for scaling small values */
/*        up is heuristic; we expect users' matrices not to be close to the */
/*        RMAX threshold. */

#line 595 "zstemr.f"
	scale = 1.;
#line 596 "zstemr.f"
	tnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 597 "zstemr.f"
	if (tnrm > 0. && tnrm < rmin) {
#line 598 "zstemr.f"
	    scale = rmin / tnrm;
#line 599 "zstemr.f"
	} else if (tnrm > rmax) {
#line 600 "zstemr.f"
	    scale = rmax / tnrm;
#line 601 "zstemr.f"
	}
#line 602 "zstemr.f"
	if (scale != 1.) {
#line 603 "zstemr.f"
	    dscal_(n, &scale, &d__[1], &c__1);
#line 604 "zstemr.f"
	    i__1 = *n - 1;
#line 604 "zstemr.f"
	    dscal_(&i__1, &scale, &e[1], &c__1);
#line 605 "zstemr.f"
	    tnrm *= scale;
#line 606 "zstemr.f"
	    if (valeig) {
/*              If eigenvalues in interval have to be found, */
/*              scale (WL, WU] accordingly */
#line 609 "zstemr.f"
		wl *= scale;
#line 610 "zstemr.f"
		wu *= scale;
#line 611 "zstemr.f"
	    }
#line 612 "zstemr.f"
	}

/*        Compute the desired eigenvalues of the tridiagonal after splitting */
/*        into smaller subblocks if the corresponding off-diagonal elements */
/*        are small */
/*        THRESH is the splitting parameter for DLARRE */
/*        A negative THRESH forces the old splitting criterion based on the */
/*        size of the off-diagonal. A positive THRESH switches to splitting */
/*        which preserves relative accuracy. */

#line 622 "zstemr.f"
	if (*tryrac) {
/*           Test whether the matrix warrants the more expensive relative approach. */
#line 624 "zstemr.f"
	    dlarrr_(n, &d__[1], &e[1], &iinfo);
#line 625 "zstemr.f"
	} else {
/*           The user does not care about relative accurately eigenvalues */
#line 627 "zstemr.f"
	    iinfo = -1;
#line 628 "zstemr.f"
	}
/*        Set the splitting criterion */
#line 630 "zstemr.f"
	if (iinfo == 0) {
#line 631 "zstemr.f"
	    thresh = eps;
#line 632 "zstemr.f"
	} else {
#line 633 "zstemr.f"
	    thresh = -eps;
/*           relative accuracy is desired but T does not guarantee it */
#line 635 "zstemr.f"
	    *tryrac = FALSE_;
#line 636 "zstemr.f"
	}

#line 638 "zstemr.f"
	if (*tryrac) {
/*           Copy original diagonal, needed to guarantee relative accuracy */
#line 640 "zstemr.f"
	    dcopy_(n, &d__[1], &c__1, &work[indd], &c__1);
#line 641 "zstemr.f"
	}
/*        Store the squares of the offdiagonal values of T */
#line 643 "zstemr.f"
	i__1 = *n - 1;
#line 643 "zstemr.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
#line 644 "zstemr.f"
	    d__1 = e[j];
#line 644 "zstemr.f"
	    work[inde2 + j - 1] = d__1 * d__1;
#line 645 "zstemr.f"
/* L5: */
#line 645 "zstemr.f"
	}
/*        Set the tolerance parameters for bisection */
#line 648 "zstemr.f"
	if (! wantz) {
/*           DLARRE computes the eigenvalues to full precision. */
#line 650 "zstemr.f"
	    rtol1 = eps * 4.;
#line 651 "zstemr.f"
	    rtol2 = eps * 4.;
#line 652 "zstemr.f"
	} else {
/*           DLARRE computes the eigenvalues to less than full precision. */
/*           ZLARRV will refine the eigenvalue approximations, and we only */
/*           need less accurate initial bisection in DLARRE. */
/*           Note: these settings do only affect the subset case and DLARRE */
#line 657 "zstemr.f"
	    rtol1 = sqrt(eps);
/* Computing MAX */
#line 658 "zstemr.f"
	    d__1 = sqrt(eps) * .005, d__2 = eps * 4.;
#line 658 "zstemr.f"
	    rtol2 = max(d__1,d__2);
#line 659 "zstemr.f"
	}
#line 660 "zstemr.f"
	dlarre_(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2], 
		&rtol1, &rtol2, &thresh, &nsplit, &iwork[iinspl], m, &w[1], &
		work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw], &
		work[indgrs], &pivmin, &work[indwrk], &iwork[iindwk], &iinfo, 
		(ftnlen)1);
#line 666 "zstemr.f"
	if (iinfo != 0) {
#line 667 "zstemr.f"
	    *info = abs(iinfo) + 10;
#line 668 "zstemr.f"
	    return 0;
#line 669 "zstemr.f"
	}
/*        Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired */
/*        part of the spectrum. All desired eigenvalues are contained in */
/*        (WL,WU] */
#line 675 "zstemr.f"
	if (wantz) {

/*           Compute the desired eigenvectors corresponding to the computed */
/*           eigenvalues */

#line 680 "zstemr.f"
	    zlarrv_(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &
		    c__1, m, &c_b18, &rtol1, &rtol2, &w[1], &work[inderr], &
		    work[indgp], &iwork[iindbl], &iwork[iindw], &work[indgrs],
		     &z__[z_offset], ldz, &isuppz[1], &work[indwrk], &iwork[
		    iindwk], &iinfo);
#line 686 "zstemr.f"
	    if (iinfo != 0) {
#line 687 "zstemr.f"
		*info = abs(iinfo) + 20;
#line 688 "zstemr.f"
		return 0;
#line 689 "zstemr.f"
	    }
#line 690 "zstemr.f"
	} else {
/*           DLARRE computes eigenvalues of the (shifted) root representation */
/*           ZLARRV returns the eigenvalues of the unshifted matrix. */
/*           However, if the eigenvectors are not desired by the user, we need */
/*           to apply the corresponding shifts from DLARRE to obtain the */
/*           eigenvalues of the original matrix. */
#line 696 "zstemr.f"
	    i__1 = *m;
#line 696 "zstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 697 "zstemr.f"
		itmp = iwork[iindbl + j - 1];
#line 698 "zstemr.f"
		w[j] += e[iwork[iinspl + itmp - 1]];
#line 699 "zstemr.f"
/* L20: */
#line 699 "zstemr.f"
	    }
#line 700 "zstemr.f"
	}

#line 703 "zstemr.f"
	if (*tryrac) {
/*           Refine computed eigenvalues so that they are relatively accurate */
/*           with respect to the original matrix T. */
#line 706 "zstemr.f"
	    ibegin = 1;
#line 707 "zstemr.f"
	    wbegin = 1;
#line 708 "zstemr.f"
	    i__1 = iwork[iindbl + *m - 1];
#line 708 "zstemr.f"
	    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 709 "zstemr.f"
		iend = iwork[iinspl + jblk - 1];
#line 710 "zstemr.f"
		in = iend - ibegin + 1;
#line 711 "zstemr.f"
		wend = wbegin - 1;
/*              check if any eigenvalues have to be refined in this block */
#line 713 "zstemr.f"
L36:
#line 714 "zstemr.f"
		if (wend < *m) {
#line 715 "zstemr.f"
		    if (iwork[iindbl + wend] == jblk) {
#line 716 "zstemr.f"
			++wend;
#line 717 "zstemr.f"
			goto L36;
#line 718 "zstemr.f"
		    }
#line 719 "zstemr.f"
		}
#line 720 "zstemr.f"
		if (wend < wbegin) {
#line 721 "zstemr.f"
		    ibegin = iend + 1;
#line 722 "zstemr.f"
		    goto L39;
#line 723 "zstemr.f"
		}
#line 725 "zstemr.f"
		offset = iwork[iindw + wbegin - 1] - 1;
#line 726 "zstemr.f"
		ifirst = iwork[iindw + wbegin - 1];
#line 727 "zstemr.f"
		ilast = iwork[iindw + wend - 1];
#line 728 "zstemr.f"
		rtol2 = eps * 4.;
#line 729 "zstemr.f"
		dlarrj_(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 
			1], &ifirst, &ilast, &rtol2, &offset, &w[wbegin], &
			work[inderr + wbegin - 1], &work[indwrk], &iwork[
			iindwk], &pivmin, &tnrm, &iinfo);
#line 735 "zstemr.f"
		ibegin = iend + 1;
#line 736 "zstemr.f"
		wbegin = wend + 1;
#line 737 "zstemr.f"
L39:
#line 737 "zstemr.f"
		;
#line 737 "zstemr.f"
	    }
#line 738 "zstemr.f"
	}

/*        If matrix was scaled, then rescale eigenvalues appropriately. */

#line 742 "zstemr.f"
	if (scale != 1.) {
#line 743 "zstemr.f"
	    d__1 = 1. / scale;
#line 743 "zstemr.f"
	    dscal_(m, &d__1, &w[1], &c__1);
#line 744 "zstemr.f"
	}
#line 745 "zstemr.f"
    }

/*     If eigenvalues are not in increasing order, then sort them, */
/*     possibly along with eigenvectors. */

#line 750 "zstemr.f"
    if (nsplit > 1 || *n == 2) {
#line 751 "zstemr.f"
	if (! wantz) {
#line 752 "zstemr.f"
	    dlasrt_("I", m, &w[1], &iinfo, (ftnlen)1);
#line 753 "zstemr.f"
	    if (iinfo != 0) {
#line 754 "zstemr.f"
		*info = 3;
#line 755 "zstemr.f"
		return 0;
#line 756 "zstemr.f"
	    }
#line 757 "zstemr.f"
	} else {
#line 758 "zstemr.f"
	    i__1 = *m - 1;
#line 758 "zstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 759 "zstemr.f"
		i__ = 0;
#line 760 "zstemr.f"
		tmp = w[j];
#line 761 "zstemr.f"
		i__2 = *m;
#line 761 "zstemr.f"
		for (jj = j + 1; jj <= i__2; ++jj) {
#line 762 "zstemr.f"
		    if (w[jj] < tmp) {
#line 763 "zstemr.f"
			i__ = jj;
#line 764 "zstemr.f"
			tmp = w[jj];
#line 765 "zstemr.f"
		    }
#line 766 "zstemr.f"
/* L50: */
#line 766 "zstemr.f"
		}
#line 767 "zstemr.f"
		if (i__ != 0) {
#line 768 "zstemr.f"
		    w[i__] = w[j];
#line 769 "zstemr.f"
		    w[j] = tmp;
#line 770 "zstemr.f"
		    if (wantz) {
#line 771 "zstemr.f"
			zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * 
				z_dim1 + 1], &c__1);
#line 772 "zstemr.f"
			itmp = isuppz[(i__ << 1) - 1];
#line 773 "zstemr.f"
			isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
#line 774 "zstemr.f"
			isuppz[(j << 1) - 1] = itmp;
#line 775 "zstemr.f"
			itmp = isuppz[i__ * 2];
#line 776 "zstemr.f"
			isuppz[i__ * 2] = isuppz[j * 2];
#line 777 "zstemr.f"
			isuppz[j * 2] = itmp;
#line 778 "zstemr.f"
		    }
#line 779 "zstemr.f"
		}
#line 780 "zstemr.f"
/* L60: */
#line 780 "zstemr.f"
	    }
#line 781 "zstemr.f"
	}
#line 782 "zstemr.f"
    }


#line 785 "zstemr.f"
    work[1] = (doublereal) lwmin;
#line 786 "zstemr.f"
    iwork[1] = liwmin;
#line 787 "zstemr.f"
    return 0;

/*     End of ZSTEMR */

} /* zstemr_ */

