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
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
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

/* > \date November 2013 */

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


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

#line 387 "zstemr.f"
    /* Parameter adjustments */
#line 387 "zstemr.f"
    --d__;
#line 387 "zstemr.f"
    --e;
#line 387 "zstemr.f"
    --w;
#line 387 "zstemr.f"
    z_dim1 = *ldz;
#line 387 "zstemr.f"
    z_offset = 1 + z_dim1;
#line 387 "zstemr.f"
    z__ -= z_offset;
#line 387 "zstemr.f"
    --isuppz;
#line 387 "zstemr.f"
    --work;
#line 387 "zstemr.f"
    --iwork;
#line 387 "zstemr.f"

#line 387 "zstemr.f"
    /* Function Body */
#line 387 "zstemr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 388 "zstemr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 389 "zstemr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 390 "zstemr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 392 "zstemr.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 393 "zstemr.f"
    zquery = *nzc == -1;
/*     DSTEMR needs WORK of size 6*N, IWORK of size 3*N. */
/*     In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N. */
/*     Furthermore, ZLARRV needs WORK of size 12*N, IWORK of size 7*N. */
#line 398 "zstemr.f"
    if (wantz) {
#line 399 "zstemr.f"
	lwmin = *n * 18;
#line 400 "zstemr.f"
	liwmin = *n * 10;
#line 401 "zstemr.f"
    } else {
/*        need less workspace if only the eigenvalues are wanted */
#line 403 "zstemr.f"
	lwmin = *n * 12;
#line 404 "zstemr.f"
	liwmin = *n << 3;
#line 405 "zstemr.f"
    }
#line 407 "zstemr.f"
    wl = 0.;
#line 408 "zstemr.f"
    wu = 0.;
#line 409 "zstemr.f"
    iil = 0;
#line 410 "zstemr.f"
    iiu = 0;
#line 411 "zstemr.f"
    nsplit = 0;
#line 413 "zstemr.f"
    if (valeig) {
/*        We do not reference VL, VU in the cases RANGE = 'I','A' */
/*        The interval (WL, WU] contains all the wanted eigenvalues. */
/*        It is either given by the user or computed in DLARRE. */
#line 417 "zstemr.f"
	wl = *vl;
#line 418 "zstemr.f"
	wu = *vu;
#line 419 "zstemr.f"
    } else if (indeig) {
/*        We do not reference IL, IU in the cases RANGE = 'V','A' */
#line 421 "zstemr.f"
	iil = *il;
#line 422 "zstemr.f"
	iiu = *iu;
#line 423 "zstemr.f"
    }

#line 425 "zstemr.f"
    *info = 0;
#line 426 "zstemr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 427 "zstemr.f"
	*info = -1;
#line 428 "zstemr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 429 "zstemr.f"
	*info = -2;
#line 430 "zstemr.f"
    } else if (*n < 0) {
#line 431 "zstemr.f"
	*info = -3;
#line 432 "zstemr.f"
    } else if (valeig && *n > 0 && wu <= wl) {
#line 433 "zstemr.f"
	*info = -7;
#line 434 "zstemr.f"
    } else if (indeig && (iil < 1 || iil > *n)) {
#line 435 "zstemr.f"
	*info = -8;
#line 436 "zstemr.f"
    } else if (indeig && (iiu < iil || iiu > *n)) {
#line 437 "zstemr.f"
	*info = -9;
#line 438 "zstemr.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 439 "zstemr.f"
	*info = -13;
#line 440 "zstemr.f"
    } else if (*lwork < lwmin && ! lquery) {
#line 441 "zstemr.f"
	*info = -17;
#line 442 "zstemr.f"
    } else if (*liwork < liwmin && ! lquery) {
#line 443 "zstemr.f"
	*info = -19;
#line 444 "zstemr.f"
    }

/*     Get machine constants. */

#line 448 "zstemr.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 449 "zstemr.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 450 "zstemr.f"
    smlnum = safmin / eps;
#line 451 "zstemr.f"
    bignum = 1. / smlnum;
#line 452 "zstemr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 453 "zstemr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 453 "zstemr.f"
    rmax = min(d__1,d__2);

#line 455 "zstemr.f"
    if (*info == 0) {
#line 456 "zstemr.f"
	work[1] = (doublereal) lwmin;
#line 457 "zstemr.f"
	iwork[1] = liwmin;

#line 459 "zstemr.f"
	if (wantz && alleig) {
#line 460 "zstemr.f"
	    nzcmin = *n;
#line 461 "zstemr.f"
	} else if (wantz && valeig) {
#line 462 "zstemr.f"
	    dlarrc_("T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp, &
		    itmp2, info, (ftnlen)1);
#line 464 "zstemr.f"
	} else if (wantz && indeig) {
#line 465 "zstemr.f"
	    nzcmin = iiu - iil + 1;
#line 466 "zstemr.f"
	} else {
/*           WANTZ .EQ. FALSE. */
#line 468 "zstemr.f"
	    nzcmin = 0;
#line 469 "zstemr.f"
	}
#line 470 "zstemr.f"
	if (zquery && *info == 0) {
#line 471 "zstemr.f"
	    i__1 = z_dim1 + 1;
#line 471 "zstemr.f"
	    z__[i__1].r = (doublereal) nzcmin, z__[i__1].i = 0.;
#line 472 "zstemr.f"
	} else if (*nzc < nzcmin && ! zquery) {
#line 473 "zstemr.f"
	    *info = -14;
#line 474 "zstemr.f"
	}
#line 475 "zstemr.f"
    }
#line 477 "zstemr.f"
    if (*info != 0) {

#line 479 "zstemr.f"
	i__1 = -(*info);
#line 479 "zstemr.f"
	xerbla_("ZSTEMR", &i__1, (ftnlen)6);

#line 481 "zstemr.f"
	return 0;
#line 482 "zstemr.f"
    } else if (lquery || zquery) {
#line 483 "zstemr.f"
	return 0;
#line 484 "zstemr.f"
    }

/*     Handle N = 0, 1, and 2 cases immediately */

#line 488 "zstemr.f"
    *m = 0;
#line 489 "zstemr.f"
    if (*n == 0) {
#line 489 "zstemr.f"
	return 0;
#line 489 "zstemr.f"
    }

#line 492 "zstemr.f"
    if (*n == 1) {
#line 493 "zstemr.f"
	if (alleig || indeig) {
#line 494 "zstemr.f"
	    *m = 1;
#line 495 "zstemr.f"
	    w[1] = d__[1];
#line 496 "zstemr.f"
	} else {
#line 497 "zstemr.f"
	    if (wl < d__[1] && wu >= d__[1]) {
#line 498 "zstemr.f"
		*m = 1;
#line 499 "zstemr.f"
		w[1] = d__[1];
#line 500 "zstemr.f"
	    }
#line 501 "zstemr.f"
	}
#line 502 "zstemr.f"
	if (wantz && ! zquery) {
#line 503 "zstemr.f"
	    i__1 = z_dim1 + 1;
#line 503 "zstemr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 504 "zstemr.f"
	    isuppz[1] = 1;
#line 505 "zstemr.f"
	    isuppz[2] = 1;
#line 506 "zstemr.f"
	}
#line 507 "zstemr.f"
	return 0;
#line 508 "zstemr.f"
    }

#line 510 "zstemr.f"
    if (*n == 2) {
#line 511 "zstemr.f"
	if (! wantz) {
#line 512 "zstemr.f"
	    dlae2_(&d__[1], &e[1], &d__[2], &r1, &r2);
#line 513 "zstemr.f"
	} else if (wantz && ! zquery) {
#line 514 "zstemr.f"
	    dlaev2_(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
#line 515 "zstemr.f"
	}
#line 516 "zstemr.f"
	if (alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1) {
#line 520 "zstemr.f"
	    ++(*m);
#line 521 "zstemr.f"
	    w[*m] = r2;
#line 522 "zstemr.f"
	    if (wantz && ! zquery) {
#line 523 "zstemr.f"
		i__1 = *m * z_dim1 + 1;
#line 523 "zstemr.f"
		d__1 = -sn;
#line 523 "zstemr.f"
		z__[i__1].r = d__1, z__[i__1].i = 0.;
#line 524 "zstemr.f"
		i__1 = *m * z_dim1 + 2;
#line 524 "zstemr.f"
		z__[i__1].r = cs, z__[i__1].i = 0.;
/*              Note: At most one of SN and CS can be zero. */
#line 526 "zstemr.f"
		if (sn != 0.) {
#line 527 "zstemr.f"
		    if (cs != 0.) {
#line 528 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 529 "zstemr.f"
			isuppz[(*m << 1) - 1] = 2;
#line 530 "zstemr.f"
		    } else {
#line 531 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 532 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 533 "zstemr.f"
		    }
#line 534 "zstemr.f"
		} else {
#line 535 "zstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 536 "zstemr.f"
		    isuppz[*m * 2] = 2;
#line 537 "zstemr.f"
		}
#line 538 "zstemr.f"
	    }
#line 539 "zstemr.f"
	}
#line 540 "zstemr.f"
	if (alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2) {
#line 544 "zstemr.f"
	    ++(*m);
#line 545 "zstemr.f"
	    w[*m] = r1;
#line 546 "zstemr.f"
	    if (wantz && ! zquery) {
#line 547 "zstemr.f"
		i__1 = *m * z_dim1 + 1;
#line 547 "zstemr.f"
		z__[i__1].r = cs, z__[i__1].i = 0.;
#line 548 "zstemr.f"
		i__1 = *m * z_dim1 + 2;
#line 548 "zstemr.f"
		z__[i__1].r = sn, z__[i__1].i = 0.;
/*              Note: At most one of SN and CS can be zero. */
#line 550 "zstemr.f"
		if (sn != 0.) {
#line 551 "zstemr.f"
		    if (cs != 0.) {
#line 552 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 553 "zstemr.f"
			isuppz[(*m << 1) - 1] = 2;
#line 554 "zstemr.f"
		    } else {
#line 555 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 556 "zstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 557 "zstemr.f"
		    }
#line 558 "zstemr.f"
		} else {
#line 559 "zstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 560 "zstemr.f"
		    isuppz[*m * 2] = 2;
#line 561 "zstemr.f"
		}
#line 562 "zstemr.f"
	    }
#line 563 "zstemr.f"
	}
#line 564 "zstemr.f"
    } else {
/*        Continue with general N */
#line 568 "zstemr.f"
	indgrs = 1;
#line 569 "zstemr.f"
	inderr = (*n << 1) + 1;
#line 570 "zstemr.f"
	indgp = *n * 3 + 1;
#line 571 "zstemr.f"
	indd = (*n << 2) + 1;
#line 572 "zstemr.f"
	inde2 = *n * 5 + 1;
#line 573 "zstemr.f"
	indwrk = *n * 6 + 1;

#line 575 "zstemr.f"
	iinspl = 1;
#line 576 "zstemr.f"
	iindbl = *n + 1;
#line 577 "zstemr.f"
	iindw = (*n << 1) + 1;
#line 578 "zstemr.f"
	iindwk = *n * 3 + 1;

/*        Scale matrix to allowable range, if necessary. */
/*        The allowable range is related to the PIVMIN parameter; see the */
/*        comments in DLARRD.  The preference for scaling small values */
/*        up is heuristic; we expect users' matrices not to be close to the */
/*        RMAX threshold. */

#line 586 "zstemr.f"
	scale = 1.;
#line 587 "zstemr.f"
	tnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 588 "zstemr.f"
	if (tnrm > 0. && tnrm < rmin) {
#line 589 "zstemr.f"
	    scale = rmin / tnrm;
#line 590 "zstemr.f"
	} else if (tnrm > rmax) {
#line 591 "zstemr.f"
	    scale = rmax / tnrm;
#line 592 "zstemr.f"
	}
#line 593 "zstemr.f"
	if (scale != 1.) {
#line 594 "zstemr.f"
	    dscal_(n, &scale, &d__[1], &c__1);
#line 595 "zstemr.f"
	    i__1 = *n - 1;
#line 595 "zstemr.f"
	    dscal_(&i__1, &scale, &e[1], &c__1);
#line 596 "zstemr.f"
	    tnrm *= scale;
#line 597 "zstemr.f"
	    if (valeig) {
/*              If eigenvalues in interval have to be found, */
/*              scale (WL, WU] accordingly */
#line 600 "zstemr.f"
		wl *= scale;
#line 601 "zstemr.f"
		wu *= scale;
#line 602 "zstemr.f"
	    }
#line 603 "zstemr.f"
	}

/*        Compute the desired eigenvalues of the tridiagonal after splitting */
/*        into smaller subblocks if the corresponding off-diagonal elements */
/*        are small */
/*        THRESH is the splitting parameter for DLARRE */
/*        A negative THRESH forces the old splitting criterion based on the */
/*        size of the off-diagonal. A positive THRESH switches to splitting */
/*        which preserves relative accuracy. */

#line 613 "zstemr.f"
	if (*tryrac) {
/*           Test whether the matrix warrants the more expensive relative approach. */
#line 615 "zstemr.f"
	    dlarrr_(n, &d__[1], &e[1], &iinfo);
#line 616 "zstemr.f"
	} else {
/*           The user does not care about relative accurately eigenvalues */
#line 618 "zstemr.f"
	    iinfo = -1;
#line 619 "zstemr.f"
	}
/*        Set the splitting criterion */
#line 621 "zstemr.f"
	if (iinfo == 0) {
#line 622 "zstemr.f"
	    thresh = eps;
#line 623 "zstemr.f"
	} else {
#line 624 "zstemr.f"
	    thresh = -eps;
/*           relative accuracy is desired but T does not guarantee it */
#line 626 "zstemr.f"
	    *tryrac = FALSE_;
#line 627 "zstemr.f"
	}

#line 629 "zstemr.f"
	if (*tryrac) {
/*           Copy original diagonal, needed to guarantee relative accuracy */
#line 631 "zstemr.f"
	    dcopy_(n, &d__[1], &c__1, &work[indd], &c__1);
#line 632 "zstemr.f"
	}
/*        Store the squares of the offdiagonal values of T */
#line 634 "zstemr.f"
	i__1 = *n - 1;
#line 634 "zstemr.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
#line 635 "zstemr.f"
	    d__1 = e[j];
#line 635 "zstemr.f"
	    work[inde2 + j - 1] = d__1 * d__1;
#line 636 "zstemr.f"
/* L5: */
#line 636 "zstemr.f"
	}
/*        Set the tolerance parameters for bisection */
#line 639 "zstemr.f"
	if (! wantz) {
/*           DLARRE computes the eigenvalues to full precision. */
#line 641 "zstemr.f"
	    rtol1 = eps * 4.;
#line 642 "zstemr.f"
	    rtol2 = eps * 4.;
#line 643 "zstemr.f"
	} else {
/*           DLARRE computes the eigenvalues to less than full precision. */
/*           ZLARRV will refine the eigenvalue approximations, and we only */
/*           need less accurate initial bisection in DLARRE. */
/*           Note: these settings do only affect the subset case and DLARRE */
#line 648 "zstemr.f"
	    rtol1 = sqrt(eps);
/* Computing MAX */
#line 649 "zstemr.f"
	    d__1 = sqrt(eps) * .005, d__2 = eps * 4.;
#line 649 "zstemr.f"
	    rtol2 = max(d__1,d__2);
#line 650 "zstemr.f"
	}
#line 651 "zstemr.f"
	dlarre_(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2], 
		&rtol1, &rtol2, &thresh, &nsplit, &iwork[iinspl], m, &w[1], &
		work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw], &
		work[indgrs], &pivmin, &work[indwrk], &iwork[iindwk], &iinfo, 
		(ftnlen)1);
#line 657 "zstemr.f"
	if (iinfo != 0) {
#line 658 "zstemr.f"
	    *info = abs(iinfo) + 10;
#line 659 "zstemr.f"
	    return 0;
#line 660 "zstemr.f"
	}
/*        Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired */
/*        part of the spectrum. All desired eigenvalues are contained in */
/*        (WL,WU] */
#line 666 "zstemr.f"
	if (wantz) {

/*           Compute the desired eigenvectors corresponding to the computed */
/*           eigenvalues */

#line 671 "zstemr.f"
	    zlarrv_(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &
		    c__1, m, &c_b18, &rtol1, &rtol2, &w[1], &work[inderr], &
		    work[indgp], &iwork[iindbl], &iwork[iindw], &work[indgrs],
		     &z__[z_offset], ldz, &isuppz[1], &work[indwrk], &iwork[
		    iindwk], &iinfo);
#line 677 "zstemr.f"
	    if (iinfo != 0) {
#line 678 "zstemr.f"
		*info = abs(iinfo) + 20;
#line 679 "zstemr.f"
		return 0;
#line 680 "zstemr.f"
	    }
#line 681 "zstemr.f"
	} else {
/*           DLARRE computes eigenvalues of the (shifted) root representation */
/*           ZLARRV returns the eigenvalues of the unshifted matrix. */
/*           However, if the eigenvectors are not desired by the user, we need */
/*           to apply the corresponding shifts from DLARRE to obtain the */
/*           eigenvalues of the original matrix. */
#line 687 "zstemr.f"
	    i__1 = *m;
#line 687 "zstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 688 "zstemr.f"
		itmp = iwork[iindbl + j - 1];
#line 689 "zstemr.f"
		w[j] += e[iwork[iinspl + itmp - 1]];
#line 690 "zstemr.f"
/* L20: */
#line 690 "zstemr.f"
	    }
#line 691 "zstemr.f"
	}

#line 694 "zstemr.f"
	if (*tryrac) {
/*           Refine computed eigenvalues so that they are relatively accurate */
/*           with respect to the original matrix T. */
#line 697 "zstemr.f"
	    ibegin = 1;
#line 698 "zstemr.f"
	    wbegin = 1;
#line 699 "zstemr.f"
	    i__1 = iwork[iindbl + *m - 1];
#line 699 "zstemr.f"
	    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 700 "zstemr.f"
		iend = iwork[iinspl + jblk - 1];
#line 701 "zstemr.f"
		in = iend - ibegin + 1;
#line 702 "zstemr.f"
		wend = wbegin - 1;
/*              check if any eigenvalues have to be refined in this block */
#line 704 "zstemr.f"
L36:
#line 705 "zstemr.f"
		if (wend < *m) {
#line 706 "zstemr.f"
		    if (iwork[iindbl + wend] == jblk) {
#line 707 "zstemr.f"
			++wend;
#line 708 "zstemr.f"
			goto L36;
#line 709 "zstemr.f"
		    }
#line 710 "zstemr.f"
		}
#line 711 "zstemr.f"
		if (wend < wbegin) {
#line 712 "zstemr.f"
		    ibegin = iend + 1;
#line 713 "zstemr.f"
		    goto L39;
#line 714 "zstemr.f"
		}
#line 716 "zstemr.f"
		offset = iwork[iindw + wbegin - 1] - 1;
#line 717 "zstemr.f"
		ifirst = iwork[iindw + wbegin - 1];
#line 718 "zstemr.f"
		ilast = iwork[iindw + wend - 1];
#line 719 "zstemr.f"
		rtol2 = eps * 4.;
#line 720 "zstemr.f"
		dlarrj_(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 
			1], &ifirst, &ilast, &rtol2, &offset, &w[wbegin], &
			work[inderr + wbegin - 1], &work[indwrk], &iwork[
			iindwk], &pivmin, &tnrm, &iinfo);
#line 726 "zstemr.f"
		ibegin = iend + 1;
#line 727 "zstemr.f"
		wbegin = wend + 1;
#line 728 "zstemr.f"
L39:
#line 728 "zstemr.f"
		;
#line 728 "zstemr.f"
	    }
#line 729 "zstemr.f"
	}

/*        If matrix was scaled, then rescale eigenvalues appropriately. */

#line 733 "zstemr.f"
	if (scale != 1.) {
#line 734 "zstemr.f"
	    d__1 = 1. / scale;
#line 734 "zstemr.f"
	    dscal_(m, &d__1, &w[1], &c__1);
#line 735 "zstemr.f"
	}
#line 736 "zstemr.f"
    }

/*     If eigenvalues are not in increasing order, then sort them, */
/*     possibly along with eigenvectors. */

#line 741 "zstemr.f"
    if (nsplit > 1 || *n == 2) {
#line 742 "zstemr.f"
	if (! wantz) {
#line 743 "zstemr.f"
	    dlasrt_("I", m, &w[1], &iinfo, (ftnlen)1);
#line 744 "zstemr.f"
	    if (iinfo != 0) {
#line 745 "zstemr.f"
		*info = 3;
#line 746 "zstemr.f"
		return 0;
#line 747 "zstemr.f"
	    }
#line 748 "zstemr.f"
	} else {
#line 749 "zstemr.f"
	    i__1 = *m - 1;
#line 749 "zstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 750 "zstemr.f"
		i__ = 0;
#line 751 "zstemr.f"
		tmp = w[j];
#line 752 "zstemr.f"
		i__2 = *m;
#line 752 "zstemr.f"
		for (jj = j + 1; jj <= i__2; ++jj) {
#line 753 "zstemr.f"
		    if (w[jj] < tmp) {
#line 754 "zstemr.f"
			i__ = jj;
#line 755 "zstemr.f"
			tmp = w[jj];
#line 756 "zstemr.f"
		    }
#line 757 "zstemr.f"
/* L50: */
#line 757 "zstemr.f"
		}
#line 758 "zstemr.f"
		if (i__ != 0) {
#line 759 "zstemr.f"
		    w[i__] = w[j];
#line 760 "zstemr.f"
		    w[j] = tmp;
#line 761 "zstemr.f"
		    if (wantz) {
#line 762 "zstemr.f"
			zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * 
				z_dim1 + 1], &c__1);
#line 763 "zstemr.f"
			itmp = isuppz[(i__ << 1) - 1];
#line 764 "zstemr.f"
			isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
#line 765 "zstemr.f"
			isuppz[(j << 1) - 1] = itmp;
#line 766 "zstemr.f"
			itmp = isuppz[i__ * 2];
#line 767 "zstemr.f"
			isuppz[i__ * 2] = isuppz[j * 2];
#line 768 "zstemr.f"
			isuppz[j * 2] = itmp;
#line 769 "zstemr.f"
		    }
#line 770 "zstemr.f"
		}
#line 771 "zstemr.f"
/* L60: */
#line 771 "zstemr.f"
	    }
#line 772 "zstemr.f"
	}
#line 773 "zstemr.f"
    }


#line 776 "zstemr.f"
    work[1] = (doublereal) lwmin;
#line 777 "zstemr.f"
    iwork[1] = liwmin;
#line 778 "zstemr.f"
    return 0;

/*     End of ZSTEMR */

} /* zstemr_ */

