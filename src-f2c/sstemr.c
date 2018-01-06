#line 1 "sstemr.f"
/* sstemr.f -- translated by f2c (version 20100827).
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

#line 1 "sstemr.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = .003;

/* > \brief \b SSTEMR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstemr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstemr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstemr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, */
/*                          M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, */
/*                          IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE */
/*       LOGICAL            TRYRAC */
/*       INTEGER            IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N */
/*       REAL               VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       REAL               D( * ), E( * ), W( * ), WORK( * ) */
/*       REAL               Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEMR computes selected eigenvalues and, optionally, eigenvectors */
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
/* > 1.SSTEMR works only on machines which follow IEEE-754 */
/* > floating-point standard in their handling of infinities and NaNs. */
/* > This permits the use of efficient inner loops avoiding a check for */
/* > zero divisors. */
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
/* > */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
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
/* >          W is REAL array, dimension (N) */
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, max(1,M) ) */
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
/* >          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) ) */
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
/* >                if INFO = 2X, internal error in SLARRV. */
/* >                Here, the digit X = ABS( IINFO ) < 10, where IINFO is */
/* >                the nonzero error code returned by SLARRE or */
/* >                SLARRV, respectively. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup realOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int sstemr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, integer *m, doublereal *w, doublereal *z__, integer *ldz,
	 integer *nzc, integer *isuppz, logical *tryrac, doublereal *work, 
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
    static integer iindw, ilast, lwmin;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
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
	     integer *, integer *, integer *, ftnlen), slarre_(char *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static doublereal thresh;
    static integer iinspl, indwrk, ifirst, liwmin, nzcmin;
    static doublereal pivmin;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int slarrj_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *), slarrr_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer nsplit;
    extern /* Subroutine */ int slarrv_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    static doublereal smlnum;
    extern /* Subroutine */ int slasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static logical lquery, zquery;


/*  -- LAPACK computational routine (version 3.7.1) -- */
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

#line 377 "sstemr.f"
    /* Parameter adjustments */
#line 377 "sstemr.f"
    --d__;
#line 377 "sstemr.f"
    --e;
#line 377 "sstemr.f"
    --w;
#line 377 "sstemr.f"
    z_dim1 = *ldz;
#line 377 "sstemr.f"
    z_offset = 1 + z_dim1;
#line 377 "sstemr.f"
    z__ -= z_offset;
#line 377 "sstemr.f"
    --isuppz;
#line 377 "sstemr.f"
    --work;
#line 377 "sstemr.f"
    --iwork;
#line 377 "sstemr.f"

#line 377 "sstemr.f"
    /* Function Body */
#line 377 "sstemr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 378 "sstemr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 379 "sstemr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 380 "sstemr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 382 "sstemr.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 383 "sstemr.f"
    zquery = *nzc == -1;
/*     SSTEMR needs WORK of size 6*N, IWORK of size 3*N. */
/*     In addition, SLARRE needs WORK of size 6*N, IWORK of size 5*N. */
/*     Furthermore, SLARRV needs WORK of size 12*N, IWORK of size 7*N. */
#line 388 "sstemr.f"
    if (wantz) {
#line 389 "sstemr.f"
	lwmin = *n * 18;
#line 390 "sstemr.f"
	liwmin = *n * 10;
#line 391 "sstemr.f"
    } else {
/*        need less workspace if only the eigenvalues are wanted */
#line 393 "sstemr.f"
	lwmin = *n * 12;
#line 394 "sstemr.f"
	liwmin = *n << 3;
#line 395 "sstemr.f"
    }
#line 397 "sstemr.f"
    wl = 0.;
#line 398 "sstemr.f"
    wu = 0.;
#line 399 "sstemr.f"
    iil = 0;
#line 400 "sstemr.f"
    iiu = 0;
#line 401 "sstemr.f"
    nsplit = 0;
#line 403 "sstemr.f"
    if (valeig) {
/*        We do not reference VL, VU in the cases RANGE = 'I','A' */
/*        The interval (WL, WU] contains all the wanted eigenvalues. */
/*        It is either given by the user or computed in SLARRE. */
#line 407 "sstemr.f"
	wl = *vl;
#line 408 "sstemr.f"
	wu = *vu;
#line 409 "sstemr.f"
    } else if (indeig) {
/*        We do not reference IL, IU in the cases RANGE = 'V','A' */
#line 411 "sstemr.f"
	iil = *il;
#line 412 "sstemr.f"
	iiu = *iu;
#line 413 "sstemr.f"
    }

#line 415 "sstemr.f"
    *info = 0;
#line 416 "sstemr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 417 "sstemr.f"
	*info = -1;
#line 418 "sstemr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 419 "sstemr.f"
	*info = -2;
#line 420 "sstemr.f"
    } else if (*n < 0) {
#line 421 "sstemr.f"
	*info = -3;
#line 422 "sstemr.f"
    } else if (valeig && *n > 0 && wu <= wl) {
#line 423 "sstemr.f"
	*info = -7;
#line 424 "sstemr.f"
    } else if (indeig && (iil < 1 || iil > *n)) {
#line 425 "sstemr.f"
	*info = -8;
#line 426 "sstemr.f"
    } else if (indeig && (iiu < iil || iiu > *n)) {
#line 427 "sstemr.f"
	*info = -9;
#line 428 "sstemr.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 429 "sstemr.f"
	*info = -13;
#line 430 "sstemr.f"
    } else if (*lwork < lwmin && ! lquery) {
#line 431 "sstemr.f"
	*info = -17;
#line 432 "sstemr.f"
    } else if (*liwork < liwmin && ! lquery) {
#line 433 "sstemr.f"
	*info = -19;
#line 434 "sstemr.f"
    }

/*     Get machine constants. */

#line 438 "sstemr.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 439 "sstemr.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 440 "sstemr.f"
    smlnum = safmin / eps;
#line 441 "sstemr.f"
    bignum = 1. / smlnum;
#line 442 "sstemr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 443 "sstemr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 443 "sstemr.f"
    rmax = min(d__1,d__2);

#line 445 "sstemr.f"
    if (*info == 0) {
#line 446 "sstemr.f"
	work[1] = (doublereal) lwmin;
#line 447 "sstemr.f"
	iwork[1] = liwmin;

#line 449 "sstemr.f"
	if (wantz && alleig) {
#line 450 "sstemr.f"
	    nzcmin = *n;
#line 451 "sstemr.f"
	} else if (wantz && valeig) {
#line 452 "sstemr.f"
	    slarrc_("T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp, &
		    itmp2, info, (ftnlen)1);
#line 454 "sstemr.f"
	} else if (wantz && indeig) {
#line 455 "sstemr.f"
	    nzcmin = iiu - iil + 1;
#line 456 "sstemr.f"
	} else {
/*           WANTZ .EQ. FALSE. */
#line 458 "sstemr.f"
	    nzcmin = 0;
#line 459 "sstemr.f"
	}
#line 460 "sstemr.f"
	if (zquery && *info == 0) {
#line 461 "sstemr.f"
	    z__[z_dim1 + 1] = (doublereal) nzcmin;
#line 462 "sstemr.f"
	} else if (*nzc < nzcmin && ! zquery) {
#line 463 "sstemr.f"
	    *info = -14;
#line 464 "sstemr.f"
	}
#line 465 "sstemr.f"
    }
#line 467 "sstemr.f"
    if (*info != 0) {

#line 469 "sstemr.f"
	i__1 = -(*info);
#line 469 "sstemr.f"
	xerbla_("SSTEMR", &i__1, (ftnlen)6);

#line 471 "sstemr.f"
	return 0;
#line 472 "sstemr.f"
    } else if (lquery || zquery) {
#line 473 "sstemr.f"
	return 0;
#line 474 "sstemr.f"
    }

/*     Handle N = 0, 1, and 2 cases immediately */

#line 478 "sstemr.f"
    *m = 0;
#line 479 "sstemr.f"
    if (*n == 0) {
#line 479 "sstemr.f"
	return 0;
#line 479 "sstemr.f"
    }

#line 482 "sstemr.f"
    if (*n == 1) {
#line 483 "sstemr.f"
	if (alleig || indeig) {
#line 484 "sstemr.f"
	    *m = 1;
#line 485 "sstemr.f"
	    w[1] = d__[1];
#line 486 "sstemr.f"
	} else {
#line 487 "sstemr.f"
	    if (wl < d__[1] && wu >= d__[1]) {
#line 488 "sstemr.f"
		*m = 1;
#line 489 "sstemr.f"
		w[1] = d__[1];
#line 490 "sstemr.f"
	    }
#line 491 "sstemr.f"
	}
#line 492 "sstemr.f"
	if (wantz && ! zquery) {
#line 493 "sstemr.f"
	    z__[z_dim1 + 1] = 1.;
#line 494 "sstemr.f"
	    isuppz[1] = 1;
#line 495 "sstemr.f"
	    isuppz[2] = 1;
#line 496 "sstemr.f"
	}
#line 497 "sstemr.f"
	return 0;
#line 498 "sstemr.f"
    }

#line 500 "sstemr.f"
    if (*n == 2) {
#line 501 "sstemr.f"
	if (! wantz) {
#line 502 "sstemr.f"
	    slae2_(&d__[1], &e[1], &d__[2], &r1, &r2);
#line 503 "sstemr.f"
	} else if (wantz && ! zquery) {
#line 504 "sstemr.f"
	    slaev2_(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
#line 505 "sstemr.f"
	}
#line 506 "sstemr.f"
	if (alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1) {
#line 510 "sstemr.f"
	    ++(*m);
#line 511 "sstemr.f"
	    w[*m] = r2;
#line 512 "sstemr.f"
	    if (wantz && ! zquery) {
#line 513 "sstemr.f"
		z__[*m * z_dim1 + 1] = -sn;
#line 514 "sstemr.f"
		z__[*m * z_dim1 + 2] = cs;
/*              Note: At most one of SN and CS can be zero. */
#line 516 "sstemr.f"
		if (sn != 0.) {
#line 517 "sstemr.f"
		    if (cs != 0.) {
#line 518 "sstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 519 "sstemr.f"
			isuppz[*m * 2] = 2;
#line 520 "sstemr.f"
		    } else {
#line 521 "sstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 522 "sstemr.f"
			isuppz[*m * 2] = 1;
#line 523 "sstemr.f"
		    }
#line 524 "sstemr.f"
		} else {
#line 525 "sstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 526 "sstemr.f"
		    isuppz[*m * 2] = 2;
#line 527 "sstemr.f"
		}
#line 528 "sstemr.f"
	    }
#line 529 "sstemr.f"
	}
#line 530 "sstemr.f"
	if (alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2) {
#line 534 "sstemr.f"
	    ++(*m);
#line 535 "sstemr.f"
	    w[*m] = r1;
#line 536 "sstemr.f"
	    if (wantz && ! zquery) {
#line 537 "sstemr.f"
		z__[*m * z_dim1 + 1] = cs;
#line 538 "sstemr.f"
		z__[*m * z_dim1 + 2] = sn;
/*              Note: At most one of SN and CS can be zero. */
#line 540 "sstemr.f"
		if (sn != 0.) {
#line 541 "sstemr.f"
		    if (cs != 0.) {
#line 542 "sstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 543 "sstemr.f"
			isuppz[*m * 2] = 2;
#line 544 "sstemr.f"
		    } else {
#line 545 "sstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 546 "sstemr.f"
			isuppz[*m * 2] = 1;
#line 547 "sstemr.f"
		    }
#line 548 "sstemr.f"
		} else {
#line 549 "sstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 550 "sstemr.f"
		    isuppz[*m * 2] = 2;
#line 551 "sstemr.f"
		}
#line 552 "sstemr.f"
	    }
#line 553 "sstemr.f"
	}
#line 554 "sstemr.f"
    } else {
/*     Continue with general N */
#line 558 "sstemr.f"
	indgrs = 1;
#line 559 "sstemr.f"
	inderr = (*n << 1) + 1;
#line 560 "sstemr.f"
	indgp = *n * 3 + 1;
#line 561 "sstemr.f"
	indd = (*n << 2) + 1;
#line 562 "sstemr.f"
	inde2 = *n * 5 + 1;
#line 563 "sstemr.f"
	indwrk = *n * 6 + 1;

#line 565 "sstemr.f"
	iinspl = 1;
#line 566 "sstemr.f"
	iindbl = *n + 1;
#line 567 "sstemr.f"
	iindw = (*n << 1) + 1;
#line 568 "sstemr.f"
	iindwk = *n * 3 + 1;

/*        Scale matrix to allowable range, if necessary. */
/*        The allowable range is related to the PIVMIN parameter; see the */
/*        comments in SLARRD.  The preference for scaling small values */
/*        up is heuristic; we expect users' matrices not to be close to the */
/*        RMAX threshold. */

#line 576 "sstemr.f"
	scale = 1.;
#line 577 "sstemr.f"
	tnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 578 "sstemr.f"
	if (tnrm > 0. && tnrm < rmin) {
#line 579 "sstemr.f"
	    scale = rmin / tnrm;
#line 580 "sstemr.f"
	} else if (tnrm > rmax) {
#line 581 "sstemr.f"
	    scale = rmax / tnrm;
#line 582 "sstemr.f"
	}
#line 583 "sstemr.f"
	if (scale != 1.) {
#line 584 "sstemr.f"
	    sscal_(n, &scale, &d__[1], &c__1);
#line 585 "sstemr.f"
	    i__1 = *n - 1;
#line 585 "sstemr.f"
	    sscal_(&i__1, &scale, &e[1], &c__1);
#line 586 "sstemr.f"
	    tnrm *= scale;
#line 587 "sstemr.f"
	    if (valeig) {
/*              If eigenvalues in interval have to be found, */
/*              scale (WL, WU] accordingly */
#line 590 "sstemr.f"
		wl *= scale;
#line 591 "sstemr.f"
		wu *= scale;
#line 592 "sstemr.f"
	    }
#line 593 "sstemr.f"
	}

/*        Compute the desired eigenvalues of the tridiagonal after splitting */
/*        into smaller subblocks if the corresponding off-diagonal elements */
/*        are small */
/*        THRESH is the splitting parameter for SLARRE */
/*        A negative THRESH forces the old splitting criterion based on the */
/*        size of the off-diagonal. A positive THRESH switches to splitting */
/*        which preserves relative accuracy. */

#line 603 "sstemr.f"
	if (*tryrac) {
/*           Test whether the matrix warrants the more expensive relative approach. */
#line 605 "sstemr.f"
	    slarrr_(n, &d__[1], &e[1], &iinfo);
#line 606 "sstemr.f"
	} else {
/*           The user does not care about relative accurately eigenvalues */
#line 608 "sstemr.f"
	    iinfo = -1;
#line 609 "sstemr.f"
	}
/*        Set the splitting criterion */
#line 611 "sstemr.f"
	if (iinfo == 0) {
#line 612 "sstemr.f"
	    thresh = eps;
#line 613 "sstemr.f"
	} else {
#line 614 "sstemr.f"
	    thresh = -eps;
/*           relative accuracy is desired but T does not guarantee it */
#line 616 "sstemr.f"
	    *tryrac = FALSE_;
#line 617 "sstemr.f"
	}

#line 619 "sstemr.f"
	if (*tryrac) {
/*           Copy original diagonal, needed to guarantee relative accuracy */
#line 621 "sstemr.f"
	    scopy_(n, &d__[1], &c__1, &work[indd], &c__1);
#line 622 "sstemr.f"
	}
/*        Store the squares of the offdiagonal values of T */
#line 624 "sstemr.f"
	i__1 = *n - 1;
#line 624 "sstemr.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
#line 625 "sstemr.f"
	    d__1 = e[j];
#line 625 "sstemr.f"
	    work[inde2 + j - 1] = d__1 * d__1;
#line 626 "sstemr.f"
/* L5: */
#line 626 "sstemr.f"
	}
/*        Set the tolerance parameters for bisection */
#line 629 "sstemr.f"
	if (! wantz) {
/*           SLARRE computes the eigenvalues to full precision. */
#line 631 "sstemr.f"
	    rtol1 = eps * 4.;
#line 632 "sstemr.f"
	    rtol2 = eps * 4.;
#line 633 "sstemr.f"
	} else {
/*           SLARRE computes the eigenvalues to less than full precision. */
/*           SLARRV will refine the eigenvalue approximations, and we can */
/*           need less accurate initial bisection in SLARRE. */
/*           Note: these settings do only affect the subset case and SLARRE */
/* Computing MAX */
#line 638 "sstemr.f"
	    d__1 = sqrt(eps) * .05, d__2 = eps * 4.;
#line 638 "sstemr.f"
	    rtol1 = max(d__1,d__2);
/* Computing MAX */
#line 639 "sstemr.f"
	    d__1 = sqrt(eps) * .005, d__2 = eps * 4.;
#line 639 "sstemr.f"
	    rtol2 = max(d__1,d__2);
#line 640 "sstemr.f"
	}
#line 641 "sstemr.f"
	slarre_(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2], 
		&rtol1, &rtol2, &thresh, &nsplit, &iwork[iinspl], m, &w[1], &
		work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw], &
		work[indgrs], &pivmin, &work[indwrk], &iwork[iindwk], &iinfo, 
		(ftnlen)1);
#line 647 "sstemr.f"
	if (iinfo != 0) {
#line 648 "sstemr.f"
	    *info = abs(iinfo) + 10;
#line 649 "sstemr.f"
	    return 0;
#line 650 "sstemr.f"
	}
/*        Note that if RANGE .NE. 'V', SLARRE computes bounds on the desired */
/*        part of the spectrum. All desired eigenvalues are contained in */
/*        (WL,WU] */
#line 656 "sstemr.f"
	if (wantz) {

/*           Compute the desired eigenvectors corresponding to the computed */
/*           eigenvalues */

#line 661 "sstemr.f"
	    slarrv_(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &
		    c__1, m, &c_b18, &rtol1, &rtol2, &w[1], &work[inderr], &
		    work[indgp], &iwork[iindbl], &iwork[iindw], &work[indgrs],
		     &z__[z_offset], ldz, &isuppz[1], &work[indwrk], &iwork[
		    iindwk], &iinfo);
#line 667 "sstemr.f"
	    if (iinfo != 0) {
#line 668 "sstemr.f"
		*info = abs(iinfo) + 20;
#line 669 "sstemr.f"
		return 0;
#line 670 "sstemr.f"
	    }
#line 671 "sstemr.f"
	} else {
/*           SLARRE computes eigenvalues of the (shifted) root representation */
/*           SLARRV returns the eigenvalues of the unshifted matrix. */
/*           However, if the eigenvectors are not desired by the user, we need */
/*           to apply the corresponding shifts from SLARRE to obtain the */
/*           eigenvalues of the original matrix. */
#line 677 "sstemr.f"
	    i__1 = *m;
#line 677 "sstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 678 "sstemr.f"
		itmp = iwork[iindbl + j - 1];
#line 679 "sstemr.f"
		w[j] += e[iwork[iinspl + itmp - 1]];
#line 680 "sstemr.f"
/* L20: */
#line 680 "sstemr.f"
	    }
#line 681 "sstemr.f"
	}

#line 684 "sstemr.f"
	if (*tryrac) {
/*           Refine computed eigenvalues so that they are relatively accurate */
/*           with respect to the original matrix T. */
#line 687 "sstemr.f"
	    ibegin = 1;
#line 688 "sstemr.f"
	    wbegin = 1;
#line 689 "sstemr.f"
	    i__1 = iwork[iindbl + *m - 1];
#line 689 "sstemr.f"
	    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 690 "sstemr.f"
		iend = iwork[iinspl + jblk - 1];
#line 691 "sstemr.f"
		in = iend - ibegin + 1;
#line 692 "sstemr.f"
		wend = wbegin - 1;
/*              check if any eigenvalues have to be refined in this block */
#line 694 "sstemr.f"
L36:
#line 695 "sstemr.f"
		if (wend < *m) {
#line 696 "sstemr.f"
		    if (iwork[iindbl + wend] == jblk) {
#line 697 "sstemr.f"
			++wend;
#line 698 "sstemr.f"
			goto L36;
#line 699 "sstemr.f"
		    }
#line 700 "sstemr.f"
		}
#line 701 "sstemr.f"
		if (wend < wbegin) {
#line 702 "sstemr.f"
		    ibegin = iend + 1;
#line 703 "sstemr.f"
		    goto L39;
#line 704 "sstemr.f"
		}
#line 706 "sstemr.f"
		offset = iwork[iindw + wbegin - 1] - 1;
#line 707 "sstemr.f"
		ifirst = iwork[iindw + wbegin - 1];
#line 708 "sstemr.f"
		ilast = iwork[iindw + wend - 1];
#line 709 "sstemr.f"
		rtol2 = eps * 4.;
#line 710 "sstemr.f"
		slarrj_(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 
			1], &ifirst, &ilast, &rtol2, &offset, &w[wbegin], &
			work[inderr + wbegin - 1], &work[indwrk], &iwork[
			iindwk], &pivmin, &tnrm, &iinfo);
#line 716 "sstemr.f"
		ibegin = iend + 1;
#line 717 "sstemr.f"
		wbegin = wend + 1;
#line 718 "sstemr.f"
L39:
#line 718 "sstemr.f"
		;
#line 718 "sstemr.f"
	    }
#line 719 "sstemr.f"
	}

/*        If matrix was scaled, then rescale eigenvalues appropriately. */

#line 723 "sstemr.f"
	if (scale != 1.) {
#line 724 "sstemr.f"
	    d__1 = 1. / scale;
#line 724 "sstemr.f"
	    sscal_(m, &d__1, &w[1], &c__1);
#line 725 "sstemr.f"
	}
#line 726 "sstemr.f"
    }

/*     If eigenvalues are not in increasing order, then sort them, */
/*     possibly along with eigenvectors. */

#line 731 "sstemr.f"
    if (nsplit > 1 || *n == 2) {
#line 732 "sstemr.f"
	if (! wantz) {
#line 733 "sstemr.f"
	    slasrt_("I", m, &w[1], &iinfo, (ftnlen)1);
#line 734 "sstemr.f"
	    if (iinfo != 0) {
#line 735 "sstemr.f"
		*info = 3;
#line 736 "sstemr.f"
		return 0;
#line 737 "sstemr.f"
	    }
#line 738 "sstemr.f"
	} else {
#line 739 "sstemr.f"
	    i__1 = *m - 1;
#line 739 "sstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 740 "sstemr.f"
		i__ = 0;
#line 741 "sstemr.f"
		tmp = w[j];
#line 742 "sstemr.f"
		i__2 = *m;
#line 742 "sstemr.f"
		for (jj = j + 1; jj <= i__2; ++jj) {
#line 743 "sstemr.f"
		    if (w[jj] < tmp) {
#line 744 "sstemr.f"
			i__ = jj;
#line 745 "sstemr.f"
			tmp = w[jj];
#line 746 "sstemr.f"
		    }
#line 747 "sstemr.f"
/* L50: */
#line 747 "sstemr.f"
		}
#line 748 "sstemr.f"
		if (i__ != 0) {
#line 749 "sstemr.f"
		    w[i__] = w[j];
#line 750 "sstemr.f"
		    w[j] = tmp;
#line 751 "sstemr.f"
		    if (wantz) {
#line 752 "sstemr.f"
			sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * 
				z_dim1 + 1], &c__1);
#line 753 "sstemr.f"
			itmp = isuppz[(i__ << 1) - 1];
#line 754 "sstemr.f"
			isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
#line 755 "sstemr.f"
			isuppz[(j << 1) - 1] = itmp;
#line 756 "sstemr.f"
			itmp = isuppz[i__ * 2];
#line 757 "sstemr.f"
			isuppz[i__ * 2] = isuppz[j * 2];
#line 758 "sstemr.f"
			isuppz[j * 2] = itmp;
#line 759 "sstemr.f"
		    }
#line 760 "sstemr.f"
		}
#line 761 "sstemr.f"
/* L60: */
#line 761 "sstemr.f"
	    }
#line 762 "sstemr.f"
	}
#line 763 "sstemr.f"
    }


#line 766 "sstemr.f"
    work[1] = (doublereal) lwmin;
#line 767 "sstemr.f"
    iwork[1] = liwmin;
#line 768 "sstemr.f"
    return 0;

/*     End of SSTEMR */

} /* sstemr_ */

