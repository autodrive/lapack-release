#line 1 "dstemr.f"
/* dstemr.f -- translated by f2c (version 20100827).
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

#line 1 "dstemr.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = .001;

/* > \brief \b DSTEMR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstemr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstemr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstemr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, */
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
/*       DOUBLE PRECISION   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEMR computes selected eigenvalues and, optionally, eigenvectors */
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
/* > 1.DSTEMR works only on machines which follow IEEE-754 */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) ) */
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
/* >                if INFO = 2X, internal error in DLARRV. */
/* >                Here, the digit X = ABS( IINFO ) < 10, where IINFO is */
/* >                the nonzero error code returned by DLARRE or */
/* >                DLARRV, respectively. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2013 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dstemr_(char *jobz, char *range, integer *n, doublereal *
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
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer lwmin;
    static logical wantz;
    extern /* Subroutine */ int dlaev2_(doublereal *, doublereal *, 
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
	     integer *), dlarrv_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), dlasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal thresh;
    static integer iinspl, ifirst, indwrk, liwmin, nzcmin;
    static doublereal pivmin;
    static integer nsplit;
    static doublereal smlnum;
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

#line 370 "dstemr.f"
    /* Parameter adjustments */
#line 370 "dstemr.f"
    --d__;
#line 370 "dstemr.f"
    --e;
#line 370 "dstemr.f"
    --w;
#line 370 "dstemr.f"
    z_dim1 = *ldz;
#line 370 "dstemr.f"
    z_offset = 1 + z_dim1;
#line 370 "dstemr.f"
    z__ -= z_offset;
#line 370 "dstemr.f"
    --isuppz;
#line 370 "dstemr.f"
    --work;
#line 370 "dstemr.f"
    --iwork;
#line 370 "dstemr.f"

#line 370 "dstemr.f"
    /* Function Body */
#line 370 "dstemr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 371 "dstemr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 372 "dstemr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 373 "dstemr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 375 "dstemr.f"
    lquery = *lwork == -1 || *liwork == -1;
#line 376 "dstemr.f"
    zquery = *nzc == -1;
/*     DSTEMR needs WORK of size 6*N, IWORK of size 3*N. */
/*     In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N. */
/*     Furthermore, DLARRV needs WORK of size 12*N, IWORK of size 7*N. */
#line 381 "dstemr.f"
    if (wantz) {
#line 382 "dstemr.f"
	lwmin = *n * 18;
#line 383 "dstemr.f"
	liwmin = *n * 10;
#line 384 "dstemr.f"
    } else {
/*        need less workspace if only the eigenvalues are wanted */
#line 386 "dstemr.f"
	lwmin = *n * 12;
#line 387 "dstemr.f"
	liwmin = *n << 3;
#line 388 "dstemr.f"
    }
#line 390 "dstemr.f"
    wl = 0.;
#line 391 "dstemr.f"
    wu = 0.;
#line 392 "dstemr.f"
    iil = 0;
#line 393 "dstemr.f"
    iiu = 0;
#line 394 "dstemr.f"
    nsplit = 0;
#line 396 "dstemr.f"
    if (valeig) {
/*        We do not reference VL, VU in the cases RANGE = 'I','A' */
/*        The interval (WL, WU] contains all the wanted eigenvalues. */
/*        It is either given by the user or computed in DLARRE. */
#line 400 "dstemr.f"
	wl = *vl;
#line 401 "dstemr.f"
	wu = *vu;
#line 402 "dstemr.f"
    } else if (indeig) {
/*        We do not reference IL, IU in the cases RANGE = 'V','A' */
#line 404 "dstemr.f"
	iil = *il;
#line 405 "dstemr.f"
	iiu = *iu;
#line 406 "dstemr.f"
    }

#line 408 "dstemr.f"
    *info = 0;
#line 409 "dstemr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 410 "dstemr.f"
	*info = -1;
#line 411 "dstemr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 412 "dstemr.f"
	*info = -2;
#line 413 "dstemr.f"
    } else if (*n < 0) {
#line 414 "dstemr.f"
	*info = -3;
#line 415 "dstemr.f"
    } else if (valeig && *n > 0 && wu <= wl) {
#line 416 "dstemr.f"
	*info = -7;
#line 417 "dstemr.f"
    } else if (indeig && (iil < 1 || iil > *n)) {
#line 418 "dstemr.f"
	*info = -8;
#line 419 "dstemr.f"
    } else if (indeig && (iiu < iil || iiu > *n)) {
#line 420 "dstemr.f"
	*info = -9;
#line 421 "dstemr.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 422 "dstemr.f"
	*info = -13;
#line 423 "dstemr.f"
    } else if (*lwork < lwmin && ! lquery) {
#line 424 "dstemr.f"
	*info = -17;
#line 425 "dstemr.f"
    } else if (*liwork < liwmin && ! lquery) {
#line 426 "dstemr.f"
	*info = -19;
#line 427 "dstemr.f"
    }

/*     Get machine constants. */

#line 431 "dstemr.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 432 "dstemr.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 433 "dstemr.f"
    smlnum = safmin / eps;
#line 434 "dstemr.f"
    bignum = 1. / smlnum;
#line 435 "dstemr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 436 "dstemr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 436 "dstemr.f"
    rmax = min(d__1,d__2);

#line 438 "dstemr.f"
    if (*info == 0) {
#line 439 "dstemr.f"
	work[1] = (doublereal) lwmin;
#line 440 "dstemr.f"
	iwork[1] = liwmin;

#line 442 "dstemr.f"
	if (wantz && alleig) {
#line 443 "dstemr.f"
	    nzcmin = *n;
#line 444 "dstemr.f"
	} else if (wantz && valeig) {
#line 445 "dstemr.f"
	    dlarrc_("T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp, &
		    itmp2, info, (ftnlen)1);
#line 447 "dstemr.f"
	} else if (wantz && indeig) {
#line 448 "dstemr.f"
	    nzcmin = iiu - iil + 1;
#line 449 "dstemr.f"
	} else {
/*           WANTZ .EQ. FALSE. */
#line 451 "dstemr.f"
	    nzcmin = 0;
#line 452 "dstemr.f"
	}
#line 453 "dstemr.f"
	if (zquery && *info == 0) {
#line 454 "dstemr.f"
	    z__[z_dim1 + 1] = (doublereal) nzcmin;
#line 455 "dstemr.f"
	} else if (*nzc < nzcmin && ! zquery) {
#line 456 "dstemr.f"
	    *info = -14;
#line 457 "dstemr.f"
	}
#line 458 "dstemr.f"
    }
#line 460 "dstemr.f"
    if (*info != 0) {

#line 462 "dstemr.f"
	i__1 = -(*info);
#line 462 "dstemr.f"
	xerbla_("DSTEMR", &i__1, (ftnlen)6);

#line 464 "dstemr.f"
	return 0;
#line 465 "dstemr.f"
    } else if (lquery || zquery) {
#line 466 "dstemr.f"
	return 0;
#line 467 "dstemr.f"
    }

/*     Handle N = 0, 1, and 2 cases immediately */

#line 471 "dstemr.f"
    *m = 0;
#line 472 "dstemr.f"
    if (*n == 0) {
#line 472 "dstemr.f"
	return 0;
#line 472 "dstemr.f"
    }

#line 475 "dstemr.f"
    if (*n == 1) {
#line 476 "dstemr.f"
	if (alleig || indeig) {
#line 477 "dstemr.f"
	    *m = 1;
#line 478 "dstemr.f"
	    w[1] = d__[1];
#line 479 "dstemr.f"
	} else {
#line 480 "dstemr.f"
	    if (wl < d__[1] && wu >= d__[1]) {
#line 481 "dstemr.f"
		*m = 1;
#line 482 "dstemr.f"
		w[1] = d__[1];
#line 483 "dstemr.f"
	    }
#line 484 "dstemr.f"
	}
#line 485 "dstemr.f"
	if (wantz && ! zquery) {
#line 486 "dstemr.f"
	    z__[z_dim1 + 1] = 1.;
#line 487 "dstemr.f"
	    isuppz[1] = 1;
#line 488 "dstemr.f"
	    isuppz[2] = 1;
#line 489 "dstemr.f"
	}
#line 490 "dstemr.f"
	return 0;
#line 491 "dstemr.f"
    }

#line 493 "dstemr.f"
    if (*n == 2) {
#line 494 "dstemr.f"
	if (! wantz) {
#line 495 "dstemr.f"
	    dlae2_(&d__[1], &e[1], &d__[2], &r1, &r2);
#line 496 "dstemr.f"
	} else if (wantz && ! zquery) {
#line 497 "dstemr.f"
	    dlaev2_(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
#line 498 "dstemr.f"
	}
#line 499 "dstemr.f"
	if (alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1) {
#line 503 "dstemr.f"
	    ++(*m);
#line 504 "dstemr.f"
	    w[*m] = r2;
#line 505 "dstemr.f"
	    if (wantz && ! zquery) {
#line 506 "dstemr.f"
		z__[*m * z_dim1 + 1] = -sn;
#line 507 "dstemr.f"
		z__[*m * z_dim1 + 2] = cs;
/*              Note: At most one of SN and CS can be zero. */
#line 509 "dstemr.f"
		if (sn != 0.) {
#line 510 "dstemr.f"
		    if (cs != 0.) {
#line 511 "dstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 512 "dstemr.f"
			isuppz[*m * 2] = 2;
#line 513 "dstemr.f"
		    } else {
#line 514 "dstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 515 "dstemr.f"
			isuppz[*m * 2] = 1;
#line 516 "dstemr.f"
		    }
#line 517 "dstemr.f"
		} else {
#line 518 "dstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 519 "dstemr.f"
		    isuppz[*m * 2] = 2;
#line 520 "dstemr.f"
		}
#line 521 "dstemr.f"
	    }
#line 522 "dstemr.f"
	}
#line 523 "dstemr.f"
	if (alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2) {
#line 527 "dstemr.f"
	    ++(*m);
#line 528 "dstemr.f"
	    w[*m] = r1;
#line 529 "dstemr.f"
	    if (wantz && ! zquery) {
#line 530 "dstemr.f"
		z__[*m * z_dim1 + 1] = cs;
#line 531 "dstemr.f"
		z__[*m * z_dim1 + 2] = sn;
/*              Note: At most one of SN and CS can be zero. */
#line 533 "dstemr.f"
		if (sn != 0.) {
#line 534 "dstemr.f"
		    if (cs != 0.) {
#line 535 "dstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 536 "dstemr.f"
			isuppz[*m * 2] = 2;
#line 537 "dstemr.f"
		    } else {
#line 538 "dstemr.f"
			isuppz[(*m << 1) - 1] = 1;
#line 539 "dstemr.f"
			isuppz[*m * 2] = 1;
#line 540 "dstemr.f"
		    }
#line 541 "dstemr.f"
		} else {
#line 542 "dstemr.f"
		    isuppz[(*m << 1) - 1] = 2;
#line 543 "dstemr.f"
		    isuppz[*m * 2] = 2;
#line 544 "dstemr.f"
		}
#line 545 "dstemr.f"
	    }
#line 546 "dstemr.f"
	}
#line 548 "dstemr.f"
    } else {
/*     Continue with general N */
#line 552 "dstemr.f"
	indgrs = 1;
#line 553 "dstemr.f"
	inderr = (*n << 1) + 1;
#line 554 "dstemr.f"
	indgp = *n * 3 + 1;
#line 555 "dstemr.f"
	indd = (*n << 2) + 1;
#line 556 "dstemr.f"
	inde2 = *n * 5 + 1;
#line 557 "dstemr.f"
	indwrk = *n * 6 + 1;

#line 559 "dstemr.f"
	iinspl = 1;
#line 560 "dstemr.f"
	iindbl = *n + 1;
#line 561 "dstemr.f"
	iindw = (*n << 1) + 1;
#line 562 "dstemr.f"
	iindwk = *n * 3 + 1;

/*        Scale matrix to allowable range, if necessary. */
/*        The allowable range is related to the PIVMIN parameter; see the */
/*        comments in DLARRD.  The preference for scaling small values */
/*        up is heuristic; we expect users' matrices not to be close to the */
/*        RMAX threshold. */

#line 570 "dstemr.f"
	scale = 1.;
#line 571 "dstemr.f"
	tnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 572 "dstemr.f"
	if (tnrm > 0. && tnrm < rmin) {
#line 573 "dstemr.f"
	    scale = rmin / tnrm;
#line 574 "dstemr.f"
	} else if (tnrm > rmax) {
#line 575 "dstemr.f"
	    scale = rmax / tnrm;
#line 576 "dstemr.f"
	}
#line 577 "dstemr.f"
	if (scale != 1.) {
#line 578 "dstemr.f"
	    dscal_(n, &scale, &d__[1], &c__1);
#line 579 "dstemr.f"
	    i__1 = *n - 1;
#line 579 "dstemr.f"
	    dscal_(&i__1, &scale, &e[1], &c__1);
#line 580 "dstemr.f"
	    tnrm *= scale;
#line 581 "dstemr.f"
	    if (valeig) {
/*              If eigenvalues in interval have to be found, */
/*              scale (WL, WU] accordingly */
#line 584 "dstemr.f"
		wl *= scale;
#line 585 "dstemr.f"
		wu *= scale;
#line 586 "dstemr.f"
	    }
#line 587 "dstemr.f"
	}

/*        Compute the desired eigenvalues of the tridiagonal after splitting */
/*        into smaller subblocks if the corresponding off-diagonal elements */
/*        are small */
/*        THRESH is the splitting parameter for DLARRE */
/*        A negative THRESH forces the old splitting criterion based on the */
/*        size of the off-diagonal. A positive THRESH switches to splitting */
/*        which preserves relative accuracy. */

#line 597 "dstemr.f"
	if (*tryrac) {
/*           Test whether the matrix warrants the more expensive relative approach. */
#line 599 "dstemr.f"
	    dlarrr_(n, &d__[1], &e[1], &iinfo);
#line 600 "dstemr.f"
	} else {
/*           The user does not care about relative accurately eigenvalues */
#line 602 "dstemr.f"
	    iinfo = -1;
#line 603 "dstemr.f"
	}
/*        Set the splitting criterion */
#line 605 "dstemr.f"
	if (iinfo == 0) {
#line 606 "dstemr.f"
	    thresh = eps;
#line 607 "dstemr.f"
	} else {
#line 608 "dstemr.f"
	    thresh = -eps;
/*           relative accuracy is desired but T does not guarantee it */
#line 610 "dstemr.f"
	    *tryrac = FALSE_;
#line 611 "dstemr.f"
	}

#line 613 "dstemr.f"
	if (*tryrac) {
/*           Copy original diagonal, needed to guarantee relative accuracy */
#line 615 "dstemr.f"
	    dcopy_(n, &d__[1], &c__1, &work[indd], &c__1);
#line 616 "dstemr.f"
	}
/*        Store the squares of the offdiagonal values of T */
#line 618 "dstemr.f"
	i__1 = *n - 1;
#line 618 "dstemr.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
#line 619 "dstemr.f"
	    d__1 = e[j];
#line 619 "dstemr.f"
	    work[inde2 + j - 1] = d__1 * d__1;
#line 620 "dstemr.f"
/* L5: */
#line 620 "dstemr.f"
	}
/*        Set the tolerance parameters for bisection */
#line 623 "dstemr.f"
	if (! wantz) {
/*           DLARRE computes the eigenvalues to full precision. */
#line 625 "dstemr.f"
	    rtol1 = eps * 4.;
#line 626 "dstemr.f"
	    rtol2 = eps * 4.;
#line 627 "dstemr.f"
	} else {
/*           DLARRE computes the eigenvalues to less than full precision. */
/*           DLARRV will refine the eigenvalue approximations, and we can */
/*           need less accurate initial bisection in DLARRE. */
/*           Note: these settings do only affect the subset case and DLARRE */
#line 632 "dstemr.f"
	    rtol1 = sqrt(eps);
/* Computing MAX */
#line 633 "dstemr.f"
	    d__1 = sqrt(eps) * .005, d__2 = eps * 4.;
#line 633 "dstemr.f"
	    rtol2 = max(d__1,d__2);
#line 634 "dstemr.f"
	}
#line 635 "dstemr.f"
	dlarre_(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2], 
		&rtol1, &rtol2, &thresh, &nsplit, &iwork[iinspl], m, &w[1], &
		work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw], &
		work[indgrs], &pivmin, &work[indwrk], &iwork[iindwk], &iinfo, 
		(ftnlen)1);
#line 641 "dstemr.f"
	if (iinfo != 0) {
#line 642 "dstemr.f"
	    *info = abs(iinfo) + 10;
#line 643 "dstemr.f"
	    return 0;
#line 644 "dstemr.f"
	}
/*        Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired */
/*        part of the spectrum. All desired eigenvalues are contained in */
/*        (WL,WU] */
#line 650 "dstemr.f"
	if (wantz) {

/*           Compute the desired eigenvectors corresponding to the computed */
/*           eigenvalues */

#line 655 "dstemr.f"
	    dlarrv_(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &
		    c__1, m, &c_b18, &rtol1, &rtol2, &w[1], &work[inderr], &
		    work[indgp], &iwork[iindbl], &iwork[iindw], &work[indgrs],
		     &z__[z_offset], ldz, &isuppz[1], &work[indwrk], &iwork[
		    iindwk], &iinfo);
#line 661 "dstemr.f"
	    if (iinfo != 0) {
#line 662 "dstemr.f"
		*info = abs(iinfo) + 20;
#line 663 "dstemr.f"
		return 0;
#line 664 "dstemr.f"
	    }
#line 665 "dstemr.f"
	} else {
/*           DLARRE computes eigenvalues of the (shifted) root representation */
/*           DLARRV returns the eigenvalues of the unshifted matrix. */
/*           However, if the eigenvectors are not desired by the user, we need */
/*           to apply the corresponding shifts from DLARRE to obtain the */
/*           eigenvalues of the original matrix. */
#line 671 "dstemr.f"
	    i__1 = *m;
#line 671 "dstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 672 "dstemr.f"
		itmp = iwork[iindbl + j - 1];
#line 673 "dstemr.f"
		w[j] += e[iwork[iinspl + itmp - 1]];
#line 674 "dstemr.f"
/* L20: */
#line 674 "dstemr.f"
	    }
#line 675 "dstemr.f"
	}

#line 678 "dstemr.f"
	if (*tryrac) {
/*           Refine computed eigenvalues so that they are relatively accurate */
/*           with respect to the original matrix T. */
#line 681 "dstemr.f"
	    ibegin = 1;
#line 682 "dstemr.f"
	    wbegin = 1;
#line 683 "dstemr.f"
	    i__1 = iwork[iindbl + *m - 1];
#line 683 "dstemr.f"
	    for (jblk = 1; jblk <= i__1; ++jblk) {
#line 684 "dstemr.f"
		iend = iwork[iinspl + jblk - 1];
#line 685 "dstemr.f"
		in = iend - ibegin + 1;
#line 686 "dstemr.f"
		wend = wbegin - 1;
/*              check if any eigenvalues have to be refined in this block */
#line 688 "dstemr.f"
L36:
#line 689 "dstemr.f"
		if (wend < *m) {
#line 690 "dstemr.f"
		    if (iwork[iindbl + wend] == jblk) {
#line 691 "dstemr.f"
			++wend;
#line 692 "dstemr.f"
			goto L36;
#line 693 "dstemr.f"
		    }
#line 694 "dstemr.f"
		}
#line 695 "dstemr.f"
		if (wend < wbegin) {
#line 696 "dstemr.f"
		    ibegin = iend + 1;
#line 697 "dstemr.f"
		    goto L39;
#line 698 "dstemr.f"
		}
#line 700 "dstemr.f"
		offset = iwork[iindw + wbegin - 1] - 1;
#line 701 "dstemr.f"
		ifirst = iwork[iindw + wbegin - 1];
#line 702 "dstemr.f"
		ilast = iwork[iindw + wend - 1];
#line 703 "dstemr.f"
		rtol2 = eps * 4.;
#line 704 "dstemr.f"
		dlarrj_(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 
			1], &ifirst, &ilast, &rtol2, &offset, &w[wbegin], &
			work[inderr + wbegin - 1], &work[indwrk], &iwork[
			iindwk], &pivmin, &tnrm, &iinfo);
#line 710 "dstemr.f"
		ibegin = iend + 1;
#line 711 "dstemr.f"
		wbegin = wend + 1;
#line 712 "dstemr.f"
L39:
#line 712 "dstemr.f"
		;
#line 712 "dstemr.f"
	    }
#line 713 "dstemr.f"
	}

/*        If matrix was scaled, then rescale eigenvalues appropriately. */

#line 717 "dstemr.f"
	if (scale != 1.) {
#line 718 "dstemr.f"
	    d__1 = 1. / scale;
#line 718 "dstemr.f"
	    dscal_(m, &d__1, &w[1], &c__1);
#line 719 "dstemr.f"
	}
#line 721 "dstemr.f"
    }

/*     If eigenvalues are not in increasing order, then sort them, */
/*     possibly along with eigenvectors. */

#line 727 "dstemr.f"
    if (nsplit > 1 || *n == 2) {
#line 728 "dstemr.f"
	if (! wantz) {
#line 729 "dstemr.f"
	    dlasrt_("I", m, &w[1], &iinfo, (ftnlen)1);
#line 730 "dstemr.f"
	    if (iinfo != 0) {
#line 731 "dstemr.f"
		*info = 3;
#line 732 "dstemr.f"
		return 0;
#line 733 "dstemr.f"
	    }
#line 734 "dstemr.f"
	} else {
#line 735 "dstemr.f"
	    i__1 = *m - 1;
#line 735 "dstemr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 736 "dstemr.f"
		i__ = 0;
#line 737 "dstemr.f"
		tmp = w[j];
#line 738 "dstemr.f"
		i__2 = *m;
#line 738 "dstemr.f"
		for (jj = j + 1; jj <= i__2; ++jj) {
#line 739 "dstemr.f"
		    if (w[jj] < tmp) {
#line 740 "dstemr.f"
			i__ = jj;
#line 741 "dstemr.f"
			tmp = w[jj];
#line 742 "dstemr.f"
		    }
#line 743 "dstemr.f"
/* L50: */
#line 743 "dstemr.f"
		}
#line 744 "dstemr.f"
		if (i__ != 0) {
#line 745 "dstemr.f"
		    w[i__] = w[j];
#line 746 "dstemr.f"
		    w[j] = tmp;
#line 747 "dstemr.f"
		    if (wantz) {
#line 748 "dstemr.f"
			dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * 
				z_dim1 + 1], &c__1);
#line 749 "dstemr.f"
			itmp = isuppz[(i__ << 1) - 1];
#line 750 "dstemr.f"
			isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
#line 751 "dstemr.f"
			isuppz[(j << 1) - 1] = itmp;
#line 752 "dstemr.f"
			itmp = isuppz[i__ * 2];
#line 753 "dstemr.f"
			isuppz[i__ * 2] = isuppz[j * 2];
#line 754 "dstemr.f"
			isuppz[j * 2] = itmp;
#line 755 "dstemr.f"
		    }
#line 756 "dstemr.f"
		}
#line 757 "dstemr.f"
/* L60: */
#line 757 "dstemr.f"
	    }
#line 758 "dstemr.f"
	}
#line 759 "dstemr.f"
    }


#line 762 "dstemr.f"
    work[1] = (doublereal) lwmin;
#line 763 "dstemr.f"
    iwork[1] = liwmin;
#line 764 "dstemr.f"
    return 0;

/*     End of DSTEMR */

} /* dstemr_ */

