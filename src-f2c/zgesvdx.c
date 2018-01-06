#line 1 "zgesvdx.f"
/* zgesvdx.f -- translated by f2c (version 20100827).
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

#line 1 "zgesvdx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> ZGESVDX computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGESVDX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesvdx
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesvdx
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesvdx
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE CGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, */
/*    $                    IL, IU, NS, S, U, LDU, VT, LDVT, WORK, */
/*    $                    LWORK, RWORK, IWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER          JOBU, JOBVT, RANGE */
/*      INTEGER            IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS */
/*      DOUBLE PRECISION   VL, VU */
/*     .. */
/*     .. Array Arguments .. */
/*      INTEGER            IWORK( * ) */
/*      DOUBLE PRECISION   S( * ), RWORK( * ) */
/*      COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*     $                   WORK( * ) */
/*     .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >  ZGESVDX computes the singular value decomposition (SVD) of a complex */
/* >  M-by-N matrix A, optionally computing the left and/or right singular */
/* >  vectors. The SVD is written */
/* > */
/* >      A = U * SIGMA * transpose(V) */
/* > */
/* >  where SIGMA is an M-by-N matrix which is zero except for its */
/* >  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and */
/* >  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA */
/* >  are the singular values of A; they are real and non-negative, and */
/* >  are returned in descending order.  The first min(m,n) columns of */
/* >  U and V are the left and right singular vectors of A. */
/* > */
/* >  ZGESVDX uses an eigenvalue problem for obtaining the SVD, which */
/* >  allows for the computation of a subset of singular values and */
/* >  vectors. See DBDSVDX for details. */
/* > */
/* >  Note that the routine returns V**T, not V. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          Specifies options for computing all or part of the matrix U: */
/* >          = 'V':  the first min(m,n) columns of U (the left singular */
/* >                  vectors) or as specified by RANGE are returned in */
/* >                  the array U; */
/* >          = 'N':  no columns of U (no left singular vectors) are */
/* >                  computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVT */
/* > \verbatim */
/* >          JOBVT is CHARACTER*1 */
/* >           Specifies options for computing all or part of the matrix */
/* >           V**T: */
/* >           = 'V':  the first min(m,n) rows of V**T (the right singular */
/* >                   vectors) or as specified by RANGE are returned in */
/* >                   the array VT; */
/* >           = 'N':  no rows of V**T (no right singular vectors) are */
/* >                   computed. */
/* > \endverbatim */
/* > */
/* > \param[in] RANGE */
/* > \verbatim */
/* >          RANGE is CHARACTER*1 */
/* >          = 'A': all singular values will be found. */
/* >          = 'V': all singular values in the half-open interval (VL,VU] */
/* >                 will be found. */
/* >          = 'I': the IL-th through IU-th singular values will be found. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the input matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the input matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the contents of A are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for singular values. VU > VL. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for singular values. VU > VL. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          smallest singular value to be returned. */
/* >          1 <= IL <= IU <= min(M,N), if min(M,N) > 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          largest singular value to be returned. */
/* >          1 <= IL <= IU <= min(M,N), if min(M,N) > 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* >          NS is INTEGER */
/* >          The total number of singular values found, */
/* >          0 <= NS <= min(M,N). */
/* >          If RANGE = 'A', NS = min(M,N); if RANGE = 'I', NS = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX*16 array, dimension (LDU,UCOL) */
/* >          If JOBU = 'V', U contains columns of U (the left singular */
/* >          vectors, stored columnwise) as specified by RANGE; if */
/* >          JOBU = 'N', U is not referenced. */
/* >          Note: The user must ensure that UCOL >= NS; if RANGE = 'V', */
/* >          the exact value of NS is not known in advance and an upper */
/* >          bound must be used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U.  LDU >= 1; if */
/* >          JOBU = 'V', LDU >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is COMPLEX*16 array, dimension (LDVT,N) */
/* >          If JOBVT = 'V', VT contains the rows of V**T (the right singular */
/* >          vectors, stored rowwise) as specified by RANGE; if JOBVT = 'N', */
/* >          VT is not referenced. */
/* >          Note: The user must ensure that LDVT >= NS; if RANGE = 'V', */
/* >          the exact value of NS is not known in advance and an upper */
/* >          bound must be used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >          The leading dimension of the array VT.  LDVT >= 1; if */
/* >          JOBVT = 'V', LDVT >= NS (see above). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK; */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          LWORK >= MAX(1,MIN(M,N)*(MIN(M,N)+4)) for the paths (see */
/* >          comments inside the code): */
/* >             - PATH 1  (M much larger than N) */
/* >             - PATH 1t (N much larger than M) */
/* >          LWORK >= MAX(1,MIN(M,N)*2+MAX(M,N)) for the other paths. */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
/* >          LRWORK >= MIN(M,N)*(MIN(M,N)*2+15*MIN(M,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (12*MIN(M,N)) */
/* >          If INFO = 0, the first NS elements of IWORK are zero. If INFO > 0, */
/* >          then IWORK contains the indices of the eigenvectors that failed */
/* >          to converge in DBDSVDX/DSTEVX. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >     INFO is INTEGER */
/* >           = 0:  successful exit */
/* >           < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >           > 0:  if INFO = i, then i eigenvectors failed to converge */
/* >                 in DBDSVDX/DSTEVX. */
/* >                 if INFO = N*2 + 1, an internal error occurred in */
/* >                 DBDSVDX */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complex16GEsing */

/*  ===================================================================== */
/* Subroutine */ int zgesvdx_(char *jobu, char *jobvt, char *range, integer *
	m, integer *n, doublecomplex *a, integer *lda, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, integer *ns, doublereal *s, 
	doublecomplex *u, integer *ldu, doublecomplex *vt, integer *ldvt, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	iwork, integer *info, ftnlen jobu_len, ftnlen jobvt_len, ftnlen 
	range_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1[2], 
	    i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, id, ie;
    static doublereal dum[1], eps;
    static integer iscl;
    static logical alls, inds;
    static integer ilqf;
    static doublereal anrm;
    static integer ierr, iqrf, itau;
    static char jobz[1];
    static logical vals;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iltgk, itemp, minmn, itaup, itauq, iutgk, itgkz, mnthr;
    static logical wantu;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), zgebrd_(integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum, abstol;
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zlascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
    static char rngtgk[1];
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static integer itempr;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), zlacpy_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int zunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen);
    static logical lquery, wantvt;
    extern /* Subroutine */ int zunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), dbdsvdx_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 327 "zgesvdx.f"
    /* Parameter adjustments */
#line 327 "zgesvdx.f"
    a_dim1 = *lda;
#line 327 "zgesvdx.f"
    a_offset = 1 + a_dim1;
#line 327 "zgesvdx.f"
    a -= a_offset;
#line 327 "zgesvdx.f"
    --s;
#line 327 "zgesvdx.f"
    u_dim1 = *ldu;
#line 327 "zgesvdx.f"
    u_offset = 1 + u_dim1;
#line 327 "zgesvdx.f"
    u -= u_offset;
#line 327 "zgesvdx.f"
    vt_dim1 = *ldvt;
#line 327 "zgesvdx.f"
    vt_offset = 1 + vt_dim1;
#line 327 "zgesvdx.f"
    vt -= vt_offset;
#line 327 "zgesvdx.f"
    --work;
#line 327 "zgesvdx.f"
    --rwork;
#line 327 "zgesvdx.f"
    --iwork;
#line 327 "zgesvdx.f"

#line 327 "zgesvdx.f"
    /* Function Body */
#line 327 "zgesvdx.f"
    *ns = 0;
#line 328 "zgesvdx.f"
    *info = 0;
#line 329 "zgesvdx.f"
    abstol = dlamch_("S", (ftnlen)1) * 2;
#line 330 "zgesvdx.f"
    lquery = *lwork == -1;
#line 331 "zgesvdx.f"
    minmn = min(*m,*n);
#line 333 "zgesvdx.f"
    wantu = lsame_(jobu, "V", (ftnlen)1, (ftnlen)1);
#line 334 "zgesvdx.f"
    wantvt = lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1);
#line 335 "zgesvdx.f"
    if (wantu || wantvt) {
#line 336 "zgesvdx.f"
	*(unsigned char *)jobz = 'V';
#line 337 "zgesvdx.f"
    } else {
#line 338 "zgesvdx.f"
	*(unsigned char *)jobz = 'N';
#line 339 "zgesvdx.f"
    }
#line 340 "zgesvdx.f"
    alls = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 341 "zgesvdx.f"
    vals = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 342 "zgesvdx.f"
    inds = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 344 "zgesvdx.f"
    *info = 0;
#line 345 "zgesvdx.f"
    if (! lsame_(jobu, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobu, "N", (
	    ftnlen)1, (ftnlen)1)) {
#line 347 "zgesvdx.f"
	*info = -1;
#line 348 "zgesvdx.f"
    } else if (! lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobvt, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 350 "zgesvdx.f"
	*info = -2;
#line 351 "zgesvdx.f"
    } else if (! (alls || vals || inds)) {
#line 352 "zgesvdx.f"
	*info = -3;
#line 353 "zgesvdx.f"
    } else if (*m < 0) {
#line 354 "zgesvdx.f"
	*info = -4;
#line 355 "zgesvdx.f"
    } else if (*n < 0) {
#line 356 "zgesvdx.f"
	*info = -5;
#line 357 "zgesvdx.f"
    } else if (*m > *lda) {
#line 358 "zgesvdx.f"
	*info = -7;
#line 359 "zgesvdx.f"
    } else if (minmn > 0) {
#line 360 "zgesvdx.f"
	if (vals) {
#line 361 "zgesvdx.f"
	    if (*vl < 0.) {
#line 362 "zgesvdx.f"
		*info = -8;
#line 363 "zgesvdx.f"
	    } else if (*vu <= *vl) {
#line 364 "zgesvdx.f"
		*info = -9;
#line 365 "zgesvdx.f"
	    }
#line 366 "zgesvdx.f"
	} else if (inds) {
#line 367 "zgesvdx.f"
	    if (*il < 1 || *il > max(1,minmn)) {
#line 368 "zgesvdx.f"
		*info = -10;
#line 369 "zgesvdx.f"
	    } else if (*iu < min(minmn,*il) || *iu > minmn) {
#line 370 "zgesvdx.f"
		*info = -11;
#line 371 "zgesvdx.f"
	    }
#line 372 "zgesvdx.f"
	}
#line 373 "zgesvdx.f"
	if (*info == 0) {
#line 374 "zgesvdx.f"
	    if (wantu && *ldu < *m) {
#line 375 "zgesvdx.f"
		*info = -15;
#line 376 "zgesvdx.f"
	    } else if (wantvt) {
#line 377 "zgesvdx.f"
		if (inds) {
#line 378 "zgesvdx.f"
		    if (*ldvt < *iu - *il + 1) {
#line 379 "zgesvdx.f"
			*info = -17;
#line 380 "zgesvdx.f"
		    }
#line 381 "zgesvdx.f"
		} else if (*ldvt < minmn) {
#line 382 "zgesvdx.f"
		    *info = -17;
#line 383 "zgesvdx.f"
		}
#line 384 "zgesvdx.f"
	    }
#line 385 "zgesvdx.f"
	}
#line 386 "zgesvdx.f"
    }

/*     Compute workspace */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 395 "zgesvdx.f"
    if (*info == 0) {
#line 396 "zgesvdx.f"
	minwrk = 1;
#line 397 "zgesvdx.f"
	maxwrk = 1;
#line 398 "zgesvdx.f"
	if (minmn > 0) {
#line 399 "zgesvdx.f"
	    if (*m >= *n) {
/* Writing concatenation */
#line 400 "zgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 400 "zgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 400 "zgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 400 "zgesvdx.f"
		mnthr = ilaenv_(&c__6, "ZGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 401 "zgesvdx.f"
		if (*m >= mnthr) {

/*                 Path 1 (M much larger than N) */

#line 405 "zgesvdx.f"
		    minwrk = *n * (*n + 5);
#line 406 "zgesvdx.f"
		    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 407 "zgesvdx.f"
		    i__2 = maxwrk, i__3 = *n * *n + (*n << 1) + (*n << 1) * 
			    ilaenv_(&c__1, "ZGEBRD", " ", n, n, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 407 "zgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 409 "zgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 410 "zgesvdx.f"
			i__2 = maxwrk, i__3 = *n * *n + (*n << 1) + *n * 
				ilaenv_(&c__1, "ZUNMQR", "LN", n, n, n, &c_n1,
				 (ftnlen)6, (ftnlen)2);
#line 410 "zgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 412 "zgesvdx.f"
		    }
#line 413 "zgesvdx.f"
		} else {

/*                 Path 2 (M at least N, but not much larger) */

#line 417 "zgesvdx.f"
		    minwrk = *n * 3 + *m;
#line 418 "zgesvdx.f"
		    maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 419 "zgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 420 "zgesvdx.f"
			i__2 = maxwrk, i__3 = (*n << 1) + *n * ilaenv_(&c__1, 
				"ZUNMQR", "LN", n, n, n, &c_n1, (ftnlen)6, (
				ftnlen)2);
#line 420 "zgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 422 "zgesvdx.f"
		    }
#line 423 "zgesvdx.f"
		}
#line 424 "zgesvdx.f"
	    } else {
/* Writing concatenation */
#line 425 "zgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 425 "zgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 425 "zgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 425 "zgesvdx.f"
		mnthr = ilaenv_(&c__6, "ZGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 426 "zgesvdx.f"
		if (*n >= mnthr) {

/*                 Path 1t (N much larger than M) */

#line 430 "zgesvdx.f"
		    minwrk = *m * (*m + 5);
#line 431 "zgesvdx.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 432 "zgesvdx.f"
		    i__2 = maxwrk, i__3 = *m * *m + (*m << 1) + (*m << 1) * 
			    ilaenv_(&c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 432 "zgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 434 "zgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 435 "zgesvdx.f"
			i__2 = maxwrk, i__3 = *m * *m + (*m << 1) + *m * 
				ilaenv_(&c__1, "ZUNMQR", "LN", m, m, m, &c_n1,
				 (ftnlen)6, (ftnlen)2);
#line 435 "zgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 437 "zgesvdx.f"
		    }
#line 438 "zgesvdx.f"
		} else {

/*                 Path 2t (N greater than M, but not much larger) */


#line 443 "zgesvdx.f"
		    minwrk = *m * 3 + *n;
#line 444 "zgesvdx.f"
		    maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "ZGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 445 "zgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 446 "zgesvdx.f"
			i__2 = maxwrk, i__3 = (*m << 1) + *m * ilaenv_(&c__1, 
				"ZUNMQR", "LN", m, m, m, &c_n1, (ftnlen)6, (
				ftnlen)2);
#line 446 "zgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 448 "zgesvdx.f"
		    }
#line 449 "zgesvdx.f"
		}
#line 450 "zgesvdx.f"
	    }
#line 451 "zgesvdx.f"
	}
#line 452 "zgesvdx.f"
	maxwrk = max(maxwrk,minwrk);
#line 453 "zgesvdx.f"
	d__1 = (doublereal) maxwrk;
#line 453 "zgesvdx.f"
	z__1.r = d__1, z__1.i = 0.;
#line 453 "zgesvdx.f"
	work[1].r = z__1.r, work[1].i = z__1.i;

#line 455 "zgesvdx.f"
	if (*lwork < minwrk && ! lquery) {
#line 456 "zgesvdx.f"
	    *info = -19;
#line 457 "zgesvdx.f"
	}
#line 458 "zgesvdx.f"
    }

#line 460 "zgesvdx.f"
    if (*info != 0) {
#line 461 "zgesvdx.f"
	i__2 = -(*info);
#line 461 "zgesvdx.f"
	xerbla_("ZGESVDX", &i__2, (ftnlen)7);
#line 462 "zgesvdx.f"
	return 0;
#line 463 "zgesvdx.f"
    } else if (lquery) {
#line 464 "zgesvdx.f"
	return 0;
#line 465 "zgesvdx.f"
    }

/*     Quick return if possible */

#line 469 "zgesvdx.f"
    if (*m == 0 || *n == 0) {
#line 470 "zgesvdx.f"
	return 0;
#line 471 "zgesvdx.f"
    }

/*     Set singular values indices accord to RANGE='A'. */

#line 475 "zgesvdx.f"
    if (alls) {
#line 476 "zgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 477 "zgesvdx.f"
	iltgk = 1;
#line 478 "zgesvdx.f"
	iutgk = min(*m,*n);
#line 479 "zgesvdx.f"
    } else if (inds) {
#line 480 "zgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 481 "zgesvdx.f"
	iltgk = *il;
#line 482 "zgesvdx.f"
	iutgk = *iu;
#line 483 "zgesvdx.f"
    } else {
#line 484 "zgesvdx.f"
	*(unsigned char *)rngtgk = 'V';
#line 485 "zgesvdx.f"
	iltgk = 0;
#line 486 "zgesvdx.f"
	iutgk = 0;
#line 487 "zgesvdx.f"
    }

/*     Get machine constants */

#line 491 "zgesvdx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 492 "zgesvdx.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 493 "zgesvdx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 497 "zgesvdx.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 498 "zgesvdx.f"
    iscl = 0;
#line 499 "zgesvdx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 500 "zgesvdx.f"
	iscl = 1;
#line 501 "zgesvdx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 502 "zgesvdx.f"
    } else if (anrm > bignum) {
#line 503 "zgesvdx.f"
	iscl = 1;
#line 504 "zgesvdx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 505 "zgesvdx.f"
    }

#line 507 "zgesvdx.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce A using the QR */
/*        decomposition. */

#line 513 "zgesvdx.f"
	if (*m >= mnthr) {

/*           Path 1 (M much larger than N): */
/*           A = Q * R = Q * ( QB * B * PB**T ) */
/*                     = Q * ( QB * ( UB * S * VB**T ) * PB**T ) */
/*           U = Q * QB * UB; V**T = VB**T * PB**T */

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 523 "zgesvdx.f"
	    itau = 1;
#line 524 "zgesvdx.f"
	    itemp = itau + *n;
#line 525 "zgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 525 "zgesvdx.f"
	    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);

/*           Copy R into WORK and bidiagonalize it: */
/*           (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB) */

#line 531 "zgesvdx.f"
	    iqrf = itemp;
#line 532 "zgesvdx.f"
	    itauq = itemp + *n * *n;
#line 533 "zgesvdx.f"
	    itaup = itauq + *n;
#line 534 "zgesvdx.f"
	    itemp = itaup + *n;
#line 535 "zgesvdx.f"
	    id = 1;
#line 536 "zgesvdx.f"
	    ie = id + *n;
#line 537 "zgesvdx.f"
	    itgkz = ie + *n;
#line 538 "zgesvdx.f"
	    zlacpy_("U", n, n, &a[a_offset], lda, &work[iqrf], n, (ftnlen)1);
#line 539 "zgesvdx.f"
	    i__2 = *n - 1;
#line 539 "zgesvdx.f"
	    i__3 = *n - 1;
#line 539 "zgesvdx.f"
	    zlaset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iqrf + 1], n, (
		    ftnlen)1);
#line 541 "zgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 541 "zgesvdx.f"
	    zgebrd_(n, n, &work[iqrf], n, &rwork[id], &rwork[ie], &work[itauq]
		    , &work[itaup], &work[itemp], &i__2, info);
#line 544 "zgesvdx.f"
	    itempr = itgkz + *n * ((*n << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*N*N+14*N) */

#line 549 "zgesvdx.f"
	    i__2 = *n << 1;
#line 549 "zgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, n, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 556 "zgesvdx.f"
	    if (wantu) {
#line 557 "zgesvdx.f"
		k = itgkz;
#line 558 "zgesvdx.f"
		i__2 = *ns;
#line 558 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 559 "zgesvdx.f"
		    i__3 = *n;
#line 559 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 560 "zgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 560 "zgesvdx.f"
			i__5 = k;
#line 560 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 560 "zgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 561 "zgesvdx.f"
			++k;
#line 562 "zgesvdx.f"
		    }
#line 563 "zgesvdx.f"
		    k += *n;
#line 564 "zgesvdx.f"
		}
#line 565 "zgesvdx.f"
		i__2 = *m - *n;
#line 565 "zgesvdx.f"
		zlaset_("A", &i__2, ns, &c_b1, &c_b1, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call ZUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 570 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 570 "zgesvdx.f"
		zunmbr_("Q", "L", "N", n, ns, n, &work[iqrf], n, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call ZUNMQR to compute Q*(QB*UB). */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 577 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 577 "zgesvdx.f"
		zunmqr_("L", "N", m, ns, n, &a[a_offset], lda, &work[itau], &
			u[u_offset], ldu, &work[itemp], &i__2, info, (ftnlen)
			1, (ftnlen)1);
#line 580 "zgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 584 "zgesvdx.f"
	    if (wantvt) {
#line 585 "zgesvdx.f"
		k = itgkz + *n;
#line 586 "zgesvdx.f"
		i__2 = *ns;
#line 586 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 587 "zgesvdx.f"
		    i__3 = *n;
#line 587 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 588 "zgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 588 "zgesvdx.f"
			i__5 = k;
#line 588 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 588 "zgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 589 "zgesvdx.f"
			++k;
#line 590 "zgesvdx.f"
		    }
#line 591 "zgesvdx.f"
		    k += *n;
#line 592 "zgesvdx.f"
		}

/*              Call ZUNMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 597 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 597 "zgesvdx.f"
		zunmbr_("P", "R", "C", ns, n, n, &work[iqrf], n, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 600 "zgesvdx.f"
	    }
#line 601 "zgesvdx.f"
	} else {

/*           Path 2 (M at least N, but not much larger) */
/*           Reduce A to bidiagonal form without QR decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB) */

#line 611 "zgesvdx.f"
	    itauq = 1;
#line 612 "zgesvdx.f"
	    itaup = itauq + *n;
#line 613 "zgesvdx.f"
	    itemp = itaup + *n;
#line 614 "zgesvdx.f"
	    id = 1;
#line 615 "zgesvdx.f"
	    ie = id + *n;
#line 616 "zgesvdx.f"
	    itgkz = ie + *n;
#line 617 "zgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 617 "zgesvdx.f"
	    zgebrd_(m, n, &a[a_offset], lda, &rwork[id], &rwork[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);
#line 620 "zgesvdx.f"
	    itempr = itgkz + *n * ((*n << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*N*N+14*N) */

#line 625 "zgesvdx.f"
	    i__2 = *n << 1;
#line 625 "zgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, n, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 632 "zgesvdx.f"
	    if (wantu) {
#line 633 "zgesvdx.f"
		k = itgkz;
#line 634 "zgesvdx.f"
		i__2 = *ns;
#line 634 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 635 "zgesvdx.f"
		    i__3 = *n;
#line 635 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 636 "zgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 636 "zgesvdx.f"
			i__5 = k;
#line 636 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 636 "zgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 637 "zgesvdx.f"
			++k;
#line 638 "zgesvdx.f"
		    }
#line 639 "zgesvdx.f"
		    k += *n;
#line 640 "zgesvdx.f"
		}
#line 641 "zgesvdx.f"
		i__2 = *m - *n;
#line 641 "zgesvdx.f"
		zlaset_("A", &i__2, ns, &c_b1, &c_b1, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call ZUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 646 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 646 "zgesvdx.f"
		zunmbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 649 "zgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 653 "zgesvdx.f"
	    if (wantvt) {
#line 654 "zgesvdx.f"
		k = itgkz + *n;
#line 655 "zgesvdx.f"
		i__2 = *ns;
#line 655 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 656 "zgesvdx.f"
		    i__3 = *n;
#line 656 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 657 "zgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 657 "zgesvdx.f"
			i__5 = k;
#line 657 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 657 "zgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 658 "zgesvdx.f"
			++k;
#line 659 "zgesvdx.f"
		    }
#line 660 "zgesvdx.f"
		    k += *n;
#line 661 "zgesvdx.f"
		}

/*              Call ZUNMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 666 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 666 "zgesvdx.f"
		zunmbr_("P", "R", "C", ns, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 669 "zgesvdx.f"
	    }
#line 670 "zgesvdx.f"
	}
#line 671 "zgesvdx.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce A using the LQ decomposition. */

#line 676 "zgesvdx.f"
	if (*n >= mnthr) {

/*           Path 1t (N much larger than M): */
/*           A = L * Q = ( QB * B * PB**T ) * Q */
/*                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q */
/*           U = QB * UB ; V**T = VB**T * PB**T * Q */

/*           Compute A=L*Q */
/*           (Workspace: need 2*M, prefer M+M*NB) */

#line 686 "zgesvdx.f"
	    itau = 1;
#line 687 "zgesvdx.f"
	    itemp = itau + *m;
#line 688 "zgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 688 "zgesvdx.f"
	    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);
/*           Copy L into WORK and bidiagonalize it: */
/*           (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB) */

#line 694 "zgesvdx.f"
	    ilqf = itemp;
#line 695 "zgesvdx.f"
	    itauq = ilqf + *m * *m;
#line 696 "zgesvdx.f"
	    itaup = itauq + *m;
#line 697 "zgesvdx.f"
	    itemp = itaup + *m;
#line 698 "zgesvdx.f"
	    id = 1;
#line 699 "zgesvdx.f"
	    ie = id + *m;
#line 700 "zgesvdx.f"
	    itgkz = ie + *m;
#line 701 "zgesvdx.f"
	    zlacpy_("L", m, m, &a[a_offset], lda, &work[ilqf], m, (ftnlen)1);
#line 702 "zgesvdx.f"
	    i__2 = *m - 1;
#line 702 "zgesvdx.f"
	    i__3 = *m - 1;
#line 702 "zgesvdx.f"
	    zlaset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ilqf + *m], m, (
		    ftnlen)1);
#line 704 "zgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 704 "zgesvdx.f"
	    zgebrd_(m, m, &work[ilqf], m, &rwork[id], &rwork[ie], &work[itauq]
		    , &work[itaup], &work[itemp], &i__2, info);
#line 707 "zgesvdx.f"
	    itempr = itgkz + *m * ((*m << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 712 "zgesvdx.f"
	    i__2 = *m << 1;
#line 712 "zgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, m, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 719 "zgesvdx.f"
	    if (wantu) {
#line 720 "zgesvdx.f"
		k = itgkz;
#line 721 "zgesvdx.f"
		i__2 = *ns;
#line 721 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 722 "zgesvdx.f"
		    i__3 = *m;
#line 722 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 723 "zgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 723 "zgesvdx.f"
			i__5 = k;
#line 723 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 723 "zgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 724 "zgesvdx.f"
			++k;
#line 725 "zgesvdx.f"
		    }
#line 726 "zgesvdx.f"
		    k += *m;
#line 727 "zgesvdx.f"
		}

/*              Call ZUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 732 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 732 "zgesvdx.f"
		zunmbr_("Q", "L", "N", m, ns, m, &work[ilqf], m, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 735 "zgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 739 "zgesvdx.f"
	    if (wantvt) {
#line 740 "zgesvdx.f"
		k = itgkz + *m;
#line 741 "zgesvdx.f"
		i__2 = *ns;
#line 741 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 742 "zgesvdx.f"
		    i__3 = *m;
#line 742 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 743 "zgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 743 "zgesvdx.f"
			i__5 = k;
#line 743 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 743 "zgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 744 "zgesvdx.f"
			++k;
#line 745 "zgesvdx.f"
		    }
#line 746 "zgesvdx.f"
		    k += *m;
#line 747 "zgesvdx.f"
		}
#line 748 "zgesvdx.f"
		i__2 = *n - *m;
#line 748 "zgesvdx.f"
		zlaset_("A", ns, &i__2, &c_b1, &c_b1, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call ZUNMBR to compute (VB**T)*(PB**T) */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 754 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 754 "zgesvdx.f"
		zunmbr_("P", "R", "C", ns, m, m, &work[ilqf], m, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call ZUNMLQ to compute ((VB**T)*(PB**T))*Q. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 761 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 761 "zgesvdx.f"
		zunmlq_("R", "N", ns, n, m, &a[a_offset], lda, &work[itau], &
			vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1);
#line 764 "zgesvdx.f"
	    }
#line 765 "zgesvdx.f"
	} else {

/*           Path 2t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB) */

#line 775 "zgesvdx.f"
	    itauq = 1;
#line 776 "zgesvdx.f"
	    itaup = itauq + *m;
#line 777 "zgesvdx.f"
	    itemp = itaup + *m;
#line 778 "zgesvdx.f"
	    id = 1;
#line 779 "zgesvdx.f"
	    ie = id + *m;
#line 780 "zgesvdx.f"
	    itgkz = ie + *m;
#line 781 "zgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 781 "zgesvdx.f"
	    zgebrd_(m, n, &a[a_offset], lda, &rwork[id], &rwork[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);
#line 784 "zgesvdx.f"
	    itempr = itgkz + *m * ((*m << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 789 "zgesvdx.f"
	    i__2 = *m << 1;
#line 789 "zgesvdx.f"
	    dbdsvdx_("L", jobz, rngtgk, m, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 796 "zgesvdx.f"
	    if (wantu) {
#line 797 "zgesvdx.f"
		k = itgkz;
#line 798 "zgesvdx.f"
		i__2 = *ns;
#line 798 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 799 "zgesvdx.f"
		    i__3 = *m;
#line 799 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 800 "zgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 800 "zgesvdx.f"
			i__5 = k;
#line 800 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 800 "zgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 801 "zgesvdx.f"
			++k;
#line 802 "zgesvdx.f"
		    }
#line 803 "zgesvdx.f"
		    k += *m;
#line 804 "zgesvdx.f"
		}

/*              Call ZUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 809 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 809 "zgesvdx.f"
		zunmbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, info, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 812 "zgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 816 "zgesvdx.f"
	    if (wantvt) {
#line 817 "zgesvdx.f"
		k = itgkz + *m;
#line 818 "zgesvdx.f"
		i__2 = *ns;
#line 818 "zgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 819 "zgesvdx.f"
		    i__3 = *m;
#line 819 "zgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 820 "zgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 820 "zgesvdx.f"
			i__5 = k;
#line 820 "zgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 820 "zgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 821 "zgesvdx.f"
			++k;
#line 822 "zgesvdx.f"
		    }
#line 823 "zgesvdx.f"
		    k += *m;
#line 824 "zgesvdx.f"
		}
#line 825 "zgesvdx.f"
		i__2 = *n - *m;
#line 825 "zgesvdx.f"
		zlaset_("A", ns, &i__2, &c_b1, &c_b1, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call ZUNMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 831 "zgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 831 "zgesvdx.f"
		zunmbr_("P", "R", "C", ns, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, 
			info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 834 "zgesvdx.f"
	    }
#line 835 "zgesvdx.f"
	}
#line 836 "zgesvdx.f"
    }

/*     Undo scaling if necessary */

#line 840 "zgesvdx.f"
    if (iscl == 1) {
#line 841 "zgesvdx.f"
	if (anrm > bignum) {
#line 841 "zgesvdx.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 841 "zgesvdx.f"
	}
#line 844 "zgesvdx.f"
	if (anrm < smlnum) {
#line 844 "zgesvdx.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 844 "zgesvdx.f"
	}
#line 847 "zgesvdx.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 851 "zgesvdx.f"
    d__1 = (doublereal) maxwrk;
#line 851 "zgesvdx.f"
    z__1.r = d__1, z__1.i = 0.;
#line 851 "zgesvdx.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

#line 853 "zgesvdx.f"
    return 0;

/*     End of ZGESVDX */

} /* zgesvdx_ */

