#line 1 "cgesvdx.f"
/* cgesvdx.f -- translated by f2c (version 20100827).
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

#line 1 "cgesvdx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> CGESVDX computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESVDX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvdx
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvdx
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvdx
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
/*      REAL               VL, VU */
/*     .. */
/*     .. Array Arguments .. */
/*     INTEGER            IWORK( * ) */
/*     REAL               S( * ), RWORK( * ) */
/*     COMPLEX            A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/*    $                   WORK( * ) */
/*     .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >  CGESVDX computes the singular value decomposition (SVD) of a complex */
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
/* >  CGESVDX uses an eigenvalue problem for obtaining the SVD, which */
/* >  allows for the computation of a subset of singular values and */
/* >  vectors. See SBDSVDX for details. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          VL is REAL */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for singular values. VU > VL. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
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
/* >          S is REAL array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX array, dimension (LDU,UCOL) */
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
/* >          VT is COMPLEX array, dimension (LDVT,N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* >          LRWORK >= MIN(M,N)*(MIN(M,N)*2+15*MIN(M,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (12*MIN(M,N)) */
/* >          If INFO = 0, the first NS elements of IWORK are zero. If INFO > 0, */
/* >          then IWORK contains the indices of the eigenvectors that failed */
/* >          to converge in SBDSVDX/SSTEVX. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >     INFO is INTEGER */
/* >           = 0:  successful exit */
/* >           < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >           > 0:  if INFO = i, then i eigenvectors failed to converge */
/* >                 in SBDSVDX/SSTEVX. */
/* >                 if INFO = N*2 + 1, an internal error occurred in */
/* >                 SBDSVDX */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexGEsing */

/*  ===================================================================== */
/* Subroutine */ int cgesvdx_(char *jobu, char *jobvt, char *range, integer *
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
    extern /* Subroutine */ int cgebrd_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), cgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), clacpy_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal abstol;
    extern /* Subroutine */ int cunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen);
    static char rngtgk[1];
    extern /* Subroutine */ int cunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer itempr;
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static logical lquery, wantvt;
    extern /* Subroutine */ int sbdsvdx_(char *, char *, char *, integer *, 
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

#line 327 "cgesvdx.f"
    /* Parameter adjustments */
#line 327 "cgesvdx.f"
    a_dim1 = *lda;
#line 327 "cgesvdx.f"
    a_offset = 1 + a_dim1;
#line 327 "cgesvdx.f"
    a -= a_offset;
#line 327 "cgesvdx.f"
    --s;
#line 327 "cgesvdx.f"
    u_dim1 = *ldu;
#line 327 "cgesvdx.f"
    u_offset = 1 + u_dim1;
#line 327 "cgesvdx.f"
    u -= u_offset;
#line 327 "cgesvdx.f"
    vt_dim1 = *ldvt;
#line 327 "cgesvdx.f"
    vt_offset = 1 + vt_dim1;
#line 327 "cgesvdx.f"
    vt -= vt_offset;
#line 327 "cgesvdx.f"
    --work;
#line 327 "cgesvdx.f"
    --rwork;
#line 327 "cgesvdx.f"
    --iwork;
#line 327 "cgesvdx.f"

#line 327 "cgesvdx.f"
    /* Function Body */
#line 327 "cgesvdx.f"
    *ns = 0;
#line 328 "cgesvdx.f"
    *info = 0;
#line 329 "cgesvdx.f"
    abstol = slamch_("S", (ftnlen)1) * 2;
#line 330 "cgesvdx.f"
    lquery = *lwork == -1;
#line 331 "cgesvdx.f"
    minmn = min(*m,*n);
#line 333 "cgesvdx.f"
    wantu = lsame_(jobu, "V", (ftnlen)1, (ftnlen)1);
#line 334 "cgesvdx.f"
    wantvt = lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1);
#line 335 "cgesvdx.f"
    if (wantu || wantvt) {
#line 336 "cgesvdx.f"
	*(unsigned char *)jobz = 'V';
#line 337 "cgesvdx.f"
    } else {
#line 338 "cgesvdx.f"
	*(unsigned char *)jobz = 'N';
#line 339 "cgesvdx.f"
    }
#line 340 "cgesvdx.f"
    alls = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 341 "cgesvdx.f"
    vals = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 342 "cgesvdx.f"
    inds = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 344 "cgesvdx.f"
    *info = 0;
#line 345 "cgesvdx.f"
    if (! lsame_(jobu, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobu, "N", (
	    ftnlen)1, (ftnlen)1)) {
#line 347 "cgesvdx.f"
	*info = -1;
#line 348 "cgesvdx.f"
    } else if (! lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobvt, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 350 "cgesvdx.f"
	*info = -2;
#line 351 "cgesvdx.f"
    } else if (! (alls || vals || inds)) {
#line 352 "cgesvdx.f"
	*info = -3;
#line 353 "cgesvdx.f"
    } else if (*m < 0) {
#line 354 "cgesvdx.f"
	*info = -4;
#line 355 "cgesvdx.f"
    } else if (*n < 0) {
#line 356 "cgesvdx.f"
	*info = -5;
#line 357 "cgesvdx.f"
    } else if (*m > *lda) {
#line 358 "cgesvdx.f"
	*info = -7;
#line 359 "cgesvdx.f"
    } else if (minmn > 0) {
#line 360 "cgesvdx.f"
	if (vals) {
#line 361 "cgesvdx.f"
	    if (*vl < 0.) {
#line 362 "cgesvdx.f"
		*info = -8;
#line 363 "cgesvdx.f"
	    } else if (*vu <= *vl) {
#line 364 "cgesvdx.f"
		*info = -9;
#line 365 "cgesvdx.f"
	    }
#line 366 "cgesvdx.f"
	} else if (inds) {
#line 367 "cgesvdx.f"
	    if (*il < 1 || *il > max(1,minmn)) {
#line 368 "cgesvdx.f"
		*info = -10;
#line 369 "cgesvdx.f"
	    } else if (*iu < min(minmn,*il) || *iu > minmn) {
#line 370 "cgesvdx.f"
		*info = -11;
#line 371 "cgesvdx.f"
	    }
#line 372 "cgesvdx.f"
	}
#line 373 "cgesvdx.f"
	if (*info == 0) {
#line 374 "cgesvdx.f"
	    if (wantu && *ldu < *m) {
#line 375 "cgesvdx.f"
		*info = -15;
#line 376 "cgesvdx.f"
	    } else if (wantvt) {
#line 377 "cgesvdx.f"
		if (inds) {
#line 378 "cgesvdx.f"
		    if (*ldvt < *iu - *il + 1) {
#line 379 "cgesvdx.f"
			*info = -17;
#line 380 "cgesvdx.f"
		    }
#line 381 "cgesvdx.f"
		} else if (*ldvt < minmn) {
#line 382 "cgesvdx.f"
		    *info = -17;
#line 383 "cgesvdx.f"
		}
#line 384 "cgesvdx.f"
	    }
#line 385 "cgesvdx.f"
	}
#line 386 "cgesvdx.f"
    }

/*     Compute workspace */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 395 "cgesvdx.f"
    if (*info == 0) {
#line 396 "cgesvdx.f"
	minwrk = 1;
#line 397 "cgesvdx.f"
	maxwrk = 1;
#line 398 "cgesvdx.f"
	if (minmn > 0) {
#line 399 "cgesvdx.f"
	    if (*m >= *n) {
/* Writing concatenation */
#line 400 "cgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 400 "cgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 400 "cgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 400 "cgesvdx.f"
		mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 401 "cgesvdx.f"
		if (*m >= mnthr) {

/*                 Path 1 (M much larger than N) */

#line 405 "cgesvdx.f"
		    minwrk = *n * (*n + 5);
#line 406 "cgesvdx.f"
		    maxwrk = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 407 "cgesvdx.f"
		    i__2 = maxwrk, i__3 = *n * *n + (*n << 1) + (*n << 1) * 
			    ilaenv_(&c__1, "CGEBRD", " ", n, n, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 407 "cgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 409 "cgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 410 "cgesvdx.f"
			i__2 = maxwrk, i__3 = *n * *n + (*n << 1) + *n * 
				ilaenv_(&c__1, "CUNMQR", "LN", n, n, n, &c_n1,
				 (ftnlen)6, (ftnlen)2);
#line 410 "cgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 412 "cgesvdx.f"
		    }
#line 413 "cgesvdx.f"
		} else {

/*                 Path 2 (M at least N, but not much larger) */

#line 417 "cgesvdx.f"
		    minwrk = *n * 3 + *m;
#line 418 "cgesvdx.f"
		    maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 419 "cgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 420 "cgesvdx.f"
			i__2 = maxwrk, i__3 = (*n << 1) + *n * ilaenv_(&c__1, 
				"CUNMQR", "LN", n, n, n, &c_n1, (ftnlen)6, (
				ftnlen)2);
#line 420 "cgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 422 "cgesvdx.f"
		    }
#line 423 "cgesvdx.f"
		}
#line 424 "cgesvdx.f"
	    } else {
/* Writing concatenation */
#line 425 "cgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 425 "cgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 425 "cgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 425 "cgesvdx.f"
		mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 426 "cgesvdx.f"
		if (*n >= mnthr) {

/*                 Path 1t (N much larger than M) */

#line 430 "cgesvdx.f"
		    minwrk = *m * (*m + 5);
#line 431 "cgesvdx.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "CGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 432 "cgesvdx.f"
		    i__2 = maxwrk, i__3 = *m * *m + (*m << 1) + (*m << 1) * 
			    ilaenv_(&c__1, "CGEBRD", " ", m, m, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 432 "cgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 434 "cgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 435 "cgesvdx.f"
			i__2 = maxwrk, i__3 = *m * *m + (*m << 1) + *m * 
				ilaenv_(&c__1, "CUNMQR", "LN", m, m, m, &c_n1,
				 (ftnlen)6, (ftnlen)2);
#line 435 "cgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 437 "cgesvdx.f"
		    }
#line 438 "cgesvdx.f"
		} else {

/*                 Path 2t (N greater than M, but not much larger) */


#line 443 "cgesvdx.f"
		    minwrk = *m * 3 + *n;
#line 444 "cgesvdx.f"
		    maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 445 "cgesvdx.f"
		    if (wantu || wantvt) {
/* Computing MAX */
#line 446 "cgesvdx.f"
			i__2 = maxwrk, i__3 = (*m << 1) + *m * ilaenv_(&c__1, 
				"CUNMQR", "LN", m, m, m, &c_n1, (ftnlen)6, (
				ftnlen)2);
#line 446 "cgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 448 "cgesvdx.f"
		    }
#line 449 "cgesvdx.f"
		}
#line 450 "cgesvdx.f"
	    }
#line 451 "cgesvdx.f"
	}
#line 452 "cgesvdx.f"
	maxwrk = max(maxwrk,minwrk);
#line 453 "cgesvdx.f"
	d__1 = (doublereal) maxwrk;
#line 453 "cgesvdx.f"
	z__1.r = d__1, z__1.i = 0.;
#line 453 "cgesvdx.f"
	work[1].r = z__1.r, work[1].i = z__1.i;

#line 455 "cgesvdx.f"
	if (*lwork < minwrk && ! lquery) {
#line 456 "cgesvdx.f"
	    *info = -19;
#line 457 "cgesvdx.f"
	}
#line 458 "cgesvdx.f"
    }

#line 460 "cgesvdx.f"
    if (*info != 0) {
#line 461 "cgesvdx.f"
	i__2 = -(*info);
#line 461 "cgesvdx.f"
	xerbla_("CGESVDX", &i__2, (ftnlen)7);
#line 462 "cgesvdx.f"
	return 0;
#line 463 "cgesvdx.f"
    } else if (lquery) {
#line 464 "cgesvdx.f"
	return 0;
#line 465 "cgesvdx.f"
    }

/*     Quick return if possible */

#line 469 "cgesvdx.f"
    if (*m == 0 || *n == 0) {
#line 470 "cgesvdx.f"
	return 0;
#line 471 "cgesvdx.f"
    }

/*     Set singular values indices accord to RANGE='A'. */

#line 475 "cgesvdx.f"
    if (alls) {
#line 476 "cgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 477 "cgesvdx.f"
	iltgk = 1;
#line 478 "cgesvdx.f"
	iutgk = min(*m,*n);
#line 479 "cgesvdx.f"
    } else if (inds) {
#line 480 "cgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 481 "cgesvdx.f"
	iltgk = *il;
#line 482 "cgesvdx.f"
	iutgk = *iu;
#line 483 "cgesvdx.f"
    } else {
#line 484 "cgesvdx.f"
	*(unsigned char *)rngtgk = 'V';
#line 485 "cgesvdx.f"
	iltgk = 0;
#line 486 "cgesvdx.f"
	iutgk = 0;
#line 487 "cgesvdx.f"
    }

/*     Get machine constants */

#line 491 "cgesvdx.f"
    eps = slamch_("P", (ftnlen)1);
#line 492 "cgesvdx.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 493 "cgesvdx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 497 "cgesvdx.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 498 "cgesvdx.f"
    iscl = 0;
#line 499 "cgesvdx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 500 "cgesvdx.f"
	iscl = 1;
#line 501 "cgesvdx.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 502 "cgesvdx.f"
    } else if (anrm > bignum) {
#line 503 "cgesvdx.f"
	iscl = 1;
#line 504 "cgesvdx.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 505 "cgesvdx.f"
    }

#line 507 "cgesvdx.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce A using the QR */
/*        decomposition. */

#line 513 "cgesvdx.f"
	if (*m >= mnthr) {

/*           Path 1 (M much larger than N): */
/*           A = Q * R = Q * ( QB * B * PB**T ) */
/*                     = Q * ( QB * ( UB * S * VB**T ) * PB**T ) */
/*           U = Q * QB * UB; V**T = VB**T * PB**T */

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 523 "cgesvdx.f"
	    itau = 1;
#line 524 "cgesvdx.f"
	    itemp = itau + *n;
#line 525 "cgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 525 "cgesvdx.f"
	    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);

/*           Copy R into WORK and bidiagonalize it: */
/*           (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB) */

#line 531 "cgesvdx.f"
	    iqrf = itemp;
#line 532 "cgesvdx.f"
	    itauq = itemp + *n * *n;
#line 533 "cgesvdx.f"
	    itaup = itauq + *n;
#line 534 "cgesvdx.f"
	    itemp = itaup + *n;
#line 535 "cgesvdx.f"
	    id = 1;
#line 536 "cgesvdx.f"
	    ie = id + *n;
#line 537 "cgesvdx.f"
	    itgkz = ie + *n;
#line 538 "cgesvdx.f"
	    clacpy_("U", n, n, &a[a_offset], lda, &work[iqrf], n, (ftnlen)1);
#line 539 "cgesvdx.f"
	    i__2 = *n - 1;
#line 539 "cgesvdx.f"
	    i__3 = *n - 1;
#line 539 "cgesvdx.f"
	    claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iqrf + 1], n, (
		    ftnlen)1);
#line 541 "cgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 541 "cgesvdx.f"
	    cgebrd_(n, n, &work[iqrf], n, &rwork[id], &rwork[ie], &work[itauq]
		    , &work[itaup], &work[itemp], &i__2, info);
#line 544 "cgesvdx.f"
	    itempr = itgkz + *n * ((*n << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*N*N+14*N) */

#line 549 "cgesvdx.f"
	    i__2 = *n << 1;
#line 549 "cgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, n, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 556 "cgesvdx.f"
	    if (wantu) {
#line 557 "cgesvdx.f"
		k = itgkz;
#line 558 "cgesvdx.f"
		i__2 = *ns;
#line 558 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 559 "cgesvdx.f"
		    i__3 = *n;
#line 559 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 560 "cgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 560 "cgesvdx.f"
			i__5 = k;
#line 560 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 560 "cgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 561 "cgesvdx.f"
			++k;
#line 562 "cgesvdx.f"
		    }
#line 563 "cgesvdx.f"
		    k += *n;
#line 564 "cgesvdx.f"
		}
#line 565 "cgesvdx.f"
		i__2 = *m - *n;
#line 565 "cgesvdx.f"
		claset_("A", &i__2, ns, &c_b1, &c_b1, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call CUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 570 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 570 "cgesvdx.f"
		cunmbr_("Q", "L", "N", n, ns, n, &work[iqrf], n, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call CUNMQR to compute Q*(QB*UB). */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 577 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 577 "cgesvdx.f"
		cunmqr_("L", "N", m, ns, n, &a[a_offset], lda, &work[itau], &
			u[u_offset], ldu, &work[itemp], &i__2, info, (ftnlen)
			1, (ftnlen)1);
#line 580 "cgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 584 "cgesvdx.f"
	    if (wantvt) {
#line 585 "cgesvdx.f"
		k = itgkz + *n;
#line 586 "cgesvdx.f"
		i__2 = *ns;
#line 586 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 587 "cgesvdx.f"
		    i__3 = *n;
#line 587 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 588 "cgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 588 "cgesvdx.f"
			i__5 = k;
#line 588 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 588 "cgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 589 "cgesvdx.f"
			++k;
#line 590 "cgesvdx.f"
		    }
#line 591 "cgesvdx.f"
		    k += *n;
#line 592 "cgesvdx.f"
		}

/*              Call CUNMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 597 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 597 "cgesvdx.f"
		cunmbr_("P", "R", "C", ns, n, n, &work[iqrf], n, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 600 "cgesvdx.f"
	    }
#line 601 "cgesvdx.f"
	} else {

/*           Path 2 (M at least N, but not much larger) */
/*           Reduce A to bidiagonal form without QR decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB) */

#line 611 "cgesvdx.f"
	    itauq = 1;
#line 612 "cgesvdx.f"
	    itaup = itauq + *n;
#line 613 "cgesvdx.f"
	    itemp = itaup + *n;
#line 614 "cgesvdx.f"
	    id = 1;
#line 615 "cgesvdx.f"
	    ie = id + *n;
#line 616 "cgesvdx.f"
	    itgkz = ie + *n;
#line 617 "cgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 617 "cgesvdx.f"
	    cgebrd_(m, n, &a[a_offset], lda, &rwork[id], &rwork[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);
#line 620 "cgesvdx.f"
	    itempr = itgkz + *n * ((*n << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*N*N+14*N) */

#line 625 "cgesvdx.f"
	    i__2 = *n << 1;
#line 625 "cgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, n, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 632 "cgesvdx.f"
	    if (wantu) {
#line 633 "cgesvdx.f"
		k = itgkz;
#line 634 "cgesvdx.f"
		i__2 = *ns;
#line 634 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 635 "cgesvdx.f"
		    i__3 = *n;
#line 635 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 636 "cgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 636 "cgesvdx.f"
			i__5 = k;
#line 636 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 636 "cgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 637 "cgesvdx.f"
			++k;
#line 638 "cgesvdx.f"
		    }
#line 639 "cgesvdx.f"
		    k += *n;
#line 640 "cgesvdx.f"
		}
#line 641 "cgesvdx.f"
		i__2 = *m - *n;
#line 641 "cgesvdx.f"
		claset_("A", &i__2, ns, &c_b1, &c_b1, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call CUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 646 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 646 "cgesvdx.f"
		cunmbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 649 "cgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 653 "cgesvdx.f"
	    if (wantvt) {
#line 654 "cgesvdx.f"
		k = itgkz + *n;
#line 655 "cgesvdx.f"
		i__2 = *ns;
#line 655 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 656 "cgesvdx.f"
		    i__3 = *n;
#line 656 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 657 "cgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 657 "cgesvdx.f"
			i__5 = k;
#line 657 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 657 "cgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 658 "cgesvdx.f"
			++k;
#line 659 "cgesvdx.f"
		    }
#line 660 "cgesvdx.f"
		    k += *n;
#line 661 "cgesvdx.f"
		}

/*              Call CUNMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 666 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 666 "cgesvdx.f"
		cunmbr_("P", "R", "C", ns, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 669 "cgesvdx.f"
	    }
#line 670 "cgesvdx.f"
	}
#line 671 "cgesvdx.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce A using the LQ decomposition. */

#line 676 "cgesvdx.f"
	if (*n >= mnthr) {

/*           Path 1t (N much larger than M): */
/*           A = L * Q = ( QB * B * PB**T ) * Q */
/*                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q */
/*           U = QB * UB ; V**T = VB**T * PB**T * Q */

/*           Compute A=L*Q */
/*           (Workspace: need 2*M, prefer M+M*NB) */

#line 686 "cgesvdx.f"
	    itau = 1;
#line 687 "cgesvdx.f"
	    itemp = itau + *m;
#line 688 "cgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 688 "cgesvdx.f"
	    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);
/*           Copy L into WORK and bidiagonalize it: */
/*           (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB) */

#line 694 "cgesvdx.f"
	    ilqf = itemp;
#line 695 "cgesvdx.f"
	    itauq = ilqf + *m * *m;
#line 696 "cgesvdx.f"
	    itaup = itauq + *m;
#line 697 "cgesvdx.f"
	    itemp = itaup + *m;
#line 698 "cgesvdx.f"
	    id = 1;
#line 699 "cgesvdx.f"
	    ie = id + *m;
#line 700 "cgesvdx.f"
	    itgkz = ie + *m;
#line 701 "cgesvdx.f"
	    clacpy_("L", m, m, &a[a_offset], lda, &work[ilqf], m, (ftnlen)1);
#line 702 "cgesvdx.f"
	    i__2 = *m - 1;
#line 702 "cgesvdx.f"
	    i__3 = *m - 1;
#line 702 "cgesvdx.f"
	    claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ilqf + *m], m, (
		    ftnlen)1);
#line 704 "cgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 704 "cgesvdx.f"
	    cgebrd_(m, m, &work[ilqf], m, &rwork[id], &rwork[ie], &work[itauq]
		    , &work[itaup], &work[itemp], &i__2, info);
#line 707 "cgesvdx.f"
	    itempr = itgkz + *m * ((*m << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 712 "cgesvdx.f"
	    i__2 = *m << 1;
#line 712 "cgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, m, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 719 "cgesvdx.f"
	    if (wantu) {
#line 720 "cgesvdx.f"
		k = itgkz;
#line 721 "cgesvdx.f"
		i__2 = *ns;
#line 721 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 722 "cgesvdx.f"
		    i__3 = *m;
#line 722 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 723 "cgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 723 "cgesvdx.f"
			i__5 = k;
#line 723 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 723 "cgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 724 "cgesvdx.f"
			++k;
#line 725 "cgesvdx.f"
		    }
#line 726 "cgesvdx.f"
		    k += *m;
#line 727 "cgesvdx.f"
		}

/*              Call CUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 732 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 732 "cgesvdx.f"
		cunmbr_("Q", "L", "N", m, ns, m, &work[ilqf], m, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 735 "cgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 739 "cgesvdx.f"
	    if (wantvt) {
#line 740 "cgesvdx.f"
		k = itgkz + *m;
#line 741 "cgesvdx.f"
		i__2 = *ns;
#line 741 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 742 "cgesvdx.f"
		    i__3 = *m;
#line 742 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 743 "cgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 743 "cgesvdx.f"
			i__5 = k;
#line 743 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 743 "cgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 744 "cgesvdx.f"
			++k;
#line 745 "cgesvdx.f"
		    }
#line 746 "cgesvdx.f"
		    k += *m;
#line 747 "cgesvdx.f"
		}
#line 748 "cgesvdx.f"
		i__2 = *n - *m;
#line 748 "cgesvdx.f"
		claset_("A", ns, &i__2, &c_b1, &c_b1, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call CUNMBR to compute (VB**T)*(PB**T) */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 754 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 754 "cgesvdx.f"
		cunmbr_("P", "R", "C", ns, m, m, &work[ilqf], m, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call CUNMLQ to compute ((VB**T)*(PB**T))*Q. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 761 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 761 "cgesvdx.f"
		cunmlq_("R", "N", ns, n, m, &a[a_offset], lda, &work[itau], &
			vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1);
#line 764 "cgesvdx.f"
	    }
#line 765 "cgesvdx.f"
	} else {

/*           Path 2t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB) */

#line 775 "cgesvdx.f"
	    itauq = 1;
#line 776 "cgesvdx.f"
	    itaup = itauq + *m;
#line 777 "cgesvdx.f"
	    itemp = itaup + *m;
#line 778 "cgesvdx.f"
	    id = 1;
#line 779 "cgesvdx.f"
	    ie = id + *m;
#line 780 "cgesvdx.f"
	    itgkz = ie + *m;
#line 781 "cgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 781 "cgesvdx.f"
	    cgebrd_(m, n, &a[a_offset], lda, &rwork[id], &rwork[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);
#line 784 "cgesvdx.f"
	    itempr = itgkz + *m * ((*m << 1) + 1);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 789 "cgesvdx.f"
	    i__2 = *m << 1;
#line 789 "cgesvdx.f"
	    sbdsvdx_("L", jobz, rngtgk, m, &rwork[id], &rwork[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[
		    itempr], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1)
		    ;

/*           If needed, compute left singular vectors. */

#line 796 "cgesvdx.f"
	    if (wantu) {
#line 797 "cgesvdx.f"
		k = itgkz;
#line 798 "cgesvdx.f"
		i__2 = *ns;
#line 798 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 799 "cgesvdx.f"
		    i__3 = *m;
#line 799 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 800 "cgesvdx.f"
			i__4 = j + i__ * u_dim1;
#line 800 "cgesvdx.f"
			i__5 = k;
#line 800 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 800 "cgesvdx.f"
			u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 801 "cgesvdx.f"
			++k;
#line 802 "cgesvdx.f"
		    }
#line 803 "cgesvdx.f"
		    k += *m;
#line 804 "cgesvdx.f"
		}

/*              Call CUNMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 809 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 809 "cgesvdx.f"
		cunmbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, info, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 812 "cgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 816 "cgesvdx.f"
	    if (wantvt) {
#line 817 "cgesvdx.f"
		k = itgkz + *m;
#line 818 "cgesvdx.f"
		i__2 = *ns;
#line 818 "cgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 819 "cgesvdx.f"
		    i__3 = *m;
#line 819 "cgesvdx.f"
		    for (j = 1; j <= i__3; ++j) {
#line 820 "cgesvdx.f"
			i__4 = i__ + j * vt_dim1;
#line 820 "cgesvdx.f"
			i__5 = k;
#line 820 "cgesvdx.f"
			z__1.r = rwork[i__5], z__1.i = 0.;
#line 820 "cgesvdx.f"
			vt[i__4].r = z__1.r, vt[i__4].i = z__1.i;
#line 821 "cgesvdx.f"
			++k;
#line 822 "cgesvdx.f"
		    }
#line 823 "cgesvdx.f"
		    k += *m;
#line 824 "cgesvdx.f"
		}
#line 825 "cgesvdx.f"
		i__2 = *n - *m;
#line 825 "cgesvdx.f"
		claset_("A", ns, &i__2, &c_b1, &c_b1, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call CUNMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 831 "cgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 831 "cgesvdx.f"
		cunmbr_("P", "R", "C", ns, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, 
			info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 834 "cgesvdx.f"
	    }
#line 835 "cgesvdx.f"
	}
#line 836 "cgesvdx.f"
    }

/*     Undo scaling if necessary */

#line 840 "cgesvdx.f"
    if (iscl == 1) {
#line 841 "cgesvdx.f"
	if (anrm > bignum) {
#line 841 "cgesvdx.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 841 "cgesvdx.f"
	}
#line 844 "cgesvdx.f"
	if (anrm < smlnum) {
#line 844 "cgesvdx.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 844 "cgesvdx.f"
	}
#line 847 "cgesvdx.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 851 "cgesvdx.f"
    d__1 = (doublereal) maxwrk;
#line 851 "cgesvdx.f"
    z__1.r = d__1, z__1.i = 0.;
#line 851 "cgesvdx.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

#line 853 "cgesvdx.f"
    return 0;

/*     End of CGESVDX */

} /* cgesvdx_ */

