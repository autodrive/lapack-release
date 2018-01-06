#line 1 "dgesvdx.f"
/* dgesvdx.f -- translated by f2c (version 20100827).
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

#line 1 "dgesvdx.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b69 = 0.;

/* > \brief <b> DGESVDX computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGESVDX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvdx
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvdx
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvdx
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE DGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, */
/*    $                    IL, IU, NS, S, U, LDU, VT, LDVT, WORK, */
/*    $                    LWORK, IWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER          JOBU, JOBVT, RANGE */
/*      INTEGER            IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS */
/*      DOUBLE PRECISION   VL, VU */
/*     .. */
/*     .. Array Arguments .. */
/*     INTEGER            IWORK( * ) */
/*     DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ), */
/*    $                   VT( LDVT, * ), WORK( * ) */
/*     .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >  DGESVDX computes the singular value decomposition (SVD) of a real */
/* >  M-by-N matrix A, optionally computing the left and/or right singular */
/* >  vectors. The SVD is written */
/* > */
/* >      A = U * SIGMA * transpose(V) */
/* > */
/* >  where SIGMA is an M-by-N matrix which is zero except for its */
/* >  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and */
/* >  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA */
/* >  are the singular values of A; they are real and non-negative, and */
/* >  are returned in descending order.  The first min(m,n) columns of */
/* >  U and V are the left and right singular vectors of A. */
/* > */
/* >  DGESVDX uses an eigenvalue problem for obtaining the SVD, which */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          VL >=0. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          If RANGE='V', the lower and upper bounds of the interval to */
/* >          be searched for singular values. VU > VL. */
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
/* >          smallest and largest singular values to be returned. */
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
/* >          U is DOUBLE PRECISION array, dimension (LDU,UCOL) */
/* >          If JOBU = 'V', U contains columns of U (the left singular */
/* >          vectors, stored columnwise) as specified by RANGE; if */
/* >          JOBU = 'N', U is not referenced. */
/* >          Note: The user must ensure that UCOL >= NS; if RANGE = 'V', */
/* >          the exact value of NS is not known ILQFin advance and an upper */
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
/* >          VT is DOUBLE PRECISION array, dimension (LDVT,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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

/* > \date November 2015 */

/* > \ingroup doubleGEsing */

/*  ===================================================================== */
/* Subroutine */ int dgesvdx_(char *jobu, char *jobvt, char *range, integer *
	m, integer *n, doublereal *a, integer *lda, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, integer *ns, doublereal *s, 
	doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, 
	doublereal *work, integer *lwork, integer *iwork, integer *info, 
	ftnlen jobu_len, ftnlen jobvt_len, ftnlen range_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1[2], 
	    i__2, i__3;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, id, ie;
    static doublereal dum[1], eps;
    static integer iscl;
    static logical alls, inds;
    static integer ilqf;
    static doublereal anrm;
    static integer ierr, iqrf, itau;
    static char jobz[1];
    static logical vals;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iltgk, itemp, minmn;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer itaup, itauq, iutgk, itgkz, mnthr;
    static logical wantu;
    extern /* Subroutine */ int dgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dlacpy_(char *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum, abstol;
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static char rngtgk[1];
    extern /* Subroutine */ int dormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dormqr_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static logical lquery, wantvt;
    extern /* Subroutine */ int dbdsvdx_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.6.0) -- */
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

#line 311 "dgesvdx.f"
    /* Parameter adjustments */
#line 311 "dgesvdx.f"
    a_dim1 = *lda;
#line 311 "dgesvdx.f"
    a_offset = 1 + a_dim1;
#line 311 "dgesvdx.f"
    a -= a_offset;
#line 311 "dgesvdx.f"
    --s;
#line 311 "dgesvdx.f"
    u_dim1 = *ldu;
#line 311 "dgesvdx.f"
    u_offset = 1 + u_dim1;
#line 311 "dgesvdx.f"
    u -= u_offset;
#line 311 "dgesvdx.f"
    vt_dim1 = *ldvt;
#line 311 "dgesvdx.f"
    vt_offset = 1 + vt_dim1;
#line 311 "dgesvdx.f"
    vt -= vt_offset;
#line 311 "dgesvdx.f"
    --work;
#line 311 "dgesvdx.f"
    --iwork;
#line 311 "dgesvdx.f"

#line 311 "dgesvdx.f"
    /* Function Body */
#line 311 "dgesvdx.f"
    *ns = 0;
#line 312 "dgesvdx.f"
    *info = 0;
#line 313 "dgesvdx.f"
    abstol = dlamch_("S", (ftnlen)1) * 2;
#line 314 "dgesvdx.f"
    lquery = *lwork == -1;
#line 315 "dgesvdx.f"
    minmn = min(*m,*n);
#line 317 "dgesvdx.f"
    wantu = lsame_(jobu, "V", (ftnlen)1, (ftnlen)1);
#line 318 "dgesvdx.f"
    wantvt = lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1);
#line 319 "dgesvdx.f"
    if (wantu || wantvt) {
#line 320 "dgesvdx.f"
	*(unsigned char *)jobz = 'V';
#line 321 "dgesvdx.f"
    } else {
#line 322 "dgesvdx.f"
	*(unsigned char *)jobz = 'N';
#line 323 "dgesvdx.f"
    }
#line 324 "dgesvdx.f"
    alls = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 325 "dgesvdx.f"
    vals = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 326 "dgesvdx.f"
    inds = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 328 "dgesvdx.f"
    *info = 0;
#line 329 "dgesvdx.f"
    if (! lsame_(jobu, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobu, "N", (
	    ftnlen)1, (ftnlen)1)) {
#line 331 "dgesvdx.f"
	*info = -1;
#line 332 "dgesvdx.f"
    } else if (! lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobvt, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 334 "dgesvdx.f"
	*info = -2;
#line 335 "dgesvdx.f"
    } else if (! (alls || vals || inds)) {
#line 336 "dgesvdx.f"
	*info = -3;
#line 337 "dgesvdx.f"
    } else if (*m < 0) {
#line 338 "dgesvdx.f"
	*info = -4;
#line 339 "dgesvdx.f"
    } else if (*n < 0) {
#line 340 "dgesvdx.f"
	*info = -5;
#line 341 "dgesvdx.f"
    } else if (*m > *lda) {
#line 342 "dgesvdx.f"
	*info = -7;
#line 343 "dgesvdx.f"
    } else if (minmn > 0) {
#line 344 "dgesvdx.f"
	if (vals) {
#line 345 "dgesvdx.f"
	    if (*vl < 0.) {
#line 346 "dgesvdx.f"
		*info = -8;
#line 347 "dgesvdx.f"
	    } else if (*vu <= *vl) {
#line 348 "dgesvdx.f"
		*info = -9;
#line 349 "dgesvdx.f"
	    }
#line 350 "dgesvdx.f"
	} else if (inds) {
#line 351 "dgesvdx.f"
	    if (*il < 1 || *il > max(1,minmn)) {
#line 352 "dgesvdx.f"
		*info = -10;
#line 353 "dgesvdx.f"
	    } else if (*iu < min(minmn,*il) || *iu > minmn) {
#line 354 "dgesvdx.f"
		*info = -11;
#line 355 "dgesvdx.f"
	    }
#line 356 "dgesvdx.f"
	}
#line 357 "dgesvdx.f"
	if (*info == 0) {
#line 358 "dgesvdx.f"
	    if (wantu && *ldu < *m) {
#line 359 "dgesvdx.f"
		*info = -15;
#line 360 "dgesvdx.f"
	    } else if (wantvt && *ldvt < minmn) {
#line 361 "dgesvdx.f"
		*info = -16;
#line 362 "dgesvdx.f"
	    }
#line 363 "dgesvdx.f"
	}
#line 364 "dgesvdx.f"
    }

/*     Compute workspace */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 373 "dgesvdx.f"
    if (*info == 0) {
#line 374 "dgesvdx.f"
	minwrk = 1;
#line 375 "dgesvdx.f"
	maxwrk = 1;
#line 376 "dgesvdx.f"
	if (minmn > 0) {
#line 377 "dgesvdx.f"
	    if (*m >= *n) {
/* Writing concatenation */
#line 378 "dgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 378 "dgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 378 "dgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 378 "dgesvdx.f"
		mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 379 "dgesvdx.f"
		if (*m >= mnthr) {

/*                 Path 1 (M much larger than N) */

#line 383 "dgesvdx.f"
		    maxwrk = *n * ((*n << 1) + 16) + *n * ilaenv_(&c__1, 
			    "DGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
/* Computing MAX */
#line 385 "dgesvdx.f"
		    i__2 = maxwrk, i__3 = *n * ((*n << 1) + 20) + (*n << 1) * 
			    ilaenv_(&c__1, "DGEBRD", " ", n, n, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 385 "dgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 387 "dgesvdx.f"
		    minwrk = *n * ((*n << 1) + 21);
#line 388 "dgesvdx.f"
		} else {

/*                 Path 2 (M at least N, but not much larger) */

#line 392 "dgesvdx.f"
		    maxwrk = *n * ((*n << 1) + 19) + (*m + *n) * ilaenv_(&
			    c__1, "DGEBRD", " ", m, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 394 "dgesvdx.f"
		    minwrk = *n * ((*n << 1) + 20) + *m;
#line 395 "dgesvdx.f"
		}
#line 396 "dgesvdx.f"
	    } else {
/* Writing concatenation */
#line 397 "dgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 397 "dgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 397 "dgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 397 "dgesvdx.f"
		mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 398 "dgesvdx.f"
		if (*n >= mnthr) {

/*                 Path 1t (N much larger than M) */

#line 402 "dgesvdx.f"
		    maxwrk = *m * ((*m << 1) + 16) + *m * ilaenv_(&c__1, 
			    "DGELQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
/* Computing MAX */
#line 404 "dgesvdx.f"
		    i__2 = maxwrk, i__3 = *m * ((*m << 1) + 20) + (*m << 1) * 
			    ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 404 "dgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 406 "dgesvdx.f"
		    minwrk = *m * ((*m << 1) + 21);
#line 407 "dgesvdx.f"
		} else {

/*                 Path 2t (N greater than M, but not much larger) */

#line 411 "dgesvdx.f"
		    maxwrk = *m * ((*m << 1) + 19) + (*m + *n) * ilaenv_(&
			    c__1, "DGEBRD", " ", m, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 413 "dgesvdx.f"
		    minwrk = *m * ((*m << 1) + 20) + *n;
#line 414 "dgesvdx.f"
		}
#line 415 "dgesvdx.f"
	    }
#line 416 "dgesvdx.f"
	}
#line 417 "dgesvdx.f"
	maxwrk = max(maxwrk,minwrk);
#line 418 "dgesvdx.f"
	work[1] = (doublereal) maxwrk;

#line 420 "dgesvdx.f"
	if (*lwork < minwrk && ! lquery) {
#line 421 "dgesvdx.f"
	    *info = -19;
#line 422 "dgesvdx.f"
	}
#line 423 "dgesvdx.f"
    }

#line 425 "dgesvdx.f"
    if (*info != 0) {
#line 426 "dgesvdx.f"
	i__2 = -(*info);
#line 426 "dgesvdx.f"
	xerbla_("DGESVDX", &i__2, (ftnlen)7);
#line 427 "dgesvdx.f"
	return 0;
#line 428 "dgesvdx.f"
    } else if (lquery) {
#line 429 "dgesvdx.f"
	return 0;
#line 430 "dgesvdx.f"
    }

/*     Quick return if possible */

#line 434 "dgesvdx.f"
    if (*m == 0 || *n == 0) {
#line 435 "dgesvdx.f"
	return 0;
#line 436 "dgesvdx.f"
    }

/*     Set singular values indices accord to RANGE. */

#line 440 "dgesvdx.f"
    if (alls) {
#line 441 "dgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 442 "dgesvdx.f"
	iltgk = 1;
#line 443 "dgesvdx.f"
	iutgk = min(*m,*n);
#line 444 "dgesvdx.f"
    } else if (inds) {
#line 445 "dgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 446 "dgesvdx.f"
	iltgk = *il;
#line 447 "dgesvdx.f"
	iutgk = *iu;
#line 448 "dgesvdx.f"
    } else {
#line 449 "dgesvdx.f"
	*(unsigned char *)rngtgk = 'V';
#line 450 "dgesvdx.f"
	iltgk = 0;
#line 451 "dgesvdx.f"
	iutgk = 0;
#line 452 "dgesvdx.f"
    }

/*     Get machine constants */

#line 456 "dgesvdx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 457 "dgesvdx.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 458 "dgesvdx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 462 "dgesvdx.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 463 "dgesvdx.f"
    iscl = 0;
#line 464 "dgesvdx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 465 "dgesvdx.f"
	iscl = 1;
#line 466 "dgesvdx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 467 "dgesvdx.f"
    } else if (anrm > bignum) {
#line 468 "dgesvdx.f"
	iscl = 1;
#line 469 "dgesvdx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 470 "dgesvdx.f"
    }

#line 472 "dgesvdx.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce A using the QR */
/*        decomposition. */

#line 478 "dgesvdx.f"
	if (*m >= mnthr) {

/*           Path 1 (M much larger than N): */
/*           A = Q * R = Q * ( QB * B * PB**T ) */
/*                     = Q * ( QB * ( UB * S * VB**T ) * PB**T ) */
/*           U = Q * QB * UB; V**T = VB**T * PB**T */

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 488 "dgesvdx.f"
	    itau = 1;
#line 489 "dgesvdx.f"
	    itemp = itau + *n;
#line 490 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 490 "dgesvdx.f"
	    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);

/*           Copy R into WORK and bidiagonalize it: */
/*           (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB) */

#line 496 "dgesvdx.f"
	    iqrf = itemp;
#line 497 "dgesvdx.f"
	    id = iqrf + *n * *n;
#line 498 "dgesvdx.f"
	    ie = id + *n;
#line 499 "dgesvdx.f"
	    itauq = ie + *n;
#line 500 "dgesvdx.f"
	    itaup = itauq + *n;
#line 501 "dgesvdx.f"
	    itemp = itaup + *n;
#line 502 "dgesvdx.f"
	    dlacpy_("U", n, n, &a[a_offset], lda, &work[iqrf], n, (ftnlen)1);
#line 503 "dgesvdx.f"
	    i__2 = *n - 1;
#line 503 "dgesvdx.f"
	    i__3 = *n - 1;
#line 503 "dgesvdx.f"
	    dlaset_("L", &i__2, &i__3, &c_b69, &c_b69, &work[iqrf + 1], n, (
		    ftnlen)1);
#line 504 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 504 "dgesvdx.f"
	    dgebrd_(n, n, &work[iqrf], n, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 511 "dgesvdx.f"
	    itgkz = itemp;
#line 512 "dgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 513 "dgesvdx.f"
	    i__2 = *n << 1;
#line 513 "dgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 519 "dgesvdx.f"
	    if (wantu) {
#line 520 "dgesvdx.f"
		j = itgkz;
#line 521 "dgesvdx.f"
		i__2 = *ns;
#line 521 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 522 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 523 "dgesvdx.f"
		    j += *n << 1;
#line 524 "dgesvdx.f"
		}
#line 525 "dgesvdx.f"
		i__2 = *m - *n;
#line 525 "dgesvdx.f"
		dlaset_("A", &i__2, n, &c_b69, &c_b69, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 530 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 530 "dgesvdx.f"
		dormbr_("Q", "L", "N", n, ns, n, &work[iqrf], n, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call DORMQR to compute Q*(QB*UB). */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 537 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 537 "dgesvdx.f"
		dormqr_("L", "N", m, ns, n, &a[a_offset], lda, &work[itau], &
			u[u_offset], ldu, &work[itemp], &i__2, info, (ftnlen)
			1, (ftnlen)1);
#line 540 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 544 "dgesvdx.f"
	    if (wantvt) {
#line 545 "dgesvdx.f"
		j = itgkz + *n;
#line 546 "dgesvdx.f"
		i__2 = *ns;
#line 546 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 547 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 548 "dgesvdx.f"
		    j += *n << 1;
#line 549 "dgesvdx.f"
		}

/*              Call DORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 554 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 554 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, n, n, &work[iqrf], n, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 557 "dgesvdx.f"
	    }
#line 558 "dgesvdx.f"
	} else {

/*           Path 2 (M at least N, but not much larger) */
/*           Reduce A to bidiagonal form without QR decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB) */

#line 568 "dgesvdx.f"
	    id = 1;
#line 569 "dgesvdx.f"
	    ie = id + *n;
#line 570 "dgesvdx.f"
	    itauq = ie + *n;
#line 571 "dgesvdx.f"
	    itaup = itauq + *n;
#line 572 "dgesvdx.f"
	    itemp = itaup + *n;
#line 573 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 573 "dgesvdx.f"
	    dgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 580 "dgesvdx.f"
	    itgkz = itemp;
#line 581 "dgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 582 "dgesvdx.f"
	    i__2 = *n << 1;
#line 582 "dgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 588 "dgesvdx.f"
	    if (wantu) {
#line 589 "dgesvdx.f"
		j = itgkz;
#line 590 "dgesvdx.f"
		i__2 = *ns;
#line 590 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 591 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 592 "dgesvdx.f"
		    j += *n << 1;
#line 593 "dgesvdx.f"
		}
#line 594 "dgesvdx.f"
		i__2 = *m - *n;
#line 594 "dgesvdx.f"
		dlaset_("A", &i__2, n, &c_b69, &c_b69, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 599 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 599 "dgesvdx.f"
		dormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 602 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 606 "dgesvdx.f"
	    if (wantvt) {
#line 607 "dgesvdx.f"
		j = itgkz + *n;
#line 608 "dgesvdx.f"
		i__2 = *ns;
#line 608 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 609 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 610 "dgesvdx.f"
		    j += *n << 1;
#line 611 "dgesvdx.f"
		}

/*              Call DORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 616 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 616 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 619 "dgesvdx.f"
	    }
#line 620 "dgesvdx.f"
	}
#line 621 "dgesvdx.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce A using the LQ decomposition. */

#line 626 "dgesvdx.f"
	if (*n >= mnthr) {

/*           Path 1t (N much larger than M): */
/*           A = L * Q = ( QB * B * PB**T ) * Q */
/*                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q */
/*           U = QB * UB ; V**T = VB**T * PB**T * Q */

/*           Compute A=L*Q */
/*           (Workspace: need 2*M, prefer M+M*NB) */

#line 636 "dgesvdx.f"
	    itau = 1;
#line 637 "dgesvdx.f"
	    itemp = itau + *m;
#line 638 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 638 "dgesvdx.f"
	    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);
/*           Copy L into WORK and bidiagonalize it: */
/*           (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB) */

#line 644 "dgesvdx.f"
	    ilqf = itemp;
#line 645 "dgesvdx.f"
	    id = ilqf + *m * *m;
#line 646 "dgesvdx.f"
	    ie = id + *m;
#line 647 "dgesvdx.f"
	    itauq = ie + *m;
#line 648 "dgesvdx.f"
	    itaup = itauq + *m;
#line 649 "dgesvdx.f"
	    itemp = itaup + *m;
#line 650 "dgesvdx.f"
	    dlacpy_("L", m, m, &a[a_offset], lda, &work[ilqf], m, (ftnlen)1);
#line 651 "dgesvdx.f"
	    i__2 = *m - 1;
#line 651 "dgesvdx.f"
	    i__3 = *m - 1;
#line 651 "dgesvdx.f"
	    dlaset_("U", &i__2, &i__3, &c_b69, &c_b69, &work[ilqf + *m], m, (
		    ftnlen)1);
#line 652 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 652 "dgesvdx.f"
	    dgebrd_(m, m, &work[ilqf], m, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 659 "dgesvdx.f"
	    itgkz = itemp;
#line 660 "dgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 661 "dgesvdx.f"
	    i__2 = *m << 1;
#line 661 "dgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 667 "dgesvdx.f"
	    if (wantu) {
#line 668 "dgesvdx.f"
		j = itgkz;
#line 669 "dgesvdx.f"
		i__2 = *ns;
#line 669 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 670 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 671 "dgesvdx.f"
		    j += *m << 1;
#line 672 "dgesvdx.f"
		}

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 677 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 677 "dgesvdx.f"
		dormbr_("Q", "L", "N", m, ns, m, &work[ilqf], m, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 680 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 684 "dgesvdx.f"
	    if (wantvt) {
#line 685 "dgesvdx.f"
		j = itgkz + *m;
#line 686 "dgesvdx.f"
		i__2 = *ns;
#line 686 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 687 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 688 "dgesvdx.f"
		    j += *m << 1;
#line 689 "dgesvdx.f"
		}
#line 690 "dgesvdx.f"
		i__2 = *n - *m;
#line 690 "dgesvdx.f"
		dlaset_("A", m, &i__2, &c_b69, &c_b69, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call DORMBR to compute (VB**T)*(PB**T) */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 695 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 695 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, m, m, &work[ilqf], m, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call DORMLQ to compute ((VB**T)*(PB**T))*Q. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 702 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 702 "dgesvdx.f"
		dormlq_("R", "N", ns, n, m, &a[a_offset], lda, &work[itau], &
			vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1);
#line 705 "dgesvdx.f"
	    }
#line 706 "dgesvdx.f"
	} else {

/*           Path 2t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB) */

#line 716 "dgesvdx.f"
	    id = 1;
#line 717 "dgesvdx.f"
	    ie = id + *m;
#line 718 "dgesvdx.f"
	    itauq = ie + *m;
#line 719 "dgesvdx.f"
	    itaup = itauq + *m;
#line 720 "dgesvdx.f"
	    itemp = itaup + *m;
#line 721 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 721 "dgesvdx.f"
	    dgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 728 "dgesvdx.f"
	    itgkz = itemp;
#line 729 "dgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 730 "dgesvdx.f"
	    i__2 = *m << 1;
#line 730 "dgesvdx.f"
	    dbdsvdx_("L", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 736 "dgesvdx.f"
	    if (wantu) {
#line 737 "dgesvdx.f"
		j = itgkz;
#line 738 "dgesvdx.f"
		i__2 = *ns;
#line 738 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 739 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 740 "dgesvdx.f"
		    j += *m << 1;
#line 741 "dgesvdx.f"
		}

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 746 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 746 "dgesvdx.f"
		dormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, info, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 749 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 753 "dgesvdx.f"
	    if (wantvt) {
#line 754 "dgesvdx.f"
		j = itgkz + *m;
#line 755 "dgesvdx.f"
		i__2 = *ns;
#line 755 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 756 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 757 "dgesvdx.f"
		    j += *m << 1;
#line 758 "dgesvdx.f"
		}
#line 759 "dgesvdx.f"
		i__2 = *n - *m;
#line 759 "dgesvdx.f"
		dlaset_("A", m, &i__2, &c_b69, &c_b69, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call DORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 764 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 764 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, 
			info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 767 "dgesvdx.f"
	    }
#line 768 "dgesvdx.f"
	}
#line 769 "dgesvdx.f"
    }

/*     Undo scaling if necessary */

#line 773 "dgesvdx.f"
    if (iscl == 1) {
#line 774 "dgesvdx.f"
	if (anrm > bignum) {
#line 774 "dgesvdx.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 774 "dgesvdx.f"
	}
#line 777 "dgesvdx.f"
	if (anrm < smlnum) {
#line 777 "dgesvdx.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 777 "dgesvdx.f"
	}
#line 780 "dgesvdx.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 784 "dgesvdx.f"
    work[1] = (doublereal) maxwrk;

#line 786 "dgesvdx.f"
    return 0;

/*     End of DGESVDX */

} /* dgesvdx_ */

