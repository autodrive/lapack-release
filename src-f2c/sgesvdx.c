#line 1 "sgesvdx.f"
/* sgesvdx.f -- translated by f2c (version 20100827).
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

#line 1 "sgesvdx.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b69 = 0.;

/* > \brief <b> SGESVDX computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGESVDX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvdx
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvdx
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvdx
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE SGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, */
/*    $                    IL, IU, NS, S, U, LDU, VT, LDVT, WORK, */
/*    $                    LWORK, IWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER          JOBU, JOBVT, RANGE */
/*      INTEGER            IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS */
/*      REAL               VL, VU */
/*     .. */
/*     .. Array Arguments .. */
/*     INTEGER            IWORK( * ) */
/*     REAL               A( LDA, * ), S( * ), U( LDU, * ), */
/*    $                   VT( LDVT, * ), WORK( * ) */
/*     .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >  SGESVDX computes the singular value decomposition (SVD) of a real */
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
/* >  SGESVDX uses an eigenvalue problem for obtaining the SVD, which */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          VL >=0. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
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
/* >          S is REAL array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array, dimension (LDU,UCOL) */
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
/* >          VT is REAL array, dimension (LDVT,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \date November 2015 */

/* > \ingroup realGEsing */

/*  ===================================================================== */
/* Subroutine */ int sgesvdx_(char *jobu, char *jobvt, char *range, integer *
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
    static integer iltgk, itemp, minmn, itaup, itauq, iutgk, itgkz, mnthr;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantu;
    extern /* Subroutine */ int sgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int sgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    slascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static doublereal abstol;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    slacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static char rngtgk[1];
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    sormbr_(char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int sormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical lquery, wantvt;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    sbdsvdx_(char *, char *, char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);


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

#line 311 "sgesvdx.f"
    /* Parameter adjustments */
#line 311 "sgesvdx.f"
    a_dim1 = *lda;
#line 311 "sgesvdx.f"
    a_offset = 1 + a_dim1;
#line 311 "sgesvdx.f"
    a -= a_offset;
#line 311 "sgesvdx.f"
    --s;
#line 311 "sgesvdx.f"
    u_dim1 = *ldu;
#line 311 "sgesvdx.f"
    u_offset = 1 + u_dim1;
#line 311 "sgesvdx.f"
    u -= u_offset;
#line 311 "sgesvdx.f"
    vt_dim1 = *ldvt;
#line 311 "sgesvdx.f"
    vt_offset = 1 + vt_dim1;
#line 311 "sgesvdx.f"
    vt -= vt_offset;
#line 311 "sgesvdx.f"
    --work;
#line 311 "sgesvdx.f"
    --iwork;
#line 311 "sgesvdx.f"

#line 311 "sgesvdx.f"
    /* Function Body */
#line 311 "sgesvdx.f"
    *ns = 0;
#line 312 "sgesvdx.f"
    *info = 0;
#line 313 "sgesvdx.f"
    abstol = slamch_("S", (ftnlen)1) * 2;
#line 314 "sgesvdx.f"
    lquery = *lwork == -1;
#line 315 "sgesvdx.f"
    minmn = min(*m,*n);
#line 317 "sgesvdx.f"
    wantu = lsame_(jobu, "V", (ftnlen)1, (ftnlen)1);
#line 318 "sgesvdx.f"
    wantvt = lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1);
#line 319 "sgesvdx.f"
    if (wantu || wantvt) {
#line 320 "sgesvdx.f"
	*(unsigned char *)jobz = 'V';
#line 321 "sgesvdx.f"
    } else {
#line 322 "sgesvdx.f"
	*(unsigned char *)jobz = 'N';
#line 323 "sgesvdx.f"
    }
#line 324 "sgesvdx.f"
    alls = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 325 "sgesvdx.f"
    vals = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 326 "sgesvdx.f"
    inds = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 328 "sgesvdx.f"
    *info = 0;
#line 329 "sgesvdx.f"
    if (! lsame_(jobu, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobu, "N", (
	    ftnlen)1, (ftnlen)1)) {
#line 331 "sgesvdx.f"
	*info = -1;
#line 332 "sgesvdx.f"
    } else if (! lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobvt, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 334 "sgesvdx.f"
	*info = -2;
#line 335 "sgesvdx.f"
    } else if (! (alls || vals || inds)) {
#line 336 "sgesvdx.f"
	*info = -3;
#line 337 "sgesvdx.f"
    } else if (*m < 0) {
#line 338 "sgesvdx.f"
	*info = -4;
#line 339 "sgesvdx.f"
    } else if (*n < 0) {
#line 340 "sgesvdx.f"
	*info = -5;
#line 341 "sgesvdx.f"
    } else if (*m > *lda) {
#line 342 "sgesvdx.f"
	*info = -7;
#line 343 "sgesvdx.f"
    } else if (minmn > 0) {
#line 344 "sgesvdx.f"
	if (vals) {
#line 345 "sgesvdx.f"
	    if (*vl < 0.) {
#line 346 "sgesvdx.f"
		*info = -8;
#line 347 "sgesvdx.f"
	    } else if (*vu <= *vl) {
#line 348 "sgesvdx.f"
		*info = -9;
#line 349 "sgesvdx.f"
	    }
#line 350 "sgesvdx.f"
	} else if (inds) {
#line 351 "sgesvdx.f"
	    if (*il < 1 || *il > max(1,minmn)) {
#line 352 "sgesvdx.f"
		*info = -10;
#line 353 "sgesvdx.f"
	    } else if (*iu < min(minmn,*il) || *iu > minmn) {
#line 354 "sgesvdx.f"
		*info = -11;
#line 355 "sgesvdx.f"
	    }
#line 356 "sgesvdx.f"
	}
#line 357 "sgesvdx.f"
	if (*info == 0) {
#line 358 "sgesvdx.f"
	    if (wantu && *ldu < *m) {
#line 359 "sgesvdx.f"
		*info = -15;
#line 360 "sgesvdx.f"
	    } else if (wantvt && *ldvt < minmn) {
#line 361 "sgesvdx.f"
		*info = -16;
#line 362 "sgesvdx.f"
	    }
#line 363 "sgesvdx.f"
	}
#line 364 "sgesvdx.f"
    }

/*     Compute workspace */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 373 "sgesvdx.f"
    if (*info == 0) {
#line 374 "sgesvdx.f"
	minwrk = 1;
#line 375 "sgesvdx.f"
	maxwrk = 1;
#line 376 "sgesvdx.f"
	if (minmn > 0) {
#line 377 "sgesvdx.f"
	    if (*m >= *n) {
/* Writing concatenation */
#line 378 "sgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 378 "sgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 378 "sgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 378 "sgesvdx.f"
		mnthr = ilaenv_(&c__6, "SGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 379 "sgesvdx.f"
		if (*m >= mnthr) {

/*                 Path 1 (M much larger than N) */

#line 383 "sgesvdx.f"
		    maxwrk = *n * ((*n << 1) + 16) + *n * ilaenv_(&c__1, 
			    "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
/* Computing MAX */
#line 385 "sgesvdx.f"
		    i__2 = maxwrk, i__3 = *n * ((*n << 1) + 20) + (*n << 1) * 
			    ilaenv_(&c__1, "SGEBRD", " ", n, n, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 385 "sgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 387 "sgesvdx.f"
		    minwrk = *n * ((*n << 1) + 21);
#line 388 "sgesvdx.f"
		} else {

/*                 Path 2 (M at least N, but not much larger) */

#line 392 "sgesvdx.f"
		    maxwrk = *n * ((*n << 1) + 19) + (*m + *n) * ilaenv_(&
			    c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 394 "sgesvdx.f"
		    minwrk = *n * ((*n << 1) + 20) + *m;
#line 395 "sgesvdx.f"
		}
#line 396 "sgesvdx.f"
	    } else {
/* Writing concatenation */
#line 397 "sgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 397 "sgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 397 "sgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 397 "sgesvdx.f"
		mnthr = ilaenv_(&c__6, "SGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 398 "sgesvdx.f"
		if (*n >= mnthr) {

/*                 Path 1t (N much larger than M) */

#line 402 "sgesvdx.f"
		    maxwrk = *m * ((*m << 1) + 16) + *m * ilaenv_(&c__1, 
			    "SGELQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (
			    ftnlen)1);
/* Computing MAX */
#line 404 "sgesvdx.f"
		    i__2 = maxwrk, i__3 = *m * ((*m << 1) + 20) + (*m << 1) * 
			    ilaenv_(&c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 404 "sgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 406 "sgesvdx.f"
		    minwrk = *m * ((*m << 1) + 21);
#line 407 "sgesvdx.f"
		} else {

/*                 Path 2t (N greater than M, but not much larger) */

#line 411 "sgesvdx.f"
		    maxwrk = *m * ((*m << 1) + 19) + (*m + *n) * ilaenv_(&
			    c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 413 "sgesvdx.f"
		    minwrk = *m * ((*m << 1) + 20) + *n;
#line 414 "sgesvdx.f"
		}
#line 415 "sgesvdx.f"
	    }
#line 416 "sgesvdx.f"
	}
#line 417 "sgesvdx.f"
	maxwrk = max(maxwrk,minwrk);
#line 418 "sgesvdx.f"
	work[1] = (doublereal) maxwrk;

#line 420 "sgesvdx.f"
	if (*lwork < minwrk && ! lquery) {
#line 421 "sgesvdx.f"
	    *info = -19;
#line 422 "sgesvdx.f"
	}
#line 423 "sgesvdx.f"
    }

#line 425 "sgesvdx.f"
    if (*info != 0) {
#line 426 "sgesvdx.f"
	i__2 = -(*info);
#line 426 "sgesvdx.f"
	xerbla_("SGESVDX", &i__2, (ftnlen)7);
#line 427 "sgesvdx.f"
	return 0;
#line 428 "sgesvdx.f"
    } else if (lquery) {
#line 429 "sgesvdx.f"
	return 0;
#line 430 "sgesvdx.f"
    }

/*     Quick return if possible */

#line 434 "sgesvdx.f"
    if (*m == 0 || *n == 0) {
#line 435 "sgesvdx.f"
	return 0;
#line 436 "sgesvdx.f"
    }

/*     Set singular values indices accord to RANGE. */

#line 440 "sgesvdx.f"
    if (alls) {
#line 441 "sgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 442 "sgesvdx.f"
	iltgk = 1;
#line 443 "sgesvdx.f"
	iutgk = min(*m,*n);
#line 444 "sgesvdx.f"
    } else if (inds) {
#line 445 "sgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 446 "sgesvdx.f"
	iltgk = *il;
#line 447 "sgesvdx.f"
	iutgk = *iu;
#line 448 "sgesvdx.f"
    } else {
#line 449 "sgesvdx.f"
	*(unsigned char *)rngtgk = 'V';
#line 450 "sgesvdx.f"
	iltgk = 0;
#line 451 "sgesvdx.f"
	iutgk = 0;
#line 452 "sgesvdx.f"
    }

/*     Get machine constants */

#line 456 "sgesvdx.f"
    eps = slamch_("P", (ftnlen)1);
#line 457 "sgesvdx.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 458 "sgesvdx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 462 "sgesvdx.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 463 "sgesvdx.f"
    iscl = 0;
#line 464 "sgesvdx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 465 "sgesvdx.f"
	iscl = 1;
#line 466 "sgesvdx.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 467 "sgesvdx.f"
    } else if (anrm > bignum) {
#line 468 "sgesvdx.f"
	iscl = 1;
#line 469 "sgesvdx.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 470 "sgesvdx.f"
    }

#line 472 "sgesvdx.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce A using the QR */
/*        decomposition. */

#line 478 "sgesvdx.f"
	if (*m >= mnthr) {

/*           Path 1 (M much larger than N): */
/*           A = Q * R = Q * ( QB * B * PB**T ) */
/*                     = Q * ( QB * ( UB * S * VB**T ) * PB**T ) */
/*           U = Q * QB * UB; V**T = VB**T * PB**T */

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 488 "sgesvdx.f"
	    itau = 1;
#line 489 "sgesvdx.f"
	    itemp = itau + *n;
#line 490 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 490 "sgesvdx.f"
	    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);

/*           Copy R into WORK and bidiagonalize it: */
/*           (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB) */

#line 496 "sgesvdx.f"
	    iqrf = itemp;
#line 497 "sgesvdx.f"
	    id = iqrf + *n * *n;
#line 498 "sgesvdx.f"
	    ie = id + *n;
#line 499 "sgesvdx.f"
	    itauq = ie + *n;
#line 500 "sgesvdx.f"
	    itaup = itauq + *n;
#line 501 "sgesvdx.f"
	    itemp = itaup + *n;
#line 502 "sgesvdx.f"
	    slacpy_("U", n, n, &a[a_offset], lda, &work[iqrf], n, (ftnlen)1);
#line 503 "sgesvdx.f"
	    i__2 = *n - 1;
#line 503 "sgesvdx.f"
	    i__3 = *n - 1;
#line 503 "sgesvdx.f"
	    slaset_("L", &i__2, &i__3, &c_b69, &c_b69, &work[iqrf + 1], n, (
		    ftnlen)1);
#line 504 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 504 "sgesvdx.f"
	    sgebrd_(n, n, &work[iqrf], n, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 511 "sgesvdx.f"
	    itgkz = itemp;
#line 512 "sgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 513 "sgesvdx.f"
	    i__2 = *n << 1;
#line 513 "sgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 519 "sgesvdx.f"
	    if (wantu) {
#line 520 "sgesvdx.f"
		j = itgkz;
#line 521 "sgesvdx.f"
		i__2 = *ns;
#line 521 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 522 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 523 "sgesvdx.f"
		    j += *n << 1;
#line 524 "sgesvdx.f"
		}
#line 525 "sgesvdx.f"
		i__2 = *m - *n;
#line 525 "sgesvdx.f"
		slaset_("A", &i__2, n, &c_b69, &c_b69, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 530 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 530 "sgesvdx.f"
		sormbr_("Q", "L", "N", n, ns, n, &work[iqrf], n, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call SORMQR to compute Q*(QB*UB). */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 537 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 537 "sgesvdx.f"
		sormqr_("L", "N", m, ns, n, &a[a_offset], lda, &work[itau], &
			u[u_offset], ldu, &work[itemp], &i__2, info, (ftnlen)
			1, (ftnlen)1);
#line 540 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 544 "sgesvdx.f"
	    if (wantvt) {
#line 545 "sgesvdx.f"
		j = itgkz + *n;
#line 546 "sgesvdx.f"
		i__2 = *ns;
#line 546 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 547 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 548 "sgesvdx.f"
		    j += *n << 1;
#line 549 "sgesvdx.f"
		}

/*              Call SORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 554 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 554 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, n, n, &work[iqrf], n, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 557 "sgesvdx.f"
	    }
#line 558 "sgesvdx.f"
	} else {

/*           Path 2 (M at least N, but not much larger) */
/*           Reduce A to bidiagonal form without QR decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB) */

#line 568 "sgesvdx.f"
	    id = 1;
#line 569 "sgesvdx.f"
	    ie = id + *n;
#line 570 "sgesvdx.f"
	    itauq = ie + *n;
#line 571 "sgesvdx.f"
	    itaup = itauq + *n;
#line 572 "sgesvdx.f"
	    itemp = itaup + *n;
#line 573 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 573 "sgesvdx.f"
	    sgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 580 "sgesvdx.f"
	    itgkz = itemp;
#line 581 "sgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 582 "sgesvdx.f"
	    i__2 = *n << 1;
#line 582 "sgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 588 "sgesvdx.f"
	    if (wantu) {
#line 589 "sgesvdx.f"
		j = itgkz;
#line 590 "sgesvdx.f"
		i__2 = *ns;
#line 590 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 591 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 592 "sgesvdx.f"
		    j += *n << 1;
#line 593 "sgesvdx.f"
		}
#line 594 "sgesvdx.f"
		i__2 = *m - *n;
#line 594 "sgesvdx.f"
		slaset_("A", &i__2, n, &c_b69, &c_b69, &u[*n + 1 + u_dim1], 
			ldu, (ftnlen)1);

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 599 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 599 "sgesvdx.f"
		sormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 602 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 606 "sgesvdx.f"
	    if (wantvt) {
#line 607 "sgesvdx.f"
		j = itgkz + *n;
#line 608 "sgesvdx.f"
		i__2 = *ns;
#line 608 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 609 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 610 "sgesvdx.f"
		    j += *n << 1;
#line 611 "sgesvdx.f"
		}

/*              Call SORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 616 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 616 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 619 "sgesvdx.f"
	    }
#line 620 "sgesvdx.f"
	}
#line 621 "sgesvdx.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce A using the LQ decomposition. */

#line 626 "sgesvdx.f"
	if (*n >= mnthr) {

/*           Path 1t (N much larger than M): */
/*           A = L * Q = ( QB * B * PB**T ) * Q */
/*                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q */
/*           U = QB * UB ; V**T = VB**T * PB**T * Q */

/*           Compute A=L*Q */
/*           (Workspace: need 2*M, prefer M+M*NB) */

#line 636 "sgesvdx.f"
	    itau = 1;
#line 637 "sgesvdx.f"
	    itemp = itau + *m;
#line 638 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 638 "sgesvdx.f"
	    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);
/*           Copy L into WORK and bidiagonalize it: */
/*           (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB) */

#line 644 "sgesvdx.f"
	    ilqf = itemp;
#line 645 "sgesvdx.f"
	    id = ilqf + *m * *m;
#line 646 "sgesvdx.f"
	    ie = id + *m;
#line 647 "sgesvdx.f"
	    itauq = ie + *m;
#line 648 "sgesvdx.f"
	    itaup = itauq + *m;
#line 649 "sgesvdx.f"
	    itemp = itaup + *m;
#line 650 "sgesvdx.f"
	    slacpy_("L", m, m, &a[a_offset], lda, &work[ilqf], m, (ftnlen)1);
#line 651 "sgesvdx.f"
	    i__2 = *m - 1;
#line 651 "sgesvdx.f"
	    i__3 = *m - 1;
#line 651 "sgesvdx.f"
	    slaset_("U", &i__2, &i__3, &c_b69, &c_b69, &work[ilqf + *m], m, (
		    ftnlen)1);
#line 652 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 652 "sgesvdx.f"
	    sgebrd_(m, m, &work[ilqf], m, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 659 "sgesvdx.f"
	    itgkz = itemp;
#line 660 "sgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 661 "sgesvdx.f"
	    i__2 = *m << 1;
#line 661 "sgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 667 "sgesvdx.f"
	    if (wantu) {
#line 668 "sgesvdx.f"
		j = itgkz;
#line 669 "sgesvdx.f"
		i__2 = *ns;
#line 669 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 670 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 671 "sgesvdx.f"
		    j += *m << 1;
#line 672 "sgesvdx.f"
		}

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 677 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 677 "sgesvdx.f"
		sormbr_("Q", "L", "N", m, ns, m, &work[ilqf], m, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 680 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 684 "sgesvdx.f"
	    if (wantvt) {
#line 685 "sgesvdx.f"
		j = itgkz + *m;
#line 686 "sgesvdx.f"
		i__2 = *ns;
#line 686 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 687 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 688 "sgesvdx.f"
		    j += *m << 1;
#line 689 "sgesvdx.f"
		}
#line 690 "sgesvdx.f"
		i__2 = *n - *m;
#line 690 "sgesvdx.f"
		slaset_("A", m, &i__2, &c_b69, &c_b69, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call SORMBR to compute (VB**T)*(PB**T) */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 695 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 695 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, m, m, &work[ilqf], m, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call SORMLQ to compute ((VB**T)*(PB**T))*Q. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 702 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 702 "sgesvdx.f"
		sormlq_("R", "N", ns, n, m, &a[a_offset], lda, &work[itau], &
			vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1);
#line 705 "sgesvdx.f"
	    }
#line 706 "sgesvdx.f"
	} else {

/*           Path 2t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB) */

#line 716 "sgesvdx.f"
	    id = 1;
#line 717 "sgesvdx.f"
	    ie = id + *m;
#line 718 "sgesvdx.f"
	    itauq = ie + *m;
#line 719 "sgesvdx.f"
	    itaup = itauq + *m;
#line 720 "sgesvdx.f"
	    itemp = itaup + *m;
#line 721 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 721 "sgesvdx.f"
	    sgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 728 "sgesvdx.f"
	    itgkz = itemp;
#line 729 "sgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 730 "sgesvdx.f"
	    i__2 = *m << 1;
#line 730 "sgesvdx.f"
	    sbdsvdx_("L", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 736 "sgesvdx.f"
	    if (wantu) {
#line 737 "sgesvdx.f"
		j = itgkz;
#line 738 "sgesvdx.f"
		i__2 = *ns;
#line 738 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 739 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 740 "sgesvdx.f"
		    j += *m << 1;
#line 741 "sgesvdx.f"
		}

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 746 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 746 "sgesvdx.f"
		sormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, info, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 749 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 753 "sgesvdx.f"
	    if (wantvt) {
#line 754 "sgesvdx.f"
		j = itgkz + *m;
#line 755 "sgesvdx.f"
		i__2 = *ns;
#line 755 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 756 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 757 "sgesvdx.f"
		    j += *m << 1;
#line 758 "sgesvdx.f"
		}
#line 759 "sgesvdx.f"
		i__2 = *n - *m;
#line 759 "sgesvdx.f"
		slaset_("A", m, &i__2, &c_b69, &c_b69, &vt[(*m + 1) * vt_dim1 
			+ 1], ldvt, (ftnlen)1);

/*              Call SORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 764 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 764 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, 
			info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 767 "sgesvdx.f"
	    }
#line 768 "sgesvdx.f"
	}
#line 769 "sgesvdx.f"
    }

/*     Undo scaling if necessary */

#line 773 "sgesvdx.f"
    if (iscl == 1) {
#line 774 "sgesvdx.f"
	if (anrm > bignum) {
#line 774 "sgesvdx.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 774 "sgesvdx.f"
	}
#line 777 "sgesvdx.f"
	if (anrm < smlnum) {
#line 777 "sgesvdx.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 777 "sgesvdx.f"
	}
#line 780 "sgesvdx.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 784 "sgesvdx.f"
    work[1] = (doublereal) maxwrk;

#line 786 "sgesvdx.f"
    return 0;

/*     End of SGESVDX */

} /* sgesvdx_ */

