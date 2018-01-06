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
static doublereal c_b109 = 0.;

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
/* >          U is REAL array, dimension (LDU,UCOL) */
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

/* > \date June 2016 */

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

#line 317 "sgesvdx.f"
    /* Parameter adjustments */
#line 317 "sgesvdx.f"
    a_dim1 = *lda;
#line 317 "sgesvdx.f"
    a_offset = 1 + a_dim1;
#line 317 "sgesvdx.f"
    a -= a_offset;
#line 317 "sgesvdx.f"
    --s;
#line 317 "sgesvdx.f"
    u_dim1 = *ldu;
#line 317 "sgesvdx.f"
    u_offset = 1 + u_dim1;
#line 317 "sgesvdx.f"
    u -= u_offset;
#line 317 "sgesvdx.f"
    vt_dim1 = *ldvt;
#line 317 "sgesvdx.f"
    vt_offset = 1 + vt_dim1;
#line 317 "sgesvdx.f"
    vt -= vt_offset;
#line 317 "sgesvdx.f"
    --work;
#line 317 "sgesvdx.f"
    --iwork;
#line 317 "sgesvdx.f"

#line 317 "sgesvdx.f"
    /* Function Body */
#line 317 "sgesvdx.f"
    *ns = 0;
#line 318 "sgesvdx.f"
    *info = 0;
#line 319 "sgesvdx.f"
    abstol = slamch_("S", (ftnlen)1) * 2;
#line 320 "sgesvdx.f"
    lquery = *lwork == -1;
#line 321 "sgesvdx.f"
    minmn = min(*m,*n);
#line 323 "sgesvdx.f"
    wantu = lsame_(jobu, "V", (ftnlen)1, (ftnlen)1);
#line 324 "sgesvdx.f"
    wantvt = lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1);
#line 325 "sgesvdx.f"
    if (wantu || wantvt) {
#line 326 "sgesvdx.f"
	*(unsigned char *)jobz = 'V';
#line 327 "sgesvdx.f"
    } else {
#line 328 "sgesvdx.f"
	*(unsigned char *)jobz = 'N';
#line 329 "sgesvdx.f"
    }
#line 330 "sgesvdx.f"
    alls = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 331 "sgesvdx.f"
    vals = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 332 "sgesvdx.f"
    inds = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 334 "sgesvdx.f"
    *info = 0;
#line 335 "sgesvdx.f"
    if (! lsame_(jobu, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobu, "N", (
	    ftnlen)1, (ftnlen)1)) {
#line 337 "sgesvdx.f"
	*info = -1;
#line 338 "sgesvdx.f"
    } else if (! lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobvt, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 340 "sgesvdx.f"
	*info = -2;
#line 341 "sgesvdx.f"
    } else if (! (alls || vals || inds)) {
#line 342 "sgesvdx.f"
	*info = -3;
#line 343 "sgesvdx.f"
    } else if (*m < 0) {
#line 344 "sgesvdx.f"
	*info = -4;
#line 345 "sgesvdx.f"
    } else if (*n < 0) {
#line 346 "sgesvdx.f"
	*info = -5;
#line 347 "sgesvdx.f"
    } else if (*m > *lda) {
#line 348 "sgesvdx.f"
	*info = -7;
#line 349 "sgesvdx.f"
    } else if (minmn > 0) {
#line 350 "sgesvdx.f"
	if (vals) {
#line 351 "sgesvdx.f"
	    if (*vl < 0.) {
#line 352 "sgesvdx.f"
		*info = -8;
#line 353 "sgesvdx.f"
	    } else if (*vu <= *vl) {
#line 354 "sgesvdx.f"
		*info = -9;
#line 355 "sgesvdx.f"
	    }
#line 356 "sgesvdx.f"
	} else if (inds) {
#line 357 "sgesvdx.f"
	    if (*il < 1 || *il > max(1,minmn)) {
#line 358 "sgesvdx.f"
		*info = -10;
#line 359 "sgesvdx.f"
	    } else if (*iu < min(minmn,*il) || *iu > minmn) {
#line 360 "sgesvdx.f"
		*info = -11;
#line 361 "sgesvdx.f"
	    }
#line 362 "sgesvdx.f"
	}
#line 363 "sgesvdx.f"
	if (*info == 0) {
#line 364 "sgesvdx.f"
	    if (wantu && *ldu < *m) {
#line 365 "sgesvdx.f"
		*info = -15;
#line 366 "sgesvdx.f"
	    } else if (wantvt) {
#line 367 "sgesvdx.f"
		if (inds) {
#line 368 "sgesvdx.f"
		    if (*ldvt < *iu - *il + 1) {
#line 369 "sgesvdx.f"
			*info = -17;
#line 370 "sgesvdx.f"
		    }
#line 371 "sgesvdx.f"
		} else if (*ldvt < minmn) {
#line 372 "sgesvdx.f"
		    *info = -17;
#line 373 "sgesvdx.f"
		}
#line 374 "sgesvdx.f"
	    }
#line 375 "sgesvdx.f"
	}
#line 376 "sgesvdx.f"
    }

/*     Compute workspace */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 385 "sgesvdx.f"
    if (*info == 0) {
#line 386 "sgesvdx.f"
	minwrk = 1;
#line 387 "sgesvdx.f"
	maxwrk = 1;
#line 388 "sgesvdx.f"
	if (minmn > 0) {
#line 389 "sgesvdx.f"
	    if (*m >= *n) {
/* Writing concatenation */
#line 390 "sgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 390 "sgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 390 "sgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 390 "sgesvdx.f"
		mnthr = ilaenv_(&c__6, "SGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 391 "sgesvdx.f"
		if (*m >= mnthr) {

/*                 Path 1 (M much larger than N) */

#line 395 "sgesvdx.f"
		    maxwrk = *n + *n * ilaenv_(&c__1, "SGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 397 "sgesvdx.f"
		    i__2 = maxwrk, i__3 = *n * (*n + 5) + (*n << 1) * ilaenv_(
			    &c__1, "SGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 397 "sgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 399 "sgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 400 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *n * (*n * 3 + 6) + *n * 
				ilaenv_(&c__1, "SORMQR", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 400 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 402 "sgesvdx.f"
		    }
#line 403 "sgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 404 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *n * (*n * 3 + 6) + *n * 
				ilaenv_(&c__1, "SORMLQ", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 404 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 406 "sgesvdx.f"
		    }
#line 407 "sgesvdx.f"
		    minwrk = *n * (*n * 3 + 20);
#line 408 "sgesvdx.f"
		} else {

/*                 Path 2 (M at least N, but not much larger) */

#line 412 "sgesvdx.f"
		    maxwrk = (*n << 2) + (*m + *n) * ilaenv_(&c__1, "SGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 414 "sgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 415 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *n * ((*n << 1) + 5) + *n * 
				ilaenv_(&c__1, "SORMQR", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 415 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 417 "sgesvdx.f"
		    }
#line 418 "sgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 419 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *n * ((*n << 1) + 5) + *n * 
				ilaenv_(&c__1, "SORMLQ", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 419 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 421 "sgesvdx.f"
		    }
/* Computing MAX */
#line 422 "sgesvdx.f"
		    i__2 = *n * ((*n << 1) + 19), i__3 = (*n << 2) + *m;
#line 422 "sgesvdx.f"
		    minwrk = max(i__2,i__3);
#line 423 "sgesvdx.f"
		}
#line 424 "sgesvdx.f"
	    } else {
/* Writing concatenation */
#line 425 "sgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 425 "sgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 425 "sgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 425 "sgesvdx.f"
		mnthr = ilaenv_(&c__6, "SGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 426 "sgesvdx.f"
		if (*n >= mnthr) {

/*                 Path 1t (N much larger than M) */

#line 430 "sgesvdx.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 432 "sgesvdx.f"
		    i__2 = maxwrk, i__3 = *m * (*m + 5) + (*m << 1) * ilaenv_(
			    &c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 432 "sgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 434 "sgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 435 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *m * (*m * 3 + 6) + *m * 
				ilaenv_(&c__1, "SORMQR", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 435 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 437 "sgesvdx.f"
		    }
#line 438 "sgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 439 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *m * (*m * 3 + 6) + *m * 
				ilaenv_(&c__1, "SORMLQ", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 439 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 441 "sgesvdx.f"
		    }
#line 442 "sgesvdx.f"
		    minwrk = *m * (*m * 3 + 20);
#line 443 "sgesvdx.f"
		} else {

/*                 Path 2t (N at least M, but not much larger) */

#line 447 "sgesvdx.f"
		    maxwrk = (*m << 2) + (*m + *n) * ilaenv_(&c__1, "SGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 449 "sgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 450 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *m * ((*m << 1) + 5) + *m * 
				ilaenv_(&c__1, "SORMQR", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 450 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 452 "sgesvdx.f"
		    }
#line 453 "sgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 454 "sgesvdx.f"
			i__2 = maxwrk, i__3 = *m * ((*m << 1) + 5) + *m * 
				ilaenv_(&c__1, "SORMLQ", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 454 "sgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 456 "sgesvdx.f"
		    }
/* Computing MAX */
#line 457 "sgesvdx.f"
		    i__2 = *m * ((*m << 1) + 19), i__3 = (*m << 2) + *n;
#line 457 "sgesvdx.f"
		    minwrk = max(i__2,i__3);
#line 458 "sgesvdx.f"
		}
#line 459 "sgesvdx.f"
	    }
#line 460 "sgesvdx.f"
	}
#line 461 "sgesvdx.f"
	maxwrk = max(maxwrk,minwrk);
#line 462 "sgesvdx.f"
	work[1] = (doublereal) maxwrk;

#line 464 "sgesvdx.f"
	if (*lwork < minwrk && ! lquery) {
#line 465 "sgesvdx.f"
	    *info = -19;
#line 466 "sgesvdx.f"
	}
#line 467 "sgesvdx.f"
    }

#line 469 "sgesvdx.f"
    if (*info != 0) {
#line 470 "sgesvdx.f"
	i__2 = -(*info);
#line 470 "sgesvdx.f"
	xerbla_("SGESVDX", &i__2, (ftnlen)7);
#line 471 "sgesvdx.f"
	return 0;
#line 472 "sgesvdx.f"
    } else if (lquery) {
#line 473 "sgesvdx.f"
	return 0;
#line 474 "sgesvdx.f"
    }

/*     Quick return if possible */

#line 478 "sgesvdx.f"
    if (*m == 0 || *n == 0) {
#line 479 "sgesvdx.f"
	return 0;
#line 480 "sgesvdx.f"
    }

/*     Set singular values indices accord to RANGE. */

#line 484 "sgesvdx.f"
    if (alls) {
#line 485 "sgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 486 "sgesvdx.f"
	iltgk = 1;
#line 487 "sgesvdx.f"
	iutgk = min(*m,*n);
#line 488 "sgesvdx.f"
    } else if (inds) {
#line 489 "sgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 490 "sgesvdx.f"
	iltgk = *il;
#line 491 "sgesvdx.f"
	iutgk = *iu;
#line 492 "sgesvdx.f"
    } else {
#line 493 "sgesvdx.f"
	*(unsigned char *)rngtgk = 'V';
#line 494 "sgesvdx.f"
	iltgk = 0;
#line 495 "sgesvdx.f"
	iutgk = 0;
#line 496 "sgesvdx.f"
    }

/*     Get machine constants */

#line 500 "sgesvdx.f"
    eps = slamch_("P", (ftnlen)1);
#line 501 "sgesvdx.f"
    smlnum = sqrt(slamch_("S", (ftnlen)1)) / eps;
#line 502 "sgesvdx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 506 "sgesvdx.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 507 "sgesvdx.f"
    iscl = 0;
#line 508 "sgesvdx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 509 "sgesvdx.f"
	iscl = 1;
#line 510 "sgesvdx.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 511 "sgesvdx.f"
    } else if (anrm > bignum) {
#line 512 "sgesvdx.f"
	iscl = 1;
#line 513 "sgesvdx.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 514 "sgesvdx.f"
    }

#line 516 "sgesvdx.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce A using the QR */
/*        decomposition. */

#line 522 "sgesvdx.f"
	if (*m >= mnthr) {

/*           Path 1 (M much larger than N): */
/*           A = Q * R = Q * ( QB * B * PB**T ) */
/*                     = Q * ( QB * ( UB * S * VB**T ) * PB**T ) */
/*           U = Q * QB * UB; V**T = VB**T * PB**T */

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 532 "sgesvdx.f"
	    itau = 1;
#line 533 "sgesvdx.f"
	    itemp = itau + *n;
#line 534 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 534 "sgesvdx.f"
	    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);

/*           Copy R into WORK and bidiagonalize it: */
/*           (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB) */

#line 540 "sgesvdx.f"
	    iqrf = itemp;
#line 541 "sgesvdx.f"
	    id = iqrf + *n * *n;
#line 542 "sgesvdx.f"
	    ie = id + *n;
#line 543 "sgesvdx.f"
	    itauq = ie + *n;
#line 544 "sgesvdx.f"
	    itaup = itauq + *n;
#line 545 "sgesvdx.f"
	    itemp = itaup + *n;
#line 546 "sgesvdx.f"
	    slacpy_("U", n, n, &a[a_offset], lda, &work[iqrf], n, (ftnlen)1);
#line 547 "sgesvdx.f"
	    i__2 = *n - 1;
#line 547 "sgesvdx.f"
	    i__3 = *n - 1;
#line 547 "sgesvdx.f"
	    slaset_("L", &i__2, &i__3, &c_b109, &c_b109, &work[iqrf + 1], n, (
		    ftnlen)1);
#line 548 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 548 "sgesvdx.f"
	    sgebrd_(n, n, &work[iqrf], n, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 555 "sgesvdx.f"
	    itgkz = itemp;
#line 556 "sgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 557 "sgesvdx.f"
	    i__2 = *n << 1;
#line 557 "sgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 563 "sgesvdx.f"
	    if (wantu) {
#line 564 "sgesvdx.f"
		j = itgkz;
#line 565 "sgesvdx.f"
		i__2 = *ns;
#line 565 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 566 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 567 "sgesvdx.f"
		    j += *n << 1;
#line 568 "sgesvdx.f"
		}
#line 569 "sgesvdx.f"
		i__2 = *m - *n;
#line 569 "sgesvdx.f"
		slaset_("A", &i__2, ns, &c_b109, &c_b109, &u[*n + 1 + u_dim1],
			 ldu, (ftnlen)1);

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 574 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 574 "sgesvdx.f"
		sormbr_("Q", "L", "N", n, ns, n, &work[iqrf], n, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call SORMQR to compute Q*(QB*UB). */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 581 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 581 "sgesvdx.f"
		sormqr_("L", "N", m, ns, n, &a[a_offset], lda, &work[itau], &
			u[u_offset], ldu, &work[itemp], &i__2, info, (ftnlen)
			1, (ftnlen)1);
#line 584 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 588 "sgesvdx.f"
	    if (wantvt) {
#line 589 "sgesvdx.f"
		j = itgkz + *n;
#line 590 "sgesvdx.f"
		i__2 = *ns;
#line 590 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 591 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 592 "sgesvdx.f"
		    j += *n << 1;
#line 593 "sgesvdx.f"
		}

/*              Call SORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 598 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 598 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, n, n, &work[iqrf], n, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 601 "sgesvdx.f"
	    }
#line 602 "sgesvdx.f"
	} else {

/*           Path 2 (M at least N, but not much larger) */
/*           Reduce A to bidiagonal form without QR decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB) */

#line 612 "sgesvdx.f"
	    id = 1;
#line 613 "sgesvdx.f"
	    ie = id + *n;
#line 614 "sgesvdx.f"
	    itauq = ie + *n;
#line 615 "sgesvdx.f"
	    itaup = itauq + *n;
#line 616 "sgesvdx.f"
	    itemp = itaup + *n;
#line 617 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 617 "sgesvdx.f"
	    sgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 624 "sgesvdx.f"
	    itgkz = itemp;
#line 625 "sgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 626 "sgesvdx.f"
	    i__2 = *n << 1;
#line 626 "sgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 632 "sgesvdx.f"
	    if (wantu) {
#line 633 "sgesvdx.f"
		j = itgkz;
#line 634 "sgesvdx.f"
		i__2 = *ns;
#line 634 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 635 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 636 "sgesvdx.f"
		    j += *n << 1;
#line 637 "sgesvdx.f"
		}
#line 638 "sgesvdx.f"
		i__2 = *m - *n;
#line 638 "sgesvdx.f"
		slaset_("A", &i__2, ns, &c_b109, &c_b109, &u[*n + 1 + u_dim1],
			 ldu, (ftnlen)1);

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 643 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 643 "sgesvdx.f"
		sormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 646 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 650 "sgesvdx.f"
	    if (wantvt) {
#line 651 "sgesvdx.f"
		j = itgkz + *n;
#line 652 "sgesvdx.f"
		i__2 = *ns;
#line 652 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 653 "sgesvdx.f"
		    scopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 654 "sgesvdx.f"
		    j += *n << 1;
#line 655 "sgesvdx.f"
		}

/*              Call SORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 660 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 660 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 663 "sgesvdx.f"
	    }
#line 664 "sgesvdx.f"
	}
#line 665 "sgesvdx.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce A using the LQ decomposition. */

#line 670 "sgesvdx.f"
	if (*n >= mnthr) {

/*           Path 1t (N much larger than M): */
/*           A = L * Q = ( QB * B * PB**T ) * Q */
/*                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q */
/*           U = QB * UB ; V**T = VB**T * PB**T * Q */

/*           Compute A=L*Q */
/*           (Workspace: need 2*M, prefer M+M*NB) */

#line 680 "sgesvdx.f"
	    itau = 1;
#line 681 "sgesvdx.f"
	    itemp = itau + *m;
#line 682 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 682 "sgesvdx.f"
	    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);
/*           Copy L into WORK and bidiagonalize it: */
/*           (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB) */

#line 688 "sgesvdx.f"
	    ilqf = itemp;
#line 689 "sgesvdx.f"
	    id = ilqf + *m * *m;
#line 690 "sgesvdx.f"
	    ie = id + *m;
#line 691 "sgesvdx.f"
	    itauq = ie + *m;
#line 692 "sgesvdx.f"
	    itaup = itauq + *m;
#line 693 "sgesvdx.f"
	    itemp = itaup + *m;
#line 694 "sgesvdx.f"
	    slacpy_("L", m, m, &a[a_offset], lda, &work[ilqf], m, (ftnlen)1);
#line 695 "sgesvdx.f"
	    i__2 = *m - 1;
#line 695 "sgesvdx.f"
	    i__3 = *m - 1;
#line 695 "sgesvdx.f"
	    slaset_("U", &i__2, &i__3, &c_b109, &c_b109, &work[ilqf + *m], m, 
		    (ftnlen)1);
#line 696 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 696 "sgesvdx.f"
	    sgebrd_(m, m, &work[ilqf], m, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 703 "sgesvdx.f"
	    itgkz = itemp;
#line 704 "sgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 705 "sgesvdx.f"
	    i__2 = *m << 1;
#line 705 "sgesvdx.f"
	    sbdsvdx_("U", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 711 "sgesvdx.f"
	    if (wantu) {
#line 712 "sgesvdx.f"
		j = itgkz;
#line 713 "sgesvdx.f"
		i__2 = *ns;
#line 713 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 714 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 715 "sgesvdx.f"
		    j += *m << 1;
#line 716 "sgesvdx.f"
		}

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 721 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 721 "sgesvdx.f"
		sormbr_("Q", "L", "N", m, ns, m, &work[ilqf], m, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 724 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 728 "sgesvdx.f"
	    if (wantvt) {
#line 729 "sgesvdx.f"
		j = itgkz + *m;
#line 730 "sgesvdx.f"
		i__2 = *ns;
#line 730 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 731 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 732 "sgesvdx.f"
		    j += *m << 1;
#line 733 "sgesvdx.f"
		}
#line 734 "sgesvdx.f"
		i__2 = *n - *m;
#line 734 "sgesvdx.f"
		slaset_("A", ns, &i__2, &c_b109, &c_b109, &vt[(*m + 1) * 
			vt_dim1 + 1], ldvt, (ftnlen)1);

/*              Call SORMBR to compute (VB**T)*(PB**T) */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 739 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 739 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, m, m, &work[ilqf], m, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call SORMLQ to compute ((VB**T)*(PB**T))*Q. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 746 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 746 "sgesvdx.f"
		sormlq_("R", "N", ns, n, m, &a[a_offset], lda, &work[itau], &
			vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1);
#line 749 "sgesvdx.f"
	    }
#line 750 "sgesvdx.f"
	} else {

/*           Path 2t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB) */

#line 760 "sgesvdx.f"
	    id = 1;
#line 761 "sgesvdx.f"
	    ie = id + *m;
#line 762 "sgesvdx.f"
	    itauq = ie + *m;
#line 763 "sgesvdx.f"
	    itaup = itauq + *m;
#line 764 "sgesvdx.f"
	    itemp = itaup + *m;
#line 765 "sgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 765 "sgesvdx.f"
	    sgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 772 "sgesvdx.f"
	    itgkz = itemp;
#line 773 "sgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 774 "sgesvdx.f"
	    i__2 = *m << 1;
#line 774 "sgesvdx.f"
	    sbdsvdx_("L", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 780 "sgesvdx.f"
	    if (wantu) {
#line 781 "sgesvdx.f"
		j = itgkz;
#line 782 "sgesvdx.f"
		i__2 = *ns;
#line 782 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 783 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 784 "sgesvdx.f"
		    j += *m << 1;
#line 785 "sgesvdx.f"
		}

/*              Call SORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 790 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 790 "sgesvdx.f"
		sormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, info, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 793 "sgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 797 "sgesvdx.f"
	    if (wantvt) {
#line 798 "sgesvdx.f"
		j = itgkz + *m;
#line 799 "sgesvdx.f"
		i__2 = *ns;
#line 799 "sgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 800 "sgesvdx.f"
		    scopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 801 "sgesvdx.f"
		    j += *m << 1;
#line 802 "sgesvdx.f"
		}
#line 803 "sgesvdx.f"
		i__2 = *n - *m;
#line 803 "sgesvdx.f"
		slaset_("A", ns, &i__2, &c_b109, &c_b109, &vt[(*m + 1) * 
			vt_dim1 + 1], ldvt, (ftnlen)1);

/*              Call SORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 808 "sgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 808 "sgesvdx.f"
		sormbr_("P", "R", "T", ns, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, 
			info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 811 "sgesvdx.f"
	    }
#line 812 "sgesvdx.f"
	}
#line 813 "sgesvdx.f"
    }

/*     Undo scaling if necessary */

#line 817 "sgesvdx.f"
    if (iscl == 1) {
#line 818 "sgesvdx.f"
	if (anrm > bignum) {
#line 818 "sgesvdx.f"
	    slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 818 "sgesvdx.f"
	}
#line 821 "sgesvdx.f"
	if (anrm < smlnum) {
#line 821 "sgesvdx.f"
	    slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 821 "sgesvdx.f"
	}
#line 824 "sgesvdx.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 828 "sgesvdx.f"
    work[1] = (doublereal) maxwrk;

#line 830 "sgesvdx.f"
    return 0;

/*     End of SGESVDX */

} /* sgesvdx_ */

