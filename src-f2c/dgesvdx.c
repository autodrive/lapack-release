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
static doublereal c_b109 = 0.;

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
/* >          U is DOUBLE PRECISION array, dimension (LDU,UCOL) */
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

/* > \date June 2016 */

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

#line 317 "dgesvdx.f"
    /* Parameter adjustments */
#line 317 "dgesvdx.f"
    a_dim1 = *lda;
#line 317 "dgesvdx.f"
    a_offset = 1 + a_dim1;
#line 317 "dgesvdx.f"
    a -= a_offset;
#line 317 "dgesvdx.f"
    --s;
#line 317 "dgesvdx.f"
    u_dim1 = *ldu;
#line 317 "dgesvdx.f"
    u_offset = 1 + u_dim1;
#line 317 "dgesvdx.f"
    u -= u_offset;
#line 317 "dgesvdx.f"
    vt_dim1 = *ldvt;
#line 317 "dgesvdx.f"
    vt_offset = 1 + vt_dim1;
#line 317 "dgesvdx.f"
    vt -= vt_offset;
#line 317 "dgesvdx.f"
    --work;
#line 317 "dgesvdx.f"
    --iwork;
#line 317 "dgesvdx.f"

#line 317 "dgesvdx.f"
    /* Function Body */
#line 317 "dgesvdx.f"
    *ns = 0;
#line 318 "dgesvdx.f"
    *info = 0;
#line 319 "dgesvdx.f"
    abstol = dlamch_("S", (ftnlen)1) * 2;
#line 320 "dgesvdx.f"
    lquery = *lwork == -1;
#line 321 "dgesvdx.f"
    minmn = min(*m,*n);
#line 323 "dgesvdx.f"
    wantu = lsame_(jobu, "V", (ftnlen)1, (ftnlen)1);
#line 324 "dgesvdx.f"
    wantvt = lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1);
#line 325 "dgesvdx.f"
    if (wantu || wantvt) {
#line 326 "dgesvdx.f"
	*(unsigned char *)jobz = 'V';
#line 327 "dgesvdx.f"
    } else {
#line 328 "dgesvdx.f"
	*(unsigned char *)jobz = 'N';
#line 329 "dgesvdx.f"
    }
#line 330 "dgesvdx.f"
    alls = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 331 "dgesvdx.f"
    vals = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 332 "dgesvdx.f"
    inds = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 334 "dgesvdx.f"
    *info = 0;
#line 335 "dgesvdx.f"
    if (! lsame_(jobu, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobu, "N", (
	    ftnlen)1, (ftnlen)1)) {
#line 337 "dgesvdx.f"
	*info = -1;
#line 338 "dgesvdx.f"
    } else if (! lsame_(jobvt, "V", (ftnlen)1, (ftnlen)1) && ! lsame_(jobvt, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 340 "dgesvdx.f"
	*info = -2;
#line 341 "dgesvdx.f"
    } else if (! (alls || vals || inds)) {
#line 342 "dgesvdx.f"
	*info = -3;
#line 343 "dgesvdx.f"
    } else if (*m < 0) {
#line 344 "dgesvdx.f"
	*info = -4;
#line 345 "dgesvdx.f"
    } else if (*n < 0) {
#line 346 "dgesvdx.f"
	*info = -5;
#line 347 "dgesvdx.f"
    } else if (*m > *lda) {
#line 348 "dgesvdx.f"
	*info = -7;
#line 349 "dgesvdx.f"
    } else if (minmn > 0) {
#line 350 "dgesvdx.f"
	if (vals) {
#line 351 "dgesvdx.f"
	    if (*vl < 0.) {
#line 352 "dgesvdx.f"
		*info = -8;
#line 353 "dgesvdx.f"
	    } else if (*vu <= *vl) {
#line 354 "dgesvdx.f"
		*info = -9;
#line 355 "dgesvdx.f"
	    }
#line 356 "dgesvdx.f"
	} else if (inds) {
#line 357 "dgesvdx.f"
	    if (*il < 1 || *il > max(1,minmn)) {
#line 358 "dgesvdx.f"
		*info = -10;
#line 359 "dgesvdx.f"
	    } else if (*iu < min(minmn,*il) || *iu > minmn) {
#line 360 "dgesvdx.f"
		*info = -11;
#line 361 "dgesvdx.f"
	    }
#line 362 "dgesvdx.f"
	}
#line 363 "dgesvdx.f"
	if (*info == 0) {
#line 364 "dgesvdx.f"
	    if (wantu && *ldu < *m) {
#line 365 "dgesvdx.f"
		*info = -15;
#line 366 "dgesvdx.f"
	    } else if (wantvt) {
#line 367 "dgesvdx.f"
		if (inds) {
#line 368 "dgesvdx.f"
		    if (*ldvt < *iu - *il + 1) {
#line 369 "dgesvdx.f"
			*info = -17;
#line 370 "dgesvdx.f"
		    }
#line 371 "dgesvdx.f"
		} else if (*ldvt < minmn) {
#line 372 "dgesvdx.f"
		    *info = -17;
#line 373 "dgesvdx.f"
		}
#line 374 "dgesvdx.f"
	    }
#line 375 "dgesvdx.f"
	}
#line 376 "dgesvdx.f"
    }

/*     Compute workspace */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 385 "dgesvdx.f"
    if (*info == 0) {
#line 386 "dgesvdx.f"
	minwrk = 1;
#line 387 "dgesvdx.f"
	maxwrk = 1;
#line 388 "dgesvdx.f"
	if (minmn > 0) {
#line 389 "dgesvdx.f"
	    if (*m >= *n) {
/* Writing concatenation */
#line 390 "dgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 390 "dgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 390 "dgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 390 "dgesvdx.f"
		mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 391 "dgesvdx.f"
		if (*m >= mnthr) {

/*                 Path 1 (M much larger than N) */

#line 395 "dgesvdx.f"
		    maxwrk = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 397 "dgesvdx.f"
		    i__2 = maxwrk, i__3 = *n * (*n + 5) + (*n << 1) * ilaenv_(
			    &c__1, "DGEBRD", " ", n, n, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 397 "dgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 399 "dgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 400 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *n * (*n * 3 + 6) + *n * 
				ilaenv_(&c__1, "DORMQR", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 400 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 402 "dgesvdx.f"
		    }
#line 403 "dgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 404 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *n * (*n * 3 + 6) + *n * 
				ilaenv_(&c__1, "DORMLQ", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 404 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 406 "dgesvdx.f"
		    }
#line 407 "dgesvdx.f"
		    minwrk = *n * (*n * 3 + 20);
#line 408 "dgesvdx.f"
		} else {

/*                 Path 2 (M at least N, but not much larger) */

#line 412 "dgesvdx.f"
		    maxwrk = (*n << 2) + (*m + *n) * ilaenv_(&c__1, "DGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 414 "dgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 415 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *n * ((*n << 1) + 5) + *n * 
				ilaenv_(&c__1, "DORMQR", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 415 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 417 "dgesvdx.f"
		    }
#line 418 "dgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 419 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *n * ((*n << 1) + 5) + *n * 
				ilaenv_(&c__1, "DORMLQ", " ", n, n, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 419 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 421 "dgesvdx.f"
		    }
/* Computing MAX */
#line 422 "dgesvdx.f"
		    i__2 = *n * ((*n << 1) + 19), i__3 = (*n << 2) + *m;
#line 422 "dgesvdx.f"
		    minwrk = max(i__2,i__3);
#line 423 "dgesvdx.f"
		}
#line 424 "dgesvdx.f"
	    } else {
/* Writing concatenation */
#line 425 "dgesvdx.f"
		i__1[0] = 1, a__1[0] = jobu;
#line 425 "dgesvdx.f"
		i__1[1] = 1, a__1[1] = jobvt;
#line 425 "dgesvdx.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 425 "dgesvdx.f"
		mnthr = ilaenv_(&c__6, "DGESVD", ch__1, m, n, &c__0, &c__0, (
			ftnlen)6, (ftnlen)2);
#line 426 "dgesvdx.f"
		if (*n >= mnthr) {

/*                 Path 1t (N much larger than M) */

#line 430 "dgesvdx.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 432 "dgesvdx.f"
		    i__2 = maxwrk, i__3 = *m * (*m + 5) + (*m << 1) * ilaenv_(
			    &c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1, (ftnlen)
			    6, (ftnlen)1);
#line 432 "dgesvdx.f"
		    maxwrk = max(i__2,i__3);
#line 434 "dgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 435 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *m * (*m * 3 + 6) + *m * 
				ilaenv_(&c__1, "DORMQR", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 435 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 437 "dgesvdx.f"
		    }
#line 438 "dgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 439 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *m * (*m * 3 + 6) + *m * 
				ilaenv_(&c__1, "DORMLQ", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 439 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 441 "dgesvdx.f"
		    }
#line 442 "dgesvdx.f"
		    minwrk = *m * (*m * 3 + 20);
#line 443 "dgesvdx.f"
		} else {

/*                 Path 2t (N at least M, but not much larger) */

#line 447 "dgesvdx.f"
		    maxwrk = (*m << 2) + (*m + *n) * ilaenv_(&c__1, "DGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 449 "dgesvdx.f"
		    if (wantu) {
/* Computing MAX */
#line 450 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *m * ((*m << 1) + 5) + *m * 
				ilaenv_(&c__1, "DORMQR", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 450 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 452 "dgesvdx.f"
		    }
#line 453 "dgesvdx.f"
		    if (wantvt) {
/* Computing MAX */
#line 454 "dgesvdx.f"
			i__2 = maxwrk, i__3 = *m * ((*m << 1) + 5) + *m * 
				ilaenv_(&c__1, "DORMLQ", " ", m, m, &c_n1, &
				c_n1, (ftnlen)6, (ftnlen)1);
#line 454 "dgesvdx.f"
			maxwrk = max(i__2,i__3);
#line 456 "dgesvdx.f"
		    }
/* Computing MAX */
#line 457 "dgesvdx.f"
		    i__2 = *m * ((*m << 1) + 19), i__3 = (*m << 2) + *n;
#line 457 "dgesvdx.f"
		    minwrk = max(i__2,i__3);
#line 458 "dgesvdx.f"
		}
#line 459 "dgesvdx.f"
	    }
#line 460 "dgesvdx.f"
	}
#line 461 "dgesvdx.f"
	maxwrk = max(maxwrk,minwrk);
#line 462 "dgesvdx.f"
	work[1] = (doublereal) maxwrk;

#line 464 "dgesvdx.f"
	if (*lwork < minwrk && ! lquery) {
#line 465 "dgesvdx.f"
	    *info = -19;
#line 466 "dgesvdx.f"
	}
#line 467 "dgesvdx.f"
    }

#line 469 "dgesvdx.f"
    if (*info != 0) {
#line 470 "dgesvdx.f"
	i__2 = -(*info);
#line 470 "dgesvdx.f"
	xerbla_("DGESVDX", &i__2, (ftnlen)7);
#line 471 "dgesvdx.f"
	return 0;
#line 472 "dgesvdx.f"
    } else if (lquery) {
#line 473 "dgesvdx.f"
	return 0;
#line 474 "dgesvdx.f"
    }

/*     Quick return if possible */

#line 478 "dgesvdx.f"
    if (*m == 0 || *n == 0) {
#line 479 "dgesvdx.f"
	return 0;
#line 480 "dgesvdx.f"
    }

/*     Set singular values indices accord to RANGE. */

#line 484 "dgesvdx.f"
    if (alls) {
#line 485 "dgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 486 "dgesvdx.f"
	iltgk = 1;
#line 487 "dgesvdx.f"
	iutgk = min(*m,*n);
#line 488 "dgesvdx.f"
    } else if (inds) {
#line 489 "dgesvdx.f"
	*(unsigned char *)rngtgk = 'I';
#line 490 "dgesvdx.f"
	iltgk = *il;
#line 491 "dgesvdx.f"
	iutgk = *iu;
#line 492 "dgesvdx.f"
    } else {
#line 493 "dgesvdx.f"
	*(unsigned char *)rngtgk = 'V';
#line 494 "dgesvdx.f"
	iltgk = 0;
#line 495 "dgesvdx.f"
	iutgk = 0;
#line 496 "dgesvdx.f"
    }

/*     Get machine constants */

#line 500 "dgesvdx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 501 "dgesvdx.f"
    smlnum = sqrt(dlamch_("S", (ftnlen)1)) / eps;
#line 502 "dgesvdx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 506 "dgesvdx.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 507 "dgesvdx.f"
    iscl = 0;
#line 508 "dgesvdx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 509 "dgesvdx.f"
	iscl = 1;
#line 510 "dgesvdx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 511 "dgesvdx.f"
    } else if (anrm > bignum) {
#line 512 "dgesvdx.f"
	iscl = 1;
#line 513 "dgesvdx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 514 "dgesvdx.f"
    }

#line 516 "dgesvdx.f"
    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce A using the QR */
/*        decomposition. */

#line 522 "dgesvdx.f"
	if (*m >= mnthr) {

/*           Path 1 (M much larger than N): */
/*           A = Q * R = Q * ( QB * B * PB**T ) */
/*                     = Q * ( QB * ( UB * S * VB**T ) * PB**T ) */
/*           U = Q * QB * UB; V**T = VB**T * PB**T */

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 532 "dgesvdx.f"
	    itau = 1;
#line 533 "dgesvdx.f"
	    itemp = itau + *n;
#line 534 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 534 "dgesvdx.f"
	    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);

/*           Copy R into WORK and bidiagonalize it: */
/*           (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB) */

#line 540 "dgesvdx.f"
	    iqrf = itemp;
#line 541 "dgesvdx.f"
	    id = iqrf + *n * *n;
#line 542 "dgesvdx.f"
	    ie = id + *n;
#line 543 "dgesvdx.f"
	    itauq = ie + *n;
#line 544 "dgesvdx.f"
	    itaup = itauq + *n;
#line 545 "dgesvdx.f"
	    itemp = itaup + *n;
#line 546 "dgesvdx.f"
	    dlacpy_("U", n, n, &a[a_offset], lda, &work[iqrf], n, (ftnlen)1);
#line 547 "dgesvdx.f"
	    i__2 = *n - 1;
#line 547 "dgesvdx.f"
	    i__3 = *n - 1;
#line 547 "dgesvdx.f"
	    dlaset_("L", &i__2, &i__3, &c_b109, &c_b109, &work[iqrf + 1], n, (
		    ftnlen)1);
#line 548 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 548 "dgesvdx.f"
	    dgebrd_(n, n, &work[iqrf], n, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 555 "dgesvdx.f"
	    itgkz = itemp;
#line 556 "dgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 557 "dgesvdx.f"
	    i__2 = *n << 1;
#line 557 "dgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 563 "dgesvdx.f"
	    if (wantu) {
#line 564 "dgesvdx.f"
		j = itgkz;
#line 565 "dgesvdx.f"
		i__2 = *ns;
#line 565 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 566 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 567 "dgesvdx.f"
		    j += *n << 1;
#line 568 "dgesvdx.f"
		}
#line 569 "dgesvdx.f"
		i__2 = *m - *n;
#line 569 "dgesvdx.f"
		dlaset_("A", &i__2, ns, &c_b109, &c_b109, &u[*n + 1 + u_dim1],
			 ldu, (ftnlen)1);

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 574 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 574 "dgesvdx.f"
		dormbr_("Q", "L", "N", n, ns, n, &work[iqrf], n, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call DORMQR to compute Q*(QB*UB). */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 581 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 581 "dgesvdx.f"
		dormqr_("L", "N", m, ns, n, &a[a_offset], lda, &work[itau], &
			u[u_offset], ldu, &work[itemp], &i__2, info, (ftnlen)
			1, (ftnlen)1);
#line 584 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 588 "dgesvdx.f"
	    if (wantvt) {
#line 589 "dgesvdx.f"
		j = itgkz + *n;
#line 590 "dgesvdx.f"
		i__2 = *ns;
#line 590 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 591 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 592 "dgesvdx.f"
		    j += *n << 1;
#line 593 "dgesvdx.f"
		}

/*              Call DORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 598 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 598 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, n, n, &work[iqrf], n, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 601 "dgesvdx.f"
	    }
#line 602 "dgesvdx.f"
	} else {

/*           Path 2 (M at least N, but not much larger) */
/*           Reduce A to bidiagonal form without QR decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB) */

#line 612 "dgesvdx.f"
	    id = 1;
#line 613 "dgesvdx.f"
	    ie = id + *n;
#line 614 "dgesvdx.f"
	    itauq = ie + *n;
#line 615 "dgesvdx.f"
	    itaup = itauq + *n;
#line 616 "dgesvdx.f"
	    itemp = itaup + *n;
#line 617 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 617 "dgesvdx.f"
	    dgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 14*N + 2*N*(N+1)) */

#line 624 "dgesvdx.f"
	    itgkz = itemp;
#line 625 "dgesvdx.f"
	    itemp = itgkz + *n * ((*n << 1) + 1);
#line 626 "dgesvdx.f"
	    i__2 = *n << 1;
#line 626 "dgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, n, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 632 "dgesvdx.f"
	    if (wantu) {
#line 633 "dgesvdx.f"
		j = itgkz;
#line 634 "dgesvdx.f"
		i__2 = *ns;
#line 634 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 635 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 636 "dgesvdx.f"
		    j += *n << 1;
#line 637 "dgesvdx.f"
		}
#line 638 "dgesvdx.f"
		i__2 = *m - *n;
#line 638 "dgesvdx.f"
		dlaset_("A", &i__2, ns, &c_b109, &c_b109, &u[*n + 1 + u_dim1],
			 ldu, (ftnlen)1);

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 643 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 643 "dgesvdx.f"
		dormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, &ierr,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 646 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 650 "dgesvdx.f"
	    if (wantvt) {
#line 651 "dgesvdx.f"
		j = itgkz + *n;
#line 652 "dgesvdx.f"
		i__2 = *ns;
#line 652 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 653 "dgesvdx.f"
		    dcopy_(n, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 654 "dgesvdx.f"
		    j += *n << 1;
#line 655 "dgesvdx.f"
		}

/*              Call DORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need N, prefer N*NB) */

#line 660 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 660 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, &
			ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 663 "dgesvdx.f"
	    }
#line 664 "dgesvdx.f"
	}
#line 665 "dgesvdx.f"
    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce A using the LQ decomposition. */

#line 670 "dgesvdx.f"
	if (*n >= mnthr) {

/*           Path 1t (N much larger than M): */
/*           A = L * Q = ( QB * B * PB**T ) * Q */
/*                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q */
/*           U = QB * UB ; V**T = VB**T * PB**T * Q */

/*           Compute A=L*Q */
/*           (Workspace: need 2*M, prefer M+M*NB) */

#line 680 "dgesvdx.f"
	    itau = 1;
#line 681 "dgesvdx.f"
	    itemp = itau + *m;
#line 682 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 682 "dgesvdx.f"
	    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2,
		     info);
/*           Copy L into WORK and bidiagonalize it: */
/*           (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB) */

#line 688 "dgesvdx.f"
	    ilqf = itemp;
#line 689 "dgesvdx.f"
	    id = ilqf + *m * *m;
#line 690 "dgesvdx.f"
	    ie = id + *m;
#line 691 "dgesvdx.f"
	    itauq = ie + *m;
#line 692 "dgesvdx.f"
	    itaup = itauq + *m;
#line 693 "dgesvdx.f"
	    itemp = itaup + *m;
#line 694 "dgesvdx.f"
	    dlacpy_("L", m, m, &a[a_offset], lda, &work[ilqf], m, (ftnlen)1);
#line 695 "dgesvdx.f"
	    i__2 = *m - 1;
#line 695 "dgesvdx.f"
	    i__3 = *m - 1;
#line 695 "dgesvdx.f"
	    dlaset_("U", &i__2, &i__3, &c_b109, &c_b109, &work[ilqf + *m], m, 
		    (ftnlen)1);
#line 696 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 696 "dgesvdx.f"
	    dgebrd_(m, m, &work[ilqf], m, &work[id], &work[ie], &work[itauq], 
		    &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 703 "dgesvdx.f"
	    itgkz = itemp;
#line 704 "dgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 705 "dgesvdx.f"
	    i__2 = *m << 1;
#line 705 "dgesvdx.f"
	    dbdsvdx_("U", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 711 "dgesvdx.f"
	    if (wantu) {
#line 712 "dgesvdx.f"
		j = itgkz;
#line 713 "dgesvdx.f"
		i__2 = *ns;
#line 713 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 714 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 715 "dgesvdx.f"
		    j += *m << 1;
#line 716 "dgesvdx.f"
		}

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 721 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 721 "dgesvdx.f"
		dormbr_("Q", "L", "N", m, ns, m, &work[ilqf], m, &work[itauq],
			 &u[u_offset], ldu, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 724 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 728 "dgesvdx.f"
	    if (wantvt) {
#line 729 "dgesvdx.f"
		j = itgkz + *m;
#line 730 "dgesvdx.f"
		i__2 = *ns;
#line 730 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 731 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 732 "dgesvdx.f"
		    j += *m << 1;
#line 733 "dgesvdx.f"
		}
#line 734 "dgesvdx.f"
		i__2 = *n - *m;
#line 734 "dgesvdx.f"
		dlaset_("A", ns, &i__2, &c_b109, &c_b109, &vt[(*m + 1) * 
			vt_dim1 + 1], ldvt, (ftnlen)1);

/*              Call DORMBR to compute (VB**T)*(PB**T) */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 739 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 739 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, m, m, &work[ilqf], m, &work[itaup],
			 &vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);

/*              Call DORMLQ to compute ((VB**T)*(PB**T))*Q. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 746 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 746 "dgesvdx.f"
		dormlq_("R", "N", ns, n, m, &a[a_offset], lda, &work[itau], &
			vt[vt_offset], ldvt, &work[itemp], &i__2, info, (
			ftnlen)1, (ftnlen)1);
#line 749 "dgesvdx.f"
	    }
#line 750 "dgesvdx.f"
	} else {

/*           Path 2t (N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */
/*           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
/*           U = QB * UB; V**T = VB**T * PB**T */

/*           Bidiagonalize A */
/*           (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB) */

#line 760 "dgesvdx.f"
	    id = 1;
#line 761 "dgesvdx.f"
	    ie = id + *m;
#line 762 "dgesvdx.f"
	    itauq = ie + *m;
#line 763 "dgesvdx.f"
	    itaup = itauq + *m;
#line 764 "dgesvdx.f"
	    itemp = itaup + *m;
#line 765 "dgesvdx.f"
	    i__2 = *lwork - itemp + 1;
#line 765 "dgesvdx.f"
	    dgebrd_(m, n, &a[a_offset], lda, &work[id], &work[ie], &work[
		    itauq], &work[itaup], &work[itemp], &i__2, info);

/*           Solve eigenvalue problem TGK*Z=Z*S. */
/*           (Workspace: need 2*M*M+14*M) */

#line 772 "dgesvdx.f"
	    itgkz = itemp;
#line 773 "dgesvdx.f"
	    itemp = itgkz + *m * ((*m << 1) + 1);
#line 774 "dgesvdx.f"
	    i__2 = *m << 1;
#line 774 "dgesvdx.f"
	    dbdsvdx_("L", jobz, rngtgk, m, &work[id], &work[ie], vl, vu, &
		    iltgk, &iutgk, ns, &s[1], &work[itgkz], &i__2, &work[
		    itemp], &iwork[1], info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           If needed, compute left singular vectors. */

#line 780 "dgesvdx.f"
	    if (wantu) {
#line 781 "dgesvdx.f"
		j = itgkz;
#line 782 "dgesvdx.f"
		i__2 = *ns;
#line 782 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 783 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &u[i__ * u_dim1 + 1], &c__1);
#line 784 "dgesvdx.f"
		    j += *m << 1;
#line 785 "dgesvdx.f"
		}

/*              Call DORMBR to compute QB*UB. */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 790 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 790 "dgesvdx.f"
		dormbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[itemp], &i__2, info, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 793 "dgesvdx.f"
	    }

/*           If needed, compute right singular vectors. */

#line 797 "dgesvdx.f"
	    if (wantvt) {
#line 798 "dgesvdx.f"
		j = itgkz + *m;
#line 799 "dgesvdx.f"
		i__2 = *ns;
#line 799 "dgesvdx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 800 "dgesvdx.f"
		    dcopy_(m, &work[j], &c__1, &vt[i__ + vt_dim1], ldvt);
#line 801 "dgesvdx.f"
		    j += *m << 1;
#line 802 "dgesvdx.f"
		}
#line 803 "dgesvdx.f"
		i__2 = *n - *m;
#line 803 "dgesvdx.f"
		dlaset_("A", ns, &i__2, &c_b109, &c_b109, &vt[(*m + 1) * 
			vt_dim1 + 1], ldvt, (ftnlen)1);

/*              Call DORMBR to compute VB**T * PB**T */
/*              (Workspace in WORK( ITEMP ): need M, prefer M*NB) */

#line 808 "dgesvdx.f"
		i__2 = *lwork - itemp + 1;
#line 808 "dgesvdx.f"
		dormbr_("P", "R", "T", ns, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, 
			info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 811 "dgesvdx.f"
	    }
#line 812 "dgesvdx.f"
	}
#line 813 "dgesvdx.f"
    }

/*     Undo scaling if necessary */

#line 817 "dgesvdx.f"
    if (iscl == 1) {
#line 818 "dgesvdx.f"
	if (anrm > bignum) {
#line 818 "dgesvdx.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 818 "dgesvdx.f"
	}
#line 821 "dgesvdx.f"
	if (anrm < smlnum) {
#line 821 "dgesvdx.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, info, (ftnlen)1);
#line 821 "dgesvdx.f"
	}
#line 824 "dgesvdx.f"
    }

/*     Return optimal workspace in WORK(1) */

#line 828 "dgesvdx.f"
    work[1] = (doublereal) maxwrk;

#line 830 "dgesvdx.f"
    return 0;

/*     End of DGESVDX */

} /* dgesvdx_ */

