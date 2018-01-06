#line 1 "ssytrd_sy2sb.f"
/* ssytrd_sy2sb.f -- translated by f2c (version 20100827).
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

#line 1 "ssytrd_sy2sb.f"
/* Table of constant values */

static integer c__20 = 20;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b17 = 0.;
static doublereal c_b23 = 1.;
static doublereal c_b39 = -.5;
static doublereal c_b42 = -1.;

/* > \brief \b SSYTRD_SY2SB */

/*  @generated from zhetrd_he2hb.f, fortran z -> s, Wed Dec  7 08:22:40 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRD_SY2SB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTRD_SY2SB( UPLO, N, KD, A, LDA, AB, LDAB, TAU, */
/*                              WORK, LWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAB, LWORK, N, KD */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), AB( LDAB, * ), */
/*                          TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRD_SY2SB reduces a real symmetric matrix A to real symmetric */
/* > band-diagonal form AB by a orthogonal similarity transformation: */
/* > Q**T * A * Q = AB. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the reduced matrix if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* >          The reduced matrix is stored in the array AB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* >          of A are overwritten by the corresponding elements of the */
/* >          tridiagonal matrix T, and the elements above the first */
/* >          superdiagonal, with the array TAU, represent the orthogonal */
/* >          matrix Q as a product of elementary reflectors; if UPLO */
/* >          = 'L', the diagonal and first subdiagonal of A are over- */
/* >          written by the corresponding elements of the tridiagonal */
/* >          matrix T, and the elements below the first subdiagonal, with */
/* >          the array TAU, represent the orthogonal matrix Q as a product */
/* >          of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          On exit, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (N-KD) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (LWORK) */
/* >          On exit, if INFO = 0, or if LWORK=-1, */
/* >          WORK(1) returns the size of LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK which should be calculated */
/* >          by a workspace query. LWORK = MAX(1, LWORK_QUERY) */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* >          LWORK_QUERY = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD */
/* >          where FACTOPTNB is the blocking used by the QR or LQ */
/* >          algorithm, usually FACTOPTNB=128 is a good choice otherwise */
/* >          putting LWORK=-1 will provide the size of WORK. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup realSYcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Implemented by Azzam Haidar. */
/* > */
/* >  All details are available on technical report, SC11, SC13 papers. */
/* > */
/* >  Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* >  Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* >  using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* >  of 2011 International Conference for High Performance Computing, */
/* >  Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* >  Article 8 , 11 pages. */
/* >  http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* >  A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* >  An improved parallel singular value algorithm and its implementation */
/* >  for multicore hardware, In Proceedings of 2013 International Conference */
/* >  for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* >  Denver, Colorado, USA, 2013. */
/* >  Article 90, 12 pages. */
/* >  http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* >  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* >  A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* >  calculations based on fine-grained memory aware tasks. */
/* >  International Journal of High Performance Computing Applications. */
/* >  Volume 28 Issue 2, Pages 196-209, May 2014. */
/* >  http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* > */
/* > \verbatim */
/* > */
/* >  If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(k)**T . . . H(2)**T H(1)**T, where k = n-kd. */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i+kd-1) = 0 and v(i+kd) = 1; conjg(v(i+kd+1:n)) is stored on exit in */
/* >  A(i,i+kd+1:n), and tau in TAU(i). */
/* > */
/* >  If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = n-kd. */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(kd+1:i) = 0 and v(i+kd+1) = 1; v(i+kd+2:n) is stored on exit in */
/* >  A(i+kd+2:n,i), and tau in TAU(i). */
/* > */
/* >  The contents of A on exit are illustrated by the following examples */
/* >  with n = 5: */
/* > */
/* >  if UPLO = 'U':                       if UPLO = 'L': */
/* > */
/* >    (  ab  ab/v1  v1      v1     v1    )              (  ab                            ) */
/* >    (      ab     ab/v2   v2     v2    )              (  ab/v1  ab                     ) */
/* >    (             ab      ab/v3  v3    )              (  v1     ab/v2  ab              ) */
/* >    (                     ab     ab/v4 )              (  v1     v2     ab/v3  ab       ) */
/* >    (                            ab    )              (  v1     v2     v3     ab/v4 ab ) */
/* > */
/* >  where d and e denote diagonal and off-diagonal elements of T, and vi */
/* >  denotes an element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ssytrd_sy2sb__(char *uplo, integer *n, integer *kd, 
	doublereal *a, integer *lda, doublereal *ab, integer *ldab, 
	doublereal *tau, doublereal *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, ab_dim1, ab_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Local variables */
    static integer i__, j, lk, pk, pn, lt, lw, ls1, ls2, ldt, ldw, lds1, lds2,
	     tpos, wpos, s1pos, s2pos;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer lwmin;
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), ssymm_(char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     ssyr2k_(char *, char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), xerbla_(char *, integer 
	    *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    sgeqrf_(integer *, integer *, doublereal *, integer *, doublereal 
	    *, doublereal *, integer *, integer *), slarft_(char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), slaset_(char *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, integer *,
	     ftnlen);
    static logical lquery;



/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Determine the minimal workspace size required */
/*     and test the input parameters */

#line 296 "ssytrd_sy2sb.f"
    /* Parameter adjustments */
#line 296 "ssytrd_sy2sb.f"
    a_dim1 = *lda;
#line 296 "ssytrd_sy2sb.f"
    a_offset = 1 + a_dim1;
#line 296 "ssytrd_sy2sb.f"
    a -= a_offset;
#line 296 "ssytrd_sy2sb.f"
    ab_dim1 = *ldab;
#line 296 "ssytrd_sy2sb.f"
    ab_offset = 1 + ab_dim1;
#line 296 "ssytrd_sy2sb.f"
    ab -= ab_offset;
#line 296 "ssytrd_sy2sb.f"
    --tau;
#line 296 "ssytrd_sy2sb.f"
    --work;
#line 296 "ssytrd_sy2sb.f"

#line 296 "ssytrd_sy2sb.f"
    /* Function Body */
#line 296 "ssytrd_sy2sb.f"
    *info = 0;
#line 297 "ssytrd_sy2sb.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 298 "ssytrd_sy2sb.f"
    lquery = *lwork == -1;
#line 299 "ssytrd_sy2sb.f"
    lwmin = ilaenv_(&c__20, "SSYTRD_SY2SB", "", n, kd, &c_n1, &c_n1, (ftnlen)
	    12, (ftnlen)0);
#line 301 "ssytrd_sy2sb.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 302 "ssytrd_sy2sb.f"
	*info = -1;
#line 303 "ssytrd_sy2sb.f"
    } else if (*n < 0) {
#line 304 "ssytrd_sy2sb.f"
	*info = -2;
#line 305 "ssytrd_sy2sb.f"
    } else if (*kd < 0) {
#line 306 "ssytrd_sy2sb.f"
	*info = -3;
#line 307 "ssytrd_sy2sb.f"
    } else if (*lda < max(1,*n)) {
#line 308 "ssytrd_sy2sb.f"
	*info = -5;
#line 309 "ssytrd_sy2sb.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 309 "ssytrd_sy2sb.f"
	i__1 = 1, i__2 = *kd + 1;
#line 309 "ssytrd_sy2sb.f"
	if (*ldab < max(i__1,i__2)) {
#line 310 "ssytrd_sy2sb.f"
	    *info = -7;
#line 311 "ssytrd_sy2sb.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 312 "ssytrd_sy2sb.f"
	    *info = -10;
#line 313 "ssytrd_sy2sb.f"
	}
#line 313 "ssytrd_sy2sb.f"
    }

#line 315 "ssytrd_sy2sb.f"
    if (*info != 0) {
#line 316 "ssytrd_sy2sb.f"
	i__1 = -(*info);
#line 316 "ssytrd_sy2sb.f"
	xerbla_("SSYTRD_SY2SB", &i__1, (ftnlen)12);
#line 317 "ssytrd_sy2sb.f"
	return 0;
#line 318 "ssytrd_sy2sb.f"
    } else if (lquery) {
#line 319 "ssytrd_sy2sb.f"
	work[1] = (doublereal) lwmin;
#line 320 "ssytrd_sy2sb.f"
	return 0;
#line 321 "ssytrd_sy2sb.f"
    }

/*     Quick return if possible */
/*     Copy the upper/lower portion of A into AB */

#line 326 "ssytrd_sy2sb.f"
    if (*n <= *kd + 1) {
#line 327 "ssytrd_sy2sb.f"
	if (upper) {
#line 328 "ssytrd_sy2sb.f"
	    i__1 = *n;
#line 328 "ssytrd_sy2sb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 329 "ssytrd_sy2sb.f"
		i__2 = *kd + 1;
#line 329 "ssytrd_sy2sb.f"
		lk = min(i__2,i__);
#line 330 "ssytrd_sy2sb.f"
		scopy_(&lk, &a[i__ - lk + 1 + i__ * a_dim1], &c__1, &ab[*kd + 
			1 - lk + 1 + i__ * ab_dim1], &c__1);
#line 332 "ssytrd_sy2sb.f"
/* L100: */
#line 332 "ssytrd_sy2sb.f"
	    }
#line 333 "ssytrd_sy2sb.f"
	} else {
#line 334 "ssytrd_sy2sb.f"
	    i__1 = *n;
#line 334 "ssytrd_sy2sb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 335 "ssytrd_sy2sb.f"
		i__2 = *kd + 1, i__3 = *n - i__ + 1;
#line 335 "ssytrd_sy2sb.f"
		lk = min(i__2,i__3);
#line 336 "ssytrd_sy2sb.f"
		scopy_(&lk, &a[i__ + i__ * a_dim1], &c__1, &ab[i__ * ab_dim1 
			+ 1], &c__1);
#line 337 "ssytrd_sy2sb.f"
/* L110: */
#line 337 "ssytrd_sy2sb.f"
	    }
#line 338 "ssytrd_sy2sb.f"
	}
#line 339 "ssytrd_sy2sb.f"
	work[1] = 1.;
#line 340 "ssytrd_sy2sb.f"
	return 0;
#line 341 "ssytrd_sy2sb.f"
    }

/*     Determine the pointer position for the workspace */

#line 345 "ssytrd_sy2sb.f"
    ldt = *kd;
#line 346 "ssytrd_sy2sb.f"
    lds1 = *kd;
#line 347 "ssytrd_sy2sb.f"
    lt = ldt * *kd;
#line 348 "ssytrd_sy2sb.f"
    lw = *n * *kd;
#line 349 "ssytrd_sy2sb.f"
    ls1 = lds1 * *kd;
#line 350 "ssytrd_sy2sb.f"
    ls2 = lwmin - lt - lw - ls1;
/*      LS2 = N*MAX(KD,FACTOPTNB) */
#line 352 "ssytrd_sy2sb.f"
    tpos = 1;
#line 353 "ssytrd_sy2sb.f"
    wpos = tpos + lt;
#line 354 "ssytrd_sy2sb.f"
    s1pos = wpos + lw;
#line 355 "ssytrd_sy2sb.f"
    s2pos = s1pos + ls1;
#line 356 "ssytrd_sy2sb.f"
    if (upper) {
#line 357 "ssytrd_sy2sb.f"
	ldw = *kd;
#line 358 "ssytrd_sy2sb.f"
	lds2 = *kd;
#line 359 "ssytrd_sy2sb.f"
    } else {
#line 360 "ssytrd_sy2sb.f"
	ldw = *n;
#line 361 "ssytrd_sy2sb.f"
	lds2 = *n;
#line 362 "ssytrd_sy2sb.f"
    }


/*     Set the workspace of the triangular matrix T to zero once such a */
/*     way everytime T is generated the upper/lower portion will be always zero */

#line 368 "ssytrd_sy2sb.f"
    slaset_("A", &ldt, kd, &c_b17, &c_b17, &work[tpos], &ldt, (ftnlen)1);

#line 370 "ssytrd_sy2sb.f"
    if (upper) {
#line 371 "ssytrd_sy2sb.f"
	i__1 = *n - *kd;
#line 371 "ssytrd_sy2sb.f"
	i__2 = *kd;
#line 371 "ssytrd_sy2sb.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 372 "ssytrd_sy2sb.f"
	    pn = *n - i__ - *kd + 1;
/* Computing MIN */
#line 373 "ssytrd_sy2sb.f"
	    i__3 = *n - i__ - *kd + 1;
#line 373 "ssytrd_sy2sb.f"
	    pk = min(i__3,*kd);

/*            Compute the LQ factorization of the current block */

#line 377 "ssytrd_sy2sb.f"
	    sgelqf_(kd, &pn, &a[i__ + (i__ + *kd) * a_dim1], lda, &tau[i__], &
		    work[s2pos], &ls2, &iinfo);

/*            Copy the upper portion of A into AB */

#line 382 "ssytrd_sy2sb.f"
	    i__3 = i__ + pk - 1;
#line 382 "ssytrd_sy2sb.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MIN */
#line 383 "ssytrd_sy2sb.f"
		i__4 = *kd, i__5 = *n - j;
#line 383 "ssytrd_sy2sb.f"
		lk = min(i__4,i__5) + 1;
#line 384 "ssytrd_sy2sb.f"
		i__4 = *ldab - 1;
#line 384 "ssytrd_sy2sb.f"
		scopy_(&lk, &a[j + j * a_dim1], lda, &ab[*kd + 1 + j * 
			ab_dim1], &i__4);
#line 385 "ssytrd_sy2sb.f"
/* L20: */
#line 385 "ssytrd_sy2sb.f"
	    }

#line 387 "ssytrd_sy2sb.f"
	    slaset_("Lower", &pk, &pk, &c_b17, &c_b23, &a[i__ + (i__ + *kd) * 
		    a_dim1], lda, (ftnlen)5);

/*            Form the matrix T */

#line 392 "ssytrd_sy2sb.f"
	    slarft_("Forward", "Rowwise", &pn, &pk, &a[i__ + (i__ + *kd) * 
		    a_dim1], lda, &tau[i__], &work[tpos], &ldt, (ftnlen)7, (
		    ftnlen)7);

/*            Compute W: */

#line 398 "ssytrd_sy2sb.f"
	    sgemm_("Conjugate", "No transpose", &pk, &pn, &pk, &c_b23, &work[
		    tpos], &ldt, &a[i__ + (i__ + *kd) * a_dim1], lda, &c_b17, 
		    &work[s2pos], &lds2, (ftnlen)9, (ftnlen)12);

#line 403 "ssytrd_sy2sb.f"
	    ssymm_("Right", uplo, &pk, &pn, &c_b23, &a[i__ + *kd + (i__ + *kd)
		     * a_dim1], lda, &work[s2pos], &lds2, &c_b17, &work[wpos],
		     &ldw, (ftnlen)5, (ftnlen)1);

#line 408 "ssytrd_sy2sb.f"
	    sgemm_("No transpose", "Conjugate", &pk, &pk, &pn, &c_b23, &work[
		    wpos], &ldw, &work[s2pos], &lds2, &c_b17, &work[s1pos], &
		    lds1, (ftnlen)12, (ftnlen)9);

#line 413 "ssytrd_sy2sb.f"
	    sgemm_("No transpose", "No transpose", &pk, &pn, &pk, &c_b39, &
		    work[s1pos], &lds1, &a[i__ + (i__ + *kd) * a_dim1], lda, &
		    c_b23, &work[wpos], &ldw, (ftnlen)12, (ftnlen)12);


/*            Update the unreduced submatrix A(i+kd:n,i+kd:n), using */
/*            an update of the form:  A := A - V'*W - W'*V */

#line 422 "ssytrd_sy2sb.f"
	    ssyr2k_(uplo, "Conjugate", &pn, &pk, &c_b42, &a[i__ + (i__ + *kd) 
		    * a_dim1], lda, &work[wpos], &ldw, &c_b23, &a[i__ + *kd + 
		    (i__ + *kd) * a_dim1], lda, (ftnlen)1, (ftnlen)9);
#line 426 "ssytrd_sy2sb.f"
/* L10: */
#line 426 "ssytrd_sy2sb.f"
	}

/*        Copy the upper band to AB which is the band storage matrix */

#line 430 "ssytrd_sy2sb.f"
	i__2 = *n;
#line 430 "ssytrd_sy2sb.f"
	for (j = *n - *kd + 1; j <= i__2; ++j) {
/* Computing MIN */
#line 431 "ssytrd_sy2sb.f"
	    i__1 = *kd, i__3 = *n - j;
#line 431 "ssytrd_sy2sb.f"
	    lk = min(i__1,i__3) + 1;
#line 432 "ssytrd_sy2sb.f"
	    i__1 = *ldab - 1;
#line 432 "ssytrd_sy2sb.f"
	    scopy_(&lk, &a[j + j * a_dim1], lda, &ab[*kd + 1 + j * ab_dim1], &
		    i__1);
#line 433 "ssytrd_sy2sb.f"
/* L30: */
#line 433 "ssytrd_sy2sb.f"
	}

#line 435 "ssytrd_sy2sb.f"
    } else {

/*         Reduce the lower triangle of A to lower band matrix */

#line 439 "ssytrd_sy2sb.f"
	i__2 = *n - *kd;
#line 439 "ssytrd_sy2sb.f"
	i__1 = *kd;
#line 439 "ssytrd_sy2sb.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 440 "ssytrd_sy2sb.f"
	    pn = *n - i__ - *kd + 1;
/* Computing MIN */
#line 441 "ssytrd_sy2sb.f"
	    i__3 = *n - i__ - *kd + 1;
#line 441 "ssytrd_sy2sb.f"
	    pk = min(i__3,*kd);

/*            Compute the QR factorization of the current block */

#line 445 "ssytrd_sy2sb.f"
	    sgeqrf_(&pn, kd, &a[i__ + *kd + i__ * a_dim1], lda, &tau[i__], &
		    work[s2pos], &ls2, &iinfo);

/*            Copy the upper portion of A into AB */

#line 450 "ssytrd_sy2sb.f"
	    i__3 = i__ + pk - 1;
#line 450 "ssytrd_sy2sb.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MIN */
#line 451 "ssytrd_sy2sb.f"
		i__4 = *kd, i__5 = *n - j;
#line 451 "ssytrd_sy2sb.f"
		lk = min(i__4,i__5) + 1;
#line 452 "ssytrd_sy2sb.f"
		scopy_(&lk, &a[j + j * a_dim1], &c__1, &ab[j * ab_dim1 + 1], &
			c__1);
#line 453 "ssytrd_sy2sb.f"
/* L50: */
#line 453 "ssytrd_sy2sb.f"
	    }

#line 455 "ssytrd_sy2sb.f"
	    slaset_("Upper", &pk, &pk, &c_b17, &c_b23, &a[i__ + *kd + i__ * 
		    a_dim1], lda, (ftnlen)5);

/*            Form the matrix T */

#line 460 "ssytrd_sy2sb.f"
	    slarft_("Forward", "Columnwise", &pn, &pk, &a[i__ + *kd + i__ * 
		    a_dim1], lda, &tau[i__], &work[tpos], &ldt, (ftnlen)7, (
		    ftnlen)10);

/*            Compute W: */

#line 466 "ssytrd_sy2sb.f"
	    sgemm_("No transpose", "No transpose", &pn, &pk, &pk, &c_b23, &a[
		    i__ + *kd + i__ * a_dim1], lda, &work[tpos], &ldt, &c_b17,
		     &work[s2pos], &lds2, (ftnlen)12, (ftnlen)12);

#line 471 "ssytrd_sy2sb.f"
	    ssymm_("Left", uplo, &pn, &pk, &c_b23, &a[i__ + *kd + (i__ + *kd) 
		    * a_dim1], lda, &work[s2pos], &lds2, &c_b17, &work[wpos], 
		    &ldw, (ftnlen)4, (ftnlen)1);

#line 476 "ssytrd_sy2sb.f"
	    sgemm_("Conjugate", "No transpose", &pk, &pk, &pn, &c_b23, &work[
		    s2pos], &lds2, &work[wpos], &ldw, &c_b17, &work[s1pos], &
		    lds1, (ftnlen)9, (ftnlen)12);

#line 481 "ssytrd_sy2sb.f"
	    sgemm_("No transpose", "No transpose", &pn, &pk, &pk, &c_b39, &a[
		    i__ + *kd + i__ * a_dim1], lda, &work[s1pos], &lds1, &
		    c_b23, &work[wpos], &ldw, (ftnlen)12, (ftnlen)12);


/*            Update the unreduced submatrix A(i+kd:n,i+kd:n), using */
/*            an update of the form:  A := A - V*W' - W*V' */

#line 490 "ssytrd_sy2sb.f"
	    ssyr2k_(uplo, "No transpose", &pn, &pk, &c_b42, &a[i__ + *kd + 
		    i__ * a_dim1], lda, &work[wpos], &ldw, &c_b23, &a[i__ + *
		    kd + (i__ + *kd) * a_dim1], lda, (ftnlen)1, (ftnlen)12);
/*            ================================================================== */
/*            RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED */
/*             DO 45 J = I, I+PK-1 */
/*                LK = MIN( KD, N-J ) + 1 */
/*                CALL SCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 ) */
/*   45        CONTINUE */
/*            ================================================================== */
#line 501 "ssytrd_sy2sb.f"
/* L40: */
#line 501 "ssytrd_sy2sb.f"
	}

/*        Copy the lower band to AB which is the band storage matrix */

#line 505 "ssytrd_sy2sb.f"
	i__1 = *n;
#line 505 "ssytrd_sy2sb.f"
	for (j = *n - *kd + 1; j <= i__1; ++j) {
/* Computing MIN */
#line 506 "ssytrd_sy2sb.f"
	    i__2 = *kd, i__3 = *n - j;
#line 506 "ssytrd_sy2sb.f"
	    lk = min(i__2,i__3) + 1;
#line 507 "ssytrd_sy2sb.f"
	    scopy_(&lk, &a[j + j * a_dim1], &c__1, &ab[j * ab_dim1 + 1], &
		    c__1);
#line 508 "ssytrd_sy2sb.f"
/* L60: */
#line 508 "ssytrd_sy2sb.f"
	}
#line 510 "ssytrd_sy2sb.f"
    }

#line 512 "ssytrd_sy2sb.f"
    work[1] = (doublereal) lwmin;
#line 513 "ssytrd_sy2sb.f"
    return 0;

/*     End of SSYTRD_SY2SB */

} /* ssytrd_sy2sb__ */

