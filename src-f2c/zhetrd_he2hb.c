#line 1 "zhetrd_he2hb.f"
/* zhetrd_he2hb.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrd_he2hb.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__20 = 20;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b33 = 1.;

/* > \brief \b ZHETRD_HE2HB */

/*  @precisions fortran z -> s d c */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRD_HE2HB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRD_HE2HB( UPLO, N, KD, A, LDA, AB, LDAB, TAU, */
/*                              WORK, LWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAB, LWORK, N, KD */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), AB( LDAB, * ), */
/*                          TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRD_HE2HB reduces a complex Hermitian matrix A to complex Hermitian */
/* > band-diagonal form AB by a unitary similarity transformation: */
/* > Q**H * A * Q = AB. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* >          of A are overwritten by the corresponding elements of the */
/* >          tridiagonal matrix T, and the elements above the first */
/* >          superdiagonal, with the array TAU, represent the unitary */
/* >          matrix Q as a product of elementary reflectors; if UPLO */
/* >          = 'L', the diagonal and first subdiagonal of A are over- */
/* >          written by the corresponding elements of the tridiagonal */
/* >          matrix T, and the elements below the first subdiagonal, with */
/* >          the array TAU, represent the unitary matrix Q as a product */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          On exit, the upper or lower triangle of the Hermitian band */
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
/* >          TAU is COMPLEX*16 array, dimension (N-KD) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (LWORK) */
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

/* > \ingroup complex16HEcomputational */

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
/* >     Q = H(k)**H . . . H(2)**H H(1)**H, where k = n-kd. */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
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
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
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
/* Subroutine */ int zhetrd_he2hb__(char *uplo, integer *n, integer *kd, 
	doublecomplex *a, integer *lda, doublecomplex *ab, integer *ldab, 
	doublecomplex *tau, doublecomplex *work, integer *lwork, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, ab_dim1, ab_offset, i__1, i__2, i__3, i__4, 
	    i__5;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, lk, pk, pn, lt, lw, ls1, ls2, ldt, ldw, lds1, lds2,
	     tpos, wpos, s1pos, s2pos;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zhemm_(char *, char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer lwmin;
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zher2k_(char *, char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zgeqrf_(integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), zlarft_(
	    char *, char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
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

#line 296 "zhetrd_he2hb.f"
    /* Parameter adjustments */
#line 296 "zhetrd_he2hb.f"
    a_dim1 = *lda;
#line 296 "zhetrd_he2hb.f"
    a_offset = 1 + a_dim1;
#line 296 "zhetrd_he2hb.f"
    a -= a_offset;
#line 296 "zhetrd_he2hb.f"
    ab_dim1 = *ldab;
#line 296 "zhetrd_he2hb.f"
    ab_offset = 1 + ab_dim1;
#line 296 "zhetrd_he2hb.f"
    ab -= ab_offset;
#line 296 "zhetrd_he2hb.f"
    --tau;
#line 296 "zhetrd_he2hb.f"
    --work;
#line 296 "zhetrd_he2hb.f"

#line 296 "zhetrd_he2hb.f"
    /* Function Body */
#line 296 "zhetrd_he2hb.f"
    *info = 0;
#line 297 "zhetrd_he2hb.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 298 "zhetrd_he2hb.f"
    lquery = *lwork == -1;
#line 299 "zhetrd_he2hb.f"
    lwmin = ilaenv_(&c__20, "ZHETRD_HE2HB", "", n, kd, &c_n1, &c_n1, (ftnlen)
	    12, (ftnlen)0);
#line 301 "zhetrd_he2hb.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 302 "zhetrd_he2hb.f"
	*info = -1;
#line 303 "zhetrd_he2hb.f"
    } else if (*n < 0) {
#line 304 "zhetrd_he2hb.f"
	*info = -2;
#line 305 "zhetrd_he2hb.f"
    } else if (*kd < 0) {
#line 306 "zhetrd_he2hb.f"
	*info = -3;
#line 307 "zhetrd_he2hb.f"
    } else if (*lda < max(1,*n)) {
#line 308 "zhetrd_he2hb.f"
	*info = -5;
#line 309 "zhetrd_he2hb.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 309 "zhetrd_he2hb.f"
	i__1 = 1, i__2 = *kd + 1;
#line 309 "zhetrd_he2hb.f"
	if (*ldab < max(i__1,i__2)) {
#line 310 "zhetrd_he2hb.f"
	    *info = -7;
#line 311 "zhetrd_he2hb.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 312 "zhetrd_he2hb.f"
	    *info = -10;
#line 313 "zhetrd_he2hb.f"
	}
#line 313 "zhetrd_he2hb.f"
    }

#line 315 "zhetrd_he2hb.f"
    if (*info != 0) {
#line 316 "zhetrd_he2hb.f"
	i__1 = -(*info);
#line 316 "zhetrd_he2hb.f"
	xerbla_("ZHETRD_HE2HB", &i__1, (ftnlen)12);
#line 317 "zhetrd_he2hb.f"
	return 0;
#line 318 "zhetrd_he2hb.f"
    } else if (lquery) {
#line 319 "zhetrd_he2hb.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 320 "zhetrd_he2hb.f"
	return 0;
#line 321 "zhetrd_he2hb.f"
    }

/*     Quick return if possible */
/*     Copy the upper/lower portion of A into AB */

#line 326 "zhetrd_he2hb.f"
    if (*n <= *kd + 1) {
#line 327 "zhetrd_he2hb.f"
	if (upper) {
#line 328 "zhetrd_he2hb.f"
	    i__1 = *n;
#line 328 "zhetrd_he2hb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 329 "zhetrd_he2hb.f"
		i__2 = *kd + 1;
#line 329 "zhetrd_he2hb.f"
		lk = min(i__2,i__);
#line 330 "zhetrd_he2hb.f"
		zcopy_(&lk, &a[i__ - lk + 1 + i__ * a_dim1], &c__1, &ab[*kd + 
			1 - lk + 1 + i__ * ab_dim1], &c__1);
#line 332 "zhetrd_he2hb.f"
/* L100: */
#line 332 "zhetrd_he2hb.f"
	    }
#line 333 "zhetrd_he2hb.f"
	} else {
#line 334 "zhetrd_he2hb.f"
	    i__1 = *n;
#line 334 "zhetrd_he2hb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 335 "zhetrd_he2hb.f"
		i__2 = *kd + 1, i__3 = *n - i__ + 1;
#line 335 "zhetrd_he2hb.f"
		lk = min(i__2,i__3);
#line 336 "zhetrd_he2hb.f"
		zcopy_(&lk, &a[i__ + i__ * a_dim1], &c__1, &ab[i__ * ab_dim1 
			+ 1], &c__1);
#line 337 "zhetrd_he2hb.f"
/* L110: */
#line 337 "zhetrd_he2hb.f"
	    }
#line 338 "zhetrd_he2hb.f"
	}
#line 339 "zhetrd_he2hb.f"
	work[1].r = 1., work[1].i = 0.;
#line 340 "zhetrd_he2hb.f"
	return 0;
#line 341 "zhetrd_he2hb.f"
    }

/*     Determine the pointer position for the workspace */

#line 345 "zhetrd_he2hb.f"
    ldt = *kd;
#line 346 "zhetrd_he2hb.f"
    lds1 = *kd;
#line 347 "zhetrd_he2hb.f"
    lt = ldt * *kd;
#line 348 "zhetrd_he2hb.f"
    lw = *n * *kd;
#line 349 "zhetrd_he2hb.f"
    ls1 = lds1 * *kd;
#line 350 "zhetrd_he2hb.f"
    ls2 = lwmin - lt - lw - ls1;
/*      LS2 = N*MAX(KD,FACTOPTNB) */
#line 352 "zhetrd_he2hb.f"
    tpos = 1;
#line 353 "zhetrd_he2hb.f"
    wpos = tpos + lt;
#line 354 "zhetrd_he2hb.f"
    s1pos = wpos + lw;
#line 355 "zhetrd_he2hb.f"
    s2pos = s1pos + ls1;
#line 356 "zhetrd_he2hb.f"
    if (upper) {
#line 357 "zhetrd_he2hb.f"
	ldw = *kd;
#line 358 "zhetrd_he2hb.f"
	lds2 = *kd;
#line 359 "zhetrd_he2hb.f"
    } else {
#line 360 "zhetrd_he2hb.f"
	ldw = *n;
#line 361 "zhetrd_he2hb.f"
	lds2 = *n;
#line 362 "zhetrd_he2hb.f"
    }


/*     Set the workspace of the triangular matrix T to zero once such a */
/*     way everytime T is generated the upper/lower portion will be always zero */

#line 368 "zhetrd_he2hb.f"
    zlaset_("A", &ldt, kd, &c_b1, &c_b1, &work[tpos], &ldt, (ftnlen)1);

#line 370 "zhetrd_he2hb.f"
    if (upper) {
#line 371 "zhetrd_he2hb.f"
	i__1 = *n - *kd;
#line 371 "zhetrd_he2hb.f"
	i__2 = *kd;
#line 371 "zhetrd_he2hb.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 372 "zhetrd_he2hb.f"
	    pn = *n - i__ - *kd + 1;
/* Computing MIN */
#line 373 "zhetrd_he2hb.f"
	    i__3 = *n - i__ - *kd + 1;
#line 373 "zhetrd_he2hb.f"
	    pk = min(i__3,*kd);

/*            Compute the LQ factorization of the current block */

#line 377 "zhetrd_he2hb.f"
	    zgelqf_(kd, &pn, &a[i__ + (i__ + *kd) * a_dim1], lda, &tau[i__], &
		    work[s2pos], &ls2, &iinfo);

/*            Copy the upper portion of A into AB */

#line 382 "zhetrd_he2hb.f"
	    i__3 = i__ + pk - 1;
#line 382 "zhetrd_he2hb.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MIN */
#line 383 "zhetrd_he2hb.f"
		i__4 = *kd, i__5 = *n - j;
#line 383 "zhetrd_he2hb.f"
		lk = min(i__4,i__5) + 1;
#line 384 "zhetrd_he2hb.f"
		i__4 = *ldab - 1;
#line 384 "zhetrd_he2hb.f"
		zcopy_(&lk, &a[j + j * a_dim1], lda, &ab[*kd + 1 + j * 
			ab_dim1], &i__4);
#line 385 "zhetrd_he2hb.f"
/* L20: */
#line 385 "zhetrd_he2hb.f"
	    }

#line 387 "zhetrd_he2hb.f"
	    zlaset_("Lower", &pk, &pk, &c_b1, &c_b2, &a[i__ + (i__ + *kd) * 
		    a_dim1], lda, (ftnlen)5);

/*            Form the matrix T */

#line 392 "zhetrd_he2hb.f"
	    zlarft_("Forward", "Rowwise", &pn, &pk, &a[i__ + (i__ + *kd) * 
		    a_dim1], lda, &tau[i__], &work[tpos], &ldt, (ftnlen)7, (
		    ftnlen)7);

/*            Compute W: */

#line 398 "zhetrd_he2hb.f"
	    zgemm_("Conjugate", "No transpose", &pk, &pn, &pk, &c_b2, &work[
		    tpos], &ldt, &a[i__ + (i__ + *kd) * a_dim1], lda, &c_b1, &
		    work[s2pos], &lds2, (ftnlen)9, (ftnlen)12);

#line 403 "zhetrd_he2hb.f"
	    zhemm_("Right", uplo, &pk, &pn, &c_b2, &a[i__ + *kd + (i__ + *kd) 
		    * a_dim1], lda, &work[s2pos], &lds2, &c_b1, &work[wpos], &
		    ldw, (ftnlen)5, (ftnlen)1);

#line 408 "zhetrd_he2hb.f"
	    zgemm_("No transpose", "Conjugate", &pk, &pk, &pn, &c_b2, &work[
		    wpos], &ldw, &work[s2pos], &lds2, &c_b1, &work[s1pos], &
		    lds1, (ftnlen)12, (ftnlen)9);

#line 413 "zhetrd_he2hb.f"
	    z__1.r = -.5, z__1.i = -0.;
#line 413 "zhetrd_he2hb.f"
	    zgemm_("No transpose", "No transpose", &pk, &pn, &pk, &z__1, &
		    work[s1pos], &lds1, &a[i__ + (i__ + *kd) * a_dim1], lda, &
		    c_b2, &work[wpos], &ldw, (ftnlen)12, (ftnlen)12);


/*            Update the unreduced submatrix A(i+kd:n,i+kd:n), using */
/*            an update of the form:  A := A - V'*W - W'*V */

#line 422 "zhetrd_he2hb.f"
	    z__1.r = -1., z__1.i = -0.;
#line 422 "zhetrd_he2hb.f"
	    zher2k_(uplo, "Conjugate", &pn, &pk, &z__1, &a[i__ + (i__ + *kd) *
		     a_dim1], lda, &work[wpos], &ldw, &c_b33, &a[i__ + *kd + (
		    i__ + *kd) * a_dim1], lda, (ftnlen)1, (ftnlen)9);
#line 426 "zhetrd_he2hb.f"
/* L10: */
#line 426 "zhetrd_he2hb.f"
	}

/*        Copy the upper band to AB which is the band storage matrix */

#line 430 "zhetrd_he2hb.f"
	i__2 = *n;
#line 430 "zhetrd_he2hb.f"
	for (j = *n - *kd + 1; j <= i__2; ++j) {
/* Computing MIN */
#line 431 "zhetrd_he2hb.f"
	    i__1 = *kd, i__3 = *n - j;
#line 431 "zhetrd_he2hb.f"
	    lk = min(i__1,i__3) + 1;
#line 432 "zhetrd_he2hb.f"
	    i__1 = *ldab - 1;
#line 432 "zhetrd_he2hb.f"
	    zcopy_(&lk, &a[j + j * a_dim1], lda, &ab[*kd + 1 + j * ab_dim1], &
		    i__1);
#line 433 "zhetrd_he2hb.f"
/* L30: */
#line 433 "zhetrd_he2hb.f"
	}

#line 435 "zhetrd_he2hb.f"
    } else {

/*         Reduce the lower triangle of A to lower band matrix */

#line 439 "zhetrd_he2hb.f"
	i__2 = *n - *kd;
#line 439 "zhetrd_he2hb.f"
	i__1 = *kd;
#line 439 "zhetrd_he2hb.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 440 "zhetrd_he2hb.f"
	    pn = *n - i__ - *kd + 1;
/* Computing MIN */
#line 441 "zhetrd_he2hb.f"
	    i__3 = *n - i__ - *kd + 1;
#line 441 "zhetrd_he2hb.f"
	    pk = min(i__3,*kd);

/*            Compute the QR factorization of the current block */

#line 445 "zhetrd_he2hb.f"
	    zgeqrf_(&pn, kd, &a[i__ + *kd + i__ * a_dim1], lda, &tau[i__], &
		    work[s2pos], &ls2, &iinfo);

/*            Copy the upper portion of A into AB */

#line 450 "zhetrd_he2hb.f"
	    i__3 = i__ + pk - 1;
#line 450 "zhetrd_he2hb.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MIN */
#line 451 "zhetrd_he2hb.f"
		i__4 = *kd, i__5 = *n - j;
#line 451 "zhetrd_he2hb.f"
		lk = min(i__4,i__5) + 1;
#line 452 "zhetrd_he2hb.f"
		zcopy_(&lk, &a[j + j * a_dim1], &c__1, &ab[j * ab_dim1 + 1], &
			c__1);
#line 453 "zhetrd_he2hb.f"
/* L50: */
#line 453 "zhetrd_he2hb.f"
	    }

#line 455 "zhetrd_he2hb.f"
	    zlaset_("Upper", &pk, &pk, &c_b1, &c_b2, &a[i__ + *kd + i__ * 
		    a_dim1], lda, (ftnlen)5);

/*            Form the matrix T */

#line 460 "zhetrd_he2hb.f"
	    zlarft_("Forward", "Columnwise", &pn, &pk, &a[i__ + *kd + i__ * 
		    a_dim1], lda, &tau[i__], &work[tpos], &ldt, (ftnlen)7, (
		    ftnlen)10);

/*            Compute W: */

#line 466 "zhetrd_he2hb.f"
	    zgemm_("No transpose", "No transpose", &pn, &pk, &pk, &c_b2, &a[
		    i__ + *kd + i__ * a_dim1], lda, &work[tpos], &ldt, &c_b1, 
		    &work[s2pos], &lds2, (ftnlen)12, (ftnlen)12);

#line 471 "zhetrd_he2hb.f"
	    zhemm_("Left", uplo, &pn, &pk, &c_b2, &a[i__ + *kd + (i__ + *kd) *
		     a_dim1], lda, &work[s2pos], &lds2, &c_b1, &work[wpos], &
		    ldw, (ftnlen)4, (ftnlen)1);

#line 476 "zhetrd_he2hb.f"
	    zgemm_("Conjugate", "No transpose", &pk, &pk, &pn, &c_b2, &work[
		    s2pos], &lds2, &work[wpos], &ldw, &c_b1, &work[s1pos], &
		    lds1, (ftnlen)9, (ftnlen)12);

#line 481 "zhetrd_he2hb.f"
	    z__1.r = -.5, z__1.i = -0.;
#line 481 "zhetrd_he2hb.f"
	    zgemm_("No transpose", "No transpose", &pn, &pk, &pk, &z__1, &a[
		    i__ + *kd + i__ * a_dim1], lda, &work[s1pos], &lds1, &
		    c_b2, &work[wpos], &ldw, (ftnlen)12, (ftnlen)12);


/*            Update the unreduced submatrix A(i+kd:n,i+kd:n), using */
/*            an update of the form:  A := A - V*W' - W*V' */

#line 490 "zhetrd_he2hb.f"
	    z__1.r = -1., z__1.i = -0.;
#line 490 "zhetrd_he2hb.f"
	    zher2k_(uplo, "No transpose", &pn, &pk, &z__1, &a[i__ + *kd + i__ 
		    * a_dim1], lda, &work[wpos], &ldw, &c_b33, &a[i__ + *kd + 
		    (i__ + *kd) * a_dim1], lda, (ftnlen)1, (ftnlen)12);
/*            ================================================================== */
/*            RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED */
/*             DO 45 J = I, I+PK-1 */
/*                LK = MIN( KD, N-J ) + 1 */
/*                CALL ZCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 ) */
/*   45        CONTINUE */
/*            ================================================================== */
#line 501 "zhetrd_he2hb.f"
/* L40: */
#line 501 "zhetrd_he2hb.f"
	}

/*        Copy the lower band to AB which is the band storage matrix */

#line 505 "zhetrd_he2hb.f"
	i__1 = *n;
#line 505 "zhetrd_he2hb.f"
	for (j = *n - *kd + 1; j <= i__1; ++j) {
/* Computing MIN */
#line 506 "zhetrd_he2hb.f"
	    i__2 = *kd, i__3 = *n - j;
#line 506 "zhetrd_he2hb.f"
	    lk = min(i__2,i__3) + 1;
#line 507 "zhetrd_he2hb.f"
	    zcopy_(&lk, &a[j + j * a_dim1], &c__1, &ab[j * ab_dim1 + 1], &
		    c__1);
#line 508 "zhetrd_he2hb.f"
/* L60: */
#line 508 "zhetrd_he2hb.f"
	}
#line 510 "zhetrd_he2hb.f"
    }

#line 512 "zhetrd_he2hb.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 513 "zhetrd_he2hb.f"
    return 0;

/*     End of ZHETRD_HE2HB */

} /* zhetrd_he2hb__ */

