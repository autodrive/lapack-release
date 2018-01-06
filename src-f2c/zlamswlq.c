#line 1 "zlamswlq.f"
/* zlamswlq.f -- translated by f2c (version 20100827).
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

#line 1 "zlamswlq.f"
/* Table of constant values */

static integer c__0 = 0;


/*  Definition: */
/*  =========== */

/*      SUBROUTINE ZLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
/*     $                LDT, C, LDC, WORK, LWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER         SIDE, TRANS */
/*      INTEGER           INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC */
/*     .. */
/*     .. Array Arguments .. */
/*      COMPLEX*16        A( LDA, * ), WORK( * ), C(LDC, * ), */
/*     $                  T( LDT, * ) */
/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLAMQRTS overwrites the general real M-by-N matrix C with */
/* > */
/* > */
/* >                    SIDE = 'L'     SIDE = 'R' */
/* >    TRANS = 'N':      Q * C          C * Q */
/* >    TRANS = 'T':      Q**T * C       C * Q**T */
/* >    where Q is a real orthogonal matrix defined as the product of blocked */
/* >    elementary reflectors computed by short wide LQ */
/* >    factorization (ZLASWLQ) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >=0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          M >= K >= 0; */
/* > */
/* > \endverbatim */
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The row block size to be used in the blocked QR. */
/* >          M >= MB >= 1 */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The column block size to be used in the blocked QR. */
/* >          NB > M. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR. */
/* >                MB > M. */
/* > */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,K) */
/* >          The i-th row must contain the vector which defines the blocked */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DLASWLQ in the first k rows of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension */
/* >          ( M * Number of blocks(CEIL(N-K/NB-K)), */
/* >          The blocked upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= MB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >         (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If SIDE = 'L', LWORK >= max(1,NB) * MB; */
/* >          if SIDE = 'R', LWORK >= max(1,M) * MB. */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > Short-Wide LQ (SWLQ) performs LQ by a sequence of orthogonal transformations, */
/* > representing Q as a product of other orthogonal matrices */
/* >   Q = Q(1) * Q(2) * . . . * Q(k) */
/* > where each Q(i) zeros out upper diagonal entries of a block of NB rows of A: */
/* >   Q(1) zeros out the upper diagonal entries of rows 1:NB of A */
/* >   Q(2) zeros out the bottom MB-N rows of rows [1:M,NB+1:2*NB-M] of A */
/* >   Q(3) zeros out the bottom MB-N rows of rows [1:M,2*NB-M+1:3*NB-2*M] of A */
/* >   . . . */
/* > */
/* > Q(1) is computed by GELQT, which represents Q(1) by Householder vectors */
/* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,1:N). */
/* > For more information see Further Details in GELQT. */
/* > */
/* > Q(i) for i>1 is computed by TPLQT, which represents Q(i) by Householder vectors */
/* > stored in columns [(i-1)*(NB-M)+M+1:i*(NB-M)+M] of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,(i-1)*M+1:i*M). */
/* > The last Q(k) may use fewer rows. */
/* > For more information see Further Details in TPQRT. */
/* > */
/* > For more details of the overall algorithm, see the description of */
/* > Sequential TSQR in Section 2.2 of [1]. */
/* > */
/* > [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations, */
/* >     J. Demmel, L. Grigori, M. Hoemmen, J. Langou, */
/* >     SIAM J. Sci. Comput, vol. 34, no. 1, 2012 */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlamswlq_(char *side, char *trans, integer *m, integer *
	n, integer *k, integer *mb, integer *nb, doublecomplex *a, integer *
	lda, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *ldc,
	 doublecomplex *work, integer *lwork, integer *info, ftnlen side_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, ii, kk, lw, ctr;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran, lquery;
    extern /* Subroutine */ int zgemlqt_(char *, char *, integer *, integer *,
	     integer *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, doublecomplex *, integer 
	    *, ftnlen, ftnlen), ztpmlqt_(char *, char *, integer *, integer *,
	     integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 234 "zlamswlq.f"
    /* Parameter adjustments */
#line 234 "zlamswlq.f"
    a_dim1 = *lda;
#line 234 "zlamswlq.f"
    a_offset = 1 + a_dim1;
#line 234 "zlamswlq.f"
    a -= a_offset;
#line 234 "zlamswlq.f"
    t_dim1 = *ldt;
#line 234 "zlamswlq.f"
    t_offset = 1 + t_dim1;
#line 234 "zlamswlq.f"
    t -= t_offset;
#line 234 "zlamswlq.f"
    c_dim1 = *ldc;
#line 234 "zlamswlq.f"
    c_offset = 1 + c_dim1;
#line 234 "zlamswlq.f"
    c__ -= c_offset;
#line 234 "zlamswlq.f"
    --work;
#line 234 "zlamswlq.f"

#line 234 "zlamswlq.f"
    /* Function Body */
#line 234 "zlamswlq.f"
    lquery = *lwork < 0;
#line 235 "zlamswlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 236 "zlamswlq.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 237 "zlamswlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 238 "zlamswlq.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 239 "zlamswlq.f"
    if (left) {
#line 240 "zlamswlq.f"
	lw = *n * *mb;
#line 241 "zlamswlq.f"
    } else {
#line 242 "zlamswlq.f"
	lw = *m * *mb;
#line 243 "zlamswlq.f"
    }

#line 245 "zlamswlq.f"
    *info = 0;
#line 246 "zlamswlq.f"
    if (! left && ! right) {
#line 247 "zlamswlq.f"
	*info = -1;
#line 248 "zlamswlq.f"
    } else if (! tran && ! notran) {
#line 249 "zlamswlq.f"
	*info = -2;
#line 250 "zlamswlq.f"
    } else if (*m < 0) {
#line 251 "zlamswlq.f"
	*info = -3;
#line 252 "zlamswlq.f"
    } else if (*n < 0) {
#line 253 "zlamswlq.f"
	*info = -4;
#line 254 "zlamswlq.f"
    } else if (*k < 0) {
#line 255 "zlamswlq.f"
	*info = -5;
#line 256 "zlamswlq.f"
    } else if (*lda < max(1,*k)) {
#line 257 "zlamswlq.f"
	*info = -9;
#line 258 "zlamswlq.f"
    } else if (*ldt < max(1,*mb)) {
#line 259 "zlamswlq.f"
	*info = -11;
#line 260 "zlamswlq.f"
    } else if (*ldc < max(1,*m)) {
#line 261 "zlamswlq.f"
	*info = -13;
#line 262 "zlamswlq.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 263 "zlamswlq.f"
	*info = -15;
#line 264 "zlamswlq.f"
    }

#line 266 "zlamswlq.f"
    if (*info != 0) {
#line 267 "zlamswlq.f"
	i__1 = -(*info);
#line 267 "zlamswlq.f"
	xerbla_("ZLAMSWLQ", &i__1, (ftnlen)8);
#line 268 "zlamswlq.f"
	work[1].r = (doublereal) lw, work[1].i = 0.;
#line 269 "zlamswlq.f"
	return 0;
#line 270 "zlamswlq.f"
    } else if (lquery) {
#line 271 "zlamswlq.f"
	work[1].r = (doublereal) lw, work[1].i = 0.;
#line 272 "zlamswlq.f"
	return 0;
#line 273 "zlamswlq.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 277 "zlamswlq.f"
    i__1 = min(*m,*n);
#line 277 "zlamswlq.f"
    if (min(i__1,*k) == 0) {
#line 278 "zlamswlq.f"
	return 0;
#line 279 "zlamswlq.f"
    }

/* Computing MAX */
#line 281 "zlamswlq.f"
    i__1 = max(*m,*n);
#line 281 "zlamswlq.f"
    if (*nb <= *k || *nb >= max(i__1,*k)) {
#line 282 "zlamswlq.f"
	zgemlqt_(side, trans, m, n, k, mb, &a[a_offset], lda, &t[t_offset], 
		ldt, &c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)
		1);
#line 284 "zlamswlq.f"
	return 0;
#line 285 "zlamswlq.f"
    }

#line 287 "zlamswlq.f"
    if (left && tran) {

/*         Multiply Q to the last block of C */

#line 291 "zlamswlq.f"
	kk = (*m - *k) % (*nb - *k);
#line 292 "zlamswlq.f"
	ctr = (*m - *k) / (*nb - *k);

#line 294 "zlamswlq.f"
	if (kk > 0) {
#line 295 "zlamswlq.f"
	    ii = *m - kk + 1;
#line 296 "zlamswlq.f"
	    ztpmlqt_("L", "C", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 299 "zlamswlq.f"
	} else {
#line 300 "zlamswlq.f"
	    ii = *m + 1;
#line 301 "zlamswlq.f"
	}

#line 303 "zlamswlq.f"
	i__1 = *nb + 1;
#line 303 "zlamswlq.f"
	i__2 = -(*nb - *k);
#line 303 "zlamswlq.f"
	for (i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+NB) */

#line 307 "zlamswlq.f"
	    --ctr;
#line 308 "zlamswlq.f"
	    i__3 = *nb - *k;
#line 308 "zlamswlq.f"
	    ztpmlqt_("L", "C", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 312 "zlamswlq.f"
	}

/*         Multiply Q to the first block of C (1:M,1:NB) */

#line 316 "zlamswlq.f"
	zgemlqt_("L", "C", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 319 "zlamswlq.f"
    } else if (left && notran) {

/*         Multiply Q to the first block of C */

#line 323 "zlamswlq.f"
	kk = (*m - *k) % (*nb - *k);
#line 324 "zlamswlq.f"
	ii = *m - kk + 1;
#line 325 "zlamswlq.f"
	ctr = 1;
#line 326 "zlamswlq.f"
	zgemlqt_("L", "N", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 329 "zlamswlq.f"
	i__2 = ii - *nb + *k;
#line 329 "zlamswlq.f"
	i__1 = *nb - *k;
#line 329 "zlamswlq.f"
	for (i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (I:I+NB,1:N) */

#line 333 "zlamswlq.f"
	    i__3 = *nb - *k;
#line 333 "zlamswlq.f"
	    ztpmlqt_("L", "N", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 336 "zlamswlq.f"
	    ++ctr;

#line 338 "zlamswlq.f"
	}
#line 339 "zlamswlq.f"
	if (ii <= *m) {

/*         Multiply Q to the last block of C */

#line 343 "zlamswlq.f"
	    ztpmlqt_("L", "N", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 347 "zlamswlq.f"
	}

#line 349 "zlamswlq.f"
    } else if (right && notran) {

/*         Multiply Q to the last block of C */

#line 353 "zlamswlq.f"
	kk = (*n - *k) % (*nb - *k);
#line 354 "zlamswlq.f"
	ctr = (*n - *k) / (*nb - *k);
#line 355 "zlamswlq.f"
	if (kk > 0) {
#line 356 "zlamswlq.f"
	    ii = *n - kk + 1;
#line 357 "zlamswlq.f"
	    ztpmlqt_("R", "N", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 360 "zlamswlq.f"
	} else {
#line 361 "zlamswlq.f"
	    ii = *n + 1;
#line 362 "zlamswlq.f"
	}

#line 364 "zlamswlq.f"
	i__1 = *nb + 1;
#line 364 "zlamswlq.f"
	i__2 = -(*nb - *k);
#line 364 "zlamswlq.f"
	for (i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 368 "zlamswlq.f"
	    --ctr;
#line 369 "zlamswlq.f"
	    i__3 = *nb - *k;
#line 369 "zlamswlq.f"
	    ztpmlqt_("R", "N", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (
		    ftnlen)1, (ftnlen)1);
#line 373 "zlamswlq.f"
	}

/*         Multiply Q to the first block of C (1:M,1:MB) */

#line 377 "zlamswlq.f"
	zgemlqt_("R", "N", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 380 "zlamswlq.f"
    } else if (right && tran) {

/*       Multiply Q to the first block of C */

#line 384 "zlamswlq.f"
	kk = (*n - *k) % (*nb - *k);
#line 385 "zlamswlq.f"
	ii = *n - kk + 1;
#line 386 "zlamswlq.f"
	zgemlqt_("R", "C", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);
#line 388 "zlamswlq.f"
	ctr = 1;

#line 390 "zlamswlq.f"
	i__2 = ii - *nb + *k;
#line 390 "zlamswlq.f"
	i__1 = *nb - *k;
#line 390 "zlamswlq.f"
	for (i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 394 "zlamswlq.f"
	    i__3 = *nb - *k;
#line 394 "zlamswlq.f"
	    ztpmlqt_("R", "C", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (
		    ftnlen)1, (ftnlen)1);
#line 397 "zlamswlq.f"
	    ++ctr;

#line 399 "zlamswlq.f"
	}
#line 400 "zlamswlq.f"
	if (ii <= *n) {

/*       Multiply Q to the last block of C */

#line 404 "zlamswlq.f"
	    ztpmlqt_("R", "C", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);

#line 408 "zlamswlq.f"
	}

#line 410 "zlamswlq.f"
    }

#line 412 "zlamswlq.f"
    work[1].r = (doublereal) lw, work[1].i = 0.;
#line 413 "zlamswlq.f"
    return 0;

/*     End of ZLAMSWLQ */

} /* zlamswlq_ */

