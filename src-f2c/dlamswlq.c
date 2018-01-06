#line 1 "dlamswlq.f"
/* dlamswlq.f -- translated by f2c (version 20100827).
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

#line 1 "dlamswlq.f"
/* Table of constant values */

static integer c__0 = 0;


/*  Definition: */
/*  =========== */

/*      SUBROUTINE DLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
/*     $                LDT, C, LDC, WORK, LWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER         SIDE, TRANS */
/*      INTEGER           INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC */
/*     .. */
/*     .. Array Arguments .. */
/*      DOUBLE        A( LDA, * ), WORK( * ), C(LDC, * ), */
/*     $                  T( LDT, * ) */
/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLAMQRTS overwrites the general real M-by-N matrix C with */
/* > */
/* > */
/* >                    SIDE = 'L'     SIDE = 'R' */
/* >    TRANS = 'N':      Q * C          C * Q */
/* >    TRANS = 'T':      Q**T * C       C * Q**T */
/* >    where Q is a real orthogonal matrix defined as the product of blocked */
/* >    elementary reflectors computed by short wide LQ */
/* >    factorization (DLASWLQ) */
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
/* >          The number of rows of the matrix C.  M >=0. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
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
/* >          T is DOUBLE PRECISION array, dimension */
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
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
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
/* >         (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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
/* Subroutine */ int dlamswlq_(char *side, char *trans, integer *m, integer *
	n, integer *k, integer *mb, integer *nb, doublereal *a, integer *lda, 
	doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	doublereal *work, integer *lwork, integer *info, ftnlen side_len, 
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
    extern /* Subroutine */ int dgemlqt_(char *, char *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), dtpmlqt_(char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 236 "dlamswlq.f"
    /* Parameter adjustments */
#line 236 "dlamswlq.f"
    a_dim1 = *lda;
#line 236 "dlamswlq.f"
    a_offset = 1 + a_dim1;
#line 236 "dlamswlq.f"
    a -= a_offset;
#line 236 "dlamswlq.f"
    t_dim1 = *ldt;
#line 236 "dlamswlq.f"
    t_offset = 1 + t_dim1;
#line 236 "dlamswlq.f"
    t -= t_offset;
#line 236 "dlamswlq.f"
    c_dim1 = *ldc;
#line 236 "dlamswlq.f"
    c_offset = 1 + c_dim1;
#line 236 "dlamswlq.f"
    c__ -= c_offset;
#line 236 "dlamswlq.f"
    --work;
#line 236 "dlamswlq.f"

#line 236 "dlamswlq.f"
    /* Function Body */
#line 236 "dlamswlq.f"
    lquery = *lwork < 0;
#line 237 "dlamswlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 238 "dlamswlq.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 239 "dlamswlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 240 "dlamswlq.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 241 "dlamswlq.f"
    if (left) {
#line 242 "dlamswlq.f"
	lw = *n * *mb;
#line 243 "dlamswlq.f"
    } else {
#line 244 "dlamswlq.f"
	lw = *m * *mb;
#line 245 "dlamswlq.f"
    }

#line 247 "dlamswlq.f"
    *info = 0;
#line 248 "dlamswlq.f"
    if (! left && ! right) {
#line 249 "dlamswlq.f"
	*info = -1;
#line 250 "dlamswlq.f"
    } else if (! tran && ! notran) {
#line 251 "dlamswlq.f"
	*info = -2;
#line 252 "dlamswlq.f"
    } else if (*m < 0) {
#line 253 "dlamswlq.f"
	*info = -3;
#line 254 "dlamswlq.f"
    } else if (*n < 0) {
#line 255 "dlamswlq.f"
	*info = -4;
#line 256 "dlamswlq.f"
    } else if (*k < 0) {
#line 257 "dlamswlq.f"
	*info = -5;
#line 258 "dlamswlq.f"
    } else if (*lda < max(1,*k)) {
#line 259 "dlamswlq.f"
	*info = -9;
#line 260 "dlamswlq.f"
    } else if (*ldt < max(1,*mb)) {
#line 261 "dlamswlq.f"
	*info = -11;
#line 262 "dlamswlq.f"
    } else if (*ldc < max(1,*m)) {
#line 263 "dlamswlq.f"
	*info = -13;
#line 264 "dlamswlq.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 265 "dlamswlq.f"
	*info = -15;
#line 266 "dlamswlq.f"
    }

#line 268 "dlamswlq.f"
    if (*info != 0) {
#line 269 "dlamswlq.f"
	i__1 = -(*info);
#line 269 "dlamswlq.f"
	xerbla_("DLAMSWLQ", &i__1, (ftnlen)8);
#line 270 "dlamswlq.f"
	work[1] = (doublereal) lw;
#line 271 "dlamswlq.f"
	return 0;
#line 272 "dlamswlq.f"
    } else if (lquery) {
#line 273 "dlamswlq.f"
	work[1] = (doublereal) lw;
#line 274 "dlamswlq.f"
	return 0;
#line 275 "dlamswlq.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 279 "dlamswlq.f"
    i__1 = min(*m,*n);
#line 279 "dlamswlq.f"
    if (min(i__1,*k) == 0) {
#line 280 "dlamswlq.f"
	return 0;
#line 281 "dlamswlq.f"
    }

/* Computing MAX */
#line 283 "dlamswlq.f"
    i__1 = max(*m,*n);
#line 283 "dlamswlq.f"
    if (*nb <= *k || *nb >= max(i__1,*k)) {
#line 284 "dlamswlq.f"
	dgemlqt_(side, trans, m, n, k, mb, &a[a_offset], lda, &t[t_offset], 
		ldt, &c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)
		1);
#line 286 "dlamswlq.f"
	return 0;
#line 287 "dlamswlq.f"
    }

#line 289 "dlamswlq.f"
    if (left && tran) {

/*         Multiply Q to the last block of C */

#line 293 "dlamswlq.f"
	kk = (*m - *k) % (*nb - *k);
#line 294 "dlamswlq.f"
	ctr = (*m - *k) / (*nb - *k);
#line 295 "dlamswlq.f"
	if (kk > 0) {
#line 296 "dlamswlq.f"
	    ii = *m - kk + 1;
#line 297 "dlamswlq.f"
	    dtpmlqt_("L", "T", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 300 "dlamswlq.f"
	} else {
#line 301 "dlamswlq.f"
	    ii = *m + 1;
#line 302 "dlamswlq.f"
	}

#line 304 "dlamswlq.f"
	i__1 = *nb + 1;
#line 304 "dlamswlq.f"
	i__2 = -(*nb - *k);
#line 304 "dlamswlq.f"
	for (i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+NB) */

#line 308 "dlamswlq.f"
	    --ctr;
#line 309 "dlamswlq.f"
	    i__3 = *nb - *k;
#line 309 "dlamswlq.f"
	    dtpmlqt_("L", "T", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 313 "dlamswlq.f"
	}

/*         Multiply Q to the first block of C (1:M,1:NB) */

#line 317 "dlamswlq.f"
	dgemlqt_("L", "T", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 320 "dlamswlq.f"
    } else if (left && notran) {

/*         Multiply Q to the first block of C */

#line 324 "dlamswlq.f"
	kk = (*m - *k) % (*nb - *k);
#line 325 "dlamswlq.f"
	ii = *m - kk + 1;
#line 326 "dlamswlq.f"
	ctr = 1;
#line 327 "dlamswlq.f"
	dgemlqt_("L", "N", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 330 "dlamswlq.f"
	i__2 = ii - *nb + *k;
#line 330 "dlamswlq.f"
	i__1 = *nb - *k;
#line 330 "dlamswlq.f"
	for (i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (I:I+NB,1:N) */

#line 334 "dlamswlq.f"
	    i__3 = *nb - *k;
#line 334 "dlamswlq.f"
	    dtpmlqt_("L", "N", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 337 "dlamswlq.f"
	    ++ctr;

#line 339 "dlamswlq.f"
	}
#line 340 "dlamswlq.f"
	if (ii <= *m) {

/*         Multiply Q to the last block of C */

#line 344 "dlamswlq.f"
	    dtpmlqt_("L", "N", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 348 "dlamswlq.f"
	}

#line 350 "dlamswlq.f"
    } else if (right && notran) {

/*         Multiply Q to the last block of C */

#line 354 "dlamswlq.f"
	kk = (*n - *k) % (*nb - *k);
#line 355 "dlamswlq.f"
	ctr = (*n - *k) / (*nb - *k);
#line 356 "dlamswlq.f"
	if (kk > 0) {
#line 357 "dlamswlq.f"
	    ii = *n - kk + 1;
#line 358 "dlamswlq.f"
	    dtpmlqt_("R", "N", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 361 "dlamswlq.f"
	} else {
#line 362 "dlamswlq.f"
	    ii = *n + 1;
#line 363 "dlamswlq.f"
	}

#line 365 "dlamswlq.f"
	i__1 = *nb + 1;
#line 365 "dlamswlq.f"
	i__2 = -(*nb - *k);
#line 365 "dlamswlq.f"
	for (i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 369 "dlamswlq.f"
	    --ctr;
#line 370 "dlamswlq.f"
	    i__3 = *nb - *k;
#line 370 "dlamswlq.f"
	    dtpmlqt_("R", "N", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (
		    ftnlen)1, (ftnlen)1);

#line 374 "dlamswlq.f"
	}

/*         Multiply Q to the first block of C (1:M,1:MB) */

#line 378 "dlamswlq.f"
	dgemlqt_("R", "N", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 381 "dlamswlq.f"
    } else if (right && tran) {

/*       Multiply Q to the first block of C */

#line 385 "dlamswlq.f"
	kk = (*n - *k) % (*nb - *k);
#line 386 "dlamswlq.f"
	ctr = 1;
#line 387 "dlamswlq.f"
	ii = *n - kk + 1;
#line 388 "dlamswlq.f"
	dgemlqt_("R", "T", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 391 "dlamswlq.f"
	i__2 = ii - *nb + *k;
#line 391 "dlamswlq.f"
	i__1 = *nb - *k;
#line 391 "dlamswlq.f"
	for (i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 395 "dlamswlq.f"
	    i__3 = *nb - *k;
#line 395 "dlamswlq.f"
	    dtpmlqt_("R", "T", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (
		    ftnlen)1, (ftnlen)1);
#line 398 "dlamswlq.f"
	    ++ctr;

#line 400 "dlamswlq.f"
	}
#line 401 "dlamswlq.f"
	if (ii <= *n) {

/*       Multiply Q to the last block of C */

#line 405 "dlamswlq.f"
	    dtpmlqt_("R", "T", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);

#line 409 "dlamswlq.f"
	}

#line 411 "dlamswlq.f"
    }

#line 413 "dlamswlq.f"
    work[1] = (doublereal) lw;
#line 414 "dlamswlq.f"
    return 0;

/*     End of DLAMSWLQ */

} /* dlamswlq_ */

