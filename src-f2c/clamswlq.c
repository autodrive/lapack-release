#line 1 "clamswlq.f"
/* clamswlq.f -- translated by f2c (version 20100827).
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

#line 1 "clamswlq.f"
/* Table of constant values */

static integer c__0 = 0;


/*  Definition: */
/*  =========== */

/*      SUBROUTINE CLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
/*     $                LDT, C, LDC, WORK, LWORK, INFO ) */


/*     .. Scalar Arguments .. */
/*      CHARACTER         SIDE, TRANS */
/*      INTEGER           INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC */
/*     .. */
/*     .. Array Arguments .. */
/*      COMPLEX        A( LDA, * ), WORK( * ), C(LDC, * ), */
/*     $                  T( LDT, * ) */
/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLAMQRTS overwrites the general real M-by-N matrix C with */
/* > */
/* > */
/* >                    SIDE = 'L'     SIDE = 'R' */
/* >    TRANS = 'N':      Q * C          C * Q */
/* >    TRANS = 'T':      Q**H * C       C * Q**H */
/* >    where Q is a real orthogonal matrix defined as the product of blocked */
/* >    elementary reflectors computed by short wide LQ */
/* >    factorization (CLASWLQ) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Transpose, apply Q**H. */
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
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the blocked */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CLASWLQ in the first k rows of its array argument A. */
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
/* >          T is COMPLEX array, dimension */
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
/* >          C is COMPLEX array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
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
/* >         (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* Subroutine */ int clamswlq_(char *side, char *trans, integer *m, integer *
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
    extern /* Subroutine */ int cgemlqt_(char *, char *, integer *, integer *,
	     integer *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, doublecomplex *, integer 
	    *, ftnlen, ftnlen), ctpmlqt_(char *, char *, integer *, integer *,
	     integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);


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

#line 236 "clamswlq.f"
    /* Parameter adjustments */
#line 236 "clamswlq.f"
    a_dim1 = *lda;
#line 236 "clamswlq.f"
    a_offset = 1 + a_dim1;
#line 236 "clamswlq.f"
    a -= a_offset;
#line 236 "clamswlq.f"
    t_dim1 = *ldt;
#line 236 "clamswlq.f"
    t_offset = 1 + t_dim1;
#line 236 "clamswlq.f"
    t -= t_offset;
#line 236 "clamswlq.f"
    c_dim1 = *ldc;
#line 236 "clamswlq.f"
    c_offset = 1 + c_dim1;
#line 236 "clamswlq.f"
    c__ -= c_offset;
#line 236 "clamswlq.f"
    --work;
#line 236 "clamswlq.f"

#line 236 "clamswlq.f"
    /* Function Body */
#line 236 "clamswlq.f"
    lquery = *lwork < 0;
#line 237 "clamswlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 238 "clamswlq.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 239 "clamswlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 240 "clamswlq.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 241 "clamswlq.f"
    if (left) {
#line 242 "clamswlq.f"
	lw = *n * *mb;
#line 243 "clamswlq.f"
    } else {
#line 244 "clamswlq.f"
	lw = *m * *mb;
#line 245 "clamswlq.f"
    }

#line 247 "clamswlq.f"
    *info = 0;
#line 248 "clamswlq.f"
    if (! left && ! right) {
#line 249 "clamswlq.f"
	*info = -1;
#line 250 "clamswlq.f"
    } else if (! tran && ! notran) {
#line 251 "clamswlq.f"
	*info = -2;
#line 252 "clamswlq.f"
    } else if (*m < 0) {
#line 253 "clamswlq.f"
	*info = -3;
#line 254 "clamswlq.f"
    } else if (*n < 0) {
#line 255 "clamswlq.f"
	*info = -4;
#line 256 "clamswlq.f"
    } else if (*k < 0) {
#line 257 "clamswlq.f"
	*info = -5;
#line 258 "clamswlq.f"
    } else if (*lda < max(1,*k)) {
#line 259 "clamswlq.f"
	*info = -9;
#line 260 "clamswlq.f"
    } else if (*ldt < max(1,*mb)) {
#line 261 "clamswlq.f"
	*info = -11;
#line 262 "clamswlq.f"
    } else if (*ldc < max(1,*m)) {
#line 263 "clamswlq.f"
	*info = -13;
#line 264 "clamswlq.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 265 "clamswlq.f"
	*info = -15;
#line 266 "clamswlq.f"
    }

#line 268 "clamswlq.f"
    if (*info != 0) {
#line 269 "clamswlq.f"
	i__1 = -(*info);
#line 269 "clamswlq.f"
	xerbla_("CLAMSWLQ", &i__1, (ftnlen)8);
#line 270 "clamswlq.f"
	work[1].r = (doublereal) lw, work[1].i = 0.;
#line 271 "clamswlq.f"
	return 0;
#line 272 "clamswlq.f"
    } else if (lquery) {
#line 273 "clamswlq.f"
	work[1].r = (doublereal) lw, work[1].i = 0.;
#line 274 "clamswlq.f"
	return 0;
#line 275 "clamswlq.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 279 "clamswlq.f"
    i__1 = min(*m,*n);
#line 279 "clamswlq.f"
    if (min(i__1,*k) == 0) {
#line 280 "clamswlq.f"
	return 0;
#line 281 "clamswlq.f"
    }

/* Computing MAX */
#line 283 "clamswlq.f"
    i__1 = max(*m,*n);
#line 283 "clamswlq.f"
    if (*nb <= *k || *nb >= max(i__1,*k)) {
#line 284 "clamswlq.f"
	cgemlqt_(side, trans, m, n, k, mb, &a[a_offset], lda, &t[t_offset], 
		ldt, &c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)
		1);
#line 286 "clamswlq.f"
	return 0;
#line 287 "clamswlq.f"
    }

#line 289 "clamswlq.f"
    if (left && tran) {

/*         Multiply Q to the last block of C */

#line 293 "clamswlq.f"
	kk = (*m - *k) % (*nb - *k);
#line 294 "clamswlq.f"
	ctr = (*m - *k) / (*nb - *k);
#line 295 "clamswlq.f"
	if (kk > 0) {
#line 296 "clamswlq.f"
	    ii = *m - kk + 1;
#line 297 "clamswlq.f"
	    ctpmlqt_("L", "C", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 300 "clamswlq.f"
	} else {
#line 301 "clamswlq.f"
	    ii = *m + 1;
#line 302 "clamswlq.f"
	}

#line 304 "clamswlq.f"
	i__1 = *nb + 1;
#line 304 "clamswlq.f"
	i__2 = -(*nb - *k);
#line 304 "clamswlq.f"
	for (i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+NB) */

#line 308 "clamswlq.f"
	    --ctr;
#line 309 "clamswlq.f"
	    i__3 = *nb - *k;
#line 309 "clamswlq.f"
	    ctpmlqt_("L", "C", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 313 "clamswlq.f"
	}

/*         Multiply Q to the first block of C (1:M,1:NB) */

#line 317 "clamswlq.f"
	cgemlqt_("L", "C", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 320 "clamswlq.f"
    } else if (left && notran) {

/*         Multiply Q to the first block of C */

#line 324 "clamswlq.f"
	kk = (*m - *k) % (*nb - *k);
#line 325 "clamswlq.f"
	ii = *m - kk + 1;
#line 326 "clamswlq.f"
	ctr = 1;
#line 327 "clamswlq.f"
	cgemlqt_("L", "N", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 330 "clamswlq.f"
	i__2 = ii - *nb + *k;
#line 330 "clamswlq.f"
	i__1 = *nb - *k;
#line 330 "clamswlq.f"
	for (i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (I:I+NB,1:N) */

#line 334 "clamswlq.f"
	    i__3 = *nb - *k;
#line 334 "clamswlq.f"
	    ctpmlqt_("L", "N", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 337 "clamswlq.f"
	    ++ctr;

#line 339 "clamswlq.f"
	}
#line 340 "clamswlq.f"
	if (ii <= *m) {

/*         Multiply Q to the last block of C */

#line 344 "clamswlq.f"
	    ctpmlqt_("L", "N", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 348 "clamswlq.f"
	}

#line 350 "clamswlq.f"
    } else if (right && notran) {

/*         Multiply Q to the last block of C */

#line 354 "clamswlq.f"
	kk = (*n - *k) % (*nb - *k);
#line 355 "clamswlq.f"
	ctr = (*n - *k) / (*nb - *k);
#line 356 "clamswlq.f"
	if (kk > 0) {
#line 357 "clamswlq.f"
	    ii = *n - kk + 1;
#line 358 "clamswlq.f"
	    ctpmlqt_("R", "N", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 361 "clamswlq.f"
	} else {
#line 362 "clamswlq.f"
	    ii = *n + 1;
#line 363 "clamswlq.f"
	}

#line 365 "clamswlq.f"
	i__1 = *nb + 1;
#line 365 "clamswlq.f"
	i__2 = -(*nb - *k);
#line 365 "clamswlq.f"
	for (i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 369 "clamswlq.f"
	    --ctr;
#line 370 "clamswlq.f"
	    i__3 = *nb - *k;
#line 370 "clamswlq.f"
	    ctpmlqt_("R", "N", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (
		    ftnlen)1, (ftnlen)1);
#line 373 "clamswlq.f"
	}

/*         Multiply Q to the first block of C (1:M,1:MB) */

#line 377 "clamswlq.f"
	cgemlqt_("R", "N", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 380 "clamswlq.f"
    } else if (right && tran) {

/*       Multiply Q to the first block of C */

#line 384 "clamswlq.f"
	kk = (*n - *k) % (*nb - *k);
#line 385 "clamswlq.f"
	ii = *n - kk + 1;
#line 386 "clamswlq.f"
	ctr = 1;
#line 387 "clamswlq.f"
	cgemlqt_("R", "C", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 390 "clamswlq.f"
	i__2 = ii - *nb + *k;
#line 390 "clamswlq.f"
	i__1 = *nb - *k;
#line 390 "clamswlq.f"
	for (i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 394 "clamswlq.f"
	    i__3 = *nb - *k;
#line 394 "clamswlq.f"
	    ctpmlqt_("R", "C", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], 
		    lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 
		    1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (
		    ftnlen)1, (ftnlen)1);
#line 397 "clamswlq.f"
	    ++ctr;

#line 399 "clamswlq.f"
	}
#line 400 "clamswlq.f"
	if (ii <= *n) {

/*       Multiply Q to the last block of C */

#line 404 "clamswlq.f"
	    ctpmlqt_("R", "C", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
		     &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);

#line 408 "clamswlq.f"
	}

#line 410 "clamswlq.f"
    }

#line 412 "clamswlq.f"
    work[1].r = (doublereal) lw, work[1].i = 0.;
#line 413 "clamswlq.f"
    return 0;

/*     End of CLAMSWLQ */

} /* clamswlq_ */

